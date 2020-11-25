#!/usr/bin/python
import numpy as np
import math
import os
import sys
import string
import time
import random
import tf
import struct



import argparse

import rospy

import baxter_interface

from baxter_interface import CHECK_VERSION

from sensor_msgs.msg import Image

from baxter_synchronisation.msg import Pos

from geometry_msgs.msg import (
    PoseStamped,
    Pose,
    Point,
    Quaternion,
)
from std_msgs.msg import Header

from baxter_core_msgs.srv import (
    SolvePositionIK,
    SolvePositionIKRequest,
)


S_Const = 100# lower limit of bandpass
I_Const = 5#upper limit of bandpass
Rate = 0.03


class slave:
	def __init__(self, start, speed, kay, scaling=0.05):# takes initial conditions and time constant
		
		self.state = start# state in x, y and z
		self.feedback = np.zeros(3)# feedback coming back from the plant
		self.power = 1# final model defaults to zero, although master feedback could force this
		
		self.state = np.array(start)# set the state
		
		self.rate = speed# set the rate

		self.DC = np.matrix(np.zeros(3))
		
		self.K = kay# this is the coupling constant

		self.scaling = scaling

		self.output = np.zeros(3)


	def update(self, mast):

		
		self.state = np.array(self.state) -self.rate*(np.array(self.state))+ self.K*(self.scaling*np.array(mast) - self.output*self.scaling)

		self.DC = self.DC + (self.rate/10)*(mast-self.output)
		self.output = self.state# + self.DC
		print self.output
		return self.output





class buff:# buffers values and performs noise filtering
	def __init__(self, lo, hi = 5):
		self.lo = lo
		self.hi = hi
		if lo > 0 and lo>=hi:
			self.wait = np.zeros((3,self.lo), dtype = np.float64)# maximum index is self.length-1
		else:
			self.wait = np.zeros((3,1), dtype = np.float64)
		self.counter = 1# counts number of iterations
		self.tell = np.zeros((3,1), dtype = np.float64)# vector of new coordinates

	def filterize(self, vPoint):
		self.tell = vPoint# pass new coordinates
		outer = np.zeros((3,1), dtype = np.float64)
		if self.lo ==0:# zero length indicates no filtering desired
			return np.array([[vPoint[0]], [vPoint[1]], [vPoint[2]]], dtype = np.float64)
		else:

			if not(np.isnan(self.tell).any()):# if the spherical transformation has not produced NaN (anything+NaN = NaN)
				self.wait[:, 1:self.lo] = self.wait[:, 0:self.lo-1]# shuffle up buffer - indices ARE CORRECT
				self.wait[:,0] = self.tell# add newest member #CHECK IF ALIGNED

				summed = np.sum(self.wait, axis=1, dtype = np.float64)/self.hi
				outer = np.reshape(summed,(3,1))/self.lo# in principle, rolling average should be updated after being used
			
			

			return outer# return high-pass filtered signal




class Master:
	def __init__(self, start, speed):# takes initial conditions and time constant
		self.state = start# state in x, y and z
		self.mast = np.zeros(3)# state of the master
		self.feedback = np.zeros(3)# feedback coming back from the plant
		self.power = 1

		self.state = start# set the state
		self.rate = speed# set the rate

		self.metro = 0# "metronome"
		self.mFeed = 0

		#self.C = 0.1

		# self.a = 0.15# Rossler constants
		# self.b = 0.2
		# self.c = 10

	def update(self):
		self.mFeed = self.mFeed + self.rate*0.2*(self.metro)
		self.metro = self.metro + self.rate*0.2*(self.power-self.mFeed)
		
		self.feedback[0] = self.feedback[0] + self.rate*(self.state[0])
		self.feedback[1] = self.feedback[1] + self.rate*(self.state[1])
		self.feedback[2] = self.feedback[2] + self.rate*(self.state[2])

		q = np.random.random()*0
		self.state[0] = self.state[0] + self.rate*(1)*(self.power+q-self.feedback[0])#+np.random.random_sample()
		self.state[1] = self.state[1] + self.rate*(1)*(self.power+q-self.feedback[1])
		self.state[2] = self.state[2] + self.rate*(1)*(self.power+q-self.feedback[2])
		

	def polarize(self,point):
		radius = np.sqrt(np.power(point[0],2) + np.power(point[1],2) + np.power(point[2],2))# radial (extension) coordinate
		theta = np.arccos(point[2]/radius)# rotation around horizontal
		azimuth = np.arctan(point[1]/point[0])# rotation around vertical
		poles = np.array([radius, theta, azimuth])# put 'em together

		return poles

	def comZip(self, limb, osc, gain, vgain, side):# packages the states of the oscillators into a velocity command for the robot
		current = limb.joint_angles()# get the angles at each joint
		vels = limb.joint_velocities()# get the velocities

		diff = gain*np.array([osc[0], osc[2]-osc[1], 1*osc[0], osc[2]-osc[1], 0, 0, 0]) - vgain*np.array([vels[side+'_s0'], vels[side+'_s1'], vels[side+'_e0'], vels[side+'_e1'], 0, 0, 0])# this has been altered to take spherical error values directly, bypassing the dynamics

		labels = [side+'_s0', side+'_s1', side+'_e0', side+'_e1', side+'_w0', side+'_w1', side+'_w2']# labels for each joint
		comm = {labels[0]:diff[0], labels[1]:diff[1], labels[2]:diff[2], labels[3]:diff[3], labels[4]:diff[4], labels[5]:diff[5], labels[6]:diff[6]}

		return comm


class rMaster:
	def __init__(self, wavelength, mag, t = 0.01):# takes initial conditions, 'wavelength' and magnitude multiplier
		self.state = np.zeros(3)# state in x, y and z
		self.R = np.zeros(3)
		self.wav = wavelength
		self.mag = mag
		self.counter = 0

		self.rate = t


	def update(self):# update state
		if self.counter < self.wav:
			self.counter +=1
		else:
			self.counter = 0
			self.R[0] = (0.5-np.random.random())*self.mag
		self.state += self.rate*(self.R-self.state)
		return self.state

		


def main():
	global S_Const
	global I_Const
	global rate

	arg_fmt = argparse.RawDescriptionHelpFormatter
	parser = argparse.ArgumentParser(formatter_class=arg_fmt,
                                  	description=main.__doc__)
	required = parser.add_argument_group('required arguments')# this section is adding the required arguments
	parser.add_argument(
    	"-W", "--Wavelength", type=int, default=0,# number of time steps each pseudo-square wave will last
    	help=("Wavelength")
	)
	parser.add_argument(
    	"-S", "--Seed", type=int, default=0,# random number seed
    	help=("Random number seed")
	)
	args = parser.parse_args(rospy.myargv()[1:])

	wav = args.Wavelength
	Seed = args.Seed


	starter = [0, 0, 0]

	master_targ = rMaster(wav,5)# create pseudorandom master process

	np.random.seed(Seed)# seed pseudorandom number generator


	rospy.init_node('master_unit')# initialise node

	Master_pub = rospy.Publisher('Master_state', Pos, queue_size=1)# master publisher
	master_buff = buff(S_Const, I_Const)# buffer for the master value

	right_pub = rospy.Publisher('right_master', Pos, queue_size=1)# I think I want to drop most missed messages, rather than let them pile up
	right_buff = buff(S_Const)# buffer for the master value

	left_pub = rospy.Publisher('left_master', Pos, queue_size=1)# I think I want to drop most missed messages, rather than let them pile up
	left_buff = buff(10)# buffer for the master value

	Right = baxter_interface.Limb('right')# Right Arm!!
	Left = baxter_interface.Limb('left')# Left Arm!!

	#offset = [0.58, -0.3, 0.1]
	pose = Left.endpoint_pose()
	offset = np.array([pose['position'].x, pose['position'].y, pose['position'].z])

	scale = 1# scale for testing

	
	counter = 0
	r = rospy.Rate(100)
	while not rospy.is_shutdown():

		master_targ.update()# update master state
		
		smoothM = master_buff.filterize(master_targ.state)# filter output


		poz = np.array([smoothM[0], smoothM[1], smoothM[2]])# place in message-compatible array
		

		posy = np.multiply(poz,[[scale], [0], [0]])# only x-axis is used for the master
		print posy[0]# print master state
		

		Master_pub.publish(Pos(posy[0]+0.7, posy[1], posy[2]))#add appropriate DC offset and publish values


		r.sleep()
	#serv.close()
	return 0

if __name__ == '__main__':
    sys.exit(main())

