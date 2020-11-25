#!/usr/bin/python
# alright, this one implements the slave target dynamics and interfaces
# with the robot simulation to produce the anticipated error
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

from baxter_synchronisation.msg import Pos

from geometry_msgs.msg import (
    PoseStamped,
    Pose,
    Point,
    Quaternion,
)
import std_msgs

from baxter_core_msgs.srv import (
    SolvePositionIK,
    SolvePositionIKRequest,
)

CARTESIAN = 0# if using cartesian control

depth = 3# number of joints recruited for movements
spread = 1# number of individual frequencies
downSample = 0# downsampling rate (2 = move every other loop)

S_Const = 0
I_Const = 0#upper limit of bandpass

SUPPRESS = 2# constant that governs how much coupling will be supressed by robot/slave discrepancy

MINS = np.array([-1.7, -2.1, -3, 0, -3, -1.5, -3])# minimum joint rotations in radians
MAXS = np.array([1.7, 1, 3, 2.6, 3, 2, 3])# maximum joint rotations

rightRot = np.array([[-np.sqrt(2)/2, np.sqrt(2)/2, 0],[-np.sqrt(2)/2, -np.sqrt(2)/2, 0],[0, 0, 1]])# rotation from central world frame to right-arm frame

leftRot = np.array([[np.sqrt(2)/2, np.sqrt(2)/2, 0],[-np.sqrt(2)/2, np.sqrt(2)/2, 0],[0, 0, 1]])
rightTran = np.array([-0.278, -0.064, 0])# translation vector from central world frame to right-arm frame
leftTran = np.array([0.278, -0.064, 0])# translation vector from central world frame to left-arm frame
	
class buff:# buffers values and performs noise filtering
	def __init__(self, lo, hi = 5, offset =  np.zeros((3,1), dtype = np.float64)):# length is the length of the buffer, pas is upper limit of passband
		self.lo = lo
		self.hi = hi
		#self.pas = np.zeros(3)# cumulative sum
		if lo > 0 and lo>=hi:
			self.wait = np.zeros((3,self.lo), dtype = np.float64)# maximum index is self.length-1
		else:
			self.wait = np.zeros((3,1), dtype = np.float64)
		self.counter = 1# counts number of iterations
		self.tell = np.zeros((3,1), dtype = np.float64)# vector of new coordinates

		self.offset = offset# instead of high-passing, which appears to damage AS, 'grab' the initial offset and assume that is the DC component
		
	def filterize(self, vPoint):
		self.tell = vPoint# pass new coordinates
		outer = np.zeros((3,1), dtype = np.float64)
		if self.lo ==0:# zero length indicates no filtering desired
			return np.array([[vPoint[0]], [vPoint[1]], [vPoint[2]]], dtype = np.float64)
		else:
			#if(self.counter < self.lo):
			self.counter += 1# increment counter
			if not(np.isnan(self.tell).any()):# if the spherical transformation has not produced NaN (anything+NaN = NaN)
				self.wait[:, 1:self.lo] = self.wait[:, 0:self.lo-1]# shuffle up buffer - indices ARE CORRECT
				self.wait[:,0] = self.tell# add newest member #CHECK IF ALIGNED
		
				#self.pas += self.tell
				pus = np.sum(self.wait[:,0:self.lo], axis=1, dtype = np.float64)#self.pas = self.pas+np.reshape(self.tell,(3,1))# add to accumulator

				summed = np.sum(self.wait[:,0:self.hi], axis=1, dtype = np.float64)/self.hi# problematic, marked for possible removal
				# using precise 
				outer = np.reshape(summed,(3,1))-self.offset#np.reshape(self.pas,(3,1))/self.counter# in principle, rolling average should be updated after being used
			
			
			return outer# return high-pass filtered signal


class armslave:
	def __init__(self, start, speed, kay, begin, side, scaling=0.05, length = 1, pas = 1, thresh = 0):# takes initial conditions and time constant

		global depth
		
		self.state = np.array(start, dtype = np.float64)# initial state in x, y and z
		self.mast = np.array([0.7,0,0], dtype = np.float64)# state of the master (position in x, y and z)
		self.feedback = np.zeros(3, dtype = np.float64)# feedback coming back from the plant
		self.power = 0# intrinsic driving term

		self.multiState = np.zeros((depth,spread), dtype = np.float64)# as many states as joints# becomes np.zeros((depth,spread)) (I think)
		self.multiFeed = np.zeros((depth,spread), dtype = np.float64)# matching number of feedbacks# becomes np.zeros((depth,spread)) (I think)

		self.Cbuff = buff(length, pas)# high pass filter for the coupling term to remove trends in the sensor noise (this is applied here because the noise sources in master and slave are independent)
		self.Dbuff = buff(length,pas)# for DC

		
		self.connect = np.ones((spread,1), dtype = np.float64)# adds up DC components

		self.mrate = np.linspace(speed, speed*spread, spread)# oscillators map directly to joints, so differing frequencies not necessary# this becomes a linspace of spread length

		self.ratio = 0# ratio between oscillation and 'DC component'/relaxation

		self.DC = np.array(start, dtype = np.float64)# DC components for each dimension

		self.mK = np.reshape(np.linspace(kay, kay, spread),(1,spread))# vector of coupling constant, currently all the same# this becomes a linspace of spread length

		self.mix = np.array([[0,0,1],[0,1,0],[0,0,0],[-1,0,0]], dtype = np.float64)#maps from dimensions to joints (JRM) - remember, keep coupling constant at equal power for each slave, including summed across sources

		self.interface = begin# arm interface allows for safe stopping
		
		self.K = kay# this is the coupling constant

		self.scaling = scaling# scale between the oscillator and real metre measurements

		self.Aoffset = np.array([[-0.5], [-1.2], [0]], dtype = np.float64)# hardcoded offset to shift zero-centred oscillators onto true joint midpoints
		self.Arange = np.array([[3], [2.4], [3.2]], dtype = np.float64)# scaling factor for joints (conversion to radians)

		self.side = side# the side being operated on

		self.outPut = np.zeros((depth,1))

		self.latch = 0# has the test begun
		
	def add_master(self, data):# grab the state of the master
		#if data.x:# this just prevents an edge case where the coupling is too large
		if self.mast.all() == 0 & self.latch == 0:# the first time a master signal is received, the offset for the master buffer is set
			self.Dbuff.offset[0] = data.x
			self.Dbuff.offset[1] = data.y
			self.Dbuff.offset[2] = data.z
			if np.sqrt(np.sum(np.power(self.mast-self.DC,2)))<2:# error handling - if the value is immediately large the experiemnt cannot progress
				self.latch =1

		self.mast[0] = data.x# set values of master/target
		
		

		self.mast[1] = data.y

		
		
		self.mast[2] = data.z


		

		print("master")# acknowledge message received

	# primary update step
	def update(self, angles, side, init, setup, smooth = np.array([0,0,0], dtype=np.float64)):# slave pose, master angles, side of master, initial angles of slave, current angles of slave, if in setup mode or not
		global depth
		global CARTESIAN

		global SUPPRESS
		
		
		otherSide = {'right':'left', 'left':'right'}# dictionary that reverses side term

		iAngles = np.array([[init[side+'_s0']], [init[side+'_s1']], [init[side+'_e1']]], dtype = np.float64)# slave initial angles
		nAngles = np.array([[angles[side+'_e1']], [-1*angles[side+'_s1']], [angles[side+'_s0']]], dtype = np.float64)# slave current angles



		discrep = np.sum(np.abs(self.outPut-nAngles))# this is a measure (of sorts) of if the robot is keeping up with its desired behaviour - if not, the strength of the coupling is altered
		newK =  self.mK/(1+SUPPRESS*discrep)# as discrepancy grows, size of coupling is decreased, giving robot time to 'catch up'

		if (self.mast.any() !=0):# if the master signal exists
			unsmoothed = polarize(smooth)
			smoothed =polarize(self.Dbuff.filterize(np.array([self.mast[0],self.mast[1],self.mast[2]], dtype = np.float64)))-polarize(truTrig(self.outPut))# spherical difference term for internal model
		else:# otherwise these tarms are 0
			unsmoothed = np.zeros((3,1), dtype = np.float64)
			smoothed = np.zeros((3,1), dtype = np.float64)
			

		if (np.isnan(smoothed).any()):# check for infinite values (this step should not be necessary, but is left in for safety)
			inter = np.zeros([3,spread], dtype = np.float64)# set all to valid (zero) values

			couple = np.zeros([depth,spread], dtype = np.float64)
			coupleR = np.zeros([depth,spread], dtype = np.float64)
		else:
			couple = np.multiply(np.dot([[0,0,1],[0,1,0],[-1,0,0]],smoothed), newK)# coupling multiplied by K vector and JRM
			print smoothed

		delta = np.multiply(self.power-self.multiFeed, self.mrate[0]) + np.reshape(couple[0:3,0],(3,1))# update state dynamics, using individual time constants and coupling terms
		self.multiState += delta#this is the difference term that must be reflected by  robot's joints
			
		felta = np.multiply(self.multiState, self.mrate[0]) - np.reshape(couple[0:3,0],(3,1))# feedback term for the oscillator - coupled, but not given to the robot
		self.multiFeed += felta
			
		dDelta = np.multiply(-self.DC, self.mrate[0]) + np.reshape(couple[0:3,0],(3,1))# update relaxation system
			
		self.DC += dDelta#DC term represents changes in the master's average position# (mAngles-iAngles)


		
		self.outPut = np.dot(self.DC,self.connect)# switch self.DC for self.multiState if using oscillators

		conventional = np.zeros(depth)# 'conventional' control output
		if (not(np.isnan(couple).any())):# LAST CONFIGURATION TESTED -CHECK IT
			conventional = np.dot(couple,self.connect)
		
		return self.outPut, dDelta, unsmoothed, conventional# return output, filtered and unfiltered difference terms, "conventional" control signal and coupling term




	def shut(self):# resets velocity and stops the robot bouncing around
		labels = [self.side+'_s0', self.side+'_s1', self.side+'_e0', self.side+'_e1', self.side+'_w0', self.side+'_w1', self.side+'_w2']# labels for each joint
		stopper = {labels[0]:0, labels[1]:0, labels[2]:0, labels[3]:0, labels[4]:0, labels[5]:0, labels[6]:0}# zero velocity command
		self.interface.set_joint_velocities(stopper)# stop

	

def translate(osc, side, gain =1):# determines precise joint angles to control robot
	global MAXS
	global MINS

	osc[0] = gain*np.maximum(MINS[0],np.minimum(MAXS[0],osc[0]))# ensure values in range
	osc[1] = gain*np.maximum(MINS[1],np.minimum(MAXS[1],osc[1]))# uses weak assumption of forearm/arm equality
	osc[2] = gain*np.maximum(MINS[3],np.minimum(MAXS[3],osc[2]))
	
	labels = [side+'_s0', side+'_s1', side+'_e0', side+'_e1', side+'_w0', side+'_w1', side+'_w2']# labels for each joint
	flexure = {labels[0]:osc[0], labels[1]:osc[1], labels[2]:0, labels[3]:osc[2], labels[4]:0, labels[5]:0, labels[6]:0}# zip control dictionary

	return flexure

def feedthru(com, side, gain, angles):# translates joint commands
	global MAXS
	global MINS

	gcom = gain*com

	hold = 4# independent large gain for fixed joints
	
	labels = [side+'_s0', side+'_s1', side+'_e0', side+'_e1', side+'_w0', side+'_w1', side+'_w2']# labels for each joint
	flexure = {labels[0]:0, labels[1]:gcom[1], labels[2]:-hold*angles[labels[2]], labels[3]:gcom[2], labels[4]:-hold*angles[labels[4]], labels[5]:-hold*angles[labels[5]], labels[6]:-hold*angles[labels[6]]}# zip control dictionary

	return flexure

def synergies(point):
	point = point.astype('float')
	radius = np.sqrt(np.power(point[0],2) + np.power(point[1],2) + np.power(point[2],2))# radial (extension) coordinate

	elbow = np.pi-np.arccos((0.5-np.power(radius,2))/(0.5))# this is based on the sensor's assumption that forearm and arm are both 0.5

	theta = np.arccos(point[2]/radius)# rotation around horizontal
	azimuth = np.arctan2(point[1],point[0])# rotation around vertical
	poles = np.array([elbow, theta, azimuth], dtype = np.float64)# put 'em together

	return poles
	
def dumbTrig(angles):# broadly, these are the 'forward kinematics' of the master humanoid that Baxter hypothesises is producing the motion to be followed

	global MAXS
	global MINS

	l = np.array([0.37,0.374])# segment lengths of locked arm

	lAngles = [np.maximum(MINS[0],np.minimum(MAXS[0],angles[0])), np.maximum(MINS[1],np.minimum(MAXS[1],angles[1])), np.maximum(MINS[3],np.minimum(MAXS[3],angles[2]))]

	rads = np.sqrt( np.power((l[0]*np.sin(lAngles[1])+l[1]*np.sin(lAngles[1]+lAngles[2])),2) + np.power((l[0]*np.cos(lAngles[1])+l[1]*np.cos(lAngles[1]+lAngles[2])),2) )

	theta = -np.arccos((l[0]*np.sin(lAngles[1])+l[1]*np.sin(lAngles[1]+lAngles[2]))/(rads))

	az = -lAngles[0]+0.8# + 1.7016/2 # azimuth is solely dependent on base joint

	
	return np.array([rads, theta, az])

def truTrig(tAngles):# kinematics of Baxter arm with fixed third angle (credit to University of Ohio mechanical school)
	global MAXS
	global MINS
	global leftRot

	# segment lengths in metres
	l0 = 0.2735
	l1 = 0.069
	l2 = 0.3645
	l3 = 0.069
	l4 = 0.37429
	l5 = 0.01
	l6 = 0.3683

	lh = np.sqrt(np.power(l2,2)+np.power(l3,2))

	angles = np.array([np.maximum(MINS[0],np.minimum(MAXS[0],tAngles[0])), np.maximum(MINS[1],np.minimum(MAXS[1],tAngles[1])), 0, np.maximum(MINS[3],np.minimum(MAXS[3],tAngles[2])), 0, 0, 0])

	c1 = np.cos(angles[0])
	c2 = np.cos(angles[1])
	c24 = np.cos(angles[1]+angles[3])
	c5 = np.cos(angles[4])

	s1 = np.sin(angles[0])
	s2 = np.sin(angles[1])
	s24 = np.sin(angles[1]+angles[3])
	s5 = np.sin(angles[4])

	eX = l1*c1 + lh*c1*c2 + l4*c1*c24 - l5*(s1*s5 + c1*s24*c5)+ l6*(c1*s5 - s1*s24*c5)
	whY = l1*s1 + lh*s1*c2 + l4*s1*c24 + l5*(c1*s5 - s1*s24*c5)- l6*(s1*s5 + c1*s24*c5)
	Zed = -lh*s2 - l4*s24 - l5*c24*c5# notably the z axis upside down compared to what I'm used to	

	return np.array([eX, whY, Zed])


def polarize(point):# spherical transformation
	point = point.astype('float')
	radius = np.sqrt(np.power(point[0],2) + np.power(point[1],2) + np.power(point[2],2))# radial (extension) coordinate
	theta = np.arccos(point[2]/radius)# rotation around horizontal
	azimuth = np.arctan2(point[1],point[0])# rotation around vertical
	poles = np.array([radius, theta, azimuth], dtype = np.float64)# put 'em together

	return poles


	
def main():

	global CARTESIAN

	global S_Const
	global I_Const
	global downSample
	global depth
#-----------------------------------Startup Stuff-----------------------------------------------------------------
	arg_fmt = argparse.RawDescriptionHelpFormatter
	parser = argparse.ArgumentParser(formatter_class=arg_fmt,
                                  	description=main.__doc__)
	required = parser.add_argument_group('required arguments')# this section is adding the required arguments
	parser.add_argument(
    	"-K", "--coupling", type=float, default=3.5,# this applies to the difference between the target and the robot's delayed state
    	help=("Strength of the master/slave coupling")
	)
	parser.add_argument(
    	"-g", "--gain", type=float, default=0.8,# 'gain' within internal model - tuned to agree with the robot's motion
    	help=("Internal controller gain")
	)
	parser.add_argument(
    	"-t", "--time", type=float, default=0.004,# the default is tuned to agree with the robot's own motion, but because it's not _a robot_ the velocity term is much smaller
    	help=("Time constant for internal model")
	)
	parser.add_argument(
	"-T", "--tau", type=int, default=1,# the time delay (defined as an integer multiple of execution steps)
    	help=("Feedback time delay")
	)
	args = parser.parse_args(rospy.myargv()[1:])

	rospy.init_node('slave_unit')# initiate node


	Left = baxter_interface.Limb('left')# arm interfaces
	Right = baxter_interface.Limb('right')

	Ktrue = args.gain# assign startup arguments to variables
	K = args.coupling
	rate = args.time
	tau = args.tau

	left_angles = Left.joint_angles()# acquire left arm joint angles
	
	slave_targ = armslave(np.zeros((3,1)), rate, Ktrue, Left, 'left', 1, S_Const, I_Const, rate)# instantiate the slave dynamics
	initer = Left.joint_angles()# initial joint angles, preserved from before the robot moves
#-----------------------------------------------------------------------------------------------------------------


	dCount = 0# downsampling counter

	Slave_pub = rospy.Publisher('Slave_state', Pos, queue_size=1)# slave state publisher

	Couple_pub = rospy.Publisher('Couple', Pos, queue_size=1)# publisher for coupling term

	Master_sub = rospy.Subscriber('Master_state', Pos, slave_targ.add_master, queue_size = 1)# master state subscriber

	#	Output block - sends time-matched messages containing the state of the master nad the robot's end effector
	m_pub = rospy.Publisher('Mast', Pos, queue_size=1)
	s_pub = rospy.Publisher('Slave', Pos, queue_size=1)


	rospy.on_shutdown(slave_targ.shut)# halt the robot's motion on shutdown


	setup_time = 200# hold in setup mode for 200 steps
	
	r = rospy.Rate(100)# set execution step at 10ms

	buff_length = tau
	delay_buffer = np.zeros((3,buff_length))# slave delay buffer

	Left.move_to_joint_positions(translate([-np.pi/4,-np.pi/4, np.pi/2], 'left'))# set up original position

	jrm = np.array([[0,0,1],[0,1,0],[-1,0,0]], dtype = np.float64)# 3d spherical jrm
	shoulder = np.array([0.069,0.278,0])# shoulder offset

	for q in range(1, setup_time):# pause before beginning
		left_pose = Left.endpoint_pose()# retrieve state of the body (combined output of oscillators)
		l_pose = np.array([left_pose['position'].x, left_pose['position'].y, left_pose['position'].z], dtype = np.float64)# extract positional coordinates
		left_angles = Left.joint_angles()# retrieve current joint angles
		com, vels, unsmoothed, conv = slave_targ.update(left_angles, 'left', initer,1, np.array([left_pose['position'].x, left_pose['position'].y, left_pose['position'].z], dtype = np.float64))# perform update step
		print 'Please stay still'

		cCom = truTrig(com)
		Slave_pub.publish(Pos(cCom[0], cCom[1], cCom[2]))

		stub = polarize(np.array([left_pose['position'].x, left_pose['position'].y, left_pose['position'].z], dtype = np.float64)-shoulder)# output robot end effector state
		s_pub.publish(Pos(stub[0],stub[1],stub[2]))

		delay_buffer[:,1:buff_length] = delay_buffer[:,0:buff_length-1]# update delay buffer
		delay_buffer[:,0] = l_pose


		r.sleep()


	

	print 'Ready to Begin'
	while not rospy.is_shutdown():# main program loop
		
		left_pose = Left.endpoint_pose()# retrieve state of the body (combined output of oscillators)
		l_pose = np.array([left_pose['position'].x, left_pose['position'].y, left_pose['position'].z], dtype = np.float64)

		left_angles = Left.joint_angles()
		right_angles = Right.joint_angles()# do the same for the master arm
		
		left_vels = Left.joint_velocities()
		
		delay_buffer[:,1:buff_length] = delay_buffer[:,0:buff_length-1]
		delay_buffer[:,0] = l_pose

		com, vels, unsmoothed, conv = slave_targ.update(left_angles, 'left', initer,1, l_pose)# update the state of the oscillators
		
		cCom = polarize(truTrig(com))#
		Slave_pub.publish(Pos(cCom[0], cCom[1], cCom[2]))

		if (slave_targ.latch == 1):# if test is live

			differ =  cCom.transpose() -polarize(delay_buffer[:,buff_length-1]-shoulder)# obtain coupling term
			Kay = np.dot(jrm,differ[0,:])#
			

			Left.set_joint_velocities(feedthru(vels*80+K*np.array([[Kay[0]],[Kay[1]],[Kay[2]]]),'left',1, left_angles))# combine parallel control inputs and send command to robot

		else:
			differ = polarize(np.array([slave_targ.mast[0],slave_targ.mast[1],slave_targ.mast[2]], dtype = np.float64))-polarize(delay_buffer[:,0]-shoulder)

			Kay = np.dot(jrm,differ)
			Left.set_joint_velocities(feedthru(K*np.array([[Kay[0]],[Kay[1]],[Kay[2]]]),'left',1, left_angles))

		Couple_pub.publish(Pos(unsmoothed[0],unsmoothed[1],unsmoothed[2]))

		#	Output block
		mMast = polarize(slave_targ.mast)# output current (from robot's perspective) master state
		m_pub.publish(Pos(mMast[0],mMast[1],mMast[2]))

		print left_angles['left_w1']

		stub = polarize(np.array([left_pose['position'].x, left_pose['position'].y, left_pose['position'].z], dtype = np.float64)-shoulder)# output current state of robot's end effector
		s_pub.publish(Pos(stub[0]+0.009, stub[1], stub[2]))
		
		print 'Running'
		r.sleep()
	return 0

if __name__ == '__main__':
    sys.exit(main())
