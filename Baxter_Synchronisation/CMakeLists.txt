cmake_minimum_required(VERSION 2.8.3)
project(baxter_synchronisation)

find_package(catkin
  REQUIRED
  COMPONENTS
  rospy
  xacro
  actionlib
  control_msgs
  trajectory_msgs
  cv_bridge
  dynamic_reconfigure
  baxter_core_msgs
  baxter_interface
  geometry_msgs
  sensor_msgs
  std_msgs
  message_generation
)

catkin_python_setup()

add_message_files(FILES
  Pos.msg
)

generate_messages(DEPENDENCIES
  geometry_msgs
  sensor_msgs
  std_msgs
)

catkin_package(
  CATKIN_DEPENDS
  rospy
  xacro
  actionlib
  sensor_msgs
  control_msgs
  trajectory_msgs
  cv_bridge
  dynamic_reconfigure
  baxter_core_msgs
  baxter_interface
  message_runtime
)

install(
  DIRECTORY scripts/
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  USE_SOURCE_PERMISSIONS
)
