cmake_minimum_required(VERSION 3.5)

project(transport_experiments)

find_package(Boost REQUIRED COMPONENTS thread system)

find_package(VISP REQUIRED)

find_package(catkin REQUIRED COMPONENTS roscpp std_msgs)

include_directories(
${Boost_INCLUDE_DIRS}
${VISP_INCLUDE_DIRS}
${catkin_INCLUDE_DIRS}
)

add_executable(leftCamera leftCamera.cpp)
target_link_libraries(leftCamera ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES})
