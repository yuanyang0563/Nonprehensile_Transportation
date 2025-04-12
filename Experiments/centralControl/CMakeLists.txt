cmake_minimum_required(VERSION 3.5)

project(transposrt_experiments)

option(CXX "enable C++ compilation" ON)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)

include_directories(include)

find_package(Boost REQUIRED COMPONENTS thread system)
include_directories(${Boost_INCLUDE_DIRS})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

find_package(VISP REQUIRED)
include_directories(${VISP_INCLUDE_DIRS})

find_package(catkin REQUIRED COMPONENTS roscpp sensor_msgs std_msgs)
include_directories(${catkin_INCLUDE_DIRS})



add_executable(pid src/pid.cpp src/common.cpp)
target_link_libraries(pid ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES})

add_executable(mpc src/mpc.cpp src/common.cpp)
target_link_libraries(mpc ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES} optimized ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
