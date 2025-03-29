cmake_minimum_required(VERSION 3.5)

project(nonprehensile_mpc)

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

find_package(catkin REQUIRED COMPONENTS roscpp std_msgs sensor_msgs)
include_directories(${catkin_INCLUDE_DIRS})



add_executable(pid_1arm src/pid_1arm.cpp src/common.cpp)
target_link_libraries(pid_1arm ${catkin_LIBRARIES})

add_executable(mpc_1arm src/mpc_1arm.cpp src/common.cpp)
target_link_libraries(mpc_1arm ${catkin_LIBRARIES} optimized ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})

add_executable(pid_2arms src/pid_2arms.cpp src/common.cpp)
target_link_libraries(pid_2arms ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES})

add_executable(mpc_2arms src/mpc_2arms.cpp src/common.cpp)
target_link_libraries(mpc_2arms ${catkin_LIBRARIES} optimized ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
