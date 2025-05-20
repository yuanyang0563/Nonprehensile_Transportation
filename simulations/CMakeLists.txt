cmake_minimum_required(VERSION 3.5)

project(nonprehensile_transportation)

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

find_package(catkin REQUIRED COMPONENTS roscpp std_msgs sensor_msgs geometry_msgs)
include_directories(${catkin_INCLUDE_DIRS})

add_executable(pbc_dual src/pbc_dual.cpp)
target_link_libraries(pbc_dual ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES})

add_executable(pbc_left src/pbc_left.cpp)
target_link_libraries(pbc_left ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES})

add_executable(pbc_right src/pbc_right.cpp)
target_link_libraries(pbc_right ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES})

add_executable(mpc_dual src/mpc_dual.cpp)
target_link_libraries(mpc_dual ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES} optimized ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})

add_executable(mpc_left src/mpc_left.cpp)
target_link_libraries(mpc_left ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES} optimized ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})

add_executable(mpc_right src/mpc_right.cpp)
target_link_libraries(mpc_right ${Boost_LIBRARIES} ${VISP_LIBRARIES} ${catkin_LIBRARIES} optimized ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
