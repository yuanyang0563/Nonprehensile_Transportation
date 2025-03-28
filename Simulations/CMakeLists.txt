cmake_minimum_required(VERSION 3.5)

project(nonprehensile_mpc)

option(CXX "enable C++ compilation" ON)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

set(COPPELIASIM_INCLUDE_DIR $ENV{COPPELIASIM_ROOT_DIR}/programming/include)
list(APPEND CMAKE_MODULE_PATH ${COPPELIASIM_INCLUDE_DIR}/cmake)
find_package(CoppeliaSim 4.1.0.0 REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_package(cppzmq REQUIRED)
include_directories(include)

add_library(RemoteAPIClient STATIC include/remote_api/RemoteAPIClient.cpp)
target_compile_definitions(RemoteAPIClient PUBLIC -DSIM_REMOTEAPICLIENT_OBJECTS)
target_link_libraries(RemoteAPIClient PUBLIC cppzmq)

add_executable(mpc_1arm src/mpc_1arm.cpp src/common.cpp)
target_link_libraries(mpc_1arm RemoteAPIClient optimized ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})

add_executable(mpc_2arms src/mpc_2arms.cpp src/common.cpp)
target_link_libraries(mpc_2arms RemoteAPIClient optimized ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
