PROJECT (ConcatArray)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#find_package(Eigen3)
#include_directories (${EIGEN3_INCLUDE_DIR})
add_definitions( -DPARDISO_AVAILABLE)
ADD_DEFINITIONS(--std=c++0x)
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../../Utilities/include")
FILE(GLOB EXE_SRC "*.cpp")

ADD_EXECUTABLE(concat ${EXE_SRC})
