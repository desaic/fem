PROJECT (MicrostructuresAnalysis)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FILE(GLOB LIB_SRC "*.cpp")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../include")

# utilities
include_directories(${PROJECT_SOURCE_DIR}/../Utilities/include)

# femlib (for Eigen)
include_directories(${PROJECT_SOURCE_DIR}/../femlib/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../femlib ${CMAKE_CURRENT_BINARY_DIR}/femlib)

# sampler utilities
include_directories(${PROJECT_SOURCE_DIR}/../SamplerUtilities/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../SamplerUtilities ${CMAKE_CURRENT_BINARY_DIR}/SamplerUtilities)

add_library(MicrostructuresAnalysis ${LIB_SRC})

set(MicrostructuresAnalysis_LIBS MicrostructuresAnalysis ${SamplerUtilities_LIBS} PARENT_SCOPE)

SET(LIB_TYPE STATIC)

