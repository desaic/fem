PROJECT (MicrostructuresAnalysis)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FILE(GLOB LIB_SRC "*.cpp")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../include")

# utilities
include_directories(${PROJECT_SOURCE_DIR}/../Utilities/include)
if (NOT TARGET Utilities)
add_subdirectory(${PROJECT_SOURCE_DIR}/../Utilities ${CMAKE_CURRENT_BINARY_DIR}/Utilities)
endif()

# femlib (for Eigen)
include_directories(${PROJECT_SOURCE_DIR}/../femlib/include)
if (NOT TARGET FEMLib)
add_subdirectory(${PROJECT_SOURCE_DIR}/../femlib ${CMAKE_CURRENT_BINARY_DIR}/femlib)
endif()

# sampler utilities
include_directories(${PROJECT_SOURCE_DIR}/../SamplerUtilities/include)
if (NOT TARGET SamplerUtilities)
add_subdirectory(${PROJECT_SOURCE_DIR}/../SamplerUtilities ${CMAKE_CURRENT_BINARY_DIR}/SamplerUtilities)
endif()

add_library(MicrostructuresAnalysis ${LIB_SRC})

set(MicrostructuresAnalysis_LIBS MicrostructuresAnalysis ${SamplerUtilities_LIBS} ${Utilities_LIBS} PARENT_SCOPE)

SET(LIB_TYPE STATIC)

