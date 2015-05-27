PROJECT (Utilities)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FILE(GLOB LIB_SRC "src/*.cpp")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/include")

set(QHULL_SOURCE_DIR $ENV{QHULL_SRC_PATH})
set(QHULL_BIN_DIR $ENV{QHULL_BUILD_PATH})

# math
include_directories(${PROJECT_SOURCE_DIR}/../vecmath/include)
if (NOT TARGET vecmath)
	add_subdirectory(${PROJECT_SOURCE_DIR}/../vecmath ${CMAKE_CURRENT_BINARY_DIR}/vecmath)
endif()

# Quickhull
include_directories(${QHULL_SOURCE_DIR} ${QHULL_SOURCE_DIR}/libqhull ${QHULL_SOURCE_DIR}/libqhullcpp)
link_directories(${QHULL_BIN_DIR})
add_definitions(-Dqh_QHpointer)

# femlib
include_directories(${PROJECT_SOURCE_DIR}/../femlib/include)
if (NOT TARGET FEMLib)
	add_subdirectory(${PROJECT_SOURCE_DIR}/../femlib ${CMAKE_CURRENT_BINARY_DIR}/femlib)
endif()

add_library(Utilities ${LIB_SRC})

set(Utilities_LIBS Utilities vecmath qhullcpp qhullstatic_p PARENT_SCOPE)

SET(LIB_TYPE STATIC)

