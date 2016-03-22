PROJECT (Utilities)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FILE(GLOB LIB_SRC "src/*.cpp")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../include")

# femlib
include_directories(${PROJECT_SOURCE_DIR}/../femlib/include)
if (NOT TARGET FEMLib)
	add_subdirectory(${PROJECT_SOURCE_DIR}/../femlib ${CMAKE_CURRENT_BINARY_DIR}/femlib)
endif()

add_library(Utilities ${LIB_SRC})

set(Utilities_LIBS Utilities PARENT_SCOPE)

SET(LIB_TYPE STATIC)

