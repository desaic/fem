cmake_minimum_required(VERSION 2.8.8) # version 2.8.8 is required for VTK

project("SamplesGenerator")

set(QHULL_SOURCE_DIR $ENV{QHULL_SRC_PATH})
set(QHULL_BIN_DIR $ENV{QHULL_BUILD_PATH})

file(GLOB src "*.cpp")
file(GLOB headers "*.h")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# utilities
include_directories(${PROJECT_SOURCE_DIR}/../Utilities/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../Utilities ${CMAKE_CURRENT_BINARY_DIR}/Utilities)
link_directories($ENV{QHULL_BUILD_PATH})

# math
include_directories(${PROJECT_SOURCE_DIR}/../vecmath/include)
if (NOT TARGET vecmath)
    add_subdirectory(${PROJECT_SOURCE_DIR}/../vecmath ${CMAKE_CURRENT_BINARY_DIR}/vecmath)
endif() 

# femlib
include_directories(${PROJECT_SOURCE_DIR}/../femlib/include)
if (NOT TARGET vecmath)
	add_subdirectory(${PROJECT_SOURCE_DIR}/../femlib ${CMAKE_CURRENT_BINARY_DIR}/femlib)
endif()

# stepper newton
include_directories(${PROJECT_SOURCE_DIR}/../StepperCPU/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../StepperCPU ${CMAKE_CURRENT_BINARY_DIR}/StepperCPU)
link_directories(${PROJECT_SOURCE_DIR}/../lib)  #for blas

# stepper cuda
#include_directories(${PROJECT_SOURCE_DIR}/../StepperCuda/include)
#add_subdirectory(${PROJECT_SOURCE_DIR}/../StepperCuda ${CMAKE_CURRENT_BINARY_DIR}/StepperCuda)

add_executable (SamplesGenerator ${src})

target_link_libraries(SamplesGenerator ${VTK_LIBRARIES} vecmath FEMLib StepperCPU libopenblas ${Utilities_LIBS})

# copy dlls to build directory
option(COPY_DLL_TO_BUILD_DIR "copy dlls to build directory" ON)
if (COPY_DLL_TO_BUILD_DIR) 
	add_custom_command(TARGET SamplesGenerator POST_BUILD
					   COMMAND set "var=${PROJECT_SOURCE_DIR}"
					   COMMAND set "var=%var:/=\\%"
					   COMMAND echo "%var%"
					   COMMAND copy "%var%\\..\\dll\\*.dll" "C:\\Users\\melina\\Documents\\Projects\\fem\\SamplesGenerator\\build\\$<CONFIGURATION>" 
					   COMMENT "copying dlls …"
				   )
endif()
