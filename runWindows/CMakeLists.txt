cmake_minimum_required(VERSION 2.8.8) # version 2.8.8 is required for VTK

project("RunWindows")

set(QHULL_SOURCE_DIR $ENV{QHULL_SRC_PATH})
set(QHULL_BIN_DIR $ENV{QHULL_BUILD_PATH})
message(${QHULL_BIN_DIR})
if(DEFINED ENV{PARDISO_LIB_PATH})
	set(IS_PARDISO_AVAILABLE ON)
endif()

file(GLOB src "*.cpp")
file(GLOB headers "*.h")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR}/../include)
# utilities
include_directories(${PROJECT_SOURCE_DIR}/../Utilities/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../Utilities ${CMAKE_CURRENT_BINARY_DIR}/Utilities)
#if(CMAKE_BUILD_TYPE STREQUAL "Release")
#  link_directories(${QHULL_BIN_DIR}/Release)
#else()
#  link_directories(${QHULL_BIN_DIR}/Debug)
#endif()
link_directories($ENV{QHULL_BUILD_PATH})

# femlib
include_directories(${PROJECT_SOURCE_DIR}/../femlib/include)
if (NOT TARGET FEMLib)
	add_subdirectory(${PROJECT_SOURCE_DIR}/../femlib ${CMAKE_CURRENT_BINARY_DIR}/femlib)
endif()

# viewer
include_directories(${PROJECT_SOURCE_DIR}/../viewer/include)
if (NOT TARGET viewer)
	add_subdirectory(${PROJECT_SOURCE_DIR}/../viewer ${CMAKE_CURRENT_BINARY_DIR}/viewer)
endif()

# stepper newton
include_directories(${PROJECT_SOURCE_DIR}/../StepperCPU/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../StepperCPU ${CMAKE_CURRENT_BINARY_DIR}/StepperCPU)
link_directories(${PROJECT_SOURCE_DIR}/../lib)  #for blas

include_directories(${PROJECT_SOURCE_DIR}/../MaterialOpt/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../MaterialOpt ${CMAKE_CURRENT_BINARY_DIR}/MaterialOpt)

if (IS_PARDISO_AVAILABLE)
	link_directories($ENV{PARDISO_LIB_PATH})
	message($ENV{PARDISO_LIB_PATH})
endif()

# stepper cuda
#include_directories(${PROJECT_SOURCE_DIR}/../StepperCuda/include)
#add_subdirectory(${PROJECT_SOURCE_DIR}/../StepperCuda ${CMAKE_CURRENT_BINARY_DIR}/StepperCuda)

add_executable (Run ${src})

target_link_libraries(Run FEMLib MaterialOpt Utilities StepperCPU libopenblas FEMViewer ${Utilities_LIBS})
if(MSVC)
  target_link_libraries(Run glfw3 opengl32 glu32)
endif()

if (IS_PARDISO_AVAILABLE)
	target_link_libraries(Run libpardiso500-WIN-X86-64)
endif()

# copy dlls to build directory
option(COPY_DLL_TO_BUILD_DIR "copy dlls to build directory" ON)
if (COPY_DLL_TO_BUILD_DIR) 
	add_custom_command(TARGET Run POST_BUILD
					   COMMAND set "var=${PROJECT_SOURCE_DIR}"
					   COMMAND set "var=%var:/=\\%"
					   COMMAND set "var1=${PROJECT_BINARY_DIR}"
					   COMMAND set "var1=%var1:/=\\%"
					   
					   COMMAND copy "%var%\\..\\dll\\*.dll" "%var1%\\$<CONFIGURATION>" 
					   COMMENT "copying dlls ..."
				   )
	if (IS_PARDISO_AVAILABLE)
		add_custom_command(TARGET Run POST_BUILD
		               COMMAND set "var1=${PROJECT_BINARY_DIR}"
					   COMMAND set "var1=%var1:/=\\%"
					   COMMAND set "var2=$ENV{PARDISO_LIB_PATH}"
					   COMMAND set "var2=%var2:/=\\%"
					   COMMAND copy "%var2%\\*.dll" "%var1%\\$<CONFIGURATION>" 
					   COMMAND copy "%var2%\\*.dll" "%var1%\\$<CONFIGURATION>" 
					   COMMENT "copying pardiso dlls ..."
				   )
	endif()
endif()
