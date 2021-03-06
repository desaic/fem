cmake_minimum_required(VERSION 2.8.8) # version 2.8.8 is required for VTK

project("SamplesGenerator")

set(QHULL_SOURCE_DIR $ENV{QHULL_SRC_PATH})
set(QHULL_BIN_DIR $ENV{QHULL_BUILD_PATH})
link_directories($ENV{QHULL_BUILD_PATH})

if(DEFINED ENV{PARDISO_LIB_PATH})
	set(IS_PARDISO_AVAILABLE ON)
endif()

file(GLOB src "*.cpp")
file(GLOB headers "*.h")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# utilities
include_directories(${PROJECT_SOURCE_DIR}/../Utilities/include)
if (NOT TARGET Utilities)
add_subdirectory(${PROJECT_SOURCE_DIR}/../Utilities ${CMAKE_CURRENT_BINARY_DIR}/Utilities)
endif()

# sampler utilities
include_directories(${PROJECT_SOURCE_DIR}/../SamplerUtilities/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../SamplerUtilities ${CMAKE_CURRENT_BINARY_DIR}/SamplerUtilities)

# microstructure family generation
include_directories(${PROJECT_SOURCE_DIR}/../MicrostructuresAnalysis)
add_subdirectory(${PROJECT_SOURCE_DIR}/../MicrostructuresAnalysis ${CMAKE_CURRENT_BINARY_DIR}/MicrostructuresAnalysis)

# math
#include_directories(${PROJECT_SOURCE_DIR}/../vecmath/include)
#if (NOT TARGET vecmath)
#    add_subdirectory(${PROJECT_SOURCE_DIR}/../vecmath ${CMAKE_CURRENT_BINARY_DIR}/vecmath)
#endif() 

# femlib
include_directories(${PROJECT_SOURCE_DIR}/../femlib/include)
if (NOT TARGET FEMLib)
	add_subdirectory(${PROJECT_SOURCE_DIR}/../femlib ${CMAKE_CURRENT_BINARY_DIR}/femlib)
endif()

# continuous material optimization
include_directories(${PROJECT_SOURCE_DIR}/../MaterialOpt/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../MaterialOpt ${CMAKE_CURRENT_BINARY_DIR}/MaterialOpt)

# stepper newton
include_directories(${PROJECT_SOURCE_DIR}/../StepperCPU/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../StepperCPU ${CMAKE_CURRENT_BINARY_DIR}/StepperCPU)
link_directories(${PROJECT_SOURCE_DIR}/../lib)  #for blas

# multigrid solver
SET(USE_AMGCL ON)
SET(USE_CUDA ON)
SET(USE_NLOPT ON)
add_subdirectory(${PROJECT_SOURCE_DIR}/../topoOptLib ${CMAKE_CURRENT_BINARY_DIR}/topoOptLib)
include_directories(${topoOptLib_INCLUDE_DIR}) 
#add_definitions(${topoOptLib_DEFS})

if (IS_PARDISO_AVAILABLE)
	link_directories($ENV{PARDISO_LIB_PATH})
endif()

# stepper cuda
#include_directories(${PROJECT_SOURCE_DIR}/../StepperCuda/include)
#add_subdirectory(${PROJECT_SOURCE_DIR}/../StepperCuda ${CMAKE_CURRENT_BINARY_DIR}/StepperCuda)

add_executable (SamplesGenerator ${src})

target_link_libraries(SamplesGenerator ${VTK_LIBRARIES} FEMLib StepperCPU libopenblas ${Utilities_LIBS} MaterialOpt ${SamplerUtilities_LIBS} ${MicrostructuresAnalysis_LIBS} ${topoOptLib_LIBS})
if (IS_PARDISO_AVAILABLE)
	target_link_libraries(SamplesGenerator libpardiso500-WIN-X86-64)
endif()

# copy dlls to build directory
option(COPY_DLL_TO_BUILD_DIR "copy dlls to build directory" OFF)
if (COPY_DLL_TO_BUILD_DIR) 
	add_custom_command(TARGET SamplesGenerator POST_BUILD
					   COMMAND set "var=${PROJECT_SOURCE_DIR}"
					   COMMAND set "var=%var:/=\\%"
					   COMMAND set "var1=${PROJECT_BINARY_DIR}"
					   COMMAND set "var1=%var1:/=\\%"
					   
					   COMMAND copy "%var%\\..\\dll\\*.dll" "%var1%\\$<CONFIGURATION>"  
					   COMMENT "copying dlls …"
				   )
	if (IS_PARDISO_AVAILABLE)
		add_custom_command(TARGET SamplesGenerator POST_BUILD
						   COMMAND set "var=$ENV{PARDISO_LIB_PATH}"
						   COMMAND set "var=%var:/=\\%"
						   COMMAND set "var1=${PROJECT_BINARY_DIR}"
					       COMMAND set "var1=%var1:/=\\%"
						   COMMAND copy "%var%\\*.dll" "%var1%\\$<CONFIGURATION>"   
						   COMMENT "copying dlls …"
				   )
	endif()
endif()
