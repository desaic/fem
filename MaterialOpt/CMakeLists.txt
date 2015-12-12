PROJECT (MaterialOpt)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#find_package(Eigen3)
#include_directories (${EIGEN3_INCLUDE_DIR})
add_definitions( -DPARDISO_AVAILABLE)
ADD_DEFINITIONS(--std=c++0x)
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../Utilities/include")

add_subdirectory(../femlib       "${CMAKE_CURRENT_BINARY_DIR}/femlib")
add_subdirectory(../viewer       "${CMAKE_CURRENT_BINARY_DIR}/viewer")
add_subdirectory(../vecmath      "${CMAKE_CURRENT_BINARY_DIR}/vecmath")
add_subdirectory(../StepperCPU   "${CMAKE_CURRENT_BINARY_DIR}/StepperCPU")
add_subdirectory(../StepperIpopt "${CMAKE_CURRENT_BINARY_DIR}/StepperIpopt")

FILE(GLOB LIB_SRC "src/*.cpp")
ADD_LIBRARY (MaterialOpt ${LIB_SRC})

FILE(GLOB EXE_SRC "*.cpp")

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../femlib/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../vecmath/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../viewer/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../StepperCPU/include")

LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/../lib-linux)

INCLUDE_DIRECTORIES(include)
ADD_EXECUTABLE(Run ${EXE_SRC})

TARGET_LINK_LIBRARIES(Run StepperCPU MaterialOpt)
TARGET_LINK_LIBRARIES(Run FEMLib)
TARGET_LINK_LIBRARIES(Run FEMViewer)
TARGET_LINK_LIBRARIES(Run vecmath)
TARGET_LINK_LIBRARIES(Run pardiso)
#TARGET_LINK_LIBRARIES(Run coinmumps)

TARGET_LINK_LIBRARIES(Run openblas)
TARGET_LINK_LIBRARIES(Run gfortran)
TARGET_LINK_LIBRARIES(Run gomp)
TARGET_LINK_LIBRARIES(Run png)
TARGET_LINK_LIBRARIES(Run glfw3 X11  Xxf86vm Xrandr Xi Xinerama Xcursor)
TARGET_LINK_LIBRARIES(Run GL GLU)
TARGET_LINK_LIBRARIES(Run pthread )
