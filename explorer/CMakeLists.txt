cmake_minimum_required(VERSION 2.8.8) # version 2.8.8 is required for VTK

project("explorer")

#set(QHULL_SOURCE_DIR $ENV{QHULL_SRC_PATH})
#set(QHULL_BIN_DIR $ENV{QHULL_BUILD_PATH})

file(GLOB src "*.cpp")
file(GLOB headers "*.h")
file(GLOB ui_sources "*.ui")

#find VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#find Qt:
find_package(Qt4 COMPONENTS QtCore QtGui REQUIRED)
include(${QT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

#QT4_ADD_RESOURCES(QtResources explorer.qrc)

QT4_WRAP_CPP(QtMocs ${headers})
QT4_WRAP_UI(QtUIs ${ui_sources})

# utilities
include_directories(${PROJECT_SOURCE_DIR}/../Utilities/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/../Utilities ${CMAKE_CURRENT_BINARY_DIR}/Utilities)
link_directories($ENV{QHULL_BUILD_PATH})

# math
#include_directories(${PROJECT_SOURCE_DIR}/../vecmath/include)
#if (NOT TARGET vecmath)
#	add_subdirectory(${PROJECT_SOURCE_DIR}/../vecmath ${CMAKE_CURRENT_BINARY_DIR}/vecmath)
#endif() 

# femlib
include_directories(${PROJECT_SOURCE_DIR}/../femlib/include)
if (NOT TARGET FEMLib)
	add_subdirectory(${PROJECT_SOURCE_DIR}/../femlib ${CMAKE_CURRENT_BINARY_DIR}/femlib)
endif()

add_executable (explorer ${src} ${QtMocs} ${QtUIs})

target_link_libraries(explorer ${QT_QTGUI_LIBRARY} ${QT_QTCORE_LIBRARY} ${VTK_LIBRARIES} FEMLib ${Utilities_LIBS})

# copy VTK dll to runtime directory
option(COPY_DLL_TO_BUILD_DIR "copy dlls to build directory" ON)
if (COPY_DLL_TO_BUILD_DIR) 
	add_custom_command(TARGET explorer POST_BUILD
					   COMMAND set "var1=${PROJECT_BINARY_DIR}"
					   COMMAND set "var1=%var1:/=\\%"
				       COMMAND COPY "$ENV{VTK_BUILD_PATH}\\bin\\$<CONFIGURATION>\\*.dll" "%var1%\\$<CONFIGURATION>" 
                       COMMENT "copying dlls …"
    )
endif()
