PROJECT (MaterialOpt)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#find_package(Eigen3)
#include_directories (${EIGEN3_INCLUDE_DIR})
ADD_DEFINITIONS(--std=c++0x)
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../Utilities/include")

FILE(GLOB LIB_SRC "src/*.cpp")
ADD_LIBRARY (MaterialOpt ${LIB_SRC})

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../femlib/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../vecmath/include")
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/../StepperCPU/include")
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include")

# if (WIN32)
    # set(Boost_USE_STATIC_LIBS ON)
# else ()
    # option(
        # USE_STATIC_BOOST_LIBS
        # "Use static versions of Boost libraries"
        # OFF
        # )

    # if (USE_STATIC_BOOST_LIBS)
        # set(Boost_USE_STATIC_LIBS ON)
    # else()
        # set(Boost_USE_STATIC_LIBS OFF)
        # add_definitions(-DBOOST_TEST_DYN_LINK)
    # endif()
# endif ()

# set(BOOST_COMPONENTS ${BOOST_COMPONENTS}
    # program_options
    # system
    # filesystem
    # serialization
    # unit_test_framework
    # )

# find_package(Boost COMPONENTS ${BOOST_COMPONENTS})
# include_directories( ${Boost_INCLUDE_DIRS} )
SET(LIB_TYPE STATIC)