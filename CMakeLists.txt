cmake_minimum_required (VERSION 3.3)

message( "build type " ${CMAKE_BUILD_TYPE} )

# remove sublime text turds
message( "removing ._* turds from " ${CMAKE_SOURCE_DIR} )
execute_process( COMMAND find ${CMAKE_SOURCE_DIR} -regex .*[.]_.* -exec rm {} \; )

set( CMAKE_INSTALL_PREFIX $ENV{SCHEME_INSTALL_PATH} )
message( "CMAKE_INSTALL_PREFIX is " ${CMAKE_INSTALL_PREFIX} )

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -w" )
add_definitions( -DCXX11 )
add_definitions( -DCXX14 )

######################### external dependencies ####################################

find_path( EIGEN3_INCLUDE_DIR eigen3/Eigen )
include_directories( ${EIGEN3_INCLUDE_DIR}/eigen3 )
MESSAGE( "add include dir for eigen3: " ${EIGEN3_INCLUDE_DIR}/eigen3 )

add_subdirectory(external/gmock)
include_directories("external/gmock/gtest/include")
include_directories("external/gmock/include")
list( APPEND EXTRA_LIBS gmock )

include_directories( "external/include" )
include_directories( "schemelib" )
include_directories( "/usr/local/include" )
link_directories( "/usr/local/lib" )
link_directories( "/usr/local/lib64" )

#####################################################################################

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")


add_subdirectory(schemelib)
add_subdirectory(apps)
