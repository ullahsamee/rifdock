cmake_minimum_required (VERSION 2.6)

execute_process( COMMAND find ${CMAKE_SOURCE_DIR} -regex .*[.]_.* -exec rm {} \; )

set( CMAKE_BUILD_TYPE Release )

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -w" )
add_definitions( -DCXX11 ) # these don't work....
add_definitions( -DCXX14 )
add_definitions( -std=c++14 )

######################### external dependencies ####################################

find_path( EIGEN3_INCLUDE_DIR eigen3/Eigen )
include_directories( ${EIGEN3_INCLUDE_DIR}/eigen3 )
MESSAGE( "add includ dir for eigen3: " ${EIGEN3_INCLUDE_DIR}/eigen3 )

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

add_subdirectory(test)
add_subdirectory(apps)
