cmake_minimum_required (VERSION 3.3)

set( PYTHON_VERSION 2 )
set( CXX_STD 11 )
message( "build type .................... " ${CMAKE_BUILD_TYPE} )
message( "python version ................ " ${PYTHON_VERSION} )
message( "c++ standard .................. " "c++${CXX_STD}" )


if( ${CXX_STD} MATCHES 11 )
  set( CXX_DIR cxx11 )
  add_definitions( -DCXX11 )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -w" )
elseif( ${CXX_STD} MATCHES 14 )
	set( CXX_DIR cxx14 )
  add_definitions( -DCXX11 )
  add_definitions( -DCXX14 )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -w" )
else()
	message( FATAL_ERROR "bad CXX_STD ${CXX_STD}" )
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)
if( ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG )
	set( ROSETTA_BIN $ENV{CMAKE_ROSETTA_PATH}/source/cmake/build_${CXX_DIR}_omp_debug )
else()
	set( ROSETTA_BIN $ENV{CMAKE_ROSETTA_PATH}/source/cmake/build_${CXX_DIR}_omp )
endif()  
message( "path to rosetta binaries ...... " ${ROSETTA_BIN} )



set( CMAKE_INSTALL_PREFIX $ENV{SCHEME_INSTALL_PATH} )
message( "install prefix ................ " ${CMAKE_INSTALL_PREFIX} )

#####################################################################################

if(${PYTHON_VERSION} MATCHES 2)
  find_path( PYTHON_INCLUDE_DIR python2.7/Python.h )
  find_library( PYTHON_LIB_FILE python2.7 )
  get_filename_component( PYTHON_LIB_DIR ${PYTHON_LIB_FILE} DIRECTORY )
  include_directories( ${PYTHON_INCLUDE_DIR}/python2.7 )
  link_directories( ${PYTHON_LIB_DIR} )
  install(DIRECTORY pysetta/ DESTINATION lib/python2/pysetta FILES_MATCHING PATTERN "*.py")
  message( "using python include .......... " ${PYTHON_INCLUDE_DIR}/python2.7 )
  message( "using python lib .............. " ${PYTHON_LIB_DIR} )
  message( "install python libs to ........ " ${CMAKE_INSTALL_PREFIX}/lib/python2 )
else()
  find_path( PYTHON_INCLUDE_DIR python3.5m/Python.h )
  find_library( PYTHON_LIB_FILE python3 )
  get_filename_component( PYTHON_LIB_DIR ${PYTHON_LIB_FILE} DIRECTORY )
  include_directories( ${PYTHON_INCLUDE_DIR}/python3.5m )
  link_directories( ${PYTHON_LIB_DIR} )
  install(DIRECTORY pysetta/ DESTINATION lib/python3/pysetta FILES_MATCHING PATTERN "*.py")
  message( "using python include .......... " ${PYTHON_INCLUDE_DIR}/python3.5m )
  message( "using python lib .............. " ${PYTHON_LIB_DIR} )
  message( "install python libs to ........ " ${CMAKE_INSTALL_PREFIX}/lib/python3 )
endif()



######################### external dependencies ####################################

find_path( EIGEN3_INCLUDE_DIR eigen3/Eigen )
if( ${EIGEN3_INCLUDE_DIR} MATCHES EIGEN3_INCLUDE_DIR-NOTFOUND )
	message( FATAL_ERROR "NO eigen3/Eigen found")
endif()
include_directories( ${EIGEN3_INCLUDE_DIR}/eigen3 )

add_subdirectory(external/gmock)
include_directories("external/gmock/gtest/include")
include_directories("external/gmock/include")
list( APPEND EXTRA_LIBS gmock )

include_directories( "external/include" )
include_directories( "schemelib" )
include_directories( "/usr/local/include" )
link_directories( "/usr/local/lib" )
link_directories( "/usr/local/lib64" )

####################################################################################

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

# # remove sublime text turds
# message( "removing ._* turds from " ${CMAKE_SOURCE_DIR} )
# execute_process( COMMAND find ${CMAKE_SOURCE_DIR} -regex .*[.]_.* -exec rm {} \; )


add_subdirectory(schemelib)
add_subdirectory(apps)

