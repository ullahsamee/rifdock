
message( "entering apps/rosetta/python")

project(pysetta)

include(CheckCXXCompilerFlag)

# Set a default build configuration if none is specified.
# 'MinSizeRel' produces the smallest binaries
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'MinSizeRel' as none was specified.")
  set(CMAKE_BUILD_TYPE MinSizeRel CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)


if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Enable link time optimization and set the default symbol
  # visibility to hidden (very important to obtain small binaries)
  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    # Default symbol visibility
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")

    # Check for Link Time Optimization support
    CHECK_CXX_COMPILER_FLAG("-flto" HAS_LTO_FLAG)
    if (HAS_LTO_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    endif()
  endif()
endif()

# Include path for Python header files
include_directories(../../../external/pybind11/include)

find_path( PYTHON_INCLUDE_DIR python3.5m/Python.h )
find_library( PYTHON_LIB_FILE python3 )
get_filename_component( PYTHON_LIB_DIR ${PYTHON_LIB_FILE} DIRECTORY )

include_directories( ${PYTHON_INCLUDE_DIR}/python3.5m )
link_directories( ${PYTHON_LIB_DIR} )

install(DIRECTORY pysetta/ DESTINATION lib/python3/pysetta FILES_MATCHING PATTERN "*.py")

message( "using python include: " ${PYTHON_INCLUDE_DIR}/python3.5m )
message( "using python lib:     " ${PYTHON_LIB_DIR} )

add_subdirectory(pysetta)
add_subdirectory(test)

