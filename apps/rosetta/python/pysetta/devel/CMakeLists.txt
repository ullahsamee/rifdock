# Create the binding library
add_library(_devel SHARED devel.cpp )
target_link_libraries( _devel python3 ${ALL_ROSETTA_LIBS} )
set_target_properties( _devel  PROPERTIES PREFIX "")

file( GLOB PY_FILES "*.py" )
file( COPY ${PY_FILES} DESTINATION . )


# Don't add a 'lib' prefix to the shared library

if (UNIX)
  # It's quite common to have multiple copies of the same Python version
  # installed on one's system. E.g.: one copy from the OS and another copy
  # that's statically linked into an application like Blender or Maya.
  # If we link our plugin library against the OS Python here and import it
  # into Blender or Maya later on, this will cause segfaults when multiple
  # conflicting Python instances are active at the same time (even when they
  # are of the same version).

  # Windows is not affected by this issue since it handles DLL imports
  # differently. The solution for Linux and Mac OS is simple: we just don't
  # link against the Python library. The resulting shared library will have
  # missing symbols, but that's perfectly fine -- they will be resolved at
  # import time.

  # .SO file extension on Linux/Mac OS
  set_target_properties(_devel  PROPERTIES SUFFIX ".so")

  # Strip unnecessary sections of the binary on Linux/Mac OS
  if(APPLE)
    set_target_properties(_devel  PROPERTIES MACOSX_RPATH ".")
    set_target_properties(_devel  PROPERTIES LINK_FLAGS "-undefined dynamic_lookup ")
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET _devel  POST_BUILD
                         COMMAND strip -u -r ${PROJECT_BINARY_DIR}/pysetta/devel/_devel.so)
    endif()
  else()
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET _devel  POST_BUILD
                         COMMAND strip ${PROJECT_BINARY_DIR}/pysetta/devel/_devel.so)
    endif()
  endif()
endif()

