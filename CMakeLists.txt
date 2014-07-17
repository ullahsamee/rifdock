cmake_minimum_required (VERSION 2.6)

add_subdirectory(external/gmock)

include_directories("external/gmock/gtest/include")
include_directories("external/gmock/include")
set (EXTRA_LIBS ${EXTRA_LIBS} gmock)

include_directories("src")

FILE(GLOB UTIL_CC "src/util/*.cc")
FILE(GLOB UTIL_MPL_CC "src/util/meta/*.cc")
FILE(GLOB NEST_CC "src/nest/*.cc")
FILE(GLOB NEST_CC "src/objective/*.cc")
FILE(GLOB NEST_MAP_CC "src/nest/maps/*.cc")

#SET( CMAKE_CXX_FLAGS         "-O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_CXX_FLAGS         "-O3 -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
SET( CMAKE_CXX_FLAGS         "-O3 -ffast-math" )
SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -ffast-math" )
#SET( CMAKE_CXX_FLAGS         "-O0 " )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O0 " )

project (main_test)
add_executable(main_test src/main/main_test.cc ${UTIL_CC} ${NEST_CC} ${NEST_MAP_CC} ${UTIL_MPL_CC})
target_link_libraries(main_test ${EXTRA_LIBS} )


project (quick_test)
add_executable(quick_test src/main/quick_test.cc)
target_link_libraries(quick_test ${EXTRA_LIBS})
