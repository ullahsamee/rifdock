cmake_minimum_required (VERSION 2.6)

add_subdirectory(external/gmock)

include_directories("external/gmock/gtest/include")
include_directories("external/gmock/include")
set (EXTRA_LIBS ${EXTRA_LIBS} gmock)

include_directories("src")

FILE(GLOB UTIL "src/util/*.cc")
FILE(GLOB UTIL_META "src/util/meta/*.cc")
FILE(GLOB UTIL_CONT "src/util/container/*.cc")
FILE(GLOB NEST "src/nest/*.cc")
FILE(GLOB OBJ "src/objective/*.cc")
FILE(GLOB NEST_MAP "src/nest/maps/*.cc")
FILE(GLOB KIN "src/kinematics/*.cc")

#SET( CMAKE_CXX_FLAGS         "-O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_CXX_FLAGS         "-O3 -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
SET( CMAKE_CXX_FLAGS         "-O0 -DDEBUG -ffast-math" )
SET( CMAKE_EXE_LINKER_FLAGS  "-O0 -DDEBUG -ffast-math" )
#SET( CMAKE_CXX_FLAGS         "-O0 " )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O0 " )

project (main_test)
add_executable(main_test src/main/main_test.cc ${UTIL} ${NEST} ${OBJ} ${NEST_MAP} ${UTIL_META} ${KIN} ${UTIL_CONT})
target_link_libraries(main_test ${EXTRA_LIBS} )


project (quick_test)
add_executable(quick_test src/main/quick_test.cc)
target_link_libraries(quick_test ${EXTRA_LIBS})
