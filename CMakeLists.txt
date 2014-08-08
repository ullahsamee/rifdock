cmake_minimum_required (VERSION 2.6)

add_subdirectory(external/gmock)

include_directories("external/gmock/gtest/include")
include_directories("external/gmock/include")
set (EXTRA_LIBS ${EXTRA_LIBS} gmock boost_system boost_filesystem boost_serialization boost_iostreams z)
#set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)

include_directories("src")

FILE(GLOB L1 "src/*.cc")
FILE(GLOB L2 "src/*/*.cc")
FILE(GLOB L3 "src/*/*/*.cc")
FILE(GLOB L4 "src/*/*/*/*.cc")
FILE(GLOB L5 "src/*/*/*/*/*.cc")

#SET( CMAKE_CXX_FLAGS         "-O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
SET( CMAKE_CXX_FLAGS         "-O3 -DDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -DDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_CXX_FLAGS         "-O0 -DDEBUG -ffast-math" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O0 -DDEBUG -ffast-math" )
#SET( CMAKE_CXX_FLAGS         "-O0 -DDEBUG" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O0 -DDEBUG" )
#SET( CMAKE_CXX_FLAGS         "-std=c++11 -ftemplate-backtrace-limit=0" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -ftemplate-backtrace-limit=0" )

project (main_test)
add_executable(main_test src/main/main_test.cc ${L1} ${L2} ${L3} ${L4} ${L5} )
target_link_libraries(main_test ${EXTRA_LIBS} )


project (quick_test)
add_executable(quick_test src/main/quick_test.cpp)
target_link_libraries(quick_test ${EXTRA_LIBS})
