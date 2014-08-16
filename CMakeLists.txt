cmake_minimum_required (VERSION 2.6)

FILE(GLOB L1 "src/*.cc")
FILE(GLOB L2 "src/*/*.cc")
FILE(GLOB L3 "src/*/*/*.cc")
FILE(GLOB L4 "src/*/*/*/*.cc")
FILE(GLOB L5 "src/*/*/*/*/*.cc")

SET (CMAKE_C_COMPILER  "/usr/bin/clang")

#SET( CMAKE_CXX_FLAGS         "-O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_CXX_FLAGS         "-O0 -DDEBUG -ffast-math" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O0 -DDEBUG -ffast-math" )
#SET( CMAKE_CXX_FLAGS         "-std=c++11 -ftemplate-backtrace-limit=0" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -ftemplate-backtrace-limit=0" )

SET( CMAKE_CXX_FLAGS         "-O3 -DDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -DDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_CXX_FLAGS         "-O3 -DDEBUG -fno-vectorize -fno-slp-vectorize -DEIGEN_DONT_VECTORIZE" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O3 -DDEBUG -fno-vectorize -fno-slp-vectorize -DEIGEN_DONT_VECTORIZE" )
#SET( CMAKE_CXX_FLAGS         "-O0 -DDEBUG" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-O0 -DDEBUG" )

add_subdirectory(external/gmock)
include_directories("external/gmock/gtest/include")
include_directories("external/gmock/include")
include_directories("/usr/include/eigen3")


# shitty "data" files to compile faster
SET_SOURCE_FILES_PROPERTIES( src/objective/rosetta/EtableParams_init.cc PROPERTIES COMPILE_FLAGS -O0 )
SET_SOURCE_FILES_PROPERTIES( src/rosetta/atype_map.cc PROPERTIES COMPILE_FLAGS -O0 )

list(APPEND EXTRA_LIBS gmock boost_system boost_filesystem boost_serialization boost_iostreams z boost_timer)
list(APPEND CMAKE_LIBRARY_PATH /usr/local/lib)


include_directories("src")
link_directories(${CMAKE_LIBRARY_PATH})

project (main_test)
add_executable(main_test apps/main_test.cc ${L1} ${L2} ${L3} ${L4} ${L5} )
target_link_libraries(main_test ${EXTRA_LIBS} )


project (quick_test)
add_executable(quick_test apps/quick_test.cc src/objective/rosetta/EtableParams_init.cc src/rosetta/atype_map.cc )
target_link_libraries(quick_test ${EXTRA_LIBS})
