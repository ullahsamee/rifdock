cmake_minimum_required (VERSION 2.6)

#SET( CMAKE_CXX_FLAGS         "-std=c++11 -w -O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -w -O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_CXX_FLAGS         "-std=c++11 -w -O0 -DDEBUG -ffast-math" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -w -O0 -DDEBUG -ffast-math" )
#SET( CMAKE_CXX_FLAGS         "-std=c++11 -w -ftemplate-backtrace-limit=0" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -w -ftemplate-backtrace-limit=0" )

SET( CMAKE_CXX_FLAGS          "-std=c++11 -w -O3 -DDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
SET( CMAKE_EXE_LINKER_FLAGS   "-std=c++11 -w -O3 -DDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )

#SET( CMAKE_CXX_FLAGS         "-std=c++11 -w -O3 -DNDEBUG -fno-vectorize -fno-slp-vectorize -DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_ALIGN" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -w -O3 -DNDEBUG -fno-vectorize -fno-slp-vectorize -DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_ALIGN" )

#SET( CMAKE_CXX_FLAGS         "-std=c++11 -w -O0 -DDEBUG" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -w -O0 -DDEBUG" )

add_subdirectory(external/gmock)
include_directories("external/gmock/gtest/include")
include_directories("external/gmock/include")
include_directories("/usr/local/include/eigen3")
include_directories("/usr/local/include")

list( APPEND LIB_CC src/rosetta/score/EtableParams_init.lib.cc src/rosetta/atype_map.lib.cc src/chemical/ligand_factory.lib.cc )

# shitty "data" files to compile faster
SET_SOURCE_FILES_PROPERTIES( ${LIB_CC} PROPERTIES COMPILE_FLAGS -O0 )

FILE(GLOB L1 "src/*.cc")
FILE(GLOB L2 "src/*/*.cc")
FILE(GLOB L3 "src/*/*/*.cc")
FILE(GLOB L4 "src/*/*/*/*.cc")
FILE(GLOB L5 "src/*/*/*/*/*.cc")

list(APPEND EXTRA_LIBS gmock boost_system boost_filesystem boost_iostreams z boost_timer)
list(APPEND CMAKE_LIBRARY_PATH /usr/local/lib)


include_directories("src")
link_directories(${CMAKE_LIBRARY_PATH})

project (main_test)
add_executable(main_test apps/main_test.cc ${L1} ${L2} ${L3} ${L4} ${L5} )
target_link_libraries(main_test ${EXTRA_LIBS} )


project (quick_test)
add_executable(quick_test apps/quick_test.cc ${LIB_CC} )
target_link_libraries(quick_test ${EXTRA_LIBS})
