cmake_minimum_required (VERSION 2.6)

#SET( CMAKE_CXX_FLAGS         "-std=c++11 -DCEREAL -w -O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -DCEREAL -w -O3 -DNDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_CXX_FLAGS         "-std=c++11 -DCEREAL -w -O0 -DDEBUG -ffast-math" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -DCEREAL -w -O0 -DDEBUG -ffast-math" )
#SET( CMAKE_CXX_FLAGS         "-std=c++11 -DCEREAL -w -ftemplate-backtrace-limit=0" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -DCEREAL -w -ftemplate-backtrace-limit=0" )

#SET( CMAKE_CXX_FLAGS          "-std=c++11 -DCEREAL -w -O3 -DDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )
#SET( CMAKE_EXE_LINKER_FLAGS   "-std=c++11 -DCEREAL -w -O3 -DDEBUG -ffast-math -mmmx -mfpmath=sse -msse -msse2 -msse3" )

#SET( CMAKE_CXX_FLAGS         "-std=c++11 -DCEREAL -w -O3 -DNDEBUG -fno-vectorize -fno-slp-vectorize -DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_ALIGN" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -DCEREAL -w -O3 -DNDEBUG -fno-vectorize -fno-slp-vectorize -DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_ALIGN" )

# SET( CMAKE_CXX_FLAGS         "-std=c++11 -DCEREAL -w -O0 -DDEBUG" )
# SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -DCEREAL -w -O0 -DDEBUG" )

#SET( CMAKE_CXX_FLAGS         "-std=c++11 -DCEREAL -DCXX11 -w -O3 -DDEBUG -ffast-math" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-std=c++11 -DCEREAL -DCXX11 -w -O3 -DDEBUG -ffast-math" )

#SET( CMAKE_CXX_FLAGS         "-w -O3 -std=c++11 -DCXX11 -DSCHEME_BENCHMARK -DDEBUG -ffast-math" )
#SET( CMAKE_EXE_LINKER_FLAGS  "-w -O3 -std=c++11 -DCXX11 -DSCHEME_BENCHMARK -DDEBUG -ffast-math" )

#SET( CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}        -std=c++14 -stdlib=libc++ -DCXX11 -DCXX14 -w -O3 -DDEBUG -ffast-math -funroll-loops -march=native" )
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lc++abi" )

# SET( CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}        -std=c++11 -DCXX11 -w -O3 -DDEBUG -ffast-math -funroll-loops -march=native" )
# SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lc++abi" )

# this removes those annoying ._* file turds you get with network fs on mac


# shitty "data" files to compile faster
list( APPEND LIB_CC ../schemelib/scheme/rosetta/score/EtableParams_init.lib.cc ../schemelib/scheme/rosetta/atype_map.lib.cc ../schemelib/scheme/chemical/ligand_factory.lib.cc )
SET_SOURCE_FILES_PROPERTIES( ${LIB_CC} PROPERTIES COMPILE_FLAGS -O0 )
add_library( scheme ${LIB_CC} )

FILE(GLOB L1 "../schemelib/scheme/*.cc")
FILE(GLOB L2 "../schemelib/scheme/*/*.cc")
FILE(GLOB L3 "../schemelib/scheme/*/*/*.cc")
FILE(GLOB L4 "../schemelib/scheme/*/*/*/*.cc")
FILE(GLOB L5 "../schemelib/scheme/*/*/*/*/*.cc")

list(APPEND EXTRA_LIBS boost_system boost_filesystem boost_iostreams z boost_timer )


include_directories(".")
link_directories(${CMAKE_LIBRARY_PATH})


project (main_test)
add_executable(main_test main_test.cc ${L1} ${L2} ${L3} ${L4} ${L5} )
target_link_libraries(main_test scheme ${EXTRA_LIBS} )


project (quick_test)
add_executable(quick_test quick_test.cc  )
target_link_libraries(quick_test scheme ${EXTRA_LIBS})
