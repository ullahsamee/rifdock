
# attemp to use a scons build...
#set( ROSETTA_SRC /work/sheffler/ssd/rosetta/master/source )
#set( ROSETTA_BIN release/linux/3.10/64/x86/gcc/6.1/cxx11-omp )
#link_directories( ${ROSETTA_SRC}/build/src/${ROSETTA_BIN} )
#link_directories( ${ROSETTA_SRC}/build/external/${ROSETTA_BIN} )


set( ROSETTA_SRC $ENV{CMAKE_ROSETTA_PATH}/source )

link_directories( ${ROSETTA_BIN} )

set( ALL_ROSETTA_LIBS 
	basic core.1 core.2 core.3 core.4 core.5 devel numeric ObjexxFCL protocols.1 protocols.3 protocols.6 protocols.7 
	protocols_a.2 protocols_a.4 protocols_a.5 protocols_b.2 protocols_b.4 protocols_b.5 protocols_c.4 protocols_c.5 
	protocols_d.4 protocols_d.5 protocols_e.4 protocols_e.5 protocols_f.4 protocols_g.4 protocols_h.4 utility 
	cifparse cppdb sqlite3 z
)
include_directories( ${ROSETTA_SRC}/src )
include_directories( ${ROSETTA_SRC}/external/include )
include_directories( ${ROSETTA_SRC}/external/dbio )
include_directories( ${ROSETTA_SRC}/external/boost_1_55_0 )
include_directories( ${ROSETTA_SRC}/src/platform/linux )
include_directories( . )


SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCXX11 -DPTR_STD -DMULTI_THREADED" )

add_subdirectory( python )

# this omp shit breaks python bindings! do it after
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_OPENMP -fopenmp" )
LIST( APPEND ALL_ROSETTA_LIBS gomp )

add_subdirectory( riflib)

add_executable( test_librosetta test_librosetta.cc  )
target_link_libraries ( test_librosetta ${ALL_ROSETTA_LIBS} )

add_executable( rifgen rifgen.cc  )
target_link_libraries ( rifgen riflib ${ALL_ROSETTA_LIBS} )
install ( TARGETS rifgen RUNTIME DESTINATION bin )

add_executable( rif_dock_test rif_dock_test.cc  )
target_link_libraries ( rif_dock_test riflib ${ALL_ROSETTA_LIBS} )
install ( TARGETS rif_dock_test RUNTIME DESTINATION bin )

