
#include <nest/maps/SphereQuad.hh>
#include <nest/NEST.hh>
#include <nest/NEST_concepts.hh>
#include <nest/NEST_test_util.hh>
#include <gtest/gtest.h>

namespace scheme {
namespace nest {
namespace maps {

using std::cout;
using std::endl;

TEST( SphereQuad, cell_centers ){
	typedef SphereQuad<2,RowVector3d> MapType;
	typedef MapType::Params PRM;
	typedef MapType::ValueType VAL;
	MapType sd;
	PRM prm;
	VAL val;
	size_t celli;

	sd.params_to_value( PRM(0.5,0.5), 0, val ); ASSERT_EQ( val, VAL( 0, 0, 1) );
	sd.params_to_value( PRM(0.5,0.5), 1, val ); ASSERT_EQ( val, VAL( 0, 1, 0) );
	sd.params_to_value( PRM(0.5,0.5), 2, val ); ASSERT_EQ( val, VAL(-1, 0, 0) );
	sd.params_to_value( PRM(0.5,0.5), 3, val ); ASSERT_EQ( val, VAL( 0,-1, 0) );
	sd.params_to_value( PRM(0.5,0.5), 4, val ); ASSERT_EQ( val, VAL( 1, 0, 0) );
	sd.params_to_value( PRM(0.5,0.5), 5, val ); ASSERT_EQ( val, VAL( 0, 0,-1) );

	sd.value_to_params( VAL( 0, 0, 1), prm, celli ); ASSERT_TRUE( (prm==PRM(0.5,0.5)).minCoeff() ); ASSERT_EQ( celli, 0 );
	sd.value_to_params( VAL( 0, 1, 0), prm, celli ); ASSERT_TRUE( (prm==PRM(0.5,0.5)).minCoeff() ); ASSERT_EQ( celli, 1 );
	sd.value_to_params( VAL(-1, 0, 0), prm, celli ); ASSERT_TRUE( (prm==PRM(0.5,0.5)).minCoeff() ); ASSERT_EQ( celli, 2 );
	sd.value_to_params( VAL( 0,-1, 0), prm, celli ); ASSERT_TRUE( (prm==PRM(0.5,0.5)).minCoeff() ); ASSERT_EQ( celli, 3 );
	sd.value_to_params( VAL( 1, 0, 0), prm, celli ); ASSERT_TRUE( (prm==PRM(0.5,0.5)).minCoeff() ); ASSERT_EQ( celli, 4 );
	sd.value_to_params( VAL( 0, 0,-1), prm, celli ); ASSERT_TRUE( (prm==PRM(0.5,0.5)).minCoeff() ); ASSERT_EQ( celli, 5 );

	// sd.params_to_value( PRM(0.75,0.75), 5, val ); cout << val << endl;
	// sd.value_to_params( VAL( 1, 1, 1), prm, celli ); cout << prm.transpose() << " " << celli << endl;
}

TEST( SphereQuad, test_index_nesting_of_bincenters ){
	generic_test_index_nesting_of_bincenters( NEST<2,Vector3d,SphereQuad>() , 6 );
	generic_test_index_nesting_of_value( NEST<2,Vector3d,SphereQuad>() , Vector3d(1,0,0), 6 );
	generic_test_index_nesting_of_value( NEST<2,Vector3d,SphereQuad>() , Vector3d(1,1,1), 6 );
}

TEST( SphereQuad, test_index_nesting_of_bincenters_various_valuetypes ){
	generic_test_index_nesting_of_bincenters( NEST<2,RowVector3d,SphereQuad>(), 3 );
	generic_test_index_nesting_of_bincenters( NEST<2,RowVector3f,SphereQuad>(), 3 );	
	generic_test_index_nesting_of_bincenters( NEST<2,   Vector3f,SphereQuad>(), 3 );		
	generic_test_index_nesting_of_bincenters( NEST<2,Array<double,3,1>,SphereQuad>(), 3 );	
	generic_test_index_nesting_of_bincenters( NEST<2,Array<double,1,3>,SphereQuad>(), 3 );	
	generic_test_index_nesting_of_bincenters( NEST<2,Array<float,1,3>,SphereQuad>(), 3 );	
	generic_test_index_nesting_of_bincenters( NEST<2,Array<float,1,3>,SphereQuad>(), 3 );	
	generic_test_index_nesting_of_bincenters( NEST<2,concept::ArrayValueArchitype<3>,SphereQuad>(), 3 );	
}


}
}
}
