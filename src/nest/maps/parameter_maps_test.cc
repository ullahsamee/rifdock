#include <nest/maps/parameter_maps.hh>
#include <nest/NEST.hh>
#include <nest/NEST_test_util.hh>
#include <gtest/gtest.h>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/foreach.hpp>

namespace scheme {
namespace nest {
namespace maps {


using std::cout;
using std::endl;

using scheme::nest::StorePointer;


TEST(NEST_scalemap,particular_values){
	{
		typedef Matrix<double,1,1> VAL;
		VAL lb, ub;
		Matrix<size_t,1,1> bs;
		lb << -16.0;
		ub <<  16.0;
		bs << 1;
		NEST<1,VAL,ScaleMap> nest1(lb,ub,bs);
		bs << 2;
		NEST<1,VAL,ScaleMap> nest2(lb,ub,bs);
		ASSERT_EQ( nest1.set_and_get(0,0)[0],  0.0 );
		ASSERT_EQ( nest1.set_and_get(0,1)[0], -8.0 );
		ASSERT_EQ( nest1.set_and_get(1,1)[0],  8.0 );

		ASSERT_EQ( nest2.set_and_get(0,0)[0], -8.0 );
		ASSERT_EQ( nest2.set_and_get(1,0)[0],  8.0 );	
		ASSERT_EQ( nest2.set_and_get(0,1)[0],-12.0 );
		ASSERT_EQ( nest2.set_and_get(1,1)[0], -4.0 );	
		ASSERT_EQ( nest2.set_and_get(2,1)[0],  4.0 );
		ASSERT_EQ( nest2.set_and_get(3,1)[0], 12.0 );	
	}
	{
		typedef Vector2d VAL;
		Matrix<double,2,1> lb, ub;
		Matrix<size_t,2,1> bs;
		lb << -16.0,-24;
		ub <<  16.0, 24;
		bs << 2,3;
		NEST<2,VAL,ScaleMap> nest(lb,ub,bs);
		ASSERT_EQ( nest.set_and_get(0,0), Vector2d(-8,-16) );
		ASSERT_EQ( nest.set_and_get(1,0), Vector2d( 8,-16) );	
		ASSERT_EQ( nest.set_and_get(2,0), Vector2d(-8,  0) );
		ASSERT_EQ( nest.set_and_get(3,0), Vector2d( 8,  0) );	
		ASSERT_EQ( nest.set_and_get(4,0), Vector2d(-8, 16) );
		ASSERT_EQ( nest.set_and_get(5,0), Vector2d( 8, 16) );
		ASSERT_EQ( nest.set_and_get(0,1), Vector2d(-12,-20) );
		ASSERT_EQ( nest.set_and_get(1,1), Vector2d( -4,-20) );
		ASSERT_EQ( nest.set_and_get(2,1), Vector2d(-12,-12) );
		ASSERT_EQ( nest.set_and_get(3,1), Vector2d( -4,-12) );
		ASSERT_EQ( nest.set_and_get(0,2), Vector2d(-14,-22) );
		ASSERT_EQ( nest.set_and_get(1,2), Vector2d(-10,-22) );
		ASSERT_EQ( nest.set_and_get(2,2), Vector2d(-14,-18) );
		ASSERT_EQ( nest.set_and_get(3,2), Vector2d(-10,-18) );
		ASSERT_EQ( nest.set_and_get(0,3), Vector2d(-15,-23) );
		ASSERT_EQ( nest.set_and_get(1,3), Vector2d(-13,-23) );
		ASSERT_EQ( nest.set_and_get(2,3), Vector2d(-15,-21) );
		ASSERT_EQ( nest.set_and_get(3,3), Vector2d(-13,-21) );
		ASSERT_EQ( nest.set_and_get(5*64+60,3), Vector2d(13,21) );
		ASSERT_EQ( nest.set_and_get(5*64+61,3), Vector2d(15,21) );
		ASSERT_EQ( nest.set_and_get(5*64+62,3), Vector2d(13,23) );
		ASSERT_EQ( nest.set_and_get(5*64+63,3), Vector2d(15,23) );
		// cout << VAL(0,0) << endl;
	}
}


TEST(UnitMap,value_to_params_for_cell){
	typedef UnitMap<2,RowVector2d> MapType;
	MapType umap(3);
	MapType::Params params;

	umap.value_to_params_for_cell( MapType::ValueType(1.5,0.5), params, 0  );
	ASSERT_EQ( params[0], 1.5 );
	ASSERT_EQ( params[1], 0.5 );

	umap.value_to_params_for_cell( MapType::ValueType(1.5,0.5), params, 1  );
	ASSERT_EQ( params[0], 0.5 );
	ASSERT_EQ( params[1], 0.5 );

	umap.value_to_params_for_cell( MapType::ValueType(1.5,0.5), params, 2  );
	ASSERT_EQ( params[0], -0.5 );
	ASSERT_EQ( params[1],  0.5 );

	umap.value_to_params_for_cell( MapType::ValueType(0.2,0.5), params, 2  );
	ASSERT_EQ( params[0], -1.8 );
	ASSERT_EQ( params[1],  0.5 );
}

TEST(ScaleMap,value_to_params_for_cell){
	{
		typedef ScaleMap<2,RowVector2d> MapType;
		MapType smap(
			MapType::Params(0,0),
			MapType::Params(4,4),
			MapType::Indices(4,4)
		);

		MapType::Params params;

		smap.value_to_params_for_cell( MapType::ValueType(1.5,0.5), params, smap.indices_to_cellindex(MapType::Indices(0,0)) );
		ASSERT_EQ( params[0], 1.5 );
		ASSERT_EQ( params[1], 0.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5,0.5), params, smap.indices_to_cellindex(MapType::Indices(1,0)) );
		ASSERT_EQ( params[0], 0.5 );
		ASSERT_EQ( params[1], 0.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5,0.5), params, smap.indices_to_cellindex(MapType::Indices(0,1)) );
		ASSERT_EQ( params[0],  1.5 );
		ASSERT_EQ( params[1], -0.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5,0.5), params, smap.indices_to_cellindex(MapType::Indices(2,0)) );
		ASSERT_EQ( params[0],-0.5 );
		ASSERT_EQ( params[1], 0.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5,0.5), params, smap.indices_to_cellindex(MapType::Indices(1,2)) );
		ASSERT_EQ( params[0],  0.5 );
		ASSERT_EQ( params[1], -1.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5,0.5), params, smap.indices_to_cellindex(MapType::Indices(3,3)) );
		ASSERT_EQ( params[0], -1.5 );
		ASSERT_EQ( params[1], -2.5 );
	}
	{
		typedef ScaleMap<2,RowVector2d> MapType;
		double scale = 2.0;
		MapType::Params shift(0.597,1.1243);
		MapType smap(
			MapType::Params(0.0,0.0)*scale-shift,
			MapType::Params(4.0,4.0)*scale-shift,
			MapType::Indices(4,4)
		);

		MapType::Params params;

		smap.value_to_params_for_cell( MapType::ValueType(1.5*scale-shift[0],0.5*scale-shift[1]), params, smap.indices_to_cellindex(MapType::Indices(0,0)) );
		ASSERT_EQ( params[0], 1.5 );
		ASSERT_EQ( params[1], 0.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5*scale-shift[0],0.5*scale-shift[1]), params, smap.indices_to_cellindex(MapType::Indices(1,0)) );
		ASSERT_EQ( params[0], 0.5 );
		ASSERT_EQ( params[1], 0.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5*scale-shift[0],0.5*scale-shift[1]), params, smap.indices_to_cellindex(MapType::Indices(0,1)) );
		ASSERT_EQ( params[0],  1.5 );
		ASSERT_EQ( params[1], -0.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5*scale-shift[0],0.5*scale-shift[1]), params, smap.indices_to_cellindex(MapType::Indices(2,0)) );
		ASSERT_EQ( params[0],-0.5 );
		ASSERT_EQ( params[1], 0.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5*scale-shift[0],0.5*scale-shift[1]), params, smap.indices_to_cellindex(MapType::Indices(1,2)) );
		ASSERT_EQ( params[0],  0.5 );
		ASSERT_EQ( params[1], -1.5 );

		smap.value_to_params_for_cell( MapType::ValueType(1.5*scale-shift[0],0.5*scale-shift[1]), params, smap.indices_to_cellindex(MapType::Indices(3,3)) );
		ASSERT_EQ( params[0], -1.5 );
		ASSERT_EQ( params[1], -2.5 );
	}
}





template<int DIM>
void test_index_lookup_scaled(){
	typedef Matrix<double,DIM,1> VAL;
	Matrix<double,10,1> lb0,ub0;
	Matrix<size_t,10,1> bs0;
	lb0 << -1.3,2.2,0,-3,-5,-9.9,1.3,44,-13.3,99;
	ub0 <<  1.3,4.2,1,10,-3, 9.9,4.3,44, 13.3,199;
	bs0 << 2,1,2,3,4,5,6,7,8,9;

	typedef Matrix<double,DIM,1> VAL;
	VAL lb, ub;
	Matrix<size_t,DIM,1> bs;
	for(size_t i = 0; i < DIM; ++i){
		lb[i] = lb0[i]; ub[i] = ub0[i]; bs[i] = bs0[i];
	}
	NEST<DIM,VAL,ScaleMap,StoreValue> nest(lb,ub,bs);
	size_t rmax = 6/DIM;
	for(size_t r = 0; r <= rmax; ++r){
		for(size_t i = 0; i < nest.size(r); ++i){
			ASSERT_TRUE( nest.set_state(i,r) );
			VAL value = nest.value();
			size_t index = nest.get_index(value,r);
			ASSERT_EQ( i, index );
			for(size_t r2 = 0; r2 <= r; ++r2){
				size_t index2 = nest.get_index(value,r2);
				ASSERT_EQ( i>>(DIM*(r-r2)), index2 );
			}
		}
	}
}

TEST(NEST_scalemap,index_lookup_scaled){
	test_index_lookup_scaled<1>();
	test_index_lookup_scaled<2>();
	test_index_lookup_scaled<3>();
	test_index_lookup_scaled<4>();
	test_index_lookup_scaled<5>();
	test_index_lookup_scaled<6>();
}


template<int DIM>
void test_map_scale_bounds(){
	BOOST_STATIC_ASSERT(DIM<10);
	Matrix<double,10,1> lb0,ub0;
	Matrix<size_t,10,1> bs0;
	lb0 << -1.3,2.2,0,-3,-5,-9.9,1.3,44,-13.3,99;
	ub0 <<  1.3,4.2,1,10,-3, 9.9,4.3,44, 13.3,199;
	bs0 << 1,2,3,4,5,6,7,8,9,10;

	typedef Matrix<double,DIM,1> VAL;
	VAL lb, ub;
	Matrix<size_t,DIM,1> bs;
	for(size_t i = 0; i < DIM; ++i){
		lb[i] = lb0[i]; ub[i] = ub0[i]; bs[i] = bs0[i];
	}
	NEST<DIM,VAL,ScaleMap,StoreValue> nest(lb,ub,bs);

	size_t resl = 8/DIM;
	for(size_t i = 0; i < nest.size(resl); ++i){
		ASSERT_TRUE( nest.set_state(i,resl) );
		for(size_t j = 0; j < DIM; ++j){ 
			ASSERT_LT( lb[j], nest.value()[j] ); ASSERT_LT( nest.value()[j] , ub[j] );
		}
		ASSERT_FALSE( nest.set_state(i+nest.size(resl),resl) );
		for(size_t j = 0; j < DIM; ++j){ 
			ASSERT_LT( lb[j], nest.value()[j] ); ASSERT_LT( nest.value()[j] , ub[j] );
		}
	}

}

TEST(NEST_scalemap,map_scale){ 
	test_map_scale_bounds<1>();
	test_map_scale_bounds<2>();
	test_map_scale_bounds<3>();
	test_map_scale_bounds<4>();
	test_map_scale_bounds<5>();
	test_map_scale_bounds<6>();
}

template<class NEST>
void test_bin_circumradius(
	NEST nest,
	typename NEST::ValueType lb,
	typename NEST::ValueType ub
){
	boost::random::mt19937 rng(time(0));
	boost::uniform_real<> uniform;
	typename NEST::ValueType randpt;
	for(size_t r = 0; r <= 10; ++r){
		double maxdis = 0;
		for(size_t iter=0; iter < 100000/NEST::DIMENSION; ++iter){
			for(size_t i = 0; i < NEST::DIMENSION; ++i) 
				randpt[i] = uniform(rng)*(ub[i]-lb[i])+lb[i];
			nest.set_state( nest.get_index(randpt,r) ,r);
			double dist = (randpt-nest.value()).norm();
			// cout << randpt.transpose() << " " << nest.value().transpose() << endl;
			maxdis = fmax(maxdis,dist);
			ASSERT_LE( dist, nest.bin_circumradius(r) );
		}
		// covering radius should be reasonably tigth
		// cout << NEST::DIMENSION << " " << r << " " << nest.bin_circumradius(r) << " " << maxdis << endl;
		if( nest.bin_circumradius(r)*(1.0-(double)NEST::DIMENSION/20.0) > maxdis ){
			cout << "WARNING(PROBABILISTIC): covering radius may be too loose DIM=" << NEST::DIMENSION << " resl=" << r << endl;
			cout << "                        covering radius " << nest.bin_circumradius(r) << " max observerd: " << maxdis << endl;
			cout << "                        if you see a handful of these, don't worry. if you see lots, then worry" << endl;
		}
		// cout << DIM << " " << maxdis << " " << nest.bin_circumradius(r) << std::endl;
	}
}

TEST(NEST_unitmap,NEST_bin_circumradius_unitmap){
	{ NEST<1>::ValueType lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<1>(), lb, ub ); }
	{ NEST<2>::ValueType lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<2>(), lb, ub ); }
	{ NEST<3>::ValueType lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<3>(), lb, ub ); }
	{ NEST<4>::ValueType lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<4>(), lb, ub ); }
	{ NEST<5>::ValueType lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<5>(), lb, ub ); }
	{ NEST<6>::ValueType lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<6>(), lb, ub ); }
}

template<int DIM>
void test_bin_circumradius_scalemap(){
	typedef NEST<DIM,Eigen::Matrix<double,DIM,1>,ScaleMap> NestType;
	typename NestType::ValueType lb, ub;
	for(size_t i = 0; i < DIM; ++i){
		lb[i] = -1;
		ub[i] = 1+2*i;
	}
	NestType nest(lb,ub);
	test_bin_circumradius< NestType >( nest, lb, ub );
}

TEST(NEST_scalemap,NEST_bin_circumradius_scalemap){
	test_bin_circumradius_scalemap<1>();
	test_bin_circumradius_scalemap<2>();
	test_bin_circumradius_scalemap<3>();
	test_bin_circumradius_scalemap<4>();
	test_bin_circumradius_scalemap<5>();
	test_bin_circumradius_scalemap<6>();
}

TEST(NEST_scalemap,test_coverage){
	boost::random::mt19937 rng(time(0));
	boost::uniform_real<> uniform;

	typedef NEST<2,RowVector2d,ScaleMap> NestType;
	NestType nest;
	std::vector<double> largest_d2_for_r(20,0.0);
	for(size_t i = 0; i < 10000; ++i){
		NestType::ValueType val;
		for(size_t j = 0; j < NestType::DIMENSION; ++j) val[j] = uniform(rng);
		generic_test_coverage_of_value( nest, val, largest_d2_for_r, 10 );
	}
	for(size_t r = 0; r <= 10; ++r) ASSERT_LT(  nest.bin_circumradius(r), 1.1*sqrt(largest_d2_for_r[r]) );
}


}
}
}

