#ifndef INCLUDED_scheme_nest_NEST_test_HH
#define INCLUDED_scheme_nest_NEST_test_HH

#include <nest/NEST.hh>
#include <gtest/gtest.h>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/unordered_set.hpp>

namespace scheme {
namespace nest {

using std::cout;
using std::endl;

TEST(NEST,particular_cases){
	NEST<1>	nest;
	EXPECT_EQ( nest.set_and_get(0,0)[0], 0.5 );
	EXPECT_EQ( nest.set_and_get(0,1)[0], 0.25 );
	EXPECT_EQ( nest.set_and_get(1,1)[0], 0.75 );
	EXPECT_EQ( nest.set_and_get(0,2)[0], 0.125 );
	EXPECT_EQ( nest.set_and_get(1,2)[0], 0.375 );
	EXPECT_EQ( nest.set_and_get(2,2)[0], 0.625 );
	EXPECT_EQ( nest.set_and_get(3,2)[0], 0.875 );
	NEST<2>	nest2;
	EXPECT_EQ( nest2.set_and_get(0,0), Vector2f(0.5,0.5)   );
	EXPECT_EQ( nest2.set_and_get(0,1), Vector2f(0.25,0.25)   );
	EXPECT_EQ( nest2.set_and_get(1,1), Vector2f(0.75,0.25)   );
	EXPECT_EQ( nest2.set_and_get(2,1), Vector2f(0.25,0.75)   );
	EXPECT_EQ( nest2.set_and_get(3,1), Vector2f(0.75,0.75)   );
}

TEST(NEST,map_discrete) {
	using std::cout;
	using std::endl;

	using namespace boost::assign;
	std::vector<double> choices;
	choices += 42.0,152.345,8049782.83402;
	NEST<0,double,DiscreteChoiceMap,StoreValue> nest0(choices);
	EXPECT_EQ(choices.size(),nest0.size());
	EXPECT_EQ(choices.size(),nest0.size(0));	
	EXPECT_EQ(choices.size(),nest0.size(3));
	EXPECT_EQ(choices.size(),nest0.size(4));			
	for(size_t i = 0; i < nest0.size(); ++i){
		EXPECT_TRUE( nest0.set_state(i) );
		EXPECT_EQ(choices[i],nest0.set_and_get(i));
		EXPECT_EQ(choices[i],nest0.set_and_get(i,0));
		EXPECT_EQ(choices[i],nest0.set_and_get(i,7));
		EXPECT_EQ(choices[i],nest0.set_and_get(i,3));
		EXPECT_EQ(choices[i],nest0.set_and_get(i,8));
	}
	EXPECT_FALSE( nest0.set_state(5) );
}


// template<int DIM>
// void test_get_index(){
// 	typedef Matrix<float,DIM,1> VAL;
// 	NEST<DIM> nest;
// 	size_t rmax = 9/DIM+1;
// 	for(size_t r = 0; r <= rmax; ++r){
// 		for(size_t i = 0; i < nest.size(r); ++i){
// 			VAL value = nest.set_and_get(i,r);
// 			size_t index = nest.get_index(value,r);
// 			EXPECT_EQ( i, index );
// 		}
// 	}
// }
// TEST(NEST,get_index){
// 	test_get_index<1>();
// 	test_get_index<2>();
// 	test_get_index<3>();
// 	test_get_index<4>();
// 	test_get_index<5>();
// 	test_get_index<6>();
// }


template<int DIM>
void test_index_nesting(){
	typedef Matrix<float,DIM,1> VAL;
	NEST<DIM> nest(2);
	size_t rmax = 9/DIM+1;
	for(size_t r = 0; r <= rmax; ++r){
		for(size_t i = 0; i < nest.size(r); ++i){
			EXPECT_TRUE( nest.set_state(i,r) );
			VAL value = nest.value();
			size_t index = nest.get_index(value,r);
			EXPECT_EQ( i, index );
			for(size_t r2 = 0; r2 <= r; ++r2){
				size_t index2 = nest.get_index(value,r2);
				EXPECT_EQ( i>>(DIM*(r-r2)), index2 );
			}
		}
	}
}
TEST(NEST,index_nesting){
	test_index_nesting<1>();
	test_index_nesting<2>();
	test_index_nesting<3>();
	test_index_nesting<4>();
	test_index_nesting<5>();
	test_index_nesting<6>();
}

template<int DIM>
void test_index_lookup_scaled(){
	typedef Matrix<float,DIM,1> VAL;
	Matrix<float,10,1> lb0,ub0;
	Matrix<size_t,10,1> bs0;
	lb0 << -1.3,2.2,0,-3,-5,-9.9,1.3,44,-13.3,99;
	ub0 <<  1.3,4.2,1,10,-3, 9.9,4.3,44, 13.3,199;
	bs0 << 2,1,2,3,4,5,6,7,8,9;

	typedef Matrix<float,DIM,1> VAL;
	VAL lb, ub;
	Matrix<size_t,DIM,1> bs;
	for(size_t i = 0; i < DIM; ++i){
		lb[i] = lb0[i]; ub[i] = ub0[i]; bs[i] = bs0[i];
	}
	NEST<DIM,VAL,ScaleMap,StoreValue> nest(lb,ub,bs);
	size_t rmax = 6/DIM;
	for(size_t r = 0; r <= rmax; ++r){
		for(size_t i = 0; i < nest.size(r); ++i){
			EXPECT_TRUE( nest.set_state(i,r) );
			VAL value = nest.value();
			size_t index = nest.get_index(value,r);
			EXPECT_EQ( i, index );
			for(size_t r2 = 0; r2 <= r; ++r2){
				size_t index2 = nest.get_index(value,r2);
				EXPECT_EQ( i>>(DIM*(r-r2)), index2 );
			}
		}
	}
}
TEST(NEST,index_lookup_scaled){
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
	Matrix<float,10,1> lb0,ub0;
	Matrix<size_t,10,1> bs0;
	lb0 << -1.3,2.2,0,-3,-5,-9.9,1.3,44,-13.3,99;
	ub0 <<  1.3,4.2,1,10,-3, 9.9,4.3,44, 13.3,199;
	bs0 << 1,2,3,4,5,6,7,8,9,10;

	typedef Matrix<float,DIM,1> VAL;
	VAL lb, ub;
	Matrix<size_t,DIM,1> bs;
	for(size_t i = 0; i < DIM; ++i){
		lb[i] = lb0[i]; ub[i] = ub0[i]; bs[i] = bs0[i];
	}
	NEST<DIM,VAL,ScaleMap,StoreValue> nest(lb,ub,bs);

	size_t resl = 6/DIM+1;
	for(size_t i = 0; i < nest.size(resl); ++i){
		EXPECT_TRUE( nest.set_state(i,resl) );
		for(size_t j = 0; j < DIM; ++j){ 
			EXPECT_LT( lb[j], nest.value()[j] ); EXPECT_LT( nest.value()[j] , ub[j] );
		}
		EXPECT_FALSE( nest.set_state(i+nest.size(resl),resl) );
		for(size_t j = 0; j < DIM; ++j){ 
			EXPECT_LT( lb[j], nest.value()[j] ); EXPECT_LT( nest.value()[j] , ub[j] );
		}
	}

}
TEST(NEST,map_scale){
	test_map_scale_bounds<1>();
	test_map_scale_bounds<2>();
	test_map_scale_bounds<3>();
	test_map_scale_bounds<4>();
	// test_map_scale_bounds<5>();
	// test_map_scale_bounds<6>();
}



}
}

#endif
