
#include <nest/NEST.hh>
#include <gtest/gtest.h>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <nest/parameter_maps.hh>

namespace scheme {
namespace nest {

using std::cout;
using std::endl;

using scheme::util::StorePointer;

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
void test_covering_radius(
	NEST nest,
	typename NEST::Params lb,
	typename NEST::Params ub
){
	boost::random::mt19937 rng(time(0));
	boost::uniform_real<> uniform;
	typename NEST::ValueType randpt;
	for(size_t r = 0; r <= 10; ++r){
		double maxdis = 0;
		for(size_t iter=0; iter < 10000/NEST::DIMENSION; ++iter){
			for(size_t i = 0; i < NEST::DIMENSION; ++i) 
				randpt[i] = uniform(rng)*(ub[i]-lb[i])+lb[i];
			nest.set_state( nest.get_index(randpt,r) ,r);
			double dist = (randpt-nest.value()).norm();
			// cout << randpt.transpose() << " " << nest.value().transpose() << endl;
			maxdis = fmax(maxdis,dist);
			ASSERT_LE( dist, nest.covering_radius(r) );
		}
		// covering radius should be reasonably tigth
		// cout << NEST::DIMENSION << " " << r << " " << nest.covering_radius(r) << " " << maxdis << endl;
		if( nest.covering_radius(r)*(1.0-(double)NEST::DIMENSION/20.0) > maxdis ){
			cout << "WARNING(PROBABILISTIC): covering radius may be too loose DIM=" << NEST::DIMENSION << " resl=" << r << endl;
			cout << "                        covering radius " << nest.covering_radius(r) << " max observerd: " << maxdis << endl;
			cout << "                        if you see a handful of these, don't worry. if you see lots, then worry" << endl;
		}
		// cout << DIM << " " << maxdis << " " << nest.covering_radius(r) << std::endl;
	}
}

TEST(NEST_unitmap,NEST_covering_radius_unitmap){
	{ NEST<1>::Params lb,ub; lb.fill(0); ub.fill(1); test_covering_radius( NEST<1>(), lb, ub ); }
	{ NEST<2>::Params lb,ub; lb.fill(0); ub.fill(1); test_covering_radius( NEST<2>(), lb, ub ); }
	{ NEST<3>::Params lb,ub; lb.fill(0); ub.fill(1); test_covering_radius( NEST<3>(), lb, ub ); }
	{ NEST<4>::Params lb,ub; lb.fill(0); ub.fill(1); test_covering_radius( NEST<4>(), lb, ub ); }
	{ NEST<5>::Params lb,ub; lb.fill(0); ub.fill(1); test_covering_radius( NEST<5>(), lb, ub ); }
	{ NEST<6>::Params lb,ub; lb.fill(0); ub.fill(1); test_covering_radius( NEST<6>(), lb, ub ); }
}

template<int DIM>
void test_covering_radius_scalemap(){
	typedef NEST<DIM,Eigen::Matrix<double,DIM,1>,ScaleMap> NestType;
	Eigen::Array<double ,DIM,1> lb, ub;
	for(size_t i = 0; i < DIM; ++i){
		lb[i] = -1;
		ub[i] = 1+2*i;
	}
	NestType nest(lb,ub);
	test_covering_radius< NestType >( nest, lb, ub );
}

TEST(NEST_scalemap,NEST_covering_radius_scalemap){
	test_covering_radius_scalemap<1>();
	test_covering_radius_scalemap<2>();
	test_covering_radius_scalemap<3>();
	test_covering_radius_scalemap<4>();
	test_covering_radius_scalemap<5>();
	test_covering_radius_scalemap<6>();
}

}
}

