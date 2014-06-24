
#include <nest/NEST.hh>
#include <gtest/gtest.h>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <nest/parameter_maps.hh>
#include <boost/foreach.hpp>

namespace scheme {
namespace nest {

using std::cout;
using std::endl;

using scheme::util::StorePointer;


void PrintTo(const RowVector2d & v, ::std::ostream* os) {
  *os << v[0] << " " << v[1] ;
}

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

		EXPECT_EQ( nest2.set_and_get(0,0)[0], -8.0 );
		EXPECT_EQ( nest2.set_and_get(1,0)[0],  8.0 );	
		EXPECT_EQ( nest2.set_and_get(0,1)[0],-12.0 );
		EXPECT_EQ( nest2.set_and_get(1,1)[0], -4.0 );	
		EXPECT_EQ( nest2.set_and_get(2,1)[0],  4.0 );
		EXPECT_EQ( nest2.set_and_get(3,1)[0], 12.0 );	
	}
	{
		typedef Vector2d VAL;
		Matrix<double,2,1> lb, ub;
		Matrix<size_t,2,1> bs;
		lb << -16.0,-24;
		ub <<  16.0, 24;
		bs << 2,3;
		NEST<2,VAL,ScaleMap> nest(lb,ub,bs);
		EXPECT_EQ( nest.set_and_get(0,0), Vector2d(-8,-16) );
		EXPECT_EQ( nest.set_and_get(1,0), Vector2d( 8,-16) );	
		EXPECT_EQ( nest.set_and_get(2,0), Vector2d(-8,  0) );
		EXPECT_EQ( nest.set_and_get(3,0), Vector2d( 8,  0) );	
		EXPECT_EQ( nest.set_and_get(4,0), Vector2d(-8, 16) );
		EXPECT_EQ( nest.set_and_get(5,0), Vector2d( 8, 16) );
		EXPECT_EQ( nest.set_and_get(0,1), Vector2d(-12,-20) );
		EXPECT_EQ( nest.set_and_get(1,1), Vector2d( -4,-20) );
		EXPECT_EQ( nest.set_and_get(2,1), Vector2d(-12,-12) );
		EXPECT_EQ( nest.set_and_get(3,1), Vector2d( -4,-12) );
		EXPECT_EQ( nest.set_and_get(0,2), Vector2d(-14,-22) );
		EXPECT_EQ( nest.set_and_get(1,2), Vector2d(-10,-22) );
		EXPECT_EQ( nest.set_and_get(2,2), Vector2d(-14,-18) );
		EXPECT_EQ( nest.set_and_get(3,2), Vector2d(-10,-18) );
		EXPECT_EQ( nest.set_and_get(0,3), Vector2d(-15,-23) );
		EXPECT_EQ( nest.set_and_get(1,3), Vector2d(-13,-23) );
		EXPECT_EQ( nest.set_and_get(2,3), Vector2d(-15,-21) );
		EXPECT_EQ( nest.set_and_get(3,3), Vector2d(-13,-21) );
		EXPECT_EQ( nest.set_and_get(5*64+60,3), Vector2d(13,21) );
		EXPECT_EQ( nest.set_and_get(5*64+61,3), Vector2d(15,21) );
		EXPECT_EQ( nest.set_and_get(5*64+62,3), Vector2d(13,23) );
		EXPECT_EQ( nest.set_and_get(5*64+63,3), Vector2d(15,23) );
		// cout << VAL(0,0) << endl;
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
	typename NEST::Params lb,
	typename NEST::Params ub
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
	{ NEST<1>::Params lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<1>(), lb, ub ); }
	{ NEST<2>::Params lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<2>(), lb, ub ); }
	{ NEST<3>::Params lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<3>(), lb, ub ); }
	{ NEST<4>::Params lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<4>(), lb, ub ); }
	{ NEST<5>::Params lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<5>(), lb, ub ); }
	{ NEST<6>::Params lb,ub; lb.fill(0); ub.fill(1); test_bin_circumradius( NEST<6>(), lb, ub ); }
}

template<int DIM>
void test_bin_circumradius_scalemap(){
	typedef NEST<DIM,Eigen::Matrix<double,DIM,1>,ScaleMap> NestType;
	Eigen::Array<double ,DIM,1> lb, ub;
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


TEST(ParamMap,UnitMap_get_neighboring_cells){
	typedef UnitMap<2,Eigen::Vector2d> MapType;
	MapType umap(99);
	std::vector<size_t> cnb;
	std::back_insert_iterator<std::vector<size_t> > biter(cnb);

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(1.5,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 1 );
	ASSERT_EQ( cnb[0], 1 );	

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(10.5,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 1 );
	ASSERT_EQ( cnb[0], 10 );	

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(10.3,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 2 );
	ASSERT_EQ( cnb[0],  9 );	
	ASSERT_EQ( cnb[1], 10 );

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(0.3,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 1 );
	ASSERT_EQ( cnb[0],  0 );	

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(0.3,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 1 );
	ASSERT_EQ( cnb[0],  0 );	

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(101.0,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 0 );

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(100.3,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 1 );
	ASSERT_EQ( cnb[0],  99 );

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(-0.5,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 0 );

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(-0.3,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 1 );
	ASSERT_EQ( cnb[0], 0 );	

}

TEST(ScaleMap,ScaleMap_get_neighboring_cells){
	typedef ScaleMap<2,Eigen::Vector2d> MapType;
	MapType umap(MapType::Indices(4,4));
	std::vector<size_t> cnb;
	std::back_insert_iterator<std::vector<size_t> > biter(cnb);

	cnb.clear();
	umap.get_neighboring_cells( MapType::ValueType(1.5,0.5), 0.4, biter );
	ASSERT_EQ( cnb.size(), 1 );
	ASSERT_EQ( cnb[0], 1 );	



	cout << cnb.size() << endl;
	BOOST_FOREACH(size_t i,cnb) cout << i << endl;
}


}
}

