
#include <nest/NEST.hh>
#include <gtest/gtest.h>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <nest/parameter_maps.hh>
#include <boost/foreach.hpp>
#include <iterator>

namespace scheme {
namespace nest {

using std::cout;
using std::endl;

using scheme::util::StorePointer;


TEST(NEST,particular_cases){
	NEST<1>	nest;
	ASSERT_EQ( nest.set_and_get(0,0)[0], 0.5 );
	ASSERT_EQ( nest.set_and_get(0,1)[0], 0.25 );
	ASSERT_EQ( nest.set_and_get(1,1)[0], 0.75 );
	ASSERT_EQ( nest.set_and_get(0,2)[0], 0.125 );
	ASSERT_EQ( nest.set_and_get(1,2)[0], 0.375 );
	ASSERT_EQ( nest.set_and_get(2,2)[0], 0.625 );
	ASSERT_EQ( nest.set_and_get(3,2)[0], 0.875 );
	NEST<2>	nest2;
	ASSERT_EQ( nest2.set_and_get(0,0), NEST<2>::ValueType(0.5,0.5)   );
	ASSERT_EQ( nest2.set_and_get(0,1), NEST<2>::ValueType(0.25,0.25)   );
	ASSERT_EQ( nest2.set_and_get(1,1), NEST<2>::ValueType(0.75,0.25)   );
	ASSERT_EQ( nest2.set_and_get(2,1), NEST<2>::ValueType(0.25,0.75)   );
	ASSERT_EQ( nest2.set_and_get(3,1), NEST<2>::ValueType(0.75,0.75)   );
}

TEST(NEST,map_discrete) {
	using std::cout;
	using std::endl;

	using namespace boost::assign;
	std::vector<double> choices;
	choices += 42.0,152.345,8049782.83402;
	NEST<0,double,DiscreteChoiceMap,StoreValue> nest0(choices);
	ASSERT_EQ(choices.size(),nest0.size());
	ASSERT_EQ(choices.size(),nest0.size(0));	
	ASSERT_EQ(choices.size(),nest0.size(3));
	ASSERT_EQ(choices.size(),nest0.size(4));			
	for(size_t i = 0; i < nest0.size(); ++i){
		ASSERT_TRUE( nest0.set_state(i) );
		ASSERT_EQ(choices[i],nest0.set_and_get(i));
		ASSERT_EQ(choices[i],nest0.set_and_get(i,0));
		ASSERT_EQ(choices[i],nest0.set_and_get(i,7));
		ASSERT_EQ(choices[i],nest0.set_and_get(i,3));
		ASSERT_EQ(choices[i],nest0.set_and_get(i,8));
	}
	ASSERT_FALSE( nest0.set_state(5) );
}


template<int DIM>
void test_coverage_unit1cell(){
	boost::random::mt19937 rng(time(0)); 
	boost::uniform_real<> uniform;
	NEST<DIM> nest;
	size_t rmax = 9/DIM+1;
	for(size_t r = 0; r <= rmax; ++r){
		double cellradius = sqrt(DIM) * 0.5 / (1<<r);
		for(size_t iter = 0; iter < 100000/DIM; ++iter){
			Matrix<double,DIM,1> tgt;
			for(size_t i = 0; i < DIM; ++i) tgt[i] = uniform(rng);
			size_t index = nest.get_index(tgt,r);
			Matrix<double,DIM,1> val = nest.set_and_get(index,r);
			ASSERT_LE( (tgt-val).norm(), cellradius );
		}
		Matrix<double,DIM,1> zeros; zeros.fill(0);
		ASSERT_LE( (zeros-nest.set_and_get(nest.get_index(zeros,r),r)).norm(), cellradius );
		Matrix<double,DIM,1> ones; ones.fill(0.999999);
		ASSERT_LE( (ones-nest.set_and_get(nest.get_index(ones,r),r)).norm(), cellradius );
	}
}

TEST(NEST,coverage_cube){
	test_coverage_unit1cell<1>();
	test_coverage_unit1cell<2>();
	test_coverage_unit1cell<3>();
	test_coverage_unit1cell<4>();
	test_coverage_unit1cell<5>();
	test_coverage_unit1cell<6>();
}



template<int DIM>
void test_index_nesting(){
	typedef Matrix<double,DIM,1> VAL;
	NEST<DIM> nest(2);
	size_t rmax = 9/DIM+1;
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
TEST(NEST,index_nesting){
	test_index_nesting<1>();
	test_index_nesting<2>();
	test_index_nesting<3>();
	test_index_nesting<4>();
	test_index_nesting<5>();
	test_index_nesting<6>();
}


template<int DIM>
void test_store_pointer_virtual(){
	typedef Matrix<double,DIM,1> VAL;
	NEST<DIM,VAL,UnitMap,StoreValue  > nest_val(2);
	NEST<DIM,VAL,UnitMap,StorePointer> nest_ptr(2);
	NestBase<> * nest_val_virtual = &nest_val;
	NestBase<> * nest_ptr_virtual = &nest_ptr;	
	StorePointer<VAL> * ptr_store      = dynamic_cast< StorePointer<VAL>* >(nest_ptr_virtual);
	StorePointer<VAL> * ptr_store_fail = dynamic_cast< StorePointer<VAL>* >(nest_val_virtual);
	ASSERT_TRUE(ptr_store);
	ASSERT_FALSE(ptr_store_fail);
	VAL val_pointed_to;
	ptr_store->set_pointer(&val_pointed_to);
	size_t rmax = 9/DIM+1;
	for(size_t r = 0; r <= rmax; ++r){
		for(size_t i = 0; i < nest_val.size(r); ++i){
			ASSERT_TRUE( nest_val_virtual->virtual_set_state(i,r) );
			ASSERT_TRUE( nest_ptr_virtual->virtual_set_state(i,r) );			
			ASSERT_EQ( nest_val.value(), val_pointed_to );
		}
	}
}
TEST(NEST,store_pointer_virtual){
	test_store_pointer_virtual<1>();
	test_store_pointer_virtual<2>();
	test_store_pointer_virtual<3>();
	test_store_pointer_virtual<4>();
	test_store_pointer_virtual<5>();
	test_store_pointer_virtual<6>();
}

template<int DIM>
void test_uniformity(){
	typedef Matrix<double,DIM,1> VAL;
	NEST<DIM,VAL,UnitMap,StoreValue  > nest(1);
	size_t rmax = 9/DIM+1;
	for(size_t r = 0; r <= rmax; ++r){
		double scale = 1<<r;
		ArrayXi counts(1<<r); counts.fill(0);
		for(size_t i = 0; i < nest.size(r); ++i){
			ASSERT_TRUE( nest.set_state(i,r) );
			for(size_t j = 0; j < DIM; ++j){
				counts[nest.value()[j]*scale]++;
			}
		}
		ASSERT_EQ( (1<<((DIM-1)*r))*DIM, counts.minCoeff() );		
		ASSERT_EQ( (1<<((DIM-1)*r))*DIM, counts.maxCoeff() );				
	}
}

TEST(NEST,uniformity){
	test_uniformity<1>();
	test_uniformity<2>();
	test_uniformity<3>();
	test_uniformity<4>();
	test_uniformity<5>();
	test_uniformity<6>();
}

TEST(NEST,bounds){
	NEST<1> nest;
	ASSERT_TRUE (nest.set_state(0,0));
	ASSERT_FALSE(nest.set_state(1,0));
	#ifndef NDEBUG
	ASSERT_DEATH( nest.set_and_get(1,0), "Assertion.*failed" );
	#endif
}



}
}

