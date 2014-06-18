
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

TEST(NEST,neighbor_lookup_single_cell){
	{
		NEST<2> nest;
		std::vector<size_t> neighbors;
	 	NEST<2>::ValueType val;
	 	val << 0.85,0.15;
	 	size_t r = 1;
	 	std::back_insert_iterator< std::vector<size_t> > back_it(neighbors);
	 	nest.get_neighbors( val, r, back_it );
	 	ASSERT_EQ(neighbors.size(),4);
	 	ASSERT_EQ(neighbors[0],0);
	 	ASSERT_EQ(neighbors[1],1);
	 	ASSERT_EQ(neighbors[2],2);
	 	ASSERT_EQ(neighbors[3],3);
	}
	{
		NEST<3> nest;
		std::vector<size_t> neighbors;
	 	NEST<3>::ValueType val;
	 	val << 0.85,0.15,0.5;
	 	size_t r = 3;
	 	std::back_insert_iterator< std::vector<size_t> > back_it(neighbors);
	 	nest.get_neighbors( val, r, back_it );
	 	ASSERT_EQ(neighbors.size(),27);
		ASSERT_EQ(neighbors[0],101);
		ASSERT_EQ(neighbors[1],108);
		ASSERT_EQ(neighbors[2],109);
		ASSERT_EQ(neighbors[3],103);
		ASSERT_EQ(neighbors[4],110);
		ASSERT_EQ(neighbors[5],111);
		ASSERT_EQ(neighbors[6],117);
		ASSERT_EQ(neighbors[7],124);
		ASSERT_EQ(neighbors[8],125);
		ASSERT_EQ(neighbors[9],321);
		ASSERT_EQ(neighbors[10],328);
		ASSERT_EQ(neighbors[11],329);
		ASSERT_EQ(neighbors[12],323);
		ASSERT_EQ(neighbors[13],330);
		ASSERT_EQ(neighbors[14],331);
		ASSERT_EQ(neighbors[15],337);
		ASSERT_EQ(neighbors[16],344);
		ASSERT_EQ(neighbors[17],345);
		ASSERT_EQ(neighbors[18],325);
		ASSERT_EQ(neighbors[19],332);
		ASSERT_EQ(neighbors[20],333);
		ASSERT_EQ(neighbors[21],327);
		ASSERT_EQ(neighbors[22],334);
		ASSERT_EQ(neighbors[23],335);
		ASSERT_EQ(neighbors[24],341);
		ASSERT_EQ(neighbors[25],348);
		ASSERT_EQ(neighbors[26],349);
	}
}

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
	ASSERT_EQ( nest2.set_and_get(0,0), Vector2f(0.5,0.5)   );
	ASSERT_EQ( nest2.set_and_get(0,1), Vector2f(0.25,0.25)   );
	ASSERT_EQ( nest2.set_and_get(1,1), Vector2f(0.75,0.25)   );
	ASSERT_EQ( nest2.set_and_get(2,1), Vector2f(0.25,0.75)   );
	ASSERT_EQ( nest2.set_and_get(3,1), Vector2f(0.75,0.75)   );
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
void test_coverage_cube(){
	boost::random::mt19937 rng; 
	boost::uniform_real<> uniform;
	NEST<DIM> nest;
	size_t rmax = 9/DIM+1;
	for(size_t r = 0; r <= rmax; ++r){
		float cellradius = sqrt(DIM) * 0.5 / (1<<r);
		for(size_t iter = 0; iter < 10000/DIM; ++iter){
			Matrix<float,DIM,1> tgt;
			for(size_t i = 0; i < DIM; ++i) tgt[i] = uniform(rng);
			size_t index = nest.get_index(tgt,r);
			Matrix<float,DIM,1> val = nest.set_and_get(index,r);
			ASSERT_LE( (tgt-val).norm(), cellradius );
		}
		Matrix<float,DIM,1> zeros; zeros.fill(0);
		ASSERT_LE( (zeros-nest.set_and_get(nest.get_index(zeros,r),r)).norm(), cellradius );
		Matrix<float,DIM,1> ones; ones.fill(0.999999);
		ASSERT_LE( (ones-nest.set_and_get(nest.get_index(ones,r),r)).norm(), cellradius );
	}
}

TEST(NEST,coverage_cube_DIM1){ test_coverage_cube<1>();    // }
/*TEST(NEST,coverage_cube_DIM2){*/ test_coverage_cube<2>();// }
/*TEST(NEST,coverage_cube_DIM3){*/ test_coverage_cube<3>();// }
/*TEST(NEST,coverage_cube_DIM4){*/ test_coverage_cube<4>();// }
/*TEST(NEST,coverage_cube_DIM5){*/ test_coverage_cube<5>();// }
/*TEST(NEST,coverage_cube_DIM6){*/ test_coverage_cube<6>(); }



template<int DIM>
void test_index_nesting(){
	typedef Matrix<float,DIM,1> VAL;
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
TEST(NEST,index_nesting_DIM1){ test_index_nesting<1>();     //}
/*TEST(NEST,index_nesting_DIM2){*/ test_index_nesting<2>(); //}
/*TEST(NEST,index_nesting_DIM3){*/ test_index_nesting<3>(); //}
/*TEST(NEST,index_nesting_DIM4){*/ test_index_nesting<4>(); //}
/*TEST(NEST,index_nesting_DIM5){*/ test_index_nesting<5>(); //}
/*TEST(NEST,index_nesting_DIM6){*/ test_index_nesting<6>(); }


template<int DIM>
void test_store_pointer_generic(){
	typedef Matrix<float,DIM,1> VAL;
	NEST<DIM,VAL,UnitMap,StoreValue  > nest_val(2);
	NEST<DIM,VAL,UnitMap,StorePointer> nest_ptr(2);
	NestBase<> * nest_val_generic = &nest_val;
	NestBase<> * nest_ptr_generic = &nest_ptr;	
	StorePointer<VAL> * ptr_store      = dynamic_cast< StorePointer<VAL>* >(nest_ptr_generic);
	StorePointer<VAL> * ptr_store_fail = dynamic_cast< StorePointer<VAL>* >(nest_val_generic);
	ASSERT_TRUE(ptr_store);
	ASSERT_FALSE(ptr_store_fail);
	VAL val_pointed_to;
	ptr_store->set_pointer(&val_pointed_to);
	size_t rmax = 9/DIM+1;
	for(size_t r = 0; r <= rmax; ++r){
		for(size_t i = 0; i < nest_val.size(r); ++i){
			ASSERT_TRUE( nest_val_generic->generic_set_state(i,r) );
			ASSERT_TRUE( nest_ptr_generic->generic_set_state(i,r) );			
			ASSERT_EQ( nest_val.value(), val_pointed_to );
		}
	}
}
TEST(NEST,store_pointer_generic_DIM1){ test_store_pointer_generic<1>();     //}
/*TEST(NEST,store_pointer_generic_DIM2){*/ test_store_pointer_generic<2>(); //}
/*TEST(NEST,store_pointer_generic_DIM3){*/ test_store_pointer_generic<3>(); //}
/*TEST(NEST,store_pointer_generic_DIM4){*/ test_store_pointer_generic<4>(); //}
/*TEST(NEST,store_pointer_generic_DIM5){*/ test_store_pointer_generic<5>(); //}
/*TEST(NEST,store_pointer_generic_DIM6){*/ test_store_pointer_generic<6>(); }

template<int DIM>
void test_uniformity(){
	typedef Matrix<float,DIM,1> VAL;
	NEST<DIM,VAL,UnitMap,StoreValue  > nest(1);
	size_t rmax = 9/DIM+1;
	for(size_t r = 0; r <= rmax; ++r){
		float scale = 1<<r;
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

TEST(NEST,uniformity_DIM1){ test_uniformity<1>();    //}
/*TEST(NEST,uniformity_DIM2){*/ test_uniformity<2>(); //}
/*TEST(NEST,uniformity_DIM3){*/ test_uniformity<3>(); //}
/*TEST(NEST,uniformity_DIM4){*/ test_uniformity<4>(); //}
/*TEST(NEST,uniformity_DIM5){*/ test_uniformity<5>(); //}
/*TEST(NEST,uniformity_DIM6){*/ test_uniformity<6>(); }




}
}

