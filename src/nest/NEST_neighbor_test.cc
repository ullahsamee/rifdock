#include <nest/NEST.hh>
#include <gtest/gtest.h>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <nest/parameter_maps.hh>
#include <boost/foreach.hpp>
#include <iterator>
#include <array>


namespace scheme {
namespace nest {

using std::cout;
using std::endl;

using scheme::util::StorePointer;

TEST(NEST_NEIGHBOR,dim2_test_case){
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
TEST(NEST_NEIGHBOR,dim3_test_case){
	NEST<3,RowVector3d> nest;
 	NEST<3,RowVector3d>::ValueType val;
 	val << 0.85,0.15,0.5;
 	// val[0] = 0.85; val[1] = 0.15; val[2] = 0.5;
 	size_t r = 3;
	std::vector<size_t> neighbors;
 	std::back_insert_iterator< std::vector<size_t> > back_it(neighbors);
 	nest.get_neighbors( val, r, back_it );
 	ASSERT_EQ(neighbors.size(),27);
	EXPECT_EQ( nest.set_and_get(neighbors[ 0],r),  RowVector3d( 0.6875, 0.0625, 0.4375 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[ 1],r),  RowVector3d( 0.8125, 0.0625, 0.4375 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[ 2],r),  RowVector3d( 0.9375, 0.0625, 0.4375 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[ 3],r),  RowVector3d( 0.6875, 0.1875, 0.4375 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[ 4],r),  RowVector3d( 0.8125, 0.1875, 0.4375 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[ 5],r),  RowVector3d( 0.9375, 0.1875, 0.4375 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[ 6],r),  RowVector3d( 0.6875, 0.3125, 0.4375 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[ 7],r),  RowVector3d( 0.8125, 0.3125, 0.4375 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[ 8],r),  RowVector3d( 0.9375, 0.3125, 0.4375 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[ 9],r),  RowVector3d( 0.6875, 0.0625, 0.5625 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[10],r),  RowVector3d( 0.8125, 0.0625, 0.5625 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[11],r),  RowVector3d( 0.9375, 0.0625, 0.5625 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[12],r),  RowVector3d( 0.6875, 0.1875, 0.5625 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[13],r),  RowVector3d( 0.8125, 0.1875, 0.5625 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[14],r),  RowVector3d( 0.9375, 0.1875, 0.5625 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[15],r),  RowVector3d( 0.6875, 0.3125, 0.5625 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[16],r),  RowVector3d( 0.8125, 0.3125, 0.5625 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[17],r),  RowVector3d( 0.9375, 0.3125, 0.5625 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[18],r),  RowVector3d( 0.6875, 0.0625, 0.6875 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[19],r),  RowVector3d( 0.8125, 0.0625, 0.6875 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[20],r),  RowVector3d( 0.9375, 0.0625, 0.6875 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[21],r),  RowVector3d( 0.6875, 0.1875, 0.6875 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[22],r),  RowVector3d( 0.8125, 0.1875, 0.6875 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[23],r),  RowVector3d( 0.9375, 0.1875, 0.6875 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[24],r),  RowVector3d( 0.6875, 0.3125, 0.6875 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[25],r),  RowVector3d( 0.8125, 0.3125, 0.6875 ) );
	EXPECT_EQ( nest.set_and_get(neighbors[26],r),  RowVector3d( 0.9375, 0.3125, 0.6875 ) );
}


TEST(NEST_NEIGHBOR,unit_1d_boundary){
	NEST<1> nest;
	std::vector<size_t> neighbors;
 	std::back_insert_iterator< std::vector<size_t> > back_it(neighbors);
 	size_t r = 2;
 	NEST<1>::ValueType val;

 	neighbors.clear();
 	val<<1.1;
 	nest.get_neighbors_unitcell( val, r, back_it );
 	// BOOST_FOREACH(size_t i,neighbors){ cout << i << " " ; } cout << endl;
 	ASSERT_EQ( neighbors.size(), 1 );
 	ASSERT_EQ( neighbors[0], 3 );

 	neighbors.clear();
 	val<<0.99;
 	nest.get_neighbors_unitcell( val, r, back_it );
 	ASSERT_EQ( neighbors.size(), 2 );
 	ASSERT_EQ( neighbors[0], 2 );
 	ASSERT_EQ( neighbors[1], 3 );

 	neighbors.clear();
 	val<<0.6;
 	nest.get_neighbors_unitcell( val, r, back_it );
 	ASSERT_EQ( neighbors.size(), 3 );
 	ASSERT_EQ( neighbors[0], 1 );
 	ASSERT_EQ( neighbors[1], 2 );
 	ASSERT_EQ( neighbors[2], 3 );
}

TEST(NEST_NEIGHBOR,unit_2d_boundary){
	NEST<2,RowVector2d> nest;
	std::vector<size_t> neighbors;
 	std::back_insert_iterator< std::vector<size_t> > back_it(neighbors);
 	size_t r = 2;

 	neighbors.clear();
 	nest.get_neighbors_unitcell( NEST<2>::ValueType(1.5,1.5), r, back_it );
 	ASSERT_EQ( neighbors.size(), 0 );
 	neighbors.clear();
 	nest.get_neighbors_unitcell( NEST<2>::ValueType(-1.5,1.5), r, back_it );
 	ASSERT_EQ( neighbors.size(), 0 );
 	neighbors.clear();
 	nest.get_neighbors_unitcell( NEST<2>::ValueType(1.5,-1.5), r, back_it );
 	ASSERT_EQ( neighbors.size(), 0 );
 	neighbors.clear();
 	nest.get_neighbors_unitcell( NEST<2>::ValueType(-1.5,-1.5), r, back_it );
 	ASSERT_EQ( neighbors.size(), 0 );

 	neighbors.clear();
 	nest.get_neighbors_unitcell( NEST<2>::ValueType(1.1,1.1), r, back_it );
 	// BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" << nest.set_and_get(i,r)<<") "; } cout << endl;
 	ASSERT_EQ( neighbors.size(), 1 );
 	ASSERT_EQ( nest.set_and_get(neighbors[0],r), RowVector2d(0.875,0.875) );

 	neighbors.clear();
 	nest.get_neighbors_unitcell( NEST<2>::ValueType(0.9,1.1), r, back_it );
 	// BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" << nest.set_and_get(i,r)<<") "; } cout << endl;
 	ASSERT_EQ( neighbors.size(), 2 );
 	ASSERT_EQ( nest.set_and_get(neighbors[0],r), RowVector2d(0.625,0.875) );
 	ASSERT_EQ( nest.set_and_get(neighbors[1],r), RowVector2d(0.875,0.875) );

 	neighbors.clear();
 	nest.get_neighbors_unitcell( NEST<2>::ValueType(1.1,0.6), r, back_it );
 	// BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" << nest.set_and_get(i,r)<<") "; } cout << endl;
 	ASSERT_EQ( neighbors.size(), 3 );
 	ASSERT_EQ( nest.set_and_get(neighbors[0],r), RowVector2d(0.875,0.375) );
 	ASSERT_EQ( nest.set_and_get(neighbors[1],r), RowVector2d(0.875,0.625) );
 	ASSERT_EQ( nest.set_and_get(neighbors[2],r), RowVector2d(0.875,0.875) );

 	neighbors.clear();
 	nest.get_neighbors_unitcell( NEST<2>::ValueType(-0.249,-0.249), r, back_it );
 	// BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" << nest.set_and_get(i,r)<<") "; } cout << endl;
 	ASSERT_EQ( neighbors.size(), 1 );
 	ASSERT_EQ( nest.set_and_get(neighbors[0],r), RowVector2d(0.125,0.125) );

 	neighbors.clear();
 	nest.get_neighbors_unitcell( NEST<2>::ValueType(-0.251,-0.251), r, back_it );
 	// BOOST_FOREACH(size_t i,neighbors){ cout << i << "(" << nest.set_and_get(i,r)<<") "; } cout << endl;
 	ASSERT_EQ( neighbors.size(), 0 );

}


/*
template<class NestType>
void generic_test_neighbors(
	NestType nest,
	typename NestType::Params lb,
	typename NestType::Params ub
){
	typedef typename NestType::IndexType Index;
	typedef std::vector<typename NestType::IndexType> IndexVec;
	typedef std::set<typename NestType::IndexType> IndexSet;
	boost::random::mt19937 rng(time(0));
	boost::uniform_real<> uniform;
	typename NestType::ValueType randpt;
	for(size_t r = 0; r <= 9; ++r){
		// typename NestType::FloatType maxdis = 0;
		for(size_t iter=0; iter < 10/NestType::DIMENSION; ++iter){
			for(size_t i = 0; i < NestType::DIMENSION; ++i){
				randpt[i] = (ub[i]-lb[i])*uniform(rng)+lb[i];
			}
			IndexSet nbset;
			Index center_index = nest.get_index(randpt,r);
	 		nest.get_neighbors( randpt, r, std::inserter(nbset,nbset.end()) );
	 		BOOST_FOREACH(Index i,nbset){
	 			cout << r << " " << center_index << " " << i << endl;
	 		}
	 		cout << endl;
	 		
			// double dist = (randpt-nest.value()).norm();
			// // cout << randpt.transpose() << " " << nest.value().transpose() << endl;
			// maxdis = fmax(maxdis,dist);
			// ASSERT_LE( dist, nest.covering_radius(r) );
		}
		// covering radius should be reasonably tigth
		// EXPECT_LT( nest.covering_radius(r)*0.75 , maxdis );
		// cout << DIM << " " << maxdis << " " << nest.covering_radius(r) << std::endl;
	}

}

template<int DIM>
void generic_test_neighbors_unit(){
	typedef NEST<DIM> NestType;
	typename NestType::Params lb,ub;
	lb.fill(0);
	ub.fill(1);
	generic_test_neighbors( NestType(), lb, ub );
}


TEST(NEST_NEIGHBOR,neighbor_radius_unit){
	generic_test_neighbors_unit<2>();
}
*/

}
}

