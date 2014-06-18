
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

TEST(NEST_NEIGHBOR,neighbor_lookup_single_cell){
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
	 	NEST<3>::ValueType val;
	 	val << 0.85,0.15,0.5;
	 	size_t r = 3;
		std::vector<size_t> neighbors;
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

}
}

