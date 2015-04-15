#include <gtest/gtest.h>

#include "scheme/objective/storage/RotamerScores.hh"

#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>

namespace scheme { namespace objective { namespace storage { namespace rstest {

using std::cout;
using std::endl;

TEST( RotamerScores, test_store_1 ){

	RotamerScores<1> rs;
	rs.add_rotamer( 0, -1.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.0 );
	rs.add_rotamer( 0, -0.9 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.0 ); // not lower
	rs.add_rotamer( 0, -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 ); // is lower

	rs.add_rotamer( 1, -1.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ),  0.0 );
	rs.add_rotamer( 1, -1.25 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ),  0.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );
	rs.add_rotamer( 0, -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ),  0.0 );

	rs.add_rotamer( 7, -2.5 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ),  0.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ),  0.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 7 ), -2.5 );	


}

TEST( RotamerScores, test_store_2 ){

	RotamerScores<2> rs;
	rs.add_rotamer( 0, -1.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.0 );
	rs.add_rotamer( 0, -0.9 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.0 ); // not lower
	rs.add_rotamer( 0, -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 ); // is lower

	rs.add_rotamer( 1, -1.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.0 );
	rs.add_rotamer( 1, -1.25 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );
	rs.add_rotamer( 0, -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );

	rs.add_rotamer( 7, -2.5 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ),  0.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 7 ), -2.5 );	


}

TEST( RotamerScores, test_store_3 ){

	RotamerScores<3> rs;
	rs.add_rotamer( 0, -1.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.0 );
	rs.add_rotamer( 0, -0.9 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.0 ); // not lower
	rs.add_rotamer( 0, -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 ); // is lower

	rs.add_rotamer( 1, -1.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.0 );
	rs.add_rotamer( 1, -1.25 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );
	rs.add_rotamer( 0, -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );

	rs.add_rotamer( 7, -2.5 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 7 ), -2.5 );	

	rs.add_rotamer( 4, -2.5 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ),  0.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 7 ), -2.5 );	
	ASSERT_FLOAT_EQ( rs.score_of_rot( 4 ), -2.5 );		

}

TEST( RotamerScores, test_store_4 ){

	RotamerScores<4> rs;

	ASSERT_EQ( rs.name(), "RotamerScores< 4, 2, 1024 >" );

	rs.add_rotamer( 0, -1.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.0 );
	rs.add_rotamer( 0, -0.9 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.0 ); // not lower
	rs.add_rotamer( 0, -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 ); // is lower

	rs.add_rotamer( 1, -1.0 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.0 );
	rs.add_rotamer( 1, -1.25 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.125 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );
	rs.add_rotamer( 0, -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );

	rs.add_rotamer( 7, -2.5 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 7 ), -2.5 );	

	rs.add_rotamer( 4, -2.5 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 0 ), -1.50 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 1 ), -1.25 );
	ASSERT_FLOAT_EQ( rs.score_of_rot( 7 ), -2.5 );	
	ASSERT_FLOAT_EQ( rs.score_of_rot( 4 ), -2.5 );	


}

TEST( RotamerScores, sort ){
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::uniform_real<> uniform;

	for( int iter = 0; iter < 100; ++iter ){

		RotamerScores<10> rs;
		for( int i = 0; i < rs.size(); ++i ){
			rs.add_rotamer( i, uniform(rng) );
		}
		rs.sort_rotamers();

		for( int i = 1; i < rs.size(); ++i ){
			ASSERT_LE( rs.scores_[i-1], rs.scores_[i] );
		}

	}

}




}}}}
