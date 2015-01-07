#include <gtest/gtest.h>

#include "scheme/objective/hash/XformMap.hh"
#include "scheme/numeric/rand_xform.hh"
#include <Eigen/Geometry>

#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>
#include <sparsehash/dense_hash_set>

#include <fstream>

namespace scheme { namespace objective { namespace hash { namespace xmtest {

using std::cout;
using std::endl;

typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;


TEST( XformMap, stores_correctly ){
	boost::random::mt19937 rng((unsigned int)time(0) + 296720384);
	boost::uniform_real<> runif;

	XformMap< Xform, double > xmap( 0.5, 10.0 );
	std::vector< std::pair<Xform,double> > dat;
	for(int i = 0; i < 1000; ++i){
		Xform x;
		numeric::rand_xform( rng, x, 256.0 );
		double val = runif(rng);
		xmap[x] = val;
		dat.push_back( std::make_pair(x,val) );
	}

	XformMap< Xform, double > const & xmap_test( xmap );
	for(int i = 0; i < dat.size(); ++i){
		Xform const & x = dat[i].first;
		double const v = dat[i].second;
		ASSERT_EQ( xmap_test[x], v );
		// cout << x.translation().transpose() << " " << v << endl;
	}

	// { // no way to check if stream in binary!
	// 	std::cout << "following failure message is expected" << std::endl;
	// 	std::ofstream out("test.sxm" );// , std::ios::binary );
	// 	ASSERT_FALSE( xmap.save( out ) );
	// 	out.close();
	// }
	std::ofstream out("test.sxm" , std::ios::binary );
	ASSERT_TRUE( xmap.save( out ) );
	out.close();

	XformMap< Xform, double > xmap_loaded;
	std::ifstream in( "test.sxm"  , std::ios::binary );
	ASSERT_TRUE( xmap_loaded.load( in ) );
	in.close();

	ASSERT_EQ( xmap.cart_resl_, xmap_loaded.cart_resl_ );
	ASSERT_EQ( xmap.ang_resl_, xmap_loaded.ang_resl_ );	
	for(int i = 0; i < dat.size(); ++i){
		Xform const & x = dat[i].first;
		double const v = dat[i].second;
		ASSERT_EQ( xmap.hasher_.get_key(x) , xmap_loaded.hasher_.get_key(x) );
		ASSERT_EQ( xmap_loaded[x], v );
		// cout << x.translation().transpose() << " " << v << endl;
	}


}


}}}}
