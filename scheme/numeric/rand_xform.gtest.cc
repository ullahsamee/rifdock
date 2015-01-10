#include <gtest/gtest.h>

#include "scheme/numeric/rand_xform.hh"
#include "scheme/numeric/bcc_lattice.hh"
#include "scheme/numeric/util.hh"
#include <boost/timer/timer.hpp>
#include <map>

namespace scheme { namespace numeric { namespace rand_xfrom_test {

using std::cout;
using std::endl;

typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;

	// boost::random::mt19937 & rng,
	// Eigen::Transform<T,3,Eigen::AffineCompact> & x,
	// double cart_bound, double quat_bound

TEST( rand_xform, cart_correctness ){
	int NSAMP = 1000;
	#ifdef NDEBUG
	NSAMP = 10000;
	#endif

	boost::random::mt19937 rng((unsigned int)time(0) + 934875);

	double cart_bound=6.0, max_cart=0.0;
	Xform x;

	for(int i = 0; i < NSAMP; ++i){
		rand_xform_quat( rng, x, cart_bound, 0.0 );
		EXPECT_LE( x.translation().norm(), cart_bound );
		max_cart = std::max( max_cart, x.translation().norm() );
	}
	EXPECT_LE( cart_bound*0.99,  max_cart );

}


TEST( rand_xform, ori_correctness ){
	int NSAMP = 10;
	#ifdef NDEBUG
	NSAMP = 1000;
	#endif

	boost::random::mt19937 rng((unsigned int)time(0) + 934875);
	boost::uniform_real<> runif;

	for(int a = 0; a < NSAMP; ++a){

		double ANG = 0.1 + runif(rng)*119.8;
		double quat_bound=numeric::deg2quat(ANG) , max_ang=0.0;

		Xform x;
		for(int i = 0; i < 1000; ++i){
			rand_xform_quat( rng, x, 0.0, quat_bound );
			double ang = Eigen::AngleAxisd( x.rotation() ).angle() * 180.0 / M_PI;
			max_ang = std::max( max_ang, ang );
			ASSERT_LE( ang, ANG );
		}
		ASSERT_LE( ANG*0.99, max_ang );

	}


}



}}}
