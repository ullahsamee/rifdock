#include <gtest/gtest.h>

#include "scheme/objective/hash/XformMap.hh"
#include <Eigen/Geometry>

#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

namespace scheme { namespace objective { namespace hash { namespace test {

using std::cout;
using std::endl;

TEST( XformMap, basic_test ){
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;


	typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;

	for(int j = 0; j < 100; ++j){

		double  ang_resl = runif(rng)*19+1.0;
		double cart_resl = runif(rng)*2.0 + 0.125;
		XformMap<Xform> xm( cart_resl, ang_resl, 512.0 );

		// cout << ang_resl << " " << cart_resl << " " << xm.approx_size()/1000000000.0 << endl;

		for(int i = 0; i < 10*1000; ++i){
			Xform x = Xform::Identity();
			x.translation()[0] = runif(rng)*512.0-256.0;
			x.translation()[1] = runif(rng)*512.0-256.0;
			x.translation()[2] = runif(rng)*512.0-256.0;
			Eigen::Quaterniond qrand( rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng) );
			qrand.normalize();
			x.rotate( qrand.matrix() );

			Xform c = xm.get_center( xm.get_key(x) );
			Xform l = c.inverse() * x;
			double dt = l.translation().norm();
			double da = Eigen::AngleAxisd(l.rotation()).angle()*180.0/M_PI;
			ASSERT_LT( dt, cart_resl );
			ASSERT_LT( da,  ang_resl );
		}

	}
}

}}}}
