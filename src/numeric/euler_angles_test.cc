#include <gtest/gtest.h>

#include "numeric/euler_angles.hh"

#include <Eigen/Geometry>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

namespace scheme { namespace numeric { namespace test {

using std::cout;
using std::endl;

TEST(euler_angles,test){
	using namespace Eigen;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> gauss;
	boost::uniform_real<> uniform;

	for(int i = 0; i < 10000; ++i){
		Quaterniond q( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
		q.normalize();
		Matrix3d m = q.matrix();
		Vector3d euler;
		euler_angles(m,euler);
		// cout << euler.transpose() << endl;
		// cout << m << endl;
		// cout << endl;
		Matrix3d m2;
		from_euler_angles(euler,m2);
		double thresh = 0.0000001;
		if( euler[2] < 0.000001 || euler[2] > M_PI-0.000001 ) thresh = 0.0002;
		if( !m2.isApprox(m,thresh) ){
			cout << euler << endl;
			cout << m << endl;
			cout << m2 << endl;
		}
		ASSERT_TRUE( m2.isApprox(m,thresh) );

	}
}



}}}

