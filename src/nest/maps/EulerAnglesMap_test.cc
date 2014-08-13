#include <gtest/gtest.h>

#include "nest/NEST.hh"
#include "nest/maps/EulerAnglesMap.hh"

#include <Eigen/Geometry>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

namespace scheme { namespace nest { namespace maps { namespace test {

using std::cout;
using std::endl;

TEST(EulerAnglesMap,test){
	using namespace Eigen;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> gauss;
	boost::uniform_real<> uniform;

		NEST<3,Matrix3d,EulerAnglesMap> nest;
	// 	for(int r = 1; r < 16; ++r){
	// 		double maxdiff = 0;
	// 		for(int i = 0; i < 1000000; ++i){
	// 			Eigen::Quaterniond q( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
	// 			q.normalize();
	// 			// cout << q << endl;
	// 			Eigen::Quaterniond qcen = nest.set_and_get( nest.get_index(q,r) , r );
	// 			maxdiff = std::max(maxdiff,q.angularDistance(qcen));
	// 		}
	// 		cout << boost::format("%2i %20i %.7d") % r % nest.size(r) % (maxdiff*180.0/M_PI) << endl;
	// 	}
	// }

}



}}}}

