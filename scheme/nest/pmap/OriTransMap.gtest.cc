
#include "scheme/nest/pmap/OriTransMap.hh"
#include "scheme/nest/NEST.hh"
#include "scheme/nest/NEST_concepts.hh"
#include "scheme/nest/NEST_test_util.hh"
#include <gtest/gtest.h>
#include <boost/random.hpp>

namespace scheme {
namespace nest {
namespace pmap {

using std::cout;
using std::endl;

TEST( OriTransMap, basic_test ){

	typedef Eigen::Transform<double,3,Eigen::AffineCompact> EigenXform;

	typedef NEST<6,EigenXform,OriTransMap> Nest;

	Nest nest( 999.0, 0.0, 100.0, 1 );

	std::cout << "THIS IS NO TEST!" << std::endl;
	for(int i = 0; i < nest.size(1)/24; ++i){
		bool valid = nest.set_state( i, 1 );
		if( !valid ) continue;
		EigenXform x = nest.value();
		Eigen::AngleAxisd aa( x.rotation() );
		std::cout << i << " " << x.translation().transpose() << " " << aa.angle() << " " << aa.axis().transpose() << std::endl;
	}


}


}
}
}
