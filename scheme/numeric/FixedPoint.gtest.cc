#include <gtest/gtest.h>

#include "scheme/numeric/FixedPoint.hh"

namespace scheme { namespace numeric { namespace test {

using std::cout;
using std::endl;

TEST( FixedPoint, test ){

	FixedPoint<-17> f;
	for(double d = 0; d > -15.0; d -= 0.1){
		f = d;
		ASSERT_EQ( (double)f, (int)(d*-17.0)/-17.0 );
	}
	f = 1.0;
	ASSERT_EQ( (double)f, 0.0 );
	f = -999.0;
	ASSERT_EQ( (double)f, -255.0/17.0 );	

}


}}}

