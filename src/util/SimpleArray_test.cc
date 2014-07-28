#include <gtest/gtest.h>
#include <util/SimpleArray.hh>
#include <boost/foreach.hpp>

namespace scheme {
namespace util {

using std::cout;
using std::endl;

TEST(SimpleArray,bounds_check){
	SimpleArray<3,int> a;
	a[3]; // non-bounds checked
	ASSERT_DEATH( a.at(3), "Assertion.*SimpleArray" );
}

TEST(SimpleArray,iteration){
	SimpleArray<3,int> a;
	int v;
	v = 0; BOOST_FOREACH(int & i, std::make_pair(a.begin(),a.end()) ) i = ++v;
	v = 0; BOOST_FOREACH(int i, std::make_pair(a.begin(),a.end()) ) ASSERT_EQ( ++v, i );
	v = 0; BOOST_FOREACH(int i, a ) ASSERT_EQ( ++v, i );
	SimpleArray<3,int> const & r = a;
	v = 0; BOOST_FOREACH(int i, std::make_pair(r.begin(),r.end()) ) ASSERT_EQ( ++v, i );	
	v = 0; BOOST_FOREACH(int const & i, r ) ASSERT_EQ( ++v, i );	
}

}
}
