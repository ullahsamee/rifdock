#include <gtest/gtest.h>

#include <types.hh>
#include <boost/any.hpp>

namespace scheme {
namespace nest {

using std::cout;
using std::endl;


TEST(MultiNest,boost_any){
	using boost::any;
	using boost::any_cast;

	double d = 123;
	double *dp = &d;
	shared_ptr<double> dsp = make_shared<double>();

	any ad = d;
	ASSERT_EQ( d, any_cast<double>(ad) );
	ASSERT_EQ( d, *any_cast<double>(&ad) );
	ASSERT_EQ( NULL, any_cast<float>(&ad) );	

	any adp = dp;

	d = 1.2345;	ASSERT_EQ( 1.2345, *any_cast<double*>(adp) );
	d = 3.2345;	ASSERT_EQ( 3.2345, *any_cast<double*>(adp) );
	d = 3.2345;	ASSERT_EQ( 3.2345, **any_cast<double*>(&adp) );
	d = 3.2345;	ASSERT_EQ( NULL, any_cast<float*>(&adp) );

	any adsp = dsp;
	*dsp = 2.34; ASSERT_EQ( 2.34, *any_cast<shared_ptr<double> >(adsp) );
	*dsp = 3.34; ASSERT_EQ( 3.34, **any_cast<shared_ptr<double> >(&adsp) );
	*dsp = 2.34; ASSERT_EQ( NULL, any_cast<double*>(&adsp) );

}

}
}

