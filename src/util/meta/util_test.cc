#include <gtest/gtest.h>
#include <util/meta/util.hh>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>

namespace scheme {
namespace util {
namespace meta {

// TEST(META,meta_product){
// 	typedef m::vector<char,int> Types;
// 	typedef m::vector<float,double> Types2;	
// 	m::for_each<Types>(PrintType());
// 	std::cout << "=========" << std::endl;
// 	m::for_each< product1<Types2>::apply<char>::type >(PrintType());
// 	m::for_each< product1<Types2>::apply<int >::type >(PrintType());
// 	std::cout << "=========" << std::endl;
// 	// typedef m::at< product<Types,Types2>::SeqOfSeq, m::int_<0> > PROD;
// 	// m::for_each< PROD >( PrintType() );

// }

}
}
}

