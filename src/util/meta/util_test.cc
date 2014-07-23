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



struct EMPTY {};
struct DOUBLE { typedef double FOO; };
struct CHAR { typedef char FOO; };

SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(FOO,int)
SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(FOO,float)

SCHEME_MEMBER_TYPE_DEFAULT_SELF_TEMPLATE(FOO)


TEST(META_UTIL,member_type_default){

	BOOST_STATIC_ASSERT( boost::is_same< int    , get_FOO_int<EMPTY  >::type >::value );
	BOOST_STATIC_ASSERT( boost::is_same< double , get_FOO_int<DOUBLE >::type >::value );
	BOOST_STATIC_ASSERT( boost::is_same< char   , get_FOO_int<CHAR   >::type >::value );

	BOOST_STATIC_ASSERT( boost::is_same< float  , get_FOO_float<EMPTY  >::type >::value );

	BOOST_STATIC_ASSERT( boost::is_same< EMPTY  , get_FOO_SELF<EMPTY  >::type >::value );
	BOOST_STATIC_ASSERT( boost::is_same< double , get_FOO_SELF<DOUBLE >::type >::value );
	BOOST_STATIC_ASSERT( boost::is_same< char   , get_FOO_SELF<CHAR   >::type >::value );

	BOOST_STATIC_ASSERT( boost::is_same< int    , get_FOO_SELF<int>::type >::value );

}

TEST(META_UTIL,PrintType){
	m::vector<int,char> m;

	PrintType()(m);

	PrintType(std::cout,"",false)(m); std::cout << std::endl;
}

}
}
}

