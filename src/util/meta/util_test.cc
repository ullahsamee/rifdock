#include <gtest/gtest.h>

#include <util/meta/util.hh>
#include <util/meta/print_type.hh>
#include <util/meta/ref_wrap.hh>

#include <boost/tuple/tuple.hpp>
#include <boost/mpl/assert.hpp>

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

struct VAL { typedef int value_type; };

SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(FOO,int)
SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(FOO,float)
SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(value_type,void)

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

	BOOST_STATIC_ASSERT( boost::is_same< void   , get_value_type_void<EMPTY>::type >::value );
	BOOST_STATIC_ASSERT( boost::is_same< int    , get_value_type_void<VAL>::type >::value );

}

TEST(META_UTIL,PrintType){
	m::vector<int,char> m;

	std::ostringstream out;
	PrintType(out).operator()(m);

	PrintType(out,"",false)(m); out << std::endl;

	// std::cout << out.str() << std::endl;
}

TEST(META_UTIL,remove_ref){
	using boost::is_same;
	using boost::tuple;
	using boost::reference_wrapper;
	using std::pair;

	BOOST_MPL_ASSERT(( is_same< int, remove_refwrap<int>::type > ));
	BOOST_MPL_ASSERT(( is_same< int, remove_refwrap<reference_wrapper<int> >::type > ));
	BOOST_MPL_ASSERT_NOT(( is_same< int, reference_wrapper<int> > ));

	BOOST_MPL_ASSERT(( is_same< int, remove_refwrap<const int>::type > ));
	BOOST_MPL_ASSERT(( is_same< int, remove_refwrap<reference_wrapper<const int> >::type > ));
	BOOST_MPL_ASSERT_NOT(( is_same< int const, reference_wrapper<int const> > ));

	BOOST_MPL_ASSERT(( is_same< int, recursive_remove_refwrap<int>::type > ));
	BOOST_MPL_ASSERT(( is_same< int, recursive_remove_refwrap<reference_wrapper<int> >::type > ));
	BOOST_MPL_ASSERT(( is_same< pair<int,char>, recursive_remove_refwrap<pair<int,char> >::type > ));
	BOOST_MPL_ASSERT(( is_same< pair<int,char>, recursive_remove_refwrap<pair<reference_wrapper<int>,char> >::type > ));
	BOOST_MPL_ASSERT(( is_same< pair<int,char>, recursive_remove_refwrap<
		pair<  reference_wrapper<int>,  reference_wrapper<char>   > >::type > ));

	BOOST_MPL_ASSERT(( is_same<
		tuple<                   int        >, recursive_remove_refwrap<
		tuple< reference_wrapper<int>       > >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                  char  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<char> >      >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                  char ,                   float  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<char>, reference_wrapper<float> >      >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                  char ,                   float ,                   int  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<char>, reference_wrapper<float>, reference_wrapper<int> >      >::type > ));

	BOOST_MPL_ASSERT(( is_same<
		tuple<                   int        >, recursive_remove_refwrap<
		tuple< reference_wrapper<int>       > >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                        char  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<const char> >      >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int       ,                  char ,                   float  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int const>,reference_wrapper<char>, reference_wrapper<float> >      >::type > ));
	BOOST_MPL_ASSERT(( is_same< 
		tuple<                   int ,                  char ,                   float       ,                   int  >  ,   recursive_remove_refwrap<
		tuple< reference_wrapper<int>,reference_wrapper<char>, reference_wrapper<float const>, reference_wrapper<int> >      >::type > ));

	// print_type< recursive_remove_refwrap< tuple< reference_wrapper<int>,reference_wrapper<char> > >::type >();

}



}
}
}

