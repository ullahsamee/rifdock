#include <gtest/gtest.h>

#include <util/meta/InstanceMap.hh>

#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/for_each.hpp>

namespace scheme {
namespace util {
namespace meta {

using std::cout;
using std::endl;

namespace bf = boost::fusion;
namespace mpl = boost::mpl;

struct PrintType {
	template<typename T>
	void operator()(T const & t){
		cout << typeid(t).name() << endl;
	}
};

struct PrintInstanceType {
	template <typename T>
	void operator()(T const & x) const {
		std::cout << typeid(x).name() << ": " << x << std::endl;
	}
};

TEST(InstanceMap,holds_types){
	using mpl::_1;
	using mpl::_2;

	typedef bf::map<
		bf::pair<int, char>
	  , bf::pair<double, std::string> >
	map_type;

	map_type m(
		bf::make_pair<int>('X')
	  , bf::make_pair<double>("Men"));

	ASSERT_EQ( 'X', bf::at_key<int>(m) );
	ASSERT_EQ( "Men", bf::at_key<double>(m) );

	typedef mpl::vector<int,char,float> Types;
	typedef mpl::vector<char,float,int> Types2;	

	{
		typedef InstanceMap<Types,_1> TEST;
		TEST imap;
		BOOST_STATIC_ASSERT( mpl::size<TEST>::value == mpl::size<Types>::value );
		BOOST_STATIC_ASSERT( bf::result_of::has_key<TEST,int>::value );
		ASSERT_TRUE( bf::has_key<int>(imap) );
		BOOST_STATIC_ASSERT( boost::is_same<int,
			bf::result_of::value_at_key<TEST,int>::type >::value );
		BOOST_STATIC_ASSERT( boost::is_same<char,
			bf::result_of::value_at_key<TEST,char>::type >::value );
		imap.get<int>() = 1;
		imap.get<char>() = 'C';	
		imap.get<float>() = 1.2345f;
		ASSERT_EQ( imap.get<int>(), 1 );
		ASSERT_EQ( imap.get<char>(), 'C' );
		ASSERT_EQ( imap.get<float>(), 1.2345f );
		// bf::for_each(imap,PrintInstanceType());
	}

	{
		InstanceMap<Types> imap;
		imap.get<int>() = 1;
		imap.get<char>() = 'C';	
		imap.get<float>() = 1.2345f;
		ASSERT_EQ( imap.get<int>(), 1 );
		ASSERT_EQ( imap.get<char>(), 'C' );
		ASSERT_EQ( imap.get<float>(), 1.2345f );
	}
	{
		InstanceMap<Types,Types> imap;
		imap.get<int>() = 1;
		imap.get<char>() = 'C';	
		imap.get<float>() = 1.2345f;
		ASSERT_EQ( imap.get<int>(), 1 );
		ASSERT_EQ( imap.get<char>(), 'C' );
		ASSERT_EQ( imap.get<float>(), 1.2345f );
	}
	{
		InstanceMap<Types,Types2> imap;
		imap.get<float>() = 1;
		imap.get<int>() = 'C';	
		imap.get<char>() = 1.2345f;
		ASSERT_EQ( imap.get<float>(), 1 );
		ASSERT_EQ( imap.get<int>(), 'C' );
		ASSERT_EQ( imap.get<char>(), 1.2345f );
	}
	{
		InstanceMap<Types,_1> imap;
		imap.get<int>() = 1;
		imap.get<char>() = 'C';	
		imap.get<float>() = 1.2345f;
		ASSERT_EQ( imap.get<int>(), 1 );
		ASSERT_EQ( imap.get<char>(), 'C' );
		ASSERT_EQ( imap.get<float>(), 1.2345f );
	}
	{
		InstanceMap<Types,std::vector<_1> > imap;
		imap.get<int>().push_back(1);
		imap.get<char>().push_back('C');
		imap.get<float>().push_back(1.2345f);
		ASSERT_EQ( imap.get<int>()[0], 1 );
		ASSERT_EQ( imap.get<char>()[0], 'C' );
		ASSERT_EQ( imap.get<float>()[0], 1.2345f );
		imap.get<float>().push_back(1.2345f);
		imap.get<float>().push_back(1.2345f);
		ASSERT_EQ( imap.get<float>().size(),3);

	}
	// bf::map< bf::pair<int,int>, bf::pair<char,char>, bf::pair<float,float> > test;
	// bf::for_each( test, PrintInstanceType() );

}



}
}
}
