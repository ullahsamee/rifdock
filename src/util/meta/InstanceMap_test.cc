#include <gtest/gtest.h>

#include <util/meta/InstanceMap.hh>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/for_each.hpp>

#include <util/SimpleArray.hh>

#include <set>

namespace scheme {
namespace util {
namespace meta {

using std::cout;
using std::endl;

namespace bf = boost::fusion;
namespace mpl = boost::mpl;

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

TEST(InstanceMap,can_use_fusion_pairs_directly){
	
	// "usual" way
	InstanceMap<
		m::vector< int , char >,
		m::vector< char, float >
	> zip_imap;
	zip_imap.get<int>() = 'a';
	zip_imap.get<char>() = 1.234f;
	ASSERT_EQ( zip_imap.get<int>(), 'a' );
	ASSERT_EQ( zip_imap.get<char>(), 1.234f );

	// make fusion map directly
	f::result_of::as_map<
		m::vector<
			f::pair<int,char>,
			f::pair<char,float>
		>
	>::type fmap;
	f::at_key<int>(fmap) = 'a';
	f::at_key<char>(fmap) = 1.234f;
	ASSERT_EQ( f::at_key<int >(fmap), 'a' );
	ASSERT_EQ( f::at_key<char>(fmap), 1.234f );

	InstanceMap<
		m::vector<
			f::pair<int,char>,
			f::pair<char,float>
		>,
		FUSION_PAIRS
	> imap;
	imap.get<int>() = 'a';
	imap.get<char>() = 1.234f;
	ASSERT_EQ( imap.get<int>(), 'a' );
	ASSERT_EQ( imap.get<char>(), 1.234f );
}

TEST(ContainerInstanceMap,basic_test){
	ContainerInstanceMap< m::vector<
		std::vector<int>,
		std::vector<char> 
	> > cmap;
	BOOST_STATIC_ASSERT( boost::is_same< std::vector<int>::value_type, int >::value );
	BOOST_STATIC_ASSERT( boost::is_same< impl::get_value_type_void<       int       >::type, void >::value );
	cmap.get<int>().push_back(1);
	cmap.get<char>().push_back('c');
	cmap.get<char>().push_back('h');
	ASSERT_EQ( cmap.get<int>().at(0), 1 );
	ASSERT_EQ( cmap.get<char>().at(0), 'c' );
	ASSERT_EQ( cmap.get<char>().at(1), 'h' );
}

TEST(ContainerInstanceMap,vector_default_test){
	ContainerInstanceMap< m::vector<
		int, // defaults to vector<int>
		char // defaults to vector<char>
	> > cmap;
	cmap.get<int>().push_back(1);
	cmap.get<char>().push_back('c');
	cmap.get<char>().push_back('h');
	ASSERT_EQ( cmap.get<int>().at(0), 1 );
	ASSERT_EQ( cmap.get<char>().at(0), 'c' );
	ASSERT_EQ( cmap.get<char>().at(1), 'h' );
}

TEST(ContainerInstanceMap,vector_default_set_test){
	ContainerInstanceMap< m::vector<
		int, // defaults to vector<int>
		std::set<char> 
	> > cmap;
	cmap.get<int>().push_back(1);
	cmap.get<char>().insert('c');
	cmap.get<char>().insert('h');
	ASSERT_EQ( cmap.get<int>().at(0), 1 );
	std::set<char> & tmp = cmap.get<char>();
	ASSERT_EQ( *tmp.find('c'), 'c' );
	ASSERT_EQ( *tmp.find('h'), 'h' );
	ASSERT_EQ( tmp.find('b'), tmp.end() );
}

TEST(ContainerInstanceMap,simple_array_test){
	ContainerInstanceMap< m::vector<
		SimpleArray<2,int>,
		SimpleArray<2,size_t>
	> > cmap;
	cmap.get<int>().at(0) = -1;
	cmap.get<size_t>()[0] = 1;
	cmap.get<size_t>()[1] = 2;
	ASSERT_EQ( cmap.get<int>().at(0), -1 );
	SimpleArray<2,size_t> & tmp = cmap.get<size_t>();
	ASSERT_EQ( tmp.at(0), (size_t)1 );
	ASSERT_EQ( tmp.at(1), (size_t)2 );
}

}
}
}
