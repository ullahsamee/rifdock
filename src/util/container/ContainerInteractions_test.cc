#include <gtest/gtest.h>

#include <util/container/ContainerInteractions.hh>


namespace scheme {
namespace util {
namespace container {

TEST(util_container_ContainerInteractions,default_on_vector){
	typedef uint32_t Index;
	typedef ContainerInteractions< std::vector<int>, std::vector<int>, Index > CI;
	typedef std::pair<Index,Index> Index2;

	BOOST_STATIC_ASSERT( boost::is_same< 
		std::pair<ContainerInteractionsIter<Index>,ContainerInteractionsIter<Index> >,
		CI::Range 
	>::value );

	std::vector<int> a(2,0);
	std::vector<int> b(3,1);

	CI::Range r = CI::get_interactions(a,b);
	// BOOST_FOREACH( Index2 p, r ) cout << p << endl;
	get_citer<CI::Range>::type beg = get_cbegin(r);
	get_citer<CI::Range>::type end = get_cend(r);
	ASSERT_EQ( *beg++, Index2(0,0) );
	ASSERT_EQ( *beg++, Index2(0,1) );
	ASSERT_EQ( *beg++, Index2(0,2) );
	ASSERT_EQ( *beg++, Index2(1,0) );
	ASSERT_EQ( *beg++, Index2(1,1) );
	ASSERT_EQ( *beg++, Index2(1,2) );
	ASSERT_EQ( beg, end );

}


}}}
namespace scheme { namespace util { namespace container {

	template<class T> struct myvector : std::vector<T> {
		myvector(size_t N,T def) : std::vector<T>(N,def){}
	};
	template<class A, class Container2, class Index>
	struct ContainerInteractions<myvector<A>,Container2,Index> {
		typedef std::vector<std::pair<Index,Index> > Range;
		static Range get_interactions(myvector<A> const & c1, Container2 const & c2){
			// std::cout << "ContainerInteractions specialization myvector10" << std::endl;
			Range r;
			for(Index i1 = 0; i1 < c1.size(); ++i1)
			for(Index i2 = 0; i2 < c2.size(); ++i2)	
				r.push_back(std::make_pair(i1,i2));			
			return r;
		}
	};

}}}


namespace scheme { 
namespace util {
namespace container {

TEST(util_container_ContainerInteractions,specialization_on_myvector){
	typedef uint16_t Index;
	typedef ContainerInteractions< myvector<int>, std::vector<int>, Index > CI;
	typedef std::pair<Index,Index> Index2;

	BOOST_STATIC_ASSERT( boost::is_same< 
		std::vector<Index2>,
		CI::Range 
	>::value );

	myvector<int> a(2,0);
	std::vector<int> b(3,1);

	CI::Range r = CI::get_interactions(a,b);
	// BOOST_FOREACH( Index2 p, r ) cout << p << endl;
	get_citer<CI::Range>::type beg = get_cbegin(r);
	get_citer<CI::Range>::type end = get_cend(r);
	ASSERT_EQ( *beg++, Index2(0,0) );
	ASSERT_EQ( *beg++, Index2(0,1) );
	ASSERT_EQ( *beg++, Index2(0,2) );
	ASSERT_EQ( *beg++, Index2(1,0) );
	ASSERT_EQ( *beg++, Index2(1,1) );
	ASSERT_EQ( *beg++, Index2(1,2) );
	ASSERT_EQ( beg, end );

	ASSERT_EQ( r.size(), 6 );
}

}}}
