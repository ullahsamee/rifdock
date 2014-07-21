#ifndef INCLUDED_util_meta_util_HH
#define INCLUDED_util_meta_util_HH

#include <iostream>
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/iter_fold.hpp>

namespace scheme {
namespace util {
namespace meta {

namespace m = boost::mpl;
namespace f = boost::fusion;


template<class T> struct type2type {};


template<class T> struct showclass;
template<int T> struct showint;

struct PrintType {
	std::ostream & out;
	std::string indent;
	PrintType(std::ostream & _out=std::cout,std::string _indent="") : out(_out),indent(_indent) {}
	template <typename T>
	void operator()(T const &) const {
		out << indent << typeid(T).name() << std::endl;
	}
	template <typename A,typename B>
	void operator()(m::pair<A,B> const &) const {
		out << indent << "mpl::pair< " << typeid(A).name() << ", " << typeid(B).name() << " > " << std::endl;		
	}
	template <typename A,typename B>
	void operator()(std::pair<A,B> const &) const {
		out << indent << "std::pair< " << typeid(A).name() << ", " << typeid(B).name() << " > " << std::endl;		
	}
	template <typename A,typename B>
	void operator()(f::pair<A,B> const &) const {
		out << indent << "fusion::pair< " << typeid(A).name() << ", " << typeid(B).name() << " > " << std::endl;		
	}
};


struct PrintInstanceType {
	std::ostream & out;
	std::string indent;
	PrintInstanceType(std::ostream & _out=std::cout,std::string _indent="") : out(_out),indent(_indent) {}
	template <typename T>
	void operator()(T const & x) const {
		out << indent << x << " (" << typeid(T).name() << ")" << std::endl;
	}
};

struct PrintBFMapofVec {
	std::ostream & out;
	std::string indent;
	PrintBFMapofVec(std::ostream & _out=std::cout,std::string _indent="") : out(_out),indent(_indent) {}
	template <typename T>
	void operator()(T const & x) const {
		out << indent << "KeyType: " << typeid(typename T::first_type).name() << std::endl;
		f::for_each( x.second, PrintInstanceType(out,indent+"    ") );
	}
};

template <typename SeqA, typename SeqB> 
struct intersect { 
	typedef typename 
		m::copy<
			SeqB ,
			m::inserter< m::set<>, m::insert< m::_1, m::_2 > > 
		>::type set_seq2;
	typedef typename 
		m::copy_if<
			SeqA,
			m::has_key< set_seq2, m::_1 > 
		>::type type;
};

/////////////// in progress.... ////////////////

// template< typename SeqOfSeq >
// struct flatten {
// 	typedef void type;
// };

// template< typename S, template<typename,typename> class pair = std::pair >
// struct product1 {
// 	template< typename T >
// 	struct apply {
// 		typedef typename
// 			m::transform< S, pair<T,m::_1> >::type
// 		type;
// 	};
// };

// template< typename S, typename T=S, template<typename,typename> class pair = std::pair >
// struct product {
// 	typedef typename m::transform< S, product1<T,pair> >::type SeqOfSeq;
// 	typedef typename flatten< SeqOfSeq >::type type;
// };



	// utility to get ResultType if it exists, else double
	// template<typename T> struct tovoid { typedef void type; };
	// template<typename T, typename Enable = void>
	// struct result_type_impl { typedef double type; };
	// template<typename T> 
	// struct result_type_impl< T, typename tovoid<typename T::ResultType>::type > {
	//     typedef typename T::ResultType  type; 
	// };


}
}
}

#endif
