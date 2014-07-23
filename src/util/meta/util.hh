#ifndef INCLUDED_util_meta_util_HH
#define INCLUDED_util_meta_util_HH

#include <iostream>
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/is_sequence.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/transform.hpp>
// #include <boost/mpl/iter_fold.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/is_sequence.hpp>

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
	bool newline;


	PrintType(
		std::ostream & _out=std::cout,
		std::string _indent="",
		bool n=true
	) : 
		out(_out),indent(_indent),newline(n) {}

	template <typename T>
	typename boost::disable_if_c<f::traits::is_sequence<T>::value>::type
	operator()( T const & ) const {
		out << indent << typeid(T).name();
		if(newline) out << std::endl; else out << " ";
	}
	template <typename T>
	typename boost::enable_if_c<f::traits::is_sequence<T>::value>::type
	operator()(T const &) const {
		out << "TypeSeq<";
		if(newline) out << std::endl; else out << " ";		
		m::for_each<T>( PrintType(out,newline?indent+"    ":"",newline) );
		out << (newline?"  ":"") << ">";
		if(newline) out << std::endl; else out << " ";		
	}
	template <typename T>
	void operator()( type2type<T> const & = type2type<T>() ) const {
		out << indent << typeid(T).name();
		if(newline) out << std::endl; else out << " ";
	}
	template <typename A,typename B>
	void operator()(m::pair<A,B> const &) const {
		out << indent << "mpl::pair< " << typeid(A).name() << ", " << typeid(B).name() << " > ";
		if(newline) out << std::endl; else out << " ";
	}
	template <typename A,typename B>
	void operator()(std::pair<A,B> const &) const {
		out << indent << "std::pair< " << typeid(A).name() << ", " << typeid(B).name() << " > ";
		if(newline) out << std::endl; else out << " ";		
	}
	template <typename A,typename B>
	void operator()(f::pair<A,B> const &) const {
		out << indent << "fusion::pair< " << typeid(A).name() << ", " << typeid(B).name() << " > ";
		if(newline) out << std::endl; else out << " ";
	}



	// template <typename A,typename B>
	// void operator()(m::vector<A,B> const &) const {
	// 	PrintType p(out,indent,false);
	// 	out << indent << "mpl::vector< ";
	// 	p(A()); out << ", "; p(B()); out << " >";
	// 	if(newline) out << std::endl;
	// }
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

// utility to get MEMBER if it exists, else double                             
struct __EMPTY_TYPE_UTILITY__ {};
template<typename T> struct __TO_EMPTY_TYPE_UTILITY__ { typedef __EMPTY_TYPE_UTILITY__ type; };

#define SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(MEMBER,DEFAULT)                                  \
	template<typename T, typename Enable = scheme::util::meta::__EMPTY_TYPE_UTILITY__>       \
		struct get_ ## MEMBER ## _ ## DEFAULT                                                \
			{ typedef DEFAULT type; };                                                       \
	template<typename T>                                                                     \
		struct get_ ## MEMBER ## _ ## DEFAULT<                                               \
		 T, typename scheme::util::meta::__TO_EMPTY_TYPE_UTILITY__<typename T::FOO>::type >  \
			{ typedef typename T::MEMBER type; };

#define SCHEME_MEMBER_TYPE_DEFAULT_SELF_TEMPLATE(MEMBER)                                     \
	template<typename T, typename Enable = scheme::util::meta::__EMPTY_TYPE_UTILITY__>       \
		struct get_ ## MEMBER ## _SELF                                                       \
			{ typedef T type; };                                                             \
	template<typename T>                                                                     \
		struct get_ ## MEMBER ## _SELF<                                                      \
		 T, typename scheme::util::meta::__TO_EMPTY_TYPE_UTILITY__<typename T::FOO>::type >  \
			{ typedef typename T::MEMBER type; };



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




}
}
}

namespace std {
	template<class A,class B>
	std::ostream & operator<<(std::ostream & out, std::pair<A,B> const & p){ 
		// return out << "std::pair< " << typeid(A).name() << ", " << typeid(B).name() << " >( " << p.first << ", " << p.second << " )";
		return out << p.first << ", " << p.second;		
	}
}



#endif
