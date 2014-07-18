#ifndef INCLUDED_util_meta_util_HH
#define INCLUDED_util_meta_util_HH

#include <iostream>
#include <boost/fusion/include/for_each.hpp>

namespace scheme {
namespace util {
namespace meta {

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
		boost::fusion::for_each( x.second, PrintInstanceType(out,indent+"    ") );
	}
};

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
