#ifndef INCLUDED_scheme_nest_NEST_HH
#define INCLUDED_scheme_nest_NEST_HH

#include <Eigen/Dense>

namespace scheme {
namespace nest {

template< int DIM, class Float >
struct DefaultFilter {

};

template< int DIM, class Value, class Float >
struct DefaultAdaptor {

};
 

template< int DIM,
	class Float   = double, 
	class Value   = Eigen::Array<Float,DIM,1>,
	class Adaptor = DefaultAdaptor<DIM,Value,Float>, 
	class Filter  = DefaultFilter<DIM,Float>,
	class Index   = size_t
>
struct NEST {
	static Index const ONE = 1;
	Eigen::Array<Index,DIM,1> base_size;
	Eigen::Array<Float,DIM,1> lower_bound;
	Eigen::Array<Float,DIM,1> upper_bound;
	NEST() :
		base_size(1),
		lower_bound(0.0),
		upper_bound(1.0)
	{}
	Index size(Index resl) const {
		return base_size.prod() * ONE<<resl;
	}
};


}
}

#endif
