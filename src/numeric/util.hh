#ifndef INCLUDED_numeric_util_HH
#define INCLUDED_numeric_util_HH


namespace scheme { namespace numeric {


template<class Position>
bool approx_eq(Position const & a, Position const & b){
	return a.isApprox(b);
}

}}

#endif
