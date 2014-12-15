#ifndef INCLUDED_numeric_FixedPoint_HH
#define INCLUDED_numeric_FixedPoint_HH

#include <limits>

namespace scheme { namespace numeric {

template<
	int Divisor,
	class Data=unsigned char
>
struct FixedPoint {
	typedef std::numeric_limits<Data> NL;
	Data data_;
	FixedPoint(){ data_ = 0; }
	FixedPoint(double d){ *this = d; }
	FixedPoint(float  d){ *this = d; }
	void operator =(double d){ data_ = std::max( (int)NL::min(), std::min( (int)NL::max(), (int)(d*(double)Divisor) ) ); }
	// void operator =(float  d){ (Data)(data_ = d*(float )Divisor); }
	operator double(){ return (double)data_/(double)Divisor; }
	// operator float (){ return (float )data_/(float )Divisor; }
	operator int(){ return (int)data_; }
};

}}


#endif
