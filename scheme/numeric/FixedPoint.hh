#ifndef INCLUDED_numeric_FixedPoint_HH
#define INCLUDED_numeric_FixedPoint_HH

namespace scheme { namespace numeric {

template<
	int Divisor,
	class Data=unsigned char
>
struct FixedPoint {
	Data data_;
	FixedPoint(){ data_ = 0; }
	FixedPoint(double d){ data_ = d*(double)Divisor; }
	FixedPoint(float  d){ data_ = d*(float )Divisor; }
	void operator =(double d){ data_ = d*(double)Divisor; }
	void operator =(float  d){ data_ = d*(float )Divisor; }	
	operator double(){ return (double)data_/(double)Divisor; }
	operator float (){ return (float )data_/(float )Divisor; }
	operator int(){ return (int)data_; }
};

}}


#endif
