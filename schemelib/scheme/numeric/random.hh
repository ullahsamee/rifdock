#ifndef INCLUDED_numeric_random_HH
#define INCLUDED_numeric_random_HH

#include "scheme/numeric/util.hh"
#include <Eigen/Geometry>
#include <random>

namespace scheme { namespace numeric {


// implements UniformRandomBitGenerator
class FastRandom {

public:
	typedef unsigned long result_type;

	result_type x, y, z;

	FastRandom( result_type seed1, result_type seed2, result_type seed3 ) 
	: x( seed1 ), y( seed2 ), z( seed3 )
	{}

	result_type min() {
		return 0;
	}

	result_type max() {
		return std::numeric_limits<unsigned long>::max();
	}
		//https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c
	unsigned long xorshf96(void) {          //period 2^96-1
	unsigned long t;
	    x ^= x << 16;
	    x ^= x >> 5;
	    x ^= x << 1;

		t = x;
		x = y;
		y = z;
		z = t ^ x ^ y;

		return z;
	}

	result_type operator()() {
		return xorshf96();
	}
};





}}


#endif
