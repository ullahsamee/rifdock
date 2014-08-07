#include <gtest/gtest.h>

#include "objective/voxel/VoxelArray.hh"

#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>

namespace scheme { namespace objective { namespace voxel { namespace test {

using std::cout;
using std::endl;

TEST(VoxelArray,simple1d){
	typedef util::SimpleArray<1,float> F1;
	// VoxelArray<0,float> a; // this should fail to compile

	VoxelArray<1,float> a1( 0, 3, 0.5 );
	ASSERT_EQ(a1.size(),7);
	a1[F1(0)] = 1;
	ASSERT_EQ(a1[F1(0)],1);
	#ifdef DEBUG
		// ASSERT_DEATH( a1[F1(-0.1)], ".*" );	
		a1[F1(0)];
		a1[F1(2.999)];		
		// ASSERT_DEATH( a1[F1(3.0)], ".*" );	
	#endif
}
TEST(VoxelArray,simple3d){
	typedef util::SimpleArray<3,float> F3;		
	VoxelArray<3,float> a3( F3(-1,-2,-3), F3(1,2,3), 0.5 );
	ASSERT_EQ(a3.shape()[0],5);
	ASSERT_EQ(a3.shape()[1],9);
	ASSERT_EQ(a3.shape()[2],13);		
	a3[F3(0,0,0)] = 1;
	a3[F3(-1,1,2.3)] = 2;
	ASSERT_EQ(a3[F3(0,0,0)],1);
	#ifdef DEBUG
		// ASSERT_DEATH( a3[F3(1,2,3)], ".*" );
		// ASSERT_DEATH( a3[F3(-1,-2,-3.1)], ".*" );			
	#endif
	ASSERT_EQ( a3[F3(0,0,0)] , 1 );
	ASSERT_EQ( a3[F3(-1,1,2.3)] , 2 );

}
TEST(VoxelArray,bounds3d){
	typedef util::SimpleArray<3,float> F3;		
	VoxelArray<3,float> a3( F3(-1,-2,-3), F3(1,2,3), 0.49 );
	ASSERT_EQ(a3.shape()[0],5);
	ASSERT_EQ(a3.shape()[1],9);
	ASSERT_EQ(a3.shape()[2],13);		
	a3[F3(0,0,0)] = 1;
	a3[F3(-1,1,2.3)] = 2;
	ASSERT_EQ(a3[F3(0,0,0)],1);
	#ifdef DEBUG
		// ASSERT_DEATH( a3[F3(1,2,3)], ".*" );
		// ASSERT_DEATH( a3[F3(-1,-2,-3.1)], ".*" );			
	#endif
	ASSERT_EQ( a3[F3(0,0,0)] , 1 );
	ASSERT_EQ( a3[F3(-1,1,2.3)] , 2 );

}

struct Ellipse3D {
	static double sqr(double f) { return f*f; }
	double c1,c2,c3,sd1,sd2,sd3;
	Ellipse3D() : c1(0),c2(0),c3(0),sd1(1),sd2(1),sd3(1) {}
	Ellipse3D(double a,double b,double c,double d,double e,double f) : c1(a),c2(b),c3(c),sd1(d),sd2(e),sd3(f) {}	
	double operator()(double f, double g, double h) const { 
		return 10.0+std::exp( -sqr((f-c1)/sd1) - sqr((g-c2)/sd2) - sqr((h-c3)/sd3) );
	}
	template<class F3> double operator()(F3 const & f3) const { return this->operator()(f3[0],f3[1],f3[2]); }
};

template<class Cache,class Field>
double test_cache_vs_field(Cache const & cache, Field field, int nsamp=10000){
	// cout << cache.num_elements()/1000000.0 <<"M" << endl;
	BOOST_STATIC_ASSERT(Cache::DIM==3);
	typedef util::SimpleArray<3,double> F3;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::uniform_real<> uniform;
	double maxerr = 0.0;
	for(int i = 0; i < nsamp; ++i){
		F3 samp = F3(uniform(rng),uniform(rng),uniform(rng)) * (cache.ub_-cache.lb_) + cache.lb_;
		double cacheval = cache[samp];
		double fieldval = field(samp);
		// std::cout << cacheval << " " << fieldval << std::endl;
		maxerr = std::max( std::abs(fieldval-cacheval), maxerr );
	}
	return maxerr;
}

TEST(FieldCache3D,test){

	Ellipse3D field(1,2,3,4,5,6);
	ASSERT_FLOAT_EQ( field(1,2,3), 11.0 );

	ASSERT_LE( test_cache_vs_field( FieldCache3D<double>(field,-11,13,1.6), field, 100000 ), 0.30 );
	ASSERT_LE( test_cache_vs_field( FieldCache3D<double>(field,-11,13,0.8), field, 100000 ), 0.15 );
	ASSERT_LE( test_cache_vs_field( FieldCache3D<double>(field,-11,13,0.4), field, 100000 ), 0.08 );
	ASSERT_LE( test_cache_vs_field( FieldCache3D<double>(field,-11,13,0.2), field, 100000 ), 0.04 );

}


}}}}
