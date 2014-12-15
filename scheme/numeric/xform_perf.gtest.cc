#include <gtest/gtest.h>

#include <Eigen/Geometry>

#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>

namespace scheme { namespace numeric { namespace pref_test {

using std::cout;
using std::endl;

typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;

template<class T>
void
rand_xform(
	boost::random::mt19937 & rng,
	boost::uniform_real<> & runif,
	boost::normal_distribution<> & rnorm,
	Eigen::Transform<T,3,Eigen::Affine> & x
){
	Eigen::Quaterniond qrand( rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng) );
	qrand.normalize();
	Eigen::Matrix3d m = qrand.matrix();

	x.data()[3] = 0; x.data()[7] = 0; x.data()[11] = 0; x.data()[15] = 1;
	for(int i = 0; i <  3; ++i) x.data()[i] = m.data()[i-0];
	for(int i = 4; i <  7; ++i) x.data()[i] = m.data()[i-1];
	for(int i = 8; i < 11; ++i) x.data()[i] = m.data()[i-2];
	x.data()[12] = runif(rng)*512.0-256.0;
	x.data()[13] = runif(rng)*512.0-256.0;
	x.data()[14] = runif(rng)*512.0-256.0;
}

template<class T>
void
rand_xform(
	boost::random::mt19937 & rng,
	boost::uniform_real<> & runif,
	boost::normal_distribution<> & rnorm,
	Eigen::Transform<T,3,Eigen::AffineCompact> & x
){
	Eigen::Quaterniond qrand( rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng) );
	qrand.normalize();
	Eigen::Matrix3d m = qrand.matrix();

	for(int i = 0; i < 9; ++i) x.data()[i] = m.data()[i];
	x.data()[ 9] = runif(rng)*512.0-256.0;
	x.data()[10] = runif(rng)*512.0-256.0;
	x.data()[11] = runif(rng)*512.0-256.0;
}

// TEST( XformMap, basic_test ){
template< class Xform >
void test_xform_perf(){
	int NSAMP = 100*1000;
	#ifdef NDEBUG
	NSAMP = 10*1000*1000;
	#endif
	boost::random::mt19937 rng;
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;
	Xform x,sum = Xform::Identity();
	rand_xform(rng,runif,rnorm,x);

	boost::timer::cpu_timer t;
	for(int i = 0; i < NSAMP; ++i){
		sum = sum * x;
	}
	double time = t.elapsed().wall;
	printf( "runtime %7.3fns nonsense: %7.3f \n", time/NSAMP, sum.translation()[0] );

}

TEST( xform_perf, preformance ){
	cout << "AffineCompact d "; test_xform_perf< Eigen::Transform<double,3,Eigen::AffineCompact> >();
	cout << "Affine        d "; test_xform_perf< Eigen::Affine3d  >();
	cout << "AffineCompact f "; test_xform_perf< Eigen::Transform<float ,3,Eigen::AffineCompact> >();
	cout << "Affine        f "; test_xform_perf< Eigen::Affine3f  >();
	// cout << "XformHash_bt24_Cubic_Zorder"; test_xform_perf< XformHash_bt24_Cubic_Zorder >();
}

}}}
