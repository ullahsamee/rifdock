#ifndef INCLUDED_numeric_rand_xform_HH
#define INCLUDED_numeric_rand_xform_HH


#include <Eigen/Geometry>

#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>

namespace scheme { namespace numeric { 


template<class T>
void
rand_xform(
	boost::random::mt19937 & rng,
	Eigen::Transform<T,3,Eigen::Affine> & x,
	T cart_bound = 512.0
){
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;
	Eigen::Quaterniond qrand( rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng) );
	qrand.normalize();
	Eigen::Matrix3d m = qrand.matrix();

	x.data()[3] = 0; x.data()[7] = 0; x.data()[11] = 0; x.data()[15] = 1;
	for(int i = 0; i <  3; ++i) x.data()[i] = m.data()[i-0];
	for(int i = 4; i <  7; ++i) x.data()[i] = m.data()[i-1];
	for(int i = 8; i < 11; ++i) x.data()[i] = m.data()[i-2];
	x.data()[12] = runif(rng) * cart_bound - cart_bound/2.0;
	x.data()[13] = runif(rng) * cart_bound - cart_bound/2.0;
	x.data()[14] = runif(rng) * cart_bound - cart_bound/2.0;
}

template<class T>
void
rand_xform(
	boost::random::mt19937 & rng,
	Eigen::Transform<T,3,Eigen::AffineCompact> & x,
	T cart_bound = 512.0
){
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;
	Eigen::Quaterniond qrand( rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng) );
	qrand.normalize();
	Eigen::Matrix3d m = qrand.matrix();
	for(int i = 0; i < 9; ++i) x.data()[i] = m.data()[i];

	x.data()[ 9] = runif(rng) * cart_bound - cart_bound/2.0;
	x.data()[10] = runif(rng) * cart_bound - cart_bound/2.0;
	x.data()[11] = runif(rng) * cart_bound - cart_bound/2.0;
}

template<class Float>
Float rad2quat(Float rad){
	return sqrt( 1.0 - cos( rad/2.0 )*cos( rad/2.0 ) );
}
template<class Float>
Float deg2quat(Float deg){
	Float rad = deg*M_PI/180.0;
	return rad2quat(rad);
}

template<class T>
void
rand_xform_quat(
	boost::random::mt19937 & rng,
	Eigen::Transform<T,3,Eigen::AffineCompact> & x,
	double cart_bound, double quat_bound
){
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;

	assert( quat_bound < sqrt(3.0)/2.0 );

	{ // ori part
		Eigen::Quaterniond qrand( 0.0, rnorm(rng), rnorm(rng), rnorm(rng) );
		double scale = 1.0 - runif(rng)*runif(rng);
		double len = qrand.norm();
		qrand.x() *= scale * quat_bound / len;
		qrand.y() *= scale * quat_bound / len;
		qrand.z() *= scale * quat_bound / len;
		qrand.w() = sqrt( 1.0 - qrand.squaredNorm() );
		assert( fabs( qrand.norm() - 1.0) < 0.00001 );

		qrand.normalize();
		Eigen::Matrix3d m = qrand.matrix();
		for(int i = 0; i < 9; ++i) x.data()[i] = m.data()[i];
	}
	
	
	{ // cart part
		x.data()[ 9] = rnorm(rng);
		x.data()[10] = rnorm(rng);
		x.data()[11] = rnorm(rng);
		double scale = 1.0 - runif(rng)*runif(rng);
		double len = sqrt( x.data()[9]*x.data()[9] + x.data()[10]*x.data()[10] + x.data()[11]*x.data()[11] );
		x.data()[ 9] *= scale * cart_bound / len;
		x.data()[10] *= scale * cart_bound / len;
		x.data()[11] *= scale * cart_bound / len;
	}
}


}}


#endif
