#ifndef INCLUDED_numeric_rand_xform_HH
#define INCLUDED_numeric_rand_xform_HH

#include "scheme/numeric/util.hh"
#include <Eigen/Geometry>
#include <random>

namespace scheme { namespace numeric {

// This returns a random, un-normalized vector from a spherical distribution
template<class RNG, class T>
void
rand_vector_sphere(
	RNG & rng,
	T & vec
){
	std::normal_distribution<> rnorm;

	vec[0] = rnorm(rng);
	vec[1] = rnorm(rng);
	vec[2] = rnorm(rng);
}


template<class RNG, class T>
void
rand_xform(
	RNG & rng,
	Eigen::Transform<T,3,Eigen::Affine> & x,
	T cart_bound = 512.0
){
	std::uniform_real_distribution<> runif;
	std::normal_distribution<> rnorm;
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

template<class RNG, class T>
void
rand_xform(
	RNG & rng,
	Eigen::Transform<T,3,Eigen::AffineCompact> & x,
	T cart_bound = 512.0
){
	std::uniform_real_distribution<> runif;
	std::normal_distribution<> rnorm;
	Eigen::Quaterniond qrand( rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng) );
	qrand.normalize();
	Eigen::Matrix3d m = qrand.matrix();
	for(int i = 0; i < 9; ++i) x.data()[i] = m.data()[i];

	x.data()[ 9] = runif(rng) * cart_bound - cart_bound/2.0;
	x.data()[10] = runif(rng) * cart_bound - cart_bound/2.0;
	x.data()[11] = runif(rng) * cart_bound - cart_bound/2.0;
}


template<class RNG, class X>
X rand_xform(
	RNG & rng,
	scalar<X> cart_bound = 512.0
){
	X x;
	rand_xform(rng,x,cart_bound);
	return x;
}

template<class RNG, class T>
void
rand_xform_cartnormal(
	RNG & rng,
	Eigen::Transform<T,3,Eigen::AffineCompact> & x,
	T const & cart_sd
){
	std::uniform_real_distribution<> runif;
	std::normal_distribution<> rnorm;
	Eigen::Quaterniond qrand( rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng) );
	qrand.normalize();
	Eigen::Matrix3d m = qrand.matrix();
	for(int i = 0; i < 9; ++i) x.data()[i] = m.data()[i];
	x.data()[ 9] = rnorm(rng) * cart_sd;
	x.data()[10] = rnorm(rng) * cart_sd;
	x.data()[11] = rnorm(rng) * cart_sd;
}

template<class RNG, class T>
void
rand_xform_quat(
	RNG & rng,
	Eigen::Transform<T,3,Eigen::AffineCompact> & x,
	double cart_bound, double quat_bound
){
	std::uniform_real_distribution<> runif;
	std::normal_distribution<> rnorm;

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

template<class RNG, class T>
void
rand_xform_sphere(
	RNG & rng,
	Eigen::Transform<T,3,Eigen::AffineCompact> & x,
	T const cart_radius,
	T const ang_radius
){
	std::normal_distribution<> rnorm;
	std::uniform_real_distribution<> runif;

	float ang = (1.0 - runif(rng)*runif(rng)) * ang_radius;
	Eigen::Matrix<float,1,3> axis( rnorm(rng), rnorm(rng), rnorm(rng) );
	axis.normalize();
	Eigen::AngleAxis<T> aa( ang, axis );
	x = Eigen::Transform<T,3,Eigen::AffineCompact> ( aa );

	// Eigen::AngleAxis<T> aa2( x.rotation() );

	// std::cout << ang << " " << ang_radius << " " << aa2.angle() << std::endl;

	Eigen::Matrix<float,1,3> delta( 9e9, 9e9, 9e9 );
	while( delta.squaredNorm() > cart_radius*cart_radius ){
		delta = 2.0 * cart_radius * Eigen::Matrix<float,1,3>( runif(rng)-0.5, runif(rng)-0.5, runif(rng)-0.5 );
	}

	x.translation() = delta;

}



}}


#endif
