#include <gtest/gtest.h>

#include "scheme/objective/hash/XformMap.hh"
#include <Eigen/Geometry>

#include <google/dense_hash_set>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>

namespace scheme { namespace objective { namespace hash { namespace test {

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
	x.data()[12] = runif(rng)*500.0-250.0;
	x.data()[13] = runif(rng)*500.0-250.0;
	x.data()[14] = runif(rng)*500.0-250.0;
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
	x.data()[ 9] = runif(rng)*500.0-250.0;
	x.data()[10] = runif(rng)*500.0-250.0;
	x.data()[11] = runif(rng)*500.0-250.0;
}

// TEST( XformMap, basic_test ){
template<
	template<class X> class XformHash
>
void test_xform_hash_perf( double cart_resl, double ang_resl, int const N2=100*1000, unsigned int seed = 0){
	boost::random::mt19937 rng((unsigned int)time(0) + seed);
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;

		double time_key = 0.0, time_cen = 0.0;
		double cart_resl2 = cart_resl*cart_resl;
		double  ang_resl2 = ang_resl* ang_resl;

		XformHash<Xform> xm( cart_resl, ang_resl, 512.0 );

		std::vector<Xform> samples(N2),centers(N2);
		
		for(int i = 0; i < N2; ++i)
			rand_xform(rng,runif,rnorm,samples[i]);

		boost::timer::cpu_timer tk;
		std::vector<uint64_t> keys(N2);
		for(int i = 0; i < N2; ++i){
			keys[i] = xm.get_key( samples[i] );
			// centers[i] = xm.get_center( keys[i] );
			// cout << endl;
		}
		time_key += (double)tk.elapsed().wall;


		boost::timer::cpu_timer tc;
		for(int i = 0; i < N2; ++i)
			centers[i] = xm.get_center( keys[i] );
		time_cen += (double)tc.elapsed().wall;

		google::dense_hash_set<size_t> idx_seen;
		idx_seen.set_empty_key(999999999999);
		for(int i = 0; i < N2; ++i) idx_seen.insert(keys[i]);

		double covrad=0, max_dt=0, max_da=0;
		for(int i = 0; i < N2; ++i){
			Xform l = centers[i].inverse() * samples[i];
			double dt = l.translation().norm();
			Eigen::Matrix3d m; for(int k=0; k<9;++k) m.data()[k] = l.data()[k];
			// cout << m << endl;
			// cout << l.rotation() << endl;
			double da = Eigen::AngleAxisd(m).angle()*180.0/M_PI;
			// double da = Eigen::AngleAxisd(l.rotation()).angle()*180.0/M_PI;
			double err = sqrt( da*da/ang_resl2*cart_resl2 + dt*dt );
			covrad = fmax(covrad,err);
			max_dt = fmax(max_dt,dt);
			max_da = fmax(max_da,da);
			ASSERT_LT( dt, cart_resl );
			ASSERT_LT( da,  ang_resl );
		}

		// double tot_cell_vol = covrad*covrad*covrad*covrad*covrad*covrad * xm.approx_size();
		printf( " %5.3f/%5.1f cr %5.3f dt %5.3f da %6.3f x2k: %7.3fns k2x: %7.3fns \n", 
			cart_resl, ang_resl, covrad/cart_resl, max_dt/cart_resl, max_da/ang_resl, time_key/N2, time_cen/N2 );

	// cout << " rate " << N1*N2/time_key << "  " << N1*N2/time_cen << endl;
}

TEST( XformHash, XformHash_Quat_BCC7_Zorder ){
	unsigned int s = 0;
	int N = 10*1000;
	#ifdef NDEBUG
	N = 1*1000*1000;
	#endif
	cout << "XformHash_Quat_BCC7_Zorder"; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder >( 4.00 , 30.0, N, ++s );
	cout << "XformHash_Quat_BCC7_Zorder"; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder >( 2.00 , 20.0, N, ++s );
	cout << "XformHash_Quat_BCC7_Zorder"; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder >( 1.00 , 15.0, N, ++s );
	cout << "XformHash_Quat_BCC7_Zorder"; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder >( 0.50 , 10.0, N, ++s );
	cout << "XformHash_Quat_BCC7_Zorder"; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder >( 0.25 ,  5.0, N, ++s );
	cout << "XformHash_Quat_BCC7_Zorder"; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder >( 0.11 ,  3.3, N, ++s );
}

TEST( XformHash, DISABLED_XformHash_bt24_BCC3 ){
	unsigned int s = 0;
	int N = 10*1000;
	#ifdef NDEBUG
	N = 1*1000*1000;
	#endif
	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 4.00 , 30.0, N, ++s );
	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 2.00 , 20.0, N, ++s );
	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 1.00 , 15.0, N, ++s );
	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 0.50 , 10.0, N, ++s );
	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 0.25 ,  5.0, N, ++s );
	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 0.11 ,  3.3, N, ++s );
}

TEST( XformHash, XformHash_bt24_BCC6 ){
	unsigned int s = 0;
	int N = 10*1000;
	#ifdef NDEBUG
	N = 1*1000*1000;
	#endif
	cout << "XformHash_bt24_BCC6"; test_xform_hash_perf< XformHash_bt24_BCC6 >( 4.00 , 30.0, N, ++s );
	cout << "XformHash_bt24_BCC6"; test_xform_hash_perf< XformHash_bt24_BCC6 >( 2.00 , 20.0, N, ++s );
	cout << "XformHash_bt24_BCC6"; test_xform_hash_perf< XformHash_bt24_BCC6 >( 1.00 , 15.0, N, ++s );
	cout << "XformHash_bt24_BCC6"; test_xform_hash_perf< XformHash_bt24_BCC6 >( 0.50 , 10.0, N, ++s );
	cout << "XformHash_bt24_BCC6"; test_xform_hash_perf< XformHash_bt24_BCC6 >( 0.25 ,  5.0, N, ++s );
	cout << "XformHash_bt24_BCC6"; test_xform_hash_perf< XformHash_bt24_BCC6 >( 0.11 ,  3.3, N, ++s );
}

TEST( XformHash, DISABLED_preformance ){

	cout << "XformHash_Quat_BCC7        "; test_xform_hash_perf< XformHash_Quat_BCC7         >( 2.0 , 1 );
	cout << "XformHash_Quat_BCC7_Zorder "; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder  >( 2.0 , 2 );
	cout << "XformHash_bt24_BCC3        "; test_xform_hash_perf< XformHash_bt24_BCC3         >( 2.0 , 3 );
	cout << "XformHash_bt24_BCC3_Zorder "; test_xform_hash_perf< XformHash_bt24_BCC3_Zorder  >( 2.0 , 4 );
	cout << "XformHash_bt24_Cubic_Zorder"; test_xform_hash_perf< XformHash_bt24_Cubic_Zorder >( 2.0 , 5 );

	cout << "XformHash_Quat_BCC7        "; test_xform_hash_perf< XformHash_Quat_BCC7         >( 0.25 , 1 );
	cout << "XformHash_Quat_BCC7_Zorder "; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder  >( 0.25 , 2 );
	cout << "XformHash_bt24_BCC3        "; test_xform_hash_perf< XformHash_bt24_BCC3         >( 0.25 , 3 );
	cout << "XformHash_bt24_BCC3_Zorder "; test_xform_hash_perf< XformHash_bt24_BCC3_Zorder  >( 0.25 , 4 );
	cout << "XformHash_bt24_Cubic_Zorder"; test_xform_hash_perf< XformHash_bt24_Cubic_Zorder >( 0.25 , 5 );
}


TEST( XformHash, XformHash_Quat_BCC7_Zorder_Neighbor_Lookup ){

}

}}}}
