#include <gtest/gtest.h>

#include "scheme/objective/hash/XformMap.hh"
#include "scheme/numeric/rand_xform.hh"
#include <Eigen/Geometry>

#include <google/dense_hash_set>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>
#include <google/dense_hash_set>

namespace scheme { namespace objective { namespace hash { namespace xhtest {

using std::cout;
using std::endl;

typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;


template< template<class X> class XformHash >
int get_num_ori_cells( int ori_nside, double & xcov ){
	boost::random::mt19937 rng((unsigned int)time(0) + 7693487);
	XformHash<Xform> xh( 1.0, ori_nside, 512.0 );
		int n_ori_bins;
	{
		google::dense_hash_set<size_t> idx_seen;
		idx_seen.set_empty_key(std::numeric_limits<uint64_t>::max());

		int NSAMP = 5000*ori_nside*ori_nside;
		Xform x;
		for(int i = 0; i < NSAMP; ++i){
			numeric::rand_xform(rng,x);
			x.translation()[0] = x.translation()[1] = x.translation()[2] = 0;						
			idx_seen.insert( xh.get_key(x) );
		}
		n_ori_bins = (int)idx_seen.size();
		xcov = (double)NSAMP / n_ori_bins;
	}
	return n_ori_bins;
}

TEST( XformHash_Quat_BCC7_Zorder, num_ori_cells ){

	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 2,512.0); ASSERT_EQ( xh.approx_nori() ,     53 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 3,512.0); ASSERT_EQ( xh.approx_nori() ,    134 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 4,512.0); ASSERT_EQ( xh.approx_nori() ,    189 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 5,512.0); ASSERT_EQ( xh.approx_nori() ,    436 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 6,512.0); ASSERT_EQ( xh.approx_nori() ,    622 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 7,512.0); ASSERT_EQ( xh.approx_nori() ,    899 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 8,512.0); ASSERT_EQ( xh.approx_nori() ,   1606 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 9,512.0); ASSERT_EQ( xh.approx_nori() ,   1996 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,10,512.0); ASSERT_EQ( xh.approx_nori() ,   2303 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,11,512.0); ASSERT_EQ( xh.approx_nori() ,   3410 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,12,512.0); ASSERT_EQ( xh.approx_nori() ,   4502 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,13,512.0); ASSERT_EQ( xh.approx_nori() ,   5510 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,14,512.0); ASSERT_EQ( xh.approx_nori() ,   6284 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,15,512.0); ASSERT_EQ( xh.approx_nori() ,   8285 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,16,512.0); ASSERT_EQ( xh.approx_nori() ,  10098 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,17,512.0); ASSERT_EQ( xh.approx_nori() ,  11634 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,18,512.0); ASSERT_EQ( xh.approx_nori() ,  13352 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,19,512.0); ASSERT_EQ( xh.approx_nori() ,  16065 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,20,512.0); ASSERT_EQ( xh.approx_nori() ,  18538 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,21,512.0); ASSERT_EQ( xh.approx_nori() ,  21205 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,22,512.0); ASSERT_EQ( xh.approx_nori() ,  23953 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,23,512.0); ASSERT_EQ( xh.approx_nori() ,  28212 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,24,512.0); ASSERT_EQ( xh.approx_nori() ,  31593 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,25,512.0); ASSERT_EQ( xh.approx_nori() ,  35653 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,26,512.0); ASSERT_EQ( xh.approx_nori() ,  38748 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,27,512.0); ASSERT_EQ( xh.approx_nori() ,  43980 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,28,512.0); ASSERT_EQ( xh.approx_nori() ,  48801 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,29,512.0); ASSERT_EQ( xh.approx_nori() ,  54661 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,30,512.0); ASSERT_EQ( xh.approx_nori() ,  58271 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,31,512.0); ASSERT_EQ( xh.approx_nori() ,  65655 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,32,512.0); ASSERT_EQ( xh.approx_nori() ,  72114 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,33,512.0); ASSERT_EQ( xh.approx_nori() ,  79038 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,34,512.0); ASSERT_EQ( xh.approx_nori() ,  84326 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,35,512.0); ASSERT_EQ( xh.approx_nori() ,  93094 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,36,512.0); ASSERT_EQ( xh.approx_nori() , 101191 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,37,512.0); ASSERT_EQ( xh.approx_nori() , 109680 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,38,512.0); ASSERT_EQ( xh.approx_nori() , 116143 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,39,512.0); ASSERT_EQ( xh.approx_nori() , 127688 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,40,512.0); ASSERT_EQ( xh.approx_nori() , 137387 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,41,512.0); ASSERT_EQ( xh.approx_nori() , 146325 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,42,512.0); ASSERT_EQ( xh.approx_nori() , 155608 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,43,512.0); ASSERT_EQ( xh.approx_nori() , 168954 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,44,512.0); ASSERT_EQ( xh.approx_nori() , 180147 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,45,512.0); ASSERT_EQ( xh.approx_nori() , 192798 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,46,512.0); ASSERT_EQ( xh.approx_nori() , 202438 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,47,512.0); ASSERT_EQ( xh.approx_nori() , 218861 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,48,512.0); ASSERT_EQ( xh.approx_nori() , 231649 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,49,512.0); ASSERT_EQ( xh.approx_nori() , 246830 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,50,512.0); ASSERT_EQ( xh.approx_nori() , 257380 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,51,512.0); ASSERT_EQ( xh.approx_nori() , 275655 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,52,512.0); ASSERT_EQ( xh.approx_nori() , 292355 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,53,512.0); ASSERT_EQ( xh.approx_nori() , 309321 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,54,512.0); ASSERT_EQ( xh.approx_nori() , 321798 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,55,512.0); ASSERT_EQ( xh.approx_nori() , 343505 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,56,512.0); ASSERT_EQ( xh.approx_nori() , 362585 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,57,512.0); ASSERT_EQ( xh.approx_nori() , 381254 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,58,512.0); ASSERT_EQ( xh.approx_nori() , 396135 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,59,512.0); ASSERT_EQ( xh.approx_nori() , 420820 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,60,512.0); ASSERT_EQ( xh.approx_nori() , 442324 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,61,512.0); ASSERT_EQ( xh.approx_nori() , 464576 ); }
	{ XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,62,512.0); ASSERT_EQ( xh.approx_nori() , 480460 ); }

	// for(int i = 2; i < 65; ++i){
	// 	double xcov;
	// 	int num_ori_cells = get_num_ori_cells<XformHash_Quat_BCC7_Zorder>( i, xcov );
	// 	printf("XformHash_Quat_BCC7_Zorder %2d %7d %8.3f\n",i, num_ori_cells, xcov );
	// 	std::cout.flush();
	// }
}

TEST( XformHash_bt24_BCC6, DISABLED_________________num_ori_cells ){
	for(int i = 19; i < 65; ++i){
		double xcov;
		int num_ori_cells = get_num_ori_cells<XformHash_bt24_BCC6>( i, xcov );
		printf("XformHash_bt24_BCC6 %2d %7d %8.3f\n",i, num_ori_cells, xcov );
		std::cout.flush();
	}
}


template< template<class X> class XformHash >
void test_xform_hash_perf( double cart_resl, double ang_resl, int const N2=100*1000, unsigned int seed = 0){
	boost::random::mt19937 rng((unsigned int)time(0) + seed);

	double time_key = 0.0, time_cen = 0.0;
	double cart_resl2 = cart_resl*cart_resl;
	double  ang_resl2 = ang_resl* ang_resl;

	XformHash<Xform> xh( cart_resl, ang_resl, 512.0 );

	std::vector<Xform> samples(N2),centers(N2);
	
	for(int i = 0; i < N2; ++i)
		numeric::rand_xform(rng,samples[i]);

	boost::timer::cpu_timer tk;
	std::vector<uint64_t> keys(N2);
	for(int i = 0; i < N2; ++i){
		keys[i] = xh.get_key( samples[i] );
		// centers[i] = xh.get_center( keys[i] );
		// cout << endl;
	}
	time_key += (double)tk.elapsed().wall;


	boost::timer::cpu_timer tc;
	for(int i = 0; i < N2; ++i)
		centers[i] = xh.get_center( keys[i] );
	time_cen += (double)tc.elapsed().wall;

	google::dense_hash_set<size_t> idx_seen;
	idx_seen.set_empty_key(std::numeric_limits<uint64_t>::max());
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

	double tot_cell_vol = covrad*covrad*covrad*covrad*covrad*covrad * xh.approx_nori() / (cart_resl*cart_resl*cart_resl);
	printf( " %5.3f/%5.1f cr %5.3f dt %5.3f da %6.3f x2k: %7.3fns k2x: %7.3fns %9.3f %7llu\n", 
		cart_resl, ang_resl, covrad, max_dt, max_da/ang_resl, time_key/N2, time_cen/N2, tot_cell_vol, xh.approx_nori() );

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

// TEST( XformHash, XformHash_bt24_BCC3 ){
// 	unsigned int s = 0;
// 	int N = 10*1000;
// 	#ifdef NDEBUG
// 	N = 1*1000*1000;
// 	#endif
// 	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 4.00 , 30.0, N, ++s );
// 	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 2.00 , 20.0, N, ++s );
// 	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 1.00 , 15.0, N, ++s );
// 	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 0.50 , 10.0, N, ++s );
// 	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 0.25 ,  5.0, N, ++s );
// 	cout << "XformHash_bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3 >( 0.11 ,  3.3, N, ++s );
// }

// TEST( XformHash, preformance ){

// 	cout << "XformHash_Quat_BCC7        "; test_xform_hash_perf< XformHash_Quat_BCC7         >( 2.0 , 1 );
// 	cout << "XformHash_Quat_BCC7_Zorder "; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder  >( 2.0 , 2 );
// 	cout << "XformHash_bt24_BCC3        "; test_xform_hash_perf< XformHash_bt24_BCC3         >( 2.0 , 3 );
// 	cout << "XformHash_bt24_BCC3_Zorder "; test_xform_hash_perf< XformHash_bt24_BCC3_Zorder  >( 2.0 , 4 );
// 	cout << "XformHash_bt24_Cubic_Zorder"; test_xform_hash_perf< XformHash_bt24_Cubic_Zorder >( 2.0 , 5 );

// 	cout << "XformHash_Quat_BCC7        "; test_xform_hash_perf< XformHash_Quat_BCC7         >( 0.25 , 1 );
// 	cout << "XformHash_Quat_BCC7_Zorder "; test_xform_hash_perf< XformHash_Quat_BCC7_Zorder  >( 0.25 , 2 );
// 	cout << "XformHash_bt24_BCC3        "; test_xform_hash_perf< XformHash_bt24_BCC3         >( 0.25 , 3 );
// 	cout << "XformHash_bt24_BCC3_Zorder "; test_xform_hash_perf< XformHash_bt24_BCC3_Zorder  >( 0.25 , 4 );
// 	cout << "XformHash_bt24_Cubic_Zorder"; test_xform_hash_perf< XformHash_bt24_Cubic_Zorder >( 0.25 , 5 );
// }


TEST( XformHash, XformHash_Quat_BCC7_Zorder_cart_shift ){
	int NSAMP = 20;
	#ifdef NDEBUG
	NSAMP = 50;
	#endif
	using namespace Eigen;
	XformHash_Quat_BCC7_Zorder<Xform> h( 1.0, 10.0 );
	Xform m = Xform::Identity();
	uint64_t mkey = h.get_key(m);
	ASSERT_EQ( m.translation(), Vector3d(0,0,0) );
	Vector3d v; 
	ASSERT_TRUE( mkey & 1 ); // middle key is odd
	for(int x = -NSAMP; x <= NSAMP; ++x){
	for(int y = -NSAMP; y <= NSAMP; ++y){
	for(int z = -NSAMP; z <= NSAMP; ++z){
		// Odd center cell
		v = h.get_center( h.cart_shift_key( mkey, x, y, z ) ).translation()/h.cart_spacing();
		ASSERT_LE( (Vector3d(x,y,z)-v).squaredNorm(), 0.00001 );
		// Even center cell
		v = h.get_center( h.cart_shift_key( mkey&~1, x, y, z ) ).translation()/h.cart_spacing();
		ASSERT_LE( (Vector3d(x-0.5,y-0.5,z-0.5)-v).squaredNorm(), 0.00001 );
	}}}
}



}}}}
