#include <gtest/gtest.h>

#include "scheme/objective/hash/XformHash.hh"
#include "scheme/objective/hash/XformHashNeighbors.hh"
#include "scheme/numeric/rand_xform.hh"

#include <Eigen/Geometry>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>

namespace scheme { namespace objective { namespace hash { namespace xhnbtest {

using std::cout;
using std::endl;

typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;

template<class XH>
void get_neighbors_ref(
	XH const & xh,
	Xform const & x,
	double cart_bound, double ang_bound,
	std::set<typename XH::Key> & nbrs,
	int NSAMP=100000
){
	std::set<typename XH::Key> keys;
	double quat_bound = numeric::deg2quat(ang_bound);
	boost::random::mt19937 rng((unsigned int)time(0) + 9384756);
	Xform xrand;
	for(int i = 0; i < NSAMP; ++i){
		numeric::rand_xform(rng,xrand,cart_bound,quat_bound);
		nbrs.insert( xh.get_key( xrand*x ) );
	}
}


TEST( XformHashNeighbors, DISABLED_Quat_BCC7_Zorder_check_rotations ){
	int NSAMP = 10;
	#ifdef NDEBUG
	NSAMP = 100;
	#endif
	boost::random::mt19937 rng((unsigned int)time(0) + 94586049);

	XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,10.0);
	cout << xh.approx_nori()/98.0 * 3.1 << endl;
	XformHashNeighbors< XformHash_Quat_BCC7_Zorder<Xform> > nb(3.0,30.0,xh);

	for(int i = 0; i < NSAMP; ++i){
		Xform x; numeric::rand_xform(rng,x);
		nb.get_ori_neighbors( xh.get_key(x) );
	}

}

TEST( XformHashNeighbors, Quat_BCC7_Zorder_can_repr_identity ){
	boost::random::mt19937 rng((unsigned int)time(0) + 94586049);
	boost::uniform_real<> runif;
	using namespace Eigen;
	for(int j = 0; j < 10000; ++j){
		for(int i = 2; i < 10; ++i){
			XformHash_Quat_BCC7_Zorder<Xform> xh( 0.1+runif(rng), i, runif(rng)*400.0+2.0 );
			util::SimpleArray<7,double> f7 = xh.grid_[xh.grid_[util::SimpleArray<7,double>(0.0)]];
			for(int k = 0; k < 7; ++k) ASSERT_NEAR( f7[k], 0.0, 0.000001 );
			Xform x = Xform::Identity();
			// x.translation()[0] = 2.0*runif(rng)-1.0; // this doesn't work... even/odd cells
			// x.translation()[1] = 2.0*runif(rng)-1.0; //
			// x.translation()[2] = 2.0*runif(rng)-1.0; //
			Quaterniond q( xh.get_center( xh.get_key( x ) ).rotation() );
			ASSERT_NEAR( q.w(), 1.0, 0.000001 );
			ASSERT_NEAR( q.x(), 0.0, 0.000001 );
			ASSERT_NEAR( q.y(), 0.0, 0.000001 );
			ASSERT_NEAR( q.z(), 0.0, 0.000001 );
		}
	}
}

TEST( XformHashNeighbors, Quat_BCC7_Zorder_quaternion_refllection_symmetry ){
	using namespace Eigen;
	typedef util::SimpleArray<7,double> F7;
	typedef util::SimpleArray<7,uint64_t> I7;	
	boost::random::mt19937 rng((unsigned int)time(0) + 293457);
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;

	for(int iter1 = 0; iter1 < 100; ++iter1){
		XformHash_Quat_BCC7_Zorder<Xform> xh( runif(rng)+0.1, 2+(int)(runif(rng)*20.0), runif(rng)*400+10.0 );
		for(int iter2 = 0; iter2 < 100; ++iter2){
			double a = runif(rng);
			double b = runif(rng);
			double c = runif(rng);		
			double x = rnorm(rng);
			double y = rnorm(rng);
			double z = rnorm(rng);
			double w = rnorm(rng);
			double l = sqrt(w*w+x*x+y*y+z*z);
			x /= l;
			w /= l;
			y /= l;
			z /= l;
			assert( fabs(w*w+x*x+y*y+z*z - 1.0) < 0.00001 );

			std::vector<uint64_t> i(16);
			i[ 0] = xh.grid_[ F7( a,b,c, +w,+x,+y,+z ) ];
			i[ 1] = xh.grid_[ F7( a,b,c, +w,+x,+y,-z ) ];
			i[ 2] = xh.grid_[ F7( a,b,c, +w,+x,-y,+z ) ];
			i[ 3] = xh.grid_[ F7( a,b,c, +w,+x,-y,-z ) ];
			i[ 4] = xh.grid_[ F7( a,b,c, +w,-x,+y,+z ) ];
			i[ 5] = xh.grid_[ F7( a,b,c, +w,-x,+y,-z ) ];
			i[ 6] = xh.grid_[ F7( a,b,c, +w,-x,-y,+z ) ];
			i[ 7] = xh.grid_[ F7( a,b,c, +w,-x,-y,-z ) ];	
			i[ 8] = xh.grid_[ F7( a,b,c, -w,+x,+y,+z ) ];
			i[ 9] = xh.grid_[ F7( a,b,c, -w,+x,+y,-z ) ];
			i[10] = xh.grid_[ F7( a,b,c, -w,+x,-y,+z ) ];
			i[11] = xh.grid_[ F7( a,b,c, -w,+x,-y,-z ) ];
			i[12] = xh.grid_[ F7( a,b,c, -w,-x,+y,+z ) ];
			i[13] = xh.grid_[ F7( a,b,c, -w,-x,+y,-z ) ];
			i[14] = xh.grid_[ F7( a,b,c, -w,-x,-y,+z ) ];
			i[15] = xh.grid_[ F7( a,b,c, -w,-x,-y,-z ) ];	


			// cout << i[0] <<" "<< i[1] <<" "<< i[2] <<" "<< i[3] <<" "<< i[4] <<" "<< i[5] <<" "<< i[6] <<" "<< i[7] << endl;	
			std::vector<F7> f(16);
			f[ 0] = xh.grid_[ xh.grid_[ F7( a,b,c, +w,+x,+y,+z ) ] ];
			f[ 1] = xh.grid_[ xh.grid_[ F7( a,b,c, +w,+x,+y,-z ) ] ];
			f[ 2] = xh.grid_[ xh.grid_[ F7( a,b,c, +w,+x,-y,+z ) ] ];
			f[ 3] = xh.grid_[ xh.grid_[ F7( a,b,c, +w,+x,-y,-z ) ] ];
			f[ 4] = xh.grid_[ xh.grid_[ F7( a,b,c, +w,-x,+y,+z ) ] ];
			f[ 5] = xh.grid_[ xh.grid_[ F7( a,b,c, +w,-x,+y,-z ) ] ];
			f[ 6] = xh.grid_[ xh.grid_[ F7( a,b,c, +w,-x,-y,+z ) ] ];
			f[ 7] = xh.grid_[ xh.grid_[ F7( a,b,c, +w,-x,-y,-z ) ] ];
			f[ 8] = xh.grid_[ xh.grid_[ F7( a,b,c, -w,+x,+y,+z ) ] ];
			f[ 9] = xh.grid_[ xh.grid_[ F7( a,b,c, -w,+x,+y,-z ) ] ];
			f[10] = xh.grid_[ xh.grid_[ F7( a,b,c, -w,+x,-y,+z ) ] ];
			f[11] = xh.grid_[ xh.grid_[ F7( a,b,c, -w,+x,-y,-z ) ] ];
			f[12] = xh.grid_[ xh.grid_[ F7( a,b,c, -w,-x,+y,+z ) ] ];
			f[13] = xh.grid_[ xh.grid_[ F7( a,b,c, -w,-x,+y,-z ) ] ];
			f[14] = xh.grid_[ xh.grid_[ F7( a,b,c, -w,-x,-y,+z ) ] ];
			f[15] = xh.grid_[ xh.grid_[ F7( a,b,c, -w,-x,-y,-z ) ] ];
			for(int j = 1; j < 16; ++j){
				ASSERT_EQ  (      i[0]%2  ,      i[j]%2    ); // all even or all odd
				ASSERT_NEAR( fabs(f[0][3]), fabs(f[j][3]) , 0.000001 );
				ASSERT_NEAR( fabs(f[0][4]), fabs(f[j][4]) , 0.000001 );
				ASSERT_NEAR( fabs(f[0][5]), fabs(f[j][5]) , 0.000001 );
				ASSERT_NEAR( fabs(f[0][6]), fabs(f[j][6]) , 0.000001 );
			}

			if( iter1+iter2 == 0 ){
				for(int j=0; j<16; ++j){
					Xform xf = Xform( nest::maps::to_half_cell( Quaterniond(f[j][3],f[j][4],f[j][5],f[j][6]) ).matrix());
					for(int k=0;k<3;++k) xf.translation()[k] = f[j][k];
					uint64_t k = xh.get_key( xf );
					uint64_t o = k & 1;
					uint64_t w =  util::undilate<7>( k>>4 ) & 63;
					uint64_t x =  util::undilate<7>( k>>5 ) & 63;
					uint64_t y =  util::undilate<7>( k>>6 ) & 63;
					uint64_t z =  util::undilate<7>( k>>7 ) & 63;
					std::cout << xh.grid_.nside_[3]-o <<" " << o << "    " << w << "\t" << x << "\t" << y << "\t" << z << std::endl;
				}
			}

		}
	}
	
}

TEST( XformHashNeighbors, Quat_BCC7_Zorder_key_symmetry ){
	using namespace Eigen;
	typedef util::SimpleArray<7,double> F7;
	typedef util::SimpleArray<7,uint64_t> I7;	
	boost::random::mt19937 rng((unsigned int)time(0) + 293457);
	boost::uniform_real<> runif;

	XformHash_Quat_BCC7_Zorder<Xform> xh( runif(rng)+0.1, 2+(int)(runif(rng)*20.0), runif(rng)*400+10.0 );

	for(int j = 0; j < 100; ++j){
		XformHash_Quat_BCC7_Zorder<Xform> xh( runif(rng)+0.1, 2+(int)(runif(rng)*20.0), runif(rng)*400+10.0 );
		for(int i = 0; i < 1000; ++i){
			Xform x; numeric::rand_xform(rng,x);
			uint64_t key = xh.get_key(x);
			uint64_t pkey,isym;
			pkey = xh.asym_key(key,isym);

				// {
				// 	uint64_t k = key;
				// 	uint64_t o = k & 1;
				// 	uint64_t w =  util::undilate<7>( k>>4 ) & 63;
				// 	uint64_t x =  util::undilate<7>( k>>5 ) & 63;
				// 	uint64_t y =  util::undilate<7>( k>>6 ) & 63;
				// 	uint64_t z =  util::undilate<7>( k>>7 ) & 63;
				// 	std::cout << "KEY  " << isym << " " << xh.grid_.nside_[3]-o <<" " << o << "    " << w << "\t" << x << "\t" << y << "\t" << z << std::endl;
				// }
				// {
				// 	uint64_t k = pkey;
				// 	uint64_t o = k & 1;
				// 	uint64_t w =  util::undilate<7>( k>>4 ) & 63;
				// 	uint64_t x =  util::undilate<7>( k>>5 ) & 63;
				// 	uint64_t y =  util::undilate<7>( k>>6 ) & 63;
				// 	uint64_t z =  util::undilate<7>( k>>7 ) & 63;
				// 	std::cout << "ASYM " << isym << " " << xh.grid_.nside_[3]-o <<" " << o << "    " << w << "\t" << x << "\t" << y << "\t" << z << std::endl;
				// }

			ASSERT_LT( isym, xh.num_key_symmetries() );
			ASSERT_EQ( key, xh.sym_key(pkey,isym) );
			if( isym==0 ) ASSERT_EQ( key, pkey );
			else          ASSERT_NE( key, pkey );

			// TODO: check geometry as in test above

		}
	}
}

}}}}



























