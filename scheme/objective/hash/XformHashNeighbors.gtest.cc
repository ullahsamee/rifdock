#include <gtest/gtest.h>

#include "scheme/objective/hash/XformHash.hh"
#include "scheme/objective/hash/XformHashNeighbors.hh"
#include "scheme/numeric/rand_xform.hh"

#include <Eigen/Geometry>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>
#include <sparsehash/dense_hash_set>

namespace scheme { namespace objective { namespace hash { namespace xhnbtest {

using std::cout;
using std::endl;

typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;

TEST( FOO, BAR ){
	ASSERT_EQ( -0.0, 0.0 );
	ASSERT_EQ( nest::maps::to_half_cell(Eigen::Quaterniond( -1, 1, 1, 1 ) ).w(), 1.0 );
	ASSERT_EQ( nest::maps::to_half_cell(Eigen::Quaterniond(  0,-1, 1, 1 ) ).x(), 1.0 );
	ASSERT_EQ( nest::maps::to_half_cell(Eigen::Quaterniond(  0, 0,-1, 1 ) ).y(), 1.0 );
	ASSERT_EQ( nest::maps::to_half_cell(Eigen::Quaterniond(  0, 0, 0,-1 ) ).z(), 1.0 );
}

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
		numeric::rand_xform_quat(rng,xrand,cart_bound,quat_bound);
		nbrs.insert( xh.get_key( xrand*x ) );
	}
}


TEST( XformHashNeighbors, Quat_BCC7_Zorder_can_repr_identity ){
	boost::random::mt19937 rng((unsigned int)time(0) + 94586049);
	boost::uniform_real<> runif;
	using namespace Eigen;
	for(int j = 0; j < 10000; ++j){
		for(int i = 2; i < 10; ++i){
			XformHash_Quat_BCC7_Zorder<Xform> xh( 0.1+runif(rng), i, runif(rng)*400.0+2.0 );
			uint64_t k0 = xh.grid_[util::SimpleArray<7,double>(0.0)];
			ASSERT_EQ( k0%2, (i+1)%2 );
			util::SimpleArray<7,double> f7 = xh.grid_[k0];
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
		XformHash_Quat_BCC7_Zorder<Xform> xh( runif(rng)+0.1, 2+(int)(runif(rng)*60.0), runif(rng)*400+10.0 );
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

			// if( iter1+iter2 == 0 ){
			// 	for(int j=0; j<8; ++j){
			// 		Xform xf = Xform( nest::maps::to_half_cell( Quaterniond(f[j][3],f[j][4],f[j][5],f[j][6]) ).matrix());
			// 		for(int k=0;k<3;++k) xf.translation()[k] = f[j][k];
			// 		uint64_t k = xh.get_key( xf );
			// 		uint64_t o = k & 1;
			// 		uint64_t w =  util::undilate<7>( k>>4 ) & 63;
			// 		uint64_t x =  util::undilate<7>( k>>5 ) & 63;
			// 		uint64_t y =  util::undilate<7>( k>>6 ) & 63;
			// 		uint64_t z =  util::undilate<7>( k>>7 ) & 63;
			// 		std::cout << xh.grid_.nside_[3]-o <<" " << o << "    " << w << "\t" << x << "\t" << y << "\t" << z << std::endl;
			// 	}
			// }

		}
	}
	
}

TEST( XformHashNeighbors, Quat_BCC7_Zorder_key_symmetry ){
	int NSAMP = 100;
	#ifdef NDEBUG
	NSAMP = 1000;
	#endif

	using namespace Eigen;
	typedef util::SimpleArray<7,double> F7;
	typedef util::SimpleArray<7,uint64_t> I7;	
	boost::random::mt19937 rng((unsigned int)time(0) + 293457);
	boost::uniform_real<> runif;

	XformHash_Quat_BCC7_Zorder<Xform> xh( runif(rng)+0.1, 2+(int)(runif(rng)*60.0), runif(rng)*400+10.0 );

	for(int j = 0; j < NSAMP; ++j){
		XformHash_Quat_BCC7_Zorder<Xform> xh( runif(rng)+0.1, 2+(int)(runif(rng)*20.0), runif(rng)*400+10.0 );
		for(int i = 0; i < NSAMP; ++i){
			Xform x; numeric::rand_xform(rng,x);
			// while( fabs(Quaterniond(x.rotation()).w()) > 0.01 ) numeric::rand_xform(rng,x);
			uint64_t key = xh.get_key(x);
			uint64_t akey,ksym;
			akey = xh.asym_key(key,ksym);
			bool fx = (ksym>>2)&1;
			bool fy = (ksym>>1)&1;
			bool fz = (ksym>>0)&1;
			// cout << "TEST FX " << fx << endl;
			// cout << "TEST FY " << fy << endl;
			// cout << "TEST FZ " << fz << endl;
			// 	{
			// 		uint64_t k = key;
			// 		uint64_t o = k & 1;
			// 		uint64_t w =  util::undilate<7>( k>>4 ) & 63;
			// 		uint64_t x =  util::undilate<7>( k>>5 ) & 63;
			// 		uint64_t y =  util::undilate<7>( k>>6 ) & 63;
			// 		uint64_t z =  util::undilate<7>( k>>7 ) & 63;
			// 		std::cout << "KEY  " << ksym << "\t" << xh.grid_.nside_[3]-o <<"\t" << o << "    " 
			// 		<< w << "\t" << x << "\t" << y << "\t" << z << std::endl;
			// 	}
			// 	{
			// 		uint64_t k = akey;
			// 		uint64_t o = k & 1;
			// 		uint64_t w =  util::undilate<7>( k>>4 ) & 63;
			// 		uint64_t x =  util::undilate<7>( k>>5 ) & 63;
			// 		uint64_t y =  util::undilate<7>( k>>6 ) & 63;
			// 		uint64_t z =  util::undilate<7>( k>>7 ) & 63;
			// 		std::cout << "ASYM " << ksym << "\t" << xh.grid_.nside_[3]-o <<"\t" << o << "    " 
			// 		<< w << "\t" << x << "\t" << y << "\t" << z << std::endl;
			// 	}

			ASSERT_LT( ksym, xh.num_key_symmetries() );
			ASSERT_EQ( key, xh.sym_key(akey,ksym) );
			ASSERT_EQ( key%2, akey%2 ); // preserves parity
			if( ksym==0 ) ASSERT_EQ( key, akey );
			else          ASSERT_NE( key, akey );

			Xform xc = xh.get_center(  key );
			Xform xa = xh.get_center( akey );
			ASSERT_EQ( xc.translation(), xa.translation() );
			Quaterniond qc( xc.rotation() );
			Quaterniond qa( xa.rotation() );
			// cout << qc.coeffs().transpose() << endl;
			// cout << qa.coeffs().transpose() << endl << endl;
			qc = nest::maps::to_half_cell(qc);
			qa = nest::maps::to_half_cell(qa);
			// cout << qc.coeffs().transpose() << endl;
			// cout << qa.coeffs().transpose() << endl << endl;
			ASSERT_NEAR( qc.w(),              qa.w()  , 0.000001 );
			if( fabs( qc.w() ) > 0.000000001  ){
				ASSERT_NEAR( qc.x(), (fx?-1.0:1.0)*qa.x() , 0.000001 );
				ASSERT_NEAR( qc.y(), (fy?-1.0:1.0)*qa.y() , 0.000001 );
				ASSERT_NEAR( qc.z(), (fz?-1.0:1.0)*qa.z() , 0.000001 );
			}
			// cout << std::showpos;

		}
	}
}

// this test is fucked.... can't get 0-registration correct...
// thinking ori + cart_shift should still work correctly.. will test that
TEST( XformHashNeighbors, DISABLED_Quat_BCC7_Zorder_check_cart_neighbors ){
	int NSAMP = 5;
	double THRESH = 0.01;
	#ifdef NDEBUG
	NSAMP = 20;
	THRESH = 0.001;
	#endif
	boost::random::mt19937 rng((unsigned int)time(0) + 4564631);
	boost::uniform_real<> runif;

	// double cart_bound=1.0+runif(rng)*4.0, ang_bound=20.0;
	double cart_bound=2.0, ang_bound=20.0;

	for(int iter1 = 0; iter1 < NSAMP; ++iter1){
		// XformHash_Quat_BCC7_Zorder<Xform> xh( 1.0, 5.0+runif(rng)*10.0, 10+runif(rng)*200.0 );
		XformHash_Quat_BCC7_Zorder<Xform> xh( 1.0, 10.0 );
		XformHashNeighbors< XformHash_Quat_BCC7_Zorder<Xform> > nb( cart_bound, ang_bound, xh, 0 );

		for(int i = 0; i < NSAMP; ++i){
			// cout << "TEST ITER " << i << endl;
			Xform x0; numeric::rand_xform(rng,x0,10.0);
			uint64_t key0 = xh.get_key(x0);
			Xform cen0 = xh.get_center(key0);
			ASSERT_EQ( xh.get_key(cen0), key0 );
			Xform x1(cen0);
			if( key0%2 ){
				x1.translation()[0] -= xh.cart_width()/2.0;
				x1.translation()[1] -= xh.cart_width()/2.0;
				x1.translation()[2] -= xh.cart_width()/2.0;								
			} else {
				x1.translation()[0] += xh.cart_width()/2.0;
				x1.translation()[1] += xh.cart_width()/2.0;
				x1.translation()[2] += xh.cart_width()/2.0;								
			}
			uint64_t key1 = xh.get_key(x1);
			cout << x0.translation().transpose() << endl;
			cout << x1.translation().transpose() << endl;			
			xh.print_key(key0);
			xh.print_key(key1);			
			ASSERT_NE( key0%2, key1%2 );
			Xform cen1 = xh.get_center(key1);

			std::vector<Eigen::Vector3d> vnbrs,vsamp;
			std::set<uint64_t> nbrs_set, nbrs_set_ori, nbrs_set_cart;
			Eigen::Vector3d nbrcen(0,0,0);
			std::vector< util::SimpleArray<3,int16_t> > shifts = nb.get_cart_shifts();
			for(int j =0; j < shifts.size(); ++j){
				uint64_t nbkey0 = xh.cart_shift_key( key0, shifts[j][0], shifts[j][1], shifts[j][2], 0 );
				uint64_t nbkey1 = xh.cart_shift_key( key1, shifts[j][0], shifts[j][1], shifts[j][2], 1 );
				nbrs_set.insert( nbkey0 );
				nbrs_set.insert( nbkey1 );				
				vnbrs.push_back( xh.get_center(nbkey0).translation() );
				vnbrs.push_back( xh.get_center(nbkey1).translation() );				
				nbrcen += xh.get_center(nbkey0).translation();
				nbrcen += xh.get_center(nbkey1).translation();

				nbrs_set_cart.insert( nbkey0 & xh.CART_MASK_NO0 );
				nbrs_set_cart.insert( nbkey1 & xh.CART_MASK_NO0 );				
				nbrs_set_ori.insert( nbkey0 & xh.ORI_MASK_NO0 );
				nbrs_set_ori.insert( nbkey1 & xh.ORI_MASK_NO0 );				
				// cout << nbkey%2 << " " <<  nbrs[j]%2 << endl;
			}

			cout << "NBRS SET ORI: " << endl;
			for(std::set<uint64_t>::const_iterator i = nbrs_set_ori.begin(); i != nbrs_set_ori.end(); ++i){
				cout << "       ";
				cout << ( util::undilate<7>( (*i)>>4 ) & 63 ) << " ";
				cout << ( util::undilate<7>( (*i)>>5 ) & 63 ) << " ";
				cout << ( util::undilate<7>( (*i)>>6 ) & 63 ) << " ";
				cout << ( util::undilate<7>( (*i)>>7 ) & 63 ) << " ";
				cout << (*i)%2;
				cout << endl;
			}
			nbrcen /= shifts.size()*2;

			int nfail=0;
			for(int j = 0; j < 1000*NSAMP; ++j){
				Xform p; numeric::rand_xform_quat(rng,p,cart_bound,0.0);
				uint64_t nk = xh.get_key( p*cen0 );
				vsamp.push_back( xh.get_center(nk).translation() );
				bool fail = nbrs_set.find(nk) == nbrs_set.end();
				nfail += fail;
				// if( fail ) cout << nk%2 << endl;
				ASSERT_TRUE( nbrs_set_cart.find( nk & xh.CART_MASK_NO0) != nbrs_set_cart.end() );
				if( fail ){
					cout << ( util::undilate<7>( nk>>4 ) & 63 ) << " ";
					cout << ( util::undilate<7>( nk>>5 ) & 63 ) << " ";
					cout << ( util::undilate<7>( nk>>6 ) & 63 ) << " ";
					cout << ( util::undilate<7>( nk>>7 ) & 63 ) << " ";
					cout << nk%2;
					cout << endl;

					// ASSERT_TRUE( nbrs_set_ori.find( nk & xh.ORI_MASK_NO0) == nbrs_set_ori.end() );
				}
			}

			std::ofstream out1("test_nbrs.pdb");
			for(int i = 0; i < vnbrs.size(); ++i) io::dump_pdb_atom( out1, "C" ,i, 10.0*vnbrs[i] );
			out1.close();
			std::ofstream out2("test_samp.pdb");
			for(int i = 0; i < vnbrs.size(); ++i) io::dump_pdb_atom( out2, "C" ,i, 10.0*vsamp[i] );
			out2.close();
			std::ofstream out3("test_cen.pdb");
			io::dump_pdb_atom( out3, "C" ,i, 10.0*cen0.translation() );
			io::dump_pdb_atom( out3, "C" ,i, 10.0*cen1.translation() );			
			out3.close();

			ASSERT_LE( nfail/1000.0/NSAMP, THRESH );
			return ;
		}
	}	
}

TEST( XformHashNeighbors, DISABLED_________Quat_BCC7_Zorder_check_ori_neighbors ){
	int NSAMP = 5;
	double THRESH = 0.01;
	#ifdef NDEBUG
	NSAMP = 20;
	THRESH = 0.001;
	#endif
	boost::random::mt19937 rng((unsigned int)time(0) + 94586049);
	boost::uniform_real<> runif;

	double cart_bound=2.0, ang_bound=10.0+runif(rng)*30.0;
	// double cart_bound=1.0, ang_bound=20.0;
	double quat_bound=numeric::deg2quat(ang_bound);

	for(int iter1 = 0; iter1 < NSAMP; ++iter1){
		XformHash_Quat_BCC7_Zorder<Xform> xh( 1.0, 5.0+runif(rng)*10.0, 10+runif(rng)*200.0 );
		// XformHash_Quat_BCC7_Zorder<Xform> xh( 1.0, 10.0 );
		XformHashNeighbors< XformHash_Quat_BCC7_Zorder<Xform> > nb( cart_bound, ang_bound, xh, NSAMP*20.0 );

		for(int i = 0; i < NSAMP; ++i){
			// cout << "TEST ITER " << i << endl;
			Xform x; numeric::rand_xform(rng,x);
			x.translation()[0] = x.translation()[1] = x.translation()[2] = 0.0;
			uint64_t key = xh.get_key(x);
			Xform c = xh.get_center( key );

			std::vector<uint64_t> const & nbrs = nb.get_ori_neighbors( key );
			google::dense_hash_set<uint64_t> nbrs_set;
			nbrs_set.set_empty_key(9999999999);
			Eigen::Vector4d qpass(0,0,0,0), qfail(0,0,0,0);
			int opass=0,ofail=0;
			for(int j =0; j < nbrs.size(); ++j) nbrs_set.insert( nbrs[j] );
			int nfail = 0;
			std::vector<uint64_t> fails,passes;
			for(int j = 0; j < 300*NSAMP; ++j){
				Xform p; numeric::rand_xform_quat(rng,p,cart_bound,quat_bound);
				p.translation()[0] = p.translation()[1] = p.translation()[2] = 0.0;
				uint64_t nk = xh.get_key( p*c ) & xh.ORI_MASK;
				// nk = nk | (uint64_t)1; // debug... this makes it correct much more of the time....
				bool fail = nbrs_set.find(nk) == nbrs_set.end();
				nfail += fail;
				uint64_t k = nk;
				uint64_t o = k & 1;
				uint64_t w =  util::undilate<7>( k>>4 ) & 63;
				uint64_t x =  util::undilate<7>( k>>5 ) & 63;
				uint64_t y =  util::undilate<7>( k>>6 ) & 63;
				uint64_t z =  util::undilate<7>( k>>7 ) & 63;

				if( fail ){ qfail += Eigen::Vector4d( w,x,y,z ); ofail += o; fails.push_back(nk); }
				else      { qpass += Eigen::Vector4d( w,x,y,z ); opass += o; passes.push_back(nk); }
			}

				// cout << fails.size() << " " << passes.size() << endl;
				// std::ofstream out("nbrs_sym_fails.pdb");
				// for(int i = 0; i < fails.size(); ++i){
				// 	Xform x = xh.get_center( fails[i] );
				// 	Eigen::Quaterniond q(x.rotation());
				// 	io::dump_pdb_atom( out, "C" ,i, 50*q.coeffs() );
				// }
				// out.close();
				// std::ofstream out2("nbrs_sym_passes.pdb");
				// for(int i = 0; i < passes.size(); ++i){
				// 	Xform x = xh.get_center( passes[i] );
				// 	Eigen::Quaterniond q(x.rotation());
				// 	io::dump_pdb_atom( out2, "C" ,i, 50*q.coeffs() );
				// }
				// out2.close();


			// cout << "PASS " << qpass.transpose() << " " << opass << endl;
			// cout << "FAIL " << qfail.transpose() << " " << ofail << endl;		
			ASSERT_LE( nfail/300.0/NSAMP, THRESH );

			// return;
		}
	}
}

TEST( XformHashNeighbors, Quat_BCC7_Zorder_check_general_neighbors ){
	int NSAMP = 5;
	double THRESH = 0.003;
	#ifdef NDEBUG
	NSAMP = 10;
	THRESH = 0.001;
	#endif
	boost::random::mt19937 rng((unsigned int)time(0) + 94586049);
	boost::uniform_real<> runif;


	for(int iter1 = 0; iter1 < NSAMP; ++iter1){
		double cart_bound=1.0+runif(rng)*2.0, ang_bound=10.0+runif(rng)*20.0;
		// double cart_bound=2.0, ang_bound=20.0;
		double quat_bound=numeric::deg2quat(ang_bound);
		// XformHash_Quat_BCC7_Zorder<Xform> xh( 1.0, 5.0+runif(rng)*10.0, 10+runif(rng)*200.0 );
		XformHash_Quat_BCC7_Zorder<Xform> xh( 1.0, 10.0 );
		XformHashNeighbors< XformHash_Quat_BCC7_Zorder<Xform> > nb( cart_bound, ang_bound, xh, NSAMP*100.0 );

		for(int iter2 = 0; iter2 < NSAMP; ++iter2){
			Xform x; numeric::rand_xform(rng,x);
			uint64_t key = xh.get_key(x);
			Xform c = xh.get_center( key );
			// xh.print_key(key);
			int ix = (int)((util::undilate<7>( key>>1 ) & 63) | ((key>>57)&127)<<6);
			int iy = (int)((util::undilate<7>( key>>2 ) & 63) | ((key>>50)&127)<<6);
			int iz = (int)((util::undilate<7>( key>>3 ) & 63) | ((key>>43)&127)<<6);

			std::set<uint64_t> nbrs_set; {
				std::vector<uint64_t> const & ori_nbrs = nb.get_ori_neighbors( key );
				std::vector< util::SimpleArray<3,int16_t> > shifts = nb.get_cart_shifts();
				for(int j =0; j < ori_nbrs.size(); ++j){
					uint64_t ori_key = ori_nbrs[j];
					ori_key = xh.cart_shift_key( ori_key, ix, iy, iz );
					// xh.print_key( ori_key );
					for(int k=0; k < shifts.size(); ++k){
						uint64_t nbkey0 = xh.cart_shift_key( ori_key, shifts[k][0], shifts[k][1], shifts[k][2], 0 );
						uint64_t nbkey1 = xh.cart_shift_key( ori_key, shifts[k][0], shifts[k][1], shifts[k][2], 1 );
						nbrs_set.insert( nbkey0 );
						nbrs_set.insert( nbkey1 );
					}
				}
				// cout << nbrs_set.size() << " " << ori_nbrs.size() << " " << shifts.size() << endl;
			}
			
			int nfail = 0;
			for(int j = 0; j < 300*NSAMP; ++j){
				Xform p; numeric::rand_xform_quat(rng,p,cart_bound,quat_bound);
				uint64_t nk = xh.get_key( c*p );
				// Xform ncen = xh.get_center( nk );
				bool fail = nbrs_set.find(nk) == nbrs_set.end();
				nfail += fail;
				// if( fail ){
				// 	xh.print_key(key);					
				// 	xh.print_key(nk);
				// 	cout << endl;
				// }
			}
			double failfrac = nfail/300.0/NSAMP;
			cout << "TEST ITER " << iter1 << " " << iter2 << " " << cart_bound << " " << ang_bound << " " << failfrac << endl;
			ASSERT_LE( failfrac, THRESH );

		}
	}
}


}}}}



























