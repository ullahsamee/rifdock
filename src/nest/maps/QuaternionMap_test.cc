#include <gtest/gtest.h>

#include <nest/NEST.hh>
#include <nest/maps/QuaternionMap.hh>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "io/dump_pdb_atom.hh"

#include <fstream>
#include <Eigen/Geometry>

namespace scheme {
namespace nest {
namespace maps {

using std::cout;
using std::endl;

TEST( QuaternionMap, test_cell_validity_check ){
	ASSERT_EQ( (int) 1.5 ,  1 ); // check round towards 0 ... IEEE?
	ASSERT_EQ( (int)-1.5 , -1 );
	typedef QuaternionMap<4> MapType;
	typedef MapType::Params P;
	MapType qmap;
	MapType::ValueType val;
	P q;

	// q = P( 0.00, 0.00, 0.00, 0.00 ); ASSERT_TRUE ( qmap.params_to_value( q/2.0+0.5, 0, 2, val ) );
	// q = P(-0.01, 0.00, 0.00, 0.00 ); ASSERT_FALSE( qmap.params_to_value( q/2.0+0.5, 0, 2, val ) );
	// q = P( 0.01,-0.01, 0.01,-0.01 ); ASSERT_TRUE ( qmap.params_to_value( q/2.0+0.5, 0, 2, val ) );
	// q = P( 0.01,-0.01, 0.01,-0.01 ); ASSERT_FALSE( qmap.params_to_value( q/2.0+0.5, 0, 3, val ) );

	NEST<4,P,QuaternionMap> nest;

	// test that a reasonable fraction of values are invalid
	int count = 0;
	for( size_t resl = 2; resl < 6; ++resl){
		int nvalid = 0;
		for( size_t i = 0; i < nest.size(resl); ++i){
			if( nest.set_state(i,resl) ){
				++nvalid;
			}
			++count;
		}
		double empirical_frac_valid = (double)nvalid/(double)nest.size(resl);
		// pi*pi is hyper area of 1/2 3sphere
		// 4.0*1<<resl is 2 * cell width
		// 16.0 is vol of 4-cube
		double analytic_frac_hypervolume_aprox = 3.14159*3.14159*(4.0/(1<<resl)) / 16.0; 
		cout << "resl " << resl << " frac valid " << empirical_frac_valid 
		     << " vs. " << analytic_frac_hypervolume_aprox << " approx cell-width hyper-shell volume fraction" << endl;
		    ASSERT_LT( empirical_frac_valid, analytic_frac_hypervolume_aprox );
		    if(resl > 2)
		    	ASSERT_GT( empirical_frac_valid, analytic_frac_hypervolume_aprox/2.0 );

	}
	cout << "Num QuaternionMap checks: " << (double)count/1000000.0 << endl;

	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> gauss;

	// test that the indices of all valid quats are valid for set_state
	for(int i = 0; i < 10000; ++i){
		q = P( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
		q = q / q.norm();
		// cout << q << endl;
		for(size_t resl = 0; resl < 15; ++resl){
			size_t i_should_be_valid = nest.get_index(q,resl);
			ASSERT_TRUE( nest.set_state( i_should_be_valid, resl ) );
		}
	}


	BOOST_STATIC_ASSERT( util::meta::has_subscript_oper<P,double&,size_t>::value );
	BOOST_STATIC_ASSERT( util::meta::has_const_subscript_oper<P,double const &,size_t>::value );

	cout << "attempt to dump visualiztion" << endl;
	{
		typedef Eigen::Vector3d V;
		NEST<4,Eigen::Quaternion<double>,QuaternionMap> nest;
		// for(size_t r = 4; r < 5; ++r){
		size_t r = 4; {
			std::ofstream out("quatmaptest.pdb");
			for(size_t i = 0; i < nest.size(r); ++i){
				if( nest.set_state(i,r) ){
					Eigen::Quaternion<double> q = nest.value();
					cout << r << " " << i << " " << q.coeffs().transpose() << endl;
					 V ximg = q.matrix() * Eigen::Vector3d(1,0,0);
					 V yimg = q.matrix() * Eigen::Vector3d(0,1.1,0);
					 V zimg = q.matrix() * Eigen::Vector3d(0,0,1.2);
					 io::dump_pdb_atom(out,"O" ,30*ximg);
					 io::dump_pdb_atom(out,"NI",30*yimg);
					 io::dump_pdb_atom(out,"N" ,30*zimg);
				}
			}
			out.close();
		}
	}


}


}
}
}
