#include <gtest/gtest.h>

#include "numeric/util.hh"
#include "io/dump_pdb_atom.hh"

#include "nest/maps/TetracontoctachoronMap.hh"
#include "nest/NEST.hh"

#include <Eigen/Dense>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <boost/timer/timer.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>

namespace scheme { namespace nest { namespace maps {

using namespace Eigen;
using std::cout;
using std::endl;


TEST(tetracontoctachoron,cell_lookup){
	// HecatonicosachoronMap<> map;
	// HecatonicosachoronMap<>::Params p;
	// for( size_t i = 0; i < 60; ++i){
	// 	size_t cell_index = 999;
	// 	map.value_to_params( cellcen<double>(i).matrix(), 0, p, cell_index ); 
	// 	// cout << p << endl;
	// 	ASSERT_EQ( i, cell_index );
	// 	ASSERT_NEAR( p[0], 0.5, 0.0000001 );
	// 	ASSERT_NEAR( p[1], 0.5, 0.0000001 );
	// 	ASSERT_NEAR( p[2], 0.5, 0.0000001 );
	// }

	NEST<3,Matrix3d,TetracontoctachoronMap> nest;
	// ASSERT_EQ( 0, nest.get_index( nest.set_and_get(0,2), 2 ));

	for(int r = 0; r <= 4; ++r){
		for(int i = 0; i < (int)nest.size(r); ++i){
			if( nest.set_state(i,r) ){
				size_t ilookup = nest.get_index( nest.value(), r );
				size_t ci0 = i >> (3*r);
				size_t ciL = ilookup >> (3*r);
				if( ci0 == ciL ) ASSERT_EQ( i , ilookup );
			}
		}
	}

}


TEST(tetracontoctachoron,covering){
	// cout << "QuaternionMap Covrad" << endl;
	int NRES = 6;
	int ITERS = 1*1000*1000;

	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> gauss;
	boost::uniform_real<> uniform;

	NEST<3,Matrix3d,TetracontoctachoronMap> nest;
	for(int r = 0; r <= NRES; ++r){
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i <= ITERS; ++i){
			Eigen::Quaterniond q( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
			q.normalize();
			// cout << q << endl;
			Eigen::Quaterniond qcen( nest.set_and_get( nest.get_index(q.matrix(),r) , r ) );
			maxdiff = std::max(maxdiff,q.angularDistance(qcen));
			avgdiff += q.angularDistance(qcen);
		}
		avgdiff /= ITERS;
		size_t count = 0; for(size_t i = 0; i < nest.size(r)/24; ++i) if(nest.set_state(i,r)) ++count;
		count *= 24;
		// size_t count = nest.size(r);
		double volfrac = (double)count*(maxdiff*maxdiff*maxdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;
		double avgfrac = (double)count*(avgdiff*avgdiff*avgdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;		
		printf("%2i %16lu %10.5f %10.5f %10.5f %10.5f %10.5f %7.5f\n", 
			r,
			count,
			maxdiff*180.0/M_PI,
			avgdiff*180.0/M_PI,
			maxdiff/avgdiff,
			volfrac,
			avgfrac,
			(double)count / nest.size(r)
		);
	}

}

TEST(tetracontoctachoron,DISABLE_visualize){

	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> gauss;
	boost::uniform_real<> uniform;

	Quaterniond qrand( gauss(rng), gauss(rng), gauss(rng), gauss(rng) );
	// Quaterniond qrand( 1,0,0,0 );
	qrand.normalize();
	Vector3d X = qrand*Vector3d(1,0  ,0  );
	Vector3d Y = qrand*Vector3d(0,1.2,0  );
	Vector3d Z = qrand*Vector3d(0,0  ,1.4);

	NEST<3,Matrix3d,TetracontoctachoronMap> nest;
	// size_t beg = 0;
	// while(!nest.set_state(beg,10)) beg = std::max<size_t>(uniform(rng)*(nest.size(10)-1000),0);
	for(size_t r = 0; r <= 8; ++r){
		int N = 4096;
		// int beg = std::max( 0, (int)nest.size(r)/12 - N/2 );
		int beg = 0;
		std::ofstream out(("tcoc_"+boost::lexical_cast<std::string>(r)+".pdb").c_str());
		size_t count1 = 0, count2 = 0;
		// cout << r << " " << nest.size(r) << " " << (beg>>(4*(10-r))) << endl;
		// continue;
		// for(size_t i = beg>>(4*(10-r)); i < nest.size(r); ++i){
		for(size_t i = beg; i < nest.size(r); ++i){		
			++count1;
			if( nest.set_state(i,r) ){
				++count2;
				if( count1 > N) break;
				Matrix3d m = nest.value();
				// cout << r << " " << i << " " << q.coeffs().transpose() << endl;
				 Vector3d ximg = m * X;
				 Vector3d yimg = m * Y;
				 Vector3d zimg = m * Z;;
				 // cout << r << " " << nest.cell_index(i,r) << endl;
				 io::dump_pdb_atom(out, "O" ,50*ximg);
				 io::dump_pdb_atom(out, "NI",50*yimg);
				 io::dump_pdb_atom(out, "N" ,50*zimg);
			}
		}
		out.close();
		cout << r << " " << count2 / (double)count1 << endl;
	}

}


}}}
