#include <gtest/gtest.h>

#include "nest/NEST.hh"
#include "nest/maps/EulerAnglesMap.hh"

#include <Eigen/Geometry>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

namespace scheme { namespace nest { namespace maps { namespace test {

using std::cout;
using std::endl;

TEST(EulerAnglesMap,DISABLED_covering){
	using namespace Eigen;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> gauss;
	boost::uniform_real<> uniform;

	cout << "EulerAnglesMap Covrad" << endl;
	int NRES = 7;
	// int const ITERS = 1000000;
	int ITERS = 1000000;	
	NEST<3,Matrix3d,EulerAnglesMap> nest;
	for(int r = 1; r <= NRES; ++r){
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < ITERS; ++i){
			Eigen::Quaterniond q( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
			q.normalize();
			Matrix3d m = nest.set_and_get( nest.get_index(q.matrix(),r) , r );
			Quaterniond qcen(m);
			// if( q.angularDistance(qcen) > maxdiff ){
			// 	RowVector3d euler; numeric::euler_angles(q.matrix(),euler);
			// 	euler[0] /= M_PI*2.0;
			// 	euler[1] /= M_PI*2.0;
			// 	euler[2] /= M_PI;
			// 	cout << r << " " << maxdiff << " " << euler << endl;
			// }
			avgdiff += q.angularDistance(qcen);
			maxdiff = std::max(maxdiff,q.angularDistance(qcen));
		}
		avgdiff /= ITERS;
		// size/2 because half samples are ignored
		double volfrac = (double)nest.size(r)/2*(maxdiff*maxdiff*maxdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;
		double avgfrac = (double)nest.size(r)/2*(avgdiff*avgdiff*avgdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;
		printf("%2i %16lu %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
			r, nest.size(r)/2, maxdiff*180.0/M_PI, avgdiff*180.0/M_PI, maxdiff/avgdiff, volfrac, avgfrac );
		// cout << boost::format("%2i %20i %.7d %.7d") % r % nest.size(r) % (maxdiff*180.0/M_PI) % volfrac << endl;
	}

}


TEST(EulerAnglesMap,shapes){
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> gauss;
	boost::uniform_real<> uniform;

	/// inspect maxdiff/avgdiff for some simple shapes
	int ITERS = 100000;
	{
		util::SimpleArray<3,double> b(1,1,1);
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < ITERS; ++i){
			util::SimpleArray<3,double> samp(uniform(rng),uniform(rng),uniform(rng));
			samp = (samp-0.5) * b;
			if(samp.norm() > 0.5){ --i; continue; }
			avgdiff += samp.norm();
			maxdiff = std::max(samp.norm(),maxdiff);
		}
		avgdiff /= ITERS;
		cout << "sphere: " << maxdiff / avgdiff << endl;
	}
	{
		util::SimpleArray<3,double> b(1,1,1);
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < ITERS; ++i){
			util::SimpleArray<3,double> samp(uniform(rng),uniform(rng),uniform(rng));
			samp = (samp-0.5) * b;
			avgdiff += samp.norm();
			maxdiff = std::max(samp.norm(),maxdiff);
		}
		avgdiff /= ITERS;
		cout << "square: " << maxdiff / avgdiff << endl;
	}
	{
		util::SimpleArray<3,double> b(2,1,1);
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < ITERS; ++i){
			util::SimpleArray<3,double> samp(uniform(rng),uniform(rng),uniform(rng));
			samp = (samp-0.5) * b;
			avgdiff += samp.norm();
			maxdiff = std::max(samp.norm(),maxdiff);
		}
		avgdiff /= ITERS;
		cout << "rect211: " << maxdiff / avgdiff << endl;
	}
	{
		util::SimpleArray<3,double> b(1,1,1), cen(0.135022,0.135022,0.135022);
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < ITERS; ++i){
			util::SimpleArray<3,double> samp(uniform(rng),uniform(rng),uniform(rng));
			samp = (samp-0.5) * b;
			if( samp.sum() < 0 ){ --i; continue; }
			avgdiff += (samp-cen).norm();
			maxdiff = std::max((samp-cen).norm(),maxdiff);
		}
		avgdiff /= ITERS;
		cout << "triang: " << maxdiff / avgdiff << endl;
	}

}



}}}}

