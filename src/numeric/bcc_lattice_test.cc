#include <gtest/gtest.h>

#include "numeric/bcc_lattice.hh"
#include "io/dump_pdb_atom.hh"
#include <fstream>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>

#include <sparsehash/dense_hash_set>

#include <Eigen/Geometry>

namespace scheme { namespace numeric { namespace test {

using std::cout;
using std::endl;

TEST(bcc_lattice,centers_map){
	typedef util::SimpleArray<3,double> V;
	typedef util::SimpleArray<3,size_t> I;

	BCC<3,double> bcc(I(3,5,7),V(0,0,0),V(6,10,14));

	std::ofstream out("test.pdb");
	for(size_t i = 0; i < bcc.size(); i+=2){
		ASSERT_EQ( i, bcc[bcc[i]] );
		// printf("%6lu %10.6f %10.6f %10.6f\n",i,bcc[i][0],bcc[i][1],bcc[i][2]);
		// io::dump_pdb_atom(out,i,bcc[i]*3.0);
		// cout << i << " " << bcc[i] << endl;
	}
	out.close();

}

template<int N, class F, class S>
F test_bcc_performance(){
	typedef util::SimpleArray<N,F> V;
	typedef util::SimpleArray<N,S> I;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::uniform_real<> uniform;
	BCC<N,F,S> bcc(I(100),V(-5),V(5));

	int NSAMP = 5*1000*1000;
	std::vector<V> samples(NSAMP);
	for(int i = 0; i < NSAMP; ++i)
		for(int j = 0; j < N; ++j)
			// samples[i][j] = uniform(rng)*9.4+0.3;
			samples[i][j] = uniform(rng)*9.0-4.5;			

	std::vector<size_t> indices(NSAMP);
	boost::timer::cpu_timer lookup_time;
	for(int i = 0; i < NSAMP; ++i)
		indices[i] = bcc[samples[i]];
	cout << N << " lookup: " << (double)NSAMP / lookup_time.elapsed().wall*1000000000.0 << " sec" << endl;

	std::vector<V> centers(NSAMP);
	boost::timer::cpu_timer getval_time;
	for(int i = 0; i < NSAMP; ++i)
		centers[i] = bcc[indices[i]];
	cout << N << " getval: " << (double)NSAMP / getval_time.elapsed().wall*1000000000.0 << " sec" << endl;


	F maxdiff = 0;
	for(int i = 0; i < NSAMP; ++i)
		maxdiff = std::max( maxdiff, (samples[i]-centers[i]).squaredNorm() );

	F frac = bcc.width_[0] * sqrt(N)/2.0 / sqrt(maxdiff);
	F improvement = 1.0; for(int i = 0; i < N; ++i) improvement *= frac;
	cout << N << " improvement over cubic: " << improvement / 2.0 << " cov: " << maxdiff << " vs. " << sqrt(N)/2.0 << endl; // 2 x num samp as cubic

	return maxdiff;


}

// TEST(bcc_lattice,performance){
// 	test_bcc_performance<3,double,size_t>();
// 	test_bcc_performance<4,double,size_t>();
// 	test_bcc_performance<5,double,size_t>();
// 	test_bcc_performance<6,double,size_t>();
// 	test_bcc_performance<7,double,size_t>();		
// }

// QuaternionMap Covrad
//  1                8  119.36447   61.22694    1.94954    3.83748    0.51790
//  2              128   59.04495   27.31148    2.16191    7.43172    0.73549
//  3             1040   27.85890   13.60802    2.04724    6.34244    0.73918
//  4             8480   13.59714    6.82994    1.99081    6.01270    0.76204
//  5            68416    6.91475    3.40731    2.02939    6.37995    0.76335
//  6           548496    3.44178    1.70610    2.01734    6.30746    0.76828
//  7          4392160    1.73778    0.85268    2.03801    6.50123    0.76802
TEST(bcc_lattice,orientatin_coverage){
	typedef double F;
	typedef size_t S;
	typedef util::SimpleArray<4,F> V;
	typedef util::SimpleArray<4,S> I;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> rnorm;
	size_t Nside = 64;
	BCC<4,F,S> bcc(I(Nside),V(-1.0-2.0/Nside),V(1.0+2.0/Nside));
	cout << bcc.lower_ << endl;
	cout << bcc.width_ << endl;
	cout << bcc.lower_+bcc.width_*Nside << endl;

	int NSAMP = 20*1000*1000;
	google::dense_hash_set<size_t> idx_seen;
	idx_seen.set_empty_key(999999999999);
	F maxdiff=0,avgdiff=0;
	for(int i = 0; i < NSAMP; ++i){
		V sample( rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng) );
		sample /= sample.norm();
		if(sample[0]<0) sample = -sample;
		ASSERT_DOUBLE_EQ( sample.norm(), 1.0 );
		size_t index = bcc[sample];
		idx_seen.insert(index);
		V center = bcc[index];
		Eigen::Quaterniond samp(sample[0],sample[1],sample[2],sample[3]);
		Eigen::Quaterniond cen (center[0],center[1],center[2],center[3]);
		cen.normalize();
		maxdiff = std::max( maxdiff, samp.angularDistance(cen) );
		avgdiff += samp.angularDistance(cen);
	}
	avgdiff /= NSAMP;
	size_t count = idx_seen.size();
	cout << "idxfrac: " << (double)count / (Nside*Nside*Nside*Nside) << " xcov " << (double)NSAMP/count << endl;
	double volfrac = (double)count*(maxdiff*maxdiff*maxdiff) * 4.0/3.0*M_PI / 8.0/M_PI/M_PI;
	double avgfrac = (double)count*(avgdiff*avgdiff*avgdiff) * 4.0/3.0*M_PI / 8.0/M_PI/M_PI;
	cout << "RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac      AvgFrac" << endl;
	printf("%2lu %16lu %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
				Nside, count, maxdiff*180.0/M_PI, avgdiff*180.0/M_PI, maxdiff/avgdiff, volfrac, avgfrac );


}

}}}

