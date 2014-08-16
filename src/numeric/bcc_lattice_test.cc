#include <gtest/gtest.h>

#include "numeric/bcc_lattice.hh"
#include "io/dump_pdb_atom.hh"
#include <fstream>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>


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
	BCC<N,F,S> bcc(I(100),V(0),V(10));

	int NSAMP = 10*1000*1000;
	std::vector<V> samples(NSAMP);
	for(int i = 0; i < NSAMP; ++i)
		for(int j = 0; j < N; ++j)
			samples[i][j] = uniform(rng)*6+2;

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

TEST(bcc_lattice,performance){
	test_bcc_performance<3,double,size_t>();
	test_bcc_performance<4,double,size_t>();
	test_bcc_performance<5,double,size_t>();
	test_bcc_performance<6,double,size_t>();
	test_bcc_performance<7,double,size_t>();		
}

}}}

