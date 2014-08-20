#include <gtest/gtest.h>

#include "numeric/bcc_lattice.hh"
#include "io/dump_pdb_atom.hh"
#include <fstream>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/foreach.hpp>
#include <boost/timer/timer.hpp>
#include <iterator>     // std::back_inserter
#include <boost/format.hpp>
#include <sparsehash/dense_hash_set>

#include <Eigen/Geometry>

namespace scheme { namespace numeric { namespace test {

using std::cout;
using std::endl;


TEST(bcc_lattice,centers_map){
	typedef util::SimpleArray<3,double> V;
	typedef util::SimpleArray<3,size_t> I;

	BCC<3,double> bcc(I(3,5,7),V(0,0,0),V(6,10,14));

	// std::ofstream out("test.pdb");
	for(size_t i = 0; i < bcc.size(); i+=2){
		ASSERT_EQ( i, bcc[bcc[i]] );
		// printf("%6lu %10.6f %10.6f %10.6f\n",i,bcc[i][0],bcc[i][1],bcc[i][2]);
		// io::dump_pdb_atom(out,i,bcc[i]*3.0);
		// cout << i << " " << bcc[i] << endl;
	}
	// out.close();

}

	// double const R3approx = 0.558099;
	// double const R4approx = 0.701687;
	// double const R5approx = 0.742306;
	// double const R6approx = 0.845359;
	// double const R7approx = 0.882879;
template<int N, class F, class S>
F
test_bcc_performance(
	int NSAMP,
	S const Nside,
	F const Width
){
	typedef util::SimpleArray<N,F> V;
	typedef util::SimpleArray<N,S> I;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::uniform_real<> runif;
	BCC<N,F,S> bcc(I(Nside),V(-Width/2),V(Width/2));

	std::vector<V> samples(NSAMP);
	for(int i = 0; i < NSAMP; ++i)
		for(int j = 0; j < N; ++j)
			// samples[i][j] = runif(rng)*9.4+0.3;
			samples[i][j] = runif(rng)*9.0-4.5;			

	std::vector<size_t> indices(NSAMP);
	boost::timer::cpu_timer lookup_time;
	for(int i = 0; i < NSAMP; ++i)
		indices[i] = bcc[samples[i]];
	cout << "BCC DIM " << N << " lookup rate: " << (double)NSAMP / lookup_time.elapsed().wall*1000000000.0 << " sec / ";

	std::vector<V> centers(NSAMP);
	boost::timer::cpu_timer getval_time;
	for(int i = 0; i < NSAMP; ++i)
		centers[i] = bcc[indices[i]];
	cout << " getval rate: " << (double)NSAMP / getval_time.elapsed().wall*1000000000.0 << " sec" << endl;


	F maxdiff = 0;
	for(int i = 0; i < NSAMP; ++i)
		maxdiff = std::max( maxdiff, (samples[i]-centers[i]).squaredNorm() );

	maxdiff = sqrt(maxdiff);
	F frac = bcc.width_[0] * sqrt(N)/2.0 / maxdiff;
	F improvement = 1.0; for(int i = 0; i < N; ++i) improvement *= frac;
	cout << "     improvement over cubic: " << improvement / 2.0 << " cov: " << maxdiff / Width * Nside
	     << " vs. " << sqrt(N)/2.0 << endl; // 2 x num samp as cubic

	return maxdiff;
}

TEST(bcc_lattice,performance){
	size_t Nside = 100;
	double Width = 10.0;
	double const R3test = test_bcc_performance<3,double,size_t>( 100000, Nside, Width );
	                      test_bcc_performance<4,double,size_t>( 100000, Nside, Width );
	                      test_bcc_performance<5,double,size_t>( 100000, Nside, Width );
	                      test_bcc_performance<6,double,size_t>( 100000, Nside, Width );
	                      test_bcc_performance<7,double,size_t>( 100000, Nside, Width );
	double const R3 = std::pow(2.0,-5.0/3.0)*sqrt(5) / std::pow(2.0,1.0/3.0) * Width / Nside;
	ASSERT_LE( R3test     , R3 );
	ASSERT_GT( R3test*1.1 , R3 );	
}


template<int N, class F, class S>
F
test_bcc_inradius(){
	typedef util::SimpleArray<N,F> V;
	typedef util::SimpleArray<N,S> I;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;	
	S const Nside = 5;
	BCC<N,F,S> bcc(I(5),V(-(F)Nside/2.0),V((F)Nside/2.0));
	BOOST_VERIFY( bcc[bcc[V(0.0)]] == V(0.0) );
	BOOST_VERIFY( bcc.width_ == V(1.0) );
	S const i0 = bcc[V(0)];

	double const RNapprox[5] = { 0.558099, 0.701687, 0.742306, 0.845359, 0.882879 };
	double const Rapprox = RNapprox[N-3];

	int NSAMP = 50*1000;

	double min_inrad = Rapprox;
	for(int i = 0; i < NSAMP; ++i){
		V samp;
		for(int j = 0; j < N; ++j) samp[j] = rnorm(rng);
		samp.normalize();
		double const radius = runif(rng)*min_inrad*0.1 + 0.9*min_inrad;
		samp *= radius;
		if( bcc[samp] != i0 ){
			min_inrad = radius;
		}

	}
	return min_inrad;
}

TEST(bcc_lattice,inradius){
	ASSERT_NEAR( (test_bcc_inradius<3,double,size_t>()), 0.433015, 0.03 );
	ASSERT_NEAR( (test_bcc_inradius<4,double,size_t>()), 0.500000, 0.03 );
	ASSERT_NEAR( (test_bcc_inradius<5,double,size_t>()), 0.500000, 0.03 );
	ASSERT_NEAR( (test_bcc_inradius<6,double,size_t>()), 0.500000, 0.03 );
	ASSERT_NEAR( (test_bcc_inradius<7,double,size_t>()), 0.500000, 0.03 );
}

template<int N, class F, class S> 
F 
test_bcc_neighbors( int NSAMP ){
	typedef util::SimpleArray<N,F> V;
	typedef util::SimpleArray<N,S> I;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;	
	S const Nside = 5;
	BCC<N,F,S> bcc(I(5),V(-(F)Nside/2.0),V((F)Nside/2.0));
	BOOST_VERIFY( bcc[bcc[V(0.0)]] == V(0.0) );
	BOOST_VERIFY( bcc.width_ == V(1.0) );
	S const i0 = bcc[V(0)];
	std::vector<size_t> nbrs,nbrs_we;		
	bcc.neighbors( i0, std::back_inserter(nbrs   ), false );
	bcc.neighbors( i0, std::back_inserter(nbrs_we), true );		

	Cubic<N,F,S> cubic(I(5),V(-(F)Nside/2.0),V((F)Nside/2.0));
	BOOST_VERIFY( cubic[cubic[V(0.0)]] == V(0.0) );
	BOOST_VERIFY( cubic.width_ == V(1.0) );
	S const i0cubic = cubic[V(0)];
	std::vector<size_t> nbrs_cubic;		
	cubic.neighbors( i0cubic, std::back_inserter(nbrs_cubic) );

	/////////////////////////////////////////////////////
	// test neighbor coverage
	///////////////////////////////////////////////////////////////

	F maxrad_99 = 0;
	// F const RNapprox[5] = { 0.558099, 0.701687, 0.742306, 0.845359, 0.882879 };
	// F const Rapprox = RNapprox[N-3];
	F inrad = 0.5;
	if(N==3) inrad = 0.433015;
	// F const radius = 1.9*inrad;
	for(F radius = 0; radius < 3.0*inrad; radius += inrad/10.0){

		int sum_in_nbrs=0, sum_in_nbrs_we=0;
		for(int i = 0; i < NSAMP; ++i){

			// pick random point in i0 cell
			V cen; for(int j = 0; j < N; ++j) cen[j] = (runif(rng)-0.5)*2.0*inrad;
			if( bcc[cen] != i0 ){ --i; continue; }

			V samp;
			for(int j = 0; j < N; ++j) samp[j] = rnorm(rng);
			samp *= (radius/samp.norm());
			S index = bcc[samp+cen];
			// if( index == i0 ){ --i; continue; }
			if( std::find( nbrs_we.begin(), nbrs_we.end(), index ) != nbrs_we.end() ){
				++sum_in_nbrs_we;
				if( std::find( nbrs.begin(), nbrs.end(), index ) != nbrs.end() ){
					++sum_in_nbrs;
				}
			}
		}

		// uncomment fto print report
		// int sum_in_nbrs_cubic=0;
		// for(int i = 0; i < NSAMP; ++i){
		// 	// pick random point in i0 cell
		// 	V cen; for(int j = 0; j < N; ++j) cen[j] = (runif(rng)-0.5);
		// 	BOOST_VERIFY( cubic[cen] == i0cubic );

		// 	V samp;
		// 	for(int j = 0; j < N; ++j) samp[j] = rnorm(rng);
		// 	samp *= (radius/inrad*0.5/samp.norm());
		// 	S index = cubic[samp+cen];
		// 	// if( index == i0 ){ --i; continue; }
		// 	if( std::find( nbrs_cubic.begin(), nbrs_cubic.end(), index ) != nbrs_cubic.end() ){
		// 		++sum_in_nbrs_cubic;
		// 	}
		// }
		// printf("%i %7.3f %9.7f %9.7f %9.7f\n",
		// 	N,
		// 	radius/inrad,
		// 	(F)sum_in_nbrs/NSAMP,
		// 	(F)sum_in_nbrs_we/NSAMP,
		// 	(F)sum_in_nbrs_cubic/NSAMP			
		// );

		if( (F)sum_in_nbrs_we/NSAMP > 0.99 ) maxrad_99 = radius;
	}

	// ///////////////////////////////////////////////////////////////////
	// // // dump neighbors
	// //////////////////////////////////////////////////////////////////
	// if(N != 3) return;
	// #define LAT bcc
	// std::ofstream out_bcc("bcc.pdb");
	// for(int i = 0; i < LAT.size(); ++i) io::dump_pdb_atom(out_bcc,i,LAT[i]*10.0);
	// out_bcc.close();
	// for(int i = 0; i < LAT.size(); ++i){
	// 	std::vector<size_t> nbrs;
	// 	LAT.neighbors( i, std::back_inserter(nbrs), true );
	// 	// BOOST_FOREACH(size_t nbr, nbrs) cout << nbr << " " << LAT[nbr] << endl;
	// 	std::string s = boost::str(boost::format("%4i") % i);
	// 	std::ofstream out("test_"+s+".pdb");
	// 	BOOST_FOREACH(size_t nbr, nbrs) io::dump_pdb_atom(out,nbr,LAT[nbr]*10.0);
	// 	out.close();
	// }

	return maxrad_99;
}

TEST(bcc_lattice,neighbors){
	// for(int i = 3; i < 8; ++i){
	// 	int nc = std::pow(3,i);
	// 	int nbccFC = 1+2*i+std::pow(2,i);
	// 	int nbccFCE = nbccFC + i*(i-1)/2 * 4;
	// 	cout << "Nnbrs: " << i <<" cubic "<< nc << " bccFC " << nbccFC << " bccFCE " << nbccFCE << endl;
	// }
	// Nnbrs: 3 cubic   27 bccFC  15 bccFCE  27
	// Nnbrs: 4 cubic   81 bccFC  25 bccFCE  49
	// Nnbrs: 5 cubic  243 bccFC  43 bccFCE  83
	// Nnbrs: 6 cubic  729 bccFC  77 bccFCE 137
	// Nnbrs: 7 cubic 2187 bccFC 143 bccFCE 227
	ASSERT_LE( 0.76, (test_bcc_neighbors<3,double,size_t>(5000)) );
	ASSERT_LE( 0.69, (test_bcc_neighbors<4,double,size_t>(5000)) );
	ASSERT_LE( 0.64, (test_bcc_neighbors<5,double,size_t>(5000)) );
	ASSERT_LE( 0.59, (test_bcc_neighbors<6,double,size_t>(5000)) );
	ASSERT_LE( 0.54, (test_bcc_neighbors<7,double,size_t>(5000)) );
	// cout << test_bcc_neighbors<3,double,size_t>(1000000) << endl; // 0.779427
	// cout << test_bcc_neighbors<4,double,size_t>(1000000) << endl; // 0.75
	// cout << test_bcc_neighbors<5,double,size_t>(1000000) << endl; // 0.7
	// cout << test_bcc_neighbors<6,double,size_t>(1000000) << endl; // 0.65
	// cout << test_bcc_neighbors<7,double,size_t>(1000000) << endl; // 0.6
}

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
	size_t Nside = 32;
	BCC<4,F,S> bcc(I(Nside),V(-1.0-2.0/Nside),V(1.0+2.0/Nside));
	// cout << bcc.lower_ << endl;
	// cout << bcc.width_ << endl;
	// cout << bcc.lower_+bcc.width_*Nside << endl;

	int NSAMP = 1*1000*1000;
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
	double xcov = (double)NSAMP/count;
	if(xcov < 50) cout << "volfrac will be low if xcov < 50 or so... (not enough samples)" << endl;
	cout << "idxfrac: " << (double)count / (Nside*Nside*Nside*Nside) << " xcov " << xcov << endl;
	double volfrac = (double)count*(maxdiff*maxdiff*maxdiff) * 4.0/3.0*M_PI / 8.0/M_PI/M_PI;
	double avgfrac = (double)count*(avgdiff*avgdiff*avgdiff) * 4.0/3.0*M_PI / 8.0/M_PI/M_PI;
	cout << "RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac      AvgFrac" << endl;
	printf("%2lu %16lu %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
				Nside, count, maxdiff*180.0/M_PI, avgdiff*180.0/M_PI, maxdiff/avgdiff, volfrac, avgfrac );

	ASSERT_LE( volfrac , 3.5 );

}

// TEST(bcc_lattice,coverage_transform_7d){
// 	using namespace Eigen;
// 	typedef Transform<double,3,AffineCompact> Xform;
// 	typedef util::SimpleArray<7,double> V;
// 	typedef util::SimpleArray<7,size_t> I;
// 	typedef Matrix<double,7,1> Vector7d;
// 	boost::random::mt19937 mt((unsigned int)time(0));
// 	boost::normal_distribution<> rnorm;
// 	size_t Nside = 64;
// 	V bounds = V(1.5,1.5,1.5,1.5,10,10,10);
// 	BCC<7,double> bcc(I(Nside),-bounds,bounds);

// 	Matrix<double,3,6> pts0;
// 	pts0 <<  1, 0, 0,-2, 0, 0,
// 	         0, 1, 0, 0,-1, 0,
// 	         0, 0, 1, 0, 0,-1;
// 	pts0.colwise() -= pts0.rowwise().sum()/pts0.cols();

// 	Xform X( AngleAxisd(2,Vector3d(1,2,3)) );

// 	// cout << pts0 << endl;
// 	// cout << X*pts0 << endl;

// 	int const NSAMP = 3;
// 	for(int i = 0; i < NSAMP; ++i){
// 		Vector4d quat( rnorm(mt), rnorm(mt), rnorm(mt), rnorm(mt) ); quat.normalize();
// 		Vector3d trans( rnorm(mt), rnorm(mt), rnorm(mt) );
// 		Vector7d samp;
// 		samp.block(0,0,4,1) = quat;
// 		samp.block(4,0,3,1) = trans;
// 		// cout << samp.transpose() << endl;
// 		V tmp = bcc[ bcc[ samp ] ];
// 		Vector7d cen;
// 		for(int i = 0; i < 7; ++i) cen[i] = tmp[i];

// 		cout << samp.transpose() << endl;
// 		cout << cen.transpose() << endl;
// 		cout << endl;

// // 	std::vector<size_t> indices(NSAMP);
// // 	boost::timer::cpu_timer lookup_time;
// // 	for(int i = 0; i < NSAMP; ++i)
// // 		indices[i] = bcc[samples[i]];
// // 	cout << N << " lookup: " << (double)NSAMP / lookup_time.elapsed().wall*1000000000.0 << " sec" << endl;

// // 	std::vector<V> centers(NSAMP);
// // 	boost::timer::cpu_timer getval_time;
// // 	for(int i = 0; i < NSAMP; ++i)
// // 		centers[i] = bcc[indices[i]];

// 	}


// }



template<int N, class F, class S> 
F 
test_bcc_children( int NSAMP ){
	typedef util::SimpleArray<N,F> V;
	typedef util::SimpleArray<N,S> I;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::uniform_real<> runif;
	boost::normal_distribution<> rnorm;	
	S const Nside = 5;
	BCC<N,F,S> bcc_parent(I(5),V(-(F)Nside),V((F)Nside));
	BCC<N,F,S> bcc(I(5),V(-(F)Nside/2.0),V((F)Nside/2.0));
	BOOST_VERIFY( bcc[bcc[V(0.0)]] == V(0.0) );
	BOOST_VERIFY( bcc.width_ == V(1.0) );
	S const i0 = bcc[V(0)];
	std::vector<size_t> nbrs,nbrs_we;		
	bcc.neighbors( i0, std::back_inserter(nbrs   ), false );
	bcc.neighbors( i0, std::back_inserter(nbrs_we), true );		

	/////////////////////////////////////////////////////
	// test neighbor coverage
	///////////////////////////////////////////////////////////////

	// F const RNapprox[5] = { 0.558099, 0.701687, 0.742306, 0.845359, 0.882879 };
	// F const Rapprox = RNapprox[N-3];
	F inrad = 0.5;
	if(N==3) inrad = 0.433015;

	int sum_in_nbrs=0, sum_in_nbrs_we=0;
	for(int i = 0; i < NSAMP; ++i){

		// pick random point in i0 parent cell
		V samp; for(int j = 0; j < N; ++j)
			samp[j] = (runif(rng)-0.5)*4.0*inrad;
		if( bcc_parent[samp] != i0 ){ --i; continue; }

		S index = bcc[samp];
		// if( index == i0 ){ --i; continue; }
		if( std::find( nbrs_we.begin(), nbrs_we.end(), index ) != nbrs_we.end() ){
			++sum_in_nbrs_we;
			if( std::find( nbrs.begin(), nbrs.end(), index ) != nbrs.end() ){
				++sum_in_nbrs;
			}
		}
	}

	int nc = std::pow(3,N);
	int nbccFC = 1+2*N+std::pow(2,N);
	int nbccFCE = nbccFC + N*(N-1)/2 * 4;
	printf("BCC child coverage: %i %3i %9.7f   %3i %9.7f\n",
		N, 
		nbccFC, (F)sum_in_nbrs/NSAMP,
		nbccFCE, (F)sum_in_nbrs_we/NSAMP
	);
	return 0;
}

TEST(bcc_lattice,children){
	cout << "BCC               DIM Nfc   frac_fc  Nfce  frac_fce" << std::endl;
	test_bcc_children<3,double,size_t>(100*1000);
	test_bcc_children<4,double,size_t>(100*1000);
	test_bcc_children<5,double,size_t>(100*1000);
	test_bcc_children<6,double,size_t>(100*1000);
	test_bcc_children<7,double,size_t>(100*1000);
}


}}}

