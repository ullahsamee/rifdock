#include <gtest/gtest.h>

#include "scheme/numeric/bcc_lattice.hh"
#include "scheme/io/dump_pdb_atom.hh"
#include "scheme/nest/maps/TetracontoctachoronMap.hh"

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/foreach.hpp>
#include <boost/timer/timer.hpp>

#include <sparsehash/dense_hash_set>

#include <Eigen/Geometry>

namespace scheme { namespace numeric { namespace test {

using std::cout;
using std::endl;



// QuaternionMap Covrad
//  1                8  119.36447   61.22694    1.94954    3.83748    0.51790
//  2              128   59.04495   27.31148    2.16191    7.43172    0.73549
//  3             1040   27.85890   13.60802    2.04724    6.34244    0.73918
//  4             8480   13.59714    6.82994    1.99081    6.01270    0.76204
//  5            68416    6.91475    3.40731    2.02939    6.37995    0.76335
//  6           548496    3.44178    1.70610    2.01734    6.30746    0.76828
//  7          4392160    1.73778    0.85268    2.03801    6.50123    0.76802

void test_orientatin_coverage_4d( size_t Nside, int NSAMP ){
	typedef double F;
	typedef size_t S;
	typedef util::SimpleArray<4,F> V;
	typedef util::SimpleArray<4,S> I;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> rnorm;

	BCC<4,F,S> bcc(I(Nside),V(-1.0-2.0/Nside),V(1.0+2.0/Nside));
	// BCC<4,F,S> bcc(I(Nside),V( -2.0 ),V( 2.0 ) );
	// cout << bcc.lower_ << endl;
	// cout << bcc.width_ << endl;
	// cout << bcc.lower_+bcc.width_*Nside << endl;

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
	if(xcov < 50){
		cout << "volfrac will be low if xcov < 50 or so... (not enough samples)" << endl;
		cout << "idxfrac: " << (double)count / (Nside*Nside*Nside*Nside) << " xcov " << xcov << endl;
	}
	double volfrac = (double)count*(maxdiff*maxdiff*maxdiff) * 4.0/3.0*M_PI / 8.0/M_PI/M_PI;
	double avgfrac = (double)count*(avgdiff*avgdiff*avgdiff) * 4.0/3.0*M_PI / 8.0/M_PI/M_PI;
	// cout << "RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac      AvgFrac" << endl;
	printf("%2lu %16lu %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
				Nside, count, maxdiff*180.0/M_PI, avgdiff*180.0/M_PI, maxdiff/avgdiff, volfrac, avgfrac );

	// ASSERT_LE( volfrac , 3.5 );

}

void test_orientatin_coverage_3d_bt24( size_t Nside, int NSAMP ){
	typedef double F;
	typedef size_t S;
	typedef util::SimpleArray<4,F> V;
	typedef util::SimpleArray<4,S> I;
	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> rnorm;
	BCC<3,F,S> bcc(I(Nside),V( -1.0/Nside ),V( 1.0+1.0/Nside ));
	// BCC<3,F,S> bcc(I(Nside),V( 0.0 ),V( 1.0));	
	// BCC<3,F,S> bcc(I(Nside),V(-1.0 ),V( 2.0));

	nest::maps::TetracontoctachoronMap<> map;
	typedef nest::maps::TetracontoctachoronMap<3,V>::Params Params;
	// cout << bcc.lower_ << endl;
	// cout << bcc.width_ << endl;
	// cout << bcc.lower_+bcc.width_*Nside << endl;

	google::dense_hash_set<size_t> idx_seen;
	idx_seen.set_empty_key(999999999999);
	F maxdiff=0,avgdiff=0;
	for(int i = 0; i < NSAMP; ++i){
		V sample( rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng) );
		sample /= sample.norm();
		if(sample[0]<0) sample = -sample;
		ASSERT_DOUBLE_EQ( sample.norm(), 1.0 );

		Eigen::Quaterniond q;
		for(int i = 0; i < 4; ++i) q.coeffs()[i] = sample[i];
		Params params; size_t cell_index;
		map.value_to_params( q.matrix(), 0, params, cell_index );
		size_t index_bcc = bcc[params];
		idx_seen.insert(index_bcc);
		// size_t index = index_bcc + bcc.size()*cell_index;
		Params pcenter = bcc[index_bcc];
		Eigen::Matrix3d mcenter;
		map.params_to_value( pcenter, cell_index, 0, mcenter );
		Eigen::Quaterniond qcenter(mcenter);
		V center;
		for(int i = 0; i < 4; ++i) center[i] = qcenter.coeffs()[i];


		Eigen::Quaterniond samp(sample[0],sample[1],sample[2],sample[3]);
		Eigen::Quaterniond cen (center[0],center[1],center[2],center[3]);
		cen.normalize();
		maxdiff = std::max( maxdiff, samp.angularDistance(cen) );
		avgdiff += samp.angularDistance(cen);
	}
	avgdiff /= NSAMP;
	// size_t count = bcc.size() * 24;
	size_t count = 24 * idx_seen.size();
	double xcov = (double)NSAMP/count*24;
	if(xcov < 50){
		cout << "volfrac will be low if xcov < 50 or so... (not enough samples)" << endl;
		cout << "idxfrac: " << (double)count / (Nside*Nside*Nside) << " xcov " << xcov << endl;
	}
	double volfrac = (double)count*(maxdiff*maxdiff*maxdiff) * 4.0/3.0*M_PI / 8.0/M_PI/M_PI;
	double avgfrac = (double)count*(avgdiff*avgdiff*avgdiff) * 4.0/3.0*M_PI / 8.0/M_PI/M_PI;
	// cout << "RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac      AvgFrac" << endl;
	printf("%2lu %16lu %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
				Nside, count, maxdiff*180.0/M_PI, avgdiff*180.0/M_PI, maxdiff/avgdiff, volfrac, avgfrac );

	// ASSERT_LE( volfrac , 3.5 );

}


int const NSAMP = 1*1000*1000;

TEST(bcc_lattice,DISABLED_orientatin_coverage_4d){
	cout << "RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac      AvgFrac" << endl;
	for(int i = 2; i < 20; ++i){
		test_orientatin_coverage_4d(  i , NSAMP );
	}
	// 	test_orientatin_coverage_4d(   8 , NSAMP );
	// 	test_orientatin_coverage_4d(  16 , NSAMP );
	// 	test_orientatin_coverage_4d(  32 , NSAMP );
	// 	// test_orientatin_coverage_4d(  64 , NSAMP );
	// 	// test_orientatin_coverage_4d( 128 , NSAMP );
}

TEST(bcc_lattice,DISABLED_orientatin_coverage_3d_bt24){
	cout << "RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac      AvgFrac" << endl;
	for(int i = 2; i < 20; ++i){
		test_orientatin_coverage_3d_bt24(  i , NSAMP );
	}
	// test_orientatin_coverage_3d_bt24(  4 , NSAMP );
	// test_orientatin_coverage_3d_bt24(  8 , NSAMP );
	// test_orientatin_coverage_3d_bt24( 16 , NSAMP );
	// test_orientatin_coverage_3d_bt24( 32 , NSAMP );
	// test_orientatin_coverage_3d_bt24( 64 , NSAMP );
}

// conclusion, use 4D for covering > ~15°, and bt24 3D for < ~15°
// or figure out how to better optimize 48-cell based lookup, esp. at low grid number
// RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac      AvgFrac
//  2                1        nan        nan        nan        nan        nan
//  3               33   74.28860   48.00784    1.54743    3.81603    1.02987
//  4               85   54.29327   30.41631    1.78500    3.83697    0.67464
//  5              169   43.45715   25.11014    1.73066    3.91202    0.75468
//  6              325   33.33101   19.38110    1.71977    3.39437    0.66734
//  7              576   29.51199   16.11477    1.83136    4.17589    0.67987
//  8              877   24.50840   14.06207    1.74287    3.64146    0.68783
//  9             1310   21.81570   12.10485    1.80223    3.83627    0.65536
// 10             1915   19.12114   10.54955    1.81251    3.77608    0.63416
// 11             2749   16.93394    9.54341    1.77441    3.76513    0.67393
// 12             3283   15.63030    8.66045    1.80479    3.53593    0.60148
// 13             4451   14.02718    7.85187    1.78647    3.46497    0.60773
// 14             5581   13.09922    7.26046    1.80419    3.53817    0.60247
// 15             7241   11.91341    6.71108    1.77519    3.45333    0.61731
// 16             8755   11.14681    6.20759    1.79567    3.42010    0.59068
// 17            10921   10.36544    5.81305    1.78313    3.43049    0.60507
// 18            13046    9.84143    5.49671    1.79042    3.50738    0.61111
// 19            15131    9.12415    5.15155    1.77114    3.24172    0.58346
// [       OK ] bcc_lattice.orientatin_coverage_4d (5256 ms)
// [ RUN      ] bcc_lattice.orientatin_coverage_3d_bt24
// RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac      AvgFrac
//  2              216   52.18960   36.09369    1.44595    8.66037    2.86470
//  3              648   30.72092   18.29302    1.67938    5.29918    1.11883
//  4             1992   21.30586   12.53629    1.69953    5.43398    1.10695
//  5             3768   15.59875    9.39145    1.66095    4.03377    0.88032
//  6             6648   12.17794    7.43833    1.63719    3.38644    0.77170
//  7            11016   10.10831    6.16428    1.63982    3.20916    0.72778
//  8            17232    8.61811    5.23222    1.64712    3.11102    0.69619
//  9            23784    7.49805    4.54899    1.64829    2.82788    0.63148
// 10            33432    6.68782    4.01809    1.66443    2.82064    0.61172
// 11            43632    5.92284    3.59687    1.64667    2.55697    0.57268
// 12            57720    5.32703    3.25537    1.63638    2.46101    0.56164
// 13            73128    4.88102    2.97284    1.64187    2.39854    0.54191
// 14            93312    4.49436    2.73480    1.64340    2.38931    0.53832
// 15           115176    4.12394    2.53137    1.62913    2.27839    0.52694
// 16           140904    3.85766    2.35500    1.63807    2.28152    0.51907
// 17           170496    3.59167    2.20211    1.63101    2.22810    0.51353
// 18           204600    3.35699    2.06766    1.62357    2.18316    0.51012
// 19           242784    3.20765    1.94856    1.64616    2.26001    0.50663


}}}

