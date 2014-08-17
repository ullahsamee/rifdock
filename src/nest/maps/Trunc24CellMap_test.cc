#include <gtest/gtest.h>

#include "numeric/util.hh"

#include <numeric/geom_4d.hh>

#include <Eigen/Dense>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <boost/timer/timer.hpp>

namespace scheme { namespace nest { namespace maps {

using namespace Eigen;
using std::cout;
using std::endl;


TEST(Trunc24CellMap,cell_lookup)
{
	typedef double Float;
	typedef uint64_t Index;
	typedef Eigen::Matrix<Float,4,1> V4;

	Map<Matrix<Float,48,4,RowMajor>const> t24( numeric::get_raw_48cell<Float>() );
	ASSERT_EQ( V4(0,0,0,-1), t24.row(7).transpose() );
	for(int i = 0; i < 48; ++i) ASSERT_FLOAT_EQ( t24.row(i).norm(), 1.0 );

	boost::random::mt19937 mt((unsigned int)time(0));
	boost::normal_distribution<> rnorm;
	boost::uniform_real<> runif;

	int NSAMP = 1*1000*1000;
	std::vector<V4> samp(NSAMP);
	std::vector<Index> cell(NSAMP),cell2(NSAMP);

	for(int i = 0; i < NSAMP; ++i){
		V4 quat(rnorm(mt),rnorm(mt),rnorm(mt),rnorm(mt));
		samp[i] = quat.normalized();
	}
		// V4 const quat_pos = quat.cwiseAbs();

	boost::timer::cpu_timer naive;
	for(int i = 0; i < NSAMP; ++i){
		(t24*samp[i]).maxCoeff(&cell[i]);
	}
	cout << "naive rate:  " << (Float)NSAMP / naive.elapsed().wall * 1000000000.0 << endl;

	boost::timer::cpu_timer clever;
	for(int i = 0; i < NSAMP; ++i){
		numeric::get_cell_48cell( samp[i], cell2[i] );
		// this is slower !?!
		// Float mx = std::max(std::max(hyperface_dist,corner_dist),edge_dist);
		// cell2[i] = hyperface_dist==mx ? facecell : (corner_dist==mx ? cornercell+8 : edgecell+24);


	}
	cout << "clever rate: " << (double)NSAMP / clever.elapsed().wall * 1000000000.0 << endl;

	for(int i = 0; i < NSAMP; ++i){
		if( cell[i] != cell2[i] ){
			ASSERT_FLOAT_EQ(
			     t24.row(cell [i]).dot(samp[i]) ,
			     t24.row(cell2[i]).dot(samp[i]) );
		} else { 
			ASSERT_EQ( cell[i], cell2[i] );
		}
	}

}

}}}
