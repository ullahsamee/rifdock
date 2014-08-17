#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <boost/timer/timer.hpp>

namespace scheme { namespace nest { namespace maps {

using namespace Eigen;
using std::cout;
using std::endl;

TEST(Trunc24CellMap,cell_lookup){
	double const r = sqrt(2)/2;
	double const h = 0.5;
	Matrix<double,48,4> t24;
	t24 <<  1, 0, 0, 0, // 0
	        0, 1, 0, 0, // 1
	        0, 0, 1, 0, // 2
	        0, 0, 0, 1, // 3
	       -1, 0, 0, 0, // 4
	        0,-1, 0, 0, // 5
	        0, 0,-1, 0, // 6
	        0, 0, 0,-1, // 7
	        h, h, h, h, // 8
	        h, h, h,-h, // 9
	        h, h,-h, h, // 10
	        h, h,-h,-h, // 11
	        h,-h, h, h, // 12
	        h,-h, h,-h, // 13
	        h,-h,-h, h, // 14
	        h,-h,-h,-h, // 15
	       -h, h, h, h, // 16
	       -h, h, h,-h, // 17
	       -h, h,-h, h, // 18
	       -h, h,-h,-h, // 19
	       -h,-h, h, h, // 20
	       -h,-h, h,-h, // 21
	       -h,-h,-h, h, // 22
	       -h,-h,-h,-h, // 23
	        r, r, 0, 0, // 24
	        r, 0, r, 0, // 25
	        r, 0, 0, r, // 26
	        0, r, r, 0, // 27
	        0, r, 0, r, // 28
	        0, 0, r, r, // 29
	       -r, r, 0, 0, // 30
	       -r, 0, r, 0, // 31
	       -r, 0, 0, r, // 32
	        0,-r, r, 0, // 33
	        0,-r, 0, r, // 34
	        0, 0,-r, r, // 35
	        r,-r, 0, 0, // 36
	        r, 0,-r, 0, // 37
	        r, 0, 0,-r, // 38
	        0, r,-r, 0, // 39
	        0, r, 0,-r, // 40
	        0, 0, r,-r, // 41
	       -r,-r, 0, 0, // 42
	       -r, 0,-r, 0, // 43
	       -r, 0, 0,-r, // 44
	        0,-r,-r, 0, // 45
	        0,-r, 0,-r, // 46
	        0, 0,-r,-r; // 47
	for(int i = 0; i < 48; ++i) ASSERT_DOUBLE_EQ( t24.row(i).norm(), 1.0 );

	boost::random::mt19937 mt((unsigned int)time(0));
	boost::normal_distribution<> rnorm;
	boost::uniform_real<> runif;

	int NSAMP = 10*1000*1000;
	std::vector<Vector4d> samp(NSAMP);
	std::vector<size_t> cell(NSAMP),cell2(NSAMP);

	for(int i = 0; i < NSAMP; ++i){
		Vector4d quat(rnorm(mt),rnorm(mt),rnorm(mt),rnorm(mt));
		samp[i] = quat.normalized();
	}
		// Vector4d const quat_pos = quat.cwiseAbs();

	boost::timer::cpu_timer naive;
	for(int i = 0; i < NSAMP; ++i){
		(t24*samp[i]).maxCoeff(&cell[i]);
	}
	cout << "naive rate:  " << (double)NSAMP / naive.elapsed().wall * 1000000000.0 << endl;

	boost::timer::cpu_timer clever;
	for(int i = 0; i < NSAMP; ++i){
		Vector4d const & quat(samp[i]);
		Vector4d const quat_pos = samp[i].cwiseAbs();
		Vector4d tmpv = quat_pos;

		size_t hyperface_axis; // closest hyperface-pair
		size_t edge_axis_1; // first axis of closest edge
		size_t edge_axis_2; // second axis of closest edge

		double const hyperface_dist = quat_pos.maxCoeff(&hyperface_axis); // dist to closest face
		double const corner_dist = quat_pos.sum()/2;
		tmpv[hyperface_axis] = -9e9;
		double const edge_dist = sqrt(2)/2 * ( hyperface_dist + tmpv.maxCoeff(&edge_axis_2) );

		edge_axis_1 = hyperface_axis<edge_axis_2 ? hyperface_axis : edge_axis_2;
		edge_axis_2 = hyperface_axis<edge_axis_2 ? edge_axis_2 : hyperface_axis;
		assert( edge_axis_1 < edge_axis_2 );

		// cell if closest if of form 1000 (add 4 if negative)
		size_t facecell = hyperface_axis | (quat[hyperface_axis]<0 ? 4 : 0);

		// cell if closest is of form 1111, bitwise by ( < 0)
		size_t const bit0 = quat[3]<0;
		size_t const bit1 = quat[2]<0;
		size_t const bit2 = quat[1]<0;
		size_t const bit3 = quat[0]<0;						
		size_t const cornercell = bit0 | bit1<<1 | bit2<<2 | bit3<<3;

		// cell if closest is of form 1100
		int perm_shift[3][4] = { {9,0,1,2}, {0,9,3,4}, {1,3,9,5} };
		size_t const sign_shift = (quat[edge_axis_1]<0)*1*6 + (quat[edge_axis_2]<0)*2*6;
		size_t const edgecell = sign_shift + perm_shift[edge_axis_1][edge_axis_2];

		// pick case 1000 1111 1100 without if statements
		size_t swtch; Vector3d(hyperface_dist,corner_dist,edge_dist).maxCoeff(&swtch);
		cell2[i] = swtch==0 ? facecell : (swtch==1 ? cornercell+8 : edgecell+24);
		// this is slower !?!
		// double const mx = std::max(std::max(hyperface_dist,corner_dist),edge_dist);
		// cell2[i] = hyperface_dist==mx ? facecell : (corner_dist==mx ? cornercell+8 : edgecell+24);


	}
	cout << "clever rate: " << (double)NSAMP / clever.elapsed().wall * 1000000000.0 << endl;

	for(int i = 0; i < NSAMP; ++i){
		ASSERT_DOUBLE_EQ( cell[i], cell2[i] );
	}

}

}}}
