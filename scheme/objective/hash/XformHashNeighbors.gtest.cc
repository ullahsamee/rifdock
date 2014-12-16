#include <gtest/gtest.h>

#include "scheme/objective/hash/XformHash.hh"
#include "scheme/objective/hash/XformHashNeighbors.hh"
#include "scheme/numeric/rand_xform.hh"

#include <Eigen/Geometry>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>

namespace scheme { namespace objective { namespace hash { namespace xhnbtest {

using std::cout;
using std::endl;

typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;

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
		numeric::rand_xform(rng,xrand,cart_bound,quat_bound);
		nbrs.insert( xh.get_key( xrand*x ) );
	}
}


TEST( XformHashNeighbors, Quat_BCC7_Zorder_check_rotations ){
	int NSAMP = 10;
	#ifdef NDEBUG
	NSAMP = 100
	#endif
	boost::random::mt19937 rng((unsigned int)time(0) + 94586049);

	XformHash_Quat_BCC7_Zorder<Xform> xh(1.0,10.0);
	cout << xh.approx_nori()/98.0 * 3.1 << endl;
	XformHashNeighbors< XformHash_Quat_BCC7_Zorder<Xform> > nb(3.0,30.0,xh);

	for(int i = 0; i < NSAMP; ++i){
		Xform x; numeric::rand_xform(rng,x);
		nb.get_ori_neighbors( xh.get_key(x) );
	}

}


}}}}
