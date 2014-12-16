#ifndef INCLUDED_objective_hash_XformHashNeighbors_HH
#define INCLUDED_objective_hash_XformHashNeighbors_HH

#include "scheme/util/SimpleArray.hh"
#include "scheme/numeric/rand_xform.hh"

#include <sparsehash/sparse_hash_map>
#include <boost/random/mersenne_twister.hpp>
#include <set>

namespace scheme { namespace objective { namespace hash {

template< class XformHash >
struct XformHashNeighbors {
	typedef typename XformHash::Key Key;
	typedef typename XformHash::Float Float;
	typedef typename XformHash::Xform Xform;

	XformHash xh_;
	google::sparse_hash_map< Key, std::vector<Key> > ori_cache_;
	double cart_bound_, ang_bound_, quat_bound_;
	int nsamp_;

	XformHashNeighbors( Float cart_bound, Float ang_bound, XformHash xh ) : xh_(xh) {
		cart_bound_ = cart_bound;
		ang_bound_ = ang_bound;
		quat_bound_ = numeric::deg2quat(ang_bound);
		Float ang_nside = 2.0 * quat_bound_ / xh_.ang_width();
		nsamp_ = ang_nside*ang_nside*ang_nside*300; // ~100 fold coverage
		std::cout << quat_bound_ << " " << ang_nside << " " << nsamp_ << std:: endl;
	}

	std::vector<typename XformHash::Key> const & get_ori_neighbors( Key key ){
		Key ori_key = key & XformHash::ORI_MASK;
		if( ori_cache_.find(ori_key) == ori_cache_.end() ){
			boost::random::mt19937 rng((unsigned int)time(0) + 23058704);
			Xform c = xh_.get_center(key);
			// c.translation()[0] = c.translation()[1] = c.translation()[2] = 0;
			std::set<Key> keys;
			for(int i = 0; i < nsamp_; ++i){
				Xform x;
				numeric::rand_xform(rng,x,cart_bound_,quat_bound_);
				x.translation()[0] = x.translation()[1] = x.translation()[2] = 0;
				Key nbkey = xh_.get_key( x * c ) & XformHash::ORI_MASK;
				// assert( ( nbkey & ~XformHash::ORI_MASK ) == 0 );
				keys.insert(nbkey);
			}
			ori_cache_[ori_key] = std::vector<Key>(keys.size());
			std::copy(keys.begin(),keys.end(),ori_cache_[ori_key].begin());
			// std::cout << "added ori_cahce_ " << ori_key << " " << (float)nsamp_/keys.size() << std::endl;
		}
		return ori_cache_[ori_key];
	}
};



}}}

#endif
