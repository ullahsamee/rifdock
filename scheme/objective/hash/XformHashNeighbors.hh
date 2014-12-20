#ifndef INCLUDED_objective_hash_XformHashNeighbors_HH
#define INCLUDED_objective_hash_XformHashNeighbors_HH

#include "scheme/util/SimpleArray.hh"
#include "scheme/numeric/rand_xform.hh"

#include <fstream>
#include "scheme/io/dump_pdb_atom.hh"

#include <sparsehash/sparse_hash_map>
#include <boost/random/mersenne_twister.hpp>
#include <set>
#include <map>

namespace scheme { namespace objective { namespace hash {

template< class XformHash >
struct XformHashNeighbors {
	typedef typename XformHash::Key Key;
	typedef typename XformHash::Float Float;
	typedef typename XformHash::Xform Xform;

	XformHash const xh_;
	// google::sparse_hash_map< Key, std::vector<Key> > ori_cache_;
	std::map< Key, std::vector<Key> > ori_cache_;
	double cart_bound_, ang_bound_, quat_bound_;
	int nsamp_;
	std::vector< util::SimpleArray<3,int16_t> > cart_shifts_;

	XformHashNeighbors( Float cart_bound, Float ang_bound, XformHash xh, double xcov_target=100.0 ) : xh_(xh) {
		cart_bound_ = cart_bound;
		ang_bound_ = ang_bound;
		quat_bound_ = numeric::deg2quat(ang_bound);
		Float ang_nside = 2.0 * quat_bound_ / xh_.ang_width();
		nsamp_ = ang_nside*ang_nside*ang_nside*3.0*xcov_target; // ~100 fold coverage
		std::cout << "XformHashNeighbors aw=" << xh_.ang_width() << ", quat_bound=" << quat_bound_ << " ang_nside=" << ang_nside << " nsamp=" << nsamp_ << std:: endl;
		fill_cart_shifts();
	}

	void fill_cart_shifts(){
		Float thresh = cart_bound_+xh_.cart_width();
		int16_t di = std::ceil( cart_bound_ / xh_.cart_width() ) ;
		for(int16_t i = -di; i <= di; ++i){
		for(int16_t j = -di; j <= di; ++j){
		for(int16_t k = -di; k <= di; ++k){					
			if( 
				Eigen::Vector3d( i-0.5, j-0.5, k-0.5 ).norm()*xh_.cart_width() < thresh || 
				Eigen::Vector3d( i    , j    , k     ).norm()*xh_.cart_width() < thresh || 
				Eigen::Vector3d( i+0.5, j+0.5, k+0.5 ).norm()*xh_.cart_width() < thresh ||
				false
			){
				cart_shifts_.push_back( util::SimpleArray<3,int16_t>(i,j,k) );
			}
		}}}
	}

	std::vector< util::SimpleArray<3,int16_t> > const & get_cart_shifts() const {
		return cart_shifts_;
	}

	std::vector<Key> const & get_ori_neighbors( Key key ){
		Key ori_key = key & XformHash::ORI_MASK;
		if( ori_cache_.find(ori_key) == ori_cache_.end() ){
			Key ksym, asym_key = xh_.asym_key( ori_key, ksym );
			// if( asym_key==ori_key ){
			if( true ){		
				// std::cout << "get_asym nbrs the hard way, store in " << ori_key << std::endl;
				boost::random::mt19937 rng((unsigned int)time(0) + 23058704);
				Xform c = xh_.get_center(key);
				// c.translation()[0] = c.translation()[1] = c.translation()[2] = 0;
				std::set<Key> keys;
				for(int i = 0; i < nsamp_; ++i){
					Xform p; numeric::rand_xform_quat(rng,p,cart_bound_,quat_bound_);
					p.translation()[0] = p.translation()[1] = p.translation()[2] = 0;
					Key nbkey = xh_.get_key( p * c ) & XformHash::ORI_MASK;
					// assert( ( nbkey & ~XformHash::ORI_MASK ) == 0 );
					keys.insert(nbkey);
				}
				ori_cache_.insert( std::make_pair( ori_key, std::vector<Key>(keys.size()) ) );
				std::copy(keys.begin(),keys.end(),ori_cache_[ori_key].begin());

				std::ofstream out("nbrs_asym.pdb");
				for(int i = 0; i < ori_cache_[ori_key].size(); ++i){
					Xform x = xh_.get_center( ori_cache_[ori_key][i] );
					Eigen::Quaterniond q(x.rotation());
					io::dump_pdb_atom( out, "C" ,i, 50*q.coeffs() );
				}
				out.close();

			} else {
				std::vector<Key> const & asym_nbrs = get_ori_neighbors( asym_key );
				// std::cout << "set sym_key nbrs from asym_nbrs " << asym_nbrs.size() << ", store in " << ori_key << " " << ori_cache_.size() << std::endl;
				ori_cache_.insert( std::make_pair( ori_key, std::vector<Key>( asym_nbrs.size() ) ) );
				for(int i = 0; i < asym_nbrs.size(); ++i){
					Key new_key = xh_.sym_key( asym_nbrs[i], ksym );
					ori_cache_[ori_key][i] = new_key;
					assert( new_key == (new_key&xh_.ORI_MASK) );
				}

				std::ofstream out("nbrs_sym.pdb");
				std::cout << "FLIP: " << ((ksym>>2)&1) << " " << ((ksym>>1)&1) << " " << ((ksym>>0)&1) << std::endl;
				for(int i = 0; i < asym_nbrs.size(); ++i){
					Xform x = xh_.get_center( ori_cache_[ori_key][i] );
					Eigen::Quaterniond q(x.rotation());
					Eigen::Vector4d v( q.coeffs() );
					// v *= q.coeffs()[3];
					// std::cout << q.w() << std::endl;
					io::dump_pdb_atom( out, "C" ,i, 50.0*v );
					// std::cout << std::showpos;
					// std::cout << Eigen::Quaterniond(xh_.get_center(ori_cache_[ori_key][i]).rotation()).coeffs().transpose() << std::endl
					//           << Eigen::Quaterniond(xh_.get_center(          asym_nbrs[i]).rotation()).coeffs().transpose() << std::endl;
				}
				out.close();

			}
			// for(int i = 0; i < ori_cache_[ori_key].size(); ++i){
			// 	std::cout << ori_cache_[ori_key][i] << std::endl;
			// }
			// std::cout << "added ori_cahce_ " << ori_key << " " << (float)nsamp_/keys.size() << std::endl;
		}
		return ori_cache_[ori_key];
	}
};



}}}

#endif
