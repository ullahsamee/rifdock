#ifndef INCLUDED_objective_hash_XformMap_HH
#define INCLUDED_objective_hash_XformMap_HH

#include "scheme/util/SimpleArray.hh"
#include "scheme/util/dilated_int.hh"
#include "scheme/numeric/FixedPoint.hh"
#include "scheme/nest/maps/TetracontoctachoronMap.hh"
#include "scheme/numeric/bcc_lattice.hh"

#include <sparsehash/dense_hash_map>

namespace scheme { namespace objective { namespace hash {


template<
	class Xform,
	class Val=numeric::FixedPoint<-17>,
	int ArrayBits=8 >
struct XformMap {
	typedef uint64_t Key;
	typedef typename Xform::Scalar Float;
    typedef util::SimpleArray< (1<<ArrayBits), Val >  ValArray;
    typedef google::dense_hash_map<Key,ValArray> Hash;
	typedef scheme::nest::maps::TetracontoctachoronMap<> OriMap;
	typedef scheme::numeric::BCC< 3, Float, uint64_t > BCC;
	typedef scheme::util::SimpleArray<3,Float> F3;
	typedef scheme::util::SimpleArray<3,uint64_t> I3;	

	Float grid_size_;
	Float grid_spacing_;
	OriMap ori_map_;
	BCC cart_bcc_, ori_bcc_;

	XformMap( Float cart_resl, Float ang_resl, Float cart_bound=512.0 )
	{
		cart_resl /= 0.56; // TODO: fix this number!
		float covrad[64] = { 49.66580,25.99805,17.48845,13.15078,10.48384, 8.76800, 7.48210, 6.56491, 5.84498, 5.27430, 4.78793, 4.35932,
		                      4.04326, 3.76735, 3.51456, 3.29493, 3.09656, 2.92407, 2.75865, 2.62890, 2.51173, 2.39665, 2.28840, 2.19235,
		                      2.09949, 2.01564, 1.94154, 1.87351, 1.80926, 1.75516, 1.69866, 1.64672, 1.59025, 1.54589, 1.50077, 1.46216,
		                      1.41758, 1.38146, 1.35363, 1.31630, 1.28212, 1.24864, 1.21919, 1.20169, 1.17003, 1.14951, 1.11853, 1.09436,
		                      1.07381, 1.05223, 1.02896, 1.00747, 0.99457, 0.97719, 0.95703, 0.93588, 0.92061, 0.90475, 0.89253, 0.87480,
		                      0.86141, 0.84846, 0.83677, 0.82164 };
		uint64_t ori_nside = 1;
		while( covrad[ori_nside-1] > ang_resl && ori_nside < 62 ) ++ori_nside;
		// std::cout << "requested ang_resl: " << ang_resl << " got " << covrad[ori_nside-1] << std::endl;
		if( 2*(int)(cart_bound/cart_resl) > 8192 ){
			throw std::out_of_range("can have at most 8192 cart cells!");
		}
		cart_bcc_.init(  I3(2.0*cart_bound/cart_resl), F3(-cart_bound), F3(cart_bound) );
		ori_bcc_.init(  I3(ori_nside+2), F3(-1.0/ori_nside), F3(1.0+1.0/ori_nside) );

		// std::cout << "cart_bcc " << cart_bcc_ << std::endl;
		// std::cout << " ori_bcc " <<  ori_bcc_ << std::endl;
		// std::cout << log2( (double)cart_bcc_.size() / ori_bcc_.size() ) << std::endl;
	}

	// Key structure
	// 5 bits bt24 cell index
	// 7 bits high order cart X bits
	// 7 bits high order cart Y bits
	// 7 bits high order cart Z bits		
	// 36 bits 6*6 zorder cart/ori
	// 2 bits cart/ori even/odd
	Key get_key( Xform const & x ) const {

		uint64_t cell_index;
		F3 params;
		ori_map_.value_to_params( x.rotation(), 0, params, cell_index );
		assert( cell_index < 24 );
		assert( 0.0 <= params[0] && params[0] <= 1.0 );
		assert( 0.0 <= params[1] && params[1] <= 1.0 );
		assert( 0.0 <= params[2] && params[2] <= 1.0 );
		
		bool ori_odd, cart_odd;
		I3 ori_indices, cart_indices;
		 ori_indices =  ori_bcc_.get_indices( params,  ori_odd );
		F3 trans(x.translation());
		cart_indices = cart_bcc_.get_indices( trans, cart_odd );

		// std::cout << "get_index  " << cell_index << " " << cart_indices << " " << cart_odd << " " << ori_indices << " " << ori_odd << std::endl;

		Key key;

		key = cell_index << 59;

		key = key | (cart_indices[0]>>6) << 52;
		key = key | (cart_indices[1]>>6) << 45;
		key = key | (cart_indices[2]>>6) << 38;

		// 6*6 zorder
		key = key>>2;
		key = key | ( util::dilate<6>(   ori_indices[0]       ) << 0 );
		key = key | ( util::dilate<6>(   ori_indices[1]       ) << 1 );
		key = key | ( util::dilate<6>(   ori_indices[2]       ) << 2 );
		key = key | ( util::dilate<6>( (cart_indices[0] & 63) ) << 3 );
		key = key | ( util::dilate<6>( (cart_indices[1] & 63) ) << 4 );
		key = key | ( util::dilate<6>( (cart_indices[2] & 63) ) << 5 );
		key = key<<2;

		// lowest two bits, even/odd
		key = key | ori_odd | cart_odd<<1;

		return key;								
	}

	Xform get_center(Key key) const {
		I3 cart_indices, ori_indices;

		uint64_t cell_index = key >> 59;

		cart_indices[0] = (((key>>52)&127) << 6) | (util::undilate<6>( (key>>5) ) & 63);
		cart_indices[1] = (((key>>45)&127) << 6) | (util::undilate<6>( (key>>6) ) & 63);
		cart_indices[2] = (((key>>38)&127) << 6) | (util::undilate<6>( (key>>7) ) & 63);			

		ori_indices[0] = util::undilate<6>( (key>>2)&((1ull<<36)-1) ) & 63;
		ori_indices[1] = util::undilate<6>( (key>>3)&((1ull<<36)-1) ) & 63;
		ori_indices[2] = util::undilate<6>( (key>>4)&((1ull<<36)-1) ) & 63;

		bool  ori_odd = key & 1;
		bool cart_odd = key & 2;

		F3 trans = cart_bcc_.get_center(cart_indices,cart_odd);
		F3 params = ori_bcc_.get_center(ori_indices,ori_odd);
		Eigen::Matrix3d m;
		ori_map_.params_to_value( params, cell_index, 0, m );
		// std::cout << "get_center " << cell_index << " " << cart_indices << " " << cart_odd << " " << ori_indices << " " << ori_odd << std::endl;

		Xform center( m );
		center.translation()[0] = trans[0];
		center.translation()[1] = trans[1];
		center.translation()[2] = trans[2];				

		return center;
	}

	Key approx_size() const { return ori_bcc_.size() * cart_bcc_.size() * 24; }

};


}}}

#endif
