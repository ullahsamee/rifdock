#ifndef INCLUDED_objective_hash_XformMap_HH
#define INCLUDED_objective_hash_XformMap_HH

#include "scheme/util/SimpleArray.hh"
#include "scheme/util/dilated_int.hh"
#include "scheme/numeric/FixedPoint.hh"
#include "scheme/nest/maps/TetracontoctachoronMap.hh"
#include "scheme/numeric/bcc_lattice.hh"
#include "scheme/objective/hash/XformHash.hh"

#include <sparsehash/dense_hash_map>

namespace scheme { namespace objective { namespace hash {






template<
	class Xform,
	class Val=numeric::FixedPoint<-17>,
	int ArrayBits=8,
	// template<class X> class _Hasher = XformHash_bt24_BCC3_Zorder >
	template<class X> class _Hasher = XformHash_bt24_Cubic_Zorder >	
struct XformMap : _Hasher<Xform> {
	typedef _Hasher<Xform> Hasher;
	typedef uint64_t Key;
	typedef typename Xform::Scalar Float;
    typedef util::SimpleArray< (1<<ArrayBits), Val >  ValArray;
    typedef google::dense_hash_map<Key,ValArray> Hash;

	XformMap( Float cart_resl, Float ang_resl, Float cart_bound=512.0 ) :
		Hasher( cart_resl, ang_resl, cart_bound )
	{
	}

};


}}}

#endif
