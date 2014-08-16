#ifndef INCLUDED_scheme_numeric_bcc_lattice_HH
#define INCLUDED_scheme_numeric_bcc_lattice_HH

#include "util/SimpleArray.hh"

namespace scheme { namespace numeric {

template<int DIM, class Float, class Index = size_t>
struct BCC {
	typedef util::SimpleArray<DIM,Index> Indices;
	typedef util::SimpleArray<DIM,Float> Floats;

	Indices sizes_, sizes_prefsum_;
	Floats lower_, width_, lower_cen_, half_width_;

	BCC(){}

	template<class Sizes>
	BCC(
		Sizes const & sizes,
		Floats lower=Floats(0), 
		Floats upper=Floats(1)
	) :
		sizes_(sizes),
		lower_(lower)
	{
		for(size_t i = 0; i < DIM; ++i)
			sizes_prefsum_[i] = sizes_.prod(i);
		width_ = (upper-lower_)/sizes_.template cast<Float>();
		half_width_ = width_ / 2.0;
		lower_cen_ = lower_ + half_width_;
	}

	Index
	size() const { 
		return sizes_.prod()*2;
	}

	Floats 
	operator[](
		Index index
	) const {
		bool odd = index & 1;
		Indices idx = ( (index>>1) / sizes_prefsum_ ) % sizes_;
		return lower_cen_ + width_ * idx.template cast<Float>() + (odd? half_width_ : 0);
	}

	Indices
	indices(
		Floats v,
		bool & odd
	) const {
		v = (v-lower_)/width_;
		Indices idx = v.template cast<Index>();
		v = v - idx.template cast<Float>() - 0.5;
		Indices idxc = idx - (v < 0);
		odd = (0.25 * DIM) < fabs( ( v.sign() * v ).sum() );
		return  odd ? idxc : idx;
	}

	Index
	operator[](
		Floats const & v
	) const {
		bool odd;
		Indices idx = indices(v,odd);
		Index index = (sizes_prefsum_*idx).sum();
		return (index<<1) + odd;
	}
};

}}

#endif
