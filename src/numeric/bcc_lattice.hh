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
		Indices indices = ( (index>>1) / sizes_prefsum_ ) % sizes_;
		return lower_cen_ + width_ * indices.template cast<Float>() + (odd? half_width_ : 0);
	}

	Indices
	indices(
		Floats value,
		bool & odd
	) const {
		value = (value-lower_)/width_;
		Indices const indices = value.template cast<Index>();
		value = value - indices.template cast<Float>() - 0.5;
		Indices const corner_indices = indices - (value < 0);
		odd = (0.25 * DIM) < fabs( ( value.sign() * value ).sum() );
		return  odd ? corner_indices : indices;
	}

	Index
	operator[](
		Floats const & value
	) const {
		bool odd;
		Indices indices = indices(value,odd);
		Index index = (sizes_prefsum_*indices).sum();
		return (index<<1) + odd;
	}
};

}}

#endif
