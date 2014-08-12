#ifndef INCLUDED_objective_voxel_FieldCache_HH
#define INCLUDED_objective_voxel_FieldCache_HH

#include "objective/voxel/VoxelArray.hh"
#include "io/cache.hh"
#include <boost/exception/all.hpp>

namespace scheme { namespace objective { namespace voxel {



template<class Float=float>
struct FieldCache3D : VoxelArray<3,Float> {
	typedef VoxelArray<3,Float> BASE;
	typedef typename BASE::Bounds Bounds;

	std::string cache_loc_;

	FieldCache3D() : cache_loc_("") {}

	template<class F1,class F2, class F3, class Field>
	FieldCache3D(
		Field const & field,
		F1 const & lb,
		F2 const & ub,
		F3 const & cs,
		std::string const & cache_loc=""
	) : BASE(lb,ub,cs), cache_loc_(cache_loc) {
		typename BASE::Indices extents;
		for(size_t i = 0; i < BASE::DIM; ++i) extents[i] = this->shape()[i];
		if( io::read_cache(cache_loc_,*this) ){
			bool extents_eq = true;
			for(size_t i = 0; i < BASE::DIM; ++i) extents_eq &= extents[i] == this->shape()[i];
			if( Bounds(lb)==this->lb_ && Bounds(ub)==this->ub_ && Bounds(cs)==this->cs_ && extents_eq ) return;
			std::cout << "Warning, FieldCache3D bad cache " << cache_loc_ << " bounds mismatch, recomputing..." << std::endl;
			this->resize(extents);
		}
		// size_t ncalls = 0;
		for(Float f = this->lb_[0]+this->cs_[0]/2.0; f < this->ub_[0]; f += this->cs_[0]){
		for(Float g = this->lb_[1]+this->cs_[1]/2.0; g < this->ub_[1]; g += this->cs_[1]){
		for(Float h = this->lb_[2]+this->cs_[2]/2.0; h < this->ub_[2]; h += this->cs_[2]){
			// ++ncalls;
			this->operator[]( util::SimpleArray<3,Float>(f,g,h) ) = field(f,g,h);
		}}}
		io::write_cache(cache_loc_,*this);

		// std::cout << this->num_elements() << " " << ncalls << std::endl;
	}

};

template<class Float=float>
struct BoundingFieldCache3D : VoxelArray<3,Float> {
	typedef VoxelArray<3,Float> BASE;
	typedef typename BASE::Bounds Bounds;
	FieldCache3D<Float> const & ref_;
	Float spread_;
	std::string cache_loc_;
	template<class F>
	BoundingFieldCache3D(
		FieldCache3D<Float> const & ref,
		Float spread,
		F const & cs,
		std::string cache_loc=""
	) : BASE(ref.lb_-spread,ref.ub_+spread,cs),
	    ref_(ref),
	    spread_(spread),
	    cache_loc_(cache_loc)
	{
		Bounds lb = ref.lb_-spread;
		Bounds ub = ref.ub_+spread;
		typename BASE::Indices extents;
		for(size_t i = 0; i < BASE::DIM; ++i) extents[i] = this->shape()[i];
		if( io::read_cache(cache_loc_,*this) ){
			bool extents_eq = true;
			for(size_t i = 0; i < BASE::DIM; ++i) extents_eq &= extents[i] == this->shape()[i];
			if( Bounds(lb)==this->lb_ && Bounds(ub)==this->ub_ && Bounds(cs)==this->cs_ && extents_eq ) return;
			std::cout << "Warning, BoundingFieldCache3D bad cache " << cache_loc_ << " bounds mismatch, recomputing..." << std::endl;
			this->resize(extents);
		}
		// size_t ncalls = 0;
		for(Float f = this->lb_[0]+this->cs_[0]/2.0; f < this->ub_[0]; f += this->cs_[0]){
		for(Float g = this->lb_[1]+this->cs_[1]/2.0; g < this->ub_[1]; g += this->cs_[1]){
		for(Float h = this->lb_[2]+this->cs_[2]/2.0; h < this->ub_[2]; h += this->cs_[2]){
			// ++ncalls;
			this->operator[]( util::SimpleArray<3,Float>(f,g,h) ) = calc_val(f,g,h);
		}}}
		io::write_cache(cache_loc_,*this);
	}
	Float calc_val(Float f, Float g, Float h) const {
		return 0.0;//ref_(Bounds(f,g,h));
	}
};


}}}

#endif
