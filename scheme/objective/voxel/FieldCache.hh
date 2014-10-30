#ifndef INCLUDED_objective_voxel_FieldCache_HH
#define INCLUDED_objective_voxel_FieldCache_HH

#include "scheme/objective/voxel/VoxelArray.hh"
#include "scheme/io/cache.hh"
// #include <boost/exception/all.hpp>

namespace scheme { namespace objective { namespace voxel {

template<class Float=float>
struct Field3D {
	virtual Float operator()(Float f, Float g, Float h) const = 0;
};

template<class Float=float>
struct FieldCache3D : VoxelArray<3,Float> {
	typedef VoxelArray<3,Float> BASE;
	typedef typename BASE::Bounds Float3;

	std::string cache_loc_;

	FieldCache3D() : cache_loc_("") {}

	template<class F1,class F2, class F3>
	FieldCache3D(
		Field3D<Float> const & field,
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
			if( Float3(lb)==this->lb_ && Float3(ub)==this->ub_ && Float3(cs)==this->cs_ && extents_eq ){
				std::cout << "EXTENTS EQ, USING CACHE" << std::endl;
				return;
			}
			std::cout << "Warning, FieldCache3D bad cache " << cache_loc_ << " bounds mismatch, recomputing..." << std::endl;
			this->resize(extents);
		} // else {
		// 	std::cout << "NO CACHE" << std::endl;
		// }
		// TODO: this should probable be based on integer indices...
		for(Float f = this->lb_[0]+this->cs_[0]/2.0; f < this->ub_[0]+this->cs_[0]/2.0; f += this->cs_[0]){
		for(Float g = this->lb_[1]+this->cs_[1]/2.0; g < this->ub_[1]+this->cs_[1]/2.0; g += this->cs_[1]){
		for(Float h = this->lb_[2]+this->cs_[2]/2.0; h < this->ub_[2]+this->cs_[2]/2.0; h += this->cs_[2]){
			this->operator[]( Float3(f,g,h) ) = field(f,g,h);
		}}}
		io::write_cache(cache_loc_,*this);
	}

};

template<class Float> struct AggMin {
	static Float initval() { return 9e9; }
	static Float aggregate( Float agg, Float newval ) { return std::min(agg,newval); }
};
template<class Float> struct AggMax {
	static Float initval() { return -9e9; }
	static Float aggregate( Float agg, Float newval ) { return std::max(agg,newval); }
};

template<class Float=float,template<class> class AGG = AggMin>
struct BoundingFieldCache3D : VoxelArray<3,Float> {
	typedef VoxelArray<3,Float> BASE;
	typedef AGG<Float> Aggregator;
	typedef typename BASE::Bounds Float3;
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
		if( ref_.cs_.norm() > spread_ )
			std::cout << "BoundingFieldCache3D warning: spread " << spread_ 
		              << " less than ref cell_size.norm() " << ref_.cs_.norm() << std::endl;
		Float3 lb = ref.lb_-spread;
		Float3 ub = ref.ub_+spread;
		typename BASE::Indices extents;
		for(size_t i = 0; i < BASE::DIM; ++i) extents[i] = this->shape()[i];
		if( io::read_cache(cache_loc_,*this) ){
			bool extents_eq = true;
			for(size_t i = 0; i < BASE::DIM; ++i) extents_eq &= extents[i] == this->shape()[i];
			if( Float3(lb)==this->lb_ && Float3(ub)==this->ub_ && Float3(cs)==this->cs_ && extents_eq ) return;
			std::cout << "Warning, BoundingFieldCache3D bad cache " << cache_loc_ << " bounds mismatch, recomputing..." << std::endl;
			this->resize(extents);
		}
		// size_t ncalls = 0;
		for(Float f = this->lb_[0]+this->cs_[0]/2.0; f < this->ub_[0]+this->cs_[0]/2.0; f += this->cs_[0]){
		for(Float g = this->lb_[1]+this->cs_[1]/2.0; g < this->ub_[1]+this->cs_[1]/2.0; g += this->cs_[1]){
		for(Float h = this->lb_[2]+this->cs_[2]/2.0; h < this->ub_[2]+this->cs_[2]/2.0; h += this->cs_[2]){
			// ++ncalls;
			this->operator[]( Float3(f,g,h) ) = calc_val(Float3(f,g,h));
		}}}
		io::write_cache(cache_loc_,*this);
	}
	Float calc_val(Float3 const & f3) const {
		Float3 beg = ref_.lb_.max(f3-spread_);
		Float3 end = ref_.ub_.min(f3+spread_);
		Float3 idx;
		Float val = Aggregator::initval();
		for(idx[0] = beg[0]; idx[0] <= end[0]; idx[0] += ref_.cs_[0]){
		for(idx[1] = beg[1]; idx[1] <= end[1]; idx[1] += ref_.cs_[1]){
		for(idx[2] = beg[2]; idx[2] <= end[2]; idx[2] += ref_.cs_[2]){
			if( (f3-idx).squaredNorm() <= spread_*spread_ )
				val = Aggregator::aggregate( val, ref_[idx] );
		}}}
		// std::cout << f3 << " | " << beg << " => " << end << std::endl;
		return val;
	}
};


}}}

#endif
