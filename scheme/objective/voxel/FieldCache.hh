#ifndef INCLUDED_objective_voxel_FieldCache_HH
#define INCLUDED_objective_voxel_FieldCache_HH

#include "scheme/objective/voxel/VoxelArray.hh"
#include "scheme/io/cache.hh"
// #include <boost/exception/all.hpp>
#include <exception>

namespace scheme { namespace objective { namespace voxel {

struct FieldException : public std::runtime_error {
	FieldException(std::string m) : std::runtime_error(m.c_str()) {}
};

template<class Float=float>
struct Field3D {
	virtual Float operator()(Float f, Float g, Float h) const = 0;
};

template<class Float=float>
struct FieldCache3D : VoxelArray<3,Float> {
	typedef VoxelArray<3,Float> BASE;
	typedef typename BASE::Bounds Float3;
	typedef typename VoxelArray<3,Float>::Indices Indices;

	std::string cache_loc_;

	FieldCache3D() : cache_loc_("") {}

	template<class F1,class F2, class F3>
	FieldCache3D(
		Field3D<Float> const & field,
		F1 const & lb,
		F2 const & ub,
		F3 const & cs,
		std::string const & cache_loc="",
		bool no_init = false
	) : BASE(lb,ub,cs), cache_loc_(cache_loc)
	{
		typename BASE::Indices extents;
		for(size_t i = 0; i < BASE::DIM; ++i) extents[i] = this->shape()[i];
		#ifdef CEREAL
			if( io::read_cache(cache_loc_,*this) ){
				bool extents_eq = true;
				for(size_t i = 0; i < BASE::DIM; ++i) extents_eq &= extents[i] == this->shape()[i];
				if( Float3(lb)==this->lb_ && Float3(ub)==this->ub_ && Float3(cs)==this->cs_ && extents_eq ){
					std::cout << "EXTENTS EQ, USING CACHE" << std::endl;
					check_against_field( field );
					return;
				}
				std::cout << "Warning, FieldCache3D bad cache " << cache_loc_ << " bounds mismatch, recomputing..." << std::endl;
				this->resize(extents);
			} // else {
		#endif
		// 	std::cout << "NO CACHE" << std::endl;
		// }
		if( !no_init ){
			// TODO: this should probable be based on integer indices...
			for(int i = 0; i < this->shape()[0]; ++i){
			for(int j = 0; j < this->shape()[1]; ++j){
			for(int k = 0; k < this->shape()[2]; ++k){							
				Float3 cen = this->indices_to_center( Indices(i,j,k) );
				this->operator[]( cen ) = field( cen[0], cen[1], cen[2] );
			}}}
		}
		#ifdef CEREAL
			io::write_cache(cache_loc_,*this);
		#endif

	}

	int check_against_field( Field3D<Float> const & field, bool permissive=false ) const {
		int nerror = 0;
		{
			// sanity check
			// std::exit(-1);
			// std::cout << this->shape()[0] << std::endl;
			// std::cout << this->shape()[1] << std::endl;
			// std::cout << this->shape()[2] << std::endl;						
			for(int i = 0; i < this->shape()[0]; i += 8){
			for(int j = 0; j < this->shape()[1]; j += 8){
			for(int k = 0; k < this->shape()[2]; k += 8){							
				Float3 cen = this->indices_to_center( Indices(i,j,k) );
				Float test1 = this->operator[]( cen );
				Float test2 = field( cen[0], cen[1], cen[2] );
				if( !(fabs(test1-test2) < 0.001 || fabs(1.0-test1/test2) < 0.001) && test2 < 9e6 ){
					std::cout << "FIELD MISMATCH stored: " << test1 << " recalculated: " << test2 << std::endl;
					++nerror;
					if(!permissive) throw FieldException("field check fails");
				} else {
					// std::cout << "check pass" << std::endl;
				}
			}}}
		}
		return nerror;
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
		#ifdef CEREAL
			if( io::read_cache(cache_loc_,*this) ){
				bool extents_eq = true;
				for(size_t i = 0; i < BASE::DIM; ++i) extents_eq &= extents[i] == this->shape()[i];
				if( Float3(lb)==this->lb_ && Float3(ub)==this->ub_ && Float3(cs)==this->cs_ && extents_eq ) return;
				std::cout << "Warning, BoundingFieldCache3D bad cache " << cache_loc_ << " bounds mismatch, recomputing..." << std::endl;
				this->resize(extents);
			}
		#endif
		// size_t ncalls = 0;
		for(Float f = this->lb_[0]+this->cs_[0]/2.0; f < this->ub_[0]+this->cs_[0]/2.0; f += this->cs_[0]){
		for(Float g = this->lb_[1]+this->cs_[1]/2.0; g < this->ub_[1]+this->cs_[1]/2.0; g += this->cs_[1]){
		for(Float h = this->lb_[2]+this->cs_[2]/2.0; h < this->ub_[2]+this->cs_[2]/2.0; h += this->cs_[2]){
			// ++ncalls;
			this->operator[]( Float3(f,g,h) ) = calc_val(Float3(f,g,h));
		}}}
		#ifdef CEREAL
			io::write_cache(cache_loc_,*this);
		#endif
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
