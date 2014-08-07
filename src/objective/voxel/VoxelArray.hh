#ifndef INCLUDED_objective_voxel_VoxelArray_HH
#define INCLUDED_objective_voxel_VoxelArray_HH

#include <boost/multi_array.hpp>
#include <boost/any.hpp>
#include "util/SimpleArray.hh"
#include <boost/type_traits.hpp>

namespace scheme { namespace objective { namespace voxel {



template<size_t _DIM,class _Float=float>
struct VoxelArray : boost::multi_array<_Float,_DIM> {
	BOOST_STATIC_ASSERT(_DIM>0);
	static size_t const DIM = _DIM;
	typedef _Float Float;
	typedef util::SimpleArray<DIM,size_t> Indices;
	util::SimpleArray<DIM,Float> lb_,ub_,cs_;

	VoxelArray(){ this->create(); }

	template<class F1,class F2, class F3>
	VoxelArray(F1 const & lb, F2 const & ub, F3 const & cs ) : lb_(lb),ub_(ub),cs_(cs) {
		Indices extents = floats_to_index(ub_);
		// std::cout << extents << std::endl;
		this->resize(extents+Indices(1)); // pad by one
	}

	template<class Floats> Indices floats_to_index(Floats const & f) const {
		Indices ind;
		for(int i = 0; i < DIM; ++i){
			Float tmp = ((f[i]-lb_[i])/cs_[i]);
			assert(tmp >= 0.0);
			ind[i] = tmp;
		}
		// std::cout << "floats_to_index " << f << " " << ind << std::endl;
		return ind;
	}

	template<class Floats>
	typename boost::disable_if< boost::is_arithmetic<Floats>, Float const & >::type
	operator[](Floats const & floats) const { return this->operator()(floats_to_index(floats)); }

	template<class Floats>
	typename boost::disable_if< boost::is_arithmetic<Floats>, Float & >::type
	operator[](Floats const & floats){ return this->operator()(floats_to_index(floats)); }	

};


template<class Float=float>
struct FieldCache3D : public VoxelArray<3,Float> {
	FieldCache3D(){}
	template<class F1,class F2, class F3, class Field>
	FieldCache3D(Field const & field, F1 const & lb, F2 const & ub, F3 const & cs) : VoxelArray<3,Float>(lb,ub,cs) {
		// size_t ncalls = 0;
		for(Float f = this->lb_[0]+this->cs_[0]/2.0; f < this->ub_[0]; f += this->cs_[0]){
		for(Float g = this->lb_[1]+this->cs_[1]/2.0; g < this->ub_[1]; g += this->cs_[1]){
		for(Float h = this->lb_[2]+this->cs_[2]/2.0; h < this->ub_[2]; h += this->cs_[2]){
			// ++ncalls;
			this->operator[]( util::SimpleArray<3,Float>(f,g,h) ) = field(f,g,h);
		}}}
		// std::cout << this->num_elements() << " " << ncalls << std::endl;
	}
};


}}}

#endif
