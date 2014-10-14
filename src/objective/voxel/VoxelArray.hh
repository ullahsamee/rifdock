#ifndef INCLUDED_objective_voxel_VoxelArray_HH
#define INCLUDED_objective_voxel_VoxelArray_HH

#include <boost/multi_array.hpp>
#include "util/SimpleArray.hh"
#include <boost/type_traits.hpp>
#include <cereal/access.hpp>
//#include <boost/serialization/access.hpp>
//#include <boost/serialization/split_member.hpp>

namespace scheme { namespace objective { namespace voxel {



template<size_t _DIM,class _Float=float>
struct VoxelArray : boost::multi_array<_Float,_DIM> {
	BOOST_STATIC_ASSERT(_DIM>0);
	typedef boost::multi_array<_Float,_DIM> BASE;
	typedef VoxelArray<_DIM,_Float> THIS;
	static size_t const DIM = _DIM;
	typedef _Float Float;
	typedef util::SimpleArray<DIM,typename BASE::size_type> Indices;
	typedef util::SimpleArray<DIM,Float> Bounds;
	Bounds lb_,ub_,cs_;

	VoxelArray() {}

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

	// void write(std::ostream & out) const { 
	// 	out.write( (char const*)&lb_, sizeof(Bounds) );
	// 	out.write( (char const*)&ub_, sizeof(Bounds) );
	// 	out.write( (char const*)&cs_, sizeof(Bounds) );
	// 	for(size_t i = 0; i < DIM; ++i) out.write( (char const*)&(this->shape()[i]), sizeof() );
	// 	out.write( (char const*)this->data(), this->num_elements()*sizeof(Float) );
	// }
	// void read(std::istream & in){ 
	// 	in.read( (char*)&lb_, sizeof(Bounds) );
	// 	in.read( (char*)&ub_, sizeof(Bounds) );
	// 	in.read( (char*)&cs_, sizeof(Bounds) );
	// 	in.read( (char*)this->data(), this->num_elements()*sizeof(Float) );
	// }
	bool operator==(THIS const & o) const { 
		return lb_==o.lb_ && ub_==o.ub_ && cs_==o.cs_ && (BASE const &)o == (BASE const &)*this;
	}

	// friend class boost::serialization::access;
    friend class cereal::access; // befriend the cereal version of access

    template<class Archive> void save(Archive & ar, const unsigned int ) const {
        ar & lb_;
        ar & ub_;
        ar & cs_;
        for(size_t i = 0; i < DIM; ++i) ar & this->shape()[i];
        for(size_t i = 0; i < this->num_elements(); ++i) ar & this->data()[i];
    }
    template<class Archive> void load(Archive & ar, const unsigned int ){
        ar & lb_;
        ar & ub_;
        ar & cs_;
        Indices extents;
        for(size_t i = 0; i < DIM; ++i) ar & extents[i];
        this->resize(extents);
        for(size_t i = 0; i < this->num_elements(); ++i) ar & this->data()[i];
    }
    // BOOST_SERIALIZATION_SPLIT_MEMBER()

};



}}}

#endif
