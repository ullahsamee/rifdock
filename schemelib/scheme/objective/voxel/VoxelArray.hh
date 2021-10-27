#ifndef INCLUDED_objective_voxel_VoxelArray_HH
#define INCLUDED_objective_voxel_VoxelArray_HH

#include <boost/multi_array.hpp>
#include "scheme/util/SimpleArray.hh"
#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <scheme/util/assert.hh>

#include <boost/format.hpp>
#include <fstream>

#include <random>
#include<boost/random/uniform_real.hpp>

#ifdef CEREAL
#include <cereal/access.hpp>
#endif
//#include <boost/serialization/access.hpp>
//#include <boost/serialization/split_member.hpp>

namespace scheme { namespace objective { namespace voxel {



template< size_t _DIM, class _Float=float, class _Value=float >
struct VoxelArray : boost::multi_array<_Value,_DIM> {
	BOOST_STATIC_ASSERT((_DIM>0));
	typedef boost::multi_array<_Value,_DIM> BASE;
	typedef VoxelArray<_DIM,_Float,_Value> THIS;
	static size_t const DIM = _DIM;
	typedef _Float Float;
    typedef _Value Value;
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
			// assert(tmp >= 0.0);
			ind[i] = tmp;
		}
		// std::cout << "floats_to_index " << f << " " << ind << std::endl;
		return ind;
	}

	Bounds indices_to_center( Indices const & idx ) const {
		Bounds c;
		for(int i = 0; i < DIM; ++i){
			c[i] = (idx[i]+0.5)*cs_[i] + lb_[i];
		}
		return c;
	}

	template<class Floats>
	typename boost::disable_if< boost::is_arithmetic<Floats>, Value const & >::type
	operator[](Floats const & floats) const { return this->operator()(floats_to_index(floats)); }

	template<class Floats>
	typename boost::disable_if< boost::is_arithmetic<Floats>, Value & >::type
	operator[](Floats const & floats){ return this->operator()(floats_to_index(floats)); }

	Value at( Float f, Float g, Float h ) const {
		Indices idx = floats_to_index( Bounds( f, g, h ) );
		if( idx[0] < this->shape()[0] && idx[1] < this->shape()[1] && idx[2] < this->shape()[2] )
			return this->operator()(idx);
		else return Value(0);
	}

	template<class V>
	Value at( V const & v ) const {
		Indices idx = floats_to_index( Bounds( v[0], v[1], v[2] ) );
		if( idx[0] < this->shape()[0] && idx[1] < this->shape()[1] && idx[2] < this->shape()[2] )
			return this->operator()(idx);
		else return Value(0);
	}

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

	#ifdef CEREAL
	// friend class boost::serialization::access;
    friend class cereal::access; // befriend the cereal version of access
    #endif

    template<class Archive> void save(Archive & ar, const unsigned int ) const {
    	BOOST_VERIFY( boost::is_pod<Float>::type::value );
        ar & lb_;
        ar & ub_;
        ar & cs_;
        for(size_t i = 0; i < DIM; ++i) ar & this->shape()[i];
        for(size_t i = 0; i < this->num_elements(); ++i) ar & this->data()[i];
    }
    template<class Archive> void load(Archive & ar, const unsigned int ){
    	BOOST_VERIFY( boost::is_pod<Float>::type::value );
        ar & lb_;
        ar & ub_;
        ar & cs_;
        Indices extents;
        for(size_t i = 0; i < DIM; ++i) ar & extents[i];
        this->resize(extents);
        for(size_t i = 0; i < this->num_elements(); ++i) ar & this->data()[i];
    }
    void save( std::ostream & out ) const {
    	BOOST_VERIFY( boost::is_pod<Float>::type::value );
  		out.write( (char*)&lb_, sizeof(Bounds) );
  		out.write( (char*)&ub_, sizeof(Bounds) );
  		out.write( (char*)&cs_, sizeof(Bounds) );
        for(size_t i = 0; i < DIM; ++i){
        	out.write( (char*)(&(this->shape()[i])), sizeof(typename BASE::size_type) );
        }
        for(size_t i = 0; i < this->num_elements(); ++i) out.write( (char*)(&(this->data()[i])), sizeof(Float) );
    }
  	void load( std::istream & in ){
    	BOOST_VERIFY( boost::is_pod<Float>::type::value );
  		ALWAYS_ASSERT( in.good() );
  		in.read( (char*)&lb_, sizeof(Bounds) );
  		ALWAYS_ASSERT( in.good() );
  		in.read( (char*)&ub_, sizeof(Bounds) );
  		ALWAYS_ASSERT( in.good() );
  		in.read( (char*)&cs_, sizeof(Bounds) );
  		ALWAYS_ASSERT( in.good() );
  		// todo: should I check these against the c'tor values? if not, should add default ctor?
        Indices extents;
        for(size_t i = 0; i < DIM; ++i){
        	in.read( (char*)(&(extents[i])), sizeof(typename BASE::size_type) );
        }
  		ALWAYS_ASSERT( in.good() );
        this->resize(extents);
        for(size_t i = 0; i < this->num_elements(); ++i) in.read( (char*)(&(this->data()[i])), sizeof(Float) );
  		ALWAYS_ASSERT( in.good() );
    }
    // BOOST_SERIALIZATION_SPLIT_MEMBER()



    void dump_pdb( std::string const & fname, Value threshold_value, bool higher_than_threshold, float fraction=1.0 ) {
        std::ofstream f( fname );
        if ( ! f ) {
            std::cerr << "Couldn't open " <<  fname << std::endl;
            ALWAYS_ASSERT( false );
        }

        std::mt19937 rng(time(0));
        boost::uniform_real<> uniform;

        size_t ires = 1;
        size_t iatom = 1;
        for ( size_t ix = 0; ix < this->shape()[0]; ix++ ) {
            for ( size_t iy = 0; iy < this->shape()[1]; iy++ ) {
                for ( size_t iz = 0; iz < this->shape()[0]; iz++ ) {
                    Indices idx(ix, iy, iz);
                    Value val = this->operator()(idx);
                    if ( (higher_than_threshold && val >= threshold_value ) ||
                        (!higher_than_threshold && val <= threshold_value ) ) {
                        if ( fraction < 1 ) {
                            if ( uniform(rng) > fraction ) continue;
                        }
                        Bounds xyz = indices_to_center( idx );
                        f << boost::str(boost::format("HETATM%5i VOXL VOX A%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s")
                            %iatom%ires%xyz[0]%xyz[1]%xyz[2]%1.0%1.0%"HB") << std::endl;
                    }
                }
            }
        }
        f.close();
    }

};

template< size_t D, class F, class V >
std::ostream & operator << ( std::ostream & out, VoxelArray<D,F,V> const & v ){
	out << "VoxelArray( lb: " << v.lb_ << " ub: " << v.ub_ << " cs: " << v.cs_ << " nelem: " << v.num_elements() << " sizeof_val: " << sizeof(V) << " )";
	return out;
}

}}}

#endif
