#ifndef INCLUDED_objective_hash_XformMap_HH
#define INCLUDED_objective_hash_XformMap_HH

#include "scheme/util/SimpleArray.hh"
#include "scheme/util/dilated_int.hh"
#include "scheme/numeric/FixedPoint.hh"
#include "scheme/nest/maps/TetracontoctachoronMap.hh"
#include "scheme/numeric/bcc_lattice.hh"
#include "scheme/objective/hash/XformHash.hh"
#include "scheme/objective/hash/XformHashNeighbors.hh"

#include <sparsehash/dense_hash_map>

namespace scheme { namespace objective { namespace hash {



template< int ArrayBits, class Key, class Val >
struct XfromMapSerializer {
    typedef util::SimpleArray< (1<<ArrayBits), Val >  ValArray;
	bool operator()( std::istream * in, std::pair<Key const,ValArray> * val ) const {
		Key & k = const_cast<Key&>( val->first );
		in->read( (char*)&k, sizeof(Key) );
		in->read( (char*)&val->second, sizeof(ValArray) );
		// std::cout << "READ " << k << " " << val->second << std::endl;
		return true;
	}
	bool operator()( std::ostream * out, std::pair<Key const,ValArray> const & val ) const {
		out->write( (char*)&val.first, sizeof(Key) );
		out->write( (char*)&val.second, sizeof(ValArray) );		
		return true;
	}
};


template<
	class Xform,
	class Val=numeric::FixedPoint<-17>,
	int ArrayBits=7,
	// template<class X> class _Hasher = XformHash_bt24_BCC6_Zorder >
	template<class X> class _Hasher = XformHash_Quat_BCC7_Zorder ,
	class ElementSerializer = XfromMapSerializer< ArrayBits, uint64_t, Val >
>
struct XformMap {
	BOOST_STATIC_ASSERT( ArrayBits >= 0 );
	typedef _Hasher<Xform> Hasher;
	typedef uint64_t Key;
	typedef typename Xform::Scalar Float;
    typedef util::SimpleArray< (1<<ArrayBits), Val >  ValArray;
    typedef google::dense_hash_map<Key,ValArray> Map;
    Hasher hasher_;
    Map map_;
	ElementSerializer element_serializer_;
    Float cart_resl_, ang_resl_, cart_bound_;

	XformMap( Float cart_resl=-1, Float ang_resl=-1, Float cart_bound=512.0 )
	{
		map_.set_empty_key( std::numeric_limits<Key>::max() );
		cart_resl_ = cart_resl;
		ang_resl_ = ang_resl;
		cart_bound_ = cart_bound;
		if( cart_resl_ != -1 && ang_resl != -1 ){
			hasher_.init( cart_resl, ang_resl, cart_bound );
		}
	}

	Val & operator[]( Xform const & x ){
		Key k = hasher_.get_key( x );
		Key k0 = k >> ArrayBits;
		Key k1 = k & (((Key)1<<ArrayBits)-1);
		typename Map::iterator iter = map_.find(k0);
		if( iter == map_.end() ){
			ValArray aval(0);
			iter = map_.insert( std::make_pair(k0,aval) ).first; // TODO: should check for failuer here
		}
		return iter->second[k1];
	}
	Val operator[]( Xform const & x ) const {
		Key k = hasher_.get_key( x );
		Key k0 = k >> ArrayBits;
		Key k1 = k & (((Key)1<<ArrayBits)-1);
		typename Map::const_iterator iter = map_.find(k0);
		if( iter == map_.end() ){ return 0.0; }
		return iter->second[k1];
	}

	void insert_sphere(
		Xform const & x,
		Val value,
		XformHashNeighbors<Hasher> & cache
	){
		// TODO
	}

	bool save( std::ostream & out ) {
		// no way to check if the stream was opened binary!
		// if( ! (out.flags() & std::ios::binary) ){
		// 	std::cerr << "XformMap::save must be binary ostream" << std::endl;
		// 	return false;
		// }
		std::ostringstream oss;
		oss << std::endl;
		oss << "=========== description ===========" << std::endl;
		oss << "Scheme Xform Map" << std::endl;
		oss << "Hasher: " << hasher_.name() << std::endl;
		oss << "Cart Resolution: " << cart_resl_ << std::endl;
		oss << "Angular Resolution: " << ang_resl_ << std::endl;
		oss << "=========== begin binary data ===========" << std::endl;
		size_t s = oss.str().size();
		out.write((char*)&s,sizeof(size_t));
		out.write(oss.str().c_str(),s);
		// begin binary data
		s = hasher_.name().size();
		// std::cout << "SIZE " << s << std::endl;
		out.write( (char*)&s, sizeof(size_t) );
		out.write( hasher_.name().c_str(), hasher_.name().size()*sizeof(char) );
		out.write( (char*)&cart_resl_, sizeof(Float) );
		out.write( (char*)&ang_resl_, sizeof(Float) );		
		out.write( (char*)&cart_bound_, sizeof(Float) );		
		// std::cout << "SIZE OUT " << map_.size() << std::endl;
		if( ! map_.serialize( element_serializer_, &out ) ){
			std::cerr << "XfromMap::load failed to unserialize sparsehash" << std::endl;
			return false;
		}
		return true;
	}
	bool load( std::istream & in ) {
		// no way to check if the stream was opened binary!
		// if( ! (in.flags() & std::ios::binary) ){
		// 	std::cerr << "XformMap::save must be binary ostream" << std::endl;
		// 	return false;
		// }
		size_t s;
		in.read((char*)&s,sizeof(size_t));
		char buf[9999];
		for(int i = 0; i < 9999; ++i) buf[i] = 0;
		in.read(buf,s);
		std::cout << "XformMap load, description: " << std::endl;
		std::cout << std::string(buf).substr(37,s-80) << std::endl;
		in.read( (char*)&s, sizeof(size_t) );
		// std::cout << "SIZE IN " << s << std::endl;
		for(int i = 0; i < 9999; ++i) buf[i] = 0;
		in.read( buf, s );
		// std::cout << "BUF: " <<  buf << " " << std::string(buf) <<  std::endl;
		if( hasher_.name() != std::string(buf) ){
			std::cerr << "XformMap::load, hasher type mismatch, expected " << hasher_.name() << " got "  << buf << std::endl;
			return false;
		}
		Float cart_resl, ang_resl, cart_bound;
		in.read((char*)&cart_resl,sizeof(Float));
		in.read((char*)&ang_resl,sizeof(Float));
		in.read((char*)&cart_bound,sizeof(Float));
		if( cart_resl_ != -1 && cart_resl_ != cart_resl ){
			std::cerr << "XformMap::load, hasher cart_resl mismatch, expected " << cart_resl_ << " got "  << cart_resl << std::endl;
			return false;			
		}
		if( ang_resl_ != -1 && ang_resl_ != ang_resl ){
			std::cerr << "XformMap::load, hasher ang_resl mismatch, expected " << ang_resl_ << " got "  << ang_resl << std::endl;
			return false;			
		}
		if( cart_resl_ == -1 ) cart_resl_ = cart_resl;
		if(  ang_resl_ == -1 )  ang_resl_ =  ang_resl;
		cart_bound_ = cart_bound;
		hasher_.init( cart_resl_, ang_resl_, cart_bound_ );

		if( ! map_.unserialize( element_serializer_, &in ) ){
			std::cerr << "XfromMap::load failed to unserialize sparsehash" << std::endl;
			return false;
		}

		// std::cout << "SIZE IN " << map_.size() << std::endl;
		return true;
	}


};


}}}

#endif
