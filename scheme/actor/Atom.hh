#ifndef INCLUDED_actor_Atom_HH
#define INCLUDED_actor_Atom_HH

#include "scheme/numeric/util.hh"
#include "scheme/types.hh"
#include "scheme/chemical/AtomData.hh"
#include "scheme/io/dump_pdb_atom.hh"
#include <iostream>

namespace scheme {
namespace actor {


template<class _Position>
struct SimpleAtom {
	typedef _Position Position;

	SimpleAtom() : 
	    position_(0,0,0),
	    type_(0)
	{}

	template<class P>
	SimpleAtom(
		P const & p,
		int32_t type = 0
	) : position_(p[0],p[1],p[2]),
	    type_(type)
	{}

	template<class Xform>
	SimpleAtom(SimpleAtom const & a,Xform const & moveby){
		position_ = moveby*a.position();
		type_ = a.type_;
	}

	int type() const { return type_; }
	void set_type(int i){ type_ = i; }

	Position const & position() const { return position_; }

	template<class P>
	void set_position( P const & pos ) { position_[0] = pos[0]; position_[1] = pos[1]; position_[2] = pos[2]; }
	
	bool operator==(SimpleAtom<Position> const & o) const {
		return numeric::approx_eq(o.position_,position_) && o.type_==type_; }

private:
	Position position_;
	int32_t type_;
};


template<class P>
std::ostream & operator<<(std::ostream & out,SimpleAtom<P> const & x){
	return out << "SimpleAtom( " << x.position() << ", " << x.type() << " )";
}

template<class P>
void write_pdb( std::ostream & out, SimpleAtom<P> const & a ){
	io::dump_pdb_atom( out, 0, a.position() );
}

using chemical::AtomData;



template<class _Position>
struct Atom {
	typedef _Position Position;

	Atom() : 
	    position_(0,0,0),
	    type_(0),
	    data_(scheme::make_shared<AtomData>()) {}
	    // data_(new AtomData) {}	    

	template<class P>
	Atom(
		P const &    p,
		int                 type     = 0,
		std::string const & atomname = AtomData::default_atomname(),        
		std::string const & resname  = AtomData::default_resname(),       
		char                chain    = AtomData::default_chain(),     
		int                 resnum   = AtomData::default_resnum(),      
		int                 atomnum  = AtomData::default_atomnum(),       
		std::string const & elem     = AtomData::default_elem(),    
		bool                ishet    = AtomData::default_ishet(),     
		float               occ      = AtomData::default_occ(),   
		float               bfac     = AtomData::default_bfac()
	) : position_(p[0],p[1],p[2]),
	    type_(type),
	    data_(scheme::make_shared<AtomData>( //TODO: why do I need scheme:: here?
	    // data_(new AtomData(
	    	atomname,resname,chain,resnum,atomnum,elem,ishet,occ,bfac)
	    )
	{}

	// // TODO: compare raw ptr for speed
	// ~Atom(){ delete data_; }

	template<class Xform>
	Atom( Atom const & a, Xform const & moveby ){
		position_ = moveby*a.position();
		type_ = a.type_;
		data_ = a.data_;
	}

	int type() const { return type_; }
	void set_type(int i){ type_ = i; }
	AtomData const & data() const { return *data_; }

	Position const & position() const { return position_; }

	template<class P>
	void set_position( P const & pos ) { position_[0] = pos[0]; position_[1] = pos[1]; position_[2] = pos[2]; }
	
	bool operator==(Atom<Position> const & o) const {
		return numeric::approx_eq(o.position_,position_) && o.type_==type_ && o.data_==data_; }

private:
	int type_;
	shared_ptr<AtomData> data_;
	// AtomData *data_;
	Position position_;
};

template<class P>
std::ostream & operator<<(std::ostream & out,Atom<P> const & x){
	return out << "Atom( " << x.position() << ", " << x.type() << ", " << x.data() << " )";
}

template<class P>
void write_pdb(std::ostream & out,Atom<P> const & a){
	io::dump_pdb_atom(out, a.position(), a.data() );
}




}
}

#endif
