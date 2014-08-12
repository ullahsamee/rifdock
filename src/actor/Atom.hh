#ifndef INCLUDED_actor_Atom_HH
#define INCLUDED_actor_Atom_HH

#include <iostream>
#include "io/dump_pdb_atom.hh"

namespace scheme {
namespace actor {

template<class _Position>
struct Atom {
	typedef _Position Position;
	Position position_;
	int type_;

	Atom() : position_(0,0,0),type_(0) {}

	Atom(Position const & p, int type) : position_(p),type_(type) {}

	template<class Xform>
	Atom(Atom const & a,Position const & moveby){
		position_ = moveby*a.position();
		type_ = a.type_;
	}

	void set_position( Position const & pos	){ position_ = pos; }

	Position const & position() const { return position_; }
	
	bool operator==(Atom<Position> const & o) const { return o.position_.isApprox(position_) && o.type_ = type_; }
};

template<class P>
std::ostream & operator<<(std::ostream & out,Atom<P> const & x){
	return out << "Atom( " << x.position() << ", " << x.type_ << " )";
}

template<class P>
void dump_pdb(std::ostream & out,Atom<P> const & a){
	//TODO: this is bad, to assume rosetta atom types
	static char const elem[26] = {' ','C','C','C','C','C','C','N','N','N','N','N','N','O','O','O','O','S','N','C','C','O','P','P','H','H'};
	io::dump_pdb_atom(out, elem[a.type_], a.position() );
}

}
}

#endif
