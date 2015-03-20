#ifndef INCLUDED_chemical_RotamerIndex_HH
#define INCLUDED_chemical_RotamerIndex_HH

#include <string>

#include <iostream>

namespace scheme { namespace chemical {

template< class Atom >
struct Rotamer {
	std::string resname_;
	std::vector<float> chi_;
	std::vector<Atom> atoms_;
};

template<class Atom, class RotamerGenerator>
struct RotamerIndex {
	std::vector<Rotamer<Atom> > rotamers_;
	RotamerGenerator rotgen_;

	RotamerIndex(){}

	void add_rotamer( std::string resname, std::vector<float> const & chi ){
		std::vector<Atom> atoms;
		rotgen_.get_atoms( resname, chi, atoms );
		Rotamer<Atom> r;
		r.resname_ = resname;
		r.chi_ = chi;
		r.atoms_ = atoms;
		rotamers_.push_back(r);
	}
};


}}

#endif
