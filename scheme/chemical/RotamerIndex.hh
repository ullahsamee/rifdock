#ifndef INCLUDED_chemical_RotamerIndex_HH
#define INCLUDED_chemical_RotamerIndex_HH

#include <string>

#include <iostream>

namespace scheme { namespace chemical {

namespace impl {

template< class _Atom >
struct Rotamer {
	typedef _Atom Atom;
	typedef struct { typename Atom::Position position; int32_t type; } AtomPT;
	std::string resname_;
	std::vector<float> chi_;
	std::vector<Atom> atoms_;
};

}

template<class _Atom, class RotamerGenerator>
struct RotamerIndex {
	typedef _Atom Atom;
	typedef impl::Rotamer<Atom> Rotamer;
	std::vector<Rotamer > rotamers_;
	RotamerGenerator rotgen_;
	typedef std::map<std::string,std::pair<size_t,size_t> > BoundsMap;
	BoundsMap bounds_map_;

	RotamerIndex(){}

	void add_rotamer( std::string resname, std::vector<float> const & chi ){
		std::vector<Atom> atoms;
		rotgen_.get_atoms( resname, chi, atoms );
		Rotamer r;
		r.resname_ = resname;
		r.chi_ = chi;
		r.atoms_ = atoms;
		rotamers_.push_back(r);
	}

	void build_index(){
		bounds_map_[ rotamers_.front().resname_ ].first = 0;
		for(int i = 1; i < rotamers_.size(); ++i){
			if( rotamers_[i].resname_ != rotamers_[i-1].resname_ ){
				bounds_map_[ rotamers_[i-1].resname_ ].second = i;
				bounds_map_[ rotamers_[i  ].resname_ ].first  = i;
			}
		}
		bounds_map_[ rotamers_.back().resname_ ].second = rotamers_.size();
	}

	std::vector<float> const & chis(int rotnum) const { return rotamers_.at(rotnum).chi_; }
	std::vector<Atom> const & atoms(int rotnum) const { return rotamers_.at(rotnum).atoms_; }

	std::pair<size_t,size_t>
	index_bounds( std::string const & resname ) const { 
		BoundsMap::const_iterator i = bounds_map_.find(resname);
		if( i == bounds_map_.end() ) return std::make_pair(0ul,0ul);
		return i->second;
	}

	size_t size() const { return rotamers_.size(); }

	void dump_pdb( std::ostream & out, std::string resname="") const {
		std::pair<size_t,size_t> b(0,rotamers_.size());
		if( resname.size() ) b = index_bounds(resname);
		for(size_t irot = b.first; irot < b.second; ++irot){
			out << "MODEL " << irot << std::endl;
			BOOST_FOREACH( Atom const & a, rotamers_[irot].atoms_ ){
				out << scheme::io::dump_pdb_atom(a) << std::endl;
			}
			out << "ENDMDL " << irot << std::endl;			
		}
	}

};


}}

#endif
