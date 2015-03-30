#ifndef INCLUDED_chemical_RotamerIndex_HH
#define INCLUDED_chemical_RotamerIndex_HH

#include "scheme/util/exit.hh"
#include "scheme/util/str.hh"

#include <boost/foreach.hpp>

#include <string>
#include <iostream>
#include <map>

namespace scheme { namespace chemical {

namespace impl {

template< class _Atom >
struct Rotamer {
	typedef _Atom Atom;
	typedef struct { typename Atom::Position position; int32_t type; } AtomPT;
	std::string resname_;
	size_t n_proton_chi_;
	std::vector<float> chi_;
	std::vector<Atom> atoms_;
	std::vector<std::pair<Atom,Atom> > hbonders_;
};

}

template<class _Atom, class RotamerGenerator>
struct RotamerIndex {
	BOOST_STATIC_ASSERT( boost::is_same< _Atom, typename RotamerGenerator::Atom >::value );
	typedef _Atom Atom;
	typedef impl::Rotamer<Atom> Rotamer;
	typedef std::map<std::string,std::pair<size_t,size_t> > BoundsMap;

	std::vector< Rotamer > rotamers_;
	std::vector< size_t > primary_rotamers_;
	std::vector< size_t > parent_rotamer_;
	RotamerGenerator rotgen_;
	BoundsMap bounds_map_;

	RotamerIndex(){}

	// assumes primary parent is added firrt among its class
	size_t
	add_rotamer(
		std::string resname,
		std::vector<float> const & chi,
		size_t n_proton_chi,
		int parent_key = -1
	){
		Rotamer r;
		rotgen_.get_atoms( resname, chi, r.atoms_, r.hbonders_ );
		r.resname_ = resname;
		r.chi_ = chi;
		r.n_proton_chi_ = n_proton_chi;
		rotamers_.push_back(r);
		size_t this_key = rotamers_.size()-1;
		if( parent_key >= 0 ) parent_rotamer_.push_back( parent_key );
		else                  parent_rotamer_.push_back( this_key );
		return this_key;
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
		for(int i = 0; i < rotamers_.size(); ++i){
			if( parent_rotamer_[i]==i ) primary_rotamers_.push_back(i);
		}
		sanity_check();
	}

	bool sanity_check() const {
		std::cout << "SANITY_CHECK" << std::endl;
		typedef std::pair<std::string,std::pair<size_t,size_t> > TMP;
		BOOST_FOREACH( TMP const & tmp, bounds_map_ ){
			// std::cout << tmp.first << " " << tmp.second.first << " " << tmp.second.second << std::endl;
			ALWAYS_ASSERT_MSG( tmp.second.first <= tmp.second.second, "RotamerIndex::sanity_check FAIL" );
		}
		for(int i = 0; i < primary_rotamers_.size(); ++i){
			ALWAYS_ASSERT( parent_rotamer_[primary_rotamers_[i]]==primary_rotamers_[i] );
		}
		for( size_t irot = 0; irot < rotamers_.size(); ++irot ){
			size_t ipri = parent_rotamer_[irot];
			ALWAYS_ASSERT_MSG( rotamers_[irot].resname_ == rotamers_[ipri].resname_, "RotamerIndex::sanity_check FAIL" );
			ALWAYS_ASSERT_MSG( rotamers_[irot].chi_.size() == rotamers_[ipri].chi_.size(), "RotamerIndex::sanity_check FAIL" );
			ALWAYS_ASSERT( rotamers_[ipri].n_proton_chi_ <= 1 );
			size_t nchi = rotamers_[ipri].chi_.size() - rotamers_[ipri].n_proton_chi_;
			for(size_t ichi = 0; ichi < nchi; ++ichi){
				float child_chi = rotamers_[irot].chi_[ichi];
				float parent_chi = rotamers_[ipri].chi_[ichi];
				float chidiff = fabs(parent_chi-child_chi);
				chidiff = std::min( chidiff, 360.0f-chidiff );
				// std::cout << rotamers_[irot].resname_ << " " << ichi << " " << parent_chi << " " << child_chi << " " << chidiff << std::endl;
				if( ichi==nchi-1 ){
					ALWAYS_ASSERT_MSG( chidiff < 30.0, "parent chi more than 30° from child! "+rotamers_[irot].resname_+", "+str(parent_chi)+" vs "+str(child_chi) );
				} else {
					ALWAYS_ASSERT_MSG( chidiff < 20.0, "parent chi more than 20° from child! "+rotamers_[irot].resname_+", "+str(parent_chi)+" vs "+str(child_chi) );					
				}
			}
		}
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
		int rescount = 0;
		for(size_t irot = b.first; irot < b.second; ++irot){
			std::cout << rotamers_[irot].resname_;
			for(int i = 0; i < rotamers_[irot].chi_.size(); ++i) std::cout << " " << rotamers_[irot].chi_[i];
			std::cout << std::endl;
			out << "MODEL " << rotamers_[irot].resname_ << " " << ++rescount << " " << irot << std::endl;
			BOOST_FOREACH( Atom const & a, rotamers_[irot].atoms_ ){
				out << scheme::io::dump_pdb_atom(a) << std::endl;
			}
			typedef std::pair<Atom,Atom> TMP;
			BOOST_FOREACH( TMP const & h, rotamers_[irot].hbonders_ ){
				out << scheme::io::dump_pdb_atom(h.first) << std::endl;
				out << scheme::io::dump_pdb_atom(h.second) << std::endl;				
			}
			out << "ENDMDL " << irot << std::endl;			
		}
	}

};


}}

#endif
