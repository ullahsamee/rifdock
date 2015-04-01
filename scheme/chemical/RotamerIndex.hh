#ifndef INCLUDED_chemical_RotamerIndex_HH
#define INCLUDED_chemical_RotamerIndex_HH

#include "scheme/util/assert.hh"
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
	typedef std::map<std::string,std::pair<int,int> > BoundsMap;
	typedef std::vector< std::pair<int,int> > ChildMap; // index is "primary" index, not rotamer number

	std::vector< Rotamer > rotamers_;
	std::vector< int > primary_rotamers_;
	std::vector< int > parent_rotamer_;
	std::vector< int > parent_primary_;
	RotamerGenerator rotgen_;
	BoundsMap bounds_map_;
	ChildMap child_map_;
	std::map<int,int> parent_to_primary_;

	RotamerIndex(){}

	// assumes primary parent is added firrt among its class
	int
	add_rotamer(
		std::string resname,
		std::vector<float> const & chi,
		int n_proton_chi,
		int parent_key = -1
	){
		Rotamer r;
		rotgen_.get_atoms( resname, chi, r.atoms_, r.hbonders_ );
		r.resname_ = resname;
		r.chi_ = chi;
		r.n_proton_chi_ = n_proton_chi;
		rotamers_.push_back(r);
		int this_key = rotamers_.size()-1;
		if( parent_key >= 0 ) parent_rotamer_.push_back( parent_key );
		else                  parent_rotamer_.push_back( this_key );
		return this_key;
	}

	void 
	build_index()
	{
		// bounds
		bounds_map_[ rotamers_.front().resname_ ].first = 0;
		for(int i = 1; i < rotamers_.size(); ++i){
			if( rotamers_[i].resname_ != rotamers_[i-1].resname_ ){
				bounds_map_[ rotamers_[i-1].resname_ ].second = i;
				bounds_map_[ rotamers_[i  ].resname_ ].first  = i;
			}
		}
		bounds_map_[ rotamers_.back().resname_ ].second = rotamers_.size();
		
		// primary rotamers
		for(int i = 0; i < rotamers_.size(); ++i){
			if( parent_rotamer_[i]==i ) primary_rotamers_.push_back(i);
			parent_to_primary_[i] = primary_rotamers_.size()-1;
		}

		child_map_.resize(primary_rotamers_.size(),std::make_pair( std::numeric_limits<int>::max(),std::numeric_limits<int>::min() ));
		for(int i = 0; i < rotamers_.size(); ++i){
			int ipri = parent_to_primary_[ parent_rotamer_[i] ];
			child_map_[ipri].first  = std::min( child_map_[ipri].first , i );
			child_map_[ipri].second = std::max( child_map_[ipri].second, i );
		}
		for(int ipri = 0; ipri < child_map_.size(); ++ipri) ++child_map_[ipri].second;

		sanity_check();
	}

	bool
	sanity_check() const 
	{
		std::cout << "SANITY_CHECK" << std::endl;
		typedef std::pair<std::string,std::pair<int,int> > TMP;
		BOOST_FOREACH( TMP const & tmp, bounds_map_ ){
			// std::cout << tmp.first << " " << tmp.second.first << " " << tmp.second.second << std::endl;
			ALWAYS_ASSERT_MSG( tmp.second.first <= tmp.second.second, "RotamerIndex::sanity_check FAIL" );
		}
		for(int i = 0; i < primary_rotamers_.size(); ++i){
			ALWAYS_ASSERT( parent_rotamer_.at(primary_rotamers_.at(i))==primary_rotamers_.at(i) );
		}
		for( int irot = 0; irot < rotamers_.size(); ++irot ){
			int ipri = parent_rotamer_[irot];
			ALWAYS_ASSERT_MSG( rotamers_.at(irot).resname_ == rotamers_.at(ipri).resname_, "RotamerIndex::sanity_check FAIL" );
			ALWAYS_ASSERT_MSG( rotamers_.at(irot).chi_.size() == rotamers_.at(ipri).chi_.size(), "RotamerIndex::sanity_check FAIL" );
			ALWAYS_ASSERT( rotamers_.at(ipri).n_proton_chi_ <= 1 );
			int nchi = rotamers_.at(ipri).chi_.size() - rotamers_.at(ipri).n_proton_chi_;
			for(int ichi = 0; ichi < nchi; ++ichi){
				float child_chi = rotamers_.at(irot).chi_.at(ichi);
				float parent_chi = rotamers_.at(ipri).chi_.at(ichi);
				float chidiff = fabs(parent_chi-child_chi);
				chidiff = std::min( chidiff, 360.0f-chidiff );
				// std::cout << rotamers_.at(irot).resname_ << " " << ichi << " " << parent_chi << " " << child_chi << " " << chidiff << std::endl;
				if( ichi==nchi-1 ){
					ALWAYS_ASSERT_MSG( chidiff < 30.0, "parent chi more than 30° from child! "
						+rotamers_.at(irot).resname_+", "+str(parent_chi)+" vs "+str(child_chi) );
				} else {
					ALWAYS_ASSERT_MSG( chidiff < 20.0, "parent chi more than 20° from child! "
						+rotamers_.at(irot).resname_+", "+str(parent_chi)+" vs "+str(child_chi) );					
				}
			}
		}

		for( int irot = 0; irot < rotamers_.size(); ++irot ){
			int ipri = parent_to_primary_.find( parent_rotamer_.at(irot) )->second;
			ALWAYS_ASSERT( child_map_.at(ipri).first <= irot );
			ALWAYS_ASSERT( child_map_.at(ipri).second > irot )
		}

		for( int ipri = 0; ipri < primary_rotamers_.size(); ++ipri ){
			int iparent = primary_rotamers_.at(ipri);
			for( int irot = child_map_[ipri].first; irot < child_map_[ipri].second; ++irot ){
				ALWAYS_ASSERT( parent_rotamer_.at(irot) == iparent );
			}
		}

	}

	std::vector<float> const & chis(int rotnum) const { return rotamers_.at(rotnum).chi_; }
	std::vector<Atom> const & atoms(int rotnum) const { return rotamers_.at(rotnum).atoms_; }

	std::pair<int,int>
	index_bounds( std::string const & resname ) const { 
		BoundsMap::const_iterator i = bounds_map_.find(resname);
		if( i == bounds_map_.end() ) return std::make_pair(0,0);
		return i->second;
	}

	std::pair<int,int>
	child_bounds_of_primary( int ipri ) const { 
		return child_map_.at(ipri);
	}


	size_t size() const { return rotamers_.size(); }

	void dump_pdb( std::ostream & out, std::string resname="") const {
		std::pair<int,int> b(0,rotamers_.size());
		if( resname.size() ) b = index_bounds(resname);
		int rescount = 0;
		for(int irot = b.first; irot < b.second; ++irot){
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
