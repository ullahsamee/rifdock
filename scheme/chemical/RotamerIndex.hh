#ifndef INCLUDED_chemical_RotamerIndex_HH
#define INCLUDED_chemical_RotamerIndex_HH

#include "scheme/util/assert.hh"
#include "scheme/util/str.hh"
#include "scheme/chemical/AtomData.hh"
#include "scheme/io/dump_pdb_atom.hh"

#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>

#include <cstdlib>
#include <string>
#include <iostream>
#include <map>
#include <cmath>

#include <Eigen/Dense>

namespace scheme { namespace chemical {

struct HBondRay {
	::Eigen::Vector3f horb_cen, direction;
};

template< class _AtomData >
struct ChemicalIndex {
	typedef _AtomData AtomData;
	std::vector<std::string> resnames_;
	std::map<std::string,int> resname2num_;
	std::vector< std::vector< AtomData > > atomdata_;
	AtomData null_atomdata_;

	ChemicalIndex() : null_atomdata_() {
	}

	bool have_res( std::string resn ) const {
		return resname2num_.find(resn) != resname2num_.end();
	}

	void add_res( std::string resn ){
		if( have_res(resn) ) return;
		resname2num_[resn] = resnames_.size();
		resnames_.push_back(resn);
	}

	void add_atomdata( int restype, int atomnum, AtomData const & data ){
		if( atomdata_.size() <= restype ) atomdata_.resize(restype+1);
		if( atomdata_[restype].size() <= atomnum ) atomdata_[restype].resize(atomnum+1);
		atomdata_[restype][atomnum] = data;
	}
	void add_atomdata( std::string resname, int atomnum, AtomData const & data ){
		add_res( resname );
		add_atomdata( resname2num_[resname], atomnum, data );
	}

	AtomData       & atom_data( int restype, int atomnum )       { if(restype<0) return null_atomdata_; return atomdata_.at(restype).at(atomnum); }
	AtomData const & atom_data( int restype, int atomnum ) const { if(restype<0) return null_atomdata_; return atomdata_.at(restype).at(atomnum); }
	AtomData       & atom_data( std::string const & resname, int atomnum )       { return atomdata_.at( resname2num_.find(resname)->second ).at(atomnum); }
	AtomData const & atom_data( std::string const & resname, int atomnum ) const { return atomdata_.at( resname2num_.find(resname)->second ).at(atomnum); }

	int resname2num( std::string const & resn ) const {
		std::map<std::string,int>::const_iterator i = resname2num_.find(resn);
		if( i == resname2num_.end() ) return -1;
		return i->second;
	}

};



namespace impl {

template< class _Atom >
struct Rotamer {
	typedef _Atom Atom;
	// typedef struct { typename Atom::Position position; int32_t type; } AtomPT;
	std::string resname_;
	size_t n_proton_chi_;
	std::vector<float> chi_;
	std::vector<Atom> atoms_;
	std::vector<std::pair<Atom,Atom> > hbonders_;
	std::vector<HBondRay> donors_;
	std::vector<HBondRay> acceptors_;
	int nheavyatoms;
	uint64_t validation_hash() const {
		uint64_t h = (uint64_t)(boost::hash<std::string>()(resname_));
		h ^= (uint64_t)boost::hash<size_t>()(n_proton_chi_);
		for(int i = 0; i < chi_.size(); ++i){
			h ^= (uint64_t)boost::hash<float>()(chi_[i]);
		}
		return h;
	}
};

}


template<class _Atom, class RotamerGenerator>
struct RotamerIndex {
	BOOST_STATIC_ASSERT(( ( boost::is_same< _Atom, typename RotamerGenerator::Atom >::value ) ));
	typedef _Atom Atom;
	typedef impl::Rotamer<Atom> Rotamer;
	typedef std::map<std::string,std::pair<int,int> > BoundsMap;
	typedef std::vector< std::pair<int,int> > ChildMap; // index is "primary" index, not rotamer number

	size_t size() const { return rotamers_.size(); }
	std::string resname(size_t i) const { return rotamers_.at(i).resname_; }
	char oneletter(size_t i) const { return oneletter_map_.at(rotamers_.at(i).resname_); }
	size_t natoms( size_t i ) const { return rotamers_.at(i).atoms_.size(); }
	size_t nheavyatoms( size_t i ) const { return rotamers_.at(i).nheavyatoms; }
	size_t nchi( size_t i ) const { return rotamers_.at(i).chi_.size(); }
	size_t nchi_noproton( size_t i ) const { return nchi(i)-nprotonchi(i); }
	size_t nprotonchi( size_t i ) const { return rotamers_.at(i).n_proton_chi_; }
	float   chi( size_t i, size_t j ) const { return rotamers_.at(i).chi_.at(j); }
	size_t nhbonds( size_t i ) const { return rotamers_.at(i).hbonders_.size(); }
	Atom const & hbond_atom1( size_t i, size_t ihb ) const { return rotamers_.at(i).hbonders_.at(ihb).first; }
	Atom const & hbond_atom2( size_t i, size_t ihb ) const { return rotamers_.at(i).hbonders_.at(ihb).second; }

	std::map<std::string,char> oneletter_map_;

	ChemicalIndex<AtomData> chem_index_;

	std::vector< Rotamer > rotamers_;
	Rotamer const & rotamer( int irot ) const { return rotamers_[irot]; }

	std::vector< int > primary_rotamers_; // rot # of non-ex rotamers
	std::vector< int > parent_rotamer_; //
	std::vector< int > parent_primary_; // pri-rot # of primary parent
	RotamerGenerator rotgen_;
	BoundsMap bounds_map_;
	ChildMap child_map_;
	std::map<int,int> parent_to_primary_;
	std::vector<int> protonchi_parent_;
	int ala_rot_;

	RotamerIndex(){
		std::cout << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << " " << "" << std::endl;
		this->fill_oneletter_map( oneletter_map_ );
		std::cout << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << " " << "" << std::endl;
	}

	// assumes primary parent is added firrt among its class
	int
	add_rotamer(
		std::string resname,
		std::vector<float> const & chi,
		int n_proton_chi,
		int parent_key = -1
	){
		Rotamer r;
		rotgen_.get_atoms( resname, chi, r.atoms_, r.hbonders_, r.nheavyatoms, r.donors_, r.acceptors_ );
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

		for(int i = 0; i < rotamers_.size(); ++i){
			if( chem_index_.have_res( rotamers_[i].resname_ ) ) continue;
			for( int ia = 0; ia < rotamers_[i].atoms_.size(); ++ia ){
				chem_index_.add_atomdata( rotamers_[i].resname_, ia, rotamers_[i].atoms_[ia].data() );
			}
		}


		protonchi_parent_.resize(this->size(),-2);
		protonchi_parent_[0] = 0;
		int last = 0;
		for( int i = 1; i < this->size(); ++i ){
			// cout << this->resname(i) << " " << this->nchi(i) << " " << this->nchi_noproton(i) << endl;
			if( this->resname(i) != this->resname(i-1) ){
				protonchi_parent_[i] = i;
				last = i;
			} else {
				for( int j = 0; j < this->nchi_noproton(i); ++j ){
					if( this->chi(i,j) != this->chi(i-1,j) ){
						protonchi_parent_[i] = i;
						last = i;
					}
				}
			}
			if( protonchi_parent_[i] == -2 ) protonchi_parent_[i] = last;
		}
		int uniq_sum = 0;
		for( int i = 0; i < this->size(); ++i ){
			// cout << I(3,i) << " " << this->resname(i);
			for( int j = 0; j < this->nchi(i); ++j ){
				// cout << " " << F(7,2,this->chi(i,j));
			}
			// cout << " " << I(3,protonchi_parent_[i]) << endl;
			uniq_sum += ( i == protonchi_parent_[i] );
		}
		std::cout << "n protonchi_parent_ " << uniq_sum << std::endl;

		for( int i = 0; i < this->size(); ++i ){
			if( resname(i) == "ALA" ){
				ala_rot_ = i;
				break;
			}
		}

		sanity_check();
	}

	bool
	sanity_check() const
	{
		// std::cout << "SANITY_CHECK" << std::endl;
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
				float chidiff = std::fabs(parent_chi-child_chi);
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

		for( int irot = 0; irot < rotamers_.size(); ++irot ){

			for( int ia = 0; ia < rotamers_[irot].atoms_.size(); ++ia ){
				AtomData const & ad1( rotamers_[irot].atoms_[ia].data() );
				// std::cout <<  rotamers_[irot].resname_ << std::endl;
				AtomData const & ad2( chem_index_.atom_data( rotamers_[irot].resname_, ia ) );
				ALWAYS_ASSERT( ad1 == ad2 );
			}
		}

		for( int i = 0; i < this->size(); ++i ){
			ALWAYS_ASSERT( 0 <= protonchi_parent_[i] && protonchi_parent_[i] < this->size() );
			ALWAYS_ASSERT( resname( protonchi_parent_[i] ) == resname(i) )
			// this is weak....
		}

		ALWAYS_ASSERT( resname(ala_rot_) == "ALA" );
	}

	std::vector<float> const & chis(int rotnum) const { return rotamers_.at(rotnum).chi_; }
	std::vector<Atom> const & atoms(int rotnum) const { return rotamers_.at(rotnum).atoms_; }

	int ala_rot() const { return ala_rot_; }

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

	bool is_primary( int irot ) const { return parent_rotamer_.at(irot)==irot; }

	void dump_pdb( std::ostream & out, std::string resname="" ) const {
		std::pair<int,int> b(0,rotamers_.size());
		if( resname.size() ) b = index_bounds(resname);
		int rescount = 0;
		for(int irot = b.first; irot < b.second; ++irot){
			// std::cout << rotamers_[irot].resname_ << " " << irot;
			// for(int i = 0; i < rotamers_[irot].chi_.size(); ++i) std::cout << " " << rotamers_[irot].chi_[i];
			// std::cout << std::endl;
			out << "MODEL " << rotamers_[irot].resname_ << " " << ++rescount << " " << irot << std::endl;
			BOOST_FOREACH( Atom const & a, rotamers_[irot].atoms_ ){
				out << scheme::io::dump_pdb_atom(a) << std::endl;
			}
			// out << "ENDMDL" << irot << std::endl;
			// out << "MODEL " << rotamers_[irot].resname_ << " " << ++rescount << " " << irot << " HBONDERS" << std::endl;
			typedef std::pair<Atom,Atom> TMP;
			BOOST_FOREACH( TMP const & h, rotamers_[irot].hbonders_ ){
				out << scheme::io::dump_pdb_atom(h.first) << std::endl;
				out << scheme::io::dump_pdb_atom(h.second) << std::endl;
			}
			BOOST_FOREACH( HBondRay const & hr, rotamers_[irot].donors_ ){
				scheme::io::dump_pdb_atom_resname_atomname( out, "DON", "CDON", hr.horb_cen );
				scheme::io::dump_pdb_atom_resname_atomname( out, "DON", "DDON", hr.horb_cen+hr.direction*0.1 );
			}
			BOOST_FOREACH( HBondRay const & hr, rotamers_[irot].acceptors_ ){
				scheme::io::dump_pdb_atom_resname_atomname( out, "ACC", "CACC", hr.horb_cen );
				scheme::io::dump_pdb_atom_resname_atomname( out, "ACC", "DACC", hr.horb_cen+hr.direction*0.1 );
			}
			out << "ENDMDL" << std::endl;
		}
	}

	template< class Xform >
	void dump_pdb( std::ostream & out, int irot, Xform x, int ires ) const {
		for( int ia = 0; ia < rotamers_[irot].atoms_.size(); ++ia ){
			io::dump_pdb_atom( out, x*rotamers_[irot].atoms_[ia].position(), rotamers_[irot].atoms_[ia].data(), ires );
		}
		typedef std::pair<Atom,Atom> TMP;
		BOOST_FOREACH( TMP const & h, rotamers_[irot].hbonders_ ){
			// io::dump_pdb_atom( out, x*h.first .position(), h.first .data(), ires );
			if( h.second.type() < 0 ){ // is acceptor orbital
				io::dump_pdb_atom( out, x*h.second.position(), h.second.data(), ires );
			}
		}
	}



	void load( std::istream & is ){
		// read resn/chi angles and rebuild?
	}

	void save( std::ostream & os ){
		// save only resn and chi angles?
	}

	uint64_t validation_hash() const {
		uint64_t h = 0;
		BOOST_FOREACH( Rotamer const & rot, rotamers_  ){
			h ^= rot.validation_hash();
		}
		return h;
	}

	int protonchi_parent( int i ) const { return protonchi_parent_[i]; }


	void fill_oneletter_map( std::map<std::string,char> & oneletter_map ){
		oneletter_map["ALA"] = 'A';
		oneletter_map["CYS"] = 'C';
		oneletter_map["ASP"] = 'D';
		oneletter_map["GLU"] = 'E';
		oneletter_map["PHE"] = 'F';
		oneletter_map["GLY"] = 'G';
		oneletter_map["HIS"] = 'H';
		oneletter_map["ILE"] = 'I';
		oneletter_map["LYS"] = 'K';
		oneletter_map["LEU"] = 'L';
		oneletter_map["MET"] = 'M';
		oneletter_map["ASN"] = 'N';
		oneletter_map["PRO"] = 'P';
		oneletter_map["GLN"] = 'Q';
		oneletter_map["ARG"] = 'R';
		oneletter_map["SER"] = 'S';
		oneletter_map["THR"] = 'T';
		oneletter_map["VAL"] = 'V';
		oneletter_map["TRP"] = 'W';
		oneletter_map["TYR"] = 'Y';
		oneletter_map["ADE"] = 'a';
		oneletter_map["CYT"] = 'c';
		oneletter_map["GUA"] = 'g';
		oneletter_map["THY"] = 't';
		oneletter_map["RAD"] = 'a';
		oneletter_map["RCY"] = 'c';
		oneletter_map["RGU"] = 'g';
		oneletter_map["URA"] = 'u';
		oneletter_map["H2O"] = 'w';
		oneletter_map["UNP"] = 'z';
		oneletter_map["UNK"] = 'Z';
		oneletter_map["VRT"] = 'X';
		oneletter_map["ala"] = 'A';
		oneletter_map["cys"] = 'C';
		oneletter_map["asp"] = 'D';
		oneletter_map["glu"] = 'E';
		oneletter_map["phe"] = 'F';
		oneletter_map["gly"] = 'G';
		oneletter_map["his"] = 'H';
		oneletter_map["ile"] = 'I';
		oneletter_map["lys"] = 'K';
		oneletter_map["leu"] = 'L';
		oneletter_map["met"] = 'M';
		oneletter_map["asn"] = 'N';
		oneletter_map["pro"] = 'P';
		oneletter_map["gln"] = 'Q';
		oneletter_map["arg"] = 'R';
		oneletter_map["ser"] = 'S';
		oneletter_map["thr"] = 'T';
		oneletter_map["val"] = 'V';
		oneletter_map["trp"] = 'W';
		oneletter_map["tyr"] = 'Y';
		oneletter_map["ade"] = 'a';
		oneletter_map["cyt"] = 'c';
		oneletter_map["gua"] = 'g';
		oneletter_map["thy"] = 't';
		oneletter_map["rad"] = 'a';
		oneletter_map["rcy"] = 'c';
		oneletter_map["rgu"] = 'g';
		oneletter_map["ura"] = 'u';
		oneletter_map["h2o"] = 'w';
		oneletter_map["unp"] = 'z';
		oneletter_map["unk"] = 'Z';
		oneletter_map["vrt"] = 'X';
	}

};


template<class A, class RG>
std::ostream & operator << ( std::ostream & out, RotamerIndex<A,RG> const & ridx ){
	out << "RotamerIndex:" << std::endl;
	std::pair<int,int> ib;
	ib=ridx.index_bounds("ALA"); out<<"    ALA "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("CYS"); out<<"    CYS "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("ASP"); out<<"    ASP "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("GLU"); out<<"    GLU "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("PHE"); out<<"    PHE "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("GLY"); out<<"    GLY "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("HIS"); out<<"    HIS "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("ILE"); out<<"    ILE "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("LYS"); out<<"    LYS "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("LEU"); out<<"    LEU "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("MET"); out<<"    MET "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("ASN"); out<<"    ASN "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("PRO"); out<<"    PRO "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("GLN"); out<<"    GLN "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("ARG"); out<<"    ARG "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("SER"); out<<"    SER "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("THR"); out<<"    THR "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("VAL"); out<<"    VAL "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("TRP"); out<<"    TRP "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("TYR"); out<<"    TYR "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
 	return out;
}



}}

#endif
