#ifndef INCLUDED_actor_BackboneHBondActor_HH
#define INCLUDED_actor_BackboneHBondActor_HH

#include <Eigen/Dense>
#include <scheme/io/dump_pdb_atom.hh>
#include <scheme/chemical/HBondRay.hh>

#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>

namespace scheme {
namespace actor {


// The hbond_rays_ holds the position information


	struct BackboneHBondActor {


		typedef ::Eigen::Transform<float,3,Eigen::AffineCompact> EigenXform;

		///@brief Position type, leave out to make actor "Fixed"
		typedef EigenXform Position;
		typedef typename Position::Scalar Float;
		typedef BackboneHBondActor THIS;
		typedef Eigen::Matrix<Float,3,1> V3;

		std::vector<::scheme::chemical::HBondRay> hbond_rays_;
		bool is_donor_;
		int index_;

		BackboneHBondActor() : is_donor_(false), index_(0)  {}

		BackboneHBondActor(std::vector<::scheme::chemical::HBondRay> const & hbond_rays, bool is_donor, int i=0) : 
			hbond_rays_(hbond_rays),
			is_donor_(is_donor),
			index_(i) 
		{}


		BackboneHBondActor(
			BackboneHBondActor const & actor0,
			Position const & to_moveby
		){
			hbond_rays_ = actor0.hbond_rays_;
			is_donor_ = actor0.is_donor_;
			index_ = actor0.index_;	

			moveby(to_moveby);
		}

		std::vector<::scheme::chemical::HBondRay> const &
		hbond_rays() const {
			return hbond_rays_;
		}

		bool
		is_donor() const {
			return is_donor_;
		}

		int
		index() const {
			return index_;
		}

		// void 
		// set_position(
		// 	Position const & pos
		// ){  }

		void 
		moveby(
			Position const & pos
		){ 
			for ( ::scheme::chemical::HBondRay & hbray : hbond_rays_ ) {
				hbray.apply_xform( pos );
			}
		}

		// Position const &
		// position() const { return Position::Identity(); }

		bool operator==(THIS const & o) const {
			return o.hbond_rays_==hbond_rays_ && 
			       o.is_donor_  ==is_donor_   &&
			       o.index_     ==index_;
		}

		///@brief necessary for testing only
		// bool operator<(THIS const & o) const { return std::make_pair(position_,aa_) < std::make_pair(o.position_,o.aa_); }

		///@brief necessary for testing only
	  	// template<class Archive> void serialize(Archive & ar, const unsigned int ){
	  	// 	ar & position_;
	  	// 	ar & aa_;
	  	// 	ar & ss_;	  		
	  	// 	ar & index_;
	  	// }

	};


template< class MetaData >
void write_pdb( std::ostream & out, BackboneHBondActor const & a, MetaData const & meta ){

	// AtomData(
	// 	std::string const & _atomname = AtomData::default_atomname(),        
	// 	std::string const & _resname  = AtomData::default_resname(),       
	// 	char                _chain    = AtomData::default_chain(),     
	// 	int                 _resnum   = AtomData::default_resnum(),      
	// 	int                 _atomnum  = AtomData::default_atomnum(),       
	// 	std::string const & _elem     = AtomData::default_elem(),    
	// 	bool                _ishet    = AtomData::default_ishet(),     
	// 	float               _occ      = AtomData::default_occ(),   
	// 	float               _bfac     = AtomData::default_bfac()
	// );

}


inline
std::ostream & operator<<(std::ostream & out,BackboneHBondActor const& a){
	return out << "BackboneHBondActor " << (a.is_donor_ ? "donor" : "acceptor" ) << " " << a.index_;
}

}
}

#endif
