#ifndef INCLUDED_actor_BackboneSasaActor_HH
#define INCLUDED_actor_BackboneSasaActor_HH

#include <Eigen/Dense>
#include <scheme/io/dump_pdb_atom.hh>

#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>

namespace scheme {
namespace actor {


// sasa_points_ holds the position information


	struct BackboneSasaActor {


		typedef ::Eigen::Transform<float,3,Eigen::AffineCompact> EigenXform;

		///@brief Position type, leave out to make actor "Fixed"
		typedef EigenXform Position;
		typedef typename Position::Scalar Float;
		typedef BackboneSasaActor THIS;
		typedef Eigen::Matrix<Float,3,1> V3;

		std::vector<Eigen::Vector3f> sasa_points_;
		int index_;

		BackboneSasaActor() : index_(0)  {}

		BackboneSasaActor(std::vector<Eigen::Vector3f> const & sasa_points, int i=0) : 
			sasa_points_(sasa_points),
			index_(i) 
		{}


		BackboneSasaActor(
			BackboneSasaActor const & actor0,
			Position const & to_moveby
		){
			sasa_points_ = actor0.sasa_points_;
			index_ = actor0.index_;	

			moveby(to_moveby);
		}

		std::vector<Eigen::Vector3f> const & 
		sasa_points() const {
			return sasa_points_;
		}

		int
		index() const {
			return index_;
		}

		void 
		moveby(
			Position const & pos
		){ 
			for ( size_t ipos = 0; ipos < sasa_points_.size(); ipos++ ) {
				sasa_points_[ipos] = pos * sasa_points_[ipos];
			}
		}

		// Position const &
		// position() const { return Position::Identity(); }

		bool operator==(THIS const & o) const {
			return o.sasa_points_==sasa_points_ && 
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
void write_pdb( std::ostream & out, BackboneSasaActor const & a, MetaData const & meta ){

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
std::ostream & operator<<(std::ostream & out,BackboneSasaActor const& a){
	return out << "BackboneSasaActor " << a.index_;
}

}
}

#endif
