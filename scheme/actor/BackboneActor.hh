#ifndef INCLUDED_actor_BackboneActor_HH
#define INCLUDED_actor_BackboneActor_HH

#include <Eigen/Dense>

namespace scheme {
namespace actor {




	template<class _Position>
	struct BackboneActor {

		///@brief Position type, leave out to make actor "Fixed"
		typedef _Position Position;
		typedef BackboneActor<Position> THIS;

		Position position_;
		char aa_,ss_;
		int index_;

		BackboneActor() : position_(), aa_('-'), ss_('-'), index_(0) {}

		BackboneActor(Position const & p, char aa, char ss, int i=0) : position_(p), aa_(aa), ss_(ss),index_(i) {}

		template <class V>
		BackboneActor(
			V n, 
			V ca, 
			V c, 
			char aa='-', 
			char ss='-', 
			int i=0
		) : aa_(aa), ss_(ss), index_(i) {
			from_n_ca_c(n,ca,c);
		}

		BackboneActor(
			BackboneActor const & actor0,
			Position const & moveby
		){
			aa_ = actor0.aa_;
			ss_ = actor0.ss_;
			index_ = actor0.index_;			
			set_position(moveby*actor0.position());
		}

		template <class Vin>
		void from_n_ca_c( Vin _n, Vin _ca, Vin _c ){
			typedef Eigen::Vector3d V;
			V n (  _n[0],  _n[1],  _n[2] );
			V ca( _ca[0], _ca[1], _ca[2] );
			V c (  _c[0],  _c[1],  _c[2] );									
			V e1 = n - ca;
			e1.normalize();
			V e3 = e1.cross( c - ca );
			e3.normalize();
			V e2 = e3.cross( e1 );
			Eigen::Matrix3d m;
			m(0,0) = e1[0];   m(0,1) = e2[0];   m(0,2) = e3[0];
			m(1,0) = e1[1];   m(1,1) = e2[1];   m(1,2) = e3[1];
			m(2,0) = e1[2];   m(2,1) = e2[2];   m(2,2) = e3[2];			

			V t = m * V(-0.865810,-1.764143,1.524857) + ca; // average 'CEN' icoor

			// std::cout << e1 << std::endl;
			// std::cout << e2 << std::endl;
			// std::cout << e3 << std::endl;						

			position_.setIdentity();
			position_.rotate( m );
			position_.translation()[0] = t[0];
			position_.translation()[1] = t[1];
			position_.translation()[2] = t[2];						
	
	// from_four_points(c,u,v,w)
	// 	Vector e1( u - v);
	// e1.normalize();
	// Vector e3( cross( e1, w - v ) );
	// e3.normalize();
	// Vector e2( cross( e3,e1) );
	// R.col_x( e1 ).col_y( e2 ).col_z( e3 );
	// t = c;	
		}

		void 
		set_position(
			Position const & pos
		){ position_ = pos; }

		Position const &
		position() const { return position_; }

		bool operator==(THIS const & o) const {
			return o.position_==position_ && 
			       o.aa_      ==aa_       &&
			       o.ss_      ==ss_       &&
			       o.index_   ==index_;
		}

		///@brief necessary for testing only
		// bool operator<(THIS const & o) const { return std::make_pair(position_,aa_) < std::make_pair(o.position_,o.aa_); }

		///@brief necessary for testing only
	  	template<class Archive> void serialize(Archive & ar, const unsigned int ){
	  		ar & position_;
	  		ar & aa_;
	  		ar & ss_;	  		
	  		ar & index_;
	  	}

	};

template<class X>
std::ostream & operator<<(std::ostream & out,BackboneActor<X> const& a){
	return out << "BackboneActor " << a.ss_ << " " << a.aa_ << " " << a.index_;
}

}
}

#endif
