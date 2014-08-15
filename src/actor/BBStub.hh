#ifndef INCLUDED_actor_ActorBBStub_HH
#define INCLUDED_actor_ActorBBStub_HH

#include <iostream>
#include "io/dump_pdb_atom.hh"
#include "chemical/stub.hh"
#include <Eigen/Geometry>

namespace scheme {
namespace actor {


// template<class Xform>
// struct ActorBBStub {
// 	typedef Xform Position;
// 	typedef Eigen::Matrix<typename Xform::Scalar,3,1> Point;
// 	Xform position_;
// 	char ss_,aa_;

// 	ActorBBStub() : position_(Xform::Identity()),ss_('*'),aa_('*') {}

// 	ActorBBStub(Position const & p, char ss='*', char aa='*') : position_(p),ss_(ss),aa_(aa) {}

// 	ActorBBStub(Point n, Point ca, Point c, char ss='*', char aa='*') : ss_(ss),aa_(aa) {
// 		chemical::backbone_stub(n,ca,c,position_);
// 	}

// 	ActorBBStub(ActorBBStub const & a,Position const & moveby){
// 		position_ = moveby*a.position();
// 		ss_ = a.ss_;
// 		aa_ = a.aa_;
// 	}

// 	void set_position( Position const & pos	){ position_ = pos; }

// 	Position const & position() const { return position_; }
	
// 	bool operator==(ActorBBStub const & o) const { return o.position_.isApprox(position_) && o.ss_==ss_; o.aa_==aa_; }
// };

// template<class Xform>
// std::ostream & operator<<(std::ostream & out,ActorBBStub<Xform> const & x){
// 	return out << "ActorBBStub( " << x.position() << " " << x.data_ << " )";
// }

// template<class Xform>
// void dump_pdb(std::ostream & out,ActorBBStub<Xform> const & a){
// 	io::dump_pdb_atom(out, "CA", a.position() * Eigen::Matrix<typename Xform::Scalar,3,1>(0,0,0) );
// 	io::dump_pdb_atom(out, "O" , a.position() * Eigen::Matrix<typename Xform::Scalar,3,1>(1.5,0,0) );	
// 	io::dump_pdb_atom(out, "NI", a.position() * Eigen::Matrix<typename Xform::Scalar,3,1>(0,1.5,0) );
// 	io::dump_pdb_atom(out, "N" , a.position() * Eigen::Matrix<typename Xform::Scalar,3,1>(0,0,1.5) );		
// }

}
}

#endif
