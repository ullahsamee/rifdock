#ifndef INCLUDED_kinematics_Body_HH
#define INCLUDED_kinematics_Body_HH

#include <util/meta/util.hh>
#include <util/meta/InstanceMap.hh>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>

namespace scheme {
namespace kinematics {

namespace m = boost::mpl;
namespace f = boost::fusion;

/* Requirements
	selective position update
	neighbor list
	parametric bodies
*/

template<class _Position, class _Data>
struct ActorConcept {

	typedef _Position Position;
	typedef _Data Data;

	Position position_;
	Data data_;

	ActorConcept() : position_(), data_() {}

	ActorConcept(Position const & p, Data const & d) : position_(p), data_(d) {}

	ActorConcept(
		ActorConcept const & actor0,
		Position const & moveby
	){
		data_ = actor0.data_;
		set_position(moveby*actor0.position());
	}

	void 
	set_position(
		Position const & pos
	){ position_ = pos; }

	Position const &
	position() const { return position_; }
};
template< class P, class D >
std::ostream & operator<<(std::ostream & out, ActorConcept<P,D> const & a){
	return out << "ActorConcept( " << a.position_ << ", " << a.data_ << " )";
}


template<class Actors>
struct Conformation : util::meta::InstanceMap<Actors,std::vector<m::_1> > {

};

template<class Actors,class Position>
struct Body {
	boost::shared_ptr<Conformation<Actors> > conformation_;
	Position position_;

	Body( boost::shared_ptr<Conformation<Actors> > cp ) : conformation_(cp) {}

	template<class Actor>
	Actor get_positioned_actor(size_t i) const {
		Actor const & actor0 = conformation_->template get<Actor>()[i];
		return Actor(actor0,position_);
	}

	void set_position(Position const & p) { position_ = p; }
	Position const & position() const { return position_; }	
};

template<class T>
struct SceneIter {

};

template<
	class _Actors,
	class _Interactions = void,
	class _Position = double
>
struct Scene {
		
	typedef _Actors Actors;
	typedef _Interactions Interactions;	
	typedef _Position Position;	

	template<class Interaction>
	struct has_interaction : m::bool_<false> {}; // f::result_of::has_key<MAP,Interaction>::value> {};

	template<class Interaction>
	struct interaction_placeholder_type { 
		typedef typename SceneIter<Interaction>::value_type type; 
	};
	template<class Interaction>
	struct interactions_type { 
		typedef std::pair<SceneIter<Interaction>,SceneIter<Interaction> > type; 
	};

	typedef std::vector<Body<Actors,Position> > Bodies;
	Bodies bodies_;
	Bodies const & bodies() const { return bodies_; }
	Bodies       & bodies()       { return bodies_; }	


	template<class Interaction>
	typename interactions_type<Interaction>::type const &
	get_interactions() const;

	template<class Interaction>
	Interaction const &
	get_interaction(
		typename interaction_placeholder_type<Interaction>::type const & placeholder
	) const;



};
template<class A,class I>
std::ostream & operator<<(std::ostream & out, Scene<A,I> const & scene){ 
	typedef Scene<A,I> Scene;
	out << "Scene";
	out << "    nbodies: " << scene.bodies().size() << std::endl;
	util::meta::PrintType(out).template operator()<typename Scene::Actors>(typename Scene::Actors());
	return out;
}

}
}

#endif
