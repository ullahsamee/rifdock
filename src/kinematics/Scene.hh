#ifndef INCLUDED_kinematics_Body_HH
#define INCLUDED_kinematics_Body_HH

#include <types.hh>
#include <util/meta/util.hh>
#include <util/meta/InstanceMap.hh>
#include <util/StoragePolicy.hh>
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
	typedef ActorConcept<Position,Data> THIS;

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

	bool operator==(THIS const & o) const { return o.position_==position_ && o.data_==data_; }
};
template< class P, class D >
std::ostream & operator<<(std::ostream & out, ActorConcept<P,D> const & a){
	return out << "ActorConcept( " << a.position_ << ", " << a.data_ << " )";
}


template< class Actors, class ActorContainer, class _Xform >
struct ConformationConcept : util::meta::InstanceMap<Actors,ActorContainer> {
	typedef _Xform Xform;

	template<class Actor>
	void add_actor(Actor const & a){
		this->template get<Actor>().insert(this->template get<Actor>().end(),a);
	}
};


///@brief holds a Conformation pointer and an Xform
template<
	class Conformation,
	template<class> class StoragePolicy = util::StoreSharedPointer
>
struct Body : public StoragePolicy<Conformation> {
	typedef StoragePolicy<Conformation> Super;
	typedef typename Conformation::Xform Xform;	

	Body() : position_(), Super() {}
	Body(Conformation & c) : position_(), Super(&c) {}

	template<class Actor> 
	Actor
	get_actor(
		size_t i
	) const { 
		return Actor(this->value()->template get<Actor>().at(i),position_);
	}
	template<class Actor> 
	Actor
	get_actor_unsafe(size_t i) const { 
		return Actor(this->value()->template get<Actor>()[i],position_); 
	}

	Xform const & xform() const { return position_; }
	Conformation const & conformation() const { return this->value(); }
 private:
	Xform position_;
 };
template< class C >
std::ostream & operator<<(std::ostream & out, Body<C> const & b){
	return out << "Body( " << b.xform() << ", " << b.conformation_ptr() << " )";
}

template<class Interaction>
struct SceneIter {

};

template<
	class _Conformation,
	class _Interactions,
	template<class> class _ConfStoragePolicy = util::StoreSharedPointer
>
struct Scene {
		
	typedef _Conformation Conformation;
	typedef _Interactions Interactions;	
	typedef Body<Conformation,_ConfStoragePolicy> Body;
	typedef typename Conformation::Keys Actors;
	// typedef typename Conformation::Xform Xform;	
	typedef std::vector<Body> Bodies;

	Bodies bodies_;

	Bodies const & bodies() const { return bodies_; }
	Body const & body       (size_t i) const { return bodies_.at(i); }
	Body const & body_unsafe(size_t i) const { return bodies_[i]; }
	void add_body(Conformation const & c){ bodies_.push_back(c); }
	void add_body(Body const & b){ bodies_.push_back(b); }

	///////// interactions

	template<class Interaction>
	struct has_interaction : m::bool_<false> {}; // f::result_of::has_key<MAP,Interaction>::value> {};

	template<class Interaction>
	struct interactions_type { 
		typedef std::pair<SceneIter<Interaction>,SceneIter<Interaction> > type; 
	};
	template<class Interaction>
	typename interactions_type<Interaction>::type const &
	get_interactions() const;

	template<class Interaction>
	struct interaction_placeholder_type { 
		typedef typename SceneIter<Interaction>::value_type type; 
	};
	template<class Interaction>
	Interaction const &
	get_interaction(
		typename interaction_placeholder_type<Interaction>::type const & placeholder
	) const;



};
template<class C,class I,template<class>class B>
std::ostream & operator<<(std::ostream & out, Scene<C,I,B> const & scene){ 
	using std::endl;
	typedef Scene<C,I,B> Scene;
	out << "Scene" << endl;
	out << "  Nbodies: " << scene.bodies().size() << endl;
	BOOST_FOREACH(typename Scene::Body const & b, scene.bodies()){ out << "    " << b << endl; }
	out << "  Types:" << endl;
	out << "    Actors:" << endl; util::meta::print_type<typename Scene::Actors>(out,"        ");
	out << "    Conformation:" << endl; util::meta::print_type<typename Scene::Conformation>(out,"        ");	
	return out;
}

}
}

#endif
