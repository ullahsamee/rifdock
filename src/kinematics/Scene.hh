#ifndef INCLUDED_kinematics_Scene_HH
#define INCLUDED_kinematics_Scene_HH

#include "types.hh"
#include "kinematics/SceneIterator.hh"

#include "util/meta/util.hh"
#include "util/meta/InstanceMap.hh"
#include "util/container/ContainerInteractions.hh"

#include <boost/foreach.hpp>

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <vector>

namespace scheme {
namespace kinematics {

namespace m = boost::mpl;
namespace f = boost::fusion;
using std::cout;
using std::endl;
/* Requirements
	selective position update
	neighbor list
	parametric bodies
	smart containers
*/


template< class ActorContainer >
struct Conformation : util::meta::ContainerInstanceMap<ActorContainer> {
	typedef typename util::meta::ContainerInstanceMap<ActorContainer>::Keys Actors;
	template<class Actor>
	void add_actor(Actor const & a){
		this->template get<Actor>().insert(this->template get<Actor>().end(),a);
	}
};

namespace impl {
	///@brief holds a Conformation pointer and an Xform
	template<
		class _Conformation,
		class _Position
	>
	struct BodyTplt {
		typedef _Conformation Conformation;
		typedef _Position Position;	

		BodyTplt() : position_(), conformation_(NULL) {}
		BodyTplt(shared_ptr<Conformation const> c) : position_(),conformation_(c) {}

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

		Position const & position() const { return position_; }
		void set_position(Position const & p){ position_ = p; }

		Conformation const & conformation() const { return *conformation_; }

	 private:
		Position position_;
		shared_ptr<Conformation const> conformation_;
	 };

}

///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Scene / //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

	template<class T,class I> struct get_placeholder_type {
		typedef typename std::pair<I,I> type; };
	template<class T1, class T2,class I> struct get_placeholder_type<std::pair<T1,T2>,I> { 
		typedef typename std::pair<std::pair<I,I>,std::pair<I,I> > type; };

	template<class Scene,class Actor>
	void
	get_scene_interaction_from_placeholder(
		Scene const & scene,
		typename get_placeholder_type<Actor,typename Scene::Index>::type const & p,
		Actor & out
	) {
		assert( scene.__conformation_unsafe__(p.first).template get<Actor>().size() > p.second );
		Actor const & a_0 = scene.__conformation_unsafe__(p.first).template get<Actor>()[p.second];
		Actor a( a_0, scene.__position_unsafe__(p.first) );
		out = a;
	}
	template<class Scene,class Actor1,class Actor2>
	void
	get_scene_interaction_from_placeholder(
		Scene const & scene,
		typename get_placeholder_type<std::pair<Actor1,Actor2>,typename Scene::Index>::type const & p,
		std::pair<Actor1,Actor2> & out
	) {
		assert( scene.__conformation_unsafe__(p.first.first  ).template get<Actor1>().size() > p.second.first );
		assert( scene.__conformation_unsafe__(p.first.second ).template get<Actor2>().size() > p.second.second );
		Actor1 const & a1_0 = scene.__conformation_unsafe__(p.first.first ).template get<Actor1>()[p.second.first ];
		Actor2 const & a2_0 = scene.__conformation_unsafe__(p.first.second).template get<Actor2>()[p.second.second];
		Actor1 a1( a1_0, scene.__position_unsafe__(p.first.first ) );
		Actor2 a2( a2_0, scene.__position_unsafe__(p.first.second) );
		// cout << "get_from_ph Actor1 " << a1_0 << " " << a1 << endl;
		// cout << "get_from_ph Actor2 " << a2_0 << " " << a2 << endl;		
		out.first  = a1;
		out.second = a2;
	}

	template<
		class _Conformation,
		class _Position,
		class _Index = size_t
	>
	struct Scene {
		
		typedef _Conformation Conformation;
		typedef _Position Position;
		typedef impl::BodyTplt<Conformation,Position> Body;
		typedef Scene<Conformation,Position> THIS;
		typedef typename Conformation::Actors Actors;
		typedef std::vector<Body> Bodies;
		typedef _Index Index;
		typedef m::true_ DefinesInteractionWeight;

		Bodies bodies_;
		std::vector<Position> symframes_; // w/o identity
		size_t n_sym_bodies_;

		void update(){
			n_sym_bodies_ = bodies_.size()*(symframes_.size()+1);
		}
		Index sym_index_map(Index & i) const {
			Index isym = i / bodies_.size();
			i = i % bodies_.size();
			return isym;
		}

		Scene(size_t nbodies=0) { for(size_t i=0; i<nbodies; ++i) add_body(); update(); }

		void set_symmetry(std::vector<Position> const & sym){ symframes_ = sym; update(); }
		void add_symframe(Position const & symframe){ symframes_.push_back(symframe); update(); }
		std::vector<Position> const & symframes() const { return symframes_; }


		/// mutators
		void add_body(){ bodies_.push_back( make_shared<Conformation const>() ); update(); }
		void add_body(shared_ptr<Conformation const> const & c){ bodies_.push_back(c); update(); }
		void add_body(Body const & b){ bodies_.push_back(b); update(); }
		void set_body(size_t i, shared_ptr<Conformation const> const & c){ bodies_[i] = c; }
		void set_body(size_t i, Body const & b){ bodies_[i] = b; }
		Conformation & mutable_conformation_asym(size_t i) { return const_cast<Conformation&>(bodies_.at(i).conformation()); }


		Conformation const & conformation(size_t i) const { return bodies_.at(i).conformation(); }
		// Body const & body       (size_t i) const { return bodies_.at(i); }
		Body const & position(size_t i) const { return bodies_.at(i).position(); }

		/// mainly internal use
		Bodies const & __bodies_unsafe__() const { return bodies_; }
		std::vector<Position> & __symframes_unsafe__(){ return symframes_; }		

		Position const & __position_unsafe_asym__(size_t i) const { return bodies_[i].position(); }
		Position       & __position_unsafe_asym__(size_t i)       { return bodies_[i].position(); }
		Position const __position_unsafe__(size_t i) const { 
			Index isym = sym_index_map(i);
			if( isym == 0 )	return bodies_[i].position();
			else return symframes_[isym-1] * bodies_[i].position();
		}
		Conformation const & __conformation_unsafe_asym__(size_t i) const { return bodies_[i].conformation(); }
		Conformation const & __conformation_unsafe__(size_t i) const { return bodies_[i%bodies_.size()].conformation(); }

		size_t nbodies() const { return n_sym_bodies_; }
		size_t nbodies_asym() const { return bodies_.size(); }

		//////////////////////////// actors //////////////////////////////////

		template<class Actor>
		std::pair<
			SceneIter1B<THIS,Actor,Symmetric,ActorCopy>,
			SceneIter1B<THIS,Actor,Symmetric,ActorCopy>
		>
		get_actors() const {
			typedef SceneIter1B<THIS,Actor,Symmetric,ActorCopy> Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}
		template<class Actor>
		std::pair<
			SceneIter1B<THIS,Actor,NotSymmetric,ActorCopy>,
			SceneIter1B<THIS,Actor,NotSymmetric,ActorCopy>
		>
		get_actors_asym() const {
			typedef SceneIter1B<THIS,Actor,NotSymmetric,ActorCopy> Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}
		template<class Actor>
		std::pair<
			SceneIter1B<THIS,Actor,Symmetric,PlaceHolder>,
			SceneIter1B<THIS,Actor,Symmetric,PlaceHolder>
		>
		get_actors_placeholder() const {
			typedef SceneIter1B<THIS,Actor,Symmetric,PlaceHolder> Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}
		template<class Actor>
		std::pair<
			SceneIter1B<THIS,Actor,NotSymmetric,PlaceHolder>,
			SceneIter1B<THIS,Actor,NotSymmetric,PlaceHolder>
		>
		get_actors_placeholder_asym() const {
			typedef SceneIter1B<THIS,Actor,NotSymmetric,PlaceHolder> Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}


		/////////////////////////// interactions //////////////////////////////////

		// TODO: fix this!
		template<class Interaction> struct
		has_interaction : m::bool_<true> {}; // f::result_of::has_key<MAP,Interaction>::value> {};

		template<class Interaction>	struct
		interaction_placeholder_type {typedef typename get_placeholder_type<Interaction,Index>::type type; };

		template<class Interaction> struct
		iter_type { typedef typename 
			m::if_< util::meta::is_pair<Interaction>,
				typename m::if_< util::meta::is_homo_pair<Interaction>,
					SceneIter2B<THIS,Interaction,CountPairUpperTriangle>,
					SceneIter2B<THIS,Interaction,CountPairNoDuplicates>
				>::type,
				SceneIter1B<THIS,Interaction,NotSymmetric,PlaceHolder>
			>::type
			type; };

		template<class Interaction> struct
		interactions_type { 
			typedef typename iter_type<Interaction>::type Iter;
			typedef std::pair<Iter,Iter> type; 
		};
		template<class Interaction>
		typename interactions_type<Interaction>::type
		get_interactions() const {
			typedef typename iter_type<Interaction>::type Iter;
			Iter beg = Iter::make_begin(*this);
			Iter end = Iter::make_end  (*this);
			return std::make_pair(beg,end);
		}

		///@brief gets a concrete interaction (Actor or pair of Actors) from a SceneIter
		template<class Iter>
		typename Iter::Interaction
		get_interaction_from_iter(
			Iter const & iter
		) const {
			typename Iter::Interaction i;
			get_scene_interaction_from_placeholder(*this,*iter,i);
			// cout << "get_interaction_from_iter " << i << endl;
			return i;
		}

		///@brief gets a concrete interaction (Actor or pair of Actors) from a placeholder
		template<class Interaction>
		Interaction
		get_interaction_from_placeholder(
			typename get_placeholder_type<Interaction,Index>::type const & p
		) const {
			Interaction i;
			get_scene_interaction_from_placeholder(*this,p,i);
			return i;
		}

		///@brief default weight is 1.0
		template<class PlaceHolder>
		double get_weight_from_placeholder(PlaceHolder const &) const {
			return 1.0;
		}
		///@brief for symmetry, if either body not in asym unit, weight should be 0.5
		double get_weight_from_placeholder(std::pair<std::pair<Index,Index>,std::pair<Index,Index> > const & ph ) const {
			return ( ph.first.first < nbodies_asym() && ph.first.second < nbodies_asym() ) ? 1.0 : 0.5;
		}

	};

}
}

#endif
