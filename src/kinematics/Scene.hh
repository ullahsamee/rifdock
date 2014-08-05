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
#include <boost/tuple/tuple.hpp>

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
		class _Position,
		class Index
	>
	struct BodyTplt {
		typedef _Conformation Conformation;
		typedef _Position Position;	

		BodyTplt() : position_(), conformation_(NULL) {}
		BodyTplt(shared_ptr<Conformation const> c) : position_(),conformation_(c) {}

		template<class Actor> 
		Actor
		get_actor(
			Index i
		) const { 
			return Actor(this->value()->template get<Actor>().at(i),position_);
		}
		template<class Actor> 
		Actor
		get_actor_unsafe(Index i) const { 
			return Actor(this->value()->template get<Actor>()[i],position_); 
		}

		Position const & position() const { return position_; }
		void set_position(Position const & p){ position_ = p; }

		Conformation const & conformation() const { return *conformation_; }

	 private:
		Position position_;
		shared_ptr<Conformation const> conformation_;
	 };

	using m::true_;
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(Symmetric,true_)

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
		
		typedef Scene<_Conformation,_Position,_Index> THIS;
		typedef _Conformation Conformation;
		typedef _Position Position;
		typedef _Index Index;
		typedef impl::BodyTplt<Conformation,Position,Index> Body;
		typedef typename Conformation::Actors Actors;
		typedef std::vector<Body> Bodies;
		typedef m::true_ DefinesInteractionWeight;
		typedef m::true_ UseVisitor;

		Bodies bodies_;
		std::vector<Position> symframes_; // w/o identity
		Index n_sym_bodies_;

		void update(){
			n_sym_bodies_ = (Index)bodies_.size()*((Index)symframes_.size()+1);
		}
		Index sym_index_map(Index & i) const {
			Index isym = i / bodies_.size();
			i = i % bodies_.size();
			return isym;
		}

		Scene(Index nbodies=0) { for(Index i=0; i<nbodies; ++i) add_body(); update(); }

		void set_symmetry(std::vector<Position> const & sym){ symframes_ = sym; update(); }
		void add_symframe(Position const & symframe){ symframes_.push_back(symframe); update(); }
		std::vector<Position> const & symframes() const { return symframes_; }


		/// mutators
		void add_body(){ bodies_.push_back( make_shared<Conformation const>() ); update(); }
		void add_body(shared_ptr<Conformation const> const & c){ bodies_.push_back(c); update(); }
		void add_body(Body const & b){ bodies_.push_back(b); update(); }
		void set_body(Index i, shared_ptr<Conformation const> const & c){ bodies_[i] = c; }
		void set_body(Index i, Body const & b){ bodies_[i] = b; }
		Conformation & mutable_conformation_asym(Index i) { return const_cast<Conformation&>(bodies_.at(i).conformation()); }


		Conformation const & conformation(Index i) const { return bodies_.at(i%bodies_.size()).conformation(); }
		// Body const & body       (Index i) const { return bodies_.at(i); }
		Position position(Index i) const {
			Index isym = sym_index_map(i);
			if( isym == 0 )	return              bodies_.at(i).position();
			else return symframes_.at(isym-1) * bodies_.at(i).position();
		}

		/// mainly internal use
		Bodies const & __bodies_unsafe__() const { return bodies_; }
		std::vector<Position> & __symframes_unsafe__(){ return symframes_; }		

		Position const & __position_unsafe_asym__(Index i) const { return bodies_[i].position(); }
		Position       & __position_unsafe_asym__(Index i)       { return bodies_[i].position(); }
		Position const __position_unsafe__(Index i) const { 
			Index isym = sym_index_map(i);
			if( isym == 0 )	return bodies_[i].position();
			else return symframes_[isym-1] * bodies_[i].position();
		}
		Conformation const & __conformation_unsafe_asym__(Index i) const { return bodies_[i].conformation(); }
		Conformation const & __conformation_unsafe__(Index i) const { return bodies_[i%bodies_.size()].conformation(); }

		Index nbodies() const { return n_sym_bodies_; }
		Index nbodies_asym() const { return bodies_.size(); }

		//////////////////////////// actor iteration //////////////////////////////////

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


		/////////////////////////// interaction iteration //////////////////////////////////

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

		///////////////////////////// Visitation ///////////////////////////////////////

		template<class Visitor>
		typename boost::disable_if<util::meta::is_pair<typename Visitor::Interaction> >::type
		visit(Visitor & visitor) const {
			typedef typename Visitor::Interaction Actor;
			typedef typename f::result_of::value_at_key<Conformation,Actor>::type Container;
			Index const NBOD = (Index)bodies_.size();
			bool const symmetric = impl::get_Symmetric_true_<Visitor>::type::value;
			Index const NSYM = symmetric ? (Index)symframes_.size()+1 : 1;
			for(Index i1 = 0; i1 < NBOD*NSYM; ++i1){
				Conformation const & c = conformation(i1);
				Position     const & p =     position(i1);
				Container const & container = c.template get<Actor>();
				BOOST_FOREACH( Actor const & a, container ){
					visitor( Actor(a,p) );
				}
			}

		}

		template<class Visitor>
		typename boost::enable_if<util::meta::is_pair<typename Visitor::Interaction> >::type
		visit(Visitor & visitor) const {
			typedef typename Visitor::Interaction::first_type Actor1;
			typedef typename Visitor::Interaction::second_type Actor2;			
			typedef typename f::result_of::value_at_key<Conformation,Actor1>::type Container1;
			typedef typename f::result_of::value_at_key<Conformation,Actor2>::type Container2;
			typedef util::container::ContainerInteractions<typename Scene::Position,Container1,Container2,Index> ContInter;
			typedef typename ContInter::Range ContRange;
			typename util::container::get_citer<ContRange>::type iter,end;
			ContRange range;

			Index const NBOD = (Index)bodies_.size();
			Index const NSYM = (Index)symframes_.size()+1;
			for(Index i1 = 0; i1 < NBOD*NSYM; ++i1){
				Conformation const & c1 = conformation(i1);
				Position     const & p1 =     position(i1);
				for(Index i2 = 0; i2 < NBOD*NSYM; ++i2){
					if( i1 >= NBOD && i2 >= NBOD ) continue;
					if( i1==i2 || (boost::is_same<Actor1,Actor2>::value && i2 <= i1) ) continue;
					Conformation const & c2 = conformation(i2);
					Position     const & p2 =     position(i2);

					// // simple version, ~10-15% faster in simple case
					// Index const NACT1 = c1.template get<Actor1>().size();
					// Index const NACT2 = c2.template get<Actor2>().size();				
					// for(Index j1 = 0; j1 < NACT1; ++j1){
					// 	Actor1 a1( c1.template get<Actor1>()[j1], p1 );
					// 	for(Index j2 = 0; j2 < NACT2; ++j2){
					// 		Actor2 a2( c2.template get<Actor2>()[j2], p2 );
					// 		visitor( std::make_pair(a1,a2), i1<NBOD&&i2<NBOD?1.0:0.5 );
					// 	}
					// }

					Container1 const & container1 = c1.template get<Actor1>();
					Container2 const & container2 = c2.template get<Actor2>();
					Position relative_position = __position_unsafe__(i2) * inverse( __position_unsafe__(i1) );
					ContInter::get_interaction_range( relative_position, container1, container2, range );
					for( iter = get_cbegin(range),end  = get_cend(range); iter != end; ++iter){
						Index j1,j2;
						boost::tie(j1,j2) = *iter;
						Actor1 a1( c1.template get<Actor1>()[j1], p1 );
						Actor2 a2( c2.template get<Actor2>()[j2], p2 );
						visitor( std::make_pair(a1,a2), i1<NBOD&&i2<NBOD?1.0:0.5 );
					}

				}
			}

		}


	};

}
}

#endif
