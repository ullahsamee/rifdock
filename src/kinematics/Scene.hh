#ifndef INCLUDED_kinematics_Scene_HH
#define INCLUDED_kinematics_Scene_HH

#include <types.hh>
#include <util/meta/util.hh>
#include <util/meta/util.hh>
#include <util/meta/print_type.hh>
#include <util/meta/InstanceMap.hh>
#include <util/container/ContainerInteractions.hh>

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>

#include <boost/iterator/iterator_facade.hpp>

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

/////////////////////// Actor Concept //////////////////////////////
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
		return out << "Actor( " << a.position_ << ", " << a.data_ << " )";
	}
	template< class P >
	std::ostream & operator<<(std::ostream & out, ActorConcept<P,int> const & a){
		return out << "ADI( " << a.position_ << ", " << a.data_ << " )";
	}
	template< class P >
	std::ostream & operator<<(std::ostream & out, ActorConcept<P,char> const & a){
		return out << "ADC( " << a.position_ << ", " << a.data_ << " )";
	}

/////////////////// Conformation ////////////////////////////

	template< class ActorContainer >
	struct Conformation : util::meta::ContainerInstanceMap<ActorContainer> {
		template<class Actor>
		void add_actor(Actor const & a){
			this->template get<Actor>().insert(this->template get<Actor>().end(),a);
		}
	};

////////////////////// Body ///////////////////////////////////

	///@brief holds a Conformation pointer and an Xform
	template<
		class _Conformation,
		class _Position
	>
	struct Body {
		typedef _Conformation Conformation;
		typedef _Position Position;	

		Body() : position_(), conformation_(NULL) {}
		Body(shared_ptr<Conformation const> c) : position_(),conformation_(c) {}

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
	template< class C, class P >
	std::ostream & operator<<(std::ostream & out, Body<C,P> const & b){
		return out << "Body( " << b.xform() << ", " << &b.conformation() << " )";
	}



///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Scene Iterator default: 1body /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

	template<class Scene, class Interaction>
	struct SceneIter1B : boost::iterator_facade<
							SceneIter1B<Scene,Interaction>,
							std::pair<typename Scene::Index,typename Scene::Index> const,
							boost::forward_traversal_tag,
							std::pair<typename Scene::Index,typename Scene::Index> const
						>					
	{
		typedef SceneIter1B<Scene,Interaction> THIS;
		typedef Interaction Actor;
		typedef typename Scene::Index Index;
		SceneIter1B(){}
		SceneIter1B(Scene const & s, Index ib, Index ia) : scene_(&s),ibody_(ib),iactor_(ia){}
		static THIS make_begin(Scene const & s){ return THIS(s,next_nonempty(s,(Index)0),(Index)0); }
		static THIS make_end  (Scene const & s){ return THIS(s,       s.nbodies()       ,(Index)0); }
		static Index next_nonempty(Scene const & s, Index i) {
			// cout << "next_nonempty " << i;
			while( i < s.nbodies_asym() && 0==s.body(i).conformation().template get<Actor>().size() ) ++i;
			// cout << " -> " << i << endl;
			return i;
		}
	private:
	    friend class boost::iterator_core_access;
		void increment(){
            assert( ibody_ < scene_->nbodies() );
			++iactor_;
			if( iactor_ == scene_->conformation_unsafe(ibody_).template get<Actor>().size() ){
				iactor_ = 0;
				ibody_ = next_nonempty(*scene_,++ibody_);
			}
		}
		std::pair<Index,Index> const dereference() const { return std::make_pair(ibody_,iactor_);	}
		bool equal(THIS const & o) const { return scene_==o.scene_ && ibody_==o.ibody_ && iactor_==o.iactor_;	}
		Scene const * scene_;
		Index ibody_,iactor_;
	};


	template<class I> struct  make_pairpair { typedef std::pair<std::pair<I,I>,std::pair<I,I> > type; };


///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// SceneIter specialization 2body hetero //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

	struct CountPairNoDuplicates {
		template<class Index>
		static
		void
		next(Index & i, Index & j, Index size){
			++j;
			j += i==j ? 1 : 0; // skip iff j- == i
			if( j >= size ){
				j = 0;
				++i;
			}			
		}
		static size_t get_i_end(size_t n){ return n; }
		static size_t get_j_end(size_t  ){ return 0; }
	};
	struct CountPairUpperTriangle {
		template<class Index>
		static
		void
		next(Index & i, Index & j, Index size){
			++j;
			if( j >= size ){
				++i;
				j = i+1;
			}			
		}
		static size_t get_i_end(size_t n){ if(n==1) return 1; return n==0?0:n-1; }
		static size_t get_j_end(size_t n){ if(n==1) return 2; return n==0?0:n; }
	};

	template<class Scene,class Interaction, class CountPair>
	struct SceneIter2B : boost::iterator_facade<
							SceneIter2B<Scene,std::pair<
								typename Interaction::first_type,
								typename Interaction::second_type>,
								CountPair >,
							typename make_pairpair<typename Scene::Index>::type const,
							boost::forward_traversal_tag,
							typename make_pairpair<typename Scene::Index>::type const
						>{
		typedef typename Interaction::first_type Actor1;
		typedef typename Interaction::second_type Actor2;
		typedef SceneIter2B<Scene,std::pair<Actor1,Actor2>,CountPair> THIS;
		typedef typename Scene::Index Index;
		typedef typename Scene::Conformation Conf;
		typedef typename f::result_of::value_at_key<Conf,Actor1>::type Container1;
		typedef typename f::result_of::value_at_key<Conf,Actor2>::type Container2;
		typedef util::container::ContainerInteractions<typename Scene::Position,Container1,Container2,Index> ContInter;
		typedef typename ContInter::Range ContRange;

		SceneIter2B(){}
		SceneIter2B(Scene const & s, Index ib1, Index ib2) : scene_(&s),ibody1_(ib1),ibody2_(ib2){}
		void
		update_range(){ // X * A = B, X = B * ~A
			ContInter::get_interaction_range(
				scene_->body_unsafe(ibody2_).position() * inverse( scene_->body_unsafe(ibody1_).position() ) ,
				scene_->conformation_unsafe(ibody1_).template get<Actor1>()  ,
				scene_->conformation_unsafe(ibody2_).template get<Actor2>()  ,
				range_  );
			iter_ = get_cbegin(range_);
			end_  = get_cend(range_);
			// cout << "UPDATE_RANGE " << ibody1_ << " " << ibody2_ << " beg " << *iter_ << " end " << *end_ << endl;
		}
		static
		THIS 
		make_begin(
			Scene const & s
		){
			Index i=0,j=s.nbodies()>1?1:0;
			next_while_empty(s,i,j); 
			THIS beg(s,i,j); 
			// cout << "MAKE BEGIN " << i << " " << j << " " << *beg << endl;
			if(i < (Index)s.nbodies_asym() && j < (Index)s.nbodies() ) beg.update_range();
			return beg;
		}
		static
		THIS
		make_end(
			Scene const & s
		){ 
			return THIS(s,
				(Index)CountPair::get_i_end(s.nbodies_asym()) ,
				(Index)CountPair::get_j_end(s.nbodies_asym()) );
		}
		static
		void
		next_while_empty(
			Scene const & s,
			Index & i,
			Index & j
		) {
			// cout << "next_while_empty " << i<<"-"<<j << endl;
 			while( i < s.nbodies_asym() && j < s.nbodies() && ( 
 				0==s.body(i).conformation().template get<Actor1>().size() ||
			    0==s.body(j).conformation().template get<Actor2>().size() )
			){
				// cout << "  NEXT " << i << " " << j;
 				CountPair::next(i,j,s.nbodies());
 				// cout << " -> " << i << " " << j << endl;
 			}
			// cout << " -> " << i << "-" << j << endl;
		}
	private:
	    friend class boost::iterator_core_access;
		void
		increment(){
			// cout << "INCREMENT iter_ " << *iter_;
			++iter_;
			// cout << " -> " << *iter_ << " end " << *end_ << endl;
			if( iter_==end_ ){ // end_ of body/body inter, must update range_
				CountPair::next(ibody1_,ibody2_,scene_->nbodies());
				// cout << "NEXT_NONEMPTY " << ibody1_ << " " << ibody2_;
				next_while_empty(*scene_,ibody1_,ibody2_);
				// cout << " -> " << ibody1_ << " " << ibody2_ << endl;
				if( ibody1_ < scene_->nbodies_asym() && ibody2_ < scene_->nbodies() ) update_range();
				// std::cout << "redo range_ " << ibody1_ << " " << ibody2_ << std::endl;
				// std::cout << "SceneIter2B NEW BODIES: " << ibody1_ << " " << ibody2_ << std::endl;
			}
			// not valid when at end
			// assert( ibody1_ < scene_->nbodies_asym() );
			// assert( ibody2_ < scene_->nbodies() );
			// assert( iter_->first  < scene_->body(ibody1_).conformation().template get<Actor1>().size() );
			// assert( iter_->second < scene_->body(ibody2_).conformation().template get<Actor2>().size() );
		}
		typename make_pairpair<Index>::type const
		dereference() const { 
			return std::make_pair(std::make_pair(ibody1_,ibody2_),*iter_);
		}
		bool
		equal(THIS const & o) const { 
			// hacky but works for checking end_
			// std::cout << "EQUAL " << ibody1_ << " " << ibody2_ << " / " << o.ibody1_ << " " << o.ibody2_ << std::endl;
			return scene_==o.scene_ && ibody1_==o.ibody1_ && ibody2_==o.ibody2_;
		}
		Scene const * scene_;
		Index ibody1_,ibody2_;
		ContRange range_;
		typename util::container::get_citer<ContRange>::type iter_,end_;
	};

///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Scene / //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

	template<class T,class I> struct get_placeholder_type {
		typedef typename std::pair<I,I> type; };
	template<class T1, class T2,class I> struct get_placeholder_type<std::pair<T1,T2>,I> { 
		typedef typename std::pair<std::pair<I,I>,std::pair<I,I> > type; };

	template<class Scene,class Actor1>
	void
	get_scene_interaction_from_placeholder(
		Scene const & scene,
		typename get_placeholder_type<Actor1,typename Scene::Index>::type const & p,
		Actor1 & out
	) {
		typedef typename Scene::Body Body;
		Body const & body1 = scene.body(p.first );
		Actor1 const & a_0 = body1.conformation().template get<Actor1>()[p.second];
		Actor1 a( a_0, body1.position() );
		out = a;
	}
	template<class Scene,class Actor1,class Actor2>
	void
	get_scene_interaction_from_placeholder(
		Scene const & scene,
		typename get_placeholder_type<std::pair<Actor1,Actor2>,typename Scene::Index>::type const & p,
		std::pair<Actor1,Actor2> & out
	) {
		typedef typename Scene::Body Body;
		Body const & body1 = scene.body(p.first.first );
		Body const & body2 = scene.body(p.first.second);
		Actor1 const & a1_0 = body1.conformation().template get<Actor1>()[p.second.first];
		Actor2 const & a2_0 = body2.conformation().template get<Actor2>()[p.second.second];		
		Actor1 a1( a1_0, body1.position() );
		Actor2 a2( a2_0, body2.position() );
		out = std::make_pair(a1,a2);
	}

	template<
		class _Conformation,
		class _Position,
		class _Index = size_t
	>
	struct Scene {
		
		typedef _Conformation Conformation;
		typedef _Position Position;
		typedef Body<Conformation,Position> Body;
		typedef Scene<Conformation,Position> THIS;
		typedef std::vector<Body> Bodies;
		typedef _Index Index;

		Bodies bodies_;
		std::vector<Position> symframes_; // w/o identity

		Bodies const & bodies() const { return bodies_; }
		Body const & body       (size_t i) const { return bodies_.at(i); }
		Body const & body_unsafe(size_t i) const { return bodies_[i]; }
		Conformation const & conformation       (size_t i) const { return bodies_.at(i).conformation(); }
		Conformation const & conformation_unsafe(size_t i) const { return bodies_[i].conformation(); }
		Conformation & mutable_conformation(size_t i) { return const_cast<Conformation&>(bodies_[i].conformation()); }
		void add_body(shared_ptr<Conformation const> const & c){ bodies_.push_back(c); }
		void add_body(Body const & b){ bodies_.push_back(b); }
		void set_body(size_t i, shared_ptr<Conformation const> const & c){ bodies_[i] = c; }
		void set_body(size_t i, Body const & b){ bodies_[i] = b; }

		size_t nbodies() const { return bodies_.size()*(symframes_.size()+1); }
		size_t nbodies_asym() const { return bodies_.size(); }

		///////// interactions

		template<class Interaction>
		struct has_interaction : m::bool_<false> {}; // f::result_of::has_key<MAP,Interaction>::value> {};

		template<class Interaction>
		struct iter_type { typedef typename 
			m::if_< util::meta::is_pair<Interaction>,
				typename m::if_< util::meta::is_homo_pair<Interaction>,
					SceneIter2B<THIS,Interaction,CountPairUpperTriangle>,
					SceneIter2B<THIS,Interaction,CountPairNoDuplicates>
				>::type,
				SceneIter1B<THIS,Interaction>
			>::type
			type; };

		template<class Interaction>
		struct interactions_type { 
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

		template<class Interaction>
		Interaction
		get_interaction(
			typename iter_type<Interaction>::type const & iter
		) const {
			Interaction i;
			get_scene_interaction_from_placeholder(*this,*iter,i);
			return i;
		}

		template<class Interaction>
		Interaction
		get_interaction(
			typename get_placeholder_type<Interaction,Index>::type const & p
		) const {
			Interaction i;
			get_scene_interaction_from_placeholder(*this,p,i);
			return i;
		}



	};
	template<class C,class P>
	std::ostream & operator<<(std::ostream & out, Scene<C,P> const & scene){ 
		using std::endl;
		typedef Scene<C,P> Scene;
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
