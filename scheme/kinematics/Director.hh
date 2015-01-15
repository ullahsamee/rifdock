#ifndef INCLUDED_kinematics_Director_HH
#define INCLUDED_kinematics_Director_HH

#include "scheme/types.hh"
#include "scheme/kinematics/Scene.hh"
#include "scheme/nest/MultiNest.hh"

#include <boost/any.hpp>

#include <vector>

/*
	NOTES

	the boost any stuff does have a performance impact...

	AS IS: set_scene rate: 1.26338e+06 / sec
	without new/del       ~1.75000e+06 / sec

	directly setting scene from NEST without boost::any or virtual func calls or multinest:
	set_scene rate: 1.86553e+07 / sec [       OK ] Director.basic_test (37 ms)


*/

namespace scheme {
namespace kinematics {

template<
	class _Position,
	class _BigIndex = uint64_t,
	class _Index = uint64_t
>
struct Director {
	
	typedef _Position Position;
	typedef _BigIndex BigIndex;
	typedef _Index Index;	
	typedef SceneBase<Position,Index> Scene;

	virtual
	void
	set_scene(
		BigIndex const & i,
		int resl,
		Scene & scene
	) const = 0;

	virtual BigIndex size(int resl) const = 0;

};


class SymElem {};

template<
	class _Position,
	class _Index = uint64_t
>
struct SceneTree {
	typedef SceneTree<_Position,_Index> THIS;
	typedef _Position Position;
	typedef _Index Index;
	typedef SceneBase<Position,Index> Scene;
	typedef nest::NestBase<Index> Nest;
	typedef shared_ptr<Nest> Nestp;
	typedef shared_ptr<SceneTree> SceneTreep;

	Position edge_xform_;
	std::vector<Nestp> position_nests_;
	std::vector<Nestp> sym_elem_nests_;	
	std::vector<Nestp> conformation_nests_;		

	std::vector<Index> bodies_;
	std::vector<SymElem> sym_elems_;

	std::vector<SceneTreep> children_;


	// mutable Position tmp_;

	SceneTree( Position const & edge_xform = Position::Identity() ) : edge_xform_( edge_xform ) {}

	void add_body(Index ibody) { bodies_.push_back(ibody); }
	void add_position_nest( Nestp nestp ) { position_nests_.push_back(nestp); }
	void add_child( SceneTreep child ) { children_.push_back(child); }

	void
	get_nests(
		std::vector<Nestp> & nests
	) const {
		nests.insert( nests.end(),     position_nests_.begin(),     position_nests_.end() );
		nests.insert( nests.end(),     sym_elem_nests_.begin(),     sym_elem_nests_.end() );		
		nests.insert( nests.end(), conformation_nests_.begin(), conformation_nests_.end() );
		BOOST_FOREACH( SceneTreep child, children_ ){
			child->get_nests( nests );
		}
	}

	void
	get_dofs( std::vector<boost::any> & dofs ) const {
		for( int ipos = 0; ipos < position_nests_.size(); ++ipos){
			// std::cerr << "create temporary Position" << std::endl;
			dofs.push_back( new Position );
			// Position *tmp = &tmp_;
			// dofs.push_back( tmp );
		}
		if(     sym_elem_nests_.size() ) throw std::logic_error("sym_elem_nests not implemented");
		if( conformation_nests_.size() ) throw std::logic_error("conformation nests not implemented");
		BOOST_FOREACH( SceneTreep child, children_ ){
			child->get_dofs( dofs );
		}

	}
	void
	clear_empty_dofs( std::vector<boost::any>::iterator idofs ) const {
		for( int ipos = 0; ipos < position_nests_.size(); ++ipos){
			// std::cerr << "delete temporary Position" << std::endl;
			delete boost::any_cast<Position*>(*idofs);
			++idofs;
		}
		if(     sym_elem_nests_.size() ) throw std::logic_error("sym_elem_nests not implemented");
		if( conformation_nests_.size() ) throw std::logic_error("conformation nests not implemented");
		BOOST_FOREACH( SceneTreep child, children_ ){
			child->clear_empty_dofs( idofs );
		}

	}

	void
	set_scene(
		Position const & parent_position,
		std::vector<boost::any>::const_iterator idof,
		Scene & scene
	) const {
		Position position( edge_xform_ * parent_position );
		for( int ipos = 0; ipos < position_nests_.size(); ++ipos){
			// std::cerr << "set_scene get nest xform" << std::endl;
			Position & nest_xform( *boost::any_cast<Position*>(*idof) );
			// std::cerr << "set_scene get nest xform DONE" << std::endl;			
			position = nest_xform * position;
			++idof;
		}
		if(     sym_elem_nests_.size() ) throw std::logic_error("sym_elem_nests not implemented");
		if( conformation_nests_.size() ) throw std::logic_error("conformation nests not implemented");
		BOOST_FOREACH( Index ibody, bodies_ ){
			scene.set_position(ibody,position);
		}
		BOOST_FOREACH( SceneTreep child, children_ ){
			child->set_scene( parent_position, idof, scene );
		}
	}

};

	// MultiNest(Nests const & nests)
	// void add_nest( Nestp nestp )

template<
	class _Position,
	class _BigIndex = uint64_t,
	class _Index = uint64_t
>
struct NestDirector : public Director<_Position,_BigIndex,_Index> {
	
	typedef _BigIndex BigIndex;
	typedef _Index Index;
	typedef _Position Position;
	// typedef SceneTree<Position,Index> SceneTree;
	typedef shared_ptr<SceneTree<Position,Index> > SceneTreep;
	typedef SceneBase<Position,Index> Scene;
	typedef nest::NestBase<Index> Nest;
	typedef shared_ptr<Nest> Nestp;
	typedef nest::MultiNest<Index,BigIndex> MultiNest;

	SceneTreep root_;
	MultiNest multinest_;

	NestDirector( SceneTreep root ) : root_(root) { init(); }

	void init(){
		std::vector<Nestp> nests;
		root_->get_nests( nests );
		multinest_.init( nests );
	}

	virtual
	void
	set_scene(
		BigIndex const & i,
		int resl,
		Scene & scene
	) const {
		std::vector<boost::any> tmp_dofs;
		// std::cerr << "get_dofs" << std::endl;
		root_->get_dofs( tmp_dofs );
		assert( tmp_dofs.size() == multinest_.size() );
		// std::cerr << "get_states" << std::endl;
		multinest_.get_states( i, resl, tmp_dofs );
		// std::cerr << *boost::any_cast<numeric::X1dim*>(tmp_dofs.front()) << std::endl;
		assert( root_ != 0 );
		// std::cerr << "set_scene" << std::endl;		
		root_->set_scene( Position::Identity(), tmp_dofs.begin(), scene );
		// std::cerr << "clear_empty_dofs" << std::endl;
		root_->clear_empty_dofs( tmp_dofs.begin() );
	}

	virtual BigIndex size(int resl) const { return multinest_.size(resl); }


private:
	NestDirector() {}
};


}
}

#endif
