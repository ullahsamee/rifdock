#include <gtest/gtest.h>

#include "kinematics/Scene_io.hh"
#include "actor/ActorConcept_io.hh"
#include "numeric/X1dim.hh"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>


#include <stdio.h>
#include <stdlib.h>

namespace scheme { namespace kinematics { namespace test {

using std::cout;
using std::endl;
using boost::tie;
using numeric::X1dim;

	typedef actor::ActorConcept<X1dim,int> ADI;
	typedef actor::ActorConcept<X1dim,char> ADC;	
	typedef m::vector< ADI, ADC > Actors;

	typedef Conformation<Actors> Conformation;
	typedef Scene<Conformation,X1dim,size_t> Scene;
	typedef size_t Index;
	typedef std::pair<size_t,size_t> Index2;
	typedef std::pair<Index2,Index2> Index4;

template<class _Interaction>
struct SetVisitor{
	typedef _Interaction Interaction;
	std::set<Interaction> set_;
	void operator()(Interaction const & i, double =1.0){
		set_.insert(i);
	}
};
template<class _Interaction>
struct AsymSetVisitor : SetVisitor<_Interaction> {
	typedef m::false_ Symmetric;
};


void test_iterator_visitor_agree(std::vector<int> nbod, int nsym){
	Scene scene;
	for(int i = 1; i < nsym; ++i) scene.add_symframe(100*i);
	for(int i = 0; i < (int)nbod.size(); ++i){
		scene.add_body();
		for(int j = 0; j < nbod[i]; ++j){
			scene.mutable_conformation_asym(i).get<ADI>().push_back(ADI(i,j));
			if( j <= i )
				scene.mutable_conformation_asym(i).get<ADC>().push_back(ADC(i,(char)j));
		}
	}

	{ // Asym ADI
		std::set<ADI> iterset; BOOST_FOREACH( ADI a, scene.get_actors_asym<ADI>() ) iterset.insert(a);
		AsymSetVisitor<ADI> v; scene.visit(v);
		ASSERT_EQ( iterset, v.set_ );
	}
	{ // SYM ADI
		std::set<ADI> iterset; BOOST_FOREACH( ADI a, scene.get_actors<ADI>() ) iterset.insert(a);
		SetVisitor<ADI> v; scene.visit(v);
		ASSERT_EQ( iterset, v.set_ );
	}
	{ // Asym ADC
		std::set<ADC> iterset; BOOST_FOREACH( ADC a, scene.get_actors_asym<ADC>() ) iterset.insert(a);
		AsymSetVisitor<ADC> v; scene.visit(v);
		ASSERT_EQ( iterset, v.set_ );
	}
	{ // SYM ADC
		std::set<ADC> iterset; BOOST_FOREACH( ADC a, scene.get_actors<ADC>() ) iterset.insert(a);
		SetVisitor<ADC> v; scene.visit(v);
		ASSERT_EQ( iterset, v.set_ );
	}

	{
		typedef std::pair<ADI,ADI> I;
		std::set<I> iterset;
		BOOST_FOREACH( Index4 i4, scene.get_interactions<I>() )
			iterset.insert( scene.get_interaction_from_placeholder<I>(i4) );
		SetVisitor<I> v; scene.visit(v);
		// ASSERT_EQ(iterset.size(),v.set_.size());
		// std::set<I>::const_iterator i=iterset.begin(),j=v.set_.begin();
		// for(;i!=iterset.end();++i,++j) cout << *i << "   /   " << *j << endl;
		ASSERT_EQ( iterset, v.set_ );
	}
	{
		typedef std::pair<ADI,ADC> I;
		std::set<I> iterset;
		BOOST_FOREACH( Index4 i4, scene.get_interactions<I>() )
			iterset.insert( scene.get_interaction_from_placeholder<I>(i4) );
		SetVisitor<I> v; scene.visit(v);
		ASSERT_EQ( iterset, v.set_ );
	}
	{
		typedef std::pair<ADC,ADC> I;
		std::set<I> iterset;
		BOOST_FOREACH( Index4 i4, scene.get_interactions<I>() )
			iterset.insert( scene.get_interaction_from_placeholder<I>(i4) );
		SetVisitor<I> v; scene.visit(v);
		ASSERT_EQ( iterset, v.set_ );
	}
	{
		typedef std::pair<ADC,ADI> I;
		std::set<I> iterset;
		BOOST_FOREACH( Index4 i4, scene.get_interactions<I>() )
			iterset.insert( scene.get_interaction_from_placeholder<I>(i4) );
		SetVisitor<I> v; scene.visit(v);
		ASSERT_EQ( iterset, v.set_ );
	}

}

TEST(Scene,test_visit_against_iteration){

	test_iterator_visitor_agree(std::vector<int>(3,2),2);
	test_iterator_visitor_agree(std::vector<int>(3,4),2);
	test_iterator_visitor_agree(std::vector<int>(13,4),2);
	test_iterator_visitor_agree(std::vector<int>(2,21),7);
	test_iterator_visitor_agree(std::vector<int>(13,11),1);

}


}
}
}
