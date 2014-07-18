#ifndef INCLUDED_objective_ObjectiveFunction_HH
#define INCLUDED_objective_ObjectiveFunction_HH

#include <util/meta/util.hh>
#include <util/meta/InstanceMap.hh>

#include <boost/foreach.hpp>

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/sort.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/filter_view.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>

namespace scheme {
namespace objective {

namespace mpl = boost::mpl;
namespace bf = boost::fusion;


namespace traits {
	struct result_type { template<class T> struct apply { typedef typename T::Result type; }; };
	struct petals_type { template<class T> struct apply { typedef typename T::Petals type; }; };
}

template<class Petals>
struct objective_petals_equal {
	template<class Objective>
	struct apply : mpl::bool_< boost::is_same<Petals,typename Objective::Petals>::value > {};
};


template<class Objectives>
struct copy_objective_if_petals_equal {
	template<class Petals>
	struct apply : mpl::copy_if<Objectives,objective_petals_equal<Petals>, mpl::back_inserter<bf::vector<> > > {};
};

template<class T> struct wrap {};


template<
	class Petals,
	class PetalsSource,
	class Results,
	class Config
>
struct EvalObjective {
	PetalsSource const & petal_source;
	Results & results;
	Config const & config;

	EvalObjective(
		PetalsSource const & p,
		Results & r,
		Config const & c
	) : petal_source(p),results(r),config(c) {}

	template<class Objective> void operator()(Objective const & objective) const {
		BOOST_STATIC_ASSERT( bf::result_of::has_key<Results,Objective>::value );
		// std::cout << "    EvalObjective source: " << typeid(petal_source.template get<Petals>()).name() << std::endl;
		// std::cout << "               objective: " << typeid(objective).name() << std::endl;
		BOOST_FOREACH( Petals const & petals, petal_source.template get<Petals>() ){
			// std::cout << "        Petals: " << petals <<
								 // " Result: " << typeid(results.template get<Objective>()).name() << std::endl;
			objective(
				petals,
				results.template get<Objective>(),
				config
			);
		}
	}
};

template<
	class PetalsSource,
	class ObjectiveMap,
	class Results,
	class Config
>
struct EvalObjectives {
	PetalsSource const & petal_source;
	ObjectiveMap const & objective_map;
	Results & results;
	Config const & config;
	
	EvalObjectives(
		PetalsSource const & p,
		ObjectiveMap const & f,
		Results & r,
		Config const & c
	) : petal_source(p),objective_map(f),results(r),config(c) {}
	
	template<class Petals> void operator()(wrap<Petals>) const {
		BOOST_STATIC_ASSERT( bf::result_of::has_key<PetalsSource,Petals>::value );
		// std::cout << "EvalObjectives Petals: " << typeid(Petals).name() << std::endl;
		bf::for_each(
			objective_map.template get<Petals>(), 
			EvalObjective<
				Petals,
				PetalsSource,
				Results,
				Config
			>(petal_source,results,config) 
		);
	}
};

template<
	typename _Objectives,
	typename _Config = int
	>
struct ObjectiveFunctionImpl {
	typedef _Objectives Objectives;
	typedef typename mpl::copy<Objectives,mpl::inserter<mpl::set<>, mpl::insert<mpl::_1,mpl::_2> > >::type SET_Objectives;
	BOOST_STATIC_ASSERT( mpl::size<SET_Objectives>::value == mpl::size<Objectives>::value );

	typedef typename mpl::transform<Objectives, traits::petals_type >::type AllPetals;
	typedef typename mpl::copy<AllPetals,mpl::inserter<mpl::set<>, mpl::insert<mpl::_1,mpl::_2> > >::type SET_Petals;
	typedef typename mpl::copy<SET_Petals,mpl::back_inserter<mpl::vector<> > >::type Petals;

	typedef typename mpl::transform< 
		Petals, 
		copy_objective_if_petals_equal<Objectives>,
		mpl::back_inserter<mpl::vector<> >
	>::type UniquePetalsObjectives;
	BOOST_STATIC_ASSERT( mpl::size<UniquePetalsObjectives>::value == mpl::size<Petals>::value );

	typedef util::meta::InstanceMap< Petals, UniquePetalsObjectives > ObjectiveMap;
	
	typedef util::meta::NumericInstanceMap<Objectives,traits::result_type> Results;

	typedef _Config Config;
	
	ObjectiveMap objective_map_;
	Config config_;

	template<class Objective>
	Objective & get_objective(){
		BOOST_STATIC_ASSERT( true ); // shold check that objective is actuall in Objectives...
		return bf::deref(bf::find<Objective>( objective_map_.template get<typename Objective::Petals>() ));
	}
	// Results weights_;
	// ObjectiveFunctionImpl(){
	// 	weights_.setall(1.0);
	// }

	template<class PetalsSource>
	void operator()(PetalsSource const & source, Results & results) const {
		BOOST_STATIC_ASSERT( util::meta::is_InstanceMap<PetalsSource>::value );
		mpl::for_each< Petals, wrap<mpl::_1> >( 
			EvalObjectives<
				PetalsSource,
				ObjectiveMap,
				Results,
				Config
			>( source, objective_map_, results, config_ )
		);
	}
	template<class PetalsSource>
	Results operator()(PetalsSource const & source) const {
		Results results(0);
		this->operator()(source,results);
		return results;
	}
};

template<class O,class C>
struct ObjectiveFunction : ObjectiveFunctionImpl<O,C> {};

template<typename O,typename C>
std::ostream & operator<<(std::ostream & out, ObjectiveFunctionImpl<O,C> const & obj){
	typedef ObjectiveFunctionImpl<O,C> OBJ;
	out << "ObjectiveFunction" << std::endl;
	out << "    Petals Combos:" << std::endl;
		mpl::for_each<typename OBJ::Petals>(util::meta::PrintType(out,"        "));
	out << "    Objectives:" << std::endl;		
	bf::for_each( obj.objective_map_, util::meta::PrintBFMapofVec(out,"        ") );
	out << "    RAW:" << std::endl;
		mpl::for_each<typename OBJ::Objectives>(util::meta::PrintType(out,"        "));
	return out;
}

}
}

#endif
