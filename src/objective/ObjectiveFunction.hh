#ifndef INCLUDED_objective_ObjectiveFunction_HH
#define INCLUDED_objective_ObjectiveFunction_HH

#include <util/meta/util.hh>
#include <util/meta/intersect.hh>
#include <util/meta/InstanceMap.hh>

#include <boost/foreach.hpp>

#include <boost/mpl/eval_if.hpp>
// #include <boost/mpl/set.hpp>
// #include <boost/mpl/inserter.hpp>
// #include <boost/mpl/insert.hpp>
#include <boost/mpl/for_each.hpp>
// #include <boost/mpl/copy.hpp>
// #include <boost/mpl/copy_if.hpp>
#include <boost/mpl/transform.hpp>
// #include <boost/mpl/int.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>

// #define DEBUG_IO 

namespace scheme {
namespace objective {

namespace m = boost::mpl;
namespace f = boost::fusion;

namespace traits {
	struct result_type { template<class T> struct apply { typedef typename T::Result type; }; };
	struct interaction_type { template<class T> struct apply { typedef typename T::Interaction type; }; };
}

namespace impl {

	// TODO: is there any way these could return references?
	template< class Interaction, class PlaceHolder, class InteractionSource, bool>
	struct get_interaction_from_placeholder_impl {
		Interaction operator()(PlaceHolder const & placeholder, InteractionSource const &){
			return placeholder;
		}
	};
	template< class Interaction, class PlaceHolder, class InteractionSource>
	struct get_interaction_from_placeholder_impl<Interaction,PlaceHolder,InteractionSource,false> {
		Interaction operator()(PlaceHolder const & placeholder, InteractionSource const & source){
			return source.template get_interaction_from_placeholder<Interaction>(placeholder);
		}
	};
	template< class Interaction, class PlaceHolder, class InteractionSource >
	Interaction get_interaction_from_placeholder(PlaceHolder const & placeholder, InteractionSource const & source){
		return get_interaction_from_placeholder_impl<
			Interaction,PlaceHolder,InteractionSource,boost::is_same<PlaceHolder,Interaction>::value>()
			(placeholder,source);
	}

	///@brief helper class tests Objective::Interaction == Interaction
	template<class Interaction>
	struct objective_interactions_equal {
		template<class Objective>
		struct apply : m::bool_< boost::is_same<Interaction,typename Objective::Interaction>::value > {};
	};

	///@brief helper class builds mpl::vector of Objectives with Objective::Interaction == Interaction
	template<class Objectives>
	struct copy_objective_if_interaction_equal {
		template<class Interaction>
		struct apply : m::copy_if<Objectives,objective_interactions_equal<Interaction>, m::back_inserter<f::vector<> > > {};
	};

	///@brief helper functor to call objective for all Interaction instances in InteractionSource
	template<
		class Interaction,
		class Results,
		class Config
	>
	struct EvalObjective {
		Interaction const & interaction;
		Results & results;
		Config const & config;
		double weight;
		EvalObjective(
			Interaction const & i,
			Results & r,
			Config const & c,
			double w
		) : interaction(i),results(r),config(c),weight(w) {}

		template<class Objective> void operator()(Objective const & objective) const {
			BOOST_STATIC_ASSERT( f::result_of::has_key<typename Results::FusionType,Objective>::value );
			#ifdef DEBUG_IO
			std::cout << "    EvalObjective:     Objective " << Objective::name() <<"( " << interaction  << " )" << std::endl;
			#endif
			typedef typename f::result_of::value_at_key<Results,Objective>::type Result;
			/// TODO: unpack args to allow for tuples of std::ref, would avoid extra copy
			Result result = objective( interaction, config );
			results.template get<Objective>() += weight*result;
		}
	};

	using m::false_;
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(DefinesInteractionWeight,false_)
	// SCHEME_HAS_MEMBER_TYPE(DefinesInteractionWeight)

	template<class InteractionSource, class Placeholder>
	typename boost::disable_if<typename get_DefinesInteractionWeight_false_<InteractionSource>::type,double>::type
	get_interaction_weight(InteractionSource const & , Placeholder const & ){ 
		return 1.0;
	}
	template<class InteractionSource, class Placeholder>
	typename boost::enable_if<typename get_DefinesInteractionWeight_false_<InteractionSource>::type,double>::type
	get_interaction_weight(InteractionSource const & is, Placeholder const & ph){ 
		return is.get_weight_from_placeholder(ph);
	}

	///@brief helper functor to call each objective for Interaction
	template<
		class InteractionSource,
		class ObjectiveMap,
		class Results,
		class Config
	>
	struct EvalObjectives {
		InteractionSource const & interaction_source;
		ObjectiveMap const & objective_map;
		Results & results;
		Config const & config;
		
		EvalObjectives(
			InteractionSource const & p,
			ObjectiveMap const & f,
			Results & r,
			Config const & c
		) : interaction_source(p),objective_map(f),results(r),config(c) {}
		
		template<class Interaction> void operator()(util::meta::type2type<Interaction>) const {
			BOOST_STATIC_ASSERT( InteractionSource::template has_interaction<Interaction>::value );		
			#ifdef DEBUG_IO
				std::cout << "    EvalObjectives Interaction: ";
				util::meta::PrintType(std::cout).operator()(Interaction());
			#endif
			typedef typename InteractionSource::template interaction_placeholder_type<Interaction>::type Placeholder;
			BOOST_FOREACH(
				Placeholder const & interaction_placeholder,
				interaction_source.template get_interactions<Interaction>()
			){
				double weight = get_interaction_weight(interaction_source,interaction_placeholder);
				// std::cout << "        Interaction: " << interaction <<
									 // " Result: " << typeid(results.template get<Objective>()).name() << std::endl;
				Interaction const & interaction =
					get_interaction_from_placeholder<
							Interaction,
							Placeholder,
							InteractionSource
						>( interaction_placeholder, interaction_source );
				f::for_each(
					objective_map.template get<Interaction>(), 
					EvalObjective< Interaction, Results, Config >
					             ( interaction, results, config, weight )
				);
			}
		}
	};

	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(InteractionTypes,void)

	// get_interaction_weight

}

///@brief Minimal Concept for Objective used in ObectiveFunction
struct ObjectiveConcept {
	///@typedef required type of result this functor generates, 
	///@note must be convertable to double but could have extra data
	typedef double Result;
	///@typedef required type of input evaluatable data, could be a pair of tuple
	typedef int Interaction;
	///@brief return name of Objective
	static std::string name(){ return "ObjectiveConcept"; }
	///@brief evaluate the Objective
	///@param Interaction main body or interaction to evaluate
	///@param Result result should be stored here
	template<class Config>
	Result operator()( Interaction const& a, Config const& ) const { return a; }
};

///@brief a generic objective function for interacting bodies
///@tparam Objectives sequence of types modeling the Objective concept
///@tparam Config a global config object passed to each Objective
template<
	typename _Objectives,
	typename _Config = int
	>
struct ObjectiveFunction {

	typedef _Objectives Objectives;
	typedef _Config Config;

	///@typedef only for error checking, make sure Objectives are non-redundant
	typedef typename
		m::copy<
			Objectives,
			m::inserter<m::set<>, m::insert<m::_1,m::_2> >
		>::type
	UniqueObjectives;

	BOOST_STATIC_ASSERT( m::size<UniqueObjectives>::value == m::size<Objectives>::value );

	///@typedef unique InteractionTypes types
	typedef typename 
		m::copy< typename // copy back into mpl::vector;
			m::copy< typename // copy to mpl::set so unique
				m::transform< // all InteractionTypes Type from Obectives
					Objectives,
					traits::interaction_type
					>::type,
				m::inserter<m::set<>,m::insert<m::_1,m::_2> >
				>::type,
			m::back_inserter<m::vector<> >
		>::type
	InteractionTypes;

	///@typedef mpl::vector of Objectives for each unique PetalsType
	typedef typename
		m::transform< 
			InteractionTypes, 
			impl::copy_objective_if_interaction_equal<Objectives>,
			m::back_inserter<m::vector<> >
		>::type
	InteractionObjectives;

	///@typedef ObjectiveMap
	typedef util::meta::InstanceMap<
		InteractionTypes,
		InteractionObjectives
	> ObjectiveMap;

	///@typedef Results
	typedef util::meta::NumericInstanceMap<
		Objectives,
		traits::result_type
	> Results;

	ObjectiveMap objective_map_;
	Config default_config_;
	Results weights_;

	///@brief default c'tor, init weights_ to 1
	ObjectiveFunction() : weights_(1.0) {}

	///@brief accessor for Objectives, may be used to configure or initialize Objective instances
	template<class Objective>
	Objective &
	get_objective(){
		BOOST_STATIC_ASSERT( true ); // shold check that objective is actuall in Objectives...
		return f::deref( f::find<Objective>( 
			objective_map_.template get<typename Objective::Interaction>() 
		));
	}

	///@brief evaluate a InteractionSource with specified Config
	///@param InteractionSource must implement the InteractionSource concept
	///@param Config override default Config
	///@return a NumericInstanceMap of ResultTypes defined by Objectives
	template<class InteractionSource>
	Results 
	operator()(
		InteractionSource const & source,
		Config const & config
	) const {
		// make sure we only operate on interactions contained in source
		typedef typename impl::get_InteractionTypes_void<InteractionSource>::type SourceInteractionTypes;
		typedef typename m::eval_if<
				boost::is_same<void,SourceInteractionTypes>,
				InteractionTypes,
				util::meta::intersect<
					InteractionTypes,
					SourceInteractionTypes
				>
			>::type
			MutualInteractionTypes;
		Results results;
		BOOST_STATIC_ASSERT( m::size<MutualInteractionTypes>::value );
		#ifdef DEBUG_IO
			std::cout << "ObjectiveFunction operator()" << std::endl;
		#endif
		m::for_each< MutualInteractionTypes, util::meta::type2type<m::_1> >( 
			impl::EvalObjectives<
				InteractionSource,
				ObjectiveMap,
				Results,
				Config
			>( source, objective_map_, results, config )
		);
		return results * weights_;
	}

	///@brief evaluate a InteractionSource with default Config
	///@param InteractionSource must implement the InteractionSource concept
	///@return a NumericInstanceMap of ResultTypes defined by Objectives
	template<class InteractionSource>
	Results 
	operator()(
		InteractionSource const & source
	) const {
		return this->operator()<InteractionSource>(source,default_config_);
	}

};


}
}

#endif
