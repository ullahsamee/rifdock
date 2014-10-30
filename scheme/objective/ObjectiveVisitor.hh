#ifndef INCLUDED_objective_ObjectiveVisitor_HH
#define INCLUDED_objective_ObjectiveVisitor_HH

#include "scheme/util/meta/util.hh"

namespace scheme {
namespace objective {

namespace objv_impl {
	SCHEME_MEMBER_TYPE_DEFAULT_TEMPLATE(Result,double)
}

template<class Objective,class Config>
struct ObjectiveVisitor {
	typedef typename Objective::Interaction Interaction;
	typedef typename objv_impl::get_Result_double<Objective>::type Result;
	Objective const & objective_;
	Config const & config_;
	Result result_;
	void clear() { result_ = 0; }
	ObjectiveVisitor( Objective const & o, Config const & c) : objective_(o),config_(c),result_() {}
	void operator()( Interaction const & i, double const & weight = 1.0 ){
		Result r = objective_.template operator()(i,config_);
		result_ += weight * r;
	}
};


}
}

#endif
