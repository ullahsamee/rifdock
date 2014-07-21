#ifndef INCLUDED_kinematics_Body_HH
#define INCLUDED_kinematics_Body_HH

#include <util/meta/util.hh>
#include <util/meta/InstanceMap.hh>
#include <boost/foreach.hpp>

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

struct ActorConcept {
	typedef double Position;
	
};


template<
	class Actors,
	class Interactions
>
struct Body {
	typedef Body<Actors,Interactions> THIS;
};

}
}

#endif
