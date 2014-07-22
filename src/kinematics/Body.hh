#ifndef INCLUDED_kinematics_Body_HH
#define INCLUDED_kinematics_Body_HH

#include <util/meta/util.hh>
#include <util/meta/InstanceMap.hh>
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

struct ActorConcept {
	typedef double Position;
	Position position_;
	void set_position(Position const & pos){ position_ = pos; }
	Position const & position() const { return position_; }

};

template<class Actors>
struct BodyZero : util::meta::InstanceMap<Actors,std::vector<m::_1> > {

};

template<
	class Actors,
	class Position
>
struct Body : BodyZero<Actors> {
	typedef Body<Actors,Position> THIS;
	boost::shared_ptr<BodyZero<Actors> > reference_;
	Position position_;

};

}
}

#endif
