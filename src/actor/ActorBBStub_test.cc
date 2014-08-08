#include <gtest/gtest.h>

#include "actor/ActorBBStub.hh"
#include "actor/ActorConcept_io.hh"
#include "objective/ObjectiveFunction.hh"
#include "objective/ObjectiveVisitor.hh"
#include "numeric/X1dim.hh"

#include <io/dump_pdb_atom.hh>

#include <boost/foreach.hpp>



#include <stdint.h>
#include <fstream>

#include <Eigen/Geometry>

namespace scheme { namespace actor { namespace test {

typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;
typedef Eigen::AngleAxis<double> AA;
typedef Eigen::Vector3d Vec;


TEST(ActorBBStub,basic){
	ActorBBStub< Xform > a(Vec(1,0,0),Vec(0,0,0),Vec(1,1,0));
	dump_pdb(std::cout,a);
}


}
}
}
