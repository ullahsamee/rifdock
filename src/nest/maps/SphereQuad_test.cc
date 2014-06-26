
#include <nest/maps/SphereQuad.hh>
#include <gtest/gtest.h>

namespace scheme {
namespace nest {
namespace maps {

using std::cout;
using std::endl;

TEST( SphereDodec, cell_centers ){
	typedef SphereQuad<2,RowVector3d> MapType;
	MapType sd;
	MapType::Params prm;
	MapType::ValueType val;
	size_t celli;

	// sd.params_to_value( MapType::Params(0.5,0.5), 0, val ); std::cout << val << std::endl;
	// sd.params_to_value( MapType::Params(0.5,0.5), 1, val ); std::cout << val << std::endl;
	// sd.params_to_value( MapType::Params(0.5,0.5), 2, val ); std::cout << val << std::endl;
	// sd.params_to_value( MapType::Params(0.5,0.5), 3, val ); std::cout << val << std::endl;
	// sd.params_to_value( MapType::Params(0.5,0.5), 4, val ); std::cout << val << std::endl;
	// sd.params_to_value( MapType::Params(0.5,0.5), 5, val ); std::cout << val << std::endl;					

	sd.value_to_params( MapType::ValueType(0,0,1), prm, celli ); std::cout << prm << " " << celli << std::endl;

}


}
}
}
