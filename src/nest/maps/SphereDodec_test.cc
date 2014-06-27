
#include <nest/maps/SphereDodec.hh>
#include <gtest/gtest.h>

namespace scheme {
namespace nest {
namespace maps {

TEST( SphereDodec, DISABLED_cell_centers ){
	typedef SphereDodec<2> MapType;
	MapType sd;
	MapType::ValueType val;
	sd.params_to_value( MapType::Params(0.5,0.5), 0, val );
}


}
}
}
