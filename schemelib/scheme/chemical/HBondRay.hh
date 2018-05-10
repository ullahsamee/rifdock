#ifndef INCLUDED_chemical_HBondRay_HH
#define INCLUDED_chemical_HBondRay_HH

#include <Eigen/Dense>

namespace scheme { namespace chemical {

struct HBondRay {
	::Eigen::Vector3f horb_cen, direction;
	int32_t id=1, group=-1;
	bool operator==(HBondRay const & other) const {
		float d1 = (horb_cen-other.horb_cen).norm();
		float d2 = (direction-other.direction).norm();
		return d1 < 0.0001 && d2 < 0.0001;
	}

	template< class Xform >
	void apply_xform( Xform const & xform ) {
		Eigen::Vector3f dirpos = horb_cen + direction;
		horb_cen = xform * horb_cen;
		direction = xform * dirpos - horb_cen;
	}
};


}}

#endif
