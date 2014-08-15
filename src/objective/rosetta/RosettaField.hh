#ifndef INCLUDED_objective_rosetta_RosettaField_HH
#define INCLUDED_objective_rosetta_RosettaField_HH

#include "objective/voxel/FieldCache.hh"
#include "objective/rosetta/AnalyticEvaluation.hh"
#include "objective/rosetta/EtableParams_init.hh"
#include "types.hh"
#include <vector>
#include <boost/foreach.hpp>

namespace scheme { namespace objective { namespace rosetta {

template<class Atom>
struct RosettaField {
	EtableParams<float> params;
	std::vector<Atom> atoms_;
	RosettaField() { init_EtableParams(params); }
	RosettaField(std::vector<Atom> const & atm) : atoms_(atm) { init_EtableParams(params); }
	void add_atom(Atom const & a) { atoms_.push_back(a); }
	float compute_rosetta_energy(float x, float y, float z, int atype) const {
		float atr=0,rep=0,sol=0;
		BOOST_FOREACH(Atom const & a,atoms_){
			EtableParamsOnePair<float> const & p = params.params_for_pair(a.type(),atype);
			float const dx = x-a.position()[0];
			float const dy = y-a.position()[1];
			float const dz = z-a.position()[2];
			float const dis2 = dx*dx+dy*dy+dz*dz;
			float const dis = std::sqrt(dis2);
			float const inv_dis2 = 1.0f/dis2;
			float atr0,rep0,sol0;
			lj_evaluation( p, dis, dis2, inv_dis2, atr0, rep0);
			lk_evaluation( p, dis, inv_dis2, sol0 ); 
			atr += atr0;
			rep += rep0;
			sol += sol0;
		}
		return 0.8*atr+0.44*rep+0.75*sol;
	}
};

template<class Atom>
struct RosettaFieldAtype : voxel::Field3D<float> {
	RosettaField<Atom> const & rf_;
	int atype_;
	RosettaFieldAtype(RosettaField<Atom> const & rf, int atype) : rf_(rf),atype_(atype) {}
	float operator()(float x, float y, float z) const {
		return rf_.compute_rosetta_energy(x,y,z,atype_);
	}
};


}}}

#endif
