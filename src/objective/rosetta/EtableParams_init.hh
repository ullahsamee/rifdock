#ifndef INCLUDED_objective_rosetta_EtableParams_init_hh
#define INCLUDED_objective_rosetta_EtableParams_init_hh

#include "objective/rosetta/EtableParams.hh"
#include <vector>

namespace scheme { namespace objective { namespace rosetta {

///@brief horrible function to fill horrible rosetta datastructure of LJ/LK params
void init_EtableParams(
	std::vector<EtableParamsOnePair<float> > & analytic_parameters
);

}}}


#endif
