#ifndef INCLUDED_rosetta_objective_EtableParams_init_hh
#define INCLUDED_rosetta_objective_EtableParams_init_hh

#include "rosetta/score/EtableParams.hh"
#include <vector>

namespace scheme { namespace rosetta { namespace score {

///@brief horrible function to fill horrible rosetta datastructure of LJ/LK params
void init_EtableParams(
	std::vector<EtableParamsOnePair<float> > & analytic_parameters
);

}}}


#endif
