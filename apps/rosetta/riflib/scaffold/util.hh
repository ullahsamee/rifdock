// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_scaffold_util_hh
#define INCLUDED_riflib_scaffold_util_hh


#include <scheme/types.hh>



namespace devel {
namespace scheme {


struct ScaffoldDataCache {
    std::vector<std::vector<float> > const * rotamer_energies_1b_;

};


typedef shared_ptr<ScaffoldDataCache> ScaffoldDataCacheOP;
typedef shared_ptr<ScaffoldDataCache const > ScaffoldDataCacheCOP;

}}



#endif