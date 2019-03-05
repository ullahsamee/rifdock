// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#ifndef INCLUDED_riflib_scaffold_ExtraScaffoldData_hh
#define INCLUDED_riflib_scaffold_ExtraScaffoldData_hh


#include <scheme/types.hh>

#include <riflib/HSearchConstraints.hh>

#include <vector>


namespace devel {
namespace scheme {

struct ExtraScaffoldData {
    ExtraScaffoldData () :
        csts(),
        force_scaffold_center( Eigen::Vector3f {std::numeric_limits<double>::quiet_NaN(), 0, 0} ),
        rotboltz_data_p( nullptr )
    {}


    std::vector<CstBaseOP> csts;
    Eigen::Vector3f force_scaffold_center;
    std::shared_ptr< std::vector< std::vector<float> > > rotboltz_data_p;
};


}
}

#endif