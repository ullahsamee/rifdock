// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/types.hh>
#include <riflib/scaffold/SingleFileScaffoldProvider.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>



namespace devel {
namespace scheme {


SingleFileScaffoldProvider::SingleFileScaffoldProvider() {}


ParametricSceneConformationCOP 
SingleFileScaffoldProvider::get_scaffold(uint64_t i) {
    if ( ! conformation_ ) {
        utility_exit_with_message("Conformation not intialized yet!!");
    }
    return conformation_;
}


uint64_t 
SingleFileScaffoldProvider::get_scaffold_index_limits() {
    return 1;
}



}}

