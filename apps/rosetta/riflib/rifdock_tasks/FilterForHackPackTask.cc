// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/FilterForHackPackTask.hh>

#include <riflib/types.hh>


#include <string>
#include <vector>



namespace devel {
namespace scheme {

std::vector<SearchPointWithRots>
FilterForHackPackTask::return_search_point_with_rotss( 
    std::vector<SearchPointWithRots> & search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) {


    size_t n_packsamp = 0;
    for( n_packsamp; n_packsamp < search_point_with_rotss.size(); ++n_packsamp ){
        if( search_point_with_rotss[n_packsamp].score > 0 ) break;
    }
    
    pd.npack = std::min( n_packsamp, (size_t)(pd.total_search_effort *
        ( rdd.opt.hack_pack_frac / (rdd.packopts.pack_n_iters*rdd.packopts.pack_iter_mult)) ) );

    search_point_with_rotss.resize( pd.npack );

    return search_point_with_rotss;
}



}}
