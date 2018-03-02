// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/SetFaModeTasks.hh>

#include <riflib/types.hh>


#include <string>
#include <vector>



namespace devel {
namespace scheme {


shared_ptr<std::vector<SearchPoint>> 
SetFaModeTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
SetFaModeTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
SetFaModeTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
SetFaModeTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    global_set_fa_mode( fa_mode_, rdd );

    return any_points;
}
    

// don't call this. Only to be used in the rif_dock_test test
void
global_set_fa_mode( bool fa_mode, RifDockData & rdd ) {

    rdd.scaffold_provider->set_fa_mode(fa_mode);

    RifDockScaffoldDirector scaffold_director(rdd.scaffold_provider, 1 );
    
    scaffold_director.set_scene( RifDockIndex(), 0, *rdd.scene_minimal );

    for ( ScenePtr & scene : rdd.scene_pt ) {
        scaffold_director.set_scene( RifDockIndex(), 0, *scene );
    }

}


}}
