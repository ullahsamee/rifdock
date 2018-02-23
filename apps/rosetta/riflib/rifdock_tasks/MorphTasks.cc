// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/MorphTasks.hh>

#include <riflib/types.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>
#include <riflib/scaffold/MorphingScaffoldProvider.hh>

#include <core/import_pose/import_pose.hh>


#include <string>
#include <vector>



namespace devel {
namespace scheme {


shared_ptr<std::vector<SearchPoint>> 
FilterByRmsdToThisPdbTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    shared_ptr<std::vector<SearchPoint>> out_points_p = make_shared<std::vector<SearchPoint>>( );


    core::pose::Pose match_this = *core::import_pose::pose_from_file( pose_filename_ );

    ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow( ScaffoldIndex() );
    EigenXform scaff2match = find_xform_from_identical_pose_to_pose( *(sdc->scaffold_centered_p ), match_this, 1 );

    float redundancy_filter_rg = sdc->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );

    int count = 0;
    for( SearchPoint const & sp : *search_points ) {

        rdd.director->set_scene( sp.index, 0, *rdd.scene_minimal );

        EigenXform x = rdd.scene_minimal->position(1);
        EigenXform xdiff = scaff2match.inverse() * x;
        float xmag =  xform_magnitude( xdiff, redundancy_filter_rg );
        if ( xmag < rmsd_ ) {
            out_points_p->push_back( sp );
        }
    }


    return out_points_p;
}



shared_ptr<std::vector<SearchPoint>> 
TestMakeChildrenTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    using ::scheme::scaffold::TreeLimits;

    shared_ptr<std::vector<SearchPoint>> out_points_p = make_shared<std::vector<SearchPoint>>( );

    shared_ptr<MorphingScaffoldProvider> morph_provider = std::dynamic_pointer_cast<MorphingScaffoldProvider>(rdd.scaffold_provider);
    morph_provider->test_make_children( ScaffoldIndex(0, 0) );

    TreeLimits limits = morph_provider->get_scaffold_index_limits();
    if (limits.size() == 1) {
        std::cout << "No scaffolds generated for morph_rules selection" << std::endl;
        return out_points_p;
    }

    uint64_t num_scaffolds = limits[1];

    out_points_p->resize( num_scaffolds * search_points->size() );

    int added = 0;
    for ( uint64_t scaffno = 0; scaffno < num_scaffolds; scaffno++ ) {
        for( SearchPoint const & sp : *search_points ) {
            SearchPoint new_sp = sp;
            new_sp.index.scaffold_index = ScaffoldIndex(1, scaffno);
            (*out_points_p)[added++] = new_sp;
        }
    }


    return out_points_p;
}



}}
