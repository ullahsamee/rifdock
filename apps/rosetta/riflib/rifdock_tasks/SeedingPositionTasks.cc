// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/SeedingPositionTasks.hh>

#include <riflib/types.hh>
#include <rif_dock_test.hh>

#include <riflib/scaffold/ScaffoldDataCache.hh>
#include <riflib/scaffold/util.hh>

#include <core/import_pose/import_pose.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {


shared_ptr<std::vector<SearchPoint>> 
DiversifyBySeedingPositionsTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
DiversifyBySeedingPositionsTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
DiversifyBySeedingPositionsTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
DiversifyBySeedingPositionsTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {


    uint64_t num_positions = rdd.director->size(0, RifDockIndex()).seeding_index;
    if ( num_positions == 0 ) {
        return any_points;
    }

    shared_ptr<std::vector<AnyPoint>> diversified = make_shared<std::vector<AnyPoint>>( num_positions * any_points->size() );

    uint64_t added = 0;
    for ( AnyPoint const & pt : *any_points ) {

        for ( uint64_t i = 0; i < num_positions; i++ ) {
            (*diversified)[added] = pt;
            (*diversified)[added].index.seeding_index = i;
            added++;
        }
    }

    any_points->clear();

    return diversified;

}
    

shared_ptr<std::vector<EigenXform>>
setup_seeding_positions( RifDockOpt & opt, ProtocolData & pd, ScaffoldProviderOP & scaffold_provider ) {

    shared_ptr<std::vector<EigenXform>> seeding_positions = make_shared<std::vector<EigenXform>>();

    if ( opt.seed_with_these_pdbs.size() > 0 ) {
        if ( opt.seed_include_input ) {
            std::cout << "Adding seeding position: " << "_SP_input" << std::endl;
            seeding_positions->push_back(EigenXform::Identity());
            pd.seeding_tags.push_back("_input");
        }

        core::pose::Pose const & input = *scaffold_provider->get_data_cache_slow(ScaffoldIndex())->scaffold_unmodified_p;

        for ( std::string const & pdb : opt.seed_with_these_pdbs ) {
            std::string tag = "_SP_" + pdb_name( pdb );
            std::cout << "Adding seeding position: " << tag << std::endl;
            core::pose::Pose pose;
            core::import_pose::pose_from_file( pose, pdb );
            std::cout << pose.size() << " " << input.size() << std::endl;
            EigenXform xform = find_xform_from_identical_pose_to_pose( input, pose, 1 );
            seeding_positions->push_back(xform);
            pd.seeding_tags.push_back( tag );
        }

        return seeding_positions;
    }


    return nullptr;

}




}}

