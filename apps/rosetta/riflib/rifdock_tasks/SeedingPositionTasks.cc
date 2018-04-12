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
#include <riflib/seeding_util.hh>

#include <core/import_pose/import_pose.hh>

#include <riflib/rifdock_tasks/HSearchTasks.hh>
#include <riflib/rifdock_tasks/SetFaModeTasks.hh>
#include <riflib/rifdock_tasks/HackPackTasks.hh>
#include <riflib/rifdock_tasks/RosettaScoreAndMinTasks.hh>
#include <riflib/rifdock_tasks/CompileAndFilterResultsTasks.hh>
#include <riflib/rifdock_tasks/OutputResultsTasks.hh>
#include <riflib/rifdock_tasks/UtilTasks.hh>

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
    


shared_ptr<std::vector<SearchPoint>> 
DiversifyByXformFileTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
DiversifyByXformFileTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
DiversifyByXformFileTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
DiversifyByXformFileTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {


    std::vector< std::pair< int64_t, EigenXform > > xform_positions;
    {
        runtime_assert_msg(parse_exhausitive_searching_file(file_name_, xform_positions /*, 10*/), "Faild to parse the xform file!!!");
    }

    uint64_t nest_size = xform_positions.size();

    shared_ptr<std::vector<AnyPoint>> diversified = make_shared<std::vector<AnyPoint>>( nest_size * any_points->size() );

    uint64_t added = 0;
    for ( AnyPoint const & pt : *any_points ) {

        for ( uint64_t i = 0; i < nest_size; i++ ) {
            (*diversified)[added] = pt;
            (*diversified)[added].index.nest_index = xform_positions[i].first;
            added++;
        }
    }

    any_points->clear();

    return diversified;
}




void
create_rifine_task( 
    std::vector<shared_ptr<Task>> & task_list, RifDockData & rdd ) {

    RifDockOpt const & opt = rdd.opt;
    int seeding_size = rdd.director->size(0, RifDockIndex()).seeding_index;
    int final_resl = rdd.rif_ptrs.size() - 1;

    task_list.push_back(make_shared<DiversifyBySeedingPositionsTask>()); 
    if ( rdd.opt.xform_fname != "ALL" ) {
        task_list.push_back(make_shared<DiversifyByXformFileTask>( rdd.opt.xform_fname ));
    } else {
        task_list.push_back(make_shared<DiversifyByNestTask>( 0 ));
    }
    
    task_list.push_back(make_shared<HSearchInit>( ));
    task_list.push_back(make_shared<HSearchScoreAtReslTask>( 0, final_resl, rdd.opt.tether_to_input_position_cut ));
    // task_list.push_back(make_shared<HSearchFinishTask>( rdd.opt.global_score_cut ));

    task_list.push_back(make_shared<FilterToBestNTask>( int( std::ceil ( opt.beam_size / seeding_size ) ), opt.filter_seeding_positions_separately, opt.filter_scaffolds_separately));
    task_list.push_back(make_shared<FilterByScoreCutTask>( opt.cluster_score_cut ));
    task_list.push_back(make_shared<FilterByFracTask>( opt.hack_pack_frac, 0, opt.filter_seeding_positions_separately, opt.filter_scaffolds_separately ));

    task_list.push_back(make_shared<SetFaModeTask>( true ));
    task_list.push_back(make_shared<FilterForHackPackTask>( 1, 0, 0 ));
    task_list.push_back(make_shared<HackPackTask>(  0, final_resl, 9e9 )); // hackpack sorts the results

    task_list.push_back(make_shared<FilterByScoreCutTask>( opt.cluster_score_cut ));
    task_list.push_back(make_shared<FilterByBiggestBlocksFracTask>( opt.keep_top_clusters_frac, opt.filter_seeding_positions_separately, opt.filter_scaffolds_separately, true ));
    task_list.push_back(make_shared<FilterByScoreCutTask>( opt.global_score_cut ));
    // task_list.push_back(make_shared<CompileAndFilterResultsTask>( 0, final_resl, 100000000, opt.redundancy_filter_mag, 0, 0, 
    //                                                                                   opt.filter_seeding_positions_separately, opt.filter_scaffolds_separately )); 
    task_list.push_back(make_shared<RemoveRedundantPointsTask>( opt.redundancy_filter_mag, 0, opt.filter_seeding_positions_separately, opt.filter_scaffolds_separately ));

    task_list.push_back(make_shared<DumpSeedingClusterScoreTask>());

    bool do_rosetta_score = opt.rosetta_score_fraction > 0;
    bool do_rosetta_min = do_rosetta_score && opt.rosetta_min_fraction > 0;

    if ( do_rosetta_score || opt.rosetta_filter_even_if_no_score ) {
                              task_list.push_back(make_shared<FilterByFracTask>( opt.rosetta_score_fraction, opt.rosetta_score_each_seeding_at_least, opt.filter_seeding_positions_separately, opt.filter_scaffolds_separately ));
        if (do_rosetta_score) task_list.push_back(make_shared<RosettaScoreTask>( 0, opt.rosetta_score_cut, do_rosetta_min, true));
    }

    if ( do_rosetta_min || opt.rosetta_filter_even_if_no_score ) {
                            task_list.push_back(make_shared<FilterForRosettaMinTask>( opt.rosetta_min_fraction, opt.rosetta_min_at_least ));
        if (do_rosetta_min) task_list.push_back(make_shared<RosettaMinTask>( 0, opt.rosetta_score_cut, false )); 
    }

    // dummy task to turn these into RifDockResults
    task_list.push_back(make_shared<CompileAndFilterResultsTask>( 0, final_resl, opt.n_pdb_out, opt.redundancy_filter_mag, 0, 0, 
                                                                                      opt.filter_seeding_positions_separately, opt.filter_scaffolds_separately )); 

    task_list.push_back(make_shared<OutputResultsTask>( 0, final_resl ));



}





}}

