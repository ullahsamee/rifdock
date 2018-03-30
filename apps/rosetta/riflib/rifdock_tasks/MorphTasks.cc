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

#include <riflib/rifdock_tasks/MorphTasks.hh>
#include <riflib/rifdock_tasks/HSearchTasks.hh>
#include <riflib/rifdock_tasks/UtilTasks.hh>
#include <riflib/rifdock_tasks/HackPackTasks.hh>
#include <riflib/rifdock_tasks/SeedingPositionTasks.hh>



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



void
create_dive_pop_hsearch_task( 
    std::vector<shared_ptr<Task>> & task_list, RifDockData & rdd ) {



    task_list.push_back(make_shared<DiversifyBySeedingPositionsTask>()); // this is a no-op if there are no seeding positions
    task_list.push_back(make_shared<DiversifyByNestTask>( 0 ));
    task_list.push_back(make_shared<HSearchInit>( ));
    for ( int i = 0; i <= rdd.opt.dive_resl-1; i++ ) {
        task_list.push_back(make_shared<HSearchScoreAtReslTask>( i, i, rdd.opt.tether_to_input_position_cut ));

        if (rdd.opt.hack_pack_during_hsearch) {
            task_list.push_back(make_shared<SortByScoreTask>( ));
            task_list.push_back(make_shared<FilterForHackPackTask>( 1, rdd.packopts.pack_n_iters, rdd.packopts.pack_iter_mult ));
            task_list.push_back(make_shared<HackPackTask>( i, i, rdd.opt.global_score_cut )); }

        task_list.push_back(make_shared<HSearchFilterSortTask>( i, rdd.opt.beam_size / rdd.opt.DIMPOW2, rdd.opt.global_score_cut, i < rdd.opt.dive_resl-1 ));

        if (rdd.opt.dump_x_frames_per_resl > 0) {
            task_list.push_back(make_shared<DumpHSearchFramesTask>( i, i, rdd.opt.dump_x_frames_per_resl, rdd.opt.dump_only_best_frames, rdd.opt.dump_only_best_stride, 
                                                                    rdd.opt.dump_prefix + "_" + rdd.scaffold_provider->get_data_cache_slow(ScaffoldIndex())->scafftag + boost::str(boost::format("_dp0_resl%i")%i) )); }
        if ( i < rdd.opt.dive_resl-1 ) {
            task_list.push_back(make_shared<HSearchScaleToReslTask>( i, i+1, rdd.opt.DIMPOW2, rdd.opt.global_score_cut )); } } 

    task_list.push_back(make_shared<HSearchFinishTask>( rdd.opt.global_score_cut ));

    task_list.push_back(make_shared<HSearchScaleToReslTask>( rdd.opt.dive_resl-1, rdd.opt.pop_resl-1, rdd.opt.DIMPOW2, rdd.opt.global_score_cut ));

    if ( rdd.opt.match_this_pdb != "") {
        task_list.push_back(make_shared<FilterByRmsdToThisPdbTask>( rdd.opt.match_this_pdb, rdd.opt.match_this_rmsd )); }

    task_list.push_back(make_shared<TestMakeChildrenTask>( ));



    task_list.push_back(make_shared<HSearchInit>( ));
    for ( int i = rdd.opt.pop_resl-1; i <= rdd.RESLS.size()-1; i++ ) {
        task_list.push_back(make_shared<HSearchScoreAtReslTask>( i, i, rdd.opt.tether_to_input_position_cut ));

        if (rdd.opt.hack_pack_during_hsearch) {
            task_list.push_back(make_shared<SortByScoreTask>( ));
            task_list.push_back(make_shared<FilterForHackPackTask>( 1, rdd.packopts.pack_n_iters, rdd.packopts.pack_iter_mult ));
            task_list.push_back(make_shared<HackPackTask>( i, i, rdd.opt.global_score_cut )); }

        task_list.push_back(make_shared<HSearchFilterSortTask>( i, rdd.opt.beam_size / rdd.opt.DIMPOW2, rdd.opt.global_score_cut, i < rdd.RESLS.size()-1 ));

        if (rdd.opt.dump_x_frames_per_resl > 0) {
            task_list.push_back(make_shared<DumpHSearchFramesTask>( i, i, rdd.opt.dump_x_frames_per_resl, rdd.opt.dump_only_best_frames, rdd.opt.dump_only_best_stride, 
                                                                    rdd.opt.dump_prefix + "_" + rdd.scaffold_provider->get_data_cache_slow(ScaffoldIndex())->scafftag + boost::str(boost::format("_dp0_resl%i")%i) )); }


        if ( i < rdd.RESLS.size()-1 ) {
            task_list.push_back(make_shared<HSearchScaleToReslTask>( i, i+1, rdd.opt.DIMPOW2, rdd.opt.global_score_cut )); } }

    task_list.push_back(make_shared<HSearchFinishTask>( rdd.opt.global_score_cut ));



}




}}
