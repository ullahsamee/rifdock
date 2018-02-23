// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_ScaffoldProviderFactory_hh
#define INCLUDED_riflib_scaffold_ScaffoldProviderFactory_hh


#include <scheme/types.hh>
#include <riflib/scaffold/Baseline9AScaffoldProvider.hh>
#include <riflib/scaffold/MorphingScaffoldProvider.hh>
#include <riflib/scaffold/SingleFileScaffoldProvider.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_tasks/MorphTasks.hh>
#include <riflib/rifdock_tasks/HSearchTasks.hh>
#include <riflib/rifdock_tasks/UtilTasks.hh>
#include <riflib/rifdock_tasks/HackPackTasks.hh>


namespace devel {
namespace scheme {


ScaffoldProviderOP
get_scaffold_provider( 
        uint64_t iscaff,
        shared_ptr< RotamerIndex > rot_index_p, 
        RifDockOpt const & opt,
        MakeTwobodyOpts const & make2bopts,
        ::devel::scheme::RotamerRFTablesManager & rotrf_table_manager,
        bool & needs_scaffold_director ) {


    if (opt.scaff_search_mode == "default") {

        needs_scaffold_director = false;
        return make_shared<SingleFileScaffoldProvider>(
                iscaff,
                rot_index_p,
                opt,
                make2bopts,
                rotrf_table_manager);

    } else if (opt.scaff_search_mode == "morph_dive_pop") {

        needs_scaffold_director = true;
        return make_shared<MorphingScaffoldProvider>(
                iscaff,
                rot_index_p,
                opt,
                make2bopts,
                rotrf_table_manager);

    } else if (opt.scaff_search_mode == "nineA_baseline") {

        needs_scaffold_director = true;
        return make_shared<Baseline9AScaffoldProvider>(
                iscaff,
                rot_index_p,
                opt,
                make2bopts,
                rotrf_table_manager);
    } else {
        utility_exit_with_message("Error!!!! -rif_dock:scaff_search_mode " + opt.scaff_search_mode + " does not name a ScaffoldProvider!");
    }

}


void
create_dive_pop_hsearch_task( 
    std::vector<shared_ptr<Task>> & task_list, RifDockData & rdd ) {



    task_list.push_back(make_shared<DiversifyByNestTask>( 0 ));
    task_list.push_back(make_shared<HSearchInit>( ));
    for ( int i = 0; i <= rdd.opt.dive_resl-1; i++ ) {
        task_list.push_back(make_shared<HSearchScoreAtReslTask>( i, rdd.opt.tether_to_input_position_cut ));

        if (rdd.opt.hack_pack_during_hsearch) {
            task_list.push_back(make_shared<SortByScoreTask>( ));
            task_list.push_back(make_shared<FilterForHackPackTask>( 1, rdd.packopts.pack_n_iters, rdd.packopts.pack_iter_mult ));
            task_list.push_back(make_shared<HackPackTask>( i,  rdd.opt.global_score_cut )); }

        task_list.push_back(make_shared<HSearchFilterSortTask>( i, rdd.opt.beam_size / rdd.opt.DIMPOW2, rdd.opt.global_score_cut, i < rdd.opt.dive_resl-1 ));

        if (rdd.opt.dump_x_frames_per_resl > 0) {
            task_list.push_back(make_shared<DumpHSearchFramesTask>( i, rdd.opt.dump_x_frames_per_resl, rdd.opt.dump_only_best_frames, rdd.opt.dump_only_best_stride, 
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
        task_list.push_back(make_shared<HSearchScoreAtReslTask>( i, rdd.opt.tether_to_input_position_cut ));

        if (rdd.opt.hack_pack_during_hsearch) {
            task_list.push_back(make_shared<SortByScoreTask>( ));
            task_list.push_back(make_shared<FilterForHackPackTask>( 1, rdd.packopts.pack_n_iters, rdd.packopts.pack_iter_mult ));
            task_list.push_back(make_shared<HackPackTask>( i, rdd.opt.global_score_cut )); }

        task_list.push_back(make_shared<HSearchFilterSortTask>( i, rdd.opt.beam_size / rdd.opt.DIMPOW2, rdd.opt.global_score_cut, i < rdd.RESLS.size()-1 ));

        if (rdd.opt.dump_x_frames_per_resl > 0) {
            task_list.push_back(make_shared<DumpHSearchFramesTask>( i, rdd.opt.dump_x_frames_per_resl, rdd.opt.dump_only_best_frames, rdd.opt.dump_only_best_stride, 
                                                                    rdd.opt.dump_prefix + "_" + rdd.scaffold_provider->get_data_cache_slow(ScaffoldIndex())->scafftag + boost::str(boost::format("_dp0_resl%i")%i) )); }


        if ( i < rdd.RESLS.size()-1 ) {
            task_list.push_back(make_shared<HSearchScaleToReslTask>( i, i+1, rdd.opt.DIMPOW2, rdd.opt.global_score_cut )); } }

    task_list.push_back(make_shared<HSearchFinishTask>( rdd.opt.global_score_cut ));



}




}}

#endif