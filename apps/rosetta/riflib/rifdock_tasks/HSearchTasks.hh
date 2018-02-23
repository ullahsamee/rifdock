// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_rifdock_tasks_HSearchTasks_hh
#define INCLUDED_riflib_rifdock_tasks_HSearchTasks_hh

#include <riflib/types.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/task/SearchPointTask.hh>
#include <riflib/task/AnyPointTask.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct DiversifyByNestTask : public AnyPointTask {

    DiversifyByNestTask(
        int resl
        ) :
        resl_( resl )
        {}

    shared_ptr<std::vector<SearchPoint>> 
    return_search_points( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

    shared_ptr<std::vector<SearchPointWithRots>> 
    return_search_point_with_rotss( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

    shared_ptr<std::vector<RifDockResult>> 
    return_rif_dock_results( 
        shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

private:
    template<class AnyPoint>
    shared_ptr<std::vector<AnyPoint>>
    return_any_points( 
        shared_ptr<std::vector<AnyPoint>> any_points, 
        RifDockData & rdd, 
        ProtocolData & pd ); // override

private:
    bool resl_;

};

struct HSearchInit : public SearchPointTask {

    HSearchInit() {}

    shared_ptr<std::vector<SearchPoint>> 
    return_search_points( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

};

struct HSearchScoreAtReslTask : public SearchPointTask {

    HSearchScoreAtReslTask(
        int resl,
        float tether_to_input_position_cut ) :
        resl_( resl ),
        tether_to_input_position_cut_( tether_to_input_position_cut )
        {}

    shared_ptr<std::vector<SearchPoint>> 
    return_search_points( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

private:
    int resl_;
    float tether_to_input_position_cut_;

};

struct HSearchFilterSortTask : public SearchPointTask {

    HSearchFilterSortTask(
        int resl,
        uint64_t num_to_keep,
        float global_score_cut,
        bool prune_extra ) :
        resl_( resl ),
        num_to_keep_( num_to_keep ),
        global_score_cut_( global_score_cut ),
        prune_extra_( prune_extra )
        {}

    shared_ptr<std::vector<SearchPoint>> 
    return_search_points( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

private:
    int resl_;
    uint64_t num_to_keep_;
    float global_score_cut_;
    bool prune_extra_;

};

struct HSearchScaleToResl : public SearchPointTask {

    HSearchScaleToResl(
        int current_resl,
        int target_resl,
        int DIMPOW2,
        float global_score_cut
         ) :
        current_resl_( current_resl ),
        target_resl_( target_resl ),
        DIMPOW2_( DIMPOW2 ),
        global_score_cut_( global_score_cut )
        {}

    shared_ptr<std::vector<SearchPoint>> 
    return_search_points( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

private:
    int current_resl_;
    int target_resl_;
    int DIMPOW2_;
    float global_score_cut_;

};

struct HSearchFinishTask : public SearchPointTask {

    HSearchFinishTask(
        float global_score_cut
         ) :
        global_score_cut_( global_score_cut )
        {}

    shared_ptr<std::vector<SearchPoint>> 
    return_search_points( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

private:
    float global_score_cut_;

};



}}

#endif
