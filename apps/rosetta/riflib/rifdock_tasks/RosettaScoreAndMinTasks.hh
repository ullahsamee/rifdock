// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_rifdock_tasks_RosettaScoreAndMinTasks_hh
#define INCLUDED_riflib_rifdock_tasks_RosettaScoreAndMinTasks_hh

#include <riflib/types.hh>
#include <riflib/task/SearchPointWithRotsTask.hh>
#include <riflib/task/AnyPointTask.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {


struct FilterForRosettaScoreTask : public AnyPointTask {

    FilterForRosettaScoreTask(
        float rosetta_score_fraction,
        float rosetta_score_then_min_below_thresh,
        int rosetta_score_at_least,
        int rosetta_score_at_most,
        bool rosetta_select_random
        ) :
        rosetta_score_fraction_( rosetta_score_fraction ),
        rosetta_score_then_min_below_thresh_( rosetta_score_then_min_below_thresh ),
        rosetta_score_at_least_( rosetta_score_at_least ),
        rosetta_score_at_most_( rosetta_score_at_most ),
        rosetta_select_random_( rosetta_select_random )
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
        ProtocolData & pd );

private:
    float rosetta_score_fraction_;
    float rosetta_score_then_min_below_thresh_;
    int rosetta_score_at_least_;
    int rosetta_score_at_most_;
    bool rosetta_select_random_;

};

struct FilterForRosettaMinTask : public AnyPointTask {

    FilterForRosettaMinTask( float rosetta_min_fraction, int rosetta_min_at_least, int rosetta_min_at_most ) :
    rosetta_min_fraction_( rosetta_min_fraction ),
    rosetta_min_at_least_( rosetta_min_at_least ),
    rosetta_min_at_most_( rosetta_min_at_most )
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
        ProtocolData & pd );

private:
    float rosetta_min_fraction_;
    int rosetta_min_at_least_;
    int rosetta_min_at_most_;

};



struct RosettaScoreTask : public SearchPointWithRotsTask {

    RosettaScoreTask( int director_resl, float rosetta_score_cut, bool will_do_min, bool store_pose ) :
        director_resl_( director_resl ),
        rosetta_score_cut_( rosetta_score_cut ),
        will_do_min_( will_do_min ),
        store_pose_( store_pose )
        {}

    shared_ptr<std::vector<SearchPointWithRots>>
    return_search_point_with_rotss( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;



private:
    int director_resl_;
    float rosetta_score_cut_;
    bool will_do_min_;
    bool store_pose_;

};

struct RosettaMinTask : public SearchPointWithRotsTask {

    RosettaMinTask( int director_resl, float rosetta_score_cut, bool store_pose ) :
    director_resl_( director_resl ),
    rosetta_score_cut_( rosetta_score_cut ),
    store_pose_( store_pose )
    {}

    shared_ptr<std::vector<SearchPointWithRots>>
    return_search_point_with_rotss( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;


private:
    int director_resl_;
    float rosetta_score_cut_;
    bool store_pose_;

};


void
rosetta_score_inner( 
    shared_ptr<std::vector< SearchPointWithRots >>  packed_results,
    RifDockData & rdd,
    ProtocolData & pd,
    int director_resl,
    float rosetta_score_cut,
    bool is_minimizing,
    bool will_do_min,
    bool store_pose
    );


}}

#endif
