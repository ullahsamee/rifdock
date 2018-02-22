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
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/task/SearchPointWithRotsTask.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct RosettaScoreTask : public SearchPointWithRotsTask {

    RosettaScoreTask( float rosetta_score_cut, bool will_do_min ) :
        rosetta_score_cut_( rosetta_score_cut ),
        will_do_min_( will_do_min )
        {}

    shared_ptr<std::vector<SearchPointWithRots>>
    return_search_point_with_rotss( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;



private:
    float rosetta_score_cut_;
    bool will_do_min_;

};

struct RosettaMinTask : public SearchPointWithRotsTask {

    RosettaMinTask( float rosetta_score_cut ) :
    rosetta_score_cut_( rosetta_score_cut )
    {}

    shared_ptr<std::vector<SearchPointWithRots>>
    return_search_point_with_rotss( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;


private:
    float rosetta_score_cut_;

};

struct FilterForRosettaScoreTask : public SearchPointWithRotsTask {

    FilterForRosettaScoreTask(
        float rosetta_score_fraction,
        float rosetta_score_then_min_below_thresh,
        int rosetta_score_at_least,
        int rosetta_score_at_most
        ) :
        rosetta_score_fraction_( rosetta_score_fraction ),
        rosetta_score_then_min_below_thresh_( rosetta_score_then_min_below_thresh ),
        rosetta_score_at_least_( rosetta_score_at_least ),
        rosetta_score_at_most_( rosetta_score_at_most )
        {}

    shared_ptr<std::vector<SearchPointWithRots>>
    return_search_point_with_rotss( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

private:
    float rosetta_score_fraction_;
    float rosetta_score_then_min_below_thresh_;
    int rosetta_score_at_least_;
    int rosetta_score_at_most_;

};

struct FilterForRosettaMinTask : public SearchPointWithRotsTask {

    FilterForRosettaMinTask( float rosetta_min_fraction ) :
    rosetta_min_fraction_( rosetta_min_fraction )
    {}

    shared_ptr<std::vector<SearchPointWithRots>>
    return_search_point_with_rotss( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

private:
    float rosetta_min_fraction_;

};

// struct RedundancyFilterForRosettaScoreTask : public SearchPointWithRotsTask {

//     RedundancyFilterForRosettaScoreTask( int n_per_block, float redundancy_mag ) :
//         n_per_block_( n_per_block ),
//         redundancy_mag_( redundancy_mag )
//     {}

//     shared_ptr<std::vector<SearchPointWithRots>>
//     return_search_point_with_rotss( 
//         shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
//         RifDockData & rdd, 
//         ProtocolData & pd );

// private:
//     int n_per_block_;
//     float redundancy_mag_;

// };

void
rosetta_score_inner( 
    shared_ptr<std::vector< SearchPointWithRots >>  packed_results,
    RifDockData & rdd,
    ProtocolData & pd,
    float rosetta_score_cut,
    bool is_minimizing,
    bool will_do_min
    );


}}

#endif
