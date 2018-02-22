// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_rifdock_tasks_CompileAndFilterResultsTask_hh
#define INCLUDED_riflib_rifdock_tasks_CompileAndFilterResultsTask_hh

#include <riflib/types.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/task/RifDockResultTask.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct CompileAndFilterResultsTask : public RifDockResultTask {

    CompileAndFilterResultsTask( 
        int resl,
        int n_per_block,
        float redundancy_mag,
        int force_output_if_close_to_input_num,
        float force_output_if_close_to_input,
        bool filter_seeding_positions_separately,
        bool filter_scaffolds_separately
         ) : 
        resl_( resl ),
        n_per_block_( n_per_block ),
        redundancy_mag_( redundancy_mag ),
        force_output_if_close_to_input_num_( force_output_if_close_to_input_num ),
        force_output_if_close_to_input_( force_output_if_close_to_input ),
        filter_seeding_positions_separately_( filter_seeding_positions_separately ),
        filter_scaffolds_separately_( filter_scaffolds_separately )
    {}

    shared_ptr<std::vector<RifDockResult>>
    return_rif_dock_results( 
        shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

    // nonstandard overrides because we actually take SeachPointWithRots

    shared_ptr<std::vector<RifDockResult>>
    return_rif_dock_results( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

    shared_ptr<std::vector<RifDockResult>> 
    return_rif_dock_results( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;


private:
    int resl_;
    int n_per_block_;
    float redundancy_mag_;
    int force_output_if_close_to_input_num_;
    float force_output_if_close_to_input_;
    bool filter_seeding_positions_separately_;
    bool filter_scaffolds_separately_;


};



}}

#endif
