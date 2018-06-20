// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_rifdock_tasks_OutputResultsTasks_hh
#define INCLUDED_riflib_rifdock_tasks_OutputResultsTasks_hh

#include <riflib/types.hh>
#include <riflib/task/RifDockResultTask.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct OutputResultsTask : public RifDockResultTask {

    OutputResultsTask( 
        int director_resl,
        int rif_resl
    ) :
        director_resl_( director_resl ),
        rif_resl_( rif_resl )
    {}

    shared_ptr<std::vector<RifDockResult>>
    return_rif_dock_results( 
        shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;



private:
    int director_resl_;
    int rif_resl_;

};

void
dump_rif_result_(
    RifDockData & rdd,
    RifDockResult const & selected_result, 
    std::string const & pdboutfile, 
    int director_resl,
    int rif_resl,
    bool quiet = true,
    std::string const & resfileoutfile = "",
    std::string const & allrifrotsoutfile = "",
    std::vector<float> const & unsat_scores = std::vector<float>(),
    std::vector<int> const & hydrophobic_counts = std::vector<int>()
    );

void
dump_search_point_(
    RifDockData & rdd,
    SearchPoint const & search_point, 
    std::string const & pdboutfile, 
    int director_resl,
    int rif_resl,
    bool quiet);


}}

#endif
