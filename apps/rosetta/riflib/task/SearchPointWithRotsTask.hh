// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_task_SearchPointWithRotsTask_hh
#define INCLUDED_riflib_task_SearchPointWithRotsTask_hh

#include <riflib/types.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/task/Task.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct SearchPointWithRotsTask : public Task {

    SearchPointWithRotsTask() {}



    SearchPoint return_search_points( std::vector<SearchPoint> const & search_points ) override;
    SearchPoint return_search_points( std::vector<SearchPointWithRots> const & search_point_with_rotss ) override;
    SearchPoint return_search_points( std::vector<RifDockResult> const & rif_dock_results ) override;

    SearchPointWithRots return_search_point_with_rotss( std::vector<SearchPoint> const & search_points ) override;
    virtual SearchPointWithRots return_search_point_with_rotss( std::vector<SearchPointWithRots> const & search_point_with_rotss ) = 0;
    SearchPointWithRots return_search_point_with_rotss( std::vector<RifDockResult> const & rif_dock_results ) override;

    RifDockResult return_rif_dock_results( std::vector<SearchPoint> const & search_points ) override;
    RifDockResult return_rif_dock_results( std::vector<SearchPointWithRots> const & search_point_with_rotss ) override;
    RifDockResult return_rif_dock_results( std::vector<RifDockResult> const & rif_dock_results ) override;




};



}}

#endif
