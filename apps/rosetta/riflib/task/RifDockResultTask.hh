// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_task_RifDockResultTask_hh
#define INCLUDED_riflib_task_RifDockResultTask_hh

#include <riflib/types.hh>
#include <riflib/task/Task.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct RifDockResultTask : public Task {

    RifDockResultTask() {}



    shared_ptr<std::vector<SearchPoint>> return_search_points( shared_ptr<std::vector<SearchPoint>> search_points, RifDockData & rdd, ProtocolData & pd ) override;
    shared_ptr<std::vector<SearchPoint>> return_search_points( shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, RifDockData & rdd, ProtocolData & pd ) override;
    shared_ptr<std::vector<SearchPoint>> return_search_points( shared_ptr<std::vector<RifDockResult>> rif_dock_results, RifDockData & rdd, ProtocolData & pd ) override;

    shared_ptr<std::vector<SearchPointWithRots>> return_search_point_with_rotss( shared_ptr<std::vector<SearchPoint>> search_points, RifDockData & rdd, ProtocolData & pd ) override;
    shared_ptr<std::vector<SearchPointWithRots>> return_search_point_with_rotss( shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, RifDockData & rdd, ProtocolData & pd ) override;
    shared_ptr<std::vector<SearchPointWithRots>> return_search_point_with_rotss( shared_ptr<std::vector<RifDockResult>> rif_dock_results, RifDockData & rdd, ProtocolData & pd ) override;

    shared_ptr<std::vector<RifDockResult>> return_rif_dock_results( shared_ptr<std::vector<SearchPoint>> search_points, RifDockData & rdd, ProtocolData & pd ) override;
    shared_ptr<std::vector<RifDockResult>> return_rif_dock_results( shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, RifDockData & rdd, ProtocolData & pd ) override;
    virtual shared_ptr<std::vector<RifDockResult>> return_rif_dock_results( shared_ptr<std::vector<RifDockResult>> rif_dock_results, RifDockData & rdd, ProtocolData & pd ) = 0;

    TaskType get_task_type() const { return RifDockResultTaskType; }



};



}}

#endif
