// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_task_TaskProtocol_hh
#define INCLUDED_riflib_task_TaskProtocol_hh

#include <riflib/types.hh>
#include <riflib/task/Task.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct ThreePointVectors {
    shared_ptr<std::vector<SearchPoint>> search_points;
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss;
    shared_ptr<std::vector<RifDockResult>> rif_dock_results;
};

struct TaskProtocol {

    TaskProtocol( std::vector<shared_ptr<Task>> const & tasks ) :
    tasks_( tasks )
    {}


    ThreePointVectors
    run( ThreePointVectors input, RifDockData & rdd, ProtocolData & pd );


private:


    std::vector<shared_ptr<Task>> tasks_;



};



}}

#endif
