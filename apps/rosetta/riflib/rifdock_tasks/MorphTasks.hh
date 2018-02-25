// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_rifdock_tasks_MorphTasks_hh
#define INCLUDED_riflib_rifdock_tasks_MorphTasks_hh

#include <riflib/types.hh>
#include <riflib/task/SearchPointTask.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct FilterByRmsdToThisPdbTask : public SearchPointTask {

    FilterByRmsdToThisPdbTask(
        std::string pose_filename,
        float rmsd
        ) :
        pose_filename_( pose_filename ),
        rmsd_( rmsd )
        {}

    shared_ptr<std::vector<SearchPoint>> 
    return_search_points( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;



private:
    std::string pose_filename_;
    float rmsd_;

};

struct TestMakeChildrenTask : public SearchPointTask {

    TestMakeChildrenTask(
        ) 
        {}

    shared_ptr<std::vector<SearchPoint>> 
    return_search_points( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;



private:


};


void
create_dive_pop_hsearch_task( 
    std::vector<shared_ptr<Task>> & task_list, RifDockData & rdd );



}}

#endif
