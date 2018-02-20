// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/task/TaskProtocol.hh>
#include <riflib/task/util.hh>

#include <riflib/types.hh>


#include <string>
#include <vector>



namespace devel {
namespace scheme {


void
TaskProtocol::run( shared_ptr<std::vector<SearchPoint>> search_points, RifDockData & rdd, ProtocolData & pd ) {

    TaskType last_task_type = SearchPointTaskType;

    shared_ptr<std::vector<SearchPoint>> working_search_points = search_points;
    shared_ptr<std::vector<SearchPointWithRots>> working_search_point_with_rotss = nullptr;
    shared_ptr<std::vector<RifDockResult>> working_rif_dock_results = nullptr;


    size_t current_taskno = 0;


    while ( current_taskno < tasks_.size() ) {

        Task & task = *tasks_[current_taskno];

        TaskType current_task_type = task.get_task_type();
        TaskType reported_task_type = current_task_type;

        switch (last_task_type) {
            case SearchPointTaskType: {
                runtime_assert( working_search_points );
                runtime_assert( ! working_search_point_with_rotss );
                runtime_assert( ! working_rif_dock_results );

                switch (current_task_type) {
                    case SearchPointTaskType: {
                        working_search_points = task.return_search_points(working_search_points, rdd, pd);
                        working_search_point_with_rotss = nullptr;
                        working_rif_dock_results = nullptr;
                        break;
                    }
                    case SearchPointWithRotsTaskType: {
                        working_search_point_with_rotss = task.return_search_point_with_rotss(working_search_points, rdd, pd);
                        working_search_points = nullptr;
                        working_rif_dock_results = nullptr;
                        break;
                    }
                    case RifDockResultTaskType: {
                        working_rif_dock_results = task.return_rif_dock_results(working_search_points, rdd, pd);
                        working_search_points = nullptr;
                        working_search_point_with_rotss = nullptr;
                        break;
                    }
                    case AnyPointTaskType: {
                        working_search_points = task.return_search_points(working_search_points, rdd, pd);
                        working_search_point_with_rotss = nullptr;
                        working_rif_dock_results = nullptr;
                        reported_task_type = SearchPointTaskType;
                        break;
                    }
                    default: { runtime_assert(false); }
                }
                break;
            }
            case SearchPointWithRotsTaskType: {
                runtime_assert( ! working_search_points );
                runtime_assert( working_search_point_with_rotss );
                runtime_assert( ! working_rif_dock_results );

                switch (current_task_type) {
                    case SearchPointTaskType: {
                        working_search_points = task.return_search_points(working_search_point_with_rotss, rdd, pd);
                        working_search_point_with_rotss = nullptr;
                        working_rif_dock_results = nullptr;
                        break;
                    }
                    case SearchPointWithRotsTaskType: {
                        working_search_point_with_rotss = task.return_search_point_with_rotss(working_search_point_with_rotss, rdd, pd);
                        working_search_points = nullptr;
                        working_rif_dock_results = nullptr;
                        break;
                    }
                    case RifDockResultTaskType: {
                        working_rif_dock_results = task.return_rif_dock_results(working_search_point_with_rotss, rdd, pd);
                        working_search_points = nullptr;
                        working_search_point_with_rotss = nullptr;
                        break;
                    }
                    case AnyPointTaskType: {
                        working_search_point_with_rotss = task.return_search_point_with_rotss(working_search_point_with_rotss, rdd, pd);
                        working_search_points = nullptr;
                        working_rif_dock_results = nullptr;
                        reported_task_type = SearchPointWithRotsTaskType;
                        break;
                    }
                    default: { runtime_assert(false); }
                }
                break;
            }
            case RifDockResultTaskType: {
                runtime_assert( ! working_search_points );
                runtime_assert( ! working_search_point_with_rotss );
                runtime_assert( working_rif_dock_results );

                switch (current_task_type) {
                    case SearchPointTaskType: {
                        working_search_points = task.return_search_points(working_rif_dock_results, rdd, pd);
                        working_search_point_with_rotss = nullptr;
                        working_rif_dock_results = nullptr;
                        break;
                    }
                    case SearchPointWithRotsTaskType: {
                        working_search_point_with_rotss = task.return_search_point_with_rotss(working_rif_dock_results, rdd, pd);
                        working_search_points = nullptr;
                        working_rif_dock_results = nullptr;
                        break;
                    }
                    case RifDockResultTaskType: {
                        working_rif_dock_results = task.return_rif_dock_results(working_rif_dock_results, rdd, pd);
                        working_search_points = nullptr;
                        working_search_point_with_rotss = nullptr;
                        break;
                    }
                    case AnyPointTaskType: {
                        working_rif_dock_results = task.return_rif_dock_results(working_rif_dock_results, rdd, pd);
                        working_search_points = nullptr;
                        working_search_point_with_rotss = nullptr;
                        reported_task_type = RifDockResultTaskType;
                        break;
                    }
                    default: { runtime_assert(false); }
                }
                break;
            }
            default: { runtime_assert(false); }
        }


        last_task_type = reported_task_type;
        current_taskno++;

    }

}



}}
