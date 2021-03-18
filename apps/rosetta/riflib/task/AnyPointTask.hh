// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_task_AnyPointTask_hh
#define INCLUDED_riflib_task_AnyPointTask_hh

#include <riflib/types.hh>
#include <riflib/task/Task.hh>

#include <string>
#include <vector>

// This class doesn't play nicely with the C++ standard
// It's only supposed to have one virtual function
// which is return_any_points
// Instead the virtualization must be done manually
// See the .cc file for implementations of the other functions


namespace devel {
namespace scheme {

struct AnyPointTask : public Task {

    AnyPointTask() {}

    // template<class AnyPoint>
    // virtual
    // shared_ptr<std::vector<AnyPoint>>
    // return_any_points( 
    //     shared_ptr<std::vector<AnyPoint>> any_points, 
    //     RifDockData & rdd, 
    //     ProtocolData & pd ) = 0;


// These are not supposed to be virtual!!! Must do them manually

    virtual shared_ptr<std::vector<SearchPoint>> return_search_points( shared_ptr<std::vector<SearchPoint>> search_points, RifDockData & rdd, ProtocolData & pd ) = 0;  
    shared_ptr<std::vector<SearchPoint>> return_search_points( shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, RifDockData & rdd, ProtocolData & pd ) override;
    shared_ptr<std::vector<SearchPoint>> return_search_points( shared_ptr<std::vector<RifDockResult>> rif_dock_results, RifDockData & rdd, ProtocolData & pd ) override;

    shared_ptr<std::vector<SearchPointWithRots>> return_search_point_with_rotss( shared_ptr<std::vector<SearchPoint>> search_points, RifDockData & rdd, ProtocolData & pd ) override;
    virtual shared_ptr<std::vector<SearchPointWithRots>> return_search_point_with_rotss( shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, RifDockData & rdd, ProtocolData & pd ) = 0;
    shared_ptr<std::vector<SearchPointWithRots>> return_search_point_with_rotss( shared_ptr<std::vector<RifDockResult>> rif_dock_results, RifDockData & rdd, ProtocolData & pd ) override;

    shared_ptr<std::vector<RifDockResult>> return_rif_dock_results( shared_ptr<std::vector<SearchPoint>> search_points, RifDockData & rdd, ProtocolData & pd ) override;
    shared_ptr<std::vector<RifDockResult>> return_rif_dock_results( shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, RifDockData & rdd, ProtocolData & pd ) override;
    virtual shared_ptr<std::vector<RifDockResult>> return_rif_dock_results( shared_ptr<std::vector<RifDockResult>> rif_dock_results, RifDockData & rdd, ProtocolData & pd ) = 0;


// boilerplate code for manual virtualization

    // shared_ptr<std::vector<SearchPoint>> 
    // return_search_points( 
    //     shared_ptr<std::vector<SearchPoint>> search_points, 
    //     RifDockData & rdd, 
    //     ProtocolData & pd ) override;

    // shared_ptr<std::vector<SearchPointWithRots>> 
    // return_search_point_with_rotss( 
    //     shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    //     RifDockData & rdd, 
    //     ProtocolData & pd ) override;

    // shared_ptr<std::vector<RifDockResult>> 
    // return_rif_dock_results( 
    //     shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    //     RifDockData & rdd, 
    //     ProtocolData & pd ) override;



    TaskType get_task_type() const { return AnyPointTaskType; }


};



}}

#endif
