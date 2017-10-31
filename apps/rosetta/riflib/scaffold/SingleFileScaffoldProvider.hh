// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_SingleFileScaffoldProvider_hh
#define INCLUDED_riflib_scaffold_SingleFileScaffoldProvider_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <scheme/scaffold/ScaffoldProviderBase.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>

#include <scheme/kinematics/Scene.hh>

#include <core/pose/Pose.hh>

// Key Assumptions of this class:
// - Only one segment is morphable


namespace devel {
namespace scheme {




struct SingleFileScaffoldProvider :
    public ::scheme::scaffold::ScaffoldProviderBase<ParametricSceneConformation, uint64_t, uint64_t> {

    SingleFileScaffoldProvider();


    ParametricSceneConformationCOP get_scaffold(uint64_t i) override;

    uint64_t get_scaffold_index_limits() override;

    ScaffoldDataCacheOP temp_function__get_writable_data_cache() {
        temp_data__data_cache_ = make_shared<ScaffoldDataCache>();
        return temp_data__data_cache_;
    }


    ParametricSceneConformationOP conformation_;
    core::pose::PoseCOP pose_;

    ScaffoldDataCacheOP temp_data__data_cache_;

};



}}

#endif
