// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_MorphingScaffoldProvider_hh
#define INCLUDED_riflib_scaffold_MorphingScaffoldProvider_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/scaffold/ScaffoldProviderBase.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>

#include <scheme/kinematics/Scene.hh>

#include <core/pose/Pose.hh>

// Key Assumptions of this class:
// - Only one segment is morphable


namespace devel {
namespace scheme {




struct MorphInfo {
    Range morph_range;
};

typedef shared_ptr<MorphInfo> MorphInfoOP;
typedef shared_ptr<MorphInfo const> MorphInfoCOP;


struct MorphMember {
    core::pose::PoseCOP pose;
    ParametricSceneConformation conformation;
    MorphInfoCOP morph_info;
};




struct MorphingScaffoldProvider :
    public TreeScaffoldProvider<ParametricSceneConformation, ScaffoldDataCache> {

    MorphingScaffoldProvider(); 


    ParametricSceneConformation get_scaffold(TreeIndex i) override;
    ScaffoldDataCache get_scaffold_cache() override;

    TreeLimits get_index_limits() override;


    // By definition, all members at depth 0 are provided by
    // the user
    std::vector< std::vector< MorphMember > > map_;

};



}}

#endif
