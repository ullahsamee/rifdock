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
#include <riflib/scaffold/ScaffoldProviderBase.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>

// Key Assumptions of this class:
// - Only one segment is morphable


namespace devel {
namespace scheme {




struct MorphInfo {
    uint64_t low_mod_point;
    uint64_t high_mod_point;
};

typedef shared_ptr<MorphInfo> MorphInfoOP;
typedef shared_ptr<MorphInfo const> MorphInfoCOP;


struct MorphMember {
    core::pose::PoseCOP pose;
    MorphInfoCOP morph_info;
};



struct MorphingScaffoldProvider :
    public ScaffoldProviderBase<MorphIndex> {

    MorphingScaffoldProvider(); 


    void get_scaffold(MorphIndex i) override;



    // By definition, all members at depth 0 are provided by
    // the user
    std::vector< std::vector< MorphMember > > map_;

};



}}

#endif
