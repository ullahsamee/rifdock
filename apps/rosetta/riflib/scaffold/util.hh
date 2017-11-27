// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_util_hh
#define INCLUDED_riflib_scaffold_util_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <rif_dock_test.hh>

#include <core/pose/Pose.hh>
#include <protocols/indexed_structure_store/movers/DirectSegmentLookupMover.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct ScaffoldDataCache;

void
get_info_for_iscaff(
    uint64_t iscaff,
    RifDockOpt const & opt, 
    std::string & scafftag,
    core::pose::Pose & scaffold,
    utility::vector1<core::Size> & scaffold_res,
    EigenXform & scaffold_perturb
    );

// historically, non_fa was used during HSearch and fa was used during hack pack
ParametricSceneConformationCOP
make_conformation_from_data_cache(shared_ptr<ScaffoldDataCache> cache, bool fa = false) ;


std::vector<core::pose::PoseOP>
apply_direct_segment_lookup_mover( 
    protocols::indexed_structure_store::movers::DirectSegmentLookupMover & dsl_mover,
    core::pose::Pose const & pose );


void
add_pdbinfo_if_missing( core::pose::Pose & pose );











}}

#endif
