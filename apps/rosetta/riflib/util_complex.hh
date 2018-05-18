// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_util_complex_hh
#define INCLUDED_riflib_util_complex_hh

#include <riflib/types.hh>


#ifdef USEGRIDSCORE
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <protocols/ligand_docking/GALigandDock/RotamerData.hh>
#endif
#include <riflib/DonorAcceptorCache.hh>

namespace devel {
namespace scheme {


void
all_by_all_rmsd( 
	std::vector<core::pose::PoseOP> const & poses,
	utility::vector1<utility::vector1<core::Real>> & table );

std::vector<std::vector<std::pair<core::pose::PoseOP, uint64_t>>>
cluster_poses_into_n_bins( 
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n,
	utility::vector1<utility::vector1<core::Real>> & rmsds );

std::vector<core::pose::PoseOP>
cluster_poses_leaving_n( 
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n );

size_t
find_cluster_center( 
	std::vector<uint64_t> indexes,
	utility::vector1<utility::vector1<core::Real>> const & rmsds );

std::vector<core::pose::PoseOP>
cluster_poses_leaving_n_representing_frac(
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n,
	float frac,
	float tol
	);

std::vector<core::pose::PoseOP>
random_selection_poses_leaving_n( 
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n );


#ifdef USEGRIDSCORE
shared_ptr<protocols::ligand_docking::ga_ligand_dock::GridScorer>
prepare_grid_scorer(
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	std::string const & atype_aas = "ACDEFGHIKLMNPQRSTVWY"
);
#endif

void
prepare_donor_acceptor_cache( 
    std::vector<HBondRay> const & target_donors,
    std::vector<HBondRay> const & target_acceptors,
    RifScoreRotamerVsTarget const & rot_tgt_scorer,
    shared_ptr<DonorAcceptorCache> & target_donor_cache,
    shared_ptr<DonorAcceptorCache> & target_acceptor_cache
);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
}

#endif
