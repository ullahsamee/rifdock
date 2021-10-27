// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols
#ifndef INCLUDED_riflib_RifFactory_hh
#define INCLUDED_riflib_RifFactory_hh

#include <riflib/types.hh>
#include <riflib/RifBase.hh>
#include <riflib/rif/RifGenerator.hh>
#include <scheme/objective/integration/SceneObjective.hh>
#include <riflib/RotamerGenerator.hh>
#include <scheme/objective/storage/TwoBodyTable.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/ScoreRotamerVsTarget.hh>
#include <riflib/BurialManager.hh>
#include <riflib/UnsatManager.hh>
#include <riflib/CBTooCloseManager.hh>
#include <riflib/HydrophobicManager.hh>
#include <riflib/AtomsCloseTogetherManager.hh>

#ifdef USEGRIDSCORE
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#endif

#include <string>
#include <vector>
#include <boost/any.hpp>
#include <boost/foreach.hpp>

#include <parallel/algorithm>

namespace scheme { namespace search { struct HackPackOpts; }}

namespace devel {
namespace scheme {

// forward declaration. This is terrible
struct UnsatManager;


/////////////////////// RifFactory ////////////////////////////

struct RifFactoryConfig
{
	std::string rif_type;
	RifFactoryConfig()
		: rif_type("")
	{}
};
struct RifSceneObjectiveConfig;

struct RifFactory
{
	RifFactoryConfig config_;

	RifFactory( RifFactoryConfig const & config ) : config_(config) {}

	virtual	RifPtr
	create_rif( float cart_resl=0, float ang_resl=0, float cart_bound=0 ) const = 0;

	virtual RifPtr
	create_rif_from_rif( RifConstPtr refrif, float cart_resl, float ang_resl, float cart_bound ) const = 0;

	virtual	RifPtr
	create_rif_from_file( std::string const & fname, std::string & description ) const = 0;

	virtual	ScenePtr
	create_scene() const = 0;

	virtual	bool
	create_objectives(
		RifSceneObjectiveConfig const & config,
		std::vector<ObjectivePtr> & objectives,
		std::vector<ObjectivePtr> & packing_objectives
	) const = 0;

	virtual	shared_ptr<rif::RifAccumulator>
	create_rif_accumulator( float cart_resl, float ang_resl, float cart_bound, size_t scratchM ) const = 0;

	RifPtr
	create_rif_from_file( std::string const & fname ) const {
		std::string tmp;
		return create_rif_from_file( fname, tmp );
	}

	RifFactoryConfig const &
	config() const { return config_; }

// DELETE THIS!@!!!!!!!!
	virtual void
    set_sasa_params(
        std::vector<ObjectivePtr> & objectives,
        shared_ptr<BurialVoxelArray> sasa_grid,
        float sasa_threshold,
        float sasa_multiplier
    ) const = 0;

	// This should not be in here. Only here because of the typedefs
	virtual std::vector<shared_ptr<UnsatManager>> &
	get_unsatperthread( ObjectivePtr & objective ) const = 0;
};


shared_ptr<RifFactory>
create_rif_factory( RifFactoryConfig const & config );

std::string get_rif_type_from_file( std::string fname );



struct HackPackOpts;
struct RifSceneObjectiveConfig
{
	bool add_native_scaffold_rots_when_packing;
	::scheme::search::HackPackOpts * packopts;
	std::vector<RifPtr> rif_ptrs;
	std::vector< std::vector< VoxelArrayPtr > > const * target_bounding_by_atype;
	RifScoreRotamerVsTarget rot_tgt_scorer;
	shared_ptr< ::devel::scheme::RotamerIndex> rot_index_p;
	int n_sat_groups;
	int require_satisfaction;
	int require_n_rifres;
    int require_hydrophobic_residue_contacts;
    float hydrophobic_ddg_cut;
	float scaff_bb_hbond_weight;
    std::vector< int > requirements;
    std::vector<std::pair<int,std::vector<int>>> requirement_groups;
    shared_ptr<BurialManager> burial_manager;
    shared_ptr<UnsatManager> unsat_manager;
    shared_ptr<BurialVoxelArray> sasa_grid;
    shared_ptr<HydrophobicManager> hydrophobic_manager;
    shared_ptr<CBTooCloseManager> CB_too_close_manager;
    shared_ptr<std::vector<AtomsCloseTogetherManager>> atoms_close_together_managers_p;
    float sasa_threshold;
    float sasa_multiplier;
    float ignore_rifres_if_worse_than;
    int num_pdbinfo_requirements_required;
    std::vector< std::vector<bool> > pdbinfo_req_active_positions;
    std::vector< std::vector<bool> > pdbinfo_req_active_requirements;
    std::vector<float> sat_bonus;
    std::vector<bool> sat_bonus_override;

};





}}

#endif
