// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/scaffold/util.hh>
#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <scheme/scaffold/ScaffoldProviderBase.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>
#include <scheme/numeric/rand_xform.hh>
#include <core/pose/PDBInfo.hh>

// includes to emulate minimize_segment.xml
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <protocols/residue_selectors/StoredResidueSubsetSelector.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>

#include <scheme/kinematics/Scene.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/file/file_sys_util.hh>



namespace devel {
namespace scheme {



void
get_info_for_iscaff(
    uint64_t iscaff,
    RifDockOpt const & opt, 
    std::string & scafftag,
    core::pose::Pose & scaffold,
    utility::vector1<core::Size> & scaffold_res,
    EigenXform & scaffold_perturb
    ) {

    std::string scaff_fname = opt.scaffold_fnames.at(iscaff);
    scafftag = utility::file_basename( utility::file::file_basename( scaff_fname ) );

    std::cout << "!!!!!!!!!!!!!!!name:: " << scaff_fname << std::endl;

    core::import_pose::pose_from_file( scaffold, scaff_fname );

    scaffold_perturb = EigenXform::Identity();
    if( opt.random_perturb_scaffold ){
        runtime_assert_msg( !opt.use_scaffold_bounding_grids,
            "opt.use_scaffold_bounding_grids incompatible with random_perturb_scaffold" );
        std::mt19937 rng( 0);// std::random_device{}() );
        ::scheme::numeric::rand_xform(rng,scaffold_perturb);
        xform_pose( scaffold, eigen2xyz(scaffold_perturb) );
    }


    scaffold_res.clear();
    std::string scaff_res_fname = "";
    if( opt.scaffold_res_fnames.size() ){
        if( opt.scaffold_res_fnames.size() == opt.scaffold_fnames.size() ){
            scaff_res_fname = opt.scaffold_res_fnames.at(iscaff);
        } else if( opt.scaffold_res_fnames.size() == 1 ){
            scaff_res_fname = opt.scaffold_res_fnames.front();
        } else {
            utility_exit_with_message( "-scaffold_res list not same length as -scaffolds list" );
        }
        if( opt.scaffold_res_use_best_guess ){
            utility_exit_with_message("should only use -scaffold_res_use_best_guess true iff not specifying scaffold_res");
        }
        scaffold_res = devel::scheme::get_res( scaff_res_fname , scaffold );
    } else if (opt.scaffold_res_use_best_guess ){
        scaffold_res = devel::scheme::get_designable_positions_best_guess( scaffold, opt.dont_use_scaffold_loops );
        std::cout << "using scaffold residues: ";
        for(auto ir:scaffold_res) std::cout << " " << ir << scaffold.residue(ir).name3();
        std::cout << std::endl;
    } else {
        for( int ir = 1; ir <= scaffold.size(); ++ir){
            if( !scaffold.residue(ir).is_protein() ) continue;
            //if( scaffold.residue(ir).name3() == "PRO" ) continue;
            //if( scaffold.residue(ir).name3() == "GLY" ) continue;
            //if( scaffold.residue(ir).name3() == "CYS" ) continue;
            scaffold_res.push_back(ir);
        }
    }

}

// historically, non_fa was used during HSearch and fa was used during hack pack
ParametricSceneConformationCOP
make_conformation_from_data_cache(ScaffoldDataCacheOP cache, bool fa /*= false*/) {
    typedef numeric::xyzVector<core::Real> Vec;
    ParametricScene scene(1);

    core::pose::Pose const & scaffold_centered = *(cache->scaffold_centered_p);
    std::vector<int> const & scaffres_g2l = *(cache->scaffres_g2l_p);
    utility::vector1<core::Size> const & scaffold_res = *(cache->scaffold_res_p);
    std::vector< SimpleAtom > const & scaffold_simple_atoms = *(cache->scaffold_simple_atoms_p);
    std::vector< SimpleAtom > const & scaffold_simple_atoms_all = *(cache->scaffold_simple_atoms_all_p);


    for( int ir = 1; ir <= scaffold_centered.size(); ++ir ){
        Vec N  = scaffold_centered.residue(ir).xyz("N" );
        Vec CA = scaffold_centered.residue(ir).xyz("CA");
        Vec C  = scaffold_centered.residue(ir).xyz("C" );

        // todo map res indices, must also edit onebody_energies
        BBActor bbactor( N, CA, C, '-', '-', scaffres_g2l[ir-1] );
        runtime_assert( bbactor.index_ == scaffres_g2l[ir-1] );


        if( std::find(scaffold_res.begin(),scaffold_res.end(),ir)!=scaffold_res.end() ){
            scene.add_actor(0,bbactor);
        }
    }

    if (fa) {
        BOOST_FOREACH( SimpleAtom const & sa, scaffold_simple_atoms_all ) scene.add_actor( 0, sa );
        runtime_assert( scene.template num_actors<SimpleAtom>(0) == scaffold_simple_atoms_all.size() );
    } else {
        BOOST_FOREACH( SimpleAtom const & sa, scaffold_simple_atoms ) scene.add_actor( 0, sa );
        runtime_assert( scene.template num_actors<SimpleAtom>(0) == scaffold_simple_atoms.size() );
    }

    ParametricSceneConformationCOP conformation = scene.conformation_ptr(0);

    ParametricSceneConformationOP conformation_mutable = std::const_pointer_cast<ParametricSceneConformation>( conformation );

    conformation_mutable->cache_data_ = cache;
    conformation_mutable->cache_data_->conformation_is_fa = fa;

    std::cout << "FA status is: " << (fa ? "True" : "False") << std::endl;

    return conformation;


}

std::vector<core::pose::PoseOP>
apply_direct_segment_lookup_mover( 
    protocols::indexed_structure_store::movers::DirectSegmentLookupMover & dsl_mover,
    core::pose::Pose const & pose ) {

    using namespace core::pack::task::operation;
    using namespace core::select::residue_selector;
    using namespace protocols::residue_selectors;

    const std::string stored_subset_name = "inserted_lookup_segment";
    dsl_mover.stored_subset_name( stored_subset_name );
    dsl_mover.structure_store_path( basic::options::option[basic::options::OptionKeys::indexed_structure_store::fragment_store]() );

    //<SCOREFXNS>

    core::scoring::ScoreFunctionOP scorefxn = make_shared<core::scoring::ScoreFunction>();
    scorefxn->set_weight(core::scoring::coordinate_constraint, 2);
    scorefxn->set_weight(core::scoring::cart_bonded, 2);
    scorefxn->set_weight(core::scoring::pro_close, 0);

    //<TASKOPERATIONS>

    RestrictAbsentCanonicalAASOP ala_only( new RestrictAbsentCanonicalAAS() );
    ala_only->include_residue( 0 );
    ala_only->keep_aas( "A" );

    ResidueSelectorCOP stored_residue_subset( new StoredResidueSubsetSelector( stored_subset_name ) );

    ResLvlTaskOperationOP prevent_repacking_rlt( new PreventRepackingRLT() );
    OperateOnResidueSubsetOP only_lookup_segment( new OperateOnResidueSubset( prevent_repacking_rlt, stored_residue_subset, true ) );


    // <MOVERS>

    core::pack::task::TaskFactoryOP to_ala_tf( new core::pack::task::TaskFactory() );
    to_ala_tf->push_back( only_lookup_segment );
    to_ala_tf->push_back( ala_only );

    protocols::simple_moves::PackRotamersMover to_ala;
    to_ala.score_function(scorefxn);
    to_ala.task_factory( to_ala_tf );


    core::pack::task::TaskFactoryOP hardmin_bb_tf( new core::pack::task::TaskFactory() );
    hardmin_bb_tf->push_back( only_lookup_segment );

    protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover() );
    min_mover->tolerance( 0.0001 );
    min_mover->min_type( "lbfgs_armijo_nonmonotone" );
    min_mover->cartesian( true );
    min_mover->score_function( scorefxn );

    protocols::simple_moves::TaskAwareMinMover hardmin_bb( min_mover, hardmin_bb_tf );
    hardmin_bb.chi( true );
    hardmin_bb.bb( true );


    //<PROTOCOLS>

    core::pose::PoseOP result = pose.clone();
    dsl_mover.apply( *result );

    std::vector<core::pose::PoseOP> results;
    do {
        if ( result->num_chains() > 1 ) {
            std::cout << "Broken pose" << std::endl;
            continue;
        }
        to_ala.apply( *result );
        hardmin_bb.apply( *result );

        results.push_back( result );
    } while ( ( result = dsl_mover.get_additional_output() ) );



    return results;
}

void
add_pdbinfo_if_missing( core::pose::Pose & pose ) {
    if ( ! pose.pdb_info() ) {
        core::pose::PDBInfoOP pdb_info = make_shared<core::pose::PDBInfo>( pose );
        pose.pdb_info(pdb_info);
    }
}



}}

