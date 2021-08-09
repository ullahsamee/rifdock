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
#include <riflib/RotamerGenerator.hh>
#include <scheme/numeric/rand_xform.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>

// includes to emulate minimize_segment.xml
#include <core/io/silent/SilentFileData.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <protocols/residue_selectors/StoredResidueSubsetSelector.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <basic/options/option.hh>
#include <protocols/moves/MoverStatus.hh>
#include <core/select/util.hh>
#include <core/sequence/SequenceProfile.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#ifdef USEHDF5
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>
#endif

#include <string>
#include <vector>
#include <boost/any.hpp>
#include <boost/format.hpp>

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
    EigenXform & scaffold_perturb,
    MorphRules & morph_rules,
    ExtraScaffoldData & extra_data,
    shared_ptr< RotamerIndex > rot_index_p,
    std::string & scaffold_filename
    ) {

    std::string scaff_fname = opt.scaffold_fnames.at(iscaff);
    scaffold_filename = scaff_fname; // NRB
    scafftag = pdb_name(scaff_fname);

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
    } else if ( !opt.scaffold_res_pdbinfo_labels.empty() ){

        utility::vector1<bool> scaffold_res_mask(scaffold.size(), false);

        std::cout << "Selecting scaffold residues using res-labels:";

        for ( std::string const & label : opt.scaffold_res_pdbinfo_labels ) {
            std::cout << " " << label;
            for ( core::Size i = 1; i <= scaffold.size(); i++ ) {
                if ( scaffold.pdb_info()->res_haslabel(i, label) ) {
                    scaffold_res_mask[i] = true;
                }
            }
        }

        scaffold_res = core::select::get_residues_from_subset( scaffold_res_mask, true );

        std::cout << std::endl;
        std::cout << "using scaffold residues: ";
        for(auto ir:scaffold_res) std::cout << " " << ir << scaffold.residue(ir).name3();
        std::cout << std::endl;

    } else if (opt.scaffold_res_use_best_guess ){
        scaffold_res = devel::scheme::get_designable_positions_best_guess( scaffold, opt.dont_use_scaffold_loops, true, 
                                                                            opt.best_guess_mutate_to_val, opt.dont_use_scaffold_helices,
                                                                            opt.dont_use_scaffold_strands );
        std::cout << "using scaffold residues: ";
        for(auto ir:scaffold_res) std::cout << " " << ir << scaffold.residue(ir).name3();
        std::cout << std::endl;
    } else {
        get_default_scaffold_res( scaffold, scaffold_res );
    }


    if( opt.cst_fnames.size() ){
        std::string cst_fname = "";
        if( opt.cst_fnames.size() == opt.scaffold_fnames.size() ){
            cst_fname = opt.cst_fnames.at(iscaff);
        } else if( opt.cst_fnames.size() == 1 ){
            cst_fname = opt.cst_fnames.front();
        } else {
            utility_exit_with_message( "-cst_files list not same length as -scaffolds list" );
        }
        runtime_assert_msg(parse_constrains_file(cst_fname, extra_data.csts), "Faild to parse the constrain file!!!");
    }

    if( opt.morph_rules_fnames.size() ){
        std::string morph_rules_fname = "";
        if( opt.morph_rules_fnames.size() == opt.scaffold_fnames.size() ){
            morph_rules_fname = opt.morph_rules_fnames.at(iscaff);
        } else if( opt.morph_rules_fnames.size() == 1 ){
            morph_rules_fname = opt.morph_rules_fnames.front();
        } else {
            utility_exit_with_message( "-morph_rules_files list not same length as -scaffolds list" );
        }
        runtime_assert_msg(parse_morph_rules_file(morph_rules_fname, morph_rules, opt), "Faild to parse the morph_rules file!!!");
    } else {
        // std::cout << "Morph rules file not specified, using command-line options" << std::endl;
        morph_rules.push_back(morph_rule_from_options(opt));
    }

    if( opt.rotamer_boltzmann_fnames.size() ){
        std::string rotboltz_fname = "";
        if( opt.rotamer_boltzmann_fnames.size() == opt.scaffold_fnames.size() ){
            rotboltz_fname = opt.rotamer_boltzmann_fnames.at(iscaff);
        } else if( opt.rotamer_boltzmann_fnames.size() == 1 ){
            rotboltz_fname = opt.rotamer_boltzmann_fnames.front();
        } else {
            utility_exit_with_message( "-rotamer_boltzmann_files list not same length as -scaffolds list" );
        }
        extra_data.accumulate_per_rotamer_custom_energies( load_rotboltz_data( rotboltz_fname, scaffold.size(), rot_index_p->size(), opt.rotboltz_ignore_missing_rots ) );
    }

    if( opt.pssm_file_fnames.size() ){
        std::string pssm_fname = "";
        if( opt.pssm_file_fnames.size() == opt.scaffold_fnames.size() ){
            pssm_fname = opt.pssm_file_fnames.at(iscaff);
        } else if( opt.pssm_file_fnames.size() == 1 ){
            pssm_fname = opt.pssm_file_fnames.front();
        } else {
            utility_exit_with_message( "-pssm_file list not same length as -scaffolds list" );
        }
        load_pssm_data( extra_data, pssm_fname, scaffold_res, rot_index_p, opt.pssm_weight, opt.pssm_cutoff, opt.pssm_higher_is_better,
                                                                                                            opt.pssm_enforce_no_ala );
    }

    if( opt.scaffold_clash_contexts.size() ){
        std::string clash_fname = "";
        if( opt.scaffold_clash_contexts.size() == opt.scaffold_fnames.size() ){
            clash_fname = opt.scaffold_clash_contexts.at(iscaff);
        } else if( opt.scaffold_clash_contexts.size() == 1 ){
            clash_fname = opt.scaffold_clash_contexts.front();
        } else {
            utility_exit_with_message( "-scaffold_clash_contexts list not same length as -scaffolds list" );
        }
        load_clash_context( extra_data, clash_fname, rot_index_p );
    }

}

core::Size
get_bb_donor_atom( core::conformation::Residue const & res ) {
    core::chemical::ResidueType const & type = res.type();

    core::chemical::AtomIndices const & h_indices = type.Hpos_polar();

    for ( core::Size hindex : h_indices ) {
        core::Size bindex = type.atom_base(hindex);
        if ( res.atom_is_backbone( bindex ) ) {
            return hindex;
        }
    }
    utility_exit_with_message( boost::str( boost::format( "Residue %i has not backbone hydrogen bond donor!!!" ) % res.seqpos() ) );
}

core::Size
get_bb_acceptor_heavy_atom( core::conformation::Residue const & res ) {
    core::chemical::ResidueType const & type = res.type();

    core::chemical::AtomIndices const & acc_indices = type.accpt_pos();

    for ( core::Size aindex : acc_indices ) {
        if ( res.atom_is_backbone( aindex ) ) {
            return aindex;
        }
    }
    utility_exit_with_message( boost::str( boost::format( "Residue %i has not backbone hydrogen bond acceptor!!!" ) % res.seqpos() ) );
}

std::vector<BBHBondActor>
get_bbhbond_actors( core::pose::Pose const & pose ) {
    std::cout << "Generating scaffold backbone hydrogen bonding actors" << std::endl;

    core::pose::Pose pose_to_score = pose;

    float ethresh = -0.01f;
    ::devel::scheme::HBRayOpts hbopt;

    std::vector<BBHBondActor> actors;

    core::scoring::ScoreFunctionOP sf = make_shared<core::scoring::ScoreFunction>();
    core::scoring::methods::EnergyMethodOptions myopt = sf->energy_method_options();
    myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
    myopt.hbond_options().bb_donor_acceptor_check(true);
    sf->set_energy_method_options(myopt);
    sf->set_weight( core::scoring::hbond_sr_bb, 1.0);
    sf->set_weight( core::scoring::hbond_lr_bb, 1.0);
    sf->score(pose_to_score);
    core::scoring::hbonds::HBondSet hbset;
    core::scoring::hbonds::fill_hbond_set(pose_to_score,false,hbset);

    hbset.resize_bb_donor_acceptor_arrays( pose.size() );   // this is so stupid. Why do they do this?

    std::set<int> n_termini;
    std::set<int> c_termini;

    for ( int ichain = 1; ichain <= pose.num_chains(); ichain++ ) {
        n_termini.insert( pose.conformation().chain_begin(ichain) );
        c_termini.insert( pose.conformation().chain_end(ichain) );
    }

    for ( int ires = 1; ires <= pose.size(); ires++ ) {
        core::conformation::Residue const & res = pose.residue(ires);


        bool allow_nh = true;
        bool allow_co = true;

        if ( n_termini.count(ires) > 0 ) allow_nh = false;
        if ( c_termini.count(ires) > 0 ) allow_co = false;

        if ( hbset.don_bbg_in_bb_bb_hbond( ires ) ) allow_nh = false;
        if ( hbset.acc_bbg_in_bb_bb_hbond( ires ) ) allow_co = false;

        if ( res.name3() == "PRO" ) allow_nh = false;

        if ( allow_nh ) {
            std::vector<HBondRay> donors;
            std::vector<std::pair<int,std::string>> anames;
            get_donor_rays( pose, ires, hbopt, donors, anames );

            std::string bb_donor_name = res.atom_name(get_bb_donor_atom( res ) );
            std::vector<HBondRay> nh_donors;
            for ( int i = 0; i < donors.size(); i++ ) {
                if ( anames[i].second == bb_donor_name ) nh_donors.push_back( donors[i] );
            }
            if ( nh_donors.size() != 1 ) {
                utility_exit_with_message(boost::str(boost::format("Found %i backbone donor atoms for residue %i")%nh_donors.size()%ires));
            }

            actors.emplace_back( nh_donors, true, ires );
        }

        if ( allow_co ) {
            std::vector<HBondRay> acceptors;
            std::vector<std::pair<int,std::string>> anames;
            get_acceptor_rays_lkball( pose, ires, hbopt, acceptors, anames );

            std::string bb_heavy_acceptor_name = res.atom_name(get_bb_acceptor_heavy_atom( res ));
            std::vector<HBondRay> co_acceptors;
            for ( int i = 0; i < acceptors.size(); i++ ) {
                if ( anames[i].second == bb_heavy_acceptor_name ) co_acceptors.push_back( acceptors[i] );
            }
            if ( co_acceptors.size() != 2 ) {
                utility_exit_with_message(boost::str(boost::format("Found %i backbone acceptor atoms for residue %i")%co_acceptors.size()%ires));
            }

            actors.emplace_back( co_acceptors, false, ires );
        }

    }

    std::cout << "Scaffold backbone donors: ";
    bool first = true;
    for ( BBHBondActor const & actor : actors ) {
        if ( actor.is_donor() ) {
            if ( !first ) std::cout << ",";
            std::cout << actor.index();
            first = false;
        }
    }
    std::cout << std::endl;

    std::cout << "Scaffold backbone acceptors: ";
    first = true;
    for ( BBHBondActor const & actor : actors ) {
        if ( ! actor.is_donor() ) {
            if ( !first ) std::cout << ",";
            std::cout << actor.index();
            first = false;
        }
    }
    std::cout << std::endl;

    return actors;
}


std::vector<BBSasaActor>
get_bbsasa_actors( core::pose::Pose const & pose ) {

    std::vector<BBSasaActor> actors;

    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {

        std::vector<Eigen::Vector3f> sasa_points;

        {
            numeric::xyzVector<core::Real> _xyz = pose.residue(resid).nbr_atom_xyz();
            Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];
            sasa_points.push_back( xyz );
        }
        {
            numeric::xyzVector<core::Real> _xyz = pose.residue(resid).xyz("N");
            Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];
            sasa_points.push_back( xyz );
        }
        {
            numeric::xyzVector<core::Real> _xyz = pose.residue(resid).xyz("CA");
            Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];
            sasa_points.push_back( xyz );
        }
        {
            numeric::xyzVector<core::Real> _xyz = pose.residue(resid).xyz("C");
            Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];
            sasa_points.push_back( xyz );
        }     
        if ( pose.residue(resid).has("CB") ) {

            numeric::xyzVector<core::Real> _xyza = pose.residue(resid).xyz("CA");
            Eigen::Vector3f xyza; xyza[0] = _xyza[0]; xyza[1] = _xyza[1]; xyza[2] = _xyza[2];

            numeric::xyzVector<core::Real> _xyzb = pose.residue(resid).xyz("CB");
            Eigen::Vector3f xyzb; xyzb[0] = _xyzb[0]; xyzb[1] = _xyzb[1]; xyzb[2] = _xyzb[2];

            Eigen::Vector3f xyz = xyzb + (xyzb - xyza);

            sasa_points.push_back( xyz );
        }

        BBSasaActor bbs( sasa_points, resid );
        actors.push_back( bbs );

    }

    return actors;

}


void
get_default_scaffold_res( core::pose::Pose const & pose,
    utility::vector1<core::Size> & scaffold_res ) {

    for( int ir = 1; ir <= pose.size(); ++ir){
        if( !pose.residue(ir).is_protein() ) continue;
        scaffold_res.push_back(ir);
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

    if ( cache->make_bbhbond_actors ) {
        std::vector<BBHBondActor> bbhbond_actors = get_bbhbond_actors( scaffold_centered );
        for ( BBHBondActor const & bbhbond_actor : bbhbond_actors ) {
            scene.add_actor(0, bbhbond_actor );
        }
    }

    if ( cache->make_bbsasa_actors ) {
        std::vector<BBSasaActor> bbsasa_actors = get_bbsasa_actors( scaffold_centered );
        for ( BBSasaActor const & bbsasa_actor : bbsasa_actors ) {
            scene.add_actor(0, bbsasa_actor);
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

    uint64_t sanity = cache->debug_sanity + 1;
    cache->debug_sanity = sanity;

    conformation_mutable->cache_data_ = cache;
    conformation_mutable->cache_data_->conformation_is_fa = fa;

    runtime_assert( conformation->cache_data_->debug_sanity == sanity );

    std::cout << "FA status is: " << (fa ? "True" : "False") << std::endl;

    return conformation;


}

#ifdef USEHDF5
std::vector<core::pose::PoseOP>
apply_direct_segment_lookup_mover( 
    protocols::indexed_structure_store::movers::DirectSegmentLookupMover & dsl_mover,
    core::pose::Pose const & in_pose,
    uint64_t low_cut_site,
    uint64_t high_cut_site,
    uint64_t minimum_loop_length,
    uint64_t max_structures,
    uint64_t max_rmsd,
    ApplyDSLMScratch & scratch ) {

    using namespace core::pack::task::operation;
    using namespace core::select::residue_selector;
    using namespace protocols::residue_selectors;


    core::pose::Pose pose;
    utility::vector1<core::Size> rmsd_range;
    core::select::residue_selector::ResidueSubset rmsd_region(in_pose.size(), false);

    for (uint64_t i = 1; i <= in_pose.size(); i++) {
        if (i <= low_cut_site || i >= high_cut_site + 1) {
            pose.append_residue_by_bond( in_pose.residue(i) );
        } else if ( i == high_cut_site ) {
            pose.append_residue_by_jump( in_pose.residue(i), 1, "N", "N", true );
        } else {
            rmsd_region[i] = true;
        }
    }




    const std::string stored_subset_name = "inserted_lookup_segment";
    dsl_mover.stored_subset_name( stored_subset_name );
    dsl_mover.structure_store_path( basic::options::option[basic::options::OptionKeys::indexed_structure_store::fragment_store]() );


    if (scratch.scorefxn_pt.size() == 0) {
        scratch.scorefxn_pt.resize(omp_max_threads());
        scratch.to_ala_pt.resize(omp_max_threads());
        scratch.hardmin_bb_pt.resize(omp_max_threads());
        scratch.rep_scorefxn_pt.resize(omp_max_threads());

        std::cout << "Generating per-thread scaffold relaxers" << std::endl;

        for ( int i = 0; i < omp_max_threads(); i++ ) {

            //<SCOREFXNS>

            core::scoring::ScoreFunctionOP scorefxn = make_shared<core::scoring::ScoreFunction>();
            // core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("beta_soft");
            scorefxn->set_weight(core::scoring::coordinate_constraint, 2);
            scorefxn->set_weight(core::scoring::cart_bonded, 2);
            scorefxn->set_weight(core::scoring::pro_close, 0);

            //<TASKOPERATIONS>

            // RestrictAbsentCanonicalAASOP ala_only( new RestrictAbsentCanonicalAAS() );
            // ala_only->include_residue( 0 );
            // ala_only->keep_aas( "A" );

            ResidueSelectorCOP stored_residue_subset( new StoredResidueSubsetSelector( stored_subset_name ) );

            ResLvlTaskOperationOP prevent_repacking_rlt( new PreventRepackingRLT() );
            OperateOnResidueSubsetOP only_lookup_segment( new OperateOnResidueSubset( prevent_repacking_rlt, stored_residue_subset, true ) );


            protocols::simple_moves::MutateResidueOP to_ala( new protocols::simple_moves::MutateResidue() );
            to_ala->set_res_name(core::chemical::AA::aa_ala);
            to_ala->set_selector(stored_residue_subset);

            // a = protocols.simple_moves.MutateResidue()

            // a.set_res_name(core.chemical.AA.aa_ala)


            // <MOVERS>

            // core::pack::task::TaskFactoryOP to_ala_tf( new core::pack::task::TaskFactory() );
            // to_ala_tf->push_back( only_lookup_segment );
            // to_ala_tf->push_back( ala_only );

            // protocols::minimization_packing::PackRotamersMoverOP to_ala ( new protocols::minimization_packing::PackRotamersMover() );
            // to_ala->score_function(scorefxn);
            // to_ala->task_factory( to_ala_tf );


            core::pack::task::TaskFactoryOP hardmin_bb_tf( new core::pack::task::TaskFactory() );
            hardmin_bb_tf->push_back( only_lookup_segment );

            protocols::minimization_packing::MinMoverOP min_mover( new protocols::minimization_packing::MinMover() );
            min_mover->tolerance( 0.0001 );
            min_mover->min_type( "lbfgs_armijo_nonmonotone" );
            min_mover->cartesian( true );
            min_mover->score_function( scorefxn );

            protocols::minimization_packing::TaskAwareMinMoverOP hardmin_bb( new protocols::minimization_packing::TaskAwareMinMover( min_mover, hardmin_bb_tf ));
            hardmin_bb->chi( true );
            hardmin_bb->bb( true );


            //<PROTOCOLS>


            core::scoring::ScoreFunctionOP rep_scorefxn = make_shared<core::scoring::ScoreFunction>();
            rep_scorefxn->set_weight(core::scoring::fa_rep, 1);

            scratch.scorefxn_pt[i] = scorefxn;
            scratch.to_ala_pt[i] = to_ala;
            scratch.hardmin_bb_pt[i] = hardmin_bb;
            scratch.rep_scorefxn_pt[i] = rep_scorefxn;

        }
    }

    core::pose::Pose clash_check_reference = pose;
    ::devel::scheme::pose_to_ala( clash_check_reference );

    (*(scratch.rep_scorefxn_pt.front()))(clash_check_reference);
    core::scoring::Energies const & clash_check_reference_energies = clash_check_reference.energies();

    std::vector<core::scoring::EnergiesOP> clash_check_reference_energies_pt(omp_max_threads());
    std::vector<core::pose::PoseOP> in_pose_pt(omp_max_threads());
    for ( int i = 0; i < omp_max_threads(); i++ ) {
        clash_check_reference_energies_pt[i] = clash_check_reference_energies.clone();
        in_pose_pt[i] = in_pose.clone();
    }

    core::Size pose_size = pose.size();

    core::pose::PoseOP result = pose.clone();
    dsl_mover.apply( *result );

    std::vector<core::pose::PoseOP> pre_results;
    pre_results.push_back(result);

    int fetch_size = std::min((uint64_t)(omp_max_threads() * 2), max_structures * 3 / 2);
    int next_fetch = 0;
    bool out_of_results = dsl_mover.get_last_move_status() == protocols::moves::MS_FAIL_DO_NOT_RETRY;
    std::vector<core::pose::PoseOP> mid_results(fetch_size);
    uint64_t found_results = 0;
    uint64_t seen_pre_results = 0;

    while ( ! out_of_results && found_results < max_structures ) {

        std::cout << "Grabbing pre-results" << std::endl;

        next_fetch += fetch_size;
        while (pre_results.size() < next_fetch) {
            result = dsl_mover.get_additional_output();
            if (! result) {
                out_of_results = true;
                break;
            }
            core::pose::PoseOP to_insert = make_shared<core::pose::Pose>();
            to_insert->detached_copy( *result );
            pre_results.push_back(to_insert);
        } 

        int pre_results_size = pre_results.size();

        if (pre_results.size() > max_structures) {
            std::cout << "Found more fragments than max_fragments: " << max_structures << "+" << std::endl;
        } else {
            std::cout << "Only found " << pre_results_size << " fragment insertions" << std::endl;
        }


        mid_results.resize(pre_results.size()); // this should always be growing it
        std::exception_ptr exception = nullptr;

        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic,1)
        #endif
        for( int ipre = seen_pre_results; ipre < pre_results_size; ++ipre )
        {
            if (found_results >= max_structures) continue;  // burn out the loop when we're done
            try {

                int const ithread = omp_get_thread_num();

                core::pose::PoseOP iresult = pre_results[ipre];

                core::scoring::ScoreFunctionOP iscorefxn = scratch.scorefxn_pt[ithread];
                protocols::simple_moves::MutateResidueOP ito_ala = scratch.to_ala_pt[ithread];
                protocols::minimization_packing::TaskAwareMinMoverOP ihardmin_bb = scratch.hardmin_bb_pt[ithread];
                core::scoring::ScoreFunctionOP irep_scorefxn = scratch.rep_scorefxn_pt[ithread];
                core::scoring::EnergiesOP iclash_check_reference_energies = clash_check_reference_energies_pt[ithread];
                core::pose::PoseOP iin_pose = in_pose_pt[ithread];

                runtime_assert(iresult);
                runtime_assert(ito_ala);
                runtime_assert(ihardmin_bb);
                runtime_assert(irep_scorefxn);
                runtime_assert(iclash_check_reference_energies);
                runtime_assert(iin_pose);
                // runtime_assert();
                // runtime_assert();

                if ( iresult->num_chains() > 1 ) {
                    std::cout << "Broken pose" << std::endl;
                    continue;
                }
                if ( iresult->size() - pose_size < minimum_loop_length ) {
                    std::cout << "Loop too short" << std::endl;
                    continue;
                }
                ito_ala->apply( *iresult );
                ihardmin_bb->apply( *iresult );

                if ( internal_comparative_clash_check( *iclash_check_reference_energies, *iresult, irep_scorefxn,
                        8, low_cut_site - 1, high_cut_site + 1 ) ) {
                    std::cout << "Internal Clash!!" << std::endl;
                    continue;
                }

        // this doesn't handle insertions or deletion!!!!!!!!
                if ( iresult->size() == pose_size ) {
                    core::Real rmsd = subset_CA_rmsd(*iresult, *iin_pose, rmsd_region, false );
                    std::cout << "RMSD: " << rmsd << std::endl;
                    // std::cout << rmsd_region << std::endl;
                    if ( rmsd > max_rmsd ) {
                        std::cout << "RMSD too high" << std::endl;
                        continue;
                    }
                }


                mid_results[ipre] = iresult;

                #pragma omp critical
                found_results ++;

                } catch(...) {
                #pragma omp critical
                exception = std::current_exception();
            }

        } // end of OMP loop
        if( exception ) std::rethrow_exception(exception);

        seen_pre_results = pre_results.size();

    }

    // by doing it this way, we make sure we keep the correct ordering of the results
    std::vector<core::pose::PoseOP> results;
    for ( int i = 0; i < mid_results.size(); i++ ) {
        if ( mid_results[i] ) {
            results.push_back(mid_results[i]);
        }
        if (results.size() >= max_structures) break;
    }


//     do {
//         if ( result->num_chains() > 1 ) {
//             std::cout << "Broken pose" << std::endl;
//             continue;
//         }
//         if ( result->size() - pose.size() < minimum_loop_length ) {
//             std::cout << "Loop too short" << std::endl;
//             continue;
//         }
//         to_ala.apply( *result );
//         hardmin_bb.apply( *result );

//         if ( internal_comparative_clash_check( clash_check_reference_energies, *result, rep_scorefxn,
//                 8, low_cut_site - 1, high_cut_site + 1 ) ) {
//             std::cout << "Internal Clash!!" << std::endl;
//             continue;
//         }

// // this doesn't handle insertions or deletion!!!!!!!!
//         core::Real rmsd = subset_CA_rmsd(*result, in_pose, rmsd_region, false );
//         std::cout << "RMSD: " << rmsd << std::endl;
//         std::cout << rmsd_region << std::endl;
//         if ( rmsd > max_rmsd ) {
//             std::cout << "RMSD too high" << std::endl;
//             continue;
//         }

//         results.push_back( result );
//     } while ( results.size() < max_structures && ( result = dsl_mover.get_additional_output() ) );


    return results;
}
#endif


void
add_pdbinfo_if_missing( core::pose::Pose & pose ) {
    if ( ! pose.pdb_info() ) {
        core::pose::PDBInfoOP pdb_info = make_shared<core::pose::PDBInfo>( pose );
        pose.pdb_info(pdb_info);
    }
}


bool
internal_comparative_clash_check( core::scoring::Energies const & original_energies,
    core::pose::Pose const & to_check,
    core::scoring::ScoreFunctionOP scorefxn,
    core::Real max_fa_rep_delta,
    core::Size low_position,
    core::Size high_position ) {

    core::pose::Pose pose = to_check;
    ::devel::scheme::pose_to_ala(pose);

    (*scorefxn)(pose);

    core::scoring::Energies const & energies = pose.energies();

    for ( core::Size i = low_position; i <= high_position; i++ ) {
        core::Real original = original_energies.residue_total_energies(i)[core::scoring::fa_rep];
        core::Real new_score = energies.residue_total_energies(i)[core::scoring::fa_rep];

        if ( new_score - original > max_fa_rep_delta ) {
            return true;
        }
    }

    return false;

}


std::vector<core::pose::PoseOP>
extract_poses_from_silent_file( std::string const & filename ) {

    std::cout << "Reading poses from silent file: " << filename << std::endl;


    core::io::silent::SilentFileOptions options;
    core::io::silent::SilentFileData sfd( options );
    sfd.read_file( filename );

    std::cout << "Turning structures into poses" << std::endl;

    std::vector<std::reference_wrapper<core::io::silent::SilentStruct const>> silent_structs;
    std::vector<std::string> tags;
    for ( std::string const & tag : sfd.tags() ) {
        silent_structs.push_back( sfd.get_structure(tag) );
        tags.push_back(tag);
    }

    std::vector<core::pose::PoseOP> poses( silent_structs.size() );

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for ( uint64_t i = 0; i < silent_structs.size(); i++ ) {
        if (exception) continue;
        try {
            core::io::silent::SilentStruct const & ss = silent_structs[i];
            core::pose::PoseOP pose( new core::pose::Pose() );
            ss.fill_pose( *pose );
            add_pdbinfo_if_missing( *pose );
            pose->pdb_info()->name( tags[i] );
            poses[i] = pose;
        } catch(...) {
            #pragma omp critical
            exception = std::current_exception();
        }
    } // end of OMP loop
    if( exception ) std::rethrow_exception(exception);


    // for ( std::string const & tag : sfd.tags() ){
    //     core::io::silent::SilentStruct const & ss = sfd.get_structure(tag);
    //     core::pose::PoseOP pose ( new core::pose::Pose() );
    //     ss.fill_pose( *pose );
    //     add_pdbinfo_if_missing( *pose );
    //     pose->pdb_info()->name( tag );
    //     poses.push_back( pose );
    // }

    return poses;
}


std::string
pdb_name( std::string const & fname ) {
    std::string name = utility::file_basename( fname );
    if (utility::endswith( name, ".gz" ) ) name = name.substr(0, name.length() - 3);
    if (utility::endswith( name, ".pdb" ) ) name = name.substr(0, name.length() - 4);
    return name;
}


std::shared_ptr< std::vector< std::vector<float> > >
load_rotboltz_data(
    std::string const & fname,
    size_t num_res,
    size_t num_rots,
    bool ignore_missing_rots
) {

    std::cout << "Loading rotamer boltzmann data from: " << fname << std::endl;

    std::shared_ptr< std::vector< std::vector<float> > > rotboltz_data_p = std::make_shared< std::vector< std::vector<float> > >(num_res);
    std::vector< std::vector<float> > & rotboltz_data = *rotboltz_data_p;
        
    runtime_assert_msg(utility::file::file_exists( fname ), "rotboltz_data file does not exist: " + fname );
    std::ifstream in;
    in.open(fname, std::ios::in);
    // utility::io::izstream in(fname);
    std::string s;
    while (std::getline(in, s)) {
        std::string save_s = s;
        if (s.empty()) continue;

        s = utility::replace_in( s, "#", " #");
        utility::vector1<std::string> comment_splt = utility::string_split_simple(s, ' ');
        utility::vector1<std::string> splt;
        for ( std::string item : comment_splt ) {
            item = utility::strip(item, " \t\n");
            if (item.empty()) continue;
            if (item[0] == '#') break;
            splt.push_back(item);
        }
        if (splt.size() == 0) continue;

        // std::cout << save_s << std::endl;
        // std::cout << splt << std::endl;
        size_t seqpos = utility::string2int(splt[1]);
        if ( seqpos > num_res ) {
            std::cout << "Error!!! rotboltz_file has seqpos above pose.size(): saw " << seqpos
            << " pose.size() " << num_res << " in file " << fname << std::endl;
            utility_exit();
        }

        if ( rotboltz_data[seqpos-1].size() > 0 ) {
            std::cout << "Error!!! rotboltz_file has seqpos listed twice: " << seqpos
            << " in file " << fname << std::endl;
            utility_exit();
        }

        rotboltz_data[seqpos-1].resize(num_rots);
        for ( size_t i = 0; i < num_rots ; i++ ) {

            float value = 0;
            if ( i + 2 > splt.size() ) {
                if ( ! ignore_missing_rots ) {
                    utility_exit_with_message("-rotamer_boltzmann_files has wrong number of rotamers for this rifdock!! " +
                        boost::str(boost::format(" saw %i expected %i")%(splt.size()-1)%num_rots));
                }
            } else {
                value = utility::string2float(splt[i+2]);
            }
            rotboltz_data[seqpos-1][i] = value;
        }

    }

    return rotboltz_data_p;
    
}


void
load_pssm_data(
    ExtraScaffoldData & extra_data,
    std::string const & pssm_fname,
    utility::vector1<core::Size> const & scaffold_res,
    shared_ptr< RotamerIndex > rot_index_p,
    float pssm_weight,
    float pssm_cutoff,
    bool pssm_higher_is_better,
    bool pssm_enforce_no_ala
) {

    using core::Size;

    core::sequence::SequenceProfile profile_obj;
    profile_obj.read_from_file( pssm_fname );

    utility::vector1< utility::vector1< core::Real > > const & profile = profile_obj.profile();

    // The PSSM code is so terribly written that the alphabet it produces is not in the right order...
    // utility::vector1< std::string > alphabet = profile_obj.alphabet();
    utility::vector1< std::string > name3_alphabet;
    for ( Size i = 1; i <= 20; i++ ) {
        name3_alphabet.push_back(core::chemical::name_from_aa( core::chemical::AA(i) ) );
        // name3_alphabet.push_back(core::chemical::name_from_aa(core::chemical::aa_from_oneletter_code( alphabet[i][0] )).substr(0, 3) );
    }

    shared_ptr< std::vector< std::vector<float> > > per_rotamer_custom_energies_p = 
                                                        std::make_shared<std::vector< std::vector<float> >>( scaffold_res.size() );
    shared_ptr<std::vector<std::vector<bool>>> allowed_irot_at_ires_p = nullptr;
    shared_ptr<std::vector<bool>> ala_disallowed_p = nullptr;
    if ( ! std::isnan( pssm_cutoff ) ) {
        allowed_irot_at_ires_p = std::make_shared<std::vector<std::vector<bool>>>( scaffold_res.size() );
        if ( pssm_enforce_no_ala ) {
            ala_disallowed_p = make_shared<std::vector<bool>>( scaffold_res.size() );
        }
    } else {
        if ( pssm_enforce_no_ala ) {
            utility_exit_with_message("If you want -pssm_enforce_no_ala to work, you must pass a value for -pssm_cutoff");
        }
    }

    for ( Size ipos = 1; ipos <= scaffold_res.size(); ipos++ ) {
        Size ir = scaffold_res[ipos];

        if ( ir > profile.size() ) {
            utility_exit_with_message("Error: No PSSM values specified for position " + utility::to_string( ir ) );
        }

        utility::vector1< core::Real > const & pssm_row = profile[ir];

        per_rotamer_custom_energies_p->at(ipos-1).resize(rot_index_p->size());
        if ( allowed_irot_at_ires_p ) allowed_irot_at_ires_p->at(ipos-1).resize(rot_index_p->size(), true);


        for ( Size iletter = 1; iletter <= name3_alphabet.size(); iletter ++ ) {
            // index_bounds is broken
            // std::pair<int, int> index_bounds = rot_index_p->index_bounds( name3_alphabet[iletter] );
            float pssm_raw_value = pssm_row[ iletter ];
            float pssm_scaled_value = ( pssm_higher_is_better ? -1 : 1 ) * pssm_weight * pssm_raw_value;
            bool is_ok = (! allowed_irot_at_ires_p) || ( pssm_higher_is_better ? pssm_raw_value >= pssm_cutoff : pssm_raw_value <= pssm_cutoff );

            // std::cout << "Pos: " << ir << " " << name3_alphabet[iletter] << " " << pssm_raw_value << " " << (is_ok ? "ok" : "not") << std::endl;
            for ( Size irot = 0; irot < rot_index_p->size(); irot++ ) {
                if ( name3_alphabet[iletter] != rot_index_p->resname( irot ) ) continue;
                
                per_rotamer_custom_energies_p->at(ipos-1)[irot] = pssm_scaled_value;
                if ( allowed_irot_at_ires_p ) allowed_irot_at_ires_p->at(ipos-1)[irot] = is_ok;
            }

            if ( ala_disallowed_p && ! is_ok && name3_alphabet[iletter] == "ALA" ) ala_disallowed_p->at( ipos-1 ) = true;


        }



    }

    extra_data.accumulate_per_rotamer_custom_energies( per_rotamer_custom_energies_p );
    if ( allowed_irot_at_ires_p ) extra_data.accumulate_allowed_irot_at_ires( allowed_irot_at_ires_p );
    if ( ala_disallowed_p ) extra_data.accumulate_ala_disallowed( ala_disallowed_p );

}


void
load_clash_context( 
    ExtraScaffoldData & extra_data,
    std::string const & clash_fname,
    shared_ptr< RotamerIndex > rot_index_p
) {
    core::pose::Pose pose;
    core::import_pose::pose_from_file( pose, clash_fname );

    std::vector<SimpleAtom> clash_simple_atoms;

    for( int ir = 1; ir <= pose.size(); ++ir ){
        utility::vector1<core::Size> resids(1,ir); // 1-index numbering
        {
            std::vector<SchemeAtom> scaff_res_atoms;
            devel::scheme::get_scheme_atoms( pose, resids, scaff_res_atoms, true ); //bb + CB
            
            int restype = rot_index_p->chem_index_.resname2num( pose.residue(ir).name3() ); // for UNK will be -1
            for( int ia = 0; ia < scaff_res_atoms.size(); ++ia){
                SchemeAtom const & a( scaff_res_atoms[ia] );
                runtime_assert( a.type() > 0 );
                if( a.type() >= 21 ) continue;
                SimpleAtom sa( a.position(), a.type(), restype, ia );
                clash_simple_atoms.push_back(sa);
            }
        }
    }
    extra_data.add_to_clash_context( clash_simple_atoms );
}


}}

