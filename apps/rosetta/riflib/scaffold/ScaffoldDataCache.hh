// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_scaffold_ScaffoldDataCache_hh
#define INCLUDED_riflib_scaffold_ScaffoldDataCache_hh


#include <scheme/types.hh>
#include <scheme/nest/NEST.hh>
#include <scheme/objective/storage/TwoBodyTable.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rosetta_field.hh>
#include <riflib/RotamerGenerator.hh>
#include <riflib/util.hh>
#include <rif_dock_test.hh>
#include <riflib/rotamer_energy_tables.hh>
#include <riflib/scaffold/MultithreadPoseCloner.hh>
#include <riflib/scaffold/util.hh>
#include <riflib/scaffold/ExtraScaffoldData.hh>
#include <riflib/HSearchConstraints.hh>
#include <riflib/BurialManager.hh>

#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <boost/format.hpp>



namespace devel {
namespace scheme {


struct ScaffoldDataCache {

// setup during constructor
    std::string scafftag;
    std::string scaff_fname;

    shared_ptr<utility::vector1<core::Size>> scaffold_res_p;                   // Seqposs of residues to design, default whole scaffold

    shared_ptr<std::vector<int>> scaffres_g2l_p;                               // maps global_seqpos -> local_seqpos  (local_seqpos.size() == scaffold_res.size())
    shared_ptr<std::vector<int>> scaffres_l2g_p;                               // maps local_seqpos  -> global_seqpos
    shared_ptr<std::vector<bool>> scaffuseres_p;                                 // maps global_seqpos -> (bool)being used

    Eigen::Vector3f scaffold_center;
    float scaff_redundancy_filter_rg;
    float scaff_radius;
    EigenXform scaffold_perturb;

    core::pose::PoseCOP scaffold_unmodified_p;                                     // original pose passed to this function
    core::pose::PoseCOP scaffold_centered_p;                                   // centered scaffold, some ALA mutations based on scaff2ala/scaff2alaselonly
    core::pose::PoseCOP scaffold_full_centered_p;                              // centered scaffold, identical to input       


    shared_ptr<std::vector< SimpleAtom >> scaffold_simple_atoms_all_p;         // all the atoms of scaffold_centered as SimpleAtoms
    shared_ptr<std::vector< SimpleAtom >> scaffold_simple_atoms_p;             // This can have a few different values
                                                                               // if opt.lowres_sterics_cbonly:
                                                                               //     all CB of scaffold as SimpleAtoms
                                                                               // else:
                                                                               //     N CA C CB for res in scaffold_res
                                                                               //     CB for others

    shared_ptr<std::vector<std::string>> scaffold_sequence_glob0_p;            // Scaffold sequence in name3 space
    shared_ptr<std::vector< std::pair<int,int> > > local_rotamers_p;           // lower and upper bounds into rotamer_index for each local_seqpos
    std::string scaff_res_hashstr;
    shared_ptr<std::vector<std::pair<core::Real,core::Real> > > scaffold_phi_psi_p; // Scaffold phi-psi
    shared_ptr<std::vector<bool>> scaffold_d_pos_p;                                  // Scaffold allow d postion base on phi-psi
    uint64_t debug_sanity;

    shared_ptr<std::vector<std::vector<bool>>> allowed_irot_at_ires_p;
    shared_ptr<std::vector<bool>> ala_disallowed_p; // More extreme than allowed_irot_at_res. Causes extra rotamers to be inserted during packing

    std::shared_ptr< std::vector< std::vector<float> > > per_rotamer_custom_energies_p;


// not setup during constructor
    shared_ptr<std::vector<std::vector<float> > > scaffold_onebody_glob0_p;    //onebodies in global numbering
    shared_ptr<std::vector<std::vector<float> > > local_onebody_p;       //onebodies in local numbering

    typedef ::scheme::objective::storage::TwoBodyTable<float> TBT;

    shared_ptr<TBT> scaffold_twobody_p;                                        // twobody_rotamer_energies using global_seqpos
    shared_ptr<TBT> local_twobody_p;                                           // twobody_rotamer_energies using local_seqpos

    std::vector<shared_ptr<TBT>> local_twobody_per_thread;                     // Used with BuriedUnsats, these can momentarily change but must be reset


    MultithreadPoseCloner mpc_both_pose;                                       // scaffold_centered_p + target
    MultithreadPoseCloner mpc_both_full_pose;                                  // scaffold_full_centered_p + target

    std::vector<CstBaseOP> csts;

    shared_ptr<BurialVoxelArray> burial_grid;

// Conformation state
    bool conformation_is_fa;


// Items that don't belong here but nontheless are here
    bool make_bbhbond_actors;                                                  // needed in make_conformation_from_data_cache
    bool make_bbsasa_actors;                                                   // needed in make_conformation_from_data_cache



    ScaffoldDataCache() {}

    // Don't add any more fields to this call
    // Consider instead that you probably want to add something to ExtraScaffoldData
    //  (you'll have to modify fewer functions if you do it this way)
    ScaffoldDataCache( core::pose::Pose & pose, 
        utility::vector1<core::Size> const & scaffold_res_in, 
        std::string const &scafftag_in,
        EigenXform const & scaffold_perturb_in,
        shared_ptr< RotamerIndex > rot_index_p,
        RifDockOpt const & opt,
        ExtraScaffoldData const & extra_data_in,
	    std::string const & scaffold_name
    ) {

        ExtraScaffoldData extra_data = extra_data_in;

        debug_sanity = 1337;

        scaffold_res_p = make_shared<utility::vector1<core::Size>>(scaffold_res_in);
        scafftag = scafftag_in;
        scaff_res_hashstr = ::devel::scheme::get_res_list_hash( *scaffold_res_p );
        scaffold_perturb = scaffold_perturb_in;

	    scaff_fname = scaffold_name;

        typedef numeric::xyzVector<core::Real> Vec;


        // setting up scaffres_g2l_p, scaffres_l2g_p, and scaffuseres_p
        int count = 0;
        scaffres_l2g_p = make_shared<std::vector<int>>();
        scaffres_g2l_p = make_shared<std::vector<int>>(pose.size(), -1);
        scaffuseres_p = make_shared<std::vector<bool>>(pose.size(), false);
        // scaffold_phi_psi_p = make_shared<std::vector< std::pair<core::Real,core::Real>>>();
        // scaffold_d_pos_p = make_shared<std::vector<int>>();
        

        // setting up the d and l maps
        std::vector<bool> d_map = rot_index_p->is_d_;
        std::vector<bool> l_map( d_map.size() );
        for ( int irot = 0; irot < l_map.size(); irot++ ) l_map[irot] = ! d_map[irot]; 

        // I did it, I commented it out. End the suffering forever -bcov
        // for ( int irot = 0; irot < l_map.size(); irot++ ) {
        //     std::cout << irot << " " << ( l_map[irot] ? "L" : " " ) << ( d_map[irot] ? "D" : " " ) << rot_index_p->oneletter(irot) << std::endl;
        // }

        if ( opt.use_dl_mix_bb ) {
            shared_ptr<std::vector<std::vector<bool>>> dl_allowed_irot_at_pos = make_shared<std::vector<std::vector<bool>>>();
            for ( core::Size ir : *scaffold_res_p) {
                if (pose.phi(ir) > 0) {
                    dl_allowed_irot_at_pos->push_back( d_map );
                } else {
                    dl_allowed_irot_at_pos->push_back( l_map );
                }
            }
            extra_data.accumulate_allowed_irot_at_ires( dl_allowed_irot_at_pos );
        }



        for ( core::Size ir : *scaffold_res_p) {
            (*scaffres_g2l_p)[ir-1] = count++;
            scaffres_l2g_p->push_back(ir-1);
            (*scaffuseres_p)[ir-1] = true;

        }




        // This is setting scaff_redundancy_filter_rg and scaff_radius
        get_rg_radius( pose, scaff_redundancy_filter_rg, scaff_radius, *scaffold_res_p, false ); 

        // Setup scaffold_unmodified_p
        scaffold_unmodified_p = make_shared<core::pose::Pose const>( pose );

        // Setup scaffold_centered_p and scaffold_full_centered_p
        scaffold_centered_p = make_shared<core::pose::Pose const>( pose );
        core::pose::Pose & scaffold_centered = const_cast<core::pose::Pose &>( *scaffold_centered_p );
        add_pdbinfo_if_missing( scaffold_centered );
        scaffold_full_centered_p = make_shared<core::pose::Pose const>( pose );
        core::pose::Pose & scaffold_full_centered = const_cast<core::pose::Pose &>( *scaffold_full_centered_p );
        add_pdbinfo_if_missing( scaffold_full_centered );

        if     ( opt.scaff2ala )        ::devel::scheme::pose_to_ala( scaffold_centered );
        else if( opt.scaff2alaselonly ) ::devel::scheme::pose_to_ala( scaffold_centered, *scaffold_res_p );

        // Setup scaffold_center
        if ( std::isnan( extra_data.force_scaffold_center[0] ) || std::isnan( extra_data.force_scaffold_center[1] )
            || std::isnan( extra_data.force_scaffold_center[2] )) {
            scaffold_center = pose_center(scaffold_centered,*scaffold_res_p);
        } else {
            scaffold_center = extra_data.force_scaffold_center;
        }

        // Move scaffold_centered_p and scaffold_full_centered_p to origin
        for( int ir = 1; ir <= scaffold_full_centered.size(); ++ir ){
            Vec tmp( scaffold_center[0], scaffold_center[1], scaffold_center[2] );
            for( int ia = 1; ia <= scaffold_centered.residue_type(ir).natoms(); ++ia ){
                core::id::AtomID aid(ia,ir);
                scaffold_centered.set_xyz( aid, scaffold_centered.xyz(aid) - tmp );
            }
            for( int ia = 1; ia <= scaffold_full_centered.residue_type(ir).natoms(); ++ia ){
                core::id::AtomID aid(ia,ir);
                scaffold_full_centered.set_xyz( aid, scaffold_full_centered.xyz(aid) - tmp );
            }
        }

        // Setup scaffold_sequence_glob0
        scaffold_sequence_glob0_p = make_shared<std::vector<std::string>>();
        for( int ir = 1; ir <= pose.size(); ++ir ){
            scaffold_sequence_glob0_p->push_back( pose.residue(ir).name3() );
        }


        // Setup local rotamers
        local_rotamers_p = make_shared<std::vector< std::pair<int,int>>>();
        for( int i = 0; i < scaffold_res_p->size(); ++i ){
            int iresglobal = scaffres_l2g_p->at(i);
            std::string name3 = scaffold_sequence_glob0_p->at(iresglobal);

            // This won't add custom rot_index_spec rotamers but I don't care
            std::pair<int,int> ib = rot_index_p->index_bounds( name3 );
            // std::cout << "local_rotamers " << i << " " << iresglobal << " " << name3 << " " << ib.first << " " << ib.second << std::endl;
            local_rotamers_p->push_back( ib );
        }


        // Setup scaffold_simple_atoms_p, and scaffold_simple_atoms_all_p

        scaffold_simple_atoms_p = make_shared<std::vector< SimpleAtom >>();
        scaffold_simple_atoms_all_p = make_shared<std::vector< SimpleAtom >>();

        for( int ir = 1; ir <= scaffold_centered.size(); ++ir ){
            utility::vector1<core::Size> resids(1,ir); // 1-index numbering
            {
                std::vector<SchemeAtom> scaff_res_atoms;
                if( !opt.lowres_sterics_cbonly && std::find( scaffold_res_p->begin(), scaffold_res_p->end(), ir ) != scaffold_res_p->end() ){
                    devel::scheme::get_scheme_atoms( scaffold_centered, resids, scaff_res_atoms, true ); //bb + CB
                } else { // is not selected residue
                    devel::scheme::get_scheme_atoms_cbonly( scaffold_centered, resids, scaff_res_atoms ); // literally only CB
                }
                int restype = rot_index_p->chem_index_.resname2num( scaffold_centered.residue(ir).name3() ); // for UNK will be -1
                for( int ia = 0; ia < scaff_res_atoms.size(); ++ia){
                    SchemeAtom const & a( scaff_res_atoms[ia] );
                    runtime_assert( a.type() > 0 );
                    if( a.type() >= 21 ) continue;
                    SimpleAtom sa( a.position(), a.type(), restype, ia );
                    scaffold_simple_atoms_p->push_back(sa);
                }
            }
            {
                std::vector<SchemeAtom> all_scaff_res_atoms;
                devel::scheme::get_scheme_atoms( scaffold_centered, resids, all_scaff_res_atoms, false );
                int restype = rot_index_p->chem_index_.resname2num( scaffold_centered.residue(ir).name3() ); // for UNK will be -1
                for( int ia = 0; ia < all_scaff_res_atoms.size(); ++ia){
                    SchemeAtom const & a( all_scaff_res_atoms[ia] );
                    runtime_assert( a.type() > 0 );
                    if( a.type() >= 21 ) continue;
                    SimpleAtom sa( a.position(), a.type(), restype, ia );
                    scaffold_simple_atoms_all_p->push_back(sa);
                }
            }
        }

        make_bbhbond_actors = opt.scaff_bb_hbond_weight > 0;
        make_bbsasa_actors = opt.need_to_calculate_sasa;

        std::cout << "scaffold selected region rg: " << scaff_redundancy_filter_rg << ", radius: " << scaff_radius << std::endl;
        std::cout << "scaffold_simple_atoms " << scaffold_simple_atoms_p->size() << std::endl;


        for(CstBaseOP cst_in : extra_data.csts ) { 
            CstBaseOP cst = cst_in->clone(); // make a copy
            cst->reset();   // clear previous scaffold related data
            csts.push_back(cst);
        }

        if( opt.native_docking ){
            std::cout << "KILLING NON-NATIVE ROTAMERS ON SCAFFOLD AND ALA!!!" << std::endl;
            shared_ptr<std::vector<std::vector<bool>>> native_allowed_irot_at_pos = make_shared<std::vector<std::vector<bool>>>();
            shared_ptr<std::vector<bool>> native_ala_disallowed = make_shared<std::vector<bool>>();
            for ( core::Size ir : *scaffold_res_p) {
                std::vector<bool> this_pos( rot_index_p->size(), true );
                for( int irot = 0; irot < rot_index_p->size(); ++irot ){
                    if( rot_index_p->resname(irot) != scaffold_sequence_glob0_p->at(ir-1) ){
                        this_pos[irot] = false;
                    }
                }
                native_allowed_irot_at_pos->push_back(this_pos);
                native_ala_disallowed->push_back( scaffold_sequence_glob0_p->at(ir-1) != "ALA" );
            }
            extra_data.accumulate_allowed_irot_at_ires( native_allowed_irot_at_pos );
            extra_data.accumulate_ala_disallowed( native_ala_disallowed );
        }


        if( opt.restrict_to_native_scaffold_res ){
            std::cout << "KILLING NON-NATIVE ROTAMERS ON SCAFFOLD!!!" << std::endl;
            shared_ptr<std::vector<std::vector<bool>>> native_allowed_irot_at_pos = make_shared<std::vector<std::vector<bool>>>();
            for ( core::Size ir : *scaffold_res_p) {
                std::vector<bool> this_pos( rot_index_p->size(), true );
                for( int irot = 0; irot < rot_index_p->size(); ++irot ){
                    if( rot_index_p->resname(irot) != scaffold_sequence_glob0_p->at(ir-1) && rot_index_p->resname(irot) != "ALA" ){
                        this_pos[irot] = true;
                    }
                }
                native_allowed_irot_at_pos->push_back(this_pos);
            }
            extra_data.accumulate_allowed_irot_at_ires( native_allowed_irot_at_pos );
        }

        if( opt.bonus_to_native_scaffold_res != 0 ){
            std::cout << "adding to native scaffold res 1BE " << opt.bonus_to_native_scaffold_res << std::endl;
            std::shared_ptr< std::vector< std::vector<float> > > native_bonuses_p = make_shared< std::vector< std::vector<float> > >();
            for ( core::Size ir : *scaffold_res_p) {
                std::vector<float> this_bonus( rot_index_p->size(), 0 );
                for( int irot = 0; irot < rot_index_p->size(); ++irot ){
                    if( rot_index_p->resname(irot) == scaffold_sequence_glob0_p->at(ir-1) ){
                        this_bonus[irot] += opt.bonus_to_native_scaffold_res;
                    }
                }
                native_bonuses_p->push_back(this_bonus);
            }
            extra_data.accumulate_per_rotamer_custom_energies( native_bonuses_p );
        }




        per_rotamer_custom_energies_p = extra_data.per_rotamer_custom_energies_p;
        allowed_irot_at_ires_p = extra_data.allowed_irot_at_ires_p;
        ala_disallowed_p = extra_data.ala_disallowed_p;

        if ( per_rotamer_custom_energies_p ) runtime_assert( scaffold_res_p->size() == per_rotamer_custom_energies_p->size() );
        if ( allowed_irot_at_ires_p ) runtime_assert( scaffold_res_p->size() == allowed_irot_at_ires_p->size() );
        if ( ala_disallowed_p ) runtime_assert( scaffold_res_p->size() == ala_disallowed_p->size() );

        if ( ala_disallowed_p ) {
            runtime_assert( allowed_irot_at_ires_p );
            for ( core::Size i = 0; i < ala_disallowed_p->size(); i++ ) {
                if ( ala_disallowed_p->at(i) ) {
                    // We assert this here because it will introduce a bug in hackpack
                    runtime_assert( ! allowed_irot_at_ires_p->at(i)[ rot_index_p->ala_rot() ] );
                    bool any = false;
                    for ( core::Size irot = 0; irot < rot_index_p->size(); irot++ ) {
                        if ( allowed_irot_at_ires_p->at(i)[irot] ) {
                            any = true;
                            break;
                        }
                    }
                    if ( ! any ) {
                        utility_exit_with_message("No rotamers allowed at scaffold position " + utility::to_string(scaffold_res_p->at(i))
                                + " and the default ALA residue has been specifically disabled.");
                    }
                }
            }
        }

        // if ( allowed_irot_at_ires_p ) {

        //     for ( int i = 1; i <= scaffold_res_p->size(); i++) {

        //         core::Size ir = scaffold_res_p->at(i);
     
        //         std::vector<bool> row = allowed_irot_at_ires_p->at(i-1);

        //         std::set<std::string> allowed;
        //         for ( int irot = 0; irot < rot_index_p->size(); irot++ ) {
        //             if ( row.at(irot) ) {
        //                 allowed.insert( rot_index_p->oneletter( irot ) );
        //             }
        //         }

        //         std::string my_str = "";
        //         for ( std::string const & str : allowed ) {
        //             my_str += str;
        //         }

        //         std::cout << "Allowed: " << ir << my_str << std::endl;
        //     }
        // }


    }

    // returns true if there are more than 0 constraints
    bool
    prepare_contraints( core::pose::Pose const & target, double resl ) {
        bool any_csts = false;
        for(CstBaseOP cst : csts ) { 
            cst->set_coordinates(target, *scaffold_centered_p);
            cst->init(target, *scaffold_centered_p, resl);
            any_csts = true;
        }
        return any_csts;
    }


    void
    setup_fake_onebody_tables(
        shared_ptr< RotamerIndex> rot_index_p,
        RifDockOpt const & opt ) {


        RotamerIndex & rot_index = *rot_index_p;

        scaffold_onebody_glob0_p = make_shared<std::vector<std::vector<float> >>(scaffold_centered_p->size());
        

        for( int ir = 1; ir <= scaffold_centered_p->size(); ++ir ){

            (*scaffold_onebody_glob0_p)[ir-1].resize( rot_index.size(), 12345.0 );

            if( std::find(scaffold_res_p->begin(), scaffold_res_p->end(), ir) == scaffold_res_p->end() ){
                continue;
            }
            
            if( ! scaffold_centered_p->residue(ir).is_protein()   ) continue;
            if(   scaffold_centered_p->residue(ir).name3()=="GLY" ) continue;
            if(   scaffold_centered_p->residue(ir).name3()=="PRO" ) continue;
            for( int jr = 0; jr < rot_index.size(); ++jr ) { 
                (*scaffold_onebody_glob0_p)[ir-1][jr] = 0;
            }

            local_onebody_p = make_shared<std::vector<std::vector<float> > >();
            for( int i = 0; i < scaffres_l2g_p->size(); ++i ){
                local_onebody_p->push_back( scaffold_onebody_glob0_p->at( scaffres_l2g_p->at(i) ) );
            }

            for( int i = 0; i < scaffres_g2l_p->size(); ++i ){
                if( (*scaffres_g2l_p)[i] < 0 ){
                    BOOST_FOREACH( float & f, (*scaffold_onebody_glob0_p)[i] ) f = 9e9;
                }
            }
        }
       
    }

    // setup scaffold_onebody_glob0_p and local_onebody_p
    void
    setup_onebody_tables(
        shared_ptr< RotamerIndex > rot_index_p,
        RifDockOpt const & opt ) {

        if (local_onebody_p) return;

        scaffold_onebody_glob0_p = make_shared<std::vector<std::vector<float> >>();

        std::string cachefile_1be = "__1BE_"+scafftag+(opt.replace_all_with_ala_1bre?"_ALLALA":"")+"_reshash"+scaff_res_hashstr+".bin.gz";
        if( ! opt.cache_scaffold_data ) cachefile_1be = "";
        std::cout << "rifdock: get_onebody_rotamer_energies" << std::endl;
        get_onebody_rotamer_energies(
                *scaffold_centered_p,
                *scaffold_res_p,           // uses 12345 as score for anything missing here
                *rot_index_p,
                *scaffold_onebody_glob0_p,
                opt.data_cache_path,
                cachefile_1be,
                opt.replace_all_with_ala_1bre,
                opt.favorable_1body_multiplier,
                opt.favorable_1body_multiplier_cutoff,
                per_rotamer_custom_energies_p
            );

        // Handled above to make things more streamlined
        // if( opt.restrict_to_native_scaffold_res ){
        //     std::cout << "KILLING NON-NATIVE ROTAMERS ON SCAFFOLD!!!" << std::endl;
        //     for( int ir = 0; ir < scaffold_onebody_glob0_p->size(); ++ir ){
        //         for( int irot = 0; irot < rot_index_p->size(); ++irot ){
        //             if( rot_index_p->resname(irot) != scaffold_sequence_glob0_p->at(ir) && rot_index_p->resname(irot) != "ALA" ){
        //                 (*scaffold_onebody_glob0_p)[ir][irot] = 9e9;
        //             }
        //         }
        //     }
        // }
        // Handled above to make things more streamlined
        // if( opt.bonus_to_native_scaffold_res != 0 ){
        //     std::cout << "adding to native scaffold res 1BE " << opt.bonus_to_native_scaffold_res << std::endl;
        //     for( int ir = 0; ir < scaffold_onebody_glob0_p->size(); ++ir ){
        //         for( int irot = 0; irot < rot_index_p->size(); ++irot ){
        //             if( rot_index_p->resname(irot) == scaffold_sequence_glob0_p->at(ir) ){
        //                 (*scaffold_onebody_glob0_p)[ir][irot] += opt.bonus_to_native_scaffold_res;
        //             }
        //         }
        //     }
        // }

        local_onebody_p = make_shared<std::vector<std::vector<float> > >();
        for( int i = 0; i < scaffres_l2g_p->size(); ++i ){
            local_onebody_p->push_back( scaffold_onebody_glob0_p->at( scaffres_l2g_p->at(i) ) );
        }

        for( int i = 0; i < scaffres_g2l_p->size(); ++i ){
            if( (*scaffres_g2l_p)[i] < 0 ){
                BOOST_FOREACH( float & f, (*scaffold_onebody_glob0_p)[i] ) f = 9e9;
            }
        }
    }

    // setup scaffold_twobody_p and local_twobody_p
    void
    setup_twobody_tables(  
        shared_ptr< RotamerIndex > rot_index_p,
        RifDockOpt const & opt ,
        MakeTwobodyOpts const & make2bopts,
        ::devel::scheme::RotamerRFTablesManager & rotrf_table_manager) {

        if (local_twobody_p) return;

        scaffold_twobody_p = make_shared<TBT>( scaffold_centered_p->size(), rot_index_p->size()  );

        std::string energy_cut = boost::str(boost::format("_ecut_%.2f")%opt.rotamer_onebody_inclusion_threshold);
        

        std::cout << "rifdock: get_twobody_tables" << std::endl;
        std::string cachefile2b = "__2BE_" + scafftag + "_reshash" + scaff_res_hashstr + energy_cut + ".bin.gz";
        if( ! opt.cache_scaffold_data || opt.extra_rotamers ) cachefile2b = "";
        std::string dscrtmp;
        get_twobody_tables(
                opt.data_cache_path,
                cachefile2b,
                dscrtmp,
                *scaffold_centered_p,
                *rot_index_p,
                *scaffold_onebody_glob0_p,
                rotrf_table_manager,
                make2bopts,
                *scaffold_twobody_p
            );


        local_twobody_p = scaffold_twobody_p->create_subtable( *scaffuseres_p, *scaffold_onebody_glob0_p, make2bopts.onebody_threshold );
    

        std::cout << "rifdock: twobody memuse: " << (float)scaffold_twobody_p->twobody_mem_use()/1000.0/1000.0 << "M" << std::endl;
        std::cout << "rifdock: onebody dimension: " << scaffold_onebody_glob0_p->size() << " " << scaffold_onebody_glob0_p->front().size() << std::endl;
        int onebody_n_allowed = 0;
        for( auto const & t : *(scaffold_onebody_glob0_p) ){
            for( auto const & v : t ){
                if( v < make2bopts.onebody_threshold ) onebody_n_allowed++;
            }
        }
        std::cout << "rifdock: onebody Nallowed: " << onebody_n_allowed << std::endl;
        std::cout << "filt_2b memuse: " << (float)local_twobody_p->twobody_mem_use()/1000.0/1000.0 << "M" << std::endl;
    }



    // setup scaffold_twobody_p and local_twobody_p
    void
    setup_twobody_tables_per_thread( ) {
        runtime_assert( local_twobody_p );

        if ( local_twobody_per_thread.size() > 0 ) return;

        local_twobody_per_thread.resize(::devel::scheme::omp_max_threads_1());

        for ( int i = 0; i < local_twobody_per_thread.size(); i++ ) {
            local_twobody_per_thread[i] = local_twobody_p->clone();
        }
    }

    bool
    check_twobody_tables_per_thread_integrity( ) {
        runtime_assert( local_twobody_per_thread.size() > 0 );

        local_twobody_per_thread.resize(::devel::scheme::omp_max_threads_1());

        for ( int i = 0; i < local_twobody_per_thread.size(); i++ ) {
            if ( ! local_twobody_per_thread[i]->check_equal( *local_twobody_p) ) {
                return false;
            }
        }
        return true;
    }


    float
    get_redundancy_filter_rg( float target_redundancy_filter_rg ) {
        return std::min( target_redundancy_filter_rg, scaff_redundancy_filter_rg );
    }



    void
    setup_both_pose( core::pose::Pose const & target ) {
        if ( mpc_both_pose.size() > 0 ) return;
        mpc_both_pose.add_pose(helper_setup_both_pose( target ));

        runtime_assert( mpc_both_pose.get_pose()->size() == scaffold_centered_p->size() + target.size() );
    }

    core::pose::PoseOP
    helper_setup_both_pose( core::pose::Pose const & target ) {
        core::pose::PoseOP __both_pose_p = make_shared<core::pose::Pose>( *scaffold_centered_p );
        ::devel::scheme::append_pose_to_pose( *__both_pose_p, target );
        add_pdbinfo_if_missing( *__both_pose_p );
        return __both_pose_p;
    }

    void
    setup_both_full_pose( core::pose::Pose const & target ) {
        if ( mpc_both_full_pose.size() > 0 ) return;
        mpc_both_full_pose.add_pose(helper_setup_both_full_pose(target));

        runtime_assert( mpc_both_full_pose.get_pose()->size() == scaffold_full_centered_p->size() + target.size() );
    }

    core::pose::PoseOP
    helper_setup_both_full_pose( core::pose::Pose const & target ) {
        core::pose::PoseOP __both_full_pose_p = make_shared<core::pose::Pose>( *scaffold_full_centered_p );
        ::devel::scheme::append_pose_to_pose( *__both_full_pose_p, target );
        add_pdbinfo_if_missing( *__both_full_pose_p );
        return __both_full_pose_p;
    }


    void
    setup_burial_grids( shared_ptr<BurialManager> const & burial_manager ) {
        burial_grid = burial_manager->get_scaffold_neighbors( *scaffold_centered_p );
    }


};


inline
void dump_bbhbond_actors_inner( shared_ptr<ScaffoldDataCache> const & cache, bool is_donor ) {
    std::string type = is_donor ? "donor" : "acceptor";
    std::string fname = cache->scafftag + "_BBHbondActor_" + type + "s.pdb.gz";

    std::cout << "Dumping BBHBondActor " + type + "s to: " + fname << std::endl;

    utility::io::ozstream out( fname );
    std::vector<HBondRay> hbrays;

    ParametricScene scene(1);
    scene.replace_body( 0, make_conformation_from_data_cache( cache, false ) );

    for ( int ihbactor = 0; ihbactor < scene.template num_actors<BBHBondActor>(0); ihbactor++ ) {
        BBHBondActor bbh = scene.template get_actor<BBHBondActor>(0,ihbactor);

        if ( bbh.is_donor() == is_donor ) {
            for ( HBondRay const & hbray : bbh.hbond_rays() ) {
                hbrays.push_back( hbray );
            }
        }   
    }

    EigenXform xform = EigenXform::Identity();
    xform.translation() = cache->scaffold_center;
    dump_hbond_rays( out, hbrays, is_donor, xform );
}

inline
void dump_bbhbond_actors( shared_ptr<ScaffoldDataCache> const & cache ) {
    dump_bbhbond_actors_inner( cache, true );
    dump_bbhbond_actors_inner( cache, false );
}





typedef shared_ptr<ScaffoldDataCache> ScaffoldDataCacheOP;
typedef shared_ptr<ScaffoldDataCache const > ScaffoldDataCacheCOP;

}}



#endif
