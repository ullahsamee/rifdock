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
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/rosetta_field.hh>
#include <riflib/RotamerGenerator.hh>
#include <riflib/util.hh>
#include <rif_dock_test.hh>
#include <riflib/rotamer_energy_tables.hh>

#include <core/pose/Pose.hh>
#include <utility/vector1.hh>



namespace devel {
namespace scheme {


struct ScaffoldDataCache {

// setup during constructor
    std::string scafftag;

    shared_ptr<utility::vector1<core::Size>> scaffold_res_p;                   // Seqposs of residues to design, default whole scaffold

    shared_ptr<std::vector<int>> scaffres_g2l_p;                               // maps global_seqpos -> local_seqpos  (local_seqpos.size() == scaffold_res.size())
    shared_ptr<std::vector<int>> scaffres_l2g_p;                               // maps local_seqpos  -> global_seqpos
    std::vector<bool> scaffuseres;                                             // maps global_seqpos -> (bool)being used

    Eigen::Vector3f scaffold_center;
    float scaff_redundancy_filter_rg;
    float scaff_radius;

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


// not setup during constructor
    shared_ptr<std::vector<std::vector<float> > > scaffold_onebody_glob0_p;    //onebodies in global numbering
    shared_ptr<std::vector<std::vector<float> > > local_onebody_p;       //onebodies in local numbering


    typedef ::scheme::objective::storage::TwoBodyTable<float> TBT;

    shared_ptr<TBT> scaffold_twobody_p;                                        // twobody_rotamer_energies using global_seqpos
    shared_ptr<TBT> local_twobody_p;                                           // twobody_rotamer_energies using local_seqpos


    ScaffoldDataCache() {}

    ScaffoldDataCache( core::pose::Pose & pose, 
        utility::vector1<core::Size> const & scaffold_res_in, 
        std::string const &scafftag_in, 
        shared_ptr< RotamerIndex > rot_index_p,
        RifDockOpt const & opt) {

        scaffold_res_p = make_shared<utility::vector1<core::Size>>(scaffold_res_in);
        scafftag = scafftag_in;
        scaff_res_hashstr = ::devel::scheme::get_res_list_hash( *scaffold_res_p );

        typedef numeric::xyzVector<core::Real> Vec;


        // setting up scaffres_g2l_p, scaffres_l2g_p, and scaffuseres
        int count = 0;
        scaffres_l2g_p = make_shared<std::vector<int>>();
        scaffres_g2l_p = make_shared<std::vector<int>>(pose.size(), -1);
        scaffuseres.resize(pose.size(), false);

        for ( core::Size ir : *scaffold_res_p) {
            (*scaffres_g2l_p)[ir-1] = count++;
            scaffres_l2g_p->push_back(ir-1);
            scaffuseres[ir-1] = true;
        }

        // This is setting scaff_redundancy_filter_rg and scaff_radius
        get_rg_radius( pose, scaff_redundancy_filter_rg, scaff_radius, *scaffold_res_p, false ); 


        // Setup scaffold_centered_p and scaffold_full_centered_p
        scaffold_centered_p = make_shared<core::pose::Pose const>( pose );
        core::pose::Pose & scaffold_centered = const_cast<core::pose::Pose &>( *scaffold_centered_p );
        scaffold_full_centered_p = make_shared<core::pose::Pose const>( pose );
        core::pose::Pose & scaffold_full_centered = const_cast<core::pose::Pose &>( *scaffold_full_centered_p );

        if     ( opt.scaff2ala )        ::devel::scheme::pose_to_ala( scaffold_centered );
        else if( opt.scaff2alaselonly ) ::devel::scheme::pose_to_ala( scaffold_centered, *scaffold_res_p );

        // Setup scaffold_center
        scaffold_center = pose_center(scaffold_centered,*scaffold_res_p);

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



    }


    // setup scaffold_onebody_glob0_p and local_onebody_p
    void
    setup_onebody_tables( 
        shared_ptr< RotamerIndex > rot_index_p,
        RifDockOpt const & opt ) {


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
                opt.replace_all_with_ala_1bre
            );

        if( opt.restrict_to_native_scaffold_res ){
            std::cout << "KILLING NON-NATIVE ROTAMERS ON SCAFFOLD!!!" << std::endl;
            for( int ir = 0; ir < scaffold_onebody_glob0_p->size(); ++ir ){
                for( int irot = 0; irot < rot_index_p->size(); ++irot ){
                    if( rot_index_p->resname(irot) != scaffold_sequence_glob0_p->at(ir) && rot_index_p->resname(irot) != "ALA" ){
                        (*scaffold_onebody_glob0_p)[ir][irot] = 9e9;
                    }
                }
            }
        }
        if( opt.bonus_to_native_scaffold_res != 0 ){
            std::cout << "adding to native scaffold res 1BE " << opt.bonus_to_native_scaffold_res << std::endl;
            for( int ir = 0; ir < scaffold_onebody_glob0_p->size(); ++ir ){
                for( int irot = 0; irot < rot_index_p->size(); ++irot ){
                    if( rot_index_p->resname(irot) == scaffold_sequence_glob0_p->at(ir) ){
                        (*scaffold_onebody_glob0_p)[ir][irot] += opt.bonus_to_native_scaffold_res;
                    }
                }
            }
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

    // setup scaffold_twobody_p and local_twobody_p
    void
    setup_twobody_tables(  
        shared_ptr< RotamerIndex > rot_index_p,
        RifDockOpt const & opt ,
        MakeTwobodyOpts const & make2bopts,
        ::devel::scheme::RotamerRFTablesManager & rotrf_table_manager) {

        scaffold_twobody_p = make_shared<TBT>( scaffold_centered_p->size(), rot_index_p->size()  );

        std::cout << "rifdock: get_twobody_tables" << std::endl;
        std::string cachefile2b = "__2BE_" + scafftag + "_reshash" + scaff_res_hashstr + ".bin.gz";
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


        local_twobody_p = scaffold_twobody_p->create_subtable( scaffuseres, *scaffold_onebody_glob0_p, make2bopts.onebody_threshold );
    }


};


typedef shared_ptr<ScaffoldDataCache> ScaffoldDataCacheOP;
typedef shared_ptr<ScaffoldDataCache const > ScaffoldDataCacheCOP;

}}



#endif