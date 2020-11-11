// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/OutputResultsTasks.hh>

#include <riflib/types.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>
#include <riflib/rifdock_tasks/HackPackTasks.hh>
#include <riflib/ScoreRotamerVsTarget.hh>
#include <riflib/RifFactory.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/kinematics/RT.hh>

#include <sys/param.h>
#include <algorithm>
#include <string>
#include <vector>

#include <ObjexxFCL/format.hh>



namespace devel {
namespace scheme {

shared_ptr<std::vector<RifDockResult>>
OutputResultsTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> selected_results_p, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::vector<RifDockResult> & selected_results = *selected_results_p;

    using std::cout;
    using std::endl;
    using ObjexxFCL::format::F;
    using ObjexxFCL::format::I;


    std::cout << " selected_results.size(): " << selected_results.size() << std::endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    print_header( "timing info" ); //////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<<"total RIF     time: "<<KMGT(pd.time_rif)<<" fraction: "<<pd.time_rif/(pd.time_rif+pd.time_pck+pd.time_ros)<<std::endl;
    std::cout<<"total Pack    time: "<<KMGT(pd.time_pck)<<" fraction: "<<pd.time_pck/(pd.time_rif+pd.time_pck+pd.time_ros)<<std::endl;
    std::cout<<"total Rosetta time: "<<KMGT(pd.time_ros)<<" fraction: "<<pd.time_ros/(pd.time_rif+pd.time_pck+pd.time_ros)<<std::endl;   


    if ( rdd.unsat_manager && rdd.opt.report_common_unsats ) {

        std::vector<shared_ptr<UnsatManager>> & unsatperthread = rdd.rif_factory->get_unsatperthread( rdd.objectives.back() );

        for ( shared_ptr<UnsatManager> const & man : unsatperthread ) {
            rdd.unsat_manager->sum_unsat_counts( *man );
        }

        rdd.unsat_manager->print_unsat_counts();
    }






    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    print_header( "output results" ); //////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    if( rdd.opt.align_to_scaffold ) std::cout << "ALIGN TO SCAFFOLD" << std::endl;
    else                        std::cout << "ALIGN TO TARGET"   << std::endl;
    utility::io::ozstream out_silent_stream;
    if ( rdd.opt.outputsilent || rdd.opt.outputlite ) {
        ScaffoldDataCacheOP example_data_cache = rdd.scaffold_provider->get_data_cache_slow( ScaffoldIndex() );
        out_silent_stream.open_append( rdd.opt.outdir + "/" + example_data_cache->scafftag + ".silent" );
    }
    for( int i_selected_result = 0; i_selected_result < selected_results.size(); ++i_selected_result ){
        RifDockResult const & selected_result = selected_results.at( i_selected_result );

// Brian Injection
        ScaffoldIndex si = selected_result.index.scaffold_index;
        ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow( si );

        std::string seeding_tag = "";
        if ( selected_result.index.seeding_index < pd.seeding_tags.size() ) seeding_tag = pd.seeding_tags[ selected_result.index.seeding_index ];

        std::string const & use_scafftag = sdc->scafftag + seeding_tag;

/////

        rdd.director->set_scene( selected_result.index, director_resl_, *rdd.scene_minimal );

        std::vector<float> unsat_scores;
        int unsats = -1;
        int buried = -1;
        if ( rdd.unsat_manager ) {
            std::vector<EigenXform> bb_positions;
            for ( int i_actor = 0; i_actor < rdd.scene_minimal->template num_actors<BBActor>(1); i_actor++ ) {
                bb_positions.push_back( rdd.scene_minimal->template get_actor<BBActor>(1,i_actor).position() );
            }
            std::vector<float> burial = rdd.burial_manager->get_burial_weights( rdd.scene_minimal->position(1), sdc->burial_grid );
            unsat_scores =  rdd.unsat_manager->get_buried_unsats( burial, selected_result.rotamers(), bb_positions, rdd.rot_tgt_scorer );

            buried = 0;
            for ( float this_burial : burial ) if ( this_burial > 0 ) buried++;
            unsats = 0;
            for ( float this_unsat_score : unsat_scores ) if ( this_unsat_score > 0 ) unsats++;
        }

        std::stringstream extra_output;

        int hydrophobic_residue_contacts;
        float hydrophobic_ddg = 0;
        if ( rdd.hydrophobic_manager ) {
            std::vector<int> hydrophobic_counts, seqposs, per_irot_counts;
            std::vector<std::pair<intRot, EigenXform>> irot_and_bbpos;
            selected_result.rotamers();
            for( int i = 0; i < selected_result.rotamers_->size(); ++i ){
                BBActor const & bb = rdd.scene_minimal->template get_actor<BBActor>( 1, selected_result.rotamers_->at(i).first );
                int seqpos = sdc->scaffres_l2g_p->at( bb.index_ ) + 1;
                int irot = selected_result.rotamers_->at(i).second;

                irot_and_bbpos.emplace_back( irot, bb.position() );
                seqposs.push_back( seqpos );
            }
            bool pass_better_than = true, pass_cation_pi = true;
            hydrophobic_residue_contacts = rdd.hydrophobic_manager->find_hydrophobic_residue_contacts( irot_and_bbpos, hydrophobic_counts, hydrophobic_ddg,
                                                                            per_irot_counts, pass_better_than, pass_cation_pi, rdd.rot_tgt_scorer );

            rdd.hydrophobic_manager->print_hydrophobic_counts( rdd.target, hydrophobic_counts, irot_and_bbpos, seqposs, per_irot_counts, 
                                                                            sdc->scaffres_g2l_p->size(), extra_output );
        }


        std::string pdboutfile = rdd.opt.outdir + "/" + use_scafftag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
        if( rdd.opt.output_tag.size() ){
            pdboutfile = rdd.opt.outdir + "/" + use_scafftag+"_" + rdd.opt.output_tag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
        }

        std::string resfileoutfile = rdd.opt.outdir + "/" + use_scafftag+"_"+devel::scheme::str(i_selected_result,9)+".resfile";
        std::string allrifrotsoutfile = rdd.opt.outdir + "/" + use_scafftag+"_allrifrots_"+devel::scheme::str(i_selected_result,9)+".pdb.gz";

        std::ostringstream oss;
        oss << "rif score: " << I(4,i_selected_result)
            << " rank "       << I(9,selected_result.isamp)
            << " dist0:    "  << F(7,2,selected_result.dist0)
            << " packscore: " << F(7,3,selected_result.score)
            << " score: "     << F(7,3,selected_result.nopackscore)
            // << " rif: "       << F(7,3,selected_result.rifscore)
            << " steric: "    << F(7,3,selected_result.stericscore);
        if (rdd.opt.scaff_bb_hbond_weight > 0) {
        oss << " bb-hbond: "  << F(7,3,selected_result.scaff_bb_hbond);
        }
        oss << " cluster: "   << I(7,selected_result.cluster_score)
            << " rifrank: "   << I(7,selected_result.prepack_rank) << " " << F(7,5,(float)selected_result.prepack_rank/(float)pd.npack);
        if ( rdd.unsat_manager ) {
        oss << " buried:" << I(4,buried);
        oss << " unsats:" << I(4, unsats);
        }
        if ( rdd.opt.need_to_calculate_sasa ) {
        oss << " sasa:" << I(5, selected_result.sasa);
        }
        if ( rdd.hydrophobic_manager ) {
        oss << " hyd-cont:" << I(3, hydrophobic_residue_contacts);
        oss << " hyd-ddg: " << F(7,3, hydrophobic_ddg);
        }
        oss << " " << pdboutfile
            << std::endl;
        std::cout << oss.str();
        rdd.dokout << oss.str(); rdd.dokout.flush();

        dump_rif_result_(rdd, selected_result, pdboutfile, director_resl_, rif_resl_, out_silent_stream, false, resfileoutfile, allrifrotsoutfile, unsat_scores);

        std::cout << extra_output.str() << std::flush;

    }

    return selected_results_p;

}


void
dump_rif_result_(
    RifDockData & rdd,
    RifDockResult const & selected_result, 
    std::string const & pdboutfile, 
    int director_resl,
    int rif_resl,
    utility::io::ozstream & out_silent_stream,
    bool quiet /* = true */,
    std::string const & resfileoutfile /* = "" */,
    std::string const & allrifrotsoutfile, /* = "" */
    std::vector<float> const & unsat_scores /* = std::vector<float>() */
    ) {

    using ObjexxFCL::format::F;
    using ObjexxFCL::format::I;

    using namespace devel::scheme;
    using std::cout;
    using std::endl;


    ScaffoldIndex si = selected_result.index.scaffold_index;
    ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow( si );
    std::vector<int> const & scaffres_g2l = *(sdc->scaffres_g2l_p);
    std::vector<int> const & scaffres_l2g = *(sdc->scaffres_l2g_p);
    std::vector<std::vector<float> > const & scaffold_onebody_glob0 = *(sdc->scaffold_onebody_glob0_p);
    uint64_t const scaffold_size = scaffres_g2l.size();


    core::pose::Pose pose_from_rif;

    if ( rdd.opt.output_full_scaffold ) {        sdc->setup_both_full_pose( rdd.target ); pose_from_rif = *(sdc->mpc_both_full_pose.get_pose());
    } else if( rdd.opt.output_scaffold_only ) {                                         pose_from_rif = *(sdc->scaffold_centered_p);
    } else if( rdd.opt.output_full_scaffold_only ) {                                    pose_from_rif = *(sdc->scaffold_full_centered_p);
    } else {                                        sdc->setup_both_pose( rdd.target ); pose_from_rif = *(sdc->mpc_both_pose.get_pose());
    }


    rdd.director->set_scene( selected_result.index, director_resl, *rdd.scene_minimal );

    EigenXform xposition1 = rdd.scene_minimal->position(1);
    EigenXform xalignout = EigenXform::Identity();
    if( rdd.opt.align_to_scaffold ){
        xalignout = xposition1.inverse();
    }

    xform_pose( pose_from_rif, eigen2xyz(xalignout)           ,        scaffold_size+1, pose_from_rif.size());
    xform_pose( pose_from_rif, eigen2xyz(xalignout*xposition1),                      1,        scaffold_size );


    std::vector< std::pair< int, std::string > > brians_infolabels;

    std::ostringstream packout, allout;
    // TYU change to vector of strings instead of string
    std::map< int, std::vector<std::string> > pikaa;
    std::map< int, std::map<std::string,float> > pssm;
    int chain_no = pose_from_rif.num_chains();   
    int res_num = pose_from_rif.size() + 1;
    const std::string chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    for( int i_actor = 0; i_actor < rdd.scene_minimal->template num_actors<BBActor>(1); ++i_actor ){
        BBActor bba = rdd.scene_minimal->template get_actor<BBActor>(1,i_actor);
        int const ires = scaffres_l2g.at( bba.index_ );


        {
            std::vector< std::pair< float, int > > rotscores;
            rdd.rif_ptrs[rif_resl]->get_rotamers_for_xform( bba.position(), rotscores );
            typedef std::pair<float,int> PairFI;
            BOOST_FOREACH( PairFI const & p, rotscores ){
                int const irot = p.second;
                float const onebody = scaffold_onebody_glob0.at( ires ).at( irot );
                float const sc = p.first + onebody;
                float const rescore = rdd.rot_tgt_scorer.score_rotamer_v_target( irot, bba.position(), 10.0, 4 );
                if( sc < 0 || rescore + onebody < 0  || p.first + onebody < 0){
                    if ( ! rdd.opt.rif_rots_as_chains) {
                        allout << "MODEL" << endl;
                    }
                    BOOST_FOREACH( SchemeAtom a, rdd.rot_index_p->rotamers_.at( irot ).atoms_ ){
                        a.set_position( xalignout * bba.position() * a.position() ); // is copy
                        a.nonconst_data().resnum = rdd.opt.rif_rots_as_chains ? res_num : ires;
                        a.nonconst_data().chain = rdd.opt.rif_rots_as_chains ? chains.at( chain_no % 52 ) : 'A';
                        ::scheme::actor::write_pdb( allout, a, nullptr );
                    }
                    if (! rdd.opt.rif_rots_as_chains) {
                        allout << "ENDMDL" << endl;
                    } else {
                        allout << "TER" << endl;
                        res_num++;
                        chain_no++;
                    }
                    // TYU change to std::string for expanded oneletter map
                    std::string oneletter = rdd.rot_index_p->oneletter(irot);
                    if( std::find( pikaa[ires+1].begin(), pikaa[ires+1].end(), oneletter ) == pikaa[ires+1].end() ){
                        pikaa[ires+1].push_back(oneletter);
                    }

                    // Brian
                    std::pair< int, int > sat1_sat2 = rdd.rif_ptrs.back()->get_sat1_sat2(bba.position(), irot);

                    // This should probably take into account more than just the sat scores
                    //  This is supposed to capture the score from RotamerFactory.cc operator()
                    float effective_score = p.first;
                    if ( rdd.opt.sat_score_bonus.size() > 0 ) {
                        if ( sat1_sat2.first != -1 ) {
                            float bonus_or_override = rdd.opt.sat_score_bonus.at( sat1_sat2.first );
                            if ( rdd.opt.sat_score_override.at( sat1_sat2.first ) ) {
                                effective_score = bonus_or_override;
                            } else {
                                effective_score += bonus_or_override;
                            }
                        }
                    }
                    float pssm_score = effective_score;
                    if ( pssm[ires+1].count( oneletter ) ) {
                        pssm_score = std::min<float>( pssm[ires+1].at( oneletter ), pssm_score );
                    }
                    pssm[ires+1][oneletter] = pssm_score;



                    bool rot_was_placed = false;
                    for ( std::pair<intRot,intRot> const & placed_rot : selected_result.rotamers() ) {
                        if ( placed_rot.first == bba.index_ && placed_rot.second == irot ) {
                            rot_was_placed = true;
                            break;
                        }
                    }

                    float rotboltz = 0;
                    if ( sdc->rotboltz_data_p ) {
                        if ( sdc->rotboltz_data_p->at(ires).size() > 0 ) {
                            rotboltz = sdc->rotboltz_data_p->at(ires)[irot];
                        }
                    }

                    if ( ! quiet ) {

                        std::cout << ( rot_was_placed ? "*" : " " );
                        std::cout << "seqpos:" << I(3, ires+1);
                        std::cout << " " << oneletter;
                        std::cout << " score:" << F(7, 2, sc);
                        std::cout << " irot:" << I(3, irot);
                        std::cout << " 1-body:" << F(7, 2, onebody );
                        if ( sdc->rotboltz_data_p ) {
                            std::cout << " rotboltz:" << F(7, 2, rotboltz );
                        }
                        std::cout << " rif score:" << F(7, 2, p.first);
                        std::cout << " rif rescore:" << F(7, 2, rescore);
                        std::cout << " sats:" << I(3, sat1_sat2.first) << " " << I(3, sat1_sat2.second) << " ";
                        std::cout << std::endl;



                    }

                    if (sat1_sat2.first > -1) {
                        std::pair< int, std::string > brian_pair;
                        brian_pair.first = ires + 1;
                        brian_pair.second = "HOT_IN:" + str(sat1_sat2.first);
                        brians_infolabels.push_back(brian_pair);
                    }

                }

            }
        }

    }

    if ( unsat_scores.size() > 0 ) {
        rdd.unsat_manager->print_buried_unsats( unsat_scores, scaffold_size );
    }


    // // TEMP debug:
    // for (auto i: scaffold_phi_psi) {
    //     std::cout << std::get<0>(i) << " " << std::get<1>(i) << std::endl;
    // }
    // for (auto i: scaffold_d_pos) {
    //     std::cout << i << " ";
    // }

    // Actually place the rotamers on the pose
    core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
    std::ostringstream resfile, expdb;
    resfile << "ALLAA" << std::endl;
    resfile << "start" << std::endl;
    expdb << "rif_residues ";

    if ( selected_result.rotamers_ ) {
        sanity_check_hackpack( rdd, selected_result.index, selected_result.rotamers_, rdd.scene_pt.front(), director_resl, rif_resl);
    }

    std::vector<int> needs_RIFRES;
    for( int ipr = 0; ipr < selected_result.numrots(); ++ipr ){
        int ires = scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
        int irot =                  selected_result.rotamers().at(ipr).second;

        std::string myResName = rdd.rot_index_p->resname(irot);
        auto myIt = rdd.rot_index_p -> d_l_map_.find(myResName);

        core::conformation::ResidueOP newrsd;
        if (myIt != rdd.rot_index_p -> d_l_map_.end()){
            core::chemical::ResidueType const & rtype = rts.lock()->name_map( myIt -> second );
            newrsd = core::conformation::ResidueFactory::create_residue( rtype );
            core::pose::Pose pose;
            pose.append_residue_by_jump(*newrsd,1);
            core::chemical::ResidueTypeSetCOP pose_rts = pose.residue_type_set_for_pose();
            core::chemical::ResidueTypeCOP pose_rt = get_restype_for_pose(pose, myIt -> second);
            core::chemical::ResidueTypeCOP d_pose_rt = pose_rts -> get_d_equivalent(pose_rt);
            newrsd = core::conformation::ResidueFactory::create_residue( *d_pose_rt );
        } else {
            newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rdd.rot_index_p->resname(irot)) );
        }
        //core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rdd.rot_index_p->resname(irot)) );

        pose_from_rif.replace_residue( ires+1, *newrsd, true );
        resfile << ires+1 << " A NATRO" << std::endl;
        expdb << ires+1 << (ipr+1<selected_result.numrots()?",":""); // skip comma on last one
        for( int ichi = 0; ichi < rdd.rot_index_p->nchi(irot); ++ichi ){
            pose_from_rif.set_chi( ichi+1, ires+1, rdd.rot_index_p->chi( irot, ichi ) );
        }
        needs_RIFRES.push_back(ires+1);
    }

    // Add PDBInfo labels if they are applicable
    bool using_rosetta_model = (selected_result.pose_ != nullptr) && !rdd.opt.override_rosetta_pose;

    core::pose::PoseOP stored_pose = selected_result.pose_;
    if ( using_rosetta_model && ( rdd.opt.output_scaffold_only || rdd.opt.output_full_scaffold_only ) ) {
        stored_pose = stored_pose->split_by_chain().front();
    }

    core::pose::Pose & pose_to_dump( using_rosetta_model ? *stored_pose : pose_from_rif );
    if( !using_rosetta_model ){
        if( rdd.opt.pdb_info_pikaa ){
            for( auto p : pikaa ){
                std::sort( p.second.begin(), p.second.end() );
                pose_to_dump.pdb_info()->add_reslabel(p.first, "PIKAA" );
                // TYU create output string for reslabel
                std::string out_string;
                for (auto i : p.second) {
                    out_string += i;
                    out_string += ",";
                }
                pose_to_dump.pdb_info()->add_reslabel(p.first, out_string );
                //pose_to_dump.pdb_info()->add_reslabel(p.first, p.second );
            }
        } else {
            std::sort(needs_RIFRES.begin(), needs_RIFRES.end());
            for( int seq_pos : needs_RIFRES ){
                pose_to_dump.pdb_info()->add_reslabel(seq_pos, "RIFRES" );
            }
        }
    }
    if ( rdd.opt.pdb_info_pssm ) {
        // For now, only output stuff at positions where we actually placed something.
        for ( int seq_pos : needs_RIFRES ) {
            utility::vector1<std::string> parts;
            for ( auto name_score : pssm[seq_pos] ) {
                std::string this_string = boost::str(boost::format("%s:%.1f")%name_score.first%name_score.second);
                parts.push_back( this_string );
            }
            if ( parts.size() > 0 ) {
                std::string full_string = "PSSM:" + utility::join(parts, ",");
                pose_to_dump.pdb_info()->add_reslabel(seq_pos, full_string );
            }
        }
    }

    for ( auto p : brians_infolabels ) {
        pose_to_dump.pdb_info()->add_reslabel(p.first, p.second);
    }

    rdd.scaffold_provider->modify_pose_for_output(si, pose_to_dump);


    // Dump the main output
    if ( !rdd.opt.outputsilent && !rdd.opt.outputlite ) {
        utility::io::ozstream out1( pdboutfile );
        out1 << expdb.str() << std::endl;
        pose_to_dump.dump_pdb(out1);
        if ( rdd.opt.dump_all_rif_rots_into_output ) {
            if ( rdd.opt.rif_rots_as_chains ) out1 << "TER" << endl;
            out1 << allout.str();
        }
        out1.close();
    }
    // Dump a resfile
    if( rdd.opt.dump_resfile ){
        utility::io::ozstream out1res( resfileoutfile );
        out1res << resfile.str();
        out1res.close();
    }

    // Dump the rif rots
    if( rdd.opt.dump_all_rif_rots ){
        utility::io::ozstream out2( allrifrotsoutfile );
        out2 << allout.str();
        out2.close();
    }
    // Dump silent file
    if (rdd.opt.outputsilent) {
        // silly thing to take care of multiple chains in PDBInfo for silentstruct
        core::Size const ch1end(1);
        //set chain ID
        utility::vector1<char> chainID_vec(pose_to_dump.conformation().size()); 
        for (auto i = 1; i <= pose_to_dump.conformation().chain_end(ch1end); i++) chainID_vec[i] = 'A';
        for (auto i = pose_to_dump.conformation().chain_end(ch1end) + 1; i <= pose_to_dump.conformation().size(); i++) chainID_vec[i] = 'B';
        pose_to_dump.pdb_info() -> set_chains(chainID_vec);
        //set chain num
        std::vector<int> v(pose_to_dump.conformation().size()); 
        std::iota (std::begin(v), std::end(v), 1);
        pose_to_dump.pdb_info() -> set_numbering(v.begin(), v.end());
        //dump to silent
        std::string model_tag = pdb_name(pdboutfile);
        core::io::silent::SilentFileOptions sf_option;
        sf_option.read_from_global_options();
        core::io::silent::SilentFileData sfd("", false, false, "binary", sf_option);
        core::io::silent::SilentStructOP ss = sfd.create_SilentStructOP();
        ss->fill_struct( pose_to_dump, model_tag );
        sfd._write_silent_struct(*ss, out_silent_stream);
    }

    // Dump silent file of PDBLite structures
    if ( rdd.opt.outputlite ) {

        core::pose::Pose lite_pose;
        core::pose::make_pose_from_sequence( lite_pose, "A", core::chemical::FA_STANDARD, false ); // Auto-termini set to false

        core::pose::PDBInfoOP pdb_info = make_shared< core::pose::PDBInfo >( lite_pose );

        pdb_info->add_reslabel( 1, "SCAFFOLD_PDB:" + absolute_path( sdc->scaff_fname ) );
        pdb_info->add_reslabel( 1, "TARGET_PDB:" + absolute_path( rdd.opt.target_pdb ) );

        EigenXform binder_xform = EigenXform( xposition1 );
        binder_xform.translate( -1 * sdc->scaffold_center );

        std::ostringstream RT_binder_stream;
        RT_binder_stream << core::kinematics::RT( eigen2xyz(binder_xform).R , eigen2xyz(binder_xform).t );

        std::string RT_binder_str = "SCAFFOLD_XFORM:" + RT_binder_stream.str().substr( 3 ); // Remove initial RT from xform string

        std::replace( RT_binder_str.begin(), RT_binder_str.end(), ' ', '_');

        pdb_info->add_reslabel( 1, RT_binder_str );

        ////////////////////////////////////////////////////////////////////////////////////////////
        // TODO: implement mutation field here - NRB
        ////////////////////////////////////////////////////////////////////////////////////////////

        // Final setting of pdb_info
        lite_pose.pdb_info( pdb_info );

        std::string model_tag = pdb_name(pdboutfile);
        core::io::silent::SilentFileOptions sf_option;
        sf_option.read_from_global_options();
        core::io::silent::SilentFileData sfd("", false, false, "binary", sf_option);
        core::io::silent::SilentStructOP ss = sfd.create_SilentStructOP();
        ss->fill_struct( lite_pose, model_tag ); // only line changed from above, would be prettier to make this it's own func - NRB
        sfd._write_silent_struct(*ss, out_silent_stream);

    }

}

std::string
absolute_path( std::string relative_path ){

    if( relative_path[0] == '/' ) return relative_path;

    // You would think that it would be easier to get the cwd but it's not
    // It's easy with C++17 but that requires gcc7 which I don't think works with rifdock - NRB
    // https://stackoverflow.com/questions/2203159/is-there-a-c-equivalent-to-getcwd
    char temp[MAXPATHLEN];
    std::string cwd = getcwd(temp, sizeof(temp)) ? std::string( temp ) : std::string("");

    // The plan is to count instances of "../" in the relative path
    // Then that number of directories from the cwd will be removed
    std::string relative_path_canonical;
    core::Size back_dir_counter = 0;

    utility::vector1< std::string > rel_path_splits = utility::string_split( relative_path, '/' );
    for ( core::Size i = 1; i <= rel_path_splits.size(); ++i ){
        if( rel_path_splits[i] == ".." ) ++back_dir_counter;
        else relative_path_canonical += "/" + rel_path_splits[i];
    }

    std::string cwd_truncated;
    utility::vector1< std::string > cwd_splits = utility::string_split( cwd, '/' );
    for ( core::Size i = 1; i <= cwd_splits.size() - back_dir_counter; ++i ){
        if( cwd_splits[i].size() <= 0 ) continue;
        cwd_truncated += "/" + cwd_splits[i];
    }
    
    return cwd_truncated + relative_path_canonical;
}

void
dump_search_point_(
    RifDockData & rdd,
    SearchPoint const & search_point, 
    std::string const & pdboutfile, 
    int director_resl,
    int rif_resl,
    bool quiet) {

    RifDockResult result;
    result = search_point;
    utility::io::ozstream trash;
    dump_rif_result_( rdd, result, pdboutfile, director_resl, rif_resl, trash,quiet, "", "" );
    trash.close();
}


}}
