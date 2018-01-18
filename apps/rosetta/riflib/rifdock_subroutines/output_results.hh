// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_output_results_hh
#define INCLUDED_riflib_rifdock_subroutines_output_results_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/rifdock_subroutines/meta.hh>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


template<class DirectorBase, class ScaffoldProvider>
struct DumpRifResultsData {
    RifDockOpt & opt;
    std::vector<float> & RESLS;
    DirectorBase & director;
    devel::scheme::ScenePtr & scene_minimal;
    std::vector<shared_ptr<devel::scheme::RifBase> > & rif_ptrs;
    devel::scheme::RotamerIndex & rot_index;
    core::pose::Pose & target;
    shared_ptr<ScaffoldProvider> & scaffold_provider;
    std::string const & resfileoutfile;
    std::string const & allrifrotsoutfile;


};



template<class DirectorBase, class ScaffoldProvider>
void
dump_search_point(
    DumpRifResultsData<DirectorBase, ScaffoldProvider> & d, 
    _SearchPoint<DirectorBase> const & search_point, 
    std::string const & pdboutfile, 
    int iresl,
    bool quiet) {

    typedef _RifDockResult<DirectorBase> RifDockResult;

    RifDockResult result;
    result.scene_index = search_point.index;

    dump_rif_result( d, result, pdboutfile, iresl, quiet );
}





template<class DirectorBase, class ScaffoldProvider>
void
dump_rif_result(
    DumpRifResultsData<DirectorBase, ScaffoldProvider> & d, 
    _RifDockResult<DirectorBase> const & selected_result, 
    std::string const & pdboutfile, 
    int iresl,
    bool quiet) {

    using namespace devel::scheme;
    using std::cout;
    using std::endl;


    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;

    ScaffoldIndex si = ::scheme::kinematics::bigindex_scaffold_index(selected_result.scene_index);
    ScaffoldDataCacheOP sdc = d.scaffold_provider->get_data_cache_slow( si );
    std::vector<int> const & scaffres_g2l = *(sdc->scaffres_g2l_p);
    std::vector<int> const & scaffres_l2g = *(sdc->scaffres_l2g_p);
    std::vector<std::vector<float> > const & scaffold_onebody_glob0 = *(sdc->scaffold_onebody_glob0_p);
    uint64_t const scaffold_size = scaffres_g2l.size();


    core::pose::Pose pose_from_rif;

    if ( d.opt.full_scaffold_output ) {        sdc->setup_both_full_pose( d.target ); pose_from_rif = *(sdc->mpc_both_full_pose.get_pose());
    } else if( d.opt.output_scaffold_only ) {                                         pose_from_rif = *(sdc->scaffold_centered_p);
    } else if( d.opt.output_full_scaffold_only ) {                                    pose_from_rif = *(sdc->scaffold_full_centered_p);
    } else {                                        sdc->setup_both_pose( d.target ); pose_from_rif = *(sdc->mpc_both_pose.get_pose());
    }


    d.director->set_scene( selected_result.scene_index, iresl, *d.scene_minimal );

    EigenXform xposition1 = d.scene_minimal->position(1);
    EigenXform xalignout = EigenXform::Identity();
    if( d.opt.align_to_scaffold ){
        xalignout = xposition1.inverse();
    }

    xform_pose( pose_from_rif, eigen2xyz(xalignout)           ,        scaffold_size+1, pose_from_rif.size());
    xform_pose( pose_from_rif, eigen2xyz(xalignout*xposition1),                      1,        scaffold_size );


    std::vector< std::pair< int, std::string > > brians_infolabels;

    std::ostringstream packout, allout;
    std::map< int, std::string > pikaa;
    int chain_no = pose_from_rif.num_chains();   
    int res_num = pose_from_rif.size() + 1;
    const std::string chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    for( int i_actor = 0; i_actor < d.scene_minimal->template num_actors<BBActor>(1); ++i_actor ){
        BBActor bba = d.scene_minimal->template get_actor<BBActor>(1,i_actor);
        int const ires = scaffres_l2g.at( bba.index_ );


        {
            std::vector< std::pair< float, int > > rotscores;
            d.rif_ptrs[iresl]->get_rotamers_for_xform( bba.position(), rotscores );
            typedef std::pair<float,int> PairFI;
            BOOST_FOREACH( PairFI const & p, rotscores ){
                int const irot = p.second;
                float const sc = p.first + scaffold_onebody_glob0.at( ires ).at( irot );
                if( sc < 0 ){
                    if (d.opt.all_rif_rots_as_models_not_chains) {
                        allout << "MODEL" << endl;
                    }
                    BOOST_FOREACH( SchemeAtom a, d.rot_index.rotamers_.at( irot ).atoms_ ){
                        a.set_position( xalignout * bba.position() * a.position() ); // is copy
                        a.nonconst_data().resnum = d.opt.all_rif_rots_as_models_not_chains ? ires : res_num;
                        a.nonconst_data().chain = d.opt.all_rif_rots_as_models_not_chains ? 'A' : chains.at( chain_no % 52 );
                        ::scheme::actor::write_pdb( allout, a, nullptr );
                    }
                    if (d.opt.all_rif_rots_as_models_not_chains) {
                        allout << "ENDMDL" << endl;
                    } else {
                        allout << "TER" << endl;
                        res_num++;
                        chain_no++;
                    }
                    char oneletter = d.rot_index.oneletter(irot);
                    if( std::find( pikaa[ires+1].begin(), pikaa[ires+1].end(), oneletter ) == pikaa[ires+1].end() ){
                        pikaa[ires+1] += oneletter;
                    }


                    // Brian
                    std::pair< int, int > sat1_sat2 = d.rif_ptrs.back()->get_sat1_sat2(bba.position(), irot);

                    if ( ! quiet ) {
                        std::cout << "Brian: " << oneletter << " " << sat1_sat2.first << " " << sat1_sat2.second << " sc: " << sc;
                        std::cout << " ires: " << ires << " irot: " << irot << " seqpos: " << ires+1;
                        std::cout << boost::str(boost::format(" rif score: %.2f 1-body: %.2f")%p.first%scaffold_onebody_glob0.at( ires ).at( irot )) << std::endl;
                    }

                    std::pair< int, std::string > brian_pair;
                    brian_pair.first = ires + 1;
                    brian_pair.second = "HOT_IN:" + str(sat1_sat2.first);
                    brians_infolabels.push_back(brian_pair);

                }

            }
        }

        // int packed_rot = -1;
        // for( int ipr = 0; ipr < selected_result.numrots(); ++ipr ){
        //     // std::cout << "checking rots " << sp.rotamers()[ipr].first << " " << d.scaffres_g2l[ires] << std::endl;
        //     if( selected_result.rotamers().at(ipr).first == scaffres_g2l.at( ires ) ){
        //         packed_rot = selected_result.rotamers().at(ipr).second;
        //     }
        // }
        // if( packed_rot >= 0 ){
        //     // packout << "MODEL" << endl;
        //     BOOST_FOREACH( SchemeAtom a, d.rot_index.rotamers_.at( packed_rot ).atoms_ ){
        //         a.set_position( xalignout * bba.position() * a.position() ); // is copy
        //         a.nonconst_data().resnum = ires;
        //         ::scheme::actor::write_pdb( packout, a, nullptr );
        //     }
        //     packout << "TER" << endl;
        // }

    }




    // place the rotamers
    core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
    std::ostringstream resfile, expdb;
    resfile << "ALLAA" << std::endl;
    resfile << "start" << std::endl;
    expdb << "rif_residues ";

    for( int ipr = 0; ipr < selected_result.numrots(); ++ipr ){
        int ires = scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
        int irot =                  selected_result.rotamers().at(ipr).second;
        core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(d.rot_index.resname(irot)) );
        pose_from_rif.replace_residue( ires+1, *newrsd, true );
        resfile << ires+1 << " A NATRO" << std::endl;
        expdb << ires+1 << (ipr+1<selected_result.numrots()?",":""); // skip comma on last one
        for( int ichi = 0; ichi < d.rot_index.nchi(irot); ++ichi ){
            pose_from_rif.set_chi( ichi+1, ires+1, d.rot_index.chi( irot, ichi ) );
        }
    }

    bool using_rosetta_model = selected_result.pose_ != nullptr;
    core::pose::Pose & pose_to_dump( *(selected_result.pose_ ? selected_result.pose_.get() : &pose_from_rif) );
    utility::io::ozstream out1( pdboutfile );
    // d.scene_full->set_position( 1, xalignout * xposition1 );
    // write_pdb( out1, dynamic_cast<Scene&>(*d.scene_full), d.rot_index.chem_index_ );
    if( !using_rosetta_model ){
        if( d.opt.pdb_info_pikaa ){
            for( auto p : pikaa ){
                std::sort( p.second.begin(), p.second.end() );
                pose_to_dump.pdb_info()->add_reslabel(p.first, "PIKAA" );
                pose_to_dump.pdb_info()->add_reslabel(p.first, p.second );
            }
        } else {
            for( auto p : pikaa ){
                pose_to_dump.pdb_info()->add_reslabel(p.first, "RIFRES" );
            }
        }
    }

    for ( auto p : brians_infolabels ) {
        pose_to_dump.pdb_info()->add_reslabel(p.first, p.second);
    }

    // if( selected_result.pose_ ){
    //  for( auto p : pikaa ){
    //      std::cout << "residue " << p.first << " " << selected_result.pose_->residue(p.first).name() << " fa_rep: "
    //                << selected_result.pose_->energies().residue_total_energies(p.first)[core::scoring::fa_rep] << std::endl;
    //  }
    // }

    out1 << expdb.str() << std::endl;
    pose_to_dump.dump_pdb(out1);
    if ( d.opt.dump_all_rif_rots_into_output ) {
        if ( ! d.opt.all_rif_rots_as_models_not_chains ) out1 << "TER" << endl;
        out1 << allout.str();
    }
    out1.close();


    // utility::io::ozstream outtmp( pdboutfile + ".orig.pdb" );
    // pose_from_rif.dump_pdb(outtmp);
    // outtmp.close();



    if( d.opt.dump_resfile ){
        utility::io::ozstream out1res( d.resfileoutfile );
        out1res << resfile.str();
        out1res.close();
    }

    if( d.opt.dump_all_rif_rots ){
        // utility_exit_with_message("this is not currently implemented, ask Will");
        utility::io::ozstream out2( d.allrifrotsoutfile );
        out2 << allout.str();
        out2.close();
    }

    // utility::io::ozstream out4(d.scafftag+"_pack_rot_"+devel::scheme::str(i_selected_result,9)+".pdb");
    // out4 << packout.str();
    // out4.close();

} // end crappy pdb io







template<class DirectorBase,class ScaffoldProvider>
struct OutputResultsData {
    RifDockOpt & opt;
    std::vector<float> & RESLS;
    DirectorBase & director;
    std::vector< _RifDockResult<DirectorBase> > & selected_results;
    int64_t & npack;
    utility::io::ozstream & dokout;
    devel::scheme::ScenePtr & scene_minimal;
    std::vector<shared_ptr<devel::scheme::RifBase> > & rif_ptrs;
    devel::scheme::RotamerIndex & rot_index;
    core::pose::Pose & target;
    shared_ptr<ScaffoldProvider> scaffold_provider;

};



template<class DirectorBase,class ScaffoldProvider>
void
output_results(
    OutputResultsData<DirectorBase,ScaffoldProvider> & d) {


    using namespace core::scoring;
        using std::cout;
        using std::endl;
        using namespace devel::scheme;
        typedef numeric::xyzVector<core::Real> Vec;
        typedef numeric::xyzMatrix<core::Real> Mat;
        // typedef numeric::xyzTransform<core::Real> Xform;
        using ObjexxFCL::format::F;
        using ObjexxFCL::format::I;
        using devel::scheme::print_header;
        using ::devel::scheme::RotamerIndex;

    typedef ::scheme::util::SimpleArray<3,float> F3;
    typedef ::scheme::util::SimpleArray<3,int> I3;

    typedef _RifDockResult<DirectorBase> RifDockResult;
    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;

    if( d.opt.align_to_scaffold ) std::cout << "ALIGN TO SCAFFOLD" << std::endl;
    else                        std::cout << "ALIGN TO TARGET"   << std::endl;
    for( int i_selected_result = 0; i_selected_result < d.selected_results.size(); ++i_selected_result ){
        RifDockResult const & selected_result = d.selected_results.at( i_selected_result );

// Brian Injection
        ScaffoldIndex si = ::scheme::kinematics::bigindex_scaffold_index(d.selected_results[i_selected_result].scene_index);
        ScaffoldDataCacheOP sdc = d.scaffold_provider->get_data_cache_slow( si );


        std::string const & scafftag = sdc->scafftag;

/////
        std::string pdboutfile = d.opt.outdir + "/" + scafftag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
        if( d.opt.output_tag.size() ){
            pdboutfile = d.opt.outdir + "/" + scafftag+"_" + d.opt.output_tag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
        }

        std::string resfileoutfile = d.opt.outdir + "/" + scafftag+"_"+devel::scheme::str(i_selected_result,9)+".resfile";
        std::string allrifrotsoutfile = d.opt.outdir + "/" + scafftag+"_allrifrots_"+devel::scheme::str(i_selected_result,9)+".pdb.gz";

        std::ostringstream oss;
        oss << "rif score: " << I(4,i_selected_result)
            << " rank "       << I(9,selected_result.isamp)
            << " dist0:    "  << F(7,2,selected_result.dist0)
            << " packscore: " << F(7,3,selected_result.packscore)
            // << " score: "     << F(7,3,selected_result.nopackscore)
            // << " rif: "       << F(7,3,selected_result.rifscore)
            << " steric: "    << F(7,3,selected_result.stericscore)
            << " cluster: "   << I(7,selected_result.cluster_score)
            << " rifrank: "   << I(7,selected_result.prepack_rank) << " " << F(7,5,(float)selected_result.prepack_rank/(float)d.npack)
            << " " << pdboutfile
            << std::endl;
        std::cout << oss.str();
        d.dokout << oss.str(); d.dokout.flush();


        {
            DumpRifResultsData<DirectorBase, ScaffoldProvider> data {
                d.opt,
                d.RESLS,
                d.director,
                d.scene_minimal,
                d.rif_ptrs,
                d.rot_index,
                d.target,
                d.scaffold_provider,
                resfileoutfile,
                allrifrotsoutfile,
            };

            dump_rif_result(data, selected_result, pdboutfile, d.RESLS.size()-1, false);
        }

         // crappy pdb io

    }



}









#endif