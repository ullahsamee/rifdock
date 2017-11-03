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

template<class DirectorBase>
struct OutputResultsData {
    RifDockOpt & opt;
    std::vector<float> & RESLS;
    DirectorBase & director;
    std::vector< _RifDockResult<DirectorBase> > & selected_results;
    std::string & scafftag;
    int64_t & npack;
    utility::io::ozstream & dokout;
    devel::scheme::ScenePtr & scene_full;
    devel::scheme::ScenePtr & scene_minimal;
    std::vector<int> & scaffres_g2l;
    std::vector<int> & scaffres_l2g;
    std::vector<shared_ptr<devel::scheme::RifBase> > & rif_ptrs;
    std::vector<std::vector<float> > & scaffold_onebody_glob0;
    devel::scheme::RotamerIndex & rot_index;
    core::pose::Pose & scaffold;
    core::pose::Pose & both_pose;
    core::pose::Pose & both_full_pose;                
    core::pose::Pose & scaffold_only_pose;            
    core::pose::Pose & scaffold_only_full_pose;
};



template<class DirectorBase>
void
output_results(
    OutputResultsData<DirectorBase> & d) {


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

    if( d.opt.align_to_scaffold ) std::cout << "ALIGN TO SCAFFOLD" << std::endl;
    else                        std::cout << "ALIGN TO TARGET"   << std::endl;
    for( int i_selected_result = 0; i_selected_result < d.selected_results.size(); ++i_selected_result ){
        RifDockResult const & selected_result = d.selected_results.at( i_selected_result );

        std::string pdboutfile = d.opt.outdir + "/" + d.scafftag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
        if( d.opt.output_tag.size() ){
            pdboutfile = d.opt.outdir + "/" + d.scafftag+"_" + d.opt.output_tag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
        }

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

        std::vector< std::pair< int, std::string > > brians_infolabels;
         // crappy pdb io
        {

            d.director->set_scene( selected_result.scene_index, d.RESLS.size()-1, *d.scene_full    );
            d.director->set_scene( selected_result.scene_index, d.RESLS.size()-1, *d.scene_minimal );

            EigenXform xposition1 = d.scene_full->position(1);
            EigenXform xalignout = EigenXform::Identity();
            if( d.opt.align_to_scaffold ){
                xalignout = xposition1.inverse();
            }

            std::ostringstream packout, allout;
            std::map< int, std::string > pikaa;
            for( int i_actor = 0; i_actor < d.scene_minimal->template num_actors<BBActor>(1); ++i_actor ){
                BBActor bba = d.scene_minimal->template get_actor<BBActor>(1,i_actor);
                int const ires = d.scaffres_l2g.at( bba.index_ );

                // if( d.opt.dump_all_rif_rots )
                {
                    std::vector< std::pair< float, int > > rotscores;
                    d.rif_ptrs.back()->get_rotamers_for_xform( bba.position(), rotscores );
                    typedef std::pair<float,int> PairFI;
                    BOOST_FOREACH( PairFI const & p, rotscores ){
                        int const irot = p.second;
                        float const sc = p.first + d.scaffold_onebody_glob0.at( ires ).at( irot );
                        if( sc < 0 ){
                            allout << "MODEL" << endl;
                            BOOST_FOREACH( SchemeAtom a, d.rot_index.rotamers_.at( irot ).atoms_ ){
                                a.set_position( xalignout * bba.position() * a.position() ); // is copy
                                a.nonconst_data().resnum = ires;
                                ::scheme::actor::write_pdb( allout, a, nullptr );
                            }
                            allout << "ENDMDL" << endl;
                            char oneletter = d.rot_index.oneletter(irot);
                            if( std::find( pikaa[ires+1].begin(), pikaa[ires+1].end(), oneletter ) == pikaa[ires+1].end() ){
                                pikaa[ires+1] += oneletter;
                            }


                            // Brian
                            std::pair< int, int > sat1_sat2 = d.rif_ptrs.back()->get_sat1_sat2(bba.position(), irot);

                            std::cout << "Brian: " << oneletter << " " << sat1_sat2.first << " " << sat1_sat2.second << " sc: " << sc;
                            std::cout << " ires: " << ires << " irot: " << irot << std::endl;

                            std::pair< int, std::string > brian_pair;
                            brian_pair.first = ires + 1;
                            brian_pair.second = "HOT_IN:" + str(sat1_sat2.first);
                            brians_infolabels.push_back(brian_pair);

                        }

                    }
                }

                int packed_rot = -1;
                for( int ipr = 0; ipr < selected_result.numrots(); ++ipr ){
                    // std::cout << "checking rots " << sp.rotamers()[ipr].first << " " << d.scaffres_g2l[ires] << std::endl;
                    if( selected_result.rotamers().at(ipr).first == d.scaffres_g2l.at( ires ) ){
                        packed_rot = selected_result.rotamers().at(ipr).second;
                    }
                }
                if( packed_rot >= 0 ){
                    // packout << "MODEL" << endl;
                    BOOST_FOREACH( SchemeAtom a, d.rot_index.rotamers_.at( packed_rot ).atoms_ ){
                        a.set_position( xalignout * bba.position() * a.position() ); // is copy
                        a.nonconst_data().resnum = ires;
                        ::scheme::actor::write_pdb( packout, a, nullptr );
                    }
                    packout << "TER" << endl;
                }

            }

            core::pose::Pose pose_from_rif;
            if     ( d.opt.full_scaffold_output ) pose_from_rif = d.both_full_pose;
            else if( d.opt.output_scaffold_only ) pose_from_rif = d.scaffold_only_pose;
            else if( d.opt.output_full_scaffold_only ) pose_from_rif = d.scaffold_only_full_pose;
            else                                pose_from_rif = d.both_pose;
            xform_pose( pose_from_rif, eigen2xyz(xalignout)           , d.scaffold.size()+1, pose_from_rif.size() );
            xform_pose( pose_from_rif, eigen2xyz(xalignout*xposition1),                      1,     d.scaffold.size() );

            // place the rotamers
            core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
            std::ostringstream resfile, expdb;
            resfile << "ALLAA" << std::endl;
            resfile << "start" << std::endl;
            expdb << "rif_residues ";

            for( int ipr = 0; ipr < selected_result.numrots(); ++ipr ){
                int ires = d.scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
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
            out1.close();


            // utility::io::ozstream outtmp( pdboutfile + ".orig.pdb" );
            // pose_from_rif.dump_pdb(outtmp);
            // outtmp.close();



            if( d.opt.dump_resfile ){
                utility::io::ozstream out1res( d.opt.outdir + "/" + d.scafftag+"_"+devel::scheme::str(i_selected_result,9)+".resfile");
                out1res << resfile.str();
                out1res.close();
            }

            if( d.opt.dump_all_rif_rots ){
                // utility_exit_with_message("this is not currently implemented, ask Will");
                utility::io::ozstream out2( d.opt.outdir + "/" + d.scafftag+"_allrifrots_"+devel::scheme::str(i_selected_result,9)+".pdb.gz");
                out2 << allout.str();
                out2.close();
            }

            // utility::io::ozstream out4(d.scafftag+"_pack_rot_"+devel::scheme::str(i_selected_result,9)+".pdb");
            // out4 << packout.str();
            // out4.close();

        } // end crappy pdb io

    }



}







#endif