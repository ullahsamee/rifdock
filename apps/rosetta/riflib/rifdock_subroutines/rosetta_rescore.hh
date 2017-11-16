// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_rosetta_rescore_hh
#define INCLUDED_riflib_rifdock_subroutines_rosetta_rescore_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/rifdock_subroutines/meta.hh>

#include <unordered_map>

using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;

template<class DirectorBase, class ScaffoldProvider >
struct RosettaRescoreData {
    RifDockOpt & opt;
    std::vector<float> & RESLS;
    DirectorBase & director;
    // std::vector<int> & scaffres_l2g;///////
    devel::scheme::RotamerIndex & rot_index;
    // core::pose::Pose & scaffold;////
    // core::pose::Pose & both_pose;////
    // core::pose::Pose & both_full_pose;//////
    // utility::vector1<core::Size> & scaffold_res;///////
    core::pose::Pose & target;
    int64_t & total_search_effort;
    std::vector< _SearchPointWithRots<DirectorBase> > & packed_results;
    std::vector< devel::scheme::ScenePtr > & scene_pt;
    std::vector< devel::scheme::VoxelArrayPtr > & target_field_by_atype;
    double & time_ros;
    shared_ptr<ScaffoldProvider> scaffold_provider;
};



template<class DirectorBase, class ScaffoldProvider>
void
rosetta_rescore(
    RosettaRescoreData<DirectorBase,ScaffoldProvider> & d) {


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

    typedef _SearchPointWithRots<DirectorBase> SearchPointWithRots;
    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;


    std::chrono::time_point<std::chrono::high_resolution_clock> start_rosetta = std::chrono::high_resolution_clock::now();

    int n_score_calculations = 0;
    int do_min = 2;
    if( d.opt.rosetta_min_fraction == 0.0 ) do_min = 1;

    // std::vector<bool> is_scaffold_fixed_res(scaffold_size+1,true);  // this is a one-indexed lookup
    // for(int designable : d.scaffold_res){
    //     is_scaffold_fixed_res[designable] = false;
    // }

    for( int minimizing = 0; minimizing < do_min; ++minimizing ){

        std::chrono::duration<double> time_copy, time_score, time_min;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( minimizing ) print_header( "rosetta min and score" ); ////////////////////////////////////////////////////////////////////
        else         print_header( "rosetta score" ); ////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector<core::kinematics::MoveMapOP> movemap_pt(omp_max_threads());
        std::vector<protocols::simple_moves::MinMoverOP> minmover_pt(omp_max_threads());
        std::vector<core::scoring::ScoreFunctionOP> scorefunc_pt(omp_max_threads());
        // std::vector<core::pose::Pose> both_full_per_thread(omp_max_threads());
        // std::vector<core::pose::Pose> both_per_thread     (omp_max_threads());
        // std::vector<core::pose::Pose> target_pt           (omp_max_threads());
        std::vector<core::pose::Pose> work_pose_pt        (omp_max_threads());
        for( int i = 0; i < omp_max_threads(); ++i){
            // both_full_per_thread[i] = d.both_full_pose;
            // if( d.opt.replace_orig_scaffold_res ){
            //     both_per_thread[i] = d.both_full_pose;
            // } else {
            //     both_per_thread[i] = d.both_pose;
            // }
            scorefunc_pt[i] = core::scoring::ScoreFunctionFactory::create_score_function(d.opt.rosetta_soft_score);
            if( minimizing ){
                if( d.opt.rosetta_hard_min ){
                    scorefunc_pt[i] = core::scoring::ScoreFunctionFactory::create_score_function(d.opt.rosetta_hard_score);
                } else {
                    scorefunc_pt[i] = core::scoring::ScoreFunctionFactory::create_score_function(d.opt.rosetta_soft_score);
                }
            } else if( do_min==2 ){
                // not minimizing, but will do minimization
                scorefunc_pt[i] = core::scoring::ScoreFunctionFactory::create_score_function(d.opt.rosetta_soft_score);
                if( !d.opt.rosetta_hard_min ){
                    scorefunc_pt[i]->set_weight( core::scoring::fa_rep, scorefunc_pt[i]->get_weight(core::scoring::fa_rep)*0.7 );
                    scorefunc_pt[i]->set_weight( core::scoring::fa_dun, scorefunc_pt[i]->get_weight(core::scoring::fa_dun)*0.7 );
                }
            } else {
                // not minimizing at all, score pass only
                scorefunc_pt[i] = core::scoring::ScoreFunctionFactory::create_score_function(d.opt.rosetta_soft_score);
                if( !d.opt.rosetta_hard_min ){
                    scorefunc_pt[i]->set_weight( core::scoring::fa_rep, scorefunc_pt[i]->get_weight(core::scoring::fa_rep)*1.0 );
                    scorefunc_pt[i]->set_weight( core::scoring::fa_dun, scorefunc_pt[i]->get_weight(core::scoring::fa_dun)*1.0 );
                }
            }
            if( d.target.size() == 1 ){
                // assume this is a ligand, so hbonding is important
                scorefunc_pt[i]->set_weight( core::scoring::fa_elec    , scorefunc_pt[i]->get_weight(core::scoring::fa_elec    )*2.0 );
                scorefunc_pt[i]->set_weight( core::scoring::hbond_sc   , scorefunc_pt[i]->get_weight(core::scoring::hbond_sc   )*2.0 );
                scorefunc_pt[i]->set_weight( core::scoring::hbond_bb_sc, scorefunc_pt[i]->get_weight(core::scoring::hbond_bb_sc)*2.0 );
            } else {
                scorefunc_pt[i]->set_weight( core::scoring::fa_elec    , scorefunc_pt[i]->get_weight(core::scoring::fa_elec    )*1.0 );
                scorefunc_pt[i]->set_weight( core::scoring::hbond_sc   , scorefunc_pt[i]->get_weight(core::scoring::hbond_sc   )*1.0 );
                scorefunc_pt[i]->set_weight( core::scoring::hbond_bb_sc, scorefunc_pt[i]->get_weight(core::scoring::hbond_bb_sc)*1.0 );
            }

            // score one pose single threaded maybe helps lkball issue
            // scorefunc_pt[i]->score(both_per_thread[i]);

            // scorefunc_pt[i]->set_weight( core::scoring::fa_rep, scorefunc_pt[i]->get_weight(core::scoring::fa_rep)*0.67 );
            // scorefunc_pt[i]->set_weight( core::scoring::fa_dun, scorefunc_pt[i]->get_weight(core::scoring::fa_dun)*0.67 );
            // core::scoring::methods::EnergyMethodOptions opts = scorefunc_pt[i]->energy_method_options();
            // core::scoring::hbonds::HBondOptions hopts = opts.hbond_options();
            // hopts.use_hb_env_dep( false );
            // opts.hbond_options( hopts );
            // scorefunc_pt[i]->set_energy_method_options( opts );
            core::kinematics::MoveMapOP movemap = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );
            movemap_pt[i] = movemap;
            minmover_pt[i] = protocols::simple_moves::MinMoverOP(
                new protocols::simple_moves::MinMover( movemap, scorefunc_pt[i], "dfpmin_armijo_nonmonotone", 0.001, true ) );
        }

        std::chrono::time_point<std::chrono::high_resolution_clock> startall = std::chrono::high_resolution_clock::now();

        size_t n_scormin = 0;
        if( minimizing ){
            // min take ~10x score time, so do on 1/10th of the scored
            n_scormin = n_score_calculations * d.opt.rosetta_min_fraction;
            n_scormin = std::ceil(1.0f*n_scormin/omp_max_threads()) * omp_max_threads();
        } else {
            // for scoring, use user cut
            n_scormin = d.opt.rosetta_score_fraction/40.0 * d.total_search_effort;
            if( d.opt.rosetta_score_then_min_below_thresh > -9e8 ){
                for( n_scormin=0; n_scormin < d.packed_results.size(); ++n_scormin ){
                    if( d.packed_results[n_scormin].score > d.opt.rosetta_score_then_min_below_thresh )
                        break;
                }
            }
            n_scormin = std::min<int>( std::max<int>( n_scormin, d.opt.rosetta_score_at_least ), d.opt.rosetta_score_at_most );
            n_scormin = std::min<int>( n_scormin, d.packed_results.size() );
            n_score_calculations = n_scormin;
        }
        d.packed_results.resize(n_scormin);


// Brian injection
        //~~~~~~~~~~~~~~~~~~~~~~~
        // 1. Determine unique scaffolds
        // 2. Populate their correct scaffold-target form

        std::unordered_map<ScaffoldIndex,bool> unique_scaffolds_dict;

        for ( int imin = 0; imin < n_scormin; ++imin ) {
            if ( unique_scaffolds_dict.count(::scheme::kinematics::bigindex_scaffold_index(d.packed_results[imin].index)) == 0) {
                unique_scaffolds_dict[::scheme::kinematics::bigindex_scaffold_index(d.packed_results[imin].index)] = true;
            }
        }
        std::vector<ScaffoldIndex> uniq_scaffolds;
        for ( std::pair<ScaffoldIndex,bool> pair : unique_scaffolds_dict ) uniq_scaffolds.push_back(pair.first);

        int n_uniq = uniq_scaffolds.size();
        std::cout << "Building " << n_uniq << " scaffold+target poses" << std::endl;
        // #ifdef USE_OPENMP
        // #pragma omp parallel for schedule(dynamic,1)
        // #endif
        for ( int iuniq = 0; iuniq < n_uniq; ++iuniq ) {
            ScaffoldIndex si = uniq_scaffolds[iuniq];
            ScaffoldDataCacheOP sdc = d.scaffold_provider->get_data_cache_slow(si);
            if( d.opt.replace_orig_scaffold_res ){
                sdc->setup_both_full_pose(*(d.target.clone()));
            } else {
                sdc->setup_both_pose(*(d.target.clone()));
            }
        }
/////////



        int64_t const out_interval = std::max<int64_t>(1,n_scormin/50);
        if( minimizing) std::cout << "rosetta min on "   << KMGT(n_scormin) << ": ";
        else            std::cout << "rosetta score on " << KMGT(n_scormin) << ": ";
        std::exception_ptr exception = nullptr;

        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic,1)
        #endif
        for( int imin = 0; imin < n_scormin; ++imin )

        {


            try
            {
                if( imin%out_interval==0 ){ cout << '*'; cout.flush();  }

                int const ithread = omp_get_thread_num();

                d.director->set_scene( d.packed_results[imin].index, d.RESLS.size()-1, *d.scene_pt[ithread] );
                EigenXform xposition1 = d.scene_pt[ithread]->position(1);
                EigenXform xalignout = EigenXform::Identity();
                if( d.opt.align_to_scaffold ) xalignout = xposition1.inverse();

                std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
                start = std::chrono::high_resolution_clock::now();

// Brian Injection
                //~~~~~~~~~~~~~~~~~~~~~~~
                // Instead of work_pose, we clone the correct scaffold-target from scene::conformation
                // Get data cache
                // clone pose
                // copy out l2g
                // copy out scaffold_res

                ScaffoldIndex si = ::scheme::kinematics::bigindex_scaffold_index(d.packed_results[imin].index);
                ScaffoldDataCacheOP sdc = d.scaffold_provider->get_data_cache_slow(si);

                core::pose::Pose & pose_to_min( work_pose_pt[ithread] );

                if( d.opt.replace_orig_scaffold_res ){
                    pose_to_min = *(sdc->mpc_both_full_pose.get_pose());
                } else {
                    pose_to_min = *(sdc->mpc_both_pose.get_pose());
                }

                // these guys are multi-thread shared. Definitely don't modify them
                shared_ptr<std::vector<int> const> scaffres_l2g_p = sdc->scaffres_l2g_p;
                shared_ptr<std::vector<bool>> scaffuseres_p = sdc->scaffuseres_p;
                int scaffold_size = scaffuseres_p->size();

//////////
                // pose_to_min = both_per_thread[ithread]; // the old way to get the pose

                // core::pose::PoseOP pose_to_min_ptr = core::pose::PoseOP(new core::pose::Pose);

                xform_pose( pose_to_min, eigen2xyz(xalignout)            , scaffold_size+1 , pose_to_min.size() );
                xform_pose( pose_to_min, eigen2xyz(xalignout*xposition1) ,                      1 ,    scaffold_size );

                // place the rotamers
                core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
                std::vector<bool> is_rif_res(pose_to_min.size(),false);
                for( int ipr = 0; ipr < d.packed_results[imin].numrots(); ++ipr ){
                    int ires = scaffres_l2g_p->at( d.packed_results[imin].rotamers().at(ipr).first );
                    int irot =                  d.packed_results[imin].rotamers().at(ipr).second;
                    core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(d.rot_index.resname(irot)) );
                    pose_to_min.replace_residue( ires+1, *newrsd, true );
                    is_rif_res[ires] = true;
                    for( int ichi = 0; ichi < d.rot_index.nchi(irot); ++ichi ){
                        pose_to_min.set_chi( ichi+1, ires+1, d.rot_index.chi( irot, ichi ) );
                    }
                }

                EigenXform Xtorifframe = xalignout.inverse();
                // clash check existing scaffold rotamers
                // #pragma omp critical
                // {
                //  pose_to_min.dump_pdb("test0.pdb");
                //  xform_pose(pose_to_min, eigen2xyz(Xtorifframe));
                //  pose_to_min.dump_pdb("test1.pdb");
                //  utility_exit_with_message("foo");
                // }

                std::vector<int> replaced_scaffold_res, rifres;
                auto alaop = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map("ALA" ) );

                // pose_to_min.dump_pdb("test_rep_scaff_rots_"+str(imin)+"_before.pdb");

                // The purpose of this loop is to figure out which residues to prune
                std::vector<int> rifatypemap = get_rif_atype_map();
                for( int ir = 1; ir <= scaffold_size; ++ir){
                    auto const & ires = pose_to_min.residue(ir);
                    if( !ires.is_protein() ) continue;
                    if( ires.aa()==core::chemical::aa_gly ||
                        ires.aa()==core::chemical::aa_ala ||
                        ires.aa()==core::chemical::aa_pro ) continue;
                    if(is_rif_res[ir-1]){
                        rifres.push_back(ir);
                        continue;
                    }
                    // if(is_scaffold_fixed_res[ir]) continue;
                    if ( ! (*scaffuseres_p)[ir-1] ) continue;   // i.e. if we weren't set to designable, continue
                    bool ir_clash = false;

                    // check against target_field_by_atype
                    float evtarget = 0.0;
                    for( int ia = 6; ia <= ires.nheavyatoms(); ++ia){
                        auto const & ixyz = ires.xyz(ia);

                        Eigen::Vector3f satm;
                        for(int k = 0; k < 3; ++k) satm[k] = ixyz[k];
                        satm = Xtorifframe * satm;
                        int const irifatype = rifatypemap[ires.atom_type_index(ia)];
                        evtarget += d.target_field_by_atype.at(irifatype)->at(satm);
                    }
                    if( evtarget > 3.0f ) ir_clash = true;

                    // check against other rif res
                    for( int ipr = 0; ipr < d.packed_results[imin].numrots(); ++ipr ){
                        int jr = 1+scaffres_l2g_p->at( d.packed_results[imin].rotamers().at(ipr).first );
                        auto const & jres = pose_to_min.residue(jr);
                        // should do rsd nbr check... but speed not critical here ATM...
                        for( int ia = 6; ia <= ires.nheavyatoms(); ++ia){
                            auto const & ixyz = ires.xyz(ia);
                            for( int ja = 6; ja <= jres.nheavyatoms(); ++ja){
                                auto const & jxyz = jres.xyz(ja);
                                if( ixyz.distance_squared(jxyz) < 9.0 ){
                                    ir_clash = true;
                                }
                            }
                        }
                    }
                    if(ir_clash){
                        pose_to_min.replace_residue(ir, *alaop, true);
                        replaced_scaffold_res.push_back(ir);
                    }
                }

                // pose_to_min.dump_pdb("test_rep_scaff_rots_"+str(imin)+"_after.pdb");
                // if(imin >= 9) utility_exit_with_message("testing...aorsenoi");


                end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed_seconds_copypose = end-start;
                #pragma omp critical
                time_copy += elapsed_seconds_copypose;
                // pose_to_min.dump_pdb("test_pre.pdb");

                if( minimizing ){

// Brian injection
                    //~~~~~~~~~~~~~~~~~~~~~~~
                    // get the movemap from the movemap list
                    // set the movemap like normal
                    // add it into the minmover_pt 

                    core::kinematics::MoveMapOP movemap = movemap_pt[ithread];
                    movemap->set_chi(true);
                    movemap->set_jump(true);
                    for(int ir = 1; ir <= pose_to_min.size(); ++ir){
                        bool is_scaffold = ir <= scaffold_size;
                        if( is_scaffold ) movemap->set_bb(ir, d.opt.rosetta_min_allbb || d.opt.rosetta_min_scaffoldbb );
                        else              movemap->set_bb(ir, d.opt.rosetta_min_allbb || d.opt.rosetta_min_targetbb );
                        if( d.opt.rosetta_min_fix_target && !is_scaffold ){
                            movemap->set_chi(ir,false);
                        }
                    }
                    minmover_pt[ithread]->set_movemap(movemap);
////////////

                    // std::cout << "MIN!" << std::endl;
                    // start = std::chrono::high_resolution_clock::now();
                    // #pragma omp critical
                    // {
                       //  pose_to_min.dump_pdb("min_"+str(imin)+"_before.pdb");
                        minmover_pt[ithread]->apply( pose_to_min );
                    //     pose_to_min.dump_pdb("min_"+str(imin)+"_after.pdb");
                    //     utility_exit_with_message("test_min");
                    // }
                    end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed_seconds_min = end-start;
                    #pragma omp critical
                    time_min += elapsed_seconds_min;
                } else {
                    // std::cout << "SCORE!" << std::endl;
                    start = std::chrono::high_resolution_clock::now();
                    scorefunc_pt[ithread]->score( pose_to_min );
                    end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed_seconds_score = end-start;
                    #pragma omp critical
                    time_score += elapsed_seconds_score;
                }

                // if( target.size() == 1 ){
                //  // ligand!
                //  float total_lj_neg  = scorefunc_pt[ithread]->get_weight(core::scoring::fa_atr) * pose_to_min.energies().total_energies()[core::scoring::fa_atr];
                //        total_lj_neg += scorefunc_pt[ithread]->get_weight(core::scoring::fa_rep) * pose_to_min.energies().total_energies()[core::scoring::fa_rep];
                //        total_lj_neg = std::max(0.0f,total_lj_neg);
                //  d.packed_results[imin].score  = 1.00*total_lj_neg;
                //  d.packed_results[imin].score += 1.00*pose_to_min.energies().total_energies()[core::scoring::hbond_sc]; // last res is ligand
                //  d.packed_results[imin].score += 1.00*pose_to_min.energies().residue_total_energy(pose_to_min.size());
                //  for( int ipr = 0; ipr < d.packed_results[imin].numrots(); ++ipr ){
                //      int ires = scaffres_l2g_p->at( d.packed_results[imin].rotamers().at(ipr).first );
                //      d.packed_results[imin].score += 0.5*pose_to_min.energies().residue_total_energy(ires);
                //  }
                // } else {
                    // not ligand! compute the fixed-everything ddg
                if( d.opt.rosetta_score_total ){
                    d.packed_results[imin].score = pose_to_min.energies().total_energy();
                } else {
                    double rosetta_score = 0.0;
                    auto const & weights = pose_to_min.energies().weights();
                    if( !d.opt.rosetta_score_ddg_only ){
                        for( int ir = 1; ir <= scaffold_size; ++ir ){
                            if( is_rif_res[ir-1] ){
                                rosetta_score += pose_to_min.energies().onebody_energies(ir).dot(weights);
                            }
                        }
                    }
                    if( !d.opt.rosetta_score_ddg_only && d.target.size()==1 ){
                        // is ligand, add it's internal energy
                        rosetta_score += pose_to_min.energies().onebody_energies(pose_to_min.size()).dot(weights);
                    }
                    auto const & egraph = pose_to_min.energies().energy_graph();
                    for(int ir = 1; ir <= egraph.num_nodes(); ++ir){
                        for ( utility::graph::Graph::EdgeListConstIter
                                iru  = egraph.get_node(ir)->const_upper_edge_list_begin(),
                                irue = egraph.get_node(ir)->const_upper_edge_list_end();
                                iru != irue; ++iru
                        ){
                            EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
                            int jr = edge.get_second_node_ind();

                            // this is DDG
                            if( ir <= scaffold_size && jr > scaffold_size ){
                                // ir in scaff, jr in target
                                rosetta_score += edge.dot(weights);
                            }
                            if( !d.opt.rosetta_score_ddg_only && jr <= scaffold_size ){
                                // ir & jr in scaffold
                                if( is_rif_res[ir-1] || is_rif_res[jr-1] ){
                                    double const edgescore = edge.dot(weights);
                                    if( edgescore > 0.0 ){
                                        // always assess full score for bad interactions
                                        rosetta_score += edgescore;
                                    } else if( is_rif_res[ir-1] && is_rif_res[jr-1] ){
                                        // both rif residues
                                        rosetta_score += d.opt.rosetta_score_rifres_rifres_weight * edgescore;
                                        // bonus for hbonds between rif residues
                                        rosetta_score += edge[core::scoring::hbond_sc];
                                    } else {
                                        // rest: one rif res, one other scaff res
                                        rosetta_score += d.opt.rosetta_score_rifres_scaffold_weight * edgescore;
                                    }
                                } else {
                                    // scaffold / scaffold ignored
                                }
                            }
                        }
                    }
                    d.packed_results[imin].score = rosetta_score;
                    // #pragma omp critical
                    // std::cout << rosetta_score << std::endl;
                }


                if( (minimizing+1 == do_min)     && d.packed_results[imin].score < d.opt.rosetta_score_cut ){
                    d.packed_results[imin].pose_ = core::pose::PoseOP( new core::pose::Pose(pose_to_min) );
                    for(int ir : rifres){
                        d.packed_results[imin].pose_->pdb_info()->add_reslabel(ir, "RIFRES" );
                    }
                    for(int ir : replaced_scaffold_res){
                        d.packed_results[imin].pose_->pdb_info()->add_reslabel(ir, "PRUNED" );
                    }
                }

                // #pragma omp critical
                // pose_to_min.energies().show(cout,pose_to_min.size());
                // pose_to_min.dump_pdb("test_post.pdb");
                // #pragma omp critical
                // std::cout << imin << " prescore: " << e_pre_min << " minscore: " << pose_to_min.energies().total_energy() << std::endl;
                // utility_exit_with_message("testing");
            } catch(...) {
                #pragma omp critical
                exception = std::current_exception();
            }

        } // end of OMP loop
        if( exception ) std::rethrow_exception(exception);


        cout << endl;
        __gnu_parallel::sort( d.packed_results.begin(), d.packed_results.end() );
        {
            size_t n_scormin = 0;
            for( n_scormin; n_scormin < d.packed_results.size(); ++n_scormin ){
                if( minimizing && d.packed_results[n_scormin].score > d.opt.rosetta_score_cut ) break;
            }
            d.packed_results.resize(n_scormin);
        }

        std::chrono::time_point<std::chrono::high_resolution_clock> stopall = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elap_sec = stopall - startall;

        if( minimizing ){
            std::cout << "total pose copy time: " << time_copy.count() << "s total min time: " << time_min.count() << "s" << " walltime total " << elap_sec.count() << std::endl;
            std::cout << "min rate: "   << (double)n_scormin / elap_sec.count()                     << " sc/rb minimizations per second" << std::endl;
            std::cout << "min rate: "   << (double)n_scormin / elap_sec.count() / omp_max_threads() << " sc/rb minimizations per second per thread" << std::endl;
        } else {
            std::cout << "total pose copy time: " << time_copy.count() << "s total score time: " << time_score.count() << "s" << " walltime total " << elap_sec.count() << std::endl;
            std::cout << "score rate: " << (double)n_scormin / elap_sec.count()                     << " rosetta scores per second" << std::endl;
            std::cout << "score rate: " << (double)n_scormin / elap_sec.count() / omp_max_threads() << " rosetta scores per second per thread" << std::endl;
        }
    }

    std::chrono::duration<double> elapsed_seconds_rosetta = std::chrono::high_resolution_clock::now()-start_rosetta;
    d.time_ros += elapsed_seconds_rosetta.count();





}




#endif