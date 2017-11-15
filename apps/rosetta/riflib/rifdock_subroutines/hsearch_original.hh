// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_hsearch_original_hh
#define INCLUDED_riflib_rifdock_subroutines_hsearch_original_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/rifdock_subroutines/meta.hh>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;

template<class DirectorBase, class ScaffoldProvider >
struct HsearchData {
    RifDockOpt & opt;
    std::vector<float> & RESLS;
    DirectorBase & director;
    int64_t & total_search_effort;
    std::vector< devel::scheme::ScenePtr > & scene_pt;
    devel::scheme::ScenePtr & scene_minimal;
    // Eigen::Vector3f & scaffold_center;
    float & target_redundancy_filter_rg;
    // core::pose::Pose & scaffold_centered;
    core::pose::Pose & target;
    // std::vector<devel::scheme::EigenXform> & symmetries_clash_check;
    // std::vector< devel::scheme::SimpleAtom > & scaffold_simple_atoms;
    devel::scheme::RotamerIndex & rot_index;
    // std::vector< std::vector< devel::scheme::VoxelArrayPtr > > & scaffold_bounding_by_atype;
    std::vector< devel::scheme::ObjectivePtr > & objectives;
    int64_t & non0_space_size;

    ScaffoldProvider & scaffold_provider;

};

// template<__DirectorBase>
// using HsearchFunctionType = typedef


template<class DirectorBase, class ScaffoldProvider >
bool
hsearch_original(
    shared_ptr<std::vector< _SearchPointWithRots<DirectorBase> > > & hsearch_results_p,
    HsearchData<DirectorBase, ScaffoldProvider > & d) {


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

    typedef _DirectorBigIndex<DirectorBase> DirectorIndex;
    typedef tmplSearchPoint<DirectorIndex> SearchPoint;

    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;

//////////////////////////////////////////////////////////////////////////////////////


            // these are in rosetta numbers
            // int rmsd_resid1 = 107;
            // int rmsd_resid2 = 54;
            // int rmsd_resid3 = 13;


            // utility::vector1<core::Size> rmsd_resids { rmsd_resid1, rmsd_resid2, rmsd_resid3 };
            // std::vector<SchemeAtom> rmsd_atoms;
            // devel::scheme::get_scheme_atoms_cbonly( scaffold_unmodified_from_file, rmsd_resids, rmsd_atoms);
            
            // SchemeAtom rmsd_cb1 = rmsd_atoms[0];
            // SchemeAtom rmsd_cb2 = rmsd_atoms[1];
            // SchemeAtom rmsd_cb3 = rmsd_atoms[2];

            // std::cout << "Checking that we have the right CBetas for rmsd" << std::endl;
            // std::cout << I(5, rmsd_resid1) << " " << rmsd_cb1.position().transpose() << std::endl;
            // std::cout << I(5, rmsd_resid2) << " " << rmsd_cb2.position().transpose() << std::endl;
            // std::cout << I(5, rmsd_resid3) << " " << rmsd_cb3.position().transpose() << std::endl;


//////////////////////////////////////////////////////////////////////////////////////


    ScaffoldDataCacheOP sdc = d.scaffold_provider.get_data_cache_slow(scaffold_index_default_value( ScaffoldIndex()));
    float redundancy_filter_rg = sdc->get_redundancy_filter_rg( d.target_redundancy_filter_rg );
    Eigen::Vector3f scaffold_center = sdc->scaffold_center;


    std::vector< std::vector< SearchPoint > > samples( d.RESLS.size() );

    bool search_failed = false;
    {
        samples[0].resize( ::scheme::kinematics::bigindex_nest_part(d.director->size(0)) );
        for( uint64_t i = 0; i < ::scheme::kinematics::bigindex_nest_part(d.director->size(0)); ++i ) samples[0][i] = SearchPoint( DirectorIndex( i, scaffold_index_default_value( ScaffoldIndex())) );
        BOOST_FOREACH( ScenePtr & s, d.scene_pt ) s = d.scene_minimal->clone_specific_deep(std::vector<uint64_t> {1});
        for( int iresl = 0; iresl < d.RESLS.size(); ++iresl )
        {
            cout << "HSearsh stage " << iresl+1 << " resl " << F(5,2,d.RESLS[iresl]) << " begin threaded sampling, " << KMGT(samples[iresl].size()) << " samples: ";
            int64_t const out_interval = samples[iresl].size()/50;
            std::exception_ptr exception = nullptr;
            std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
            start = std::chrono::high_resolution_clock::now();
            d.total_search_effort += samples[iresl].size();

            std::vector< double > rmsds( samples[iresl].size(), 1000 );

            bool answer_exists = false;

            #ifdef USE_OPENMP
            #pragma omp parallel for schedule(dynamic,64)
            #endif
            for( int64_t i = 0; i < samples[iresl].size(); ++i ){
                if( exception ) continue;
                try {
                    if( i%out_interval==0 ){ cout << '*'; cout.flush(); }
                    DirectorIndex const isamp = samples[iresl][i].index;

///////////////////////////////////////////////////////////////////////////////////////////////

                    // if ( isamp < 3000000000000000 && isamp != 47729600827993 && 
                    //   isamp != 745775012937 && isamp != 11652734577 &&
                    //   isamp != 182073977 && isamp != 2844905 ) {
                    //  continue;
                    // }
                    // if ( isamp > 3000000000000000 && isamp != 3054694452991568 ) {
                    //  continue;
                    // }


                    // money 3054694452991568
///////////////////////////////////////////////////////////////////////////////////////////////


                    ScenePtr tscene( d.scene_pt[omp_get_thread_num()] );
                    d.director->set_scene( isamp, iresl, *tscene );

                    if( d.opt.tether_to_input_position ){
                        EigenXform x = tscene->position(1);
                        x.translation() -= scaffold_center;
                        float xmag =  xform_magnitude( x, redundancy_filter_rg );
                        if( xmag > d.opt.tether_to_input_position_cut + d.RESLS[iresl] ){
                            samples[iresl][i].score = 9e9;
                            continue;
                        } else {
                            // std::cout << "inbounds " << iresl << " " << xform_magnitude( tscene->position(1), d.redundancy_filter_rg ) << std::endl;
                        }
                    }

    // This is not debug code, this is the old symmetry clash check 
    // This needs to be reimplemented after the ScaffoldProvider merge
                    // float tot_sym_score = 0;
                    // if( d.opt.nfold_symmetry > 1 ){
                    //     bool dump = false;
                    //     if( iresl == 2 ) dump = true;
                    //     if(dump){
                    //         d.scaffold_centered.dump_pdb("test_scaff.pdb");
                    //         d.target.dump_pdb("test_target.pdb");
                    //     }
                    //     EigenXform x = tscene->position(1);
                    //     EigenXform xinv = x.inverse();
                    //     for(int isym = 0; isym < d.symmetries_clash_check.size(); ++isym){
                    //         EigenXform const & xsym = d.symmetries_clash_check[isym];
                    //     // for( EigenXform const & xsym : d.symmetries_clash_check ){
                    //         utility::io::ozstream * outp = nullptr;
                    //         if(dump){
                    //             outp = new utility::io::ozstream("test_sym_"+str(isym)+".pdb");
                    //         }
                    //         EigenXform x_to_internal = xinv * xsym * x;
                    //         for( SimpleAtom a : d.scaffold_simple_atoms ){
                    //             a.set_position( x_to_internal * a.position() );
                    //             if(dump) ::scheme::actor::write_pdb( *outp, a, d.rot_index.chem_index_ );
                    //             float atom_score = d.scaffold_bounding_by_atype.at(iresl).at(5)->at(a.position());
                    //             if( atom_score > 0.0 ){ // clash only
                    //                 tot_sym_score += atom_score;
                    //             }
                    //         }
                    //         if(dump){
                    //             outp->close();
                    //             delete outp;
                    //         }
                    //     }
                    //     if(dump){
                    //         std::cout << "tot_sym_score " << tot_sym_score << std::endl;
                    //         utility_exit_with_message("testing symmetric clash check");
                    //     }
                    // }


                    // the real rif score!!!!!!
                    samples[iresl][i].score = d.objectives[iresl]->score( *tscene );// + tot_sym_score;




///////////////////////////////////////////////////////////////////////////////////////////

                    // if ( isamp == 3054694452991568 || isamp == 47729600827993 ||
                    //   isamp == 745775012937 || isamp == 11652734577 ||
                    //   isamp == 182073977 || isamp == 2844905 ) {


//                          // samples[iresl][i].score = d.objectives[iresl]->score( *tscene ) + tot_sym_score;
                    
//                          // float score = d.objectives[iresl]->score( *tscene ) + tot_sym_score;
                    //  // answer_exists = true;

                    //  #pragma omp critical
                    //  std::cout << "Score for the one: " << F(6, 2, samples[iresl][i].score) << std::endl;
                    // }



//                         const SimpleAtom scene_cb1 = tscene->template get_actor<SimpleAtom>(1, rmsd_resid1-1);
//                         const SimpleAtom scene_cb2 = tscene->template get_actor<SimpleAtom>(1, rmsd_resid2-1);
//                         const SimpleAtom scene_cb3 = tscene->template get_actor<SimpleAtom>(1, rmsd_resid3-1);

//                         const float cb1_dist2 = ( scene_cb1.position() - rmsd_cb1.position() ).squaredNorm();
//                         const float cb2_dist2 = ( scene_cb2.position() - rmsd_cb2.position() ).squaredNorm();
//                         const float cb3_dist2 = ( scene_cb3.position() - rmsd_cb3.position() ).squaredNorm();

//                         const float rmsd_squared = cb1_dist2 + cb2_dist2 + cb3_dist2;
//                         // const float rmsd = std::sqrt( (cb1_dist * cb1_dist + cb2_dist * cb2_dist + cb3_dist * cb3_dist) / 3.00 );

//                         rmsds[i] = rmsd_squared;

//                         samples[iresl][i].score = rmsd_squared - 200.0;

//                         if (isamp > 3000000000000000) {
//                          // samples[iresl][i].score = d.objectives[iresl]->score( *tscene ) + tot_sym_score;
//                         }


                    // if ( isamp == 47780569615988) {
                    //  double score = d.objectives[iresl]->score( *tscene ) + tot_sym_score;
                    //  std::cout << "Score for the one: " << score << std::endl;
                    // }

                    // #pragma omp critical
                    // std::cout << I(20, isamp) << F(7, 1, cb1_dist) << F(7, 1, cb2_dist) << F(7, 1, cb3_dist) << std::endl;

                    // SimpleAtom sa = tscene->template get_actor<SimpleAtom>(1,0);
                    // SimpleAtom sa2 = tscene->template get_actor<SimpleAtom>(1,13);
                    // bool success = true;
                    // // // bool success = tscene->get_actor( 1, 0, sa ); 
                    // if ( success ) {
                    //  #pragma omp critical
                    //  std::cout << "Succss: " << I(8, isamp) << " SA: " << sa << sa2 << std::endl;
                    // } else {
                    //  #pragma omp critical
                    //  std::cout << "Failure" << std::endl;
                    // }


////////////////////////////////////////////////////////////////////////////////////////////



                } catch( std::exception const & ex ) {
                    #ifdef USE_OPENMP
                    #pragma omp critical
                    #endif
                    exception = std::current_exception();
                }
            }
            if( exception ) std::rethrow_exception(exception);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_seconds_rif = end-start;
            float rate = (double)samples[iresl].size()/ elapsed_seconds_rif.count()/omp_max_threads();
            cout << endl;// << "done threaded sampling, partitioning data..." << endl;

//////////////////////////////////////////////////////////////////////////////////////

            // shared_ptr<DirectorOriTrans6D> nest_director = std::dynamic_pointer_cast<DirectorOriTrans6D>(d.director);
            // NestOriTrans6D const & nest( nest_director->nest() );

            // std::cout << "Sorting rmsds... " << std::endl;
            // std::vector< uint64_t > sorted( rmsds.size() );
            // uint64_t n = 0;
            // std::generate( std::begin(sorted), std::end(sorted), [&]{ return n++; });
            // std::sort( std::begin( sorted), std::end(sorted),
            //  [&](int i1, int i2) { return rmsds[i1] < rmsds[i2]; } );

            // for ( uint64_t i = 0; i < 64; i++ ) {
            //  uint64_t index = sorted[i];
            //  SearchPoint sp = samples[iresl][index];
            //  uint64_t isamp = sp.index;


            //  EigenXform trans;
            //  bool success = nest.get_state( isamp, iresl, trans );

            //  Eigen::Matrix<float,3,3> rot;
            //  scheme::objective::hash::get_transform_rotation<float>( trans, rot );

            //  // Eigen::Matrix<double,3,3> rot2 = rot.cast<double>();

            //  Eigen::Vector3f axis;
            //  axis.setRandom();
            //  axis.normalize();
            //  Eigen::AngleAxisf angle_axis( rot );


            //  std::cout << I(20, isamp) << F(7, 1, rmsds[index]) << "  Trans: " 
            //  << F(7, 1, trans.translation()[0]) 
            //  << F(7, 1, trans.translation()[1]) 
            //  << F(7, 1, trans.translation()[2]) 
            //  << "  Angle: " << F(7, 1, angle_axis.angle()*180.0/M_PI) 
            //  << "  Score: " << F(7, 2, sp.score)

            //  << std::endl;
            // }
////////////////////////////////////////////////////////////////////////////////////////


            SearchPoint max_pt, min_pt;
            int64_t len = samples[iresl].size();
            if( samples[iresl].size() > d.opt.beam_size/d.opt.DIMPOW2 ){
                __gnu_parallel::nth_element( samples[iresl].begin(), samples[iresl].begin()+d.opt.beam_size/d.opt.DIMPOW2, samples[iresl].end() );
                len = d.opt.beam_size/d.opt.DIMPOW2;
                min_pt = *__gnu_parallel::min_element( samples[iresl].begin(), samples[iresl].begin()+len );
                max_pt = *(samples[iresl].begin()+d.opt.beam_size/d.opt.DIMPOW2);
            } else {
                min_pt = *__gnu_parallel::min_element( samples[iresl].begin(), samples[iresl].end() );
                max_pt = *__gnu_parallel::max_element( samples[iresl].begin(), samples[iresl].end() );
            }

            cout << "HSearsh stage " << iresl+1 << " complete, resl. " << F(7,3,d.RESLS[iresl]) << ", "
                  << " " << KMGT(samples[iresl].size()) << ", promote: " << F(9,6,min_pt.score) << " to "
                  << F(9,6, std::min(d.opt.global_score_cut,max_pt.score)) << " rate " << KMGT(rate) << "/s/t " << std::endl;

            // cout << "Answer: " << ( answer_exists ? "exists" : "doesn't exist" ) << std::endl;

            if( iresl+1 == samples.size() ) break;

            for( int64_t i = 0; i < len; ++i ){
                uint64_t isamp0 = ::scheme::kinematics::bigindex_nest_part(samples[iresl][i].index);
                if( samples[iresl][i].score >= d.opt.global_score_cut ) continue;
                if( iresl == 0 ) ++d.non0_space_size;
                for( uint64_t j = 0; j < d.opt.DIMPOW2; ++j ){
                    uint64_t isamp = isamp0 * d.opt.DIMPOW2 + j;
                    samples[iresl+1].push_back( SearchPoint(DirectorIndex(isamp, scaffold_index_default_value( ScaffoldIndex()))) );
                }
            }
            if( 0 == samples[iresl+1].size() ){
                search_failed = true;
                std::cout << "search fail, no valid samples!" << std::endl;
                break;
            }
            samples[iresl].clear();

        }
        if( search_failed ) return false;
        std::cout << "full sort of final samples" << std::endl;
        __gnu_parallel::sort( samples.back().begin(), samples.back().end() );
    }
    if( search_failed ) return false;

    std::cout << "total non-0 space size was approx " << float(d.non0_space_size)*1024.0*1024.0*1024.0 << " grid points" << std::endl;
    std::cout << "total search effort " << KMGT(d.total_search_effort) << std::endl;


    hsearch_results_p = make_shared<std::vector< SearchPointWithRots >>();
    std::vector< SearchPointWithRots > & hsearch_results = *hsearch_results_p;


    hsearch_results.resize( samples.back().size() );
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1024)
    #endif
    for( int ipack = 0; ipack < hsearch_results.size(); ++ipack ){
        hsearch_results[ipack].score = samples.back()[ipack].score;
        hsearch_results[ipack].index = samples.back()[ipack].index;
    }



    return true;



}




#endif