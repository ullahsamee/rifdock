// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/CompileAndFilterResultsTasks.hh>

#include <riflib/types.hh>
#include <riflib/task/util.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>

#include <string>
#include <vector>
#include <unordered_map>



namespace devel {
namespace scheme {


shared_ptr<std::vector<RifDockResult>>
CompileAndFilterResultsTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss = search_point_with_rotss_from_rif_dock_results( rif_dock_results );
    return return_rif_dock_results( search_point_with_rotss, rdd, pd );
}


shared_ptr<std::vector<RifDockResult>>
CompileAndFilterResultsTask::return_rif_dock_results( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss = search_point_with_rotss_from_search_points( search_points );
    return return_rif_dock_results( search_point_with_rotss, rdd, pd );
}


typedef int32_t intRot;

template<class EigenXform, class ScaffoldIndex>
struct tmplXRtriple {
    EigenXform xform;
    ScaffoldIndex scaffold_index;
    uint64_t result_num;
};


// how can I fix this??? make the whole prototype into a class maybe???
// what does it do?
//  set and rescore scene with nopackscore, record more score detail
//  compute dist0
//  select results with some redundancy filtering
template<
    class EigenXform,
    // class Scene,
    class ScenePtr,
    class ObjectivePtr,
    class ScaffoldIndex
>
void
awful_compile_output_helper_(
    int64_t isamp,
    int director_resl,
    std::vector< SearchPointWithRots > const & packed_results,
    std::vector< ScenePtr > & scene_pt,
    DirectorBase director,
    float redundancy_filter_rg,
    float redundancy_filter_mag,
    Eigen::Vector3f scaffold_center,
    std::vector< std::vector< RifDockResult > > & allresults_pt,
                 std::vector< RifDockResult >   & selected_results,
    std::vector< tmplXRtriple<EigenXform, ScaffoldIndex> > & selected_xforms,
    int n_pdb_out,
    #ifdef USE_OPENMP
        omp_lock_t & dump_lock,
    #endif
    ObjectivePtr objective,
    int & nclose,
    int nclosemax,
    float nclosethresh,
    EigenXform scaffold_perturb
) {

    typedef tmplXRtriple<EigenXform, ScaffoldIndex> XRtriple;

    SearchPointWithRots const & sp = packed_results[isamp];
    if( sp.score >= 0.0f ) return;
    ScenePtr scene_minimal( scene_pt[omp_get_thread_num()] );
    director->set_scene( sp.index, director_resl, *scene_minimal );
    std::vector<float> sc = objective->scores(*scene_minimal);
    float const nopackscore = sc[0]+sc[1]; //result.sum();
    float const rifscore = sc[0]; //result.template get<MyScoreBBActorRIF>();
    float const stericscore = sc[1]; //result.template get<MyClashScore>();

    // dist0 is only important to the nclose* options
    float dist0; {
        EigenXform x = scene_minimal->position(1);
        x = scaffold_perturb * x;
        x.translation() -= scaffold_center;
        dist0 = ::devel::scheme::xform_magnitude( x, redundancy_filter_rg );
    }

    RifDockResult r; // float dist0, packscore, nopackscore, rifscore, stericscore;
    r.isamp = isamp;
    r.prepack_rank = sp.prepack_rank;
    r.index = sp.index;
    r.score = sp.score;
    r.nopackscore = nopackscore;
    r.rifscore = rifscore;
    r.stericscore = stericscore;
    r.dist0 = dist0;
    r.cluster_score = 0.0;
    r.pose_ = sp.pose_;
    allresults_pt.at( omp_get_thread_num() ).push_back( r ); // recorded w/o rotamers here

    bool force_selected = ( dist0 < nclosethresh && ++nclose < nclosemax ); // not thread-safe... is this important?

    if( selected_xforms.size() < n_pdb_out || force_selected ){

        EigenXform xposition1 = scene_minimal->position(1);
        EigenXform xposition1inv = xposition1.inverse();

        float mindiff_candidate = 9e9;
        int64_t i_closest_result;
        BOOST_FOREACH( XRtriple const & xrp, selected_xforms ){
            EigenXform const & xsel = xrp.xform;
            EigenXform const xdiff = xposition1inv * xsel;
            float diff = devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg );
            if( diff < mindiff_candidate ){
                mindiff_candidate = diff;
                i_closest_result = xrp.result_num;
            }
            // todo: also compare AA composition of rotamers
        }

        if( mindiff_candidate < redundancy_filter_mag ){ // redundant result
            selected_results[i_closest_result].cluster_score += 1.0; //sp.score==0.0 ? nopackscore : sp.score;
        }

        if( mindiff_candidate > redundancy_filter_mag || force_selected ){

            #ifdef USE_OPENMP
            omp_set_lock( &dump_lock );
            #endif
            {
                // std::cout << "checking again to add selected " << selected_xforms.size() << " " << omp_get_thread_num() << std::endl;
                float mindiff_actual = 9e9;
                BOOST_FOREACH( XRtriple const & xrp, selected_xforms ){
                    EigenXform const & xsel = xrp.xform;
                    EigenXform const xdiff = xposition1inv * xsel;
                    mindiff_actual = std::min( mindiff_actual, devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg ) );
                }
                if( mindiff_actual > redundancy_filter_mag || force_selected ){
                    if( redundancy_filter_mag > 0.0001 ) {
                        selected_xforms.push_back( XRtriple {
                            xposition1, 
                            sp.index,
                            (int64_t)selected_results.size()
                        } );
                    }
                    r.rotamers_ = sp.rotamers_;
                    selected_results.push_back( r ); // recorded with rotamers here
                } else {
                    // std::cout << " second check failed" << std::endl;
                }
            }
            #ifdef USE_OPENMP
            omp_unset_lock( &dump_lock );
            #endif

        } // end if( mindiff > redundancy_filter_mag ){

    } // end    if( selected_xforms.size() < n_pdb_out || force_selected )

}



shared_ptr<std::vector<RifDockResult>> 
CompileAndFilterResultsTask::return_rif_dock_results( 
    shared_ptr<std::vector<SearchPointWithRots>> packed_results_p, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::vector<SearchPointWithRots> & packed_results = *packed_results_p;
    std::vector< RifDockResult > allresults;

    shared_ptr<std::vector<RifDockResult>> selected_results_p = make_shared<std::vector<RifDockResult>>();
    std::vector< RifDockResult > & selected_results = *selected_results_p;

    using std::cout;
    using std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    print_header( "compile and filter results" ); ///////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    typedef tmplXRtriple<EigenXform, RifDockIndex> XRtriple;

    int64_t Nout = packed_results.size(); 

    std::vector< std::vector< RifDockResult > > allresults_pt( omp_max_threads() );

    SelectiveRifDockIndexHasher   hasher( false, filter_seeding_positions_separately_, filter_scaffolds_separately_ );
    SelectiveRifDockIndexEquater equater( false, filter_seeding_positions_separately_, filter_scaffolds_separately_ );

    std::unordered_map< RifDockIndex, std::vector< XRtriple >, SelectiveRifDockIndexHasher, SelectiveRifDockIndexEquater > 
        selected_xforms_map(1000, hasher, equater);  // default value here needs to be .reserve(65536)
    std::unordered_map< RifDockIndex, int, SelectiveRifDockIndexHasher, SelectiveRifDockIndexEquater > 
        nclose_map(1000, hasher, equater); // default value here needs to be 0

    for ( uint64_t isamp = 0; isamp < Nout; isamp++ ) {
        RifDockIndex rdi = packed_results[isamp].index;
        if ( selected_xforms_map.count(rdi) == 0 ) {
            selected_xforms_map[ rdi ].reserve(65536); // init big to reduce liklihood of resizes
            nclose_map[ rdi ] = 0;
        }
    }


    ////////////////////


    int nclosemax      = force_output_if_close_to_input_num_;
    float nclosethresh = force_output_if_close_to_input_;

    std::cout << "redundancy_filter_mag " << redundancy_mag_ << "A \"rmsd\"" << std::endl;
    int64_t Nout_singlethread = std::min( (int64_t)10000, Nout );

    std::cout << "going throuth 10K results (1 thread): ";
    int64_t out_interval = 10000/81;
    for( int64_t isamp = 0; isamp < Nout_singlethread; ++isamp ){
        if( isamp%out_interval==0 ){ cout << '*'; cout.flush(); }

        RifDockIndex rdi = packed_results[isamp].index;
        ScaffoldIndex si = packed_results[isamp].index.scaffold_index;
        std::vector< XRtriple > & selected_xforms = selected_xforms_map.at( rdi );
        int & nclose = nclose_map.at( rdi );
        ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(si);
        float redundancy_filter_rg = sdc->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );
        EigenXform scaffold_perturb = sdc->scaffold_perturb;
        Eigen::Vector3f scaffold_center = sdc->scaffold_center;

                                        
        awful_compile_output_helper_< EigenXform, ScenePtr, ObjectivePtr >(
            isamp, director_resl_, packed_results, rdd.scene_pt, rdd.director,
            redundancy_filter_rg, redundancy_mag_, scaffold_center,
            allresults_pt, selected_results, selected_xforms, n_per_block_,
            #ifdef USE_OPENMP
                rdd.dump_lock,
            #endif
            rdd.objectives.at(rif_resl_), nclose, nclosemax, nclosethresh,
            scaffold_perturb
        );

            std::cout << selected_xforms.size() << std::endl;
    }
    std::cout << std::endl;

    std::cout << "going throuth all results (threaded): ";
    out_interval = Nout / 82;
    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,8)
    #endif
    for( int64_t isamp = Nout_singlethread; isamp < Nout; ++isamp ){
        if( exception ) continue;
        try{
            if( isamp%out_interval==0 ){ cout << '*'; cout.flush(); }

            RifDockIndex rdi = packed_results[isamp].index;
            ScaffoldIndex si = packed_results[isamp].index.scaffold_index;
            std::vector< XRtriple > & selected_xforms = selected_xforms_map.at( rdi );
            int & nclose = nclose_map.at( rdi );
            ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(si);
            float redundancy_filter_rg = sdc->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );
            EigenXform scaffold_perturb = sdc->scaffold_perturb;
            Eigen::Vector3f scaffold_center = sdc->scaffold_center;


            awful_compile_output_helper_< EigenXform, ScenePtr, ObjectivePtr >(
                isamp, director_resl_, packed_results, rdd.scene_pt, rdd.director,
                redundancy_filter_rg, redundancy_mag_, scaffold_center,
                allresults_pt, selected_results, selected_xforms, n_per_block_,
                #ifdef USE_OPENMP
                    rdd.dump_lock,
                #endif
                rdd.objectives.at(rif_resl_), nclose, nclosemax, nclosethresh,
                scaffold_perturb
            );
        } catch(...) {
            #pragma omp critical
            exception = std::current_exception();
        }
    }
    if( exception ) std::rethrow_exception(exception);
    std::cout << std::endl;

    std::cout << "sort compiled results" << std::endl;
    BOOST_FOREACH( std::vector<RifDockResult> const & rs, allresults_pt ){
        BOOST_FOREACH( RifDockResult const & r, rs ){
            allresults.push_back( r );
        }
    }
    __gnu_parallel::sort( allresults.begin(), allresults.end() );


    return selected_results_p;

}



}}
