// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/HackPackTasks.hh>

#include <riflib/types.hh>
#include <riflib/util.hh>
#include <riflib/ScoreRotamerVsTarget.hh>


#include <string>
#include <vector>



namespace devel {
namespace scheme {


shared_ptr<std::vector<SearchPoint>> 
FilterForHackPackTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
FilterForHackPackTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
FilterForHackPackTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
FilterForHackPackTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    size_t n_packsamp = 0;
    for( n_packsamp; n_packsamp < any_points->size(); ++n_packsamp ){
<<<<<<< HEAD
        if( (*any_points)[n_packsamp].score > global_score_cut_ ) break;
=======
        if( (*any_points)[n_packsamp].score > hack_pack_score_cut_ ) break;
>>>>>>> rifdock/master
    }
    
    pd.npack = std::min( n_packsamp, (size_t)(pd.total_search_effort *
        ( hack_pack_frac_ / (pack_n_iters_*pack_iter_mult_)) ) );

    if ( hack_pack_frac_ >= 1 ) pd.npack = n_packsamp;

    any_points->resize( pd.npack );

    return any_points;
}
    

shared_ptr<std::vector<SearchPointWithRots>>
HackPackTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> packed_results_p, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::chrono::time_point<std::chrono::high_resolution_clock> start_pack = std::chrono::high_resolution_clock::now();

    using namespace devel::scheme;
    using std::cout;
    using std::endl;

    std::vector<SearchPointWithRots> & packed_results = *packed_results_p;

    std::cout << "Building twobody tables before hack-pack" << std::endl;
    for( int ipack = 0; ipack < pd.npack; ++ipack ) {
        ScaffoldIndex si = packed_results[ipack].index.scaffold_index;
        rdd.scaffold_provider->setup_twobody_tables( si );
    }

    if ( rdd.unsat_manager ) {
        std::cout << "Building twobody tables per thread for unsats" << std::endl;
        for( int ipack = 0; ipack < pd.npack; ++ipack ) {
            ScaffoldIndex si = packed_results[ipack].index.scaffold_index;
            rdd.scaffold_provider->setup_twobody_tables_per_thread( si );
        }
    }

    print_header( "hack-packing top " + KMGT(pd.npack) );

    std::cout << "packing options: " << rdd.packopts << std::endl;
    std::cout << "packing w/rif rofts ";


    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    int64_t const out_interval = std::max<int64_t>(1,pd.npack/100);
    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for( int ipack = 0; ipack < pd.npack; ++ipack ){
        if( exception ) continue;
        try {
            if( ipack%out_interval==0 ){ cout << '*'; cout.flush(); }

            bool bad_score = packed_results[ipack].score > global_score_cut_;
            
            bool director_success = false;
            ScenePtr tscene;

            if ( ! bad_score ) {
                RifDockIndex isamp = packed_results[ipack].index;
                packed_results[ ipack ].index = isamp;
                packed_results[ ipack ].prepack_rank = ipack;
                tscene = ( rdd.scene_pt[omp_get_thread_num()] );
                director_success = rdd.director->set_scene( isamp, director_resl_, *tscene );
            }

            if ( ! director_success ) {
                packed_results[ ipack ].rotamers(); // this initializes it to blank
                packed_results[ ipack ].score = 9e9;
                continue;
            }

            std::vector<float> scores;
            packed_results[ ipack ].score = rdd.packing_objectives[rif_resl_]->score_with_rotamers( *tscene, scores, packed_results[ ipack ].rotamers() );
            packed_results[ ipack ].sasa = (uint16_t) ( scores[3] / SASA_SUBVERT_MULTIPLIER );


        } catch( std::exception const & ex ) {
            #ifdef USE_OPENMP
            #pragma omp critical
            #endif
            exception = std::current_exception();
        }
    }
    if( exception ) std::rethrow_exception(exception);
    end = std::chrono::high_resolution_clock::now();
    std::cout << std::endl;

    std::chrono::duration<double> elapsed_seconds_pack = end-start;
    std::cout << "packing rate: " << (double)pd.npack/elapsed_seconds_pack.count()                   << " iface packs per second" << std::endl;
    std::cout << "packing rate: " << (double)pd.npack/elapsed_seconds_pack.count()/omp_max_threads() << " iface packs per second per thread" << std::endl;



    std::cout << "full sort of packed samples" << std::endl;
    __gnu_parallel::sort( packed_results.begin(), packed_results.end() );

    int to_check = std::min(1000, (int)packed_results.size());
    std::cout << "Check " << to_check << " results after hackpack" << std::endl;
    for ( int i = 0; i < to_check; i++ ) {
        SearchPointWithRots const & packed_result = packed_results[i];
        if (packed_result.rotamers().size() == 0) continue;
        ScenePtr tscene( rdd.scene_pt[omp_get_thread_num()] );
        sanity_check_hackpack( rdd, packed_result.index, packed_result.rotamers_, tscene, director_resl_, rif_resl_);
    }


    std::chrono::duration<double> elapsed_seconds_all_pack = std::chrono::high_resolution_clock::now()-start_pack;
    pd.time_pck += elapsed_seconds_all_pack.count();


    return packed_results_p;
}




void
sanity_check_rots(
    RifDockData & rdd, 
    RifDockIndex i,
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers,
    ScenePtr scene,
    bool original,
    int /*director_resl*/,
    int rif_resl 
) {



    bool only_bad = true;
    bool all_missing = true;
    bool all_ala = true;

    for( int ipr = 0; ipr < rotamers->size(); ++ipr ){
        int irot = rotamers->at(ipr).second;

        BBActor bba = scene->template get_actor<BBActor>(1,rotamers->at(ipr).first);

        float rescore = rdd.rot_tgt_scorer.score_rotamer_v_target( irot, bba.position(), 10.0, 4 );
        if (rescore >= 0) {
        } else {
        }

        std::vector< std::pair< float, int > > rotscores;
        rdd.rif_ptrs[rif_resl]->get_rotamers_for_xform( bba.position(), rotscores );

        bool exists = false;
        for ( std::pair<float,int> const & p : rotscores ) {
            if (p.second == irot) {
                exists = true;
                break;
            } else {
            }
        }
        if (irot == 0) {
        } else {
            if (exists) {
                all_missing=false;
            } else {
                if (original) {
                    std::cout << "ZMISSING ROT";
                } else {
                    std::cout << "RECALCMISSING ROT";
                }
            }
            all_ala = false;
        }


    }

}

void 
sanity_check_hackpack(
    RifDockData & rdd, 
    RifDockIndex i,
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers,
    ScenePtr scene,
    int director_resl,
    int rif_resl ) {

    bool success = rdd.director->set_scene( i, director_resl, *scene );
    if ( ! success ) {
        std::cout << "Bad index" << std::endl;
        return;
    }
    sanity_check_rots(rdd, i, rotamers, scene, true, director_resl, rif_resl);

    rdd.director->set_scene( i, director_resl, *scene );
    SearchPointWithRots temp;

    rdd.packing_objectives[rif_resl]->score_with_rotamers( *scene, temp.rotamers() );

    sanity_check_rots(rdd, i, temp.rotamers_, scene, false, director_resl, rif_resl);


}



}}
