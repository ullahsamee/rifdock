// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/HackPackTask.hh>

#include <riflib/types.hh>


#include <string>
#include <vector>



namespace devel {
namespace scheme {

shared_ptr<std::vector<SearchPointWithRots>>
HackPackTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> packed_results_p, 
    RifDockData & rdd, 
    ProtocolData & pd ) {


    using namespace devel::scheme;
    using std::cout;
    using std::endl;

    std::vector<SearchPointWithRots> & packed_results = *packed_results_p;

    std::cout << "Building twobody tables before hack-pack" << std::endl;
    for( int ipack = 0; ipack < pd.npack; ++ipack ) {
        ScaffoldIndex si = packed_results[ipack].index.scaffold_index;
        rdd.scaffold_provider->setup_twobody_tables( si );
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
                director_success = rdd.director->set_scene( isamp, resl_, *tscene );
            }

            if ( ! director_success ) {
                packed_results[ ipack ].rotamers(); // this initializes it to blank
                packed_results[ ipack ].score = 9e9;
                continue;
            }

            packed_results[ ipack ].score = rdd.packing_objectives[resl_]->score_with_rotamers( *tscene, packed_results[ ipack ].rotamers() );

        } catch( std::exception const & ex ) {
            #ifdef USE_OPENMP
            #pragma omp critical
            #endif
            exception = std::current_exception();
        }
    }
    if( exception ) std::rethrow_exception(exception);
    end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_seconds_pack = end-start;
    std::cout << "packing rate: " << (double)pd.npack/elapsed_seconds_pack.count()                   << " iface packs per second" << std::endl;
    std::cout << "packing rate: " << (double)pd.npack/elapsed_seconds_pack.count()/omp_max_threads() << " iface packs per second per thread" << std::endl;



    std::cout << std::endl;
    std::cout << "full sort of packed samples" << std::endl;
    __gnu_parallel::sort( packed_results.begin(), packed_results.end() );

    int to_check = std::min(1000, (int)packed_results.size());
    std::cout << "Check " << to_check << " results after hackpack" << std::endl;
    for ( int i = 0; i < to_check; i++ ) {
        SearchPointWithRots const & packed_result = packed_results[i];
        if (packed_result.rotamers().size() == 0) continue;
        ScenePtr tscene( rdd.scene_pt[omp_get_thread_num()] );
        sanity_check_hackpack( rdd, packed_result.index, packed_result.rotamers_, tscene, resl_);
    }



    return packed_results_p;
}



}}
