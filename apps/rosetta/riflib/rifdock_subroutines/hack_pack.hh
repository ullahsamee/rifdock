// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_hack_pack_hh
#define INCLUDED_riflib_rifdock_subroutines_hack_pack_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


namespace devel {
namespace scheme {


void
hack_pack(
    shared_ptr< std::vector< SearchPointWithRots > > & hsearch_results_p,
    std::vector< SearchPointWithRots > & packed_results,
    RifDockData & rdd,
    int64_t total_search_effort, int64_t & npack) {


    using namespace devel::scheme;
    using std::cout;
    using std::endl;



    std::vector< SearchPointWithRots > & hsearch_results = *hsearch_results_p;


    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    size_t n_packsamp = 0;
    for( n_packsamp; n_packsamp < hsearch_results.size(); ++n_packsamp ){
        if( hsearch_results[n_packsamp].score > 0 ) break;
    }
    int const config = rdd.RESLS.size()-1;
    npack = std::min( n_packsamp, (size_t)(total_search_effort *
        ( rdd.opt.hack_pack_frac / (rdd.packopts.pack_n_iters*rdd.packopts.pack_iter_mult)) ) );

    packed_results.resize( npack );


    print_header( "hack-packing top " + KMGT(npack) );

    std::cout << "Building twobody tables before hack-pack" << std::endl;
    for( int ipack = 0; ipack < npack; ++ipack ) {
        ScaffoldIndex si = hsearch_results[ipack].index.scaffold_index;
        rdd.scaffold_provider->setup_twobody_tables( si );
    }


    std::cout << "packing options: " << rdd.packopts << std::endl;
    std::cout << "packing w/rif rofts ";
    int64_t const out_interval = std::max<int64_t>(1,npack/100);
    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for( int ipack = 0; ipack < npack; ++ipack ){
        if( exception ) continue;
        try {
            if( ipack%out_interval==0 ){ cout << '*'; cout.flush(); }
            RifDockIndex isamp = hsearch_results[ipack].index;
            if( hsearch_results[ipack].score > rdd.opt.global_score_cut ) continue;
            packed_results[ ipack ].index = isamp;
            packed_results[ ipack ].prepack_rank = ipack;
            ScenePtr tscene( rdd.scene_pt[omp_get_thread_num()] );
            rdd.director->set_scene( isamp, rdd.RESLS.size()-1, *tscene );
            packed_results[ ipack ].score = rdd.packing_objective->score_with_rotamers( *tscene, packed_results[ ipack ].rotamers() );

            sanity_check_hackpack( rdd, isamp, packed_results[ ipack ].rotamers_, tscene);
        } catch( std::exception const & ex ) {
            #ifdef USE_OPENMP
            #pragma omp critical
            #endif
            exception = std::current_exception();
        }
    }
    if( exception ) std::rethrow_exception(exception);
    end = std::chrono::high_resolution_clock::now();

        std::cout << "Check after hackpack111" << std::endl;
    for (SearchPointWithRots const & packed_result : packed_results) {
        ScenePtr tscene( rdd.scene_pt[omp_get_thread_num()] );
        rdd.director->set_scene( packed_result.index, rdd.RESLS.size()-1, *tscene );
        sanity_check_hackpack( rdd, packed_result.index, packed_result.rotamers_, tscene);
    }



    // std::vector< RifDockIndex > rdis(packed_results.size());
    // std::vector< std::vector< std::pair<intRot,intRot> >> rotamerss(packed_results.size());
    // for (int i = 0; i < packed_results.size(); i++) {
    //     rdis[i] = packed_results[i].index;
    //     rotamerss[i] = *(packed_results[i].rotamers_);
    //     // old_copies[i].rotamers_ = make_shared<std::vector< std::pair<intRot,intRot> >>(*old_copies[i].rotamers_);
    // }






    // std::vector< SearchPointWithRots > old_copies(packed_results.size());
    // for (int i = 0; i < old_copies.size(); i++) {
    //     old_copies[i] = packed_results[i];
    //     old_copies[i].rotamers_ = make_shared<std::vector< std::pair<intRot,intRot> >>(*old_copies[i].rotamers_);
    // }



    std::cout << std::endl;
    std::cout << "full sort of packed samples" << std::endl;
    // __gnu_parallel::sort( packed_results.begin(), packed_results.end() );


    
    // std::sort( packed_results.begin(), packed_results.end() );



    // for ( uint64_t i = 0; i < packed_results.size() - 10000; i++ ) {
    //     // packed_results[i] = packed_results[i+1];
    //     swap(packed_results[i], packed_results[i+10000]);
    // }

    // std::cout << "Checking old copies" << std::endl;
    // for (int i = 0; i < old_copies.size(); i++) {
    //     if (i%100 == 0) {
    //         std::cout << i << std::endl;
    //     }
    //     RifDockIndex index = old_copies[i].index;

    //     bool found_it = false;
    //     for (SearchPointWithRots const & sp : packed_results) {
    //         if (!(index == sp.index)) {
    //             found_it = true;
    //             if (*(sp.rotamers_) != *(old_copies[i].rotamers_)) {
    //                 std::cout << "DIFF " << std::endl;
    //             }
    //         }
    //     }
    //     if (!found_it) {
    //         std::cout << "NOTFOUND " << std::endl;
    //     }
    // }







    //     std::cout << "Checking old copies" << std::endl;
    // for (int i = 0; i < rdis.size(); i++) {
    //     if (i%100 == 0) {
    //         std::cout << i << std::endl;
    //     }
    //     RifDockIndex index = rdis[i];

    //     bool found_it = false;
    //     for (SearchPointWithRots const & sp : packed_results) {
    //         if ((index == sp.index)) {
    //             found_it = true;
    //             if (*(sp.rotamers_) != rotamerss[i]) {
    //                 std::cout << "DIFF " << std::endl;
    //             }
    //         }
    //     }
    //     if (!found_it) {
    //         std::cout << "NOTFOUND " << std::endl;
    //     }
    // }











    std::chrono::duration<double> elapsed_seconds_pack = end-start;
    std::cout << "packing rate: " << (double)npack/elapsed_seconds_pack.count()                   << " iface packs per second" << std::endl;
    std::cout << "packing rate: " << (double)npack/elapsed_seconds_pack.count()/omp_max_threads() << " iface packs per second per thread" << std::endl;



    //         std::cout << "Check after hackpack" << std::endl;
    // for (int i = 0; i < packed_results.size(); i++) {
    //     if (i%100 == 0) {
    //         std::cout << i << std::endl;
    //     }
    //     RifDockIndex rdi = rdis[i];
    //     shared_ptr<std::vector< std::pair<intRot,intRot> >> rotamers = make_shared<std::vector< std::pair<intRot,intRot> >>(rotamerss[i]);

    //     ScenePtr tscene( rdd.scene_pt[omp_get_thread_num()] );
    //     rdd.director->set_scene( rdi, rdd.RESLS.size()-1, *tscene );
    //     sanity_check_hackpack( rdd, rdi, rotamers, tscene);
    // }


}

}}


#endif