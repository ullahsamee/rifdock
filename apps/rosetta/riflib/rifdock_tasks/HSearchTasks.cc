// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/HSearchTasks.hh>

#include <riflib/types.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>
#include <riflib/rifdock_tasks/OutputResultsTasks.hh>


#include <string>
#include <vector>
#include <unordered_map>


#include <ObjexxFCL/format.hh>



namespace devel {
namespace scheme {


shared_ptr<std::vector<SearchPoint>> 
DiversifyByNestTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
DiversifyByNestTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
DiversifyByNestTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
DiversifyByNestTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {


    uint64_t nest_size = rdd.director->size(resl_, RifDockIndex()).nest_index;

    shared_ptr<std::vector<AnyPoint>> diversified = make_shared<std::vector<AnyPoint>>( nest_size * any_points->size() );

    uint64_t added = 0;
    for ( AnyPoint const & pt : *any_points ) {

        for ( uint64_t i = 0; i < nest_size; i++ ) {
            (*diversified)[added] = pt;
            (*diversified)[added].index.nest_index = i;
            added++;
        }
    }

    any_points->clear();

    return diversified;
}



shared_ptr<std::vector<SearchPoint>> 
HSearchInit::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {


    pd.start_rif = std::chrono::high_resolution_clock::now();

    SelectiveRifDockIndexHasher   hasher( false, false, true );
    SelectiveRifDockIndexEquater equater( false, false, true );
    std::unordered_map<RifDockIndex, bool, SelectiveRifDockIndexHasher, SelectiveRifDockIndexEquater> uniq_scaffolds (
        1000, hasher, equater);

    for ( SearchPoint const & sp : *search_points ) {
        uniq_scaffolds[ sp.index ] = true;
    }

    pd.unique_scaffolds.clear();
    for ( std::pair<RifDockIndex, bool> pair : uniq_scaffolds ) {
        ScaffoldIndex si = pair.first.scaffold_index;
        rdd.scaffold_provider->get_data_cache_slow(si)->setup_onebody_tables( rdd.rot_index_p, rdd.opt);
        pd.unique_scaffolds.push_back(si);
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    print_header( "perform hierarchical search" ); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////


    if (rdd.opt.multiply_beam_by_scaffolds) pd.beam_multiplier *= pd.unique_scaffolds.size();
    if (rdd.opt.multiply_beam_by_seeding_positions) pd.beam_multiplier *= rdd.director->size(0, RifDockIndex()).seeding_index;

    if ( pd.beam_multiplier > rdd.opt.max_beam_multiplier ) {
        std::cout << "Beam multiplier exceeds max_beam_multiplier: " << pd.beam_multiplier << std::endl;
        pd.beam_multiplier = rdd.opt.max_beam_multiplier;
    }

    std::cout << "Beam size multiplier: " << pd.beam_multiplier << std::endl;

    return search_points;
}


shared_ptr<std::vector<SearchPoint>> 
HSearchScoreAtReslTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points_p, 
    RifDockData & rdd, 
    ProtocolData & pd ) {


    using ObjexxFCL::format::F;
    using ObjexxFCL::format::I;
    using std::cout;
    using std::endl;

    std::vector<SearchPoint> & search_points = *search_points_p;

    bool using_csts = false;
    for ( ScaffoldIndex si : pd.unique_scaffolds ) {
        ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(si);
        using_csts |= sdc->prepare_contraints( rdd.target, rdd.RESLS[rif_resl_] );
    }

    bool need_sdc = using_csts || tether_to_input_position_cut_ != 0;


    cout << "HSearsh stage " << rif_resl_+1 << " resl " << F(5,2,rdd.RESLS[rif_resl_]) << " begin threaded sampling, " << KMGT(search_points.size()) << " samples: ";
    int64_t const out_interval = search_points.size()/50;
    std::exception_ptr exception = nullptr;
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();
    pd.total_search_effort += search_points.size();

    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for( int64_t i = 0; i < search_points.size(); ++i ){
        if( exception ) continue;
        try {
            if( i%out_interval==0 ){ cout << '*'; cout.flush(); }
            RifDockIndex const isamp = search_points[i].index;

            ScenePtr tscene( rdd.scene_pt[omp_get_thread_num()] );
            bool director_success = rdd.director->set_scene( isamp, director_resl_, *tscene );
            if ( ! director_success ) {
                search_points[i].score = 9e9;
                continue;
            }

            if ( need_sdc ) {
                ScaffoldIndex si = isamp.scaffold_index;
                ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(si);

                if( tether_to_input_position_cut_ > 0 ){
                    float redundancy_filter_rg = sdc->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );

                    EigenXform x;// = tscene->position(1);
                    rdd.nest.get_state( isamp.nest_index, director_resl_, x );
                    x.translation() -= sdc->scaffold_center;
                    float xmag =  xform_magnitude( x, redundancy_filter_rg );
                    if( xmag > tether_to_input_position_cut_ + rdd.RESLS[rif_resl_] ){
                        search_points[i].score = 9e9;
                        continue;
                    } 
                }

                /////////////////////////////////////////////////////
                /////// Longxing' code  ////////////////////////////
                ////////////////////////////////////////////////////
                if (using_csts) {
                    EigenXform x = tscene->position(1);
                    bool pass_all = true;
                    for(CstBaseOP p : sdc->csts) {
                        if (!p->apply( x )) {
                            pass_all = false;
                            break;
                        }
                    }
                    if (!pass_all) {
                        search_points[i].score = 9e9;
                        continue;
                    }
                }
            }

            // the real rif score!!!!!!
            search_points[i].score = rdd.objectives[rif_resl_]->score( *tscene );// + tot_sym_score;


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
    pd.hsearch_rate = (double)search_points.size()/ elapsed_seconds_rif.count()/omp_max_threads();
    cout << endl;// << "done threaded sampling, partitioning data..." << endl;


    return search_points_p;
}


shared_ptr<std::vector<SearchPoint>> 
HSearchFilterSortTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points_p, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    using ObjexxFCL::format::F;
    using ObjexxFCL::format::I;

    std::vector<SearchPoint> & search_points = *search_points_p;

    SearchPoint max_pt, min_pt;
    int64_t len = search_points.size();
    uint64_t keeping = num_to_keep_ * pd.beam_multiplier;
    if( search_points.size() > keeping ){
        __gnu_parallel::nth_element( search_points.begin(), search_points.begin()+ keeping, search_points.end() );
        len = keeping;
        min_pt = *__gnu_parallel::min_element( search_points.begin(), search_points.begin()+len );
        max_pt = *(search_points.begin()+keeping);
    } else {
        min_pt = *__gnu_parallel::min_element( search_points.begin(), search_points.end() );
        max_pt = *__gnu_parallel::max_element( search_points.begin(), search_points.end() );
    }

    std::cout << "HSearsh stage " << resl_+1 << " complete, resl. " << F(7,3,rdd.RESLS[resl_]) << ", "
          << " " << KMGT(search_points.size()) << ", promote: " << F(9,6,min_pt.score) << " to "
          << F(9,6, std::min(global_score_cut_,max_pt.score)) << " rate " << KMGT(pd.hsearch_rate) << "/s/t " << std::endl;

    if ( prune_extra_ ) {
        search_points.resize(len);
    }

    return search_points_p;
}

shared_ptr<std::vector<SearchPoint>> 
HSearchScaleToReslTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points_p, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::vector<SearchPoint> & search_points = *search_points_p;

    shared_ptr<std::vector<SearchPoint>> out_points_p = make_shared<std::vector<SearchPoint>>( );
    std::vector<SearchPoint> & out_points = *out_points_p;

    // this is the normal path
    if ( target_resl_ >= current_resl_ ) {
        int num_resls = target_resl_ - current_resl_;

        int use_pow2 = 1;
        for ( int i = 0; i < num_resls; i++) {
            use_pow2 *= DIMPOW2_;
        }

        __gnu_parallel::sort( search_points.begin(), search_points.end() );
        size_t good_points = 0;
        for ( good_points = 0; good_points < search_points.size(); good_points++ ) {
            if ( search_points[good_points].score >= global_score_cut_ ) break;
        }

        if( current_resl_ == 0 ) pd.non0_space_size += good_points;

        out_points.resize( use_pow2 * good_points );

        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic,64)
        #endif
        for( int64_t i = 0; i < good_points; ++i ){
            RifDockIndex rdi0 = search_points[i].index;
            uint64_t isamp0 = use_pow2 * rdi0.nest_index;
            uint64_t array_offset = i*use_pow2;

            for( uint64_t j = 0; j < use_pow2; ++j ){
                SearchPoint & rdi = out_points[array_offset+j];
                rdi = rdi0;
                rdi.index.nest_index = isamp0 + j;
            }
        }

    } else {

        int dropped_resls = current_resl_ - target_resl_;
        int shift_factor = dropped_resls * 6;

        SelectiveRifDockIndexHasher   hasher( true, true, true );
        SelectiveRifDockIndexEquater equater( true, true, true );
        std::unordered_map<RifDockIndex, bool, SelectiveRifDockIndexHasher, SelectiveRifDockIndexEquater> uniq_positions(
            search_points.size()/100 + 1000, hasher, equater );

        for ( SearchPoint const & sp : search_points ) {
            RifDockIndex rdi = sp.index;
            rdi.nest_index = rdi.nest_index >> shift_factor;
            uniq_positions[rdi] = true;
        }

        out_points.resize(uniq_positions.size());

        int added = 0;
        for ( std::pair<RifDockIndex, bool> pair : uniq_positions ) {
            out_points[added++] = pair.first;
        }

        std::cout << "Num uniq positions: " << out_points.size() << std::endl;

    }


    search_points.clear();

    return out_points_p;

}

shared_ptr<std::vector<SearchPoint>> 
HSearchFinishTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points_p, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::vector<SearchPoint> & search_points = *search_points_p;

    std::cout << "full sort of final samples" << std::endl;
    __gnu_parallel::sort( search_points.begin(), search_points.end() );

    size_t good_points = 0;
    for( good_points; good_points < search_points.size(); ++good_points ){
        if( search_points[good_points].score > 0 ) break;
    }

    search_points.resize(good_points);

    pd.beam_multiplier = 1;

    std::cout << "total non-0 space size was approx " << float(pd.non0_space_size)*1024.0*1024.0*1024.0 << " grid points" << std::endl;
    std::cout << "total search effort " << KMGT(pd.total_search_effort) << std::endl;


    std::chrono::duration<double> elapsed_seconds_rif = std::chrono::high_resolution_clock::now()-pd.start_rif;
    pd.time_rif += elapsed_seconds_rif.count();

    return search_points_p;
}

shared_ptr<std::vector<SearchPoint>> 
DumpHSearchFramesTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points_p, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::cout << "Dumping frames with prefix: " << prefix_ << std::endl;

    uint64_t dump_every = std::floor( search_points_p->size() / dump_x_frames_per_resl_ );

    if ( dump_only_best_frames_ ) {
        dump_every = std::max( 1, dump_only_best_stride_ );
        __gnu_parallel::sort( search_points_p->begin(), search_points_p->end() );
    }

    for ( uint64_t i = 0; i < search_points_p->size(); i++ ) {
        if ( (   dump_only_best_frames_ && i < dump_x_frames_per_resl_ ) ||
             ( ! dump_only_best_frames_ && ( i % dump_every ) == 0 )) {
            std::string filename = prefix_ + boost::str( boost::format( "_%06i.pdb.gz" ) % (i/dump_every));
            dump_search_point_( rdd, (*search_points_p)[i], filename, director_resl_, rif_resl_, true );
        }
    }


    return search_points_p;
}



}}
