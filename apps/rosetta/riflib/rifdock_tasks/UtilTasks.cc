// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/UtilTasks.hh>

#include <riflib/types.hh>
#include <riflib/util.hh>

#include <riflib/ScoreRotamerVsTarget.hh>

#include <riflib/scaffold/ScaffoldDataCache.hh>

#include <string>
#include <vector>
#include <riflib/seeding_util.hh>

#include <boost/format.hpp>



namespace devel {
namespace scheme {


shared_ptr<std::vector<SearchPoint>> 
SortByScoreTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
SortByScoreTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
SortByScoreTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
SortByScoreTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    __gnu_parallel::sort( any_points->begin(), any_points->end() );

    return any_points;
}
    

shared_ptr<std::vector<SearchPoint>> 
FilterToBestNTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
FilterToBestNTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
FilterToBestNTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
FilterToBestNTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    typedef _AnyPointVectorsMap<AnyPoint> AnyPointVectorsMap;

    AnyPointVectorsMap map = sort_into_blocks( any_points, false, 
                                        filter_seeding_positions_separately_, 
                                        filter_scaffolds_separately_ );

    std::vector<RifDockIndex> keys;
    keys.reserve(map.size());
    for ( auto pair : map ) {
        keys.push_back(pair.first);
    }

    int num_keys = keys.size();
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for ( int i = 0; i < num_keys; i++ ) {
        std::nth_element( map[keys[i]].begin(), map[keys[i]].begin()+n_, map[keys[i]].end());
    }

    any_points->resize(0);
    any_points->reserve(map.size() * n_ );

    for ( int i = 0; i < num_keys; i++ ) {
        std::vector<AnyPoint> const & vec = map[keys[i]];
        int count = std::min<int>( n_, vec.size());
        for ( int j = 0; j < count; j++ ) {
            any_points->push_back( vec[j] );
        }
    }

    return any_points;
}
    


shared_ptr<std::vector<SearchPoint>> 
FilterByScoreCutTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
FilterByScoreCutTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
FilterByScoreCutTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
FilterByScoreCutTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    shared_ptr<std::vector<AnyPoint>> out = make_shared<std::vector<AnyPoint>>();
    out->reserve( any_points->size());

    for ( AnyPoint const & pt : *any_points ) {
        if ( pt.score > score_cut_ ) continue;
        out->push_back( pt );
    }

    any_points->resize(0);
    return out;
}
    



shared_ptr<std::vector<SearchPoint>> 
FilterByFracTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
FilterByFracTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
FilterByFracTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
FilterByFracTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    typedef _AnyPointVectorsMap<AnyPoint> AnyPointVectorsMap;

    AnyPointVectorsMap map = sort_into_blocks( any_points, false, 
                                        filter_seeding_positions_separately_, 
                                        filter_scaffolds_separately_ );

    std::vector<RifDockIndex> keys;
    keys.reserve(map.size());
    for ( auto pair : map ) {
        keys.push_back(pair.first);
    }

    int num_keys = keys.size();
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for ( int i = 0; i < num_keys; i++ ) {
        std::sort( map[keys[i]].begin(), map[keys[i]].end());
    }

    int trial_size = any_points->size() * frac_;
    any_points->resize(0);
    any_points->reserve(trial_size );

    for ( int i = 0; i < num_keys; i++ ) {
        std::vector<AnyPoint> const & vec = map[keys[i]];
        int count = int(std::ceil(vec.size() * frac_));
        count = std::max<int>( keep_at_least_, count );
        count = std::min<int>( count, vec.size() );
        for ( int j = 0; j < count; j++ ) {
            any_points->push_back( vec[j] );
        }
    }

    return any_points;
}
    
shared_ptr<std::vector<SearchPoint>> 
FilterByBiggestBlocksFracTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
FilterByBiggestBlocksFracTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
FilterByBiggestBlocksFracTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
FilterByBiggestBlocksFracTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    typedef _AnyPointVectorsMap<AnyPoint> AnyPointVectorsMap;

    AnyPointVectorsMap map = sort_into_blocks( any_points, false, 
                                        filter_seeding_positions_separately_, 
                                        filter_scaffolds_separately_ );

    if ( print_seeds_ ) {
        std::cout << "Cluster score of each seeding pos:" << std::endl;
        int seeding_size = rdd.director->size(0, RifDockIndex()).seeding_index;
        for ( int i = 0; i < seeding_size; i++ ) {
            int size = 0;
            RifDockIndex rdi = RifDockIndex(0, i, ScaffoldIndex());
            if (map.count(rdi) > 0) size = map[rdi].size();
            std::cout << i << ":" << size << " ";
        }
        std::cout << std::endl << std::endl;
    }

    std::vector<RifDockIndex> keys;
    keys.reserve(map.size());
    std::vector<size_t> sizes;
    sizes.reserve(map.size());
    for ( auto pair : map ) {
        keys.push_back(pair.first);
        sizes.push_back(pair.second.size());
    }

    __gnu_parallel::sort( sizes.begin(), sizes.end() );

    size_t last_index = int( std::floor ( (1 - frac_) * sizes.size() ) );
    size_t cut_size = sizes.at(last_index);

    int trial_size = any_points->size() * frac_;
    any_points->resize(0);
    any_points->reserve(trial_size );

    int num_keys = keys.size();
    for ( int i = 0; i < num_keys; i++ ) {
        std::vector<AnyPoint> const & vec = map[keys[i]];
        if (vec.size() < cut_size) continue;

        int count = vec.size();
        for ( int j = 0; j < count; j++ ) {
            any_points->push_back( vec[j] );
        }
    }

 
    return any_points;
}
    
shared_ptr<std::vector<SearchPoint>> 
RemoveRedundantPointsTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
RemoveRedundantPointsTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
RemoveRedundantPointsTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
RemoveRedundantPointsTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    typedef _AnyPointVectorsMap<AnyPoint> AnyPointVectorsMap;

    AnyPointVectorsMap map = sort_into_blocks( any_points, false, 
                                        filter_seeding_positions_separately_, 
                                        filter_scaffolds_separately_ );

    std::vector<RifDockIndex> keys;
    keys.reserve(map.size());
    for ( auto pair : map ) {
        keys.push_back(pair.first);
    }

    int trial_size = any_points->size();
    any_points->resize(0);
    any_points->reserve(trial_size );

    float redundancy_filter_rg = rdd.scaffold_provider->get_data_cache_slow(ScaffoldIndex())
                                        ->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );

    // for ( RifDockIndex const & rdi : keys ) {

    int seeding_size = rdd.director->size(0, RifDockIndex()).seeding_index;
    for ( int seed = 0; seed < seeding_size; seed++ ) {
        int size = 0;
        RifDockIndex rdi = RifDockIndex(0, seed, ScaffoldIndex());
        if (map.count(rdi) == 0) continue;
        std::vector<AnyPoint> const & vec = map[rdi];
        runtime_assert( vec.size() > 0 );

        int first_one = any_points->size();
        any_points->push_back(vec[0]);
        int end = any_points->size();


        for ( int i = 1; i < vec.size(); i++ ) {
            rdd.director->set_scene( vec[i].index, director_resl_, *rdd.scene_minimal );
            EigenXform p1 = rdd.scene_minimal->position(1);

            bool is_redundant = false;
            for ( int j = first_one; j < end; j++ ) {
                rdd.director->set_scene( any_points->at(j).index, director_resl_, *rdd.scene_minimal );
                EigenXform p2 = rdd.scene_minimal->position(1);

                float mag = devel::scheme::xform_magnitude( p2 * p1.inverse(), redundancy_filter_rg );
                if ( mag <= redundancy_mag_ ) {
                    is_redundant = true;
                    break;
                }
            }
            if ( ! is_redundant ) {
                any_points->push_back( vec[i] );
                end = any_points->size();
            }

        }


    }


    std::cout << "Number of total searching points left: " << any_points->size() << std::endl << std::endl;

    return any_points;
}
    

shared_ptr<std::vector<SearchPoint>> 
DumpSeedingClusterScoreTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
DumpSeedingClusterScoreTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
DumpSeedingClusterScoreTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
DumpSeedingClusterScoreTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    typedef _AnyPointVectorsMap<AnyPoint> AnyPointVectorsMap;

    AnyPointVectorsMap map = sort_into_blocks( any_points, false, 
                                        true, 
                                        true );



    std::cout << "Cluster score of each seeding pos:" << std::endl;
    int seeding_size = rdd.director->size(0, RifDockIndex()).seeding_index;
    for ( int i = 0; i < seeding_size; i++ ) {
        int size = 0;
        RifDockIndex rdi = RifDockIndex(0, i, ScaffoldIndex());
        if (map.count(rdi) > 0) size = map[rdi].size();
        if (size > 0 ) std::cout << i << ":" << size << " ";
    }
    std::cout << std::endl << std::endl;

    return any_points;
}
    


shared_ptr<std::vector<SearchPoint>> 
DumpIndividualScoresTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
DumpIndividualScoresTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
DumpIndividualScoresTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
DumpIndividualScoresTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    typedef _AnyPointVectorsMap<AnyPoint> AnyPointVectorsMap;

    AnyPointVectorsMap map = sort_into_blocks( any_points, false, 
                                        true, 
                                        true );



    std::cout << "Cluster score of each seeding pos:" << std::endl;
    int seeding_size = rdd.director->size(0, RifDockIndex()).seeding_index;
    for ( int i = 0; i < seeding_size; i++ ) {
        int size = 0;
        RifDockIndex rdi = RifDockIndex(0, i, ScaffoldIndex());
        if (map.count(rdi) == 0) continue;
        std::vector<AnyPoint> const & vec = map[rdi];

        for ( int j = 0; j < vec.size(); j++ ) {
            std::cout << i << ":" << j << " " << vec[0].score << " " << vec[0].index.nest_index << std::endl;
        }
    }


    std::cout << std::endl << std::endl;

    return any_points;
}
    




shared_ptr<std::vector<SearchPoint>> 
DumpScoresTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
DumpScoresTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
DumpScoresTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
DumpScoresTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::ofstream out;
    out.open( file_name_ );

    out << "RifDockIndex\tscore" << std::endl;
    for ( AnyPoint const & pt : *any_points ) {
        out << pt.index << "\t" << pt.score << std::endl;
    }

    out.close();

    return any_points;
}
    
shared_ptr<std::vector<SearchPointWithRots>>
DumpRotScoresTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::ofstream out;
    out.open( file_name_ );

    if ( use_rif_ ) {
        utility_exit_with_message("Not implemented!!!");
    } else {

        devel::scheme::ScoreRotamerVsTarget<
            VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
        > rot_tgt_scorer;
        rot_tgt_scorer.rot_index_p_ = rdd.rot_index_p;
        rot_tgt_scorer.target_field_by_atype_ = rdd.target_field_by_atype;
        rot_tgt_scorer.target_donors_ = *rdd.target_donors;
        rot_tgt_scorer.target_acceptors_ = *rdd.target_acceptors;
        rot_tgt_scorer.hbond_weight_ = rdd.packopts.hbond_weight;
        rot_tgt_scorer.upweight_iface_ = rdd.packopts.upweight_iface;
        rot_tgt_scorer.upweight_multi_hbond_ = rdd.packopts.upweight_multi_hbond;
#ifdef USEGRIDSCORE
        rot_tgt_scorer.grid_scorer_ = rdd.grid_scorer;
        rot_tgt_scorer.soft_grid_energies_ = rdd.opt.soft_rosetta_grid_energies;
#endif

        out << "RifDockIndex\ttotal";


        bool first_line = true;
        for ( SearchPointWithRots const & sp : *search_point_with_rotss ) {
            rdd.director->set_scene( sp.index, resl_, *rdd.scene_minimal );

            int num_actors = rdd.scene_minimal->template num_actors<BBActor>(1);

            if (first_line) {
                out << "RifDockIndex";
                for( int i_actor = 0; i_actor < num_actors; ++i_actor ) {
                    out << boost::str(boost::format("\tscore%i\t1body%i")%i_actor%i_actor);
                }
                out << std::endl;
            }
            first_line = false;



            ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(sp.index.scaffold_index);
            runtime_assert( sdc );
            std::vector<std::vector<float> > const * rotamer_energies_1b = sdc->local_onebody_p.get();


            out << sp.index;

            std::vector<int> rot_at_res(num_actors, -1);
            for( int ipr = 0; ipr < sp.numrots(); ++ipr ) {
                int ires = sp.rotamers().at(ipr).first;
                int irot = sp.rotamers().at(ipr).second;

                rot_at_res[ires] = irot;
            }

            out << boost::str(boost::format("\t%.3f")%sp.score);

            for( int i_actor = 0; i_actor < num_actors; ++i_actor ) {
                
                int irot = rot_at_res[i_actor];
                if (irot == -1) {
                    out << "\tNaN\tNaN";
                    continue;
                }

                BBActor bba = rdd.scene_minimal->template get_actor<BBActor>(1,i_actor);
                runtime_assert( bba.index_ == i_actor );

                float const onebody = (*rotamer_energies_1b).at(i_actor).at(irot);
                float score = rot_tgt_scorer.score_rotamer_v_target( irot, bba.position(), 10.0, 4 );

                out << boost::str(boost::format("\t%.3f\t%.3f")%score%onebody);

            }

            out << std::endl;


        }


    }

    out.close();

    return search_point_with_rotss;


}


}}
