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

#include <riflib/scaffold/ScaffoldDataCache.hh>

#include <string>
#include <vector>

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
