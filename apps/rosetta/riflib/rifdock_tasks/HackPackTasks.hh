// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_rifdock_tasks_HackPackTasks_hh
#define INCLUDED_riflib_rifdock_tasks_HackPackTasks_hh

#include <riflib/types.hh>
#include <riflib/task/SearchPointWithRotsTask.hh>
#include <riflib/task/AnyPointTask.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {

struct FilterForHackPackTask : public AnyPointTask {

    FilterForHackPackTask(
        float hack_pack_frac,
        int pack_n_iters,
        float pack_iter_mult
        ) :
        hack_pack_frac_( hack_pack_frac ),
        pack_n_iters_( pack_n_iters ),
        pack_iter_mult_( pack_iter_mult )
        {}

    shared_ptr<std::vector<SearchPoint>> 
    return_search_points( 
        shared_ptr<std::vector<SearchPoint>> search_points, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

    shared_ptr<std::vector<SearchPointWithRots>> 
    return_search_point_with_rotss( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;

    shared_ptr<std::vector<RifDockResult>> 
    return_rif_dock_results( 
        shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;


private:
    template<class AnyPoint>
    shared_ptr<std::vector<AnyPoint>>
    return_any_points( 
        shared_ptr<std::vector<AnyPoint>> any_points, 
        RifDockData & rdd, 
        ProtocolData & pd ); // override

private:
    float hack_pack_frac_;
    int pack_n_iters_;
    float pack_iter_mult_;

};


struct HackPackTask : public SearchPointWithRotsTask {

    HackPackTask( int resl, float global_score_cut ) : 
        resl_(resl),
        global_score_cut_(global_score_cut)
    {}

    shared_ptr<std::vector<SearchPointWithRots>>
    return_search_point_with_rotss( 
        shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
        RifDockData & rdd, 
        ProtocolData & pd ) override;


private:
    int resl_;
    float global_score_cut_;


};



void
sanity_check_rots(
    RifDockData & rdd, 
    RifDockIndex i,
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers,
    ScenePtr scene,
    bool original,
    int resl 
); 


void 
sanity_check_hackpack(
    RifDockData & rdd, 
    RifDockIndex i,
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers,
    ScenePtr scene,
    int resl 
);




}}

#endif
