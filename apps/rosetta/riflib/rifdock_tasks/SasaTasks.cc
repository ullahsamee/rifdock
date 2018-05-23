// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/rifdock_tasks/SasaTasks.hh>

#include <riflib/types.hh>


#include <string>
#include <vector>



namespace devel {
namespace scheme {


shared_ptr<std::vector<SearchPoint>> 
FilterBySasaTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
FilterBySasaTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
FilterBySasaTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
FilterBySasaTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::cout << "Sasa cut: " << sasa_cut_ << std::endl;

    shared_ptr<std::vector<AnyPoint>> out = make_shared<std::vector<AnyPoint>>();
    out->reserve( any_points->size());

    for ( AnyPoint const & pt : *any_points ) {
        if ( (float)pt.sasa < sasa_cut_ ) continue;
        out->push_back( pt );
    }

    any_points->resize(0);
    return out;

}
    


shared_ptr<std::vector<SearchPoint>> 
FilterByScorePer1000SasaTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
FilterByScorePer1000SasaTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
FilterByScorePer1000SasaTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}

template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
FilterByScorePer1000SasaTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::cout << "ScorePer1000Sasa cut: " << score_per_1000_sasa_cut_ << std::endl;

    shared_ptr<std::vector<AnyPoint>> out = make_shared<std::vector<AnyPoint>>();
    out->reserve( any_points->size());

    for ( AnyPoint const & pt : *any_points ) {
        if ( pt.score * 1000.0f / (float)pt.sasa > score_per_1000_sasa_cut_ ) continue;
        out->push_back( pt );
    }

    any_points->resize(0);
    return out;

}
    


shared_ptr<std::vector<SearchPoint>> 
SortByScorePer1000SasaTask::return_search_points( 
    shared_ptr<std::vector<SearchPoint>> search_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {
    return return_any_points( search_points, rdd, pd );
}
shared_ptr<std::vector<SearchPointWithRots>> 
SortByScorePer1000SasaTask::return_search_point_with_rotss( 
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( search_point_with_rotss, rdd, pd );
}
shared_ptr<std::vector<RifDockResult>> 
SortByScorePer1000SasaTask::return_rif_dock_results( 
    shared_ptr<std::vector<RifDockResult>> rif_dock_results, 
    RifDockData & rdd, 
    ProtocolData & pd ) { 
    return return_any_points( rif_dock_results, rdd, pd );
}



template<class AnyPoint>
shared_ptr<std::vector<AnyPoint>>
SortByScorePer1000SasaTask::return_any_points( 
    shared_ptr<std::vector<AnyPoint>> any_points, 
    RifDockData & rdd, 
    ProtocolData & pd ) {

    std::cout << "SortByScorePer1000SasaTask " << std::endl;

    ScorePer1000SasaComparator comparator;

    __gnu_parallel::sort( any_points->begin(), any_points->end(), comparator );


    return any_points;

}
    

}}
