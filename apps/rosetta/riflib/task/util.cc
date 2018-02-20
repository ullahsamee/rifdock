// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:



#include <riflib/types.hh>

#include <riflib/task/util.hh>



namespace devel {
namespace scheme {


shared_ptr<std::vector<SearchPointWithRots>>
search_point_with_rotss_from_search_points(shared_ptr<std::vector<SearchPoint>> search_points) {

    size_t length = search_points->size();
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss = make_shared<std::vector<SearchPointWithRots>>( length );

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        (*search_point_with_rotss)[i] = (*search_points)[i];
    }

    return search_point_with_rotss;
}
shared_ptr<std::vector<RifDockResult>>
rif_dock_results_from_search_points(shared_ptr<std::vector<SearchPoint>> search_points) {

    size_t length = search_points->size();
    shared_ptr<std::vector<RifDockResult>> rif_dock_results = make_shared<std::vector<RifDockResult>>( length);

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        (*rif_dock_results)[i] = (*search_points)[i];
    }

    return rif_dock_results;
}


shared_ptr<std::vector<SearchPoint>>
search_points_from_search_point_with_rotss(shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss) {

    size_t length = search_point_with_rotss->size();
    shared_ptr<std::vector<SearchPoint>> search_points = make_shared<std::vector<SearchPoint>>( length);

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        (*search_points)[i] = (*search_point_with_rotss)[i];
    }

    return search_points;
}
shared_ptr<std::vector<RifDockResult>>
rif_dock_results_from_search_point_with_rotss(shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss) {

    size_t length = search_point_with_rotss->size();
    shared_ptr<std::vector<RifDockResult>> rif_dock_results = make_shared<std::vector<RifDockResult>>( length );

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        (*rif_dock_results)[i] = (*search_point_with_rotss)[i];
    }

    return rif_dock_results;
}


shared_ptr<std::vector<SearchPoint>>
search_points_from_rif_dock_results(shared_ptr<std::vector<RifDockResult>> rif_dock_results) {

    size_t length = rif_dock_results->size();
    shared_ptr<std::vector<SearchPoint>> search_points = make_shared<std::vector<SearchPoint>>( length );

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        (*search_points)[i] = (*rif_dock_results)[i];
    }

    return search_points;
}
shared_ptr<std::vector<SearchPointWithRots>>
search_point_with_rotss_from_rif_dock_results(shared_ptr<std::vector<RifDockResult>> rif_dock_results) {

    size_t length = rif_dock_results->size();
    shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss = make_shared<std::vector<SearchPointWithRots>>( length );

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        (*search_point_with_rotss)[i] = (*rif_dock_results)[i];
    }

    return search_point_with_rotss;
}



}}


