// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:



#include <riflib/types.hh>

#include <riflib/task/util.hh>



using ::scheme::make_shared;
using ::scheme::shared_ptr;


namespace devel {
namespace scheme {


std::vector<SearchPointWithRots>
search_point_with_rotss_from_search_points(std::vector<SearchPoint> const & search_points) {

    size_t length = search_points.size();
    std::vector<SearchPointWithRots> search_point_with_rotss( length );

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        search_point_with_rotss[i] = search_points[i];
    }

    return search_point_with_rotss;
}
std::vector<RifDockResult>
rif_dock_results_from_search_points(std::vector<SearchPoint> const & search_points) {

    size_t length = search_points.size();
    std::vector<RifDockResult> rif_dock_results( length);

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        rif_dock_results[i] = search_points[i];
    }

    return rif_dock_results;
}


std::vector<SearchPoint>
search_points_from_search_point_with_rotss(std::vector<SearchPointWithRots> const & search_point_with_rotss) {

    size_t length = search_point_with_rotss.size();
    std::vector<SearchPoint> search_points( length);

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        search_points[i] = search_point_with_rotss[i];
    }

    return search_points;
}
std::vector<RifDockResult>
rif_dock_results_from_search_point_with_rotss(std::vector<SearchPointWithRots> const & search_point_with_rotss) {

    size_t length = search_point_with_rotss.size();
    std::vector<RifDockResult> rif_dock_results( length );

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        rif_dock_results[i] = search_point_with_rotss[i];
    }

    return rif_dock_results;
}


std::vector<SearchPoint>
search_points_from_rif_dock_results(std::vector<RifDockResult> const & rif_dock_results) {

    size_t length = rif_dock_results.size();
    std::vector<SearchPoint> search_points( length );

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        search_points[i] = rif_dock_results[i];
    }

    return search_points;
}
std::vector<SearchPointWithRots>
search_point_with_rotss_from_rif_dock_results(std::vector<RifDockResult> const & rif_dock_results) {

    size_t length = rif_dock_results.size();
    std::vector<SearchPointWithRots> search_point_with_rotss( length );

    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for ( size_t i = 0; i < length; i++) {
        search_point_with_rotss[i] = rif_dock_results[i];
    }

    return search_point_with_rotss;
}



}}


