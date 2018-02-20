// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_task_util_hh
#define INCLUDED_riflib_rifdock_task_util_hh


#include <riflib/types.hh>

#include <riflib/rifdock_subroutines/util.hh>



using ::scheme::make_shared;
using ::scheme::shared_ptr;


namespace devel {
namespace scheme {


shared_ptr<std::vector<SearchPointWithRots>>
search_point_with_rotss_from_search_points(shared_ptr<std::vector<SearchPoint>> search_points);
shared_ptr<std::vector<RifDockResult>>
rif_dock_results_from_search_points(shared_ptr<std::vector<SearchPoint>> search_points);


shared_ptr<std::vector<SearchPoint>>
search_points_from_search_point_with_rotss(shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss);
shared_ptr<std::vector<RifDockResult>>
rif_dock_results_from_search_point_with_rotss(shared_ptr<std::vector<SearchPointWithRots>> search_point_with_rotss);


shared_ptr<std::vector<SearchPoint>>
search_points_from_rif_dock_results(shared_ptr<std::vector<RifDockResult>> rif_dock_results);
shared_ptr<std::vector<SearchPointWithRots>>
search_point_with_rotss_from_rif_dock_results(shared_ptr<std::vector<RifDockResult>> rif_dock_results);



}}



#endif