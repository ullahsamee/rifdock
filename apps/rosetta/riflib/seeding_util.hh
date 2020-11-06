// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_seeding_util_hh
#define INCLUDED_riflib_seeding_util_hh

#include <riflib/types.hh>
#include <core/types.hh>

#include <rif_dock_test.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/task/types.hh>



namespace devel {
namespace scheme {


shared_ptr<std::vector<EigenXform>>
setup_seeding_positions( RifDockOpt & opt, ProtocolData & pd, ScaffoldProviderOP & scaffold_provider, int iscaff );

bool 
parse_exhausitive_searching_file(
    std::string fname, 
    std::vector<std::pair< int64_t, devel::scheme::EigenXform > > & searching_positions, 
    double maximum_ang = 999
);

bool 
parse_seeding_file(
    std::string fname, 
    std::vector<devel::scheme::EigenXform> & seeding_positions, 
    bool seeding_by_patchdock,
		float patchdock_min_sasa = -1000.0,
		int patchdock_top_ranks = 99999
);

template<class AnyPoint>
_AnyPointVectorsMap<AnyPoint>
sort_into_blocks( 
    shared_ptr<std::vector<AnyPoint>> const & any_points,
    bool treat_nests_differently,
    bool treat_seeds_differently,
    bool treat_scaffolds_differently ) {

    typedef _AnyPointVectorsMap<AnyPoint> AnyPointVectorsMap;

    SelectiveRifDockIndexHasher   hasher( treat_nests_differently, treat_seeds_differently, treat_scaffolds_differently );
    SelectiveRifDockIndexEquater equater( treat_nests_differently, treat_seeds_differently, treat_scaffolds_differently );

    AnyPointVectorsMap map (10000, hasher, equater );

    for ( AnyPoint pt : *any_points) {
        map[pt.index].push_back( pt );
    }

    return map;

}


void dump_xform_file(
    DirectorBase const & director,
    ScenePtr const & scene_minimal,
    double cart_radius,
    double cart_resl,
    double angle_radius,
    double angle_resl
);

}
}

#endif
