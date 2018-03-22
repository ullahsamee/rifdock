// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_UnsatManager_hh
#define INCLUDED_riflib_UnsatManager_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>

#include <core/pose/Pose.hh>

#include <string>
#include <vector>

#include <boost/any.hpp>



namespace devel {
namespace scheme {

namespace hbond {

enum AType {
    Hydroxyl,
    TyrHydroxyl,
    HydroxylH,
    TyrHydroxylH,
    BBAmide,
    BBCarbonyl,
    Amide,
    AmideCarbonyl,
    ArgNH,
    ArgNE,
    LysN,
    Carboxylate,
    TrypN,
    HisN,
    HisH,
    Unknown
};

std::vector<int> NumOrbitals {
    2,  // "Hydroxyl",
    2,  // "TyrHydroxyl",
    1,  // "HydroxylH",
    1,  // "TyrHydroxylH",
    1,  // "BBAmide",
    2,  // "BBCarbonyl",
    2,  // "Amide",
    2,  // "AmideCarbonyl",
    2,  // "ArgNH",
    1,  // "ArgNE",
    3,  // "LysN",
    2,  // "Carboxylate",
    1,  // "TrypN",
    1,  // "HisN",
    1,  // "HisH",
    0  // "Unknown"
};

std::vector<std::string> ATypeNames {
    "Hydroxyl",
    "TyrHydroxyl",
    "HydroxylH",
    "TyrHydroxylH",
    "BBAmide",
    "BBCarbonyl",
    "Amide",
    "AmideCarbonyl",
    "ArgNH",
    "ArgNE",
    "LysN",
    "Carboxylate",
    "TrypN",
    "HisN",
    "HisH",
    "Unknown"
};

// This is really arbitrary - brian
// Loosely based on this http://www.biochem.ucl.ac.uk/bsm/atlas/index.html

// penalty__first_atom = DefaultUnsatPenalties[X][0]
// penalty_second_atom = DefaultUnsatPenalties[X][0] * DefaultUnsatPenalties[X][1]
// penalty__third_atom = DefaultUnsatPenalties[X][0] * ( 1  - 2*( 1 - DefaultUnsatPenalties[X][1] ) )

// first number is penalty for last unsat
// second number is ratio of that for second one
//    the third one is 1 - ( 2*(1 - ratio ) )
std::vector<std::vector<float>> DefaultUnsatPenalties {
    std::vector<float> {0, 0},    // "Hydroxyl",
    std::vector<float> {0, 0},    // "TyrHydroxyl",
    std::vector<float> {4, 0},    // "HydroxylH",
    std::vector<float> {4, 0},    // "TyrHydroxylH",
    std::vector<float> {4, 0},    // "BBAmide",
    std::vector<float> {4, 0},    // "BBCarbonyl",
    std::vector<float> {4, 0.5},    // "Amide",
    std::vector<float> {4, 0},    // "AmideCarbonyl",
    std::vector<float> {3, 1},    // "ArgNH",
    std::vector<float> {3, 0},    // "ArgNE",
    std::vector<float> {4, 1},    // "LysN",
    std::vector<float> {3, 1},    // "Carboxylate",
    std::vector<float> {4, 0},    // "TrypN",
    std::vector<float> {0, 0},    // "HisN",
    std::vector<float> {4, 0},    // "HisH"
    std::vector<float> {0, 0}    // "Unknown"
};


// std::vector<std::vector<float>> DefaultUnsatPenalties {
//     std::vector<float> {6, 2, 0, 0},
//     std::vector<float> {6, 2, 0, 0},
//     std::vector<float> {6, 2, 0, 0},
//     std::vector<float> {6, 2, 0, 0},
//     std::vector<float> {4, 0},
//     std::vector<float> {6, 2, 0},
//     std::vector<float> {6, 2, 0},
//     std::vector<float> {6, 2, 0},
//     std::vector<float> {8, 4, 0},
//     std::vector<float> {4, 0},
//     std::vector<float> {12, 8, 4, 0},
//     std::vector<float> {8, 4, 0},
//     std::vector<float> {4, 0},
//     std::vector<float> {4, 0}
// };

struct HeavyAtom {
    int AType = 0;
    int resid = 0;
    std::vector<int> sat_groups;
    std::string name = "";
    Eigen::Vector3f xyz;
};


AType
identify_donor( std::string const & dname, std::string & heavy_atom_name );

AType
identify_acceptor( std::string const & aname, char one_letter, std::string & heavy_atom_name );



}

struct UnsatManager {

    UnsatManager() {} // used by clone()

    // UnsatManager(
    // ) 
    // {

    // }


    void
    set_target_donors_acceptors( 
        core::pose::Pose const & target, 
        std::vector<HBondRay> const & target_donors, 
        std::vector<HBondRay> const & target_acceptors,
        std::vector<std::pair<int,std::string>> const & donor_anames,
        std::vector<std::pair<int,std::string>> const & acceptor_anames );

    void
    find_target_presatisfied(
        core::pose::Pose const & target
    );

    void
    dump_presatisfied();

    std::vector<Eigen::Vector3f>
    get_heavy_atom_xyzs();


private:
    int
    find_heavy_atom(
        int resid,
        std::string const & heavy_atom_name
    );

    void
    add_target_sat( 
        int resid,
        std::string const & heavy_atom_name,
        hbond::AType heavy_atom_type,
        Eigen::Vector3f const & xyz,
        int sat
    );


    bool
    validate_heavy_atoms();

private:

    int num_donors_;
    std::vector<HBondRay> target_donors_acceptors_;
    std::vector<hbond::HeavyAtom> target_heavy_atoms_;
    std::vector<bool> target_presatisfied_;


};




}}

#endif
