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

#include <scheme/search/HackPack.hh>
#include <riflib/BurialManager.hh>

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

const std::vector<int> NumOrbitals {
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

const std::vector<std::string> ATypeNames {
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
// See UnsatManager::prepare_packer for an explanation of why the scores are formatted like this

// penalty__first_atom = DefaultUnsatPenalties[X][0]
// penalty_second_atom = DefaultUnsatPenalties[X][0] * DefaultUnsatPenalties[X][1]
// penalty__third_atom = DefaultUnsatPenalties[X][0] * ( 1  - 2*( 1 - DefaultUnsatPenalties[X][1] ) )

// first number is penalty for last unsat
// second number is ratio of that for second one
//    the third one is 1 - ( 2*(1 - ratio ) )
const std::vector<std::vector<float>> DefaultUnsatPenalties {
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
const std::vector<std::vector<float>> ScottUnsatPenalties {
    std::vector<float> {0, 0},    // "Hydroxyl",
    std::vector<float> {0, 0},    // "TyrHydroxyl",
    std::vector<float> {4, 0},    // "HydroxylH",
    std::vector<float> {4, 0},    // "TyrHydroxylH",
    std::vector<float> {4, 0},    // "BBAmide",
    std::vector<float> {4, 0},    // "BBCarbonyl",
    std::vector<float> {4, 0},    // "Amide",
    std::vector<float> {4, 0},    // "AmideCarbonyl",
    std::vector<float> {4, 0},    // "ArgNH",
    std::vector<float> {4, 0},    // "ArgNE",
    std::vector<float> {4, 0.5},    // "LysN",  // making this 0 is glitchy
    std::vector<float> {4, 0},    // "Carboxylate",
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


enum Patch {
    NONE,
    NOT_BURIED,
    ACTUALLY_HBONDING
};


}

struct ToPackRot {
    int ires;
    int irot;
    float score;
    int sat1;
    int sat2;
    ToPackRot( int _ires, int _irot, float _score, int _sat1, int _sat2 ) :
        ires(_ires), irot(_irot), score(_score), sat1(_sat1), sat2(_sat2) {}
};

struct UnsatManager {

    UnsatManager() {} // used by clone()

    UnsatManager( 
        std::vector<std::vector<float>> const & unsat_penalties, 
        shared_ptr< RotamerIndex > rot_index_p,
        float unsat_score_scalar,
        int require_burial,
        float score_offset,
        bool debug,
        bool store_common_unsats
    );

    shared_ptr<UnsatManager>
    clone() const;

    void
    reset();

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

    std::vector<float>
    get_initial_unsats( std::vector<float> const & burial_weights ) const;

    std::vector<float>
    get_buried_unsats( 
        std::vector<float> const & burial_weights,
        std::vector< std::pair<intRot,intRot> > const & rotamers,
        std::vector< EigenXform > const & bb_positions,
        RifScoreRotamerVsTarget const & rot_tgt_scorer
    ) const;

    void
    print_buried_unsats( std::vector<float> const & unsat_penalties, int scaff_size=-1) const;

    void
    print_unsats_help( std::vector<float> const & unsat_penalties) const;

    std::vector<Eigen::Vector3f>
    get_heavy_atom_xyzs();

    std::vector<bool> const &
    get_presatisfied() { return target_presatisfied_; }

    float
    calculate_nonpack_score( 
        std::vector<float> const & burial_weights,
        std::vector<bool> const & is_satisfied
    );

    float
    calculate_unsat_score( float P0, float P1, int wants, int satisfied, float weight ) const;


    void
    add_to_pack_rot( 
        int ires,
        int irot,
        float score,
        int sat1,
        int sat2
    );

    float
    prepare_packer( 
        ::scheme::search::HackPack & packer, 
        std::vector<float> const & burial_weights,
        std::vector<bool> const & pre_and_bb_satisfied
    );

    void
    insert_to_pack_rots_into_packer( ::scheme::search::HackPack & packer );

    void
    fix_packer( 
        ::scheme::search::HackPack & packer, 
        shared_ptr<::scheme::objective::storage::TwoBodyTable<float> const> reference_twobody
    );

    bool
    patch_heavy_atoms( 
        int resid,
        std::string const & heavy_atom,
        hbond::Patch patch,
        shared_ptr<BurialManager> const & burial_manager
    );

    void
    apply_unsat_helper( 
        std::string const & unsat_helper_fname,
        shared_ptr<BurialManager> const & burial_manager
    );

    void
    sum_unsat_counts( UnsatManager const & other );

    void
    print_unsat_counts() const;

// private:
    int
    find_heavy_atom(
        int resid,
        std::string const & heavy_atom_name,
        bool strict = false
    );

    void
    add_target_sat( 
        int resid,
        std::string const & heavy_atom_name,
        hbond::AType heavy_atom_type,
        Eigen::Vector3f const & xyz,
        int sat
    );

    float
    handle_twobody( 
        int satisfier1,
        int satisfier2,
        float penalty,
        ::scheme::search::HackPack & packer
    );

    bool
    validate_heavy_atoms();

// private:

    int num_donors_;
    std::vector<HBondRay> target_donors_acceptors_;
    std::vector<hbond::HeavyAtom> target_heavy_atoms_;
    std::vector<bool> target_presatisfied_;
    std::vector<std::vector<float>> unsat_penalties_;
    std::vector<std::vector<float>> total_first_twob_;  // total penalty, P0, P0 * (1 - P1)
    shared_ptr< RotamerIndex > rot_index_p;
    float unsat_score_scalar_;
    int require_burial_;
    float score_offset_;
    bool debug_;

    bool store_common_unsats_;
    std::vector<uint64_t> unsat_counts_;

// things that are resetable
    std::vector<ToPackRot> to_pack_rots_;
    std::vector<std::array<int, 4>> modified_edges_;

};




}}

#endif
