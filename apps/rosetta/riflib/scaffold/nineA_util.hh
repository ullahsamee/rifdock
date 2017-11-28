// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_nineA_util_hh
#define INCLUDED_riflib_scaffold_nineA_util_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <rif_dock_test.hh>

#include <core/pose/Pose.hh>

#include <utility/io/izstream.hh>

#include <string>
#include <vector>



namespace devel {
namespace scheme {


struct NineAMember {
    ParametricSceneConformationCOP conformation;
    uint64_t cd_index;
    uint64_t clust;
    uint64_t gindex;
    std::vector<MorphAction> morph_history;
    ::scheme::scaffold::TreeRelation tree_relation;
};


namespace DanielAtomIndex {
    const uint64_t N = 3; 
    const uint64_t CA = 0; 
    const uint64_t C = 1; 
    const uint64_t O = 2; 
}

struct StatRow {
    std::vector<uint32_t> int_fields;
    std::vector<float> float_fields;
    std::string filename;
};

typedef std::vector<StatRow> StatTable;

enum IntFields {
    Clust = 0,
    NumCon,
    PrevAssig,
    last_int_field
};

enum FloatFields {
    AvgRad = 0,
    MaxRad,
    X_1,
    X_2,
    X_3,
    X_4,
    X_5,
    X_6,
    X_7,
    X_8,
    X_9,
    X_10,
    X_11,
    X_12,
    X_13,
    X_14,
    X_15,
    X_16,
    X_17,
    X_18,
    X_19,
    X_20,
    X_21,
    X_22,
    X_23,
    X_24,
    X_25,
    X_26,
    X_27,
    X_28,
    X_29,
    X_30,
    X_31,
    X_32,
    X_33,
    X_34,
    X_35,
    X_36,
    Y_1,
    Y_2,
    Y_3,
    Y_4,
    Y_5,
    Y_6,
    Y_7,
    Y_8,
    Y_9,
    Y_10,
    Y_11,
    Y_12,
    Y_13,
    Y_14,
    Y_15,
    Y_16,
    Y_17,
    Y_18,
    Y_19,
    Y_20,
    Y_21,
    Y_22,
    Y_23,
    Y_24,
    Y_25,
    Y_26,
    Y_27,
    Y_28,
    Y_29,
    Y_30,
    Y_31,
    Y_32,
    Y_33,
    Y_34,
    Y_35,
    Y_36,
    Z_1,
    Z_2,
    Z_3,
    Z_4,
    Z_5,
    Z_6,
    Z_7,
    Z_8,
    Z_9,
    Z_10,
    Z_11,
    Z_12,
    Z_13,
    Z_14,
    Z_15,
    Z_16,
    Z_17,
    Z_18,
    Z_19,
    Z_20,
    Z_21,
    Z_22,
    Z_23,
    Z_24,
    Z_25,
    Z_26,
    Z_27,
    Z_28,
    Z_29,
    Z_30,
    Z_31,
    Z_32,
    Z_33,
    Z_34,
    Z_35,
    Z_36,
    Phi_1,
    Psi_1,
    Ome_1,
    Phi_2,
    Psi_2,
    Ome_2,
    Phi_3,
    Psi_3,
    Ome_3,
    Phi_4,
    Psi_4,
    Ome_4,
    Phi_5,
    Psi_5,
    Ome_5,
    Phi_6,
    Psi_6,
    Ome_6,
    Phi_7,
    Psi_7,
    Ome_7,
    Phi_8,
    Psi_8,
    Ome_8,
    Phi_9,
    Psi_9,
    Ome_9,
    last_float_field
};


const std::vector<std::string> cluster_data_names {
    "_res5.22",
    "_res5.26",
    "_res5.31",
    "_res5.24",
    "_res5.32",
    "_res2.23",
    "_res1.90",
    "_res1.86",
    "_res1.65",
    "_res1.58",
    "final_res1.41"
};

const uint64_t num_clusters = cluster_data_names.size();

inline
StatTable
load_stat_table( uint64_t cdindex, RifDockOpt opt ) {

    std::string filename = opt.nineA_cluster_path 
                            + "/kcenters_stats_al1.dat" 
                            + cluster_data_names.at( cdindex );

    utility::izstream f( filename );

}






}}

#endif
