// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rif_requirements_util_hh
#define INCLUDED_riflib_rif_requirements_util_hh

#include <string>
#include <vector>
#include <utility/string_util.hh>

// TODO:: change this tuning file manager to a factory ..... Singleton ?

namespace devel {
    namespace scheme {
        namespace rif {
            
            struct DonorDefinition {
                int res_num = -1;
                std::vector< std::string > allowed_donor_res;
            };
            
            struct AcceptorDefinition {
                int res_num = -1;
                std::vector< std::string > allowed_acceptor_res;
            };
            
            struct HBondDefinition {
                std::string atom_name = "";
                int			    res_num = -1;
                std::vector< std::string > allowed_rot_names;
            };
            
            struct BidentateDefinition {
                std::string atom1_name = "";
                int res1_num = -1;
                std::string atom2_name = "";
                int res2_num = -1;
            };
            
            struct HotspotRequirement {
                int req_num = -1;
                int hotspot_num = -1;
            };
            
            struct HbondRequirement {
                int req_num = -1;
                std::string atom_name = "";
                int res_num = -1;
            };
            
            struct BidentateRequirement {
                int req_num = -1;
                std::string atom1_name = "";
                int res1_num = -1;
                std::string atom2_name = "";
                int res2_num = -1;
            };
            
            
            struct ApoReqTerm {
                std::string atom_name = "";
                int res_num = -1;
                float distance = 0;
            };
            struct ApoRequirement {
                int req_num = -1;
                std::vector<std::string> allowed_rot_names;
                std::vector<ApoReqTerm> terms;
                HbondRequirement apo_hbond;
            };
            
            struct PiPiStackingTerm {
                std::string atom_name = "";
                int res_num = -1;
            };
            struct PiPiStackingRequirement {
                int req_num = -1;
                std::vector<std::string> allowed_rot_names;
                std::vector<PiPiStackingTerm> terms;
                HbondRequirement apo_hbond;
            };
            
            struct CationPiRequirement {
                int req_num = -1;
                int res_num = -1;
                std::vector<std::string> allowed_rot_names;
                HbondRequirement apo_hbond;
            };
            
            std::vector < DonorDefinition > get_donor_definitions( std::string tuning_file );
            std::vector < AcceptorDefinition > get_acceptor_definitions( std::string tuning_file );
            std::vector< HBondDefinition > get_hbond_definitions( std::string tuning_file );
            std::vector< BidentateDefinition > get_bidentate_definitions( std::string tuning_file );
            std::vector< HbondRequirement > get_hbond_requirement_definitions( std::string tuning_file );
            std::vector< BidentateRequirement > get_bidentate_requirement_definitions( std::string tuning_file );
            std::vector< HotspotRequirement > get_hotspot_requirement_definitions( std::string tuning_file );
            std::vector< ApoRequirement > get_Apo_requirement_definitions( std::string tuning_file );
            std::vector< PiPiStackingRequirement > get_pipi_stacking_requirement_definitions( std::string tuning_file );
            std::vector< CationPiRequirement > get_cationpi_requirement_definitions( std::string tuning_file );
            bool check_requirement_definition_exists( std::string tuning_file );

            template< class Requirement >
            void
            clear_apo_hbond( Requirement & req ) {
                req.apo_hbond.atom_name = "";
                req.apo_hbond.res_num = -1;
            }

            template< class Requirement >
            utility::vector1<std::string>
            parse_apo_hbond( utility::vector1<std::string> const & split, size_t start, Requirement & req ) {
                size_t hbond_start = 0;
                for ( size_t i = start; i <= split.size(); i++ ) {
                    if ( split[i] == "HBOND" ) {
                        hbond_start = i;
                        break;
                    }
                }
                if ( hbond_start == 0 ) return split;

                if ( hbond_start + 2 > split.size() ) {
                    utility_exit_with_message("Bad apo-hbond definition");
                }
                req.apo_hbond.atom_name = split[hbond_start+1];
                req.apo_hbond.res_num = utility::string2int(split[hbond_start+2]);

                utility::vector1<std::string> new_split;
                for ( size_t i = 1; i < hbond_start; i++ ) new_split.push_back(split[i]);
                return new_split;
            }

        }
    }
}

#endif