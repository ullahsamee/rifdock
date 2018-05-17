

#include <riflib/UnsatManager.hh>

#include <riflib/RotamerGenerator.hh>
#include <riflib/util.hh>
#include <riflib/ScoreRotamerVsTarget.hh>

#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/format.hh>
#include <boost/format.hpp>


using Eigen::Vector3f;

namespace devel {
namespace scheme {



namespace hbond {



AType
identify_donor( std::string const & dname, std::string & heavy_atom_name ) {

    AType heavy_atom_type = Unknown;
    if ( " H  " == dname ) {
        heavy_atom_name = " N  ";
        heavy_atom_type = BBAmide;

    } else if ( " HE2" == dname ) {
        heavy_atom_name = " NE2";
        heavy_atom_type = HisN;

    } else if ( " HD1" == dname ) {
        heavy_atom_name = " ND1";
        heavy_atom_type = HisN;

    } else if ( "1HZ " == dname || "2HZ " == dname || "3HZ " == dname ) {
        heavy_atom_name = " NZ ";
        heavy_atom_type = LysN;

    } else if ( "1HD2" == dname || "2HD2" == dname ) {
        heavy_atom_name = " ND2";
        heavy_atom_type = Amide;

    } else if ( "1HE2" == dname || "2HE2" == dname ) {
        heavy_atom_name = " NE2";
        heavy_atom_type = Amide;

    } else if ( " HE " == dname ) {
        heavy_atom_name = " NE ";
        heavy_atom_type = ArgNE;

    } else if ( "1HH2" == dname || "2HH2" == dname ) {
        heavy_atom_name = " NH2";
        heavy_atom_type = ArgNH;

    } else if ( "1HH1" == dname || "2HH1" == dname ) {
        heavy_atom_name = " NH1";
        heavy_atom_type = ArgNH;

    } else if ( " HG " == dname ) {
        heavy_atom_name = " HG ";
        heavy_atom_type = HydroxylH;

    } else if ( " HG1" == dname ) {
        heavy_atom_name = " HG1";
        heavy_atom_type = HydroxylH;

    } else if ( " HE1"  == dname) {
        heavy_atom_name = " NE1";
        heavy_atom_type = TrypN;

    } else if ( " HH " == dname ) {
        heavy_atom_name = " HH ";
        heavy_atom_type = TyrHydroxylH;

    }
    
    return heavy_atom_type;
}

AType
identify_acceptor( std::string const & aname, char one_letter, std::string & heavy_atom_name ) {

    AType heavy_atom_type = Unknown;
    if ( " O  " == aname ) {
        heavy_atom_name = " O  ";
        heavy_atom_type = BBCarbonyl;

    } else if ( " OD1" == aname ) {
        if ( 'D' == one_letter ) {
            heavy_atom_name = " OD1";
            heavy_atom_type = Carboxylate;

        } 
        if ( 'N' == one_letter ) {
            heavy_atom_name = " OD1";
            heavy_atom_type = AmideCarbonyl;

        }
    } else if ( " OD2" == aname ) {
        if ( 'D' == one_letter ) {
            heavy_atom_name = " OD2";
            heavy_atom_type = Carboxylate;

        } 
        if ( 'N' == one_letter ) {
            heavy_atom_name = " OD2";
            heavy_atom_type = AmideCarbonyl;

        }
    } else if ( " OE1" == aname ) {
        if ( 'E' == one_letter ) {
            heavy_atom_name = " OE1";
            heavy_atom_type = Carboxylate;

        }
        if ( 'Q' == one_letter ) {
            heavy_atom_name = " OE1";
            heavy_atom_type = AmideCarbonyl;

        }
    } else if ( " OE2" == aname ) {
        if ( 'E' == one_letter ) {
            heavy_atom_name = " OE2";
            heavy_atom_type = Carboxylate;

        }
        if ( 'Q' == one_letter ) {
            heavy_atom_name = " OE2";
            heavy_atom_type = AmideCarbonyl;

        }
    } else if ( " ND1" == aname ) {
        heavy_atom_name = " ND1";
        heavy_atom_type = HisN;

    } else if ( " NE2" == aname ) {
        heavy_atom_name = " NE2";
        heavy_atom_type = HisN;

    } else if ( " OG " == aname ) {
        heavy_atom_name = " OG ";
        heavy_atom_type = Hydroxyl;

    } else if ( " OG1" == aname ) {
        heavy_atom_name = " OG1";
        heavy_atom_type = Hydroxyl;

    } else if ( " OH " == aname ) {
        heavy_atom_name = " OH ";
        heavy_atom_type = TyrHydroxyl;

    }
    return heavy_atom_type;
}




}

UnsatManager::UnsatManager( 
    std::vector<std::vector<float>> const & unsat_penalties,
    shared_ptr< RotamerIndex > rot_index_p_in,
    float unsat_score_scalar,
    int require_burial,
    float score_offset,
    bool debug,
    bool store_common_unsats
) :
    unsat_penalties_( unsat_penalties ),
    rot_index_p( rot_index_p_in ),
    require_burial_( require_burial ),
    unsat_score_scalar_( unsat_score_scalar ),
    score_offset_( score_offset ),
    debug_( debug ),
    store_common_unsats_( store_common_unsats )
{

    using ObjexxFCL::format::F;
    std::cout << "Buried orbital/hydrogen penalties: " << std::endl;
    std::cout << "  # of unsats:   0    1    2    3  " << std::endl;
    runtime_assert( unsat_penalties_.size() == hbond::DefaultUnsatPenalties.size() );
    total_first_twob_.resize(unsat_penalties_.size(), std::vector<float>(3));

    for ( int i = 0; i < total_first_twob_.size(); i++ ) {
        std::vector<float> const & penalties = unsat_penalties_[i];
        
        const float P0 = penalties[0];
        const float P1 = penalties[1];
        const int wants = hbond::NumOrbitals[i];

        total_first_twob_[i][0] = calculate_unsat_score( P0, P1, wants, 0, 1 ) / unsat_score_scalar;
        total_first_twob_[i][1] = P0;
        total_first_twob_[i][2] = P0 * ( 1 - P1 );

        for (int j = 0; j < 13 - hbond::ATypeNames[i].length(); j++ ) std::cout << " ";

        std::cout << hbond::ATypeNames[i] << ": ";
        for ( int j = 0; j <= wants; j++ ) {
            std::cout << F(5,1, calculate_unsat_score( P0, P1, wants, wants-j, 1 ) );
        }
        std::cout << std::endl;
    }


}

shared_ptr<UnsatManager>
UnsatManager::clone() const {
    shared_ptr<UnsatManager> ot = make_shared<UnsatManager>();

    *ot = *this;

    return ot;
}

void
UnsatManager::reset() {
    to_pack_rots_.clear();  // this supposedly doesn't mess with the memory
    to_pack_rots_.reserve(512);

    modified_edges_.clear();
    modified_edges_.reserve(1024);
}


void
UnsatManager::set_target_donors_acceptors( 
    core::pose::Pose const & target, 
    std::vector<HBondRay> const & target_donors, 
    std::vector<HBondRay> const & target_acceptors,
    std::vector<std::pair<int,std::string>> const & donor_anames,
    std::vector<std::pair<int,std::string>> const & acceptor_anames ) {

    using hbond::AType;

    std::string sequence = target.sequence();

    target_heavy_atoms_.resize(0);

    int so_far = 0;
    for ( int isat = 0; isat < target_donors.size(); isat++ ) {

        HBondRay const & ray = target_donors[isat];

        Vector3f xyz_; 
        while ( so_far < donor_anames.size() ) {
            std::pair<int,std::string> const & atom_pair = donor_anames[so_far];
            numeric::xyzVector<core::Real> xyz = target.residue(atom_pair.first).xyz(atom_pair.second);
            xyz_[0] = xyz[0]; xyz_[1] = xyz[1]; xyz_[2] = xyz[2];

            float dist = ( xyz_ - ray.horb_cen ).norm();

            if (dist < 0.1) break;

            so_far++;
        }

        if ( so_far == donor_anames.size() ) {
            utility_exit_with_message("Error with hydrogen bond donors!!!");
        }
        std::pair<int,std::string> const & atom_pair = donor_anames[so_far];

        std::string heavy_atom_name = "";
        AType heavy_atom_type = hbond::identify_donor( atom_pair.second, heavy_atom_name );
        
        if ( heavy_atom_name == "" ) {
            utility_exit_with_message("Unknown hydrogen bond donor type: " + atom_pair.second);
        }

        add_target_sat( atom_pair.first, heavy_atom_name, heavy_atom_type, xyz_, isat );
        so_far++;
    }



    so_far = 0;
    for ( int isat = 0; isat < target_acceptors.size(); isat++ ) {

        HBondRay const & ray = target_acceptors[isat];

        Vector3f xyz_; 
        while ( so_far < acceptor_anames.size() ) {
            std::pair<int,std::string> const & atom_pair = acceptor_anames[so_far];
            numeric::xyzVector<core::Real> xyz = target.residue(atom_pair.first).xyz(atom_pair.second);
            xyz_[0] = xyz[0]; xyz_[1] = xyz[1]; xyz_[2] = xyz[2];

            float dist = ( xyz_ - ( ray.horb_cen - ray.direction*ORBLEN ) ).norm();

            if (dist < 0.1) break;
            so_far++;
        }

        if (so_far == acceptor_anames.size() ) {
            utility_exit_with_message("Error with hydrogen bond acceptors!!!");
        }
        std::pair<int,std::string> const & atom_pair = acceptor_anames[so_far];


        std::string heavy_atom_name = "";
        AType heavy_atom_type = hbond::identify_acceptor( atom_pair.second, sequence[atom_pair.first-1], heavy_atom_name );
        
        if ( heavy_atom_name == "" ) {
            std::cout << atom_pair.first << std::endl;
            utility_exit_with_message("Unknown hydrogen bond acceptor type: " + atom_pair.second + " " + sequence[atom_pair.first-1]);
        }

        add_target_sat( atom_pair.first, heavy_atom_name, heavy_atom_type, xyz_, isat + target_donors.size() );

        so_far ++;
    }

    target_donors_acceptors_ = target_donors;
    target_donors_acceptors_.insert( target_donors_acceptors_.end(), target_acceptors.begin(), target_acceptors.end());
    num_donors_ = target_donors.size();


    if (debug_) {
        for ( int i = 0; i < target_heavy_atoms_.size(); i++ ) {
            hbond::HeavyAtom const & ha = target_heavy_atoms_[i];
            std::cout << "Heavy Atom: " << i << " res: " << ha.resid << " name: " << ha.name 
                      << " type: " << hbond::ATypeNames[ha.AType] << std::endl;
        }
    }


    if ( store_common_unsats_ ) {
        unsat_counts_.resize( target_heavy_atoms_.size(), 0 );
    }

    runtime_assert( validate_heavy_atoms() );


}

int
UnsatManager::find_heavy_atom( int resid, std::string const & heavy_atom_name, bool strict /* =false */ ) {
    int heavy_atom_no = 0;
    for ( heavy_atom_no = 0; heavy_atom_no < target_heavy_atoms_.size(); heavy_atom_no++ ) {
        hbond::HeavyAtom const & ha = target_heavy_atoms_[heavy_atom_no];
        if ( ha.resid != resid ) continue;
        if ( strict ) {
            if ( ha.name != heavy_atom_name ) continue;
        } else {
            if ( utility::strip(ha.name) != utility::strip(heavy_atom_name) ) continue;
        }
        break;
    }

    if (heavy_atom_no == target_heavy_atoms_.size()) {
        return -1;
    }

    return heavy_atom_no;
}

void
UnsatManager::add_target_sat( 
    int resid,
    std::string const & heavy_atom_name,
    hbond::AType heavy_atom_type,
    Eigen::Vector3f const & xyz,
    int sat
) {
    int heavy_atom_no = find_heavy_atom( resid, heavy_atom_name );

    if (heavy_atom_no == -1) {
        hbond::HeavyAtom ha;
        ha.resid = resid;
        ha.name = heavy_atom_name;
        ha.AType = heavy_atom_type;
        ha.xyz = xyz;
        target_heavy_atoms_.push_back(ha);
        heavy_atom_no = target_heavy_atoms_.size()-1;
    }

    hbond::HeavyAtom & ha = target_heavy_atoms_[heavy_atom_no];
    ha.sat_groups.push_back(sat);   

}


bool
UnsatManager::validate_heavy_atoms() {
    bool good = true;
    for ( hbond::HeavyAtom const & ha : target_heavy_atoms_ ) {
        if ( ha.sat_groups.size() != hbond::NumOrbitals[ha.AType] ) {
            good = false;
            std::cout << "Error: Residue " << ha.resid << " atom " << ha.name << " wrong number of acceptor/donors. Has " << ha.sat_groups.size()
                    << " wants " << hbond::NumOrbitals[ha.AType] << std::endl;
        }
    }
    if ( store_common_unsats_ ) {
        runtime_assert( unsat_counts_.size() == target_heavy_atoms_.size() );
    }
    return good;
}

void
UnsatManager::find_target_presatisfied(
    core::pose::Pose const & target
) {
    std::cout << "UnsatManager: finding target pre-satisfied donors/acceptors" << std::endl;

    runtime_assert( target_heavy_atoms_.size() > 0 );

    target_presatisfied_.resize(0);
    target_presatisfied_.resize(target_donors_acceptors_.size());

    std::string sequence = target.sequence();
    std::vector<std::pair<int,int>> presat_res_atoms = ::devel::scheme::get_satisfied_atoms_vector(target, -0.00000001f);

    runtime_assert( presat_res_atoms.size() % 2 == 0 );

    Vector3f xyz;
    Vector3f last_xyz;
    for ( int ipresat = 0; ipresat < presat_res_atoms.size(); ipresat++ ) {
        last_xyz = xyz;

        bool is_donor = ipresat % 2 == 0;

        std::pair<int, int> const & res_atomno = presat_res_atoms[ipresat];
        int atomno = res_atomno.first;
        int resid = res_atomno.second;
        std::string atom_name = target.residue(resid).atom_name(atomno);
        numeric::xyzVector<core::Real> _xyz = target.residue(resid).xyz(atomno);
        xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

        std::string heavy_atom_name = "";
        if ( is_donor ) {
            hbond::identify_donor( atom_name, heavy_atom_name );
        } else {
            hbond::identify_acceptor( atom_name, sequence[resid-1], heavy_atom_name );
        }        
        if ( heavy_atom_name == "" ) {
            utility_exit_with_message("Unknown hydrogen bond acceptor/donor type: " + atom_name);
        }


        int heavy_atom_no = find_heavy_atom( resid, heavy_atom_name );
        if ( heavy_atom_no == -1 ) continue;

        hbond::HeavyAtom const & ha = target_heavy_atoms_[heavy_atom_no];
        bool fully_sat = true;
        int satno = -1;

        if ( ha.sat_groups.size() == 1 ) {
            int found_sat = ha.sat_groups.front();
            if ( target_presatisfied_[found_sat] ) fully_sat = true;
            else satno = found_sat;
        } else if ( is_donor ) {
            // if it's a donor, we just check hydrogen distances

            int found_sat = -1;
            for ( int sat : ha.sat_groups ) {
                if ( sat >= num_donors_ ) continue;
                if ( target_presatisfied_[sat] ) continue;
                fully_sat = false;

                HBondRay const & ray = target_donors_acceptors_[sat];
                const float dist = ( ray.horb_cen - xyz ).norm();
                if ( dist > 0.1 ) continue;

                if ( found_sat != -1 ) std::cout << "Error: Overlapping hydrogens??? " << std::endl;
                found_sat = sat;
            }
            satno = found_sat;

        } else {
            // if it's an acceptor, we need to check directionality.
            //  Use closest by angle, but use next closest if that's already satisfied

            Vector3f direction = xyz - last_xyz;

            float lowest_angle = 200;
            int lowest_sat = -1;

            for ( int sat : ha.sat_groups ) {
                if ( sat < num_donors_ ) continue;
                if ( target_presatisfied_[sat] ) continue;
                fully_sat = false;

                const float angle_deg = angle_between_vectors( direction, target_donors_acceptors_[sat].direction ) * M_PI / 180;
                if ( angle_deg < lowest_angle ) {
                    lowest_angle = angle_deg;
                    lowest_sat = sat;
                }
            }

            satno = lowest_sat;
            if (lowest_angle > 120 && lowest_sat != -1) {
                std::cout << "High angle acceptor sat: " << lowest_angle << " " << resid << " " << atom_name << std::endl;
            }

        }

        if ( satno == -1 ) {
            if ( fully_sat ) std::cout << "Super satisfied atom: " << resid << " " << atom_name << std::endl;
            else std::cout << "Error!!! Could not find matching orbital: " << resid << " " << atom_name << std::endl;
            continue;
        }

        target_presatisfied_[satno] = true;

    }

    int presats = 0;
    for ( int i = 0; i < target_presatisfied_.size(); i++ ) {
        if ( target_presatisfied_[i] ) presats ++;
    }
    std::cout << "Found " << presats << " presatisfied donors/acceptors" << std::endl;

}

void
UnsatManager::dump_presatisfied() {

    std::vector<HBondRay> donor_dump;

    for ( int i = 0; i < num_donors_; i++ ) {
        if ( target_presatisfied_[i] ) {
            donor_dump.push_back(target_donors_acceptors_[i]);
        }
    }

    std::vector<HBondRay> acceptor_dump;

    for ( int i = num_donors_; i < target_donors_acceptors_.size(); i++ ) {
        if ( target_presatisfied_[i] ) {
            acceptor_dump.push_back(target_donors_acceptors_[i]);
        }
    }

    utility::io::ozstream donout("presatisfied_donors.pdb.gz");
    ::devel::scheme::dump_hbond_rays( donout, donor_dump, true );
    donout.close();

    utility::io::ozstream accout("presatisfied_acceptors.pdb.gz");
    ::devel::scheme::dump_hbond_rays( accout, acceptor_dump, false );
    accout.close();


}


bool
UnsatManager::patch_heavy_atoms( int resid, std::string const & heavy_atom, hbond::Patch patch, shared_ptr<BurialManager> const & burial_manager ) {

    int heavy_atom_no = find_heavy_atom( resid, heavy_atom, false );
    if ( heavy_atom_no == -1 ) return false;
    
    switch ( patch ) {
        case hbond::ACTUALLY_HBONDING: {
            target_heavy_atoms_.erase( target_heavy_atoms_.begin() + heavy_atom_no );
            if ( store_common_unsats_ ) {
                unsat_counts_.erase( unsat_counts_.begin() + heavy_atom_no );
            }
            int remaining = burial_manager->remove_heavy_atom( heavy_atom_no );
            runtime_assert( target_heavy_atoms_.size() == remaining );
            break;
        }
        case hbond::NOT_BURIED: {
            burial_manager->unbury_heavy_atom( heavy_atom_no );
            break;
        }
    }

    return true;

}

void
UnsatManager::apply_unsat_helper( 
    std::string const & unsat_helper_fname,
    shared_ptr<BurialManager> const & burial_manager
) {
    if ( unsat_helper_fname == "" ) return;

    runtime_assert_msg(utility::file::file_exists( unsat_helper_fname ), "unsat helper file does not exits: " + unsat_helper_fname );

    std::ifstream in;
    in.open(unsat_helper_fname, std::ios::in);

    std::string s;

    while ( std::getline(in, s)){
        if (s.empty()) continue;

        utility::vector1<std::string> splt = utility::quoted_split( s );

        if ( splt[1].find("#") == 0 ) continue;

        runtime_assert_msg( splt.size() >= 3, "Bad line in unsat helper file: " + s );

        int resid = stoi( splt[1] );
        std::string heavy_atom = splt[2];
        std::string flag = splt[3];

        hbond::Patch patch = hbond::NONE;
        if ( flag == "HB" ) patch = hbond::ACTUALLY_HBONDING;
        if ( flag == "NB" ) patch = hbond::NOT_BURIED;
        if ( patch == hbond::NONE ) {
            utility_exit_with_message("Invalid patch operation: " + s );
        }
            std::cout << resid << " " << heavy_atom << " " << patch << std::endl;

        if ( ! patch_heavy_atoms( resid, heavy_atom, patch, burial_manager ) ) {
            utility_exit_with_message( "Heavy atom not found: " + s );
        }

    }
}


std::vector<float>
UnsatManager::get_initial_unsats( 
    std::vector<float> const & burial_weights
) const {
    std::vector< std::pair<intRot,intRot> > rotamers;
    std::vector< EigenXform > bb_positions;
    RifScoreRotamerVsTarget rot_tgt_scorer;

    return get_buried_unsats( burial_weights, rotamers, bb_positions, rot_tgt_scorer );
}

std::vector<float>
UnsatManager::get_buried_unsats( 
    std::vector<float> const & burial_weights,
    std::vector< std::pair<intRot,intRot> > const & rotamers,
    std::vector< EigenXform > const & bb_positions,
    RifScoreRotamerVsTarget const & rot_tgt_scorer
) const {

    std::vector<bool> satisfied = target_presatisfied_;

    // runtime_assert( rotamers.size() == bb_positions.size() );   // this is false

    for ( int i = 0; i < rotamers.size(); i++ ) {
        int ires = rotamers[i].first;
        int irot = rotamers[i].second;

        int sat1 = -1, sat2 = -1, hbcount = 0;
        rot_tgt_scorer.score_rotamer_v_target_sat( irot, bb_positions.at(ires), sat1, sat2, true, hbcount, 10.0, 4 );

        if ( debug_ ) {
            std::cout << rot_index_p->oneletter(irot)  << " rotamer: " << i << " ires: " << ires << " irot: " << irot 
                << " sats: " << sat1 << " " << sat2 << std::endl;
        }

        if ( sat1 > -1 ) satisfied[sat1] = true;
        if ( sat2 > -1 ) satisfied[sat2] = true;
    }

    std::vector<float> unsat_scores( target_heavy_atoms_.size() );
    for ( int iheavy = 0; iheavy < target_heavy_atoms_.size(); iheavy++ ) {

        const float weight = burial_weights[iheavy];
        if ( weight == 0) continue;

        hbond::HeavyAtom const & ha = target_heavy_atoms_[iheavy];
        int sat_count = 0;
        for ( int isat : ha.sat_groups ) {
            if ( satisfied[isat] ) {
                sat_count += 1;
            }
        }

        const int wants = hbond::NumOrbitals[ha.AType];
        std::vector<float> const & penalties = unsat_penalties_[ha.AType];
        float score = calculate_unsat_score( penalties[0], penalties[1], wants, sat_count, weight );
        unsat_scores[iheavy] = score * weight * unsat_score_scalar_;

        if ( debug_ ) {
            std::cout << "BuriedHeavy: " << iheavy << " res: " << ha.resid << " name: " << ha.name 
                          << " type: " << hbond::ATypeNames[ha.AType] << " burial: " << weight 
                          << " wants: " << wants << " has: " << sat_count << " raw_unsat_score: " << score << std::endl;
        }
    }

    return unsat_scores;
}

void
UnsatManager::print_buried_unsats( std::vector<float> const & unsat_penalties) const {

    for ( int iheavy = 0; iheavy < target_heavy_atoms_.size(); iheavy++ ) {

        float score = unsat_penalties[iheavy];
         if ( score > 0 ) {

            hbond::HeavyAtom const & ha = target_heavy_atoms_[iheavy];
            std::cout << " Buried unsat: " << boost::str(boost::format(" resid: %i name: %s score: %.3f ")%ha.resid%ha.name%score);
            std::cout << " satnos: ";
            for ( int isat : ha.sat_groups ) {
                std::cout << isat << ",";
            }
            std::cout << std::endl;
         }
    }

}

void
UnsatManager::print_unsats_help( std::vector<float> const & unsat_penalties) const {

    bool any = false;
    for ( float penalty : unsat_penalties ) {
        any |= penalty > 0;
    }

    if ( !any ) return;

    std::cout << "These initial unsats are going to cause problems with scoring." << std::endl;
    std::cout << "You need to create a -unsat_helper file to address these." << std::endl;
    std::cout << "HB - this atom is actually making an H-bond but rosetta missed it." << std::endl;
    std::cout << "NB - this atom isn't actually buried (but could be)." << std::endl;
    std::cout << "Template:" << std::endl << std::endl;

    for ( int iheavy = 0; iheavy < target_heavy_atoms_.size(); iheavy++ ) {

        float score = unsat_penalties[iheavy];
         if ( score > 0 ) {

            hbond::HeavyAtom const & ha = target_heavy_atoms_[iheavy];
            std::cout << boost::str(boost::format("%i %s HB/NB ")%ha.resid%ha.name) <<  std::endl;
         }
    }

    std::cout << std::endl;

}


void
UnsatManager::sum_unsat_counts( UnsatManager const & other ) {
    runtime_assert( other.unsat_counts_.size() == unsat_counts_.size() );

    for ( int i = 0; i < unsat_counts_.size(); i++ ) {
        unsat_counts_[i] += other.unsat_counts_[i];
    }
}

void
UnsatManager::print_unsat_counts() const {
    using ObjexxFCL::format::I;

    runtime_assert( unsat_counts_.size() == target_heavy_atoms_.size() );

    std::cout << std::endl;
    std::cout << "Commonly buried unsatisfied polars" << std::endl;

    uint64_t max = 0;
    for ( int i = 0; i < unsat_counts_.size(); i++ ) max = std::max( max, unsat_counts_[i] );

    const int width = 80;

    uint64_t unsats_per_star = max / width;

    for ( int iheavy = 0; iheavy < unsat_counts_.size(); iheavy++ ) {

        uint64_t count = unsat_counts_[iheavy];
        if ( count == 0 ) continue;

        hbond::HeavyAtom const & ha = target_heavy_atoms_[iheavy];

        // 13
        //                                                 6  + 4 + 7  + 4  + 7  + 13 = 41
        std::string to_print = boost::str(boost::format("resid:%4i name: %s type: %s")%ha.resid%ha.name%hbond::ATypeNames[ha.AType]);

        while ( to_print.length() < 41 + 2 ) {
            to_print = to_print + " ";
        }

        to_print += I(12, count );

        uint64_t stars = count / unsats_per_star;

        for ( uint64_t i = 0; i < stars; i++ ) {
            to_print += "*";
        }

        std::cout << to_print << std::endl;

    }


}





std::vector<Eigen::Vector3f>
UnsatManager::get_heavy_atom_xyzs() {
    std::vector<Eigen::Vector3f> xyzs( target_heavy_atoms_.size() );
    for ( int i = 0; i < target_heavy_atoms_.size(); i++ ) {
        xyzs[i] = target_heavy_atoms_[i].xyz;
    }
    return xyzs;
}


// I encourage you to check my math. But this calculates the score based on the P0 P1 scheme.
float
UnsatManager::calculate_unsat_score( float P0, float P1, int wants, int satisfied, float weight ) const {

    float penalty = P0;
    float score = 0;
    const float subtracter = P0 * ( 1.0f - P1 );

    while ( wants > 0 ) {
        satisfied --;
        if ( satisfied < 0 ) score += penalty * weight * unsat_score_scalar_;
        penalty -= subtracter;
        wants --;
    }
    return score;
}


float
UnsatManager::calculate_nonpack_score( 
    std::vector<float> const & burial_weights,
    std::vector<bool> const & is_satisfied
) {
    runtime_assert( burial_weights.size() == target_heavy_atoms_.size() );
    runtime_assert( is_satisfied.size() == target_donors_acceptors_.size() );

    float score = 0;
    score += score_offset_;

    int buried = 0;

    for ( int iheavy = 0; iheavy < target_heavy_atoms_.size(); iheavy++ ) {

        if ( burial_weights[iheavy] == 0 ) continue;
        buried ++;
        float weight = burial_weights[iheavy];

        hbond::HeavyAtom const & ha = target_heavy_atoms_[iheavy];

        const int wants = hbond::NumOrbitals[ha.AType];
        int satisfied = 0;
        for ( int sat : ha.sat_groups ) {
            if ( is_satisfied[sat] ) satisfied ++;
        }

        if (satisfied == wants) continue;

        std::vector<float> const & penalties = unsat_penalties_[ha.AType];

        const float adding = calculate_unsat_score( penalties[0], penalties[1], wants, satisfied, weight );
        score += adding;

        if ( debug_ ) std::cout << "iheavy: " << iheavy << " adding: " << adding << " score: " << score << std::endl;

        if ( store_common_unsats_ && adding > 0 ) {
            unsat_counts_[iheavy] ++;
        }
    }

    if (buried < require_burial_) {
        score = 9e9;
    }


    return score;
}

void
UnsatManager::add_to_pack_rot( 
    int ires,
    int irot,
    float score,
    int sat1,
    int sat2
) {
    to_pack_rots_.push_back(ToPackRot(ires, irot, score, sat1, sat2));
}


// If only the target has unsats
// 
// Believe it or not, but it's possible to add the buried unsat penalty as a pairwise decomposable term
//                     AB -- acceptors A and B
//                     H
// Example: Amide    N
//                     H
//                     CD -- acceptors C and D
// 2 Unsats = 6 kcal
// 1 Unsat  = 2 kcal
// 0 Unsat  = 0 kcal
//
// Penalty definition: {4, 0.5}
//
// Step 1:
//   Penalize this buried amide by adding 6 to the total score (function return value)
//
// Step 2:
//   Give a bonus to any rotamer that satisfies the first unsat ( add -4 to the onebody of A, B, C, D )
//
// Step 3:
//   Penalize rotamers that both satisfy the same orbital ( add +4 to the twobody between A-B and C-D )
//
// Step 4:
//   Penalize rotamers that satisfy the different orbitals ( add +2 to the twobody between A-C, A-D, B-C, B-D )
//
// And this will magically apply the penalty definition
//
// Calculated scores:
// 0     6
// A     2
// B     2
// C     2
// D     2
// AB    2
// AC    0
// CD    2
// ABC   2  // not ideal but oh well
// ABCD  6  // not ideal but oh well


// Stages for this calculation
//
// 1. Identify all buried heavy atoms
// 2.   Assign max penalty to each heavy atom
// 3. Find all orbitals of said heavy atoms
// 4. Identify all of their satisfiers (including scaffold backbone and target)
// 5.   Assign P0 as a bonus to all satisfiers
// 6. For all pairs of satisfiers at one orbital
// 7.   Assign P0 as a twobody penalty
// 8. For all pairs on one heavy atom
// 9.   Assign P0 * ( 1 - P1 ) as a twobody penalty

// 10. Insert rotamers into HackPack

// Assignment of onebody energies
//   If object is a rotamer, add to score
//   If object is scaffold bb or target, add to total score

// Assignment of twobody energies
//   If both are rotamers, add to twobody
//   If both are same rotamer, add to onebody
//   If one is rotamer and other is scaffold bb or target, add to rotamer onebody
//   If both are target, add to total score

float
UnsatManager::prepare_packer( 
    ::scheme::search::HackPack & packer, 
    std::vector<float> const & burial_weights,
    std::vector<bool> const & pre_and_bb_satisfied
) {
    runtime_assert( burial_weights.size() == target_heavy_atoms_.size() );
    runtime_assert( pre_and_bb_satisfied.size() == target_presatisfied_.size());

    if (debug_) {

        for ( int ih = 0; ih < burial_weights.size(); ih++ ) {

            const float weight = burial_weights[ih];
            hbond::HeavyAtom const & ha = target_heavy_atoms_[ih];

            std::cout << "Heavy atom: " << ih << " resid: " << ha.resid << " heavy: " << ha.name << " type: " << hbond::ATypeNames[ha.AType] 
                << " burial: " << weight << " sats: ";

            for ( int sat : ha.sat_groups ) {
                std::cout << sat << ",";
            }
            std::cout << std::endl;

        }



        for ( int i = 0; i < to_pack_rots_.size(); i++ ) {
            std::cout << "ToPackRot: " << i << " " << rot_index_p->oneletter(to_pack_rots_[i].irot) 
                << " ires: " << to_pack_rots_[i].ires << " irot: " << to_pack_rots_[i].irot 
                << " score: " << to_pack_rots_[i].score 
                << " sat1: " << to_pack_rots_[i].sat1 << " sat2: " << to_pack_rots_[i].sat2 
                << std::endl;
        }
    }

// data structures

    float zerobody_penalty = 0;
    zerobody_penalty += score_offset_;

    // satisfiers are the index of a rotamer in to_pack_rots_
    // -1 is the target
    std::vector<std::vector<std::vector<int>>> per_heavy_per_orb_satisfiers( burial_weights.size() );

    std::vector<int> heavy_atom_per_sat( target_donors_acceptors_.size(), -1 );


////////////////////////////////////////////////////////////////////////


    for ( int ih = 0; ih < burial_weights.size(); ih++ ) {
        const float weight = burial_weights[ih];

    // 1. Identify all buried heavy atoms
        if ( weight == 0 ) continue;

        hbond::HeavyAtom const & ha = target_heavy_atoms_[ih];

    // 2.   Assign max penalty to each heavy atom
        zerobody_penalty += total_first_twob_[ha.AType][0] * weight * unsat_score_scalar_;

        if (debug_) std::cout << "0body unsat: " << total_first_twob_[ha.AType][0] * weight * unsat_score_scalar_ 
                    << " heavy atom: " << ih << std::endl;

    // 3. Find all orbitals of said heavy atoms
        for ( int sat : ha.sat_groups ) heavy_atom_per_sat[sat] = ih;

    }

        


    // 4. Identify all of their satisfiers (including scaffold backbone and target)
    // 5.   Assign P0 as a bonus to all satisfiers

    std::vector<std::vector<int>> sat_satsifiers( target_donors_acceptors_.size() );

    //       first find target presatisfiers

    for ( int isat = 0; isat < pre_and_bb_satisfied.size(); isat++ ) {
        if ( pre_and_bb_satisfied[isat] ) {
            sat_satsifiers[isat].push_back(-1);

            const int heavy_atom_no = heavy_atom_per_sat[isat];
            if ( heavy_atom_no > -1 ) {
                const int heavy_atom_type = target_heavy_atoms_[heavy_atom_no].AType;
                const float weight = burial_weights[heavy_atom_no];
                zerobody_penalty += - total_first_twob_[heavy_atom_type][1] * weight * unsat_score_scalar_;

                if (debug_) std::cout << "0body self-satisfy: " << - total_first_twob_[heavy_atom_type][1] * weight * unsat_score_scalar_
                                   << " heavy atom: " << heavy_atom_no << std::endl;
            }
        }
    }

    //       second find rotamer satisfiers

    for ( int ipack = 0; ipack < to_pack_rots_.size(); ipack++ ) {
        ToPackRot & to_pack_rot = to_pack_rots_[ipack];

        if ( to_pack_rot.sat1 > -1 ) {
            sat_satsifiers[to_pack_rot.sat1].push_back(ipack);

            const int heavy_atom_no = heavy_atom_per_sat[to_pack_rot.sat1];
            if ( heavy_atom_no > -1 ) {
                const int heavy_atom_type = target_heavy_atoms_[heavy_atom_no].AType;
                const float weight = burial_weights[heavy_atom_no];
                to_pack_rot.score += - total_first_twob_[heavy_atom_type][1] * weight * unsat_score_scalar_;

                if (debug_) std::cout << "1body rot-satisfy: " << - total_first_twob_[heavy_atom_type][1] * weight * unsat_score_scalar_
                                  << " heavy atom: " << heavy_atom_no << " ipack: " << ipack << std::endl;
            }
        }

        if ( to_pack_rot.sat2 > -1 ) {
            sat_satsifiers[to_pack_rot.sat2].push_back(ipack);

            const int heavy_atom_no = heavy_atom_per_sat[to_pack_rot.sat2];
            if ( heavy_atom_no > -1 ) {
                const int heavy_atom_type = target_heavy_atoms_[heavy_atom_no].AType;
                const float weight = burial_weights[heavy_atom_no];
                to_pack_rot.score += - total_first_twob_[heavy_atom_type][1] * weight * unsat_score_scalar_;

                if (debug_) std::cout << "1body rot-satisfy: " << - total_first_twob_[heavy_atom_type][1] * weight * unsat_score_scalar_
                                      << " heavy atom: " << heavy_atom_no << " ipack: " << ipack << std::endl;
            }
        }
    }



    for ( int ih = 0; ih < burial_weights.size(); ih++ ) {
        const float weight = burial_weights[ih];

        if ( weight == 0 ) continue;

        hbond::HeavyAtom const & ha = target_heavy_atoms_[ih];

    // 6. For all pairs of satisfiers at one orbital
    // 7.   Assign P0 as a twobody penalty

        const float P0 = total_first_twob_[ha.AType][1] * weight * unsat_score_scalar_;

        for ( int isat : ha.sat_groups ) {
            std::vector<int> const & local_sat_satsifiers = sat_satsifiers[isat];

            // upper triangle for loop
            for ( int ilocal = 0; ilocal < (int)local_sat_satsifiers.size() - 1; ilocal ++ ) {
                for ( int jlocal = ilocal + 1; jlocal < local_sat_satsifiers.size(); jlocal++ ) {
                    zerobody_penalty += handle_twobody( local_sat_satsifiers[ilocal], local_sat_satsifiers[jlocal], P0, packer );


                    if (debug_) std::cout << "orbital clash: " << P0 << " heavy_atom: " << ih
                                << " ipack1: " << local_sat_satsifiers[ilocal] << " ipack2: " << local_sat_satsifiers[jlocal] << std::endl;
                }
            }
        }



    // 8. For all pairs on one heavy atom
    // 9.   Assign P0 * ( 1 - P1 ) as a twobody penalty


        const float twob_penalty = total_first_twob_[ha.AType][2] * weight * unsat_score_scalar_;

        // upper triangle of orbitals
        for ( int iorb = 0; iorb < (int)ha.sat_groups.size() - 1; iorb ++) {

            const int isat = ha.sat_groups[iorb];
            std::vector<int> const & iorb_satisfiers = sat_satsifiers[isat];

            for ( int jorb = iorb + 1; jorb < ha.sat_groups.size(); jorb++ ) {

                const int jsat = ha.sat_groups[jorb];
                std::vector<int> const & jorb_satisfiers = sat_satsifiers[jsat];

                // all x all of iorb_satisfiers against jorb_satisfiers
                for ( int ilocal = 0; ilocal < iorb_satisfiers.size(); ilocal++ ) {
                    for ( int jlocal = 0; jlocal < jorb_satisfiers.size(); jlocal++ ) {
                        zerobody_penalty += handle_twobody( iorb_satisfiers[ilocal], jorb_satisfiers[jlocal], twob_penalty, packer );

                        if (debug_) std::cout << "heavy atom clash: " << twob_penalty << " heavy_atom: " << ih
                                    << " ipack1: " << iorb_satisfiers[ilocal] << " ipack2: " << jorb_satisfiers[jlocal] << std::endl;
                    }
                }




            }
        }


    }
    insert_to_pack_rots_into_packer(packer);

    return zerobody_penalty;

}

void
UnsatManager::insert_to_pack_rots_into_packer(::scheme::search::HackPack & packer) {
    if ( debug_ ) std::cout << "Inserting into packer:" << std::endl;
    for ( int i = 0; i < to_pack_rots_.size(); i++ ) {
        ToPackRot const & rot = to_pack_rots_[i];
        packer.add_tmp_rot( rot.ires, rot.irot, rot.score );

        if ( debug_ ) {
            std::cout << "ToPackRot: " << i << " " << rot_index_p->oneletter(rot.irot) 
                << " ires: " << rot.ires << " irot: " << rot.irot 
                << " score: " << rot.score 
                << " sat1: " << rot.sat1 << " sat2: " << rot.sat2 
                << std::endl;
        }
    }
}

float
UnsatManager::handle_twobody( int satisfier1, int satisfier2, float penalty, ::scheme::search::HackPack & packer ) {

    if ( satisfier1 == -1 && satisfier2 == -1 ) return penalty;

    if ( satisfier1 == -1 ) {
        to_pack_rots_[ satisfier2 ].score += penalty;
        return 0;
    }

    if ( satisfier2 == -1 ) {
        to_pack_rots_[ satisfier1 ].score += penalty;
        return 0;
    }

    if ( satisfier1 == satisfier2 ) {
        to_pack_rots_[ satisfier1 ].score += penalty;
        return 0;
    }

    ToPackRot const & pack1 = to_pack_rots_[ satisfier1 ];
    ToPackRot const & pack2 = to_pack_rots_[ satisfier2 ];

    packer.twob_->upweight_edge( pack1.ires, pack2.ires, pack1.irot, pack2.irot, penalty );
    modified_edges_.push_back( std::array<int, 4> { pack1.ires, pack2.ires, pack1.irot, pack2.irot } );

    return 0;

}

void
UnsatManager::fix_packer( 
    ::scheme::search::HackPack & packer, 
    shared_ptr<::scheme::objective::storage::TwoBodyTable<float> const> reference_twobody
) {
    for ( std::array<int, 4> const & arr : modified_edges_ ) {
        packer.twob_->restore_edge( arr[0], arr[1], arr[2], arr[3], reference_twobody );
    }
}



}}
