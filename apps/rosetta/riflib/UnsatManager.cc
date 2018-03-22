

#include <riflib/UnsatManager.hh>

#include <riflib/RotamerGenerator.hh>
#include <riflib/util.hh>

#include <utility/io/ozstream.hh>

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

            float dist = ( xyz_ - ( ray.horb_cen + ray.direction ) ).norm();

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

            float dist = ( xyz_ - ( ray.horb_cen + ray.direction*(1+ORBLEN) ) ).norm();

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

    runtime_assert( validate_heavy_atoms() );

}

int
UnsatManager::find_heavy_atom( int resid, std::string const & heavy_atom_name ) {
    int heavy_atom_no = 0;
    for ( heavy_atom_no = 0; heavy_atom_no < target_heavy_atoms_.size(); heavy_atom_no++ ) {
        hbond::HeavyAtom const & ha = target_heavy_atoms_[heavy_atom_no];
        if ( ha.resid != resid ) continue;
        if ( ha.name != heavy_atom_name ) continue;
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
    std::vector<std::pair<int,int>> presat_res_atoms = ::devel::scheme::get_satisfied_atoms_vector(target, -0.1f);

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
                const float dist = ( (ray.horb_cen + ray.direction) - xyz ).norm();
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

std::vector<Eigen::Vector3f>
UnsatManager::get_heavy_atom_xyzs() {
    std::vector<Eigen::Vector3f> xyzs( target_heavy_atoms_.size() );
    for ( int i = 0; i < target_heavy_atoms_.size(); i++ ) {
        xyzs[i] = target_heavy_atoms_[i].xyz;
    }
    return xyzs;
}





}}
