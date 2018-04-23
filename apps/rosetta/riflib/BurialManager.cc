

#include <riflib/BurialManager.hh>

#include <ObjexxFCL/format.hh>

#include <fstream>

namespace devel {
namespace scheme {



shared_ptr<BurialManager>
BurialManager::clone() const {
    shared_ptr<BurialManager> ot = make_shared<BurialManager>();
    *ot = *this;
    return ot;
}

void
BurialManager::reset() {
    // std::fill(other_neighbor_counts_.begin(), other_neighbor_counts_.end(), 0);
}


void
fill_voxel_near_xyz( BurialVoxelArray & grid, Eigen::Vector3f const & xyz, float radius ) {


    Eigen::Vector3f lbs( xyz[0] - radius, xyz[1] - radius, xyz[2] - radius );
    Eigen::Vector3f ubs( xyz[0] + radius, xyz[1] + radius, xyz[2] + radius  );

    bool overlap = true;

    for ( int i = 0; i < 3; i++ ) {
        if ( lbs[i] > grid.ub_[i] ) overlap = false;
        if ( ubs[i] < grid.lb_[i] ) overlap = false;
    }

    if ( !overlap ) return;

    const float radius_sq = radius * radius;

    const float step = grid.cs_[0];

    Eigen::Vector3f worker;

    for ( float x = lbs[0] + step/2; x < ubs[0]; x += step ) {
        if ( x < grid.lb_[0] || x > grid.ub_[0] ) continue;
        worker[0] = x;

        for ( float y = lbs[1] + step/2; y < ubs[1]; y += step ) {
                if ( y < grid.lb_[1] || y > grid.ub_[1] ) continue;
                worker[1] = y;

            for ( float z = lbs[2] + step/2; z < ubs[2]; z += step ) {
                if ( z < grid.lb_[2] || z > grid.ub_[2] ) continue;
                worker[2] = z;

                if ( ( xyz - worker ).squaredNorm() < radius_sq ) {
                    grid[worker] += 1;
                }
            }
        }
    }    

}

/////// !!!!!!!!!!!!!!!!!!!!!  STOLEN FROM ROSETTA  !!!!!!!!!!!!!!! /////////////////////////
/// @brief Given a point in 3D space, and a vector and floats defining a cone, determine the extent to which the point
/// is in the cone.
/// @details The return value ranges from 0 (not in the cone) to 1 (fully in the cone).  The cone has fuzzy boundaries, so
/// non-integer return values are possible.
/// @param [in] point_coordinates The coordinates of a test point that may or may not be in the cone.  For the layer selector, this is the beta- or alpha-carbon
/// of another residue.
/// @param[in] conevect A vector defining the direction in which the test cone opens.  For the layer selector, this is the CA-CB vector of an amino acid.
/// @param[in] conevect_coordinate_2 The coordinate in space of the base of the cone.  For the layer selector, this is the CB atom.
/// @param[in] angle_exponent A value defining how rapidly the cone falls off in the angular direction.
/// @param[in] angle_shift_factor A value that shifts the angluar falloff.
/// @param[in] dist_exponent A value defining how rapidly the cone falls off with distance.
/// @param[in] dist_midpoint A value defining the midpoint of the distance falloff.
inline core::Real calculate_point_in_cone(
    Eigen::Vector3f const &point_coordinates,
    Eigen::Vector3f const &conevect,
    Eigen::Vector3f const &conevect_coordinate_2,
    float const angle_exponent,
    float const angle_shift_factor,
    float const dist_exponent,
    float const dist_midpoint
) {
    Eigen::Vector3f vect = point_coordinates - conevect_coordinate_2;
    float const dist_term(1.0f / (1.0f + (float)std::exp( dist_exponent*(vect.norm() - dist_midpoint)  ))); 
    core::Real angle_term( ( conevect.dot(vect.normalized()) + angle_shift_factor ) / (1 + angle_shift_factor ) );
    if ( angle_term < 0 ) {
        angle_term = 0.0;
    }
    return (dist_term * std::pow<float>(angle_term, angle_exponent) );
}




void
fill_voxel_near_xyz_cone( BurialVoxelArray & grid, Eigen::Vector3f const & ca, Eigen::Vector3f const & cb, float test1, float test2 ) {


    Eigen::Vector3f conevect = (cb - ca).normalized();

    float radius = 10;

    Eigen::Vector3f lbs( cb[0] - radius, cb[1] - radius, cb[2] - radius );
    Eigen::Vector3f ubs( cb[0] + radius, cb[1] + radius, cb[2] + radius  );

    bool overlap = true;

    for ( int i = 0; i < 3; i++ ) {
        if ( lbs[i] > grid.ub_[i] ) overlap = false;
        if ( ubs[i] < grid.lb_[i] ) overlap = false;
    }

    if ( !overlap ) return;


    const float step = grid.cs_[0];

    Eigen::Vector3f worker;

    for ( float x = lbs[0] + step/2; x < ubs[0]; x += step ) {
        if ( x < grid.lb_[0] || x > grid.ub_[0] ) continue;
        worker[0] = x;

        for ( float y = lbs[1] + step/2; y < ubs[1]; y += step ) {
                if ( y < grid.lb_[1] || y > grid.ub_[1] ) continue;
                worker[1] = y;

            for ( float z = lbs[2] + step/2; z < ubs[2]; z += step ) {
                if ( z < grid.lb_[2] || z > grid.ub_[2] ) continue;
                worker[2] = z;


                float in_cone = calculate_point_in_cone( worker, conevect, ca, 1.0f, -0.4f, 4.0, 4.0f);
                
                if (in_cone < 0.3) {
                    in_cone = 0;
                } else {
                    in_cone = 1;
                }

                grid[worker] += in_cone*test1;
                
            }
        }
    }    

}


void
BurialManager::set_target_neighbors( core::pose::Pose const & pose ) {
    target_burial_grid_ = generate_burial_grid( pose );
}


shared_ptr<BurialVoxelArray>
BurialManager::generate_burial_grid( core::pose::Pose const & pose ) {

    // float cutoff_sq = opts_.neighbor_distance_cutoff*opts_.neighbor_distance_cutoff;

    // for ( int i_res = 1; i_res <= pose.size(); i_res++ ) {
    //     numeric::xyzVector<core::Real> xyz = pose.residue(i_res).nbr_atom_xyz();
    //     Eigen::Vector3f CB_xyz;
    //     CB_xyz[0] = xyz[0]; CB_xyz[1] = xyz[1]; CB_xyz[2] = xyz[2];

    //     for ( int i_pt = 0; i_pt < target_burial_points_.size(); i_pt++ ) {
    //         const float dist_sq = (target_burial_points_[i_pt] - CB_xyz).squaredNorm();
    //         if ( dist_sq > cutoff_sq ) continue;
    //         target_neighbor_counts_[i_pt]++;
    //     }
    // }

    std::cout << "Building burial grid" << std::endl;


    using ObjexxFCL::format::F;

    Eigen::Vector3f lbs( 9e9, 9e9, 9e9 );
    Eigen::Vector3f ubs( -9e9, -9e9, -9e9  );

    std::cout << target_burial_points_.size() << std::endl;
    for ( Eigen::Vector3f const & xyz : target_burial_points_ ) {
        for ( int i = 0; i < 3; i++ ) {
            lbs[i] = std::min<float>( lbs[i], xyz[i] );
            ubs[i] = std::max<float>( ubs[i], xyz[i] );
        }
    }

    for ( int i = 0; i < 3; i++ ) {
        lbs[i] -= /*opts.neighbor_distance_cutoff +*/ 3;
        ubs[i] += /*opts.neighbor_distance_cutoff +*/ 3;
    }

    Eigen::Vector3f resolution( opts_.burial_grid_spacing, opts_.burial_grid_spacing, opts_.burial_grid_spacing );

    std::cout << "Grid bounds: ";
    for ( int i = 0; i < 3; i++ ) {
        if (i > 0) std::cout << " x ";
        std::cout << "( " << F(5,1,lbs[i]) << ", " << F(5,1,ubs[i]) << " )";
    }
    std::cout << std::endl;

    // let's pretend for now that this gets initialized to 0
    shared_ptr<BurialVoxelArray> burial_grid = make_shared<BurialVoxelArray>( lbs, ubs, resolution );


    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
        numeric::xyzVector<core::Real> _xyz = pose.residue(resid).nbr_atom_xyz();
        Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

        fill_voxel_near_xyz( *burial_grid, xyz, opts_.neighbor_distance_cutoff );
    }
    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
        numeric::xyzVector<core::Real> _xyz = pose.residue(resid).xyz("N");
        Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

        fill_voxel_near_xyz( *burial_grid, xyz, opts_.neighbor_distance_cutoff );
    }
    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
        numeric::xyzVector<core::Real> _xyz = pose.residue(resid).xyz("CA");
        Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

        fill_voxel_near_xyz( *burial_grid, xyz, opts_.neighbor_distance_cutoff );
    }
    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
        numeric::xyzVector<core::Real> _xyz = pose.residue(resid).xyz("C");
        Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

        fill_voxel_near_xyz( *burial_grid, xyz, opts_.neighbor_distance_cutoff );
    }
    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
        if ( ! pose.residue(resid).has("CB") ) continue;

        numeric::xyzVector<core::Real> _xyza = pose.residue(resid).xyz("CA");
        Eigen::Vector3f xyza; xyza[0] = _xyza[0]; xyza[1] = _xyza[1]; xyza[2] = _xyza[2];

        numeric::xyzVector<core::Real> _xyzb = pose.residue(resid).xyz("CB");
        Eigen::Vector3f xyzb; xyzb[0] = _xyzb[0]; xyzb[1] = _xyzb[1]; xyzb[2] = _xyzb[2];

        Eigen::Vector3f xyz = xyzb + (xyzb - xyza);

        fill_voxel_near_xyz( *burial_grid, xyz, opts_.neighbor_distance_cutoff );
    }

    // for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
    //     if ( ! pose.residue(resid).has("CB")) continue;

    //     numeric::xyzVector<core::Real> _xyzb = pose.residue(resid).xyz("CB");
    //     Eigen::Vector3f xyzb; xyzb[0] = _xyzb[0]; xyzb[1] = _xyzb[1]; xyzb[2] = _xyzb[2];

    //     numeric::xyzVector<core::Real> _xyza = pose.residue(resid).xyz("CA");
    //     Eigen::Vector3f xyza; xyza[0] = _xyza[0]; xyza[1] = _xyza[1]; xyza[2] = _xyza[2];

    //     fill_voxel_near_xyz_cone( *target_burial_grid_, xyza, xyzb, test1, test2  );
    // }

    return burial_grid;
}



float
burial_lookup( 
    Eigen::Vector3f const & pt,
    shared_ptr<BurialVoxelArray> const & target_grid, 
    EigenXform const & scaff_inv_transform,
    shared_ptr<BurialVoxelArray> const & scaff_grid ) {

    float burial = target_grid->at( pt );

    if ( !scaff_grid ) return burial;

    Eigen::Vector3f scaff_pt = scaff_inv_transform * pt;

    burial += scaff_grid->at( scaff_pt );

    return burial;
}


void
BurialManager::dump_burial_grid( 
    std::string const & fname,  
    EigenXform const & scaff_transform, 
    shared_ptr<BurialVoxelArray> const & scaff_grid 
) {
    runtime_assert( target_burial_grid_ );

    std::cout << "Dumping burial grid to " + fname << std::endl;
    std::ofstream out( fname );

    typedef typename BurialVoxelArray::Bounds Bounds;

    BurialVoxelArray & grid = *target_burial_grid_;
    EigenXform const & scaff_inv_transform = scaff_transform.inverse();

    Bounds lb = grid.lb_;
    Bounds ub = grid.ub_;

    float step = grid.cs_[0];

    int anum=0, rnum = 0;

    Eigen::Vector3f worker;
    for ( float x = lb[0] + step/2; x < ub[0]; x += step ) { worker[0] = x;
    for ( float y = lb[1] + step/2; y < ub[1]; y += step ) { worker[1] = y;
    for ( float z = lb[2] + step/2; z < ub[2]; z += step ) { worker[2] = z;

        float burial = burial_lookup( worker, target_burial_grid_, scaff_inv_transform, scaff_grid );
        int value = int( burial + 0.5f );
        char buf[128];

        if ( opts_.neighbor_count_weights[value] ) {
            std::string aname;
            snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
                "HETATM",
                (++anum)%100000,
                " BUR",
                "BUR",
                'A',
                (++rnum)%10000,
                x,y,z,
                1.0,
                1.0,
                "B"
            );
            out << buf;
        }

    }
    }
    }    

    out.close();


}


float
BurialManager::get_burial_count( 
    Eigen::Vector3f const & xyz,
    EigenXform const & scaff_inv_transform,
    shared_ptr<BurialVoxelArray> const & scaff_grid
) const {
    Eigen::Vector3f work = xyz;

    // std::cout << work << " og" << std::endl;

    float burial_count = 9e9;
    for ( int dim = 0; dim < 3; dim++ ) {
        for ( float dir = -1; dir < 2; dir += 2 ) {
            work[dim] += dir * 2.0f;
            const float this_burial_count = burial_lookup( work, target_burial_grid_, scaff_inv_transform, scaff_grid);
            burial_count = std::min<float>(burial_count, this_burial_count);
            // if (debug_) std::cout << this_burial_count;
            // std::cout << work << std::endl;
            work[dim] -= dir * 2.0f;
        }
    }

    return burial_count;

}

std::vector<float>
BurialManager::get_burial_weights( EigenXform const & scaff_transform, shared_ptr<BurialVoxelArray> const & scaff_grid) const {

    EigenXform const & scaff_inv_transform = scaff_transform.inverse();
    std::vector<float> weights( target_burial_points_.size() );

    for ( int i_pt = 0; i_pt < target_burial_points_.size(); i_pt ++ ) {

        const float burial_count = get_burial_count( target_burial_points_[i_pt], scaff_inv_transform, scaff_grid ) 
                                    - unburial_adjust_[i_pt];

        // if (debug_) std::cout << "Burial: iheavy: " << i_pt << " count: " << burial_count << std::endl;
        const float burial = opts_.neighbor_count_weights[int(burial_count + 0.5f)];
        weights[i_pt] = burial;
    }


    return weights;
}

// void
// BurialManager::accumulate_neighbors( BBActor const & bb ) {

//     float cutoff_sq = opts_.neighbor_distance_cutoff*opts_.neighbor_distance_cutoff;

//     Eigen::Vector3f CB_xyz = bb.position().translation();

//     for ( int i_pt = 0; i_pt < target_burial_points_.size(); i_pt++ ) {
//         const float dist_sq = (target_burial_points_[i_pt] - CB_xyz).squaredNorm();
//         if ( dist_sq > cutoff_sq ) continue;
//         other_neighbor_counts_[i_pt]++;
//     }
// }

int
BurialManager::remove_heavy_atom( int heavy_atom_no ) {
    target_burial_points_.erase( target_burial_points_.begin() + heavy_atom_no );
    unburial_adjust_.erase( unburial_adjust_.begin() + heavy_atom_no );

    runtime_assert( target_burial_points_.size() == unburial_adjust_.size() );
    return target_burial_points_.size();
}


void
BurialManager::unbury_heavy_atom( int heavy_atom_no ) {
    float count = get_burial_count( target_burial_points_[heavy_atom_no], EigenXform::Identity(), nullptr );

    float unburial_amount = 0;

    while ( count - unburial_amount >= 0 && opts_.neighbor_count_weights[int(count - unburial_amount + 0.5f)] > 0 ) {
        unburial_amount += 1;
    }

    if ( count - unburial_amount < 0 ) {
        utility_exit_with_message( "Weird error: You can't use a unsat_helper file if you make all atoms buried.");
    }

    unburial_adjust_[heavy_atom_no] = unburial_amount;

}



}}
