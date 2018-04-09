

#include <riflib/BurialManager.hh>

#include <ObjexxFCL/format.hh>

#include <fstream>

namespace devel {
namespace scheme {



shared_ptr<BurialManager>
BurialManager::clone() const {
    shared_ptr<BurialManager> ot = make_shared<BurialManager>();
    ot->opts_ = opts_;
    // ot->num_donors_ = num_donors_;
    // ot->donor_acceptors_ = donor_acceptors_;
    ot->target_burial_points_ = target_burial_points_;
    ot->target_neighbor_counts_ = target_neighbor_counts_;
    ot->other_neighbor_counts_ = other_neighbor_counts_;
    return ot;
}

void
BurialManager::reset() {
    std::fill(other_neighbor_counts_.begin(), other_neighbor_counts_.end(), 0);
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

    float radius_sq = radius * radius;

    float step = grid.cs_[0];

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


void
BurialManager::set_target_neighbors( core::pose::Pose const & pose ) {

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

    std::cout << "Building target burial grid" << std::endl;


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
        lbs[i] -= /*opts.neighbor_distance_cutoff +*/ 1;
        ubs[i] += /*opts.neighbor_distance_cutoff +*/ 1;
    }

    Eigen::Vector3f resolution( opts_.burial_grid_spacing, opts_.burial_grid_spacing, opts_.burial_grid_spacing );

    std::cout << "Grid bounds: ";
    for ( int i = 0; i < 3; i++ ) {
        if (i > 0) std::cout << " x ";
        std::cout << "( " << F(5,1,lbs[i]) << ", " << F(5,1,ubs[i]) << " )";
    }
    std::cout << std::endl;

    // let's pretend for now that this gets initialized to 0
    target_burial_grid_ = make_shared<BurialVoxelArray>( lbs, ubs, resolution );


    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
        numeric::xyzVector<core::Real> _xyz = pose.residue(resid).nbr_atom_xyz();
        Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

        fill_voxel_near_xyz( *target_burial_grid_, xyz, opts_.neighbor_distance_cutoff );
    }
    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
        numeric::xyzVector<core::Real> _xyz = pose.residue(resid).xyz("N");
        Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

        fill_voxel_near_xyz( *target_burial_grid_, xyz, opts_.neighbor_distance_cutoff );
    }
    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
        numeric::xyzVector<core::Real> _xyz = pose.residue(resid).xyz("CA");
        Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

        fill_voxel_near_xyz( *target_burial_grid_, xyz, opts_.neighbor_distance_cutoff );
    }
    for ( core::Size resid = 1; resid <= pose.size(); resid ++ ) {
        numeric::xyzVector<core::Real> _xyz = pose.residue(resid).xyz("C");
        Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

        fill_voxel_near_xyz( *target_burial_grid_, xyz, opts_.neighbor_distance_cutoff );
    }



}


void
BurialManager::dump_burial_grid( std::string const & fname ) {


    std::cout << "Dumping burial grid to " + fname << std::endl;
    std::ofstream out( fname );

    typedef typename BurialVoxelArray::Bounds Bounds;

    BurialVoxelArray & grid = *target_burial_grid_;

    Bounds lb = grid.lb_;
    Bounds ub = grid.ub_;

    float step = grid.cs_[0];

    int anum=0, rnum = 0;

    Eigen::Vector3f worker;
    for ( float x = lb[0] + step/2; x < ub[0]; x += step ) { worker[0] = x;
    for ( float y = lb[1] + step/2; y < ub[1]; y += step ) { worker[1] = y;
    for ( float z = lb[2] + step/2; z < ub[2]; z += step ) { worker[2] = z;

        int value = int( grid.at(worker) + 0.1f );
        char buf[128];

        if ( opts_.neighbor_count_weights[value] ) {
            std::string aname;
            snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
                "HETATM",
                ++anum,
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



std::vector<float>
BurialManager::get_burial_weights( ) const {

    std::vector<float> weights( target_neighbor_counts_.size() );
    for ( int i_pt = 0; i_pt < target_neighbor_counts_.size(); i_pt++ ) {
        weights[i_pt] = opts_.neighbor_count_weights[ target_neighbor_counts_[i_pt] 
                                                    + other_neighbor_counts_[i_pt]];
    }
    return weights;
}

void
BurialManager::accumulate_neighbors( BBActor const & bb ) {

    float cutoff_sq = opts_.neighbor_distance_cutoff*opts_.neighbor_distance_cutoff;

    Eigen::Vector3f CB_xyz = bb.position().translation();

    for ( int i_pt = 0; i_pt < target_burial_points_.size(); i_pt++ ) {
        const float dist_sq = (target_burial_points_[i_pt] - CB_xyz).squaredNorm();
        if ( dist_sq > cutoff_sq ) continue;
        other_neighbor_counts_[i_pt]++;
    }
}






}}
