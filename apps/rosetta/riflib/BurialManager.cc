

#include <riflib/BurialManager.hh>


namespace devel {
namespace scheme {




void
BurialManager::set_target_neighbors( core::pose::Pose const & pose ) {

    float cutoff_sq = opts_.neighbor_distance_cutoff*opts_.neighbor_distance_cutoff;

    for ( int i_res = 1; i_res <= pose.size(); i_res++ ) {
        numeric::xyzVector<core::Real> xyz = pose.residue(i_res).nbr_atom_xyz();
        Eigen::Vector3f CB_xyz;
        CB_xyz[0] = xyz[0]; CB_xyz[1] = xyz[1]; CB_xyz[2] = xyz[2];

        for ( int i_ray = 0; i_ray < donor_acceptors_.size(); i_ray++ ) {
            const float dist_sq = (donor_acceptors_[i_ray].horb_cen - CB_xyz).squaredNorm();
            if ( dist_sq > cutoff_sq ) continue;
            target_neighbor_counts_[i_ray]++;
        }
    }

}


void
BurialManager::get_burial_weights( std::vector<float> & weights ) const {
    runtime_assert( weights.size() == target_neighbor_counts_.size());

    for ( int i_ray = 0; i_ray < target_neighbor_counts_.size(); i_ray++ ) {
        weights[i_ray] = opts_.neighbor_count_weights[ target_neighbor_counts_[i_ray] 
                                                    + other_neighbor_counts_[i_ray]];
    }
}

void
BurialManager::accumulate_neighbors( BBActor const & bb ) {

    float cutoff_sq = opts_.neighbor_distance_cutoff*opts_.neighbor_distance_cutoff;

    Eigen::Vector3f CB_xyz = bb.position().translation();

    for ( int i_ray = 0; i_ray < donor_acceptors_.size(); i_ray++ ) {
        const float dist_sq = (donor_acceptors_[i_ray].horb_cen - CB_xyz).squaredNorm();
        if ( dist_sq > cutoff_sq ) continue;
        other_neighbor_counts_[i_ray]++;
    }
}







}}
