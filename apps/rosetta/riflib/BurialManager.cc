

#include <riflib/BurialManager.hh>


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
BurialManager::set_target_neighbors( core::pose::Pose const & pose ) {

    float cutoff_sq = opts_.neighbor_distance_cutoff*opts_.neighbor_distance_cutoff;

    for ( int i_res = 1; i_res <= pose.size(); i_res++ ) {
        numeric::xyzVector<core::Real> xyz = pose.residue(i_res).nbr_atom_xyz();
        Eigen::Vector3f CB_xyz;
        CB_xyz[0] = xyz[0]; CB_xyz[1] = xyz[1]; CB_xyz[2] = xyz[2];

        for ( int i_pt = 0; i_pt < target_burial_points_.size(); i_pt++ ) {
            const float dist_sq = (target_burial_points_[i_pt] - CB_xyz).squaredNorm();
            if ( dist_sq > cutoff_sq ) continue;
            target_neighbor_counts_[i_pt]++;
        }
    }

}


std::vector<float>
BurialManager::get_burial_weights( ) const {

    std::vector<float> weights( target_neighbor_counts_.size() );
    for ( int i_pt = 0; i_pt < target_neighbor_counts_.size(); i_pt++ ) {
        weights[i_pt] = opts_.neighbor_count_weights[ target_neighbor_counts_[i_pt] 
                                                    + other_neighbor_counts_[i_pt]];
    }
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
