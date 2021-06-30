// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <scheme/numeric/rand_xform.hh>

#include <riflib/util.hh>
#include <riflib/util_complex.hh>
#include <riflib/scaffold/MultithreadPoseCloner.hh>
#include <riflib/ScoreRotamerVsTarget.hh>
#include <riflib/RotamerGenerator.hh>

#include <numeric/agglomerative_hierarchical_clustering.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>

#include <ObjexxFCL/format.hh>


namespace devel {
namespace scheme {


void
all_by_all_rmsd( 
	std::vector<core::pose::PoseOP> const & poses,
	utility::vector1<utility::vector1<core::Real>> & table ) {
	
	table.clear();
	table.resize( poses.size() );

	for ( core::Size i = 1; i <= table.size(); i++ ) {
		table[i].resize(poses.size(), 0);	// initialize the diagonal to 0
	}

	// Make it threadsafe
	utility::vector1<shared_ptr<MultithreadPoseCloner>> mpcs;
	for ( core::pose::PoseOP const & pose : poses ) {
		mpcs.push_back(make_shared<MultithreadPoseCloner>(pose));
	}

	std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1)
    #endif
	for ( core::Size i = 1; i <= mpcs.size(); i++ ) {
		if (exception) continue;
		try {
			core::pose::PoseCOP outer_pose = mpcs[i]->get_pose();

			for ( core::Size j = i + 1; j <= mpcs.size(); j++) {
				core::pose::PoseCOP inner_pose = mpcs[j]->get_pose();

				core::Real rmsd = core::scoring::CA_rmsd(*outer_pose, *inner_pose);

				runtime_assert( table[i][j] == 0 );
				runtime_assert( table[j][i] == 0 );

				table[i][j] = rmsd;
				table[j][i] = rmsd;
			}

		} catch(...) {
            #pragma omp critical
            exception = std::current_exception();
        }

	} // end of OMP loop
    if( exception ) std::rethrow_exception(exception);

}

std::vector<std::vector<std::pair<core::pose::PoseOP, uint64_t>>>
cluster_poses_into_n_bins( 
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n,
	utility::vector1<utility::vector1<core::Real>> & rmsds ) {

	runtime_assert( poses.size() >= n );

	if (rmsds.size() != poses.size()) {
		std::cout << "Calculating n^2 rmsd table" << std::endl;
		all_by_all_rmsd( poses, rmsds );
	}

	numeric::AverageLinkClusterer alc;
	utility::vector1<numeric::ClusteringTreeNodeOP> roots = alc.cluster(rmsds, n);

	utility::vector1<uint64_t> input_indices( poses.size() );
	for ( core::Size i = 1; i <= poses.size(); i++ ) {
		input_indices[i] = i - 1;
	}

	std::vector<bool> pose_got_used(poses.size(), false);

	std::vector<std::vector<std::pair<core::pose::PoseOP, uint64_t>>> bins;

	for ( core::Size i = 1; i <= roots.size(); i++ ) {
		std::vector<std::pair<core::pose::PoseOP, uint64_t>> this_bin;

		utility::vector1<uint64_t> this_bin_indices;
		numeric::get_cluster_data( input_indices, roots[i], this_bin_indices );

		runtime_assert( this_bin_indices.size() > 0 );

		for ( uint64_t index : this_bin_indices ) {
			runtime_assert( ! pose_got_used[index] );
			pose_got_used[index] = true;
			this_bin.push_back( std::pair<core::pose::PoseOP, uint64_t>(poses[index], index) );
		}

		bins.push_back( this_bin );
	}

	for ( uint64_t i = 0; i < pose_got_used.size(); i++ ) {
		runtime_assert(pose_got_used[i]);
	}

	return bins;
}

std::vector<core::pose::PoseOP>
cluster_poses_leaving_n( 
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n ) {

	if ( n >= poses.size() ) return poses;

	utility::vector1<utility::vector1<core::Real>> rmsds;
	std::vector<std::vector<std::pair<core::pose::PoseOP, uint64_t>>> bins = cluster_poses_into_n_bins( poses, n, rmsds );

	std::vector<core::pose::PoseOP> output_poses;

	for (std::vector<std::pair<core::pose::PoseOP, uint64_t>> const & bin : bins) {
		output_poses.push_back(bin.front().first);
	}

	return output_poses;
}

// this returns the index into the indexes list
// indexes should be 0 based even though rmsds is 1 based
size_t
find_cluster_center( 
	std::vector<uint64_t> indexes,
	utility::vector1<utility::vector1<core::Real>> const & rmsds ) {

	std::vector<core::Real> rmsd_sums( indexes.size() );

	for ( size_t i = 0; i < indexes.size(); i++ ) {
		uint64_t index = indexes[i];

		core::Real this_sum = 0;
		for ( uint64_t ind : indexes ) {
			this_sum += rmsds.at(ind+1).at(index+1);
		}
		// std::cout << this_sum << std::endl;
		rmsd_sums[i] = this_sum;
	}
	ptrdiff_t center_index = std::min_element(rmsd_sums.begin(), rmsd_sums.end()) - rmsd_sums.begin();

	std::cout << "Cluster center is " << center_index << " with rmsd sum " << rmsd_sums[center_index] << std::endl;

	return center_index;
}

// uses a binary search to try to find a number
// of bins to cluster poses into such than the top n
// represent frac of poses within tol
std::vector<core::pose::PoseOP>
cluster_poses_leaving_n_representing_frac(
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n,
	float frac,
	float tol
	) {

	using ObjexxFCL::format::I;

	runtime_assert( frac > 0 && frac <= 1);

	uint64_t trial = n * (uint64_t)( 1.0 / frac );
	uint64_t initial_trial = trial;
	std::vector<uint64_t> trial_history;

	std::vector<size_t> idx ;
	utility::vector1<utility::vector1<core::Real>> rmsds;

	std::vector<std::vector<std::pair<core::pose::PoseOP, uint64_t>>> bins;

	while (true) {

		trial_history.push_back( trial );
		std::cout << "Round " << trial_history.size() << " : trial size " << trial << std::endl;

		bins = cluster_poses_into_n_bins( poses, trial, rmsds );

// sort indexes based on size of bin
		idx.clear();
		idx.resize(bins.size());
  		std::iota(idx.begin(), idx.end(), 0);
		std::sort(idx.begin(), idx.end(),
       		[&bins](size_t i1, size_t i2) {return bins[i1].size() < bins[i2].size();});


		size_t biggest_bin = bins[idx.back()].size();

// display block //////////////////////////
		{
			// find the number of digits
			int digits = 0; 
			{
				size_t temp = biggest_bin;
				while (temp != 0) { temp /= 10; digits++; }
			}

			size_t star_every = std::max<size_t>(biggest_bin / 50, 1);

			for ( size_t i = 0; i < bins.size(); i++ ) {
				if ( i == n ) std::cout << "============================ cut =========================" << std::endl;
				// going in reverse order and using indexes
				size_t bin = idx[ bins.size() - i - 1 ];

				std::cout << I( digits + 1, bins[bin].size() ) << " ";

				{
					int temp = bins[bin].size() - star_every;
					while ( temp > 0 ) { std::cout << "*"; temp -= star_every; }
				}
				std::cout << std::endl;
			}
		}

/////////////////////////////////////////////

// decide if we have met the tolerance criteria
		size_t current_count = 0;
		for ( size_t i = 0; i < n; i++ ) {
			size_t bin = idx[ bins.size() - i - 1 ];
			current_count += bins[bin].size();
		}

		float current_frac = (float)current_count / (float)poses.size();
		float error = current_frac - frac;
		std::cout << "Error = " << error << std::endl;
		if ( std::abs(error) < tol ) {
			std::cout << "Tolerance satisfied, clustering complete" << std::endl;
			break;
		}

// find position in array
		std::sort(trial_history.begin(), trial_history.end());
		ptrdiff_t pos = std::find(trial_history.begin(), trial_history.end(), trial ) - trial_history.begin();

// either one plus or one minus
		int other_relevant_pos;
		if ( error > 0 ) {
			other_relevant_pos = pos + 1;
		} else {
			other_relevant_pos = pos - 1;
		}

// naively set the next trial size
		if ( other_relevant_pos < 0 ) {
			trial /= 2;
		} else if ( other_relevant_pos >= trial_history.size() ) {
			trial *= 4;
		} else {
			uint64_t other_trial = trial_history[other_relevant_pos];
			trial = ( other_trial + trial ) / 2;
		}

// cap it to meaningful values
		trial = std::max<uint64_t>( 1, trial );
		trial = std::min<uint64_t>( poses.size(), trial );

// make sure next trial will be unique
		if ( std::find(trial_history.begin(), trial_history.end(), trial ) != trial_history.end() ) {
			std::cout << "Repeating trial, clustering complete" << std::endl;
			break;
		}
	}

	std::vector<core::pose::PoseOP> output_poses;

	for ( size_t i = 0; i < n; i++ ) {
		size_t bin = idx.at( bins.size() - i - 1 );
		std::vector<uint64_t> indexes;
		for ( std::pair<core::pose::PoseOP, uint64_t> pair : bins[bin] ) {
			indexes.push_back( pair.second );
		}
		size_t cluster_center = find_cluster_center( indexes, rmsds );
		output_poses.push_back( bins.at(bin).at(cluster_center).first );
	}

	return output_poses;
}

std::vector<core::pose::PoseOP>
random_selection_poses_leaving_n( 
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n ) {

	if ( n >= poses.size() ) return poses;

	std::vector<uint64_t> indexes( poses.size());
	for ( uint64_t i = 0; i < indexes.size(); i++ ) {
		indexes[i] = i;
	}

	std::random_shuffle(indexes.begin(), indexes.end());

	std::vector<core::pose::PoseOP> output_poses(n);

	for ( uint64_t i = 0; i < n; i++ ) {
		output_poses[i] = poses[indexes[i]];
	}

	return output_poses;
}



#ifdef USEGRIDSCORE
shared_ptr<protocols::ligand_docking::ga_ligand_dock::GridScorer>
prepare_grid_scorer(
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	std::string const & atype_aas /*= "ACDEFGHIKLMNPQRSTVWY" */
) {

	using ObjexxFCL::format::F;

	std::vector<std::pair<core::Real,core::Real>> xyz_min_max(3);
	for ( std::pair<core::Real,core::Real> & pair : xyz_min_max ) { pair.first = 9e9; pair.second = -9e9; }

	for ( core::Size ir : target_res ) {
		for( int ia = 1; ia <= target.residue(ir).nheavyatoms(); ++ia ){
			numeric::xyzVector<core::Real> const & xyz = target.residue(ir).xyz(ia);
			
			for ( int j = 0; j < 3; j++ ) {
				xyz_min_max[j].first = std::min(xyz_min_max[j].first, xyz[j]);
				xyz_min_max[j].second = std::max(xyz_min_max[j].second, xyz[j]);
			}
		}
	}

	numeric::xyzVector< core::Real > centers(3);
	numeric::xyzVector< core::Real > radii(3);

	for ( int j = 0; j < 3; j++ ) {
		centers[j] = ( xyz_min_max[j].first + xyz_min_max[j].second ) / 2.0;
		radii[j] = ( xyz_min_max[j].second - centers[j] );
		runtime_assert( radii[j] > 0 );
	}


	std::cout << "Grid bounds: ";
	for ( int j = 0; j < 3; j++ ) {
		if (j > 0) std::cout << " x ";
		std::cout << "( " << F(5,1,xyz_min_max[j].first) << ", " << F(5,1,xyz_min_max[j].second) << " )";
	}
	std::cout << std::endl;

	radii += 5.0;	// max interaction cutoff

	std::cout << "Padded grid dimensions: ";
	for ( int j = 0; j < 3; j++ ) {
		if (j > 0) std::cout << " x ";
		std::cout << F(5,1,radii[j]); 
	}
	std::cout << std::endl;

	std::cout << "Grid volume: " << F(7,1,radii[0]*radii[1]*radii[2]*8) << std::endl;

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("beta_nov16");
	shared_ptr<protocols::ligand_docking::ga_ligand_dock::GridScorer> grid_scorer 
		= make_shared<protocols::ligand_docking::ga_ligand_dock::GridScorer>( scorefxn );
	grid_scorer->set_grid_com_radii( centers, radii );

	// first A is the "ligand"
	std::string atype_aas_wlig = "A" + atype_aas;
	numeric::xyzVector< core::Real> far_away( 10000, 10000, 10000);
	numeric::xyzMatrix< core::Real> identity = numeric::xyzMatrix< core::Real>::identity();

	std::cout << "Using residue types: " << atype_aas << " for atom types" << std::endl;

	core::pose::PoseOP grid_pose = target.clone();
	utility::vector1< core::Size > fake_movingSCs;
	int ligand_resid = grid_pose->size() + 1;
	for ( int i = 0; i < atype_aas_wlig.size(); i++ ) {
		core::conformation::ResidueOP res = core::conformation::get_residue_from_name1( atype_aas_wlig[i] );
		res->apply_transform_Rx_plus_v( identity, far_away );
		grid_pose->append_residue_by_jump( *res, 1 );
		if ( i > 0 ) fake_movingSCs.push_back( grid_pose->size() );
	}

	grid_scorer->calculate_grid( *grid_pose, ligand_resid, fake_movingSCs );

	return grid_scorer;

}

#endif


// True on success
bool
test_donor_acceptor_caches(
    RifScoreRotamerVsTarget const & rot_tgt_scorer_in,
    shared_ptr<DonorAcceptorCache> const & target_donor_cache,
    shared_ptr<DonorAcceptorCache> const & target_acceptor_cache
) {

	// Prepare the correct and test scoreres
	RifScoreRotamerVsTarget rot_tgt_scorer = rot_tgt_scorer_in;
	rot_tgt_scorer.hbond_weight_ = 1.0;
    RifScoreRotamerVsTarget rot_tgt_scorer_test = rot_tgt_scorer;
	rot_tgt_scorer_test.target_donor_cache_ = target_donor_cache;
	rot_tgt_scorer_test.target_acceptor_cache_ = target_acceptor_cache;

	// get center of hbonds to be more efficient
	Eigen::Vector3f center( 0, 0, 0);
	for ( HBondRay const & ray : rot_tgt_scorer.target_donors_) {
		center += ray.horb_cen;
	}
	for ( HBondRay const & ray : rot_tgt_scorer.target_acceptors_) {
		center += ray.horb_cen;
	}
	center /= rot_tgt_scorer.target_donors_.size() + rot_tgt_scorer.target_acceptors_.size();

	float max_distance = 0;
	for ( HBondRay const & ray : rot_tgt_scorer.target_donors_) {
		max_distance = std::max<float>( max_distance, ( ray.horb_cen - center ).norm() );
	}
	for ( HBondRay const & ray : rot_tgt_scorer.target_acceptors_) {
		max_distance = std::max<float>( max_distance, ( ray.horb_cen - center ).norm() );
	}

	max_distance += 5.0f;


	// Randomly try hbonds until we accumulate 1000 or attempt 1000000
    std::vector<std::mt19937> rng_pt; 
	for( int i  = 0; i < ::devel::scheme::omp_max_threads_1(); ++i ){
		rng_pt.emplace_back(time(0) + i);
	}

	size_t max_iter = 1000000;
	size_t required = 1000;
	bool all_pass = true;
	size_t so_far = 0;

	#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
	for ( size_t i = 0; i < max_iter; i++ ) {

		if ( ! all_pass || so_far >= required ) continue;

		std::mt19937 & rng = rng_pt.at(::devel::scheme::omp_thread_num());
		std::uniform_real_distribution<> runif;

		Eigen::Vector3f position;
		::scheme::numeric::rand_vector_sphere(rng, position);
		position.normalize();
		position *= max_distance * runif(rng);
		position += center;

		Eigen::Vector3f direction;
		::scheme::numeric::rand_vector_sphere(rng, direction);
		direction.normalize();

		HBondRay ray;
		ray.horb_cen[0] = position[0]; ray.horb_cen[1] = position[1]; ray.horb_cen[2] = position[2];
		ray.direction[0] = direction[0]; ray.direction[1] = direction[1]; ray.direction[2] = direction[2];

		std::vector<HBondRay> rays { ray };

		// Test Acceptors
		{
			int sat1 = -1, sat2 = -1, hbcount = 0;
			float score = rot_tgt_scorer.score_acceptor_rays_v_target( rays, sat1, sat2, hbcount );

			int og_sat1 = -1, og_sat2 = -1, og_hbcount = 0;
			float og_score = rot_tgt_scorer_test.score_acceptor_rays_v_target( rays, og_sat1, og_sat2, og_hbcount );

			if ( score < -0.05) {
				#ifdef USE_OPENMP
				#pragma omp critical
				#endif
				{
					so_far += 1;
				}
				if ( sat1 != og_sat1 || sat2 != og_sat2 || score != og_score ) {
					// std::cout << "acceptors: " << score << " " << sat1 << " " << sat2 << " " << og_score << " " << og_sat1 << " " << og_sat2 << std::endl;
					all_pass = false;
				}
			}

		}

	}
	if ( ! all_pass ) {
		std::cout << "Hbond error" << std::endl;
		return false;
	}
	std::cout << "Tested " << so_far << " hbonds without errors" << std::endl;

	return true;

}

// True on success
bool
inner_prepare_donor_acceptor_cache( 
    std::vector<HBondRay> const & target_donors,
    std::vector<HBondRay> const & target_acceptors,
    RifScoreRotamerVsTarget const & rot_tgt_scorer,
    float safe_distance,
    shared_ptr<DonorAcceptorCache> & target_donor_cache,
    shared_ptr<DonorAcceptorCache> & target_acceptor_cache
) {
							 // ideal length - orbital position + max extension + extra extension + resl + safety
	const float max_hbond_interaction = ( 2.00 - 0.61 + 0.8 + rot_tgt_scorer.long_hbond_fudge_distance_ + 1.0 + safe_distance );
	std::cout << "Donors: ";
	target_donor_cache = make_shared<DonorAcceptorCache>( target_donors, max_hbond_interaction );
	std::cout << "Acceptors: ";
	target_acceptor_cache = make_shared<DonorAcceptorCache>( target_acceptors, max_hbond_interaction );

	return test_donor_acceptor_caches( rot_tgt_scorer, target_donor_cache, target_acceptor_cache );
}

void
prepare_donor_acceptor_cache( 
    std::vector<HBondRay> const & target_donors,
    std::vector<HBondRay> const & target_acceptors,
    RifScoreRotamerVsTarget const & rot_tgt_scorer,
    shared_ptr<DonorAcceptorCache> & target_donor_cache,
    shared_ptr<DonorAcceptorCache> & target_acceptor_cache
) {
	if (target_donors.size() == 0 && target_acceptors.size() == 0) return;

	float safe = 0.1;
	while ( safe < 5.0f ) {
		std::cout << boost::str(boost::format("Building DonorAcceptorCache with %.1f safe distance")%safe) << std::endl;
		if ( inner_prepare_donor_acceptor_cache( target_donors, target_acceptors, rot_tgt_scorer, safe, target_donor_cache, target_acceptor_cache ) ) {
			return;
		}
		safe += 0.5f;
	}

	utility_exit_with_message("Something is horribly wrong. Tell Brian. (bcov@uw.edu)");

}



}
}


