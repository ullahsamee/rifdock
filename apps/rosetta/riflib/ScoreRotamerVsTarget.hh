// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_ScoreRotamerVsTarget_hh
#define INCLUDED_riflib_ScoreRotamerVsTarget_hh

#include <riflib/util.hh>
#include <scheme/numeric/util.hh>
#include <scheme/chemical/HBondRay.hh>
#include <riflib/DonorAcceptorCache.hh>

#ifdef USEGRIDSCORE
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <protocols/ligand_docking/GALigandDock/RotamerData.hh>
#endif


namespace devel {
namespace scheme {


inline
float
dot( Eigen::Vector3f const & v1, Eigen::Vector3f const & v2 ) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

// on the donor:
// horb_cen = H xyz
// horb_cen - direction = O/N xyz

// on the acceptor:
// horb_cen = orbital by rosetta definition
// horb_cen - direction * ORBLEN = N/O xyz

// ORBLEN = 0.61

// Notes from brian and longxing
// This used to be based on the H -> orbital distance
// This gave a pearsonr=0.1
// After changing it to use the H -> acceptor distance
//   the pearsonr=0.64
//
// Further improvements:
//   Dirscore used to be based on the acceptor donor ray dot product
//   A better way is to dot both against the A -> H ray
//   Doing this and either averaging or multiplying the dots gives r=0.88
//   Squaring the H term gives r=0.93

// This code is *INCREDIBLY* hot. Be very careful here
template< class HBondRay >
float score_hbond_rays(
    HBondRay const & don,
    HBondRay const & acc,
    float non_directional_fraction = 0.0, // 0.0 - 1.0, // this is disabled by a comment below
    float long_hbond_fudge_distance = 0.0
){
    using ::scheme::numeric::sqr;

    const Eigen::Vector3f accep_O = acc.horb_cen - acc.direction*::scheme::chemical::ORBLEN;


    const Eigen::Vector3f h_to_a = ( accep_O - don.horb_cen ).normalized();
    // float const h_dirscore = don.direction.dot( h_to_a );
    float const h_dirscore = dot( don.direction, h_to_a );
    if ( h_dirscore <= 0 ) return 0;

    // float const a_dirscore = acc.direction.dot( h_to_a ) * -1.0;
    float const a_dirscore = dot( acc.direction, h_to_a ) * -1.0;
    if ( a_dirscore <= 0 ) return 0;


    float diff = ( (don.horb_cen-accep_O).norm() - 2.00 );
    diff = diff < 0 ? diff*1.5 : std::max<float>( 0.0, diff - long_hbond_fudge_distance ); // increase dis pen if too close
    float const max_diff = 0.8;
    if ( diff >= max_diff || diff <= -max_diff ) return 0;

    // sigmoid -like shape on distance score
    float score = sqr( 1.0 - sqr( diff/max_diff ) ) * -1.0;
    assert( score <= 0.0 );

    float dirscore = h_dirscore * h_dirscore * a_dirscore;

// Uncomment this if you want to use non_directional_fraction
    // float const nds = non_directional_fraction;
    // score = ( 1.0 - nds )*( score * dirscore ) + nds * score;

    return score * dirscore;
}

template< class VoxelArrayPtr, class HBondRay, class RotamerIndex >
struct ScoreRotamerVsTarget {
    ::scheme::shared_ptr< RotamerIndex const > rot_index_p_ = nullptr;
    std::vector<VoxelArrayPtr> target_field_by_atype_;
    std::vector< HBondRay > target_donors_, target_acceptors_;

    // Ok, the names shouldn't be here, but it's convenient
    std::vector<std::pair<int, std::string> > target_donor_names;
    std::vector<std::pair<int, std::string> > target_acceptor_names;

    float hbond_weight_ = 2.0;
    float upweight_iface_ = 1.0;
    float upweight_multi_hbond_ = 0.0;
    float min_hb_quality_for_multi_ = -0.5;
    float min_hb_quality_for_satisfaction_ = -0.6;
    float long_hbond_fudge_distance_ = 0.0;
#ifdef USEGRIDSCORE
    shared_ptr<protocols::ligand_docking::ga_ligand_dock::GridScorer> grid_scorer_;
    bool soft_grid_energies_;
#endif
    shared_ptr<DonorAcceptorCache> target_donor_cache_;
    shared_ptr<DonorAcceptorCache> target_acceptor_cache_;

    ScoreRotamerVsTarget(){}

    template< class Xform, class Int >
    float
    score_rotamer_v_target(
        Int const & irot,
        Xform const & rbpos,
        float bad_score_thresh = 10.0, // hbonds won't get computed if grid score is above this
        int start_atom = 0 // to score only SC, use 4... N,CA,C,CB (?)
    ) const {
        int tmp1=-12345, tmp2=-12345, hbcount=0;
        return score_rotamer_v_target_sat( irot, rbpos, tmp1, tmp2, false, hbcount, bad_score_thresh, start_atom );
    }

    template< class Xform, class Int >
    float
    score_rotamer_v_target_sat(
        Int const & irot,
        Xform const & rbpos,
        int & sat1,
        int & sat2,
        bool want_sats, // need to do double score if we want_sats and grid scoring
        int & hbcount, // how many hbonds? Requires (want_sat or ! grid_scorer_)
        float bad_score_thresh = 10.0, // hbonds won't get computed if grid score is above this
        int start_atom = 0 // to score only SC, use 4... N,CA,C,CB (?)
    ) const {
        using devel::scheme::score_hbond_rays;
        assert( rot_index_p_ );
        assert( target_field_by_atype_.size() == 22 );
        float score = 0;
        typedef typename RotamerIndex::Atom Atom;

        bool use_grid_scorer = false;
#ifdef USEGRIDSCORE
        use_grid_scorer = (bool)grid_scorer_;
#endif

        if ( use_grid_scorer ) {
#ifdef USEGRIDSCORE
            core::conformation::ResidueOP residue = rot_index_p_->get_per_thread_rotamer_at_identity(omp_get_thread_num(), irot);
            apply_xform_to_residue( *residue, rbpos );
            core::scoring::lkball::LKB_ResidueInfoOP lkbrinfo = rot_index_p_->get_per_thread_lkbrinfo(omp_get_thread_num(), irot);
            protocols::ligand_docking::ga_ligand_dock::ReweightableRepEnergy rerep_energy 
                = grid_scorer_->get_1b_energy( *residue, lkbrinfo, soft_grid_energies_, true );
            score += rerep_energy.score(1.0);
#endif
        } else {
            for( int iatom = start_atom; iatom < rot_index_p_->nheavyatoms(irot); ++iatom )
            {
                Atom const & atom = rot_index_p_->rotamer(irot).atoms_.at(iatom);
                typename Atom::Position pos = rbpos * atom.position();
                score += target_field_by_atype_.at(atom.type())->at( pos );
            }
        }


        bool calculate_hbonds = ( score < bad_score_thresh ) && ( ! use_grid_scorer || want_sats );
        // in one test: 244m with this, 182m without... need to optimize...
        if( calculate_hbonds ){
            float hbscore = 0;
            // int hbcount = 0;
            if( rot_index_p_->rotamer(irot).acceptors_.size() > 0 ||
                rot_index_p_->rotamer(irot).donors_   .size() > 0 )
            {

                std::vector<HBondRay> acceptor_rays = rot_index_p_->rotamer(irot).acceptors_;
                for ( HBondRay & hb_ray : acceptor_rays ) hb_ray.apply_xform( rbpos );
                hbscore += score_acceptor_rays_v_target( acceptor_rays, sat1, sat2, hbcount );

                std::vector<HBondRay> donor_rays = rot_index_p_->rotamer(irot).donors_;
                for ( HBondRay & hb_ray : donor_rays ) hb_ray.apply_xform( rbpos );
                hbscore += score_donor_rays_v_target( donor_rays, sat1, sat2, hbcount );
                
            }

            // oh god, fix me..... what should the logic be??? probably "softer" thresh on thishb to count
            if( upweight_multi_hbond_ != 0 ){
                if( rot_index_p_->resname(irot)!="TYR" && // hack to aleviate my OH problems...
                    rot_index_p_->resname(irot)!="SER" &&
                    rot_index_p_->resname(irot)!="THR"
                ){
                    float multihb = 0;
                    int nchi = rot_index_p_->nchi(irot) - rot_index_p_->nprotonchi(irot);
                    switch( nchi ){
                        case 0:
                        case 1:
                            if( hbcount <= 1 ) hbscore *= 1.0;
                            if( hbcount > 1 ) multihb = 0.7*(hbcount-1);
                            break;
                        case 2:
                            if( hbcount <= 1 ) hbscore *= 0.9;
                            if( hbcount > 1 ) multihb = 1.0*(hbcount-1);
                            break;
                        case 3:
                            if( hbcount <= 1 ) hbscore *= 0.7;
                            if( hbcount > 1 ) multihb = 0.8*(hbcount-1);
                            break;
                        default:
                            if( hbcount <= 1 ) hbscore *= 0.5;
                            if( hbcount > 1 ) multihb = 0.7*(hbcount-1);
                            break;
                    }
                    multihb = std::max( 0.0f, multihb );
                    hbscore += hbscore * multihb * upweight_multi_hbond_; // should multihb be additive or multiplicative?
                }
            }
            if( hbscore < 0 && ! use_grid_scorer ) score += hbscore;
        }
        return score * upweight_iface_;
    }

    float
    score_acceptor_rays_v_target( std::vector<HBondRay> const & acceptor_rays, int & sat1, int & sat2, int & hbcount ) const {
        float hbscore = 0;
        const size_t target_donors_size = target_donors_.size();

        // this is faster than std::vector
        float used_tgt_donor   [target_donors_size];
        for( int i = 0; i < target_donors_size; ++i ) used_tgt_donor   [i] = 9e9;


        for( HBondRay const & hr_rot_acc : acceptor_rays ) {
            
            float best_score = 100;
            int best_sat = -1;

            if ( target_donor_cache_ ) {

                typedef typename DonorAcceptorCache::Sat Sat;
                std::vector<Sat>::const_iterator sats_iter = target_donor_cache_->at( hr_rot_acc.horb_cen );

                Sat i_hr_tgt_don = 0;
                while ( (i_hr_tgt_don = *(sats_iter++)) != DonorAcceptorCache::CACHE_MAX_SAT ) {

                    /////////////// DUPLICATE CODE ///////////////////////////////////////////
                    HBondRay const & hr_tgt_don = target_donors_.at(i_hr_tgt_don);
                    float const thishb = score_hbond_rays( hr_tgt_don, hr_rot_acc, 0.0, long_hbond_fudge_distance_ );
                    if ( thishb < best_score ) {
                        best_score = thishb;
                        best_sat = i_hr_tgt_don;
                    }
                    /////////////// DUPLICATE CODE ///////////////////////////////////////////

                }

            } else {
                for( int i_hr_tgt_don = 0; i_hr_tgt_don < target_donors_.size(); ++i_hr_tgt_don )
                {
                    /////////////// DUPLICATE CODE ///////////////////////////////////////////
                    HBondRay const & hr_tgt_don = target_donors_.at(i_hr_tgt_don);
                    float const thishb = score_hbond_rays( hr_tgt_don, hr_rot_acc, 0.0, long_hbond_fudge_distance_ );
                    if ( thishb < best_score ) {
                        best_score = thishb;
                        best_sat = i_hr_tgt_don;
                    }
                    /////////////// DUPLICATE CODE ///////////////////////////////////////////
                }
            }

            if ( best_sat > -1 && used_tgt_donor[best_sat] > best_score ) {
                // I don't think there is any need to use if ... else ..., but to make things more clear.
                if ( used_tgt_donor[ best_sat ] < this->min_hb_quality_for_satisfaction_ ){
                    if ( sat1 == -1 || sat1 == best_sat ){
                        sat1 = best_sat;
                    } else if ( sat2 == -1 || sat2 == best_sat) {
                        sat2 = best_sat;
                    }
                } else if ( best_score < this->min_hb_quality_for_satisfaction_ ) {
                    if(      sat1==-1 ) sat1 = best_sat;
                    else if( sat2==-1 ) sat2 = best_sat;
                }
                if( best_score < min_hb_quality_for_multi_ ){
                    if( used_tgt_donor[best_sat] >= min_hb_quality_for_multi_ ) ++hbcount;
                }
                // remove the double counting of hbond??? Do I need to do this, or ..........
                if ( used_tgt_donor[ best_sat ] <= 9e5 ){
                    hbscore += ( best_score - used_tgt_donor[best_sat] ) * hbond_weight_;
                } else {
                    hbscore += best_score * hbond_weight_;
                }
                used_tgt_donor[ best_sat ] = best_score;
            }
        }

        return hbscore;

    }



    float
    score_donor_rays_v_target( std::vector<HBondRay> const & donor_rays, int & sat1, int & sat2, int & hbcount ) const {
        float hbscore = 0;

        const size_t target_acceptors_size = target_acceptors_.size();

        // This is faster than std::vector
        float used_tgt_acceptor[target_acceptors_size];
        for( int i = 0; i < target_acceptors_size; ++i ) used_tgt_acceptor[i] = 9e9;

        for( HBondRay const & hr_rot_don : donor_rays ) {
            
            float best_score = 100;
            int best_sat = -1;


            if ( target_acceptor_cache_ ) {

                typedef typename DonorAcceptorCache::Sat Sat;
                std::vector<Sat>::const_iterator sats_iter =  target_acceptor_cache_->at( hr_rot_don.horb_cen );

                Sat i_hr_tgt_acc = 0;
                while ( (i_hr_tgt_acc = *(sats_iter++)) != DonorAcceptorCache::CACHE_MAX_SAT ) {

                    /////////////// DUPLICATE CODE ///////////////////////////////////////////
                    HBondRay const & hr_tgt_acc = target_acceptors_.at(i_hr_tgt_acc);
                    float const thishb = score_hbond_rays( hr_rot_don, hr_tgt_acc, 0.0, long_hbond_fudge_distance_ );
                    if ( thishb < best_score ) {
                        best_score = thishb;
                        best_sat = i_hr_tgt_acc + target_donors_.size();
                    }
                    /////////////// DUPLICATE CODE ///////////////////////////////////////////
                }

            } else {
                for( int i_hr_tgt_acc = 0; i_hr_tgt_acc < target_acceptors_.size(); ++i_hr_tgt_acc )
                {
                    /////////////// DUPLICATE CODE ///////////////////////////////////////////
                    HBondRay const & hr_tgt_acc = target_acceptors_.at(i_hr_tgt_acc);
                    float const thishb = score_hbond_rays( hr_rot_don, hr_tgt_acc, 0.0, long_hbond_fudge_distance_ );
                    if ( thishb < best_score ) {
                        best_score = thishb;
                        best_sat = i_hr_tgt_acc + target_donors_.size();
                    }
                    /////////////// DUPLICATE CODE ///////////////////////////////////////////
                }
            }


            if ( best_sat > -1 && used_tgt_acceptor[best_sat - target_donors_.size()] > best_score ) {
                // I don't think there is any need to use if ... else ..., but to make things more clear.
                if ( used_tgt_acceptor[ best_sat - target_donors_.size() ] < this->min_hb_quality_for_satisfaction_ ){
                    if ( sat1 == -1 || sat1 == best_sat ){
                        sat1 = best_sat;
                    } else if ( sat2 == -1 || sat2 == best_sat) {
                        sat2 = best_sat;
                    }
                } else if ( best_score < this->min_hb_quality_for_satisfaction_ ) {
                    if(      sat1==-1 ) sat1 = best_sat;
                    else if( sat2==-1 ) sat2 = best_sat;
                }
                if( best_score < min_hb_quality_for_multi_ ){
                    if( used_tgt_acceptor[best_sat - target_donors_.size()] >= min_hb_quality_for_multi_ ) ++hbcount;
                }
                // remove the double counting of hbond??? Do I need to do this, or ..........
                if ( used_tgt_acceptor[ best_sat - target_donors_.size() ] <= 9e5 ){
                    hbscore += ( best_score - used_tgt_acceptor[best_sat - target_donors_.size()] ) * hbond_weight_;
                } else {
                    hbscore += best_score * hbond_weight_;
                }
                used_tgt_acceptor[ best_sat - target_donors_.size() ] = best_score;
            }
            
        }

        return hbscore;

    }



};







}
}

#endif
