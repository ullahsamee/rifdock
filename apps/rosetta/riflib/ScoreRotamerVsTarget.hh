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


#ifdef USEGRIDSCORE
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <protocols/ligand_docking/GALigandDock/RotamerData.hh>
#endif


namespace devel {
namespace scheme {

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
// The plot of rif energies vs rosetta eneriges now looks like a triangle
//   Everything rif think is bad, rosetta agrees
//   But rif accepts some things that rosetta does not    

template< class HBondRay >
float score_hbond_rays(
    HBondRay const & don,
    HBondRay const & acc,
    float non_directional_fraction = 0.0 // 0.0 - 1.0
){

    const Eigen::Vector3f accep_O = acc.horb_cen - acc.direction*ORBLEN;
    float diff = ( (don.horb_cen-accep_O).norm() - 2.00 );
    diff = diff < 0 ? diff*1.5 : diff; // increase dis pen if too close
    float const max_diff = 0.8;
    diff = diff >  max_diff ?  max_diff : diff;
    diff = diff < -max_diff ? -max_diff : diff;
    // if( diff > max_diff ) return 0.0;
    // sigmoid -like shape on distance score
    float score = sqr( 1.0 - sqr( diff/max_diff ) ) * -1.0;
    assert( score <= 0.0 );
    float dirscore = -don.direction.dot( acc.direction );
    dirscore = dirscore < 0 ? 0.0 : dirscore; // is positive
    float const nds = non_directional_fraction;
    score = ( 1.0 - nds )*( score * dirscore ) + nds * score;

    return score;
}

template< class VoxelArrayPtr, class HBondRay, class RotamerIndex >
struct ScoreRotamerVsTarget {
    ::scheme::shared_ptr< RotamerIndex const > rot_index_p_ = nullptr;
    std::vector<VoxelArrayPtr> target_field_by_atype_;
    std::vector< HBondRay > target_donors_, target_acceptors_;
    float hbond_weight_ = 2.0;
    float upweight_iface_ = 1.0;
    float upweight_multi_hbond_ = 0.0;
    float min_hb_quality_for_multi_ = -0.5;
    float min_hb_quality_for_satisfaction_ = -0.6;
#ifdef USEGRIDSCORE
    shared_ptr<protocols::ligand_docking::ga_ligand_dock::GridScorer> grid_scorer_;
    bool soft_grid_energies_;
#endif
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
        int & hbcount, // how many hbonds? Requires (want_sat or ! grid_scorer_) and upweight_multi_hbond_
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
                // alloca style stack bump... dangerous... don't piss memory here...
                bool used_tgt_donor   [target_donors_   .size()];
                bool used_tgt_acceptor[target_acceptors_.size()];
                bool used_rot_donor   [rot_index_p_->rotamer(irot).donors_   .size()];
                bool used_rot_acceptor[rot_index_p_->rotamer(irot).acceptors_.size()];
                for( int i = 0; i < target_donors_   .size(); ++i ) used_tgt_donor   [i] = false;
                for( int i = 0; i < target_acceptors_.size(); ++i ) used_tgt_acceptor[i] = false;
                for( int i = 0; i < rot_index_p_->rotamer(irot).donors_   .size(); ++i ) used_rot_donor   [i] = false;
                for( int i = 0; i < rot_index_p_->rotamer(irot).acceptors_.size(); ++i ) used_rot_acceptor[i] = false;

                for( int i_hr_rot_acc = 0; i_hr_rot_acc < rot_index_p_->rotamer(irot).acceptors_.size(); ++i_hr_rot_acc )
                {
                    HBondRay hr_rot_acc = rot_index_p_->rotamer(irot).acceptors_.at(i_hr_rot_acc);
                    Eigen::Vector3f dirpos = hr_rot_acc.horb_cen + hr_rot_acc.direction;
                    hr_rot_acc.horb_cen  = rbpos * hr_rot_acc.horb_cen;
                    hr_rot_acc.direction = rbpos * dirpos - hr_rot_acc.horb_cen;
                    for( int i_hr_tgt_don = 0; i_hr_tgt_don < target_donors_.size(); ++i_hr_tgt_don )
                    {
                        HBondRay const & hr_tgt_don = target_donors_.at(i_hr_tgt_don);
                        float const thishb = score_hbond_rays( hr_tgt_don, hr_rot_acc );
                        hbscore += thishb * hbond_weight_;
                        if( thishb < this->min_hb_quality_for_satisfaction_ ){
                            if(      sat1==-1 ) sat1 = i_hr_tgt_don;
                            else if( sat2==-1 ) sat2 = i_hr_tgt_don;
                        }
                        if( upweight_multi_hbond_ && thishb < min_hb_quality_for_multi_ ){
                            if( !used_tgt_donor[i_hr_tgt_don] && !used_rot_acceptor[i_hr_rot_acc] ) ++hbcount;
                            used_tgt_donor   [i_hr_tgt_don] = true;
                            used_rot_acceptor[i_hr_rot_acc] = true;
                        }
                    }
                }
                for( int i_hr_rot_don = 0; i_hr_rot_don < rot_index_p_->rotamer(irot).donors_.size(); ++i_hr_rot_don )
                {
                    HBondRay hr_rot_don = rot_index_p_->rotamer(irot).donors_.at(i_hr_rot_don);
                    Eigen::Vector3f dirpos = hr_rot_don.horb_cen + hr_rot_don.direction;
                    hr_rot_don.horb_cen  = rbpos * hr_rot_don.horb_cen;
                    hr_rot_don.direction = rbpos * dirpos - hr_rot_don.horb_cen;
                    for( int i_hr_tgt_acc = 0; i_hr_tgt_acc < target_acceptors_.size(); ++i_hr_tgt_acc )
                    {
                        HBondRay const & hr_tgt_acc = target_acceptors_.at(i_hr_tgt_acc);
                        float const thishb = score_hbond_rays( hr_rot_don, hr_tgt_acc );
                        hbscore += thishb * hbond_weight_;
                        if( thishb < this->min_hb_quality_for_satisfaction_ ){
                            if(      sat1==-1 ) sat1 = i_hr_tgt_acc + target_donors_.size();
                            else if( sat2==-1 ) sat2 = i_hr_tgt_acc + target_donors_.size();
                        }
                        if( upweight_multi_hbond_ && thishb < min_hb_quality_for_multi_ ){
                            if( !used_rot_donor[i_hr_rot_don] && !used_tgt_acceptor[i_hr_tgt_acc] ) ++hbcount;
                            used_tgt_acceptor[i_hr_tgt_acc] = true;
                            used_rot_donor[i_hr_rot_don] = true;
                        }
                    }
                }
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


};







}
}

#endif
