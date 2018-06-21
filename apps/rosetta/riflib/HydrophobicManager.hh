// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_HydrophobicManager_hh
#define INCLUDED_riflib_HydrophobicManager_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>

#include "scheme/util/SimpleArray.hh"
#include <riflib/RotamerGenerator.hh>

#include <core/pose/Pose.hh>

#include <string>
#include <vector>

#include <boost/any.hpp>
#include <boost/format.hpp>

#include <riflib/util.hh>


namespace devel {
namespace scheme {


// This class is pretty sketch
// We're trying to figure out whether the scaffold interacts with each target hydrophobic residue with more than
//   0.5 kcal fa_sol, fa_atr, fa_rep
//
// The plan:
// Make a 0.5 A grid
// Give each hydrophobic a Hyd number
// Start filling in the grid for each hydrophobic carbon
// From 3.0 - 4.5: fill the grid twice with the Hyd // not doing this for now
// From 4.5 - 6.0: fill the grid once with the Hyd

// Then after packing:
// Scan through all rifres carbon atoms and start tallying the Hyds
// Any Hyd that gets more than 2 is tallied
// This count is the hydrophobic_residue_contacts



struct HydrophobicManager {

    typedef uint16_t Hyd;
    static int const CACHE_MAX_HYD = std::numeric_limits<Hyd>::max();

    typedef ::scheme::util::SimpleArray<3,size_t> Indices;
    typedef Eigen::Vector3f Bounds;
    Bounds lb_,ub_,cs_;

    const float max_interaction_range_ = 7.0;   // It's actually 6.0, but just to be safe
    size_t max_hyds_;

    const std::set<char> hydrophobic_name1s_ {'A', 'C', 'F', 'I', 'L', 'M', 'P', 'T', 'V', 'W', 'Y'}; 
    std::vector<core::Size> hydrophobic_res_;

    // Cache these for speed
    std::vector<int> rif_atype_map_;
    std::vector<bool> rif_hydrophobic_map_;


    std::vector<Hyd> voxel_map_;

    // This only actually needs 1. But include a few extra in case it accidentally gets written to
    std::vector<Hyd> OOB_LIST = {CACHE_MAX_HYD, CACHE_MAX_HYD, CACHE_MAX_HYD, CACHE_MAX_HYD, CACHE_MAX_HYD, CACHE_MAX_HYD};

    Indices shape_;

    shared_ptr< RotamerIndex > rot_index_p;

    HydrophobicManager(
        core::pose::Pose const & target,
        utility::vector1<core::Size> const & target_res,
        shared_ptr< RotamerIndex > rot_index_p_in
    ) : rot_index_p( rot_index_p_in )
    {

        rif_atype_map_ = get_rif_atype_map();
        rif_hydrophobic_map_ = get_rif_hydrophobic_map();

        identify_hydrophobic_residues( target, target_res );

        prepare_bounds( target );
        std::vector<std::vector<Hyd>> early_map = first_pass_fill( target );

        create_and_fill_voxel_map( early_map );
    }

    void
    identify_hydrophobic_residues( 
        core::pose::Pose const & target,
        utility::vector1<core::Size> const & target_res
    ) {
        hydrophobic_res_.clear();

        for ( core::Size seqpos : target_res ) {
            if ( hydrophobic_name1s_.count( target.residue(seqpos).name1() ) != 0 ) {
                hydrophobic_res_.push_back(seqpos);
            }
        }
    }

    bool
    is_rosetta_atom_hydrophobic( core::conformation::Residue const & res, core::Size atno ) {
        int iatype = rif_atype_map_.at(res.atom_type_index(atno));
        if (iatype >= rif_hydrophobic_map_.size()) {
            return false;
        }
        return rif_hydrophobic_map_.at(iatype);
    }


    void
    prepare_bounds( core::pose::Pose const & target ) {

        Eigen::Vector3f lbs( 9e9, 9e9, 9e9 );
        Eigen::Vector3f ubs( -9e9, -9e9, -9e9 );

        for ( core::Size seqpos : hydrophobic_res_ ) {
            core::conformation::Residue const & res = target.residue(seqpos);
            for ( core::Size atno = 4; atno < res.nheavyatoms(); atno++ ) {  // Don't include backbone
                if ( ! is_rosetta_atom_hydrophobic( res, atno ) ) continue;

                numeric::xyzVector<core::Real> xyz = res.xyz(atno);
                for ( int i = 0; i < 3; i++ ) {
                    lbs[i] = std::min<float>( lbs[i], xyz[i] );
                    ubs[i] = std::max<float>( ubs[i], xyz[i] );
                }
            }
            
        }

        for ( int i = 0; i < 3; i++ ) {
            lbs[i] -= max_interaction_range_;
            ubs[i] += max_interaction_range_;
        }

        lb_ = lbs;
        ub_ = ubs;
        cs_ = Eigen::Vector3f( 0.5, 0.5, 0.5 );

        Indices extents = floats_to_index( ub_ );
        shape_ = extents + Indices(1);

    }

    std::vector<std::vector<Hyd>>
    first_pass_fill( core::pose::Pose const & target ) {

        std::vector<std::vector<Hyd>> early_map;

        size_t elements = shape_[0] * shape_[1] * shape_[2];
        // std::cout << shape_ << std::endl;
        // std::cout << floats_to_index( ub_ ) << std::endl;
        // std::cout << elements << " " << index_to_map_index( floats_to_index( ub_ ) ) << std::endl;
        runtime_assert( elements - 1 == index_to_map_index( floats_to_index( ub_ ) ) );

        early_map.resize(elements);

        for ( Hyd ihyd = 0; ihyd < hydrophobic_res_.size(); ihyd++ ) {
            core::conformation::Residue const & res = target.residue(hydrophobic_res_[ihyd]);

            for ( core::Size atno = 4; atno < res.nheavyatoms(); atno++ ) {  // Don't include backbone
                if ( ! is_rosetta_atom_hydrophobic( res, atno ) ) continue;

                numeric::xyzVector<core::Real> _xyz = res.xyz( atno );
                Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

                Eigen::Vector3f lbs( xyz[0] - max_interaction_range_, xyz[1] - max_interaction_range_, xyz[2] - max_interaction_range_ );
                Eigen::Vector3f ubs( xyz[0] + max_interaction_range_, xyz[1] + max_interaction_range_, xyz[2] + max_interaction_range_ );

                const float low_rad_sq = 3.0f*3.0f;
                const float med_rad_sq = 4.5f*4.5f;
                const float long_rad_sq = 6.0f*6.0f;

                const float step = cs_[0];

                Eigen::Vector3f worker;

                for ( float x = lbs[0] - step/2; x < ubs[0] + step; x += step ) {
                    if ( x < lb_[0] || x > ub_[0] ) continue;
                    worker[0] = x;

                    for ( float y = lbs[1] - step/2; y < ubs[1] + step; y += step ) {
                        if ( y < lb_[1] || y > ub_[1] ) continue;
                        worker[1] = y;

                        for ( float z = lbs[2] - step/2; z < ubs[2] + step; z += step ) {
                            if ( z < lb_[2] || z > ub_[2] ) continue;
                            worker[2] = z;

                            const float squared_dist = ( xyz - worker ).squaredNorm();

                            if ( squared_dist < low_rad_sq ) {
                                // do nothing
                            } else if ( squared_dist < med_rad_sq ) {
                                size_t offset = index_to_map_index( floats_to_index( worker ) );
                                early_map.at(offset).push_back( ihyd );
                                // early_map.at(offset).push_back( ihyd );
                            } else if ( squared_dist < long_rad_sq ) {
                                size_t offset = index_to_map_index( floats_to_index( worker ) );
                                early_map.at(offset).push_back( ihyd );
                            } else {
                                // do nothing
                            }
                        }
                    }
                }
            }
        }
        return early_map;
    }



    void
    create_and_fill_voxel_map( std::vector<std::vector<Hyd>> const & early_map ) {

        max_hyds_ = 0;

        for ( std::vector<Hyd> these_hyds : early_map ) {
            max_hyds_ = std::max<size_t>( max_hyds_, these_hyds.size() );
        }

        std::cout << "Max hydrophobics at one voxel: " << max_hyds_ << std::endl;

        // max_hyds_ is 1 too big. This means that the receiver can always loop until they hit a CACHE_MAX_HYD
        max_hyds_ += 1;

        voxel_map_.resize( max_hyds_ * shape_[0] * shape_[1] * shape_[2], CACHE_MAX_HYD );

        // We already asserted the two maps have the same size, so just array fill
        for ( size_t imap = 0; imap < early_map.size(); imap++ ) {

            std::vector<Hyd> these_hyds = early_map[imap];
            bool wrote_a_null = false;
            for ( size_t idx = 0; idx < max_hyds_; idx++ ) {
                size_t offset = imap * max_hyds_ + idx;
                if ( idx < these_hyds.size() ) {
                    voxel_map_.at(offset) = these_hyds.at(idx);
                } else {
                    voxel_map_.at(offset) = CACHE_MAX_HYD;
                    wrote_a_null = true;
                }
            }
            runtime_assert( wrote_a_null );     // Redundancy. Make sure each list is terminated
        }
    }


    int
    find_hydrophobic_residue_contacts( 
        std::vector<std::pair<intRot, EigenXform>> const & irot_and_bbpos,
        std::vector<int> & hyd_counts,
        float & hydrophobic_ddg
    ) const {
        hyd_counts.clear();
        hyd_counts.resize( hydrophobic_res_.size(), 0 );

        typedef typename RotamerIndex::Atom Atom;

        for ( std::pair<intRot, EigenXform> pair : irot_and_bbpos ) {
            int irot = pair.first;
            EigenXform const & bbpos = pair.second;

            for( int iatom = 4; iatom < rot_index_p->nheavyatoms(irot); ++iatom )
            {
                Atom const & atom = rot_index_p->rotamer(irot).atoms_.at(iatom);

                if ( ! rif_hydrophobic_map_.at(atom.type()) ) continue;

                typename Atom::Position pos = bbpos * atom.position();

                std::vector<Hyd>::const_iterator hyds_iter = this->at( pos );

                Hyd this_hyd = 0;
                while ( (this_hyd = *(hyds_iter++)) != CACHE_MAX_HYD ) {
                    hyd_counts[this_hyd] ++;
                }
            }
        }

        int total_sum = 0;
        int hydrophobic_residue_contacts = 0;
        for ( int count : hyd_counts ) {
            if (count >= 2) {
                hydrophobic_residue_contacts++;
            }
            total_sum += count;
        }

        hydrophobic_ddg = -0.17f * total_sum; // This is how it was parameterized, but holy crap this is a hack.
        return hydrophobic_residue_contacts;
    }

    void
    print_hydrophobic_counts( core::pose::Pose const & target, std::vector<int> const & hyd_counts, int scaff_size ) const {
        std::cout << "Hydrophobic counts:" << std::endl;

        assert( hyd_counts.size() == hydrophobic_res_.size() );

        for ( int i = 0; i < hyd_counts.size(); i++ ) {
            int count = hyd_counts[i];
            if ( count >= 2 ) {
                core::Size seqpos = hydrophobic_res_[i];
                core::conformation::Residue const & res = target.residue(seqpos);
                std::cout << " " << seqpos << "/" << seqpos + scaff_size << " " << res.name3() << ": " << count << std::endl;
            }
        }
    }




    template<class Floats> Indices floats_to_index(Floats const & f) const {
        Indices ind;
        for(int i = 0; i < 3; ++i){
            float tmp = ((f[i]-lb_[i])/cs_[i]);
            ind[i] = tmp;
        }
        return ind;
    }

    size_t index_to_map_index( Indices const & ind ) const {

        size_t accum = ind[0];
        accum = accum * shape_[1] + ind[1];
        accum = accum * shape_[2] + ind[2];

        return accum;

    }

    size_t index_to_offset( Indices const & ind ) const {

        return index_to_map_index( ind ) * max_hyds_;

    }

    template<class Floats>
    std::vector<Hyd>::const_iterator
    operator[](Floats const & floats) const { 
        Indices ind = floats_to_index(floats);
        return voxel_map_.cbegin() + index_to_offset(ind);
    }

    template<class Floats>
    std::vector<Hyd>::iterator
    operator[](Floats const & floats){
        Indices ind = floats_to_index(floats);
        return voxel_map_.begin() + index_to_offset(ind);
    }

    std::vector<Hyd>::const_iterator 
    at( float f, float g, float h ) const {
        Indices idx = floats_to_index( Bounds( f, g, h ) );
        if( idx[0] < shape_[0] && idx[1] < shape_[1] && idx[2] < shape_[2] )
            return voxel_map_.cbegin() + index_to_offset(idx);
        else return OOB_LIST.cbegin();
    }

    template<class V>
    std::vector<Hyd>::const_iterator
    at( V const & v ) const {
        Indices idx = floats_to_index( Bounds( v[0], v[1], v[2] ) );
        if( idx[0] < shape_[0] && idx[1] < shape_[1] && idx[2] < shape_[2] )
            return voxel_map_.cbegin() + index_to_offset(idx);
        else return OOB_LIST.cbegin();
    }




};




}}

#endif
