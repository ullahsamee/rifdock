// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_CBTooCloseManager_hh
#define INCLUDED_riflib_CBTooCloseManager_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>

#include "scheme/util/SimpleArray.hh"
#include <riflib/RotamerGenerator.hh>
#include <riflib/ScoreRotamerVsTarget.hh>
#include <scheme/objective/voxel/VoxelArray.hh>

#include <core/pose/Pose.hh>

#include <string>
#include <vector>

#include <boost/any.hpp>
#include <boost/format.hpp>

#include <riflib/util.hh>

#include <ObjexxFCL/format.hh>


namespace devel {
namespace scheme {



struct CBTooCloseManager {

    float resl_;

    shared_ptr<::scheme::objective::voxel::VoxelArray<3>> voxel_array_;




    CBTooCloseManager(
        core::pose::Pose const & target,
        float resl,
        float too_close_dist,
        float penalty,
        size_t max_target_res_atom_idx
    ) : resl_( resl )
    {
        std::cout << "Creating CB too close grid" << std::endl;
        prepare_bounds( target, too_close_dist, max_target_res_atom_idx );
        create_and_fill_voxel_map( target, too_close_dist, penalty, max_target_res_atom_idx );
    

        // voxel_array_->dump_pdb("CB_too_close.pdb", 0.1, true, 0.1);
    }


    void
    prepare_bounds( core::pose::Pose const & target, float too_close_dist, size_t max_target_res_atom_idx ) {

        Eigen::Vector3f lbs( 9e9, 9e9, 9e9 );
        Eigen::Vector3f ubs( -9e9, -9e9, -9e9 );

        for ( core::Size seqpos = 1; seqpos <= target.size(); seqpos++ ) {
            core::conformation::Residue const & res = target.residue(seqpos);
            size_t loop_ub = std::min<size_t>( res.nheavyatoms(), max_target_res_atom_idx );
            for ( core::Size atno = 1; atno <= loop_ub; atno++ ) { 

                numeric::xyzVector<core::Real> xyz = res.xyz(atno);
                for ( int i = 0; i < 3; i++ ) {
                    lbs[i] = std::min<float>( lbs[i], xyz[i] );
                    ubs[i] = std::max<float>( ubs[i], xyz[i] );
                }
            }
            
        }

        
        for ( int i = 0; i < 3; i++ ) {
            lbs[i] -= too_close_dist + resl_ * 2;
            ubs[i] += too_close_dist + resl_ * 2;
        }

        // lb_ = lbs;
        // ub_ = ubs;
        // cs_ = Eigen::Vector3f( resl_, resl_, resl_ );

        voxel_array_ = make_shared<::scheme::objective::voxel::VoxelArray<3>>( lbs, ubs, Eigen::Vector3f( resl_, resl_, resl_ ) );

        // Indices extents = floats_to_index( ub_ );
        // shape_ = extents + Indices(1);

    }

    void
    create_and_fill_voxel_map(
        core::pose::Pose const & target, 
        float too_close_dist,
        float penalty,
        size_t max_target_res_atom_idx
        ) {

        // size_t elements = shape_[0] * shape_[1] * shape_[2];
        // std::cout << shape_ << std::endl;
        // std::cout << floats_to_index( ub_ ) << std::endl;
        // std::cout << elements << " " << index_to_map_index( floats_to_index( ub_ ) ) << std::endl;
        // runtime_assert( elements - 1 == index_to_map_index( floats_to_index( ub_ ) ) );

        // voxel_map_.resize( elements, 0 );


        ::scheme::util::SimpleArray<3,float> lb = voxel_array_->lb_;
        ::scheme::util::SimpleArray<3,float> ub = voxel_array_->ub_;

        float const too_close_dist2 = too_close_dist * too_close_dist;


        for ( size_t seqpos = 1; seqpos <= target.size(); seqpos ++  ) {
            core::conformation::Residue const & res = target.residue(seqpos);

            size_t loop_ub = std::min<size_t>( res.nheavyatoms(), max_target_res_atom_idx );
            for ( core::Size atno = 1; atno <= loop_ub; atno++ ) {  

                numeric::xyzVector<core::Real> _xyz = res.xyz( atno );
                Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

                Eigen::Vector3f lbs( xyz[0] - too_close_dist, xyz[1] - too_close_dist, xyz[2] - too_close_dist );
                Eigen::Vector3f ubs( xyz[0] + too_close_dist, xyz[1] + too_close_dist, xyz[2] + too_close_dist );

                const float step = resl_;

                Eigen::Vector3f worker;

                for ( float x = lbs[0] - step/2; x < ubs[0] + step; x += step ) {
                    if ( x < lb[0] || x > ub[0] ) continue;
                    worker[0] = x;

                    for ( float y = lbs[1] - step/2; y < ubs[1] + step; y += step ) {
                        if ( y < lb[1] || y > ub[1] ) continue;
                        worker[1] = y;

                        for ( float z = lbs[2] - step/2; z < ubs[2] + step; z += step ) {
                            if ( z < lb[2] || z > ub[2] ) continue;
                            worker[2] = z;

                            const float squared_dist = ( xyz - worker ).squaredNorm();

                            if ( squared_dist < too_close_dist2 ) {
                                (*voxel_array_)[worker] = penalty;
                                // size_t offset = index_to_map_index( floats_to_index( worker ) );
                                // voxel_map_.at(offset) = penalty;
                            } 
                        }
                    }
                }
            }
        }
    }




/////////////////////////////////////////////////////////////////////////
    float
    get_CB_penalty(
        EigenXform const & bbpos
    ) const {

        Eigen::Matrix<float,3,1> CB = bbpos * Eigen::Matrix<float,3,1>( 1.0264273 ,  0.25245885, -0.308907 );

        return voxel_array_->at( CB );

    }




    // template<class Floats> Indices floats_to_index(Floats const & f) const {
    //     Indices ind;
    //     for(int i = 0; i < 3; ++i){
    //         float tmp = ((f[i]-lb_[i])/cs_[i]);
    //         ind[i] = tmp;
    //     }
    //     return ind;
    // }

    // size_t index_to_map_index( Indices const & ind ) const {

    //     size_t accum = ind[0];
    //     accum = accum * shape_[1] + ind[1];
    //     accum = accum * shape_[2] + ind[2];

    //     return accum;

    // }

    // size_t index_to_offset( Indices const & ind ) const {

    //     return index_to_map_index( ind );

    // }



    // float const &
    // at( float f, float g, float h ) const {
    //     Indices idx = floats_to_index( Bounds( f, g, h ) );
    //     if( idx[0] < shape_[0] && idx[1] < shape_[1] && idx[2] < shape_[2] )
    //         return voxel_map_.at( index_to_offset(idx) );
    //     else return 0;
    // }

    // template<class V>
    // float const &
    // at( V const & v ) const {
    //     Indices idx = floats_to_index( Bounds( v[0], v[1], v[2] ) );
    //     if( idx[0] < shape_[0] && idx[1] < shape_[1] && idx[2] < shape_[2] )
    //         return voxel_map_.at( index_to_offset(idx) );
    //     else return 0;
    // }


};




}}

#endif
