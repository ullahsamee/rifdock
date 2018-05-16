// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_DonorAcceptorCache_hh
#define INCLUDED_riflib_DonorAcceptorCache_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>

#include "scheme/util/SimpleArray.hh"

#include <core/pose/Pose.hh>

#include <string>
#include <vector>

#include <boost/any.hpp>
#include <boost/format.hpp>


namespace devel {
namespace scheme {


struct DonorAcceptorCache {

    typedef uint16_t Sat;
    static int const CACHE_MAX_SAT = std::numeric_limits<Sat>::max();

    typedef ::scheme::util::SimpleArray<3,size_t> Indices;
    typedef Eigen::Vector3f Bounds;
    Bounds lb_,ub_,cs_;

    float max_interaction_range_;
    size_t max_sats_;

    std::vector<Sat> voxel_map_;

    // This only actually needs 1. But include a few extra in case it accidentally gets written to
    std::vector<Sat> OOB_LIST = {CACHE_MAX_SAT, CACHE_MAX_SAT, CACHE_MAX_SAT, CACHE_MAX_SAT, CACHE_MAX_SAT, CACHE_MAX_SAT};

    Indices shape_;

    DonorAcceptorCache(
        std::vector<HBondRay> const & rays,
        float max_interaction_range
    ) : max_interaction_range_(max_interaction_range)
    {

        std::cout << "Preparing DonorAcceptorCache with interaction range " 
                  << boost::str(boost::format("%.1f...")%max_interaction_range) << std::endl;

        prepare_bounds( rays );
        std::vector<std::vector<Sat>> early_map = first_pass_fill( rays );

        create_and_fill_voxel_map( early_map );

    }

    void
    create_and_fill_voxel_map( std::vector<std::vector<Sat>> const & early_map ) {

        max_sats_ = 0;

        for ( std::vector<Sat> these_sats : early_map ) {
            max_sats_ = std::max<size_t>( max_sats_, these_sats.size() );
        }

        std::cout << "Max sats at one voxel: " << max_sats_ << std::endl;

        // max_sats_ is 1 too big. This means that the receiver can always loop until they hit a CACHE_MAX_SAT
        max_sats_ += 1;

        voxel_map_.resize( max_sats_ * shape_[0] * shape_[1] * shape_[2], CACHE_MAX_SAT );

        // We already asserted the two maps have the same size, so just array fill
        for ( size_t imap = 0; imap < early_map.size(); imap++ ) {

            std::vector<Sat> these_sats = early_map[imap];
            bool wrote_a_null = false;
            for ( size_t idx = 0; idx < max_sats_; idx++ ) {
                size_t offset = imap * max_sats_ + idx;
                if ( idx < these_sats.size() ) {
                    voxel_map_.at(offset) = these_sats.at(idx);
                } else {
                    voxel_map_.at(offset) = CACHE_MAX_SAT;
                    wrote_a_null = true;
                }
            }
            runtime_assert( wrote_a_null );     // Redundancy. Make sure each list is terminated
        }

    }


    void
    prepare_bounds( std::vector<HBondRay> const & rays ) {

        Eigen::Vector3f lbs( 9e9, 9e9, 9e9 );
        Eigen::Vector3f ubs( -9e9, -9e9, -9e9  );

        for ( HBondRay const & ray : rays ) {
            Eigen::Vector3f xyz = ray.horb_cen;
            for ( int i = 0; i < 3; i++ ) {
                lbs[i] = std::min<float>( lbs[i], xyz[i] );
                ubs[i] = std::max<float>( ubs[i], xyz[i] );
            }
        }

        for ( int i = 0; i < 3; i++ ) {
            lbs[i] -= max_interaction_range_;
            ubs[i] += max_interaction_range_;
        }

        lb_ = lbs;
        ub_ = ubs;
        cs_ = Eigen::Vector3f( 1.0, 1.0, 1.0 );

        Indices extents = floats_to_index( ub_ );
        shape_ = extents + Indices(1);

    }



    std::vector<std::vector<Sat>>
    first_pass_fill( std::vector<HBondRay> const & rays ) {

        std::vector<std::vector<Sat>> early_map;

        size_t elements = shape_[0] * shape_[1] * shape_[2];
        std::cout << shape_ << std::endl;
        std::cout << floats_to_index( ub_ ) << std::endl;
        std::cout << elements << " " << index_to_map_index( floats_to_index( ub_ ) ) << std::endl;
        runtime_assert( elements - 1 == index_to_map_index( floats_to_index( ub_ ) ) );

        early_map.resize(elements);

        for ( Sat isat = 0; isat < rays.size(); isat++ ) {
            HBondRay const & ray = rays[isat];

            Eigen::Vector3f xyz = ray.horb_cen + ray.direction * 0.5;

            Eigen::Vector3f lbs( xyz[0] - max_interaction_range_, xyz[1] - max_interaction_range_, xyz[2] - max_interaction_range_ );
            Eigen::Vector3f ubs( xyz[0] + max_interaction_range_, xyz[1] + max_interaction_range_, xyz[2] + max_interaction_range_ );

            const float radius_sq = max_interaction_range_*max_interaction_range_;
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

                        if ( ( xyz - worker ).squaredNorm() <= radius_sq ) {

                            size_t offset = index_to_map_index( floats_to_index( worker ) );
                            early_map.at(offset).push_back( isat );
                        }
                    }
                }
            }
        }

        return early_map;

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

        return index_to_map_index( ind ) * max_sats_;

    }

    template<class Floats>
    std::vector<Sat>::const_iterator
    operator[](Floats const & floats) const { 
        Indices ind = floats_to_index(floats);
        return voxel_map_.cbegin() + index_to_offset(ind);
    }

    template<class Floats>
    std::vector<Sat>::iterator
    operator[](Floats const & floats){
        Indices ind = floats_to_index(floats);
        return voxel_map_.begin() + index_to_offset(ind);
    }

    std::vector<Sat>::const_iterator 
    at( float f, float g, float h ) const {
        Indices idx = floats_to_index( Bounds( f, g, h ) );
        if( idx[0] < shape_[0] && idx[1] < shape_[1] && idx[2] < shape_[2] )
            return voxel_map_.cbegin() + index_to_offset(idx);
        else return OOB_LIST.cbegin();
    }

    template<class V>
    std::vector<Sat>::const_iterator
    at( V const & v ) const {
        Indices idx = floats_to_index( Bounds( v[0], v[1], v[2] ) );
        if( idx[0] < shape_[0] && idx[1] < shape_[1] && idx[2] < shape_[2] )
            return voxel_map_.cbegin() + index_to_offset(idx);
        else return OOB_LIST.cbegin();
    }




};




}}

#endif
