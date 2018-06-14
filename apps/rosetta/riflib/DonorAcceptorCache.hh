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
    );

    void
    create_and_fill_voxel_map( std::vector<std::vector<Sat>> const & early_map );


    void
    prepare_bounds( std::vector<HBondRay> const & rays );


    std::vector<std::vector<Sat>>
    first_pass_fill( std::vector<HBondRay> const & rays );



    template<class Floats> Indices floats_to_index(Floats const & f) const {
        Indices ind;
        for(int i = 0; i < 3; ++i){
            float tmp = ((f[i]-lb_[i])/cs_[i]);
            ind[i] = tmp;
        }
        return ind;
    }


    size_t index_to_map_index( Indices const & ind ) const;

    size_t index_to_offset( Indices const & ind ) const;

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
