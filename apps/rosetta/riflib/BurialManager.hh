// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_BurialManager_hh
#define INCLUDED_riflib_BurialManager_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>

#include <core/pose/Pose.hh>

#include <string>
#include <vector>

#include <boost/any.hpp>



namespace devel {
namespace scheme {


struct BurialOpts {
    float neighbor_distance_cutoff = 6;
    std::vector<float> neighbor_count_weights;
};


struct BurialManager {

    BurialManager() {} // used by clone()

    BurialManager( 
        BurialOpts const & opts,
        std::vector< HBondRay > target_donors,
        std::vector< HBondRay > target_acceptors
    ) :
        opts_( opts )
    {
        donor_acceptors_ = target_donors;
        donor_acceptors_.insert( donor_acceptors_.end(), target_acceptors.begin(), target_acceptors.end());
        num_donors_ = target_donors.size();

        target_neighbor_counts_.resize( donor_acceptors_.size(), 0 );
        other_neighbor_counts_.resize( donor_acceptors_.size(), 0 );

        runtime_assert( opts_.neighbor_count_weights.size() >= 20 );
    }

    shared_ptr<BurialManager>
    clone() const;

    void
    reset();

    void
    set_target_neighbors( core::pose::Pose const & pose );

    std::vector<float>
    get_burial_weights( ) const;

    void
    accumulate_neighbors( BBActor const & bb );

private:

    BurialOpts opts_;

    int num_donors_;
    std::vector< HBondRay > donor_acceptors_;
    std::vector<int> target_neighbor_counts_;
    std::vector<int> other_neighbor_counts_;

};




}}

#endif
