// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_AtomsCloseTogetherManager_hh
#define INCLUDED_riflib_AtomsCloseTogetherManager_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>

#include "scheme/util/SimpleArray.hh"
#include <riflib/RotamerGenerator.hh>
#include <riflib/ScoreRotamerVsTarget.hh>
#include <scheme/objective/voxel/VoxelArray.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pose/Pose.hh>

#include <string>
#include <vector>

#include <boost/any.hpp>
#include <boost/format.hpp>

#include <riflib/util.hh>

#include <ObjexxFCL/format.hh>


namespace devel {
namespace scheme {



struct AtomsCloseTogetherManager {


    shared_ptr<::scheme::objective::voxel::VoxelArray<3>> voxel_array_;

    std::string primary_info_label_;
    float bonus_;
    float resl_;
    float min_dist_;
    float max_dist_;
    std::vector<std::pair<int, int>> resnum_atomnums_;

    AtomsCloseTogetherManager(
        core::pose::Pose const & target,
        std::string const & spec_string
    )
    {
        std::cout << "Creating AtomsCloseTogether grid with spec: " << spec_string << std::endl;
        parse_spec_string( target, spec_string );
        prepare_bounds( target );
        create_and_fill_voxel_map( target );
    

        // voxel_array_->dump_pdb("CB_too_close.pdb", 0.1, true, 0.1);
    }


    void
    parse_spec_string(
        core::pose::Pose const & target,
        std::string const & spec_string
    ) {
        utility::vector1<std::string> split = utility::string_split( spec_string, ',' );
        if ( split.size() < 6 ) {
            utility_exit_with_message("Not enough fields for -specific_atoms_close_bonus: " + spec_string );
        }
        primary_info_label_ = split[1];
        try {
            bonus_ = utility::from_string( split[2], float(0) );
            resl_ = utility::from_string( split[3], float(0) );
            min_dist_ = utility::from_string( split[4], float(0) );
            max_dist_ = utility::from_string( split[5], float(0) );
        } catch ( std::exception const & ) {
            utility_exit_with_message( "Failed to parse a number in -specific_atoms_close_bonus: " + spec_string);
        }
        runtime_assert( resl_ > 0 );
        runtime_assert( max_dist_ >= min_dist_ );

        for ( size_t i = 6; i <= split.size(); i++ ) {
            utility::vector1<std::string> sp = utility::string_split( split[i], ':' );
            if ( sp.size() != 2 ) {
                utility_exit_with_message("Badly formatted atom selector in -specific_atoms_close_bonus: " + split[i] );
            }
            core::Size resnum = 0;
            try {
                resnum = utility::from_string( sp[1], core::Size(0 ) );
            } catch ( std::exception const & ) {
                utility_exit_with_message( "Failed to parse a number in atom selector in -specific_atoms_close_bonus: " + split[i]);
            }
            if ( resnum <= 0 || resnum > target.size() ) {
                utility_exit_with_message( "Resnum out of bounds in -specific_atoms_close_bonus: " + split[i]);
            }
            core::conformation::Residue const & res = target.residue(resnum);
            if ( ! res.has( sp[2] ) ) {
                utility_exit_with_message( "Atom name doesn't exist on residue in -specific_atoms_close_bonus: " + split[i]);   
            }
            resnum_atomnums_.emplace_back( resnum, res.atom_index( sp[2] ) );
        }
    }

    void
    prepare_bounds( core::pose::Pose const & target  ) {

        Eigen::Vector3f lbs( 9e9, 9e9, 9e9 );
        Eigen::Vector3f ubs( -9e9, -9e9, -9e9 );

        for ( std::pair<int, int> const & resnum_atomnum : resnum_atomnums_ ) {
            core::conformation::Residue const & res = target.residue(resnum_atomnum.first);

            numeric::xyzVector<core::Real> xyz = res.xyz(resnum_atomnum.second);
            for ( int i = 0; i < 3; i++ ) {
                lbs[i] = std::min<float>( lbs[i], xyz[i] );
                ubs[i] = std::max<float>( ubs[i], xyz[i] );
            }
            
        }

        
        for ( int i = 0; i < 3; i++ ) {
            lbs[i] -= max_dist_ + resl_ * 2;
            ubs[i] += max_dist_ + resl_ * 2;
        }

        voxel_array_ = make_shared<::scheme::objective::voxel::VoxelArray<3>>( lbs, ubs, Eigen::Vector3f( resl_, resl_, resl_ ) );

    }

    void
    create_and_fill_voxel_map(
        core::pose::Pose const & target
    ) {

        ::scheme::util::SimpleArray<3,float> lb = voxel_array_->lb_;
        ::scheme::util::SimpleArray<3,float> ub = voxel_array_->ub_;

        float const max_dist2 = max_dist_ * max_dist_;
        float const max_dist_min_min_dist = max_dist_ - min_dist_;

        
        for ( std::pair<int, int> const & resnum_atomnum : resnum_atomnums_ ) {
            core::conformation::Residue const & res = target.residue(resnum_atomnum.first);

            numeric::xyzVector<core::Real> _xyz = res.xyz( resnum_atomnum.second );
            Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

            Eigen::Vector3f lbs( xyz[0] - max_dist_, xyz[1] - max_dist_, xyz[2] - max_dist_ );
            Eigen::Vector3f ubs( xyz[0] + max_dist_, xyz[1] + max_dist_, xyz[2] + max_dist_ );

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

                        if ( squared_dist < max_dist2 ) {
                            float dist = std::sqrt(squared_dist);
                            float score = bonus_;
                            if ( dist > min_dist_ ) {
                                score = bonus_ * ( max_dist_ - dist ) / max_dist_min_min_dist;
                            }
                            if ( std::abs<float>( score ) > std::abs<float>( (*voxel_array_)[worker] ) ) {
                                (*voxel_array_)[worker] = score;
                            }
                        } 
                    }
                }
            }
        }
    }


    std::vector<Eigen::Matrix<float,3,1>>
    prepare_scaffold_atoms(
        core::pose::Pose const & pose
    ) const {
        core::pose::PDBInfoCOP pdb_info = pose.pdb_info();

        std::string search_string = primary_info_label_ + ":";

        std::vector<Eigen::Matrix<float,3,1>> atoms;

        for ( size_t seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
            utility::vector1< std::string > reslabels = pdb_info->get_reslabels( seqpos );
            for ( std::string const & label : reslabels ) {
                if ( utility::startswith( label, search_string ) ) {
                    std::string atom_name = label.substr( search_string.size() );
                    if ( ! pose.residue(seqpos).has( atom_name ) ) {
                        utility_exit_with_message("AtomsCloseTogetherManager with pdbinfolabel " + primary_info_label_ + " can't find atom at"
                            + " position " + utility::to_string(seqpos) + " from label " + label );
                    }
                    numeric::xyzVector<core::Real> xyz = pose.residue(seqpos).xyz( atom_name );
                    atoms.emplace_back( xyz.x(), xyz.y(), xyz.z() );
                }
            }
        }
        return atoms;
    }



/////////////////////////////////////////////////////////////////////////
    float
    get_bonus(
        EigenXform const & scaff_xform,
        std::vector<Eigen::Matrix<float,3,1>> const & atoms
    ) const {
        float bonus = 0;

        for ( Eigen::Matrix<float,3,1> const & atom : atoms ) {
            Eigen::Matrix<float,3,1> pos = scaff_xform * atom;
            bonus += voxel_array_->at( pos );
        }

        return bonus;
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
