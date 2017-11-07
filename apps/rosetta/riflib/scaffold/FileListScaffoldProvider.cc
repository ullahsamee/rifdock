// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/scaffold/FileListScaffoldProvider.hh>

#include <riflib/types.hh>
#include <scheme/numeric/rand_xform.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/file/file_sys_util.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>



namespace devel {
namespace scheme {


// FileListScaffoldProvider::FileListScaffoldProvider() {}

FileListScaffoldProvider::FileListScaffoldProvider( 
        std::vector<std::string> const & scaff_fnames, 
        std::vector<utility::vector1<core::Size>> const & scaffold_ress,
        shared_ptr< RotamerIndex > rot_index_p_in, 
        RifDockOpt const & opt_in ) :
        scaff_fnames_( scaff_fnames ),
        scaffold_ress_( scaffold_ress ),
        rot_index_p( rot_index_p_in), 
        opt(opt_in) {
            assert( scaff_fnames.size() == scaffold_ress.size() );
            temp__data_caches_.resize( scaff_fnames.size() );
        }


ScaffoldDataCacheOP 
FileListScaffoldProvider::temp_function__get_data_cache(uint64_t i) {
    assert(i < get_scaffold_index_limits());

    std::string scafftag = utility::file_basename( utility::file::file_basename( scaff_fnames_[i] ) );
    core::pose::Pose scaffold;
    core::import_pose::pose_from_file( scaffold, scaff_fnames_[i] );

    return make_shared<ScaffoldDataCache>( 
        scaffold, 
        scaffold_ress_[i],
        scafftag,
        rot_index_p,
        opt);
}


// FileListScaffoldProvider::FileListScaffoldProvider(
//         std::string const & scaff_fname, 
//         std::string const & scaff_res_fname,
//         shared_ptr< RotamerIndex > rot_index_p_in, 
//         RifDockOpt const & opt_in ) :
//         rot_index_p( rot_index_p_in), opt(opt_in) {

//     std::string scafftag = utility::file_basename( utility::file::file_basename( scaff_fname ) );
//     core::pose::Pose scaffold;
//     core::import_pose::pose_from_file( scaffold, scaff_fname );

//     if( opt.random_perturb_scaffold ){
//         EigenXform scaffold_perturb = EigenXform::Identity();
//         runtime_assert_msg( !opt.use_scaffold_bounding_grids,
//             "opt.use_scaffold_bounding_grids incompatible with random_perturb_scaffold" );

//         std::mt19937 rng( 0);// std::random_device{}() );
//         ::scheme::numeric::rand_xform(rng,scaffold_perturb);
//         xform_pose( scaffold, eigen2xyz(scaffold_perturb) );
//     }

//     // setup scaffold_res
//     utility::vector1<core::Size> scaffold_res;
//     if( opt.scaffold_res_fnames.size() ){
//         if( opt.scaffold_res_use_best_guess ){
//             utility_exit_with_message("should only use -scaffold_res_use_best_guess true iff not specifying scaffold_res");
//         }
//         scaffold_res = devel::scheme::get_res( scaff_res_fname , scaffold );
//     } else if (opt.scaffold_res_use_best_guess ){
//         scaffold_res = devel::scheme::get_designable_positions_best_guess( scaffold, opt.dont_use_scaffold_loops );
//         std::cout << "using scaffold residues: ";
//         for(auto ir:scaffold_res) std::cout << " " << ir << scaffold.residue(ir).name3();
//         std::cout << std::endl;
//     } else {
//         for( int ir = 1; ir <= scaffold.size(); ++ir){
//             if( !scaffold.residue(ir).is_protein() ) continue;
//             //if( scaffold.residue(ir).name3() == "PRO" ) continue;
//             //if( scaffold.residue(ir).name3() == "GLY" ) continue;
//             //if( scaffold.residue(ir).name3() == "CYS" ) continue;
//             scaffold_res.push_back(ir);
//         }
//     }

//     data_cache_ = make_shared<ScaffoldDataCache>( 
//         scaffold, 
//         scaffold_res,
//         scafftag,
//         rot_index_p,
//         opt);
// }


ParametricSceneConformationCOP 
FileListScaffoldProvider::get_scaffold(uint64_t i) {
    if ( ! conformation_ ) {
        utility_exit_with_message("Conformation not intialized yet!!");
    }
    return conformation_;
}


uint64_t 
FileListScaffoldProvider::get_scaffold_index_limits() {
    return scaff_fnames_.size();
}



}}

