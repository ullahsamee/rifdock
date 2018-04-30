// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_util_hh
#define INCLUDED_riflib_util_hh

#include <riflib/types.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <core/pose/Pose.hh>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <numeric/xyzTransform.hh>
#include <exception>

#include <utility/io/ozstream.fwd.hh>
#include <utility/io/izstream.fwd.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/id/AtomID.hh>


#ifdef USEGRIDSCORE
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <protocols/ligand_docking/GALigandDock/RotamerData.hh>
#endif

namespace devel {
namespace scheme {



#ifdef USE_OPENMP
 #include <omp.h>
 #endif
 static core::Size omp_max_threads_1(){
	#ifdef USE_OPENMP
		return omp_get_max_threads();
	#else
		return 1;
	#endif
 }
 static core::Size omp_thread_num_1(){
	#ifdef USE_OPENMP
		return omp_get_thread_num() + 1;
	#else
		return 1;
	#endif
 }
 static core::Size omp_max_threads(){
	#ifdef USE_OPENMP
		return omp_get_max_threads();
	#else
		return 1;
	#endif
 }
 static core::Size omp_thread_num(){
	#ifdef USE_OPENMP
		return omp_get_thread_num();
	#else
		return 0;
	#endif
 }

utility::vector1<core::Size> get_res(
	std::string fname,
	core::pose::Pose const & pose,
	bool nocgp = true
 );

utility::vector1<core::Size> get_designable_positions_best_guess(
	  core::pose::Pose pose
	, bool noloops
	, bool nocpg = true
 );


std::string
get_res_list_hash( utility::vector1<core::Size> const & reslist);


template<class T>
std::string str(T const & t, core::Size N=0){
	// std::ostringstream oss;
	// oss << t;
	std::string s = boost::lexical_cast<std::string>(t);
	while(s.size()<N) s = "0"+s;
	return s;
}


void
apply_xform_to_pose( core::pose::Pose & pose, EigenXform const & xform);

void
apply_xform_to_residue( core::conformation::Residue & residue, EigenXform const & xform);



void
append_pose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	bool new_chain = true
);

/// @brief Append specified residues of pose2 to pose1.
void
append_subpose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	core::Size start_res,
	core::Size end_res,
	bool new_chain = true
);


std::vector<int> get_rif_atype_map();



template<class Float>
Eigen::Matrix<Float,3,3>
xyz2eigen( numeric::xyzMatrix<Float> const & m ){
	Eigen::Matrix<Float,3,3> rot;
	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			rot(i,j) = m(i+1,j+1);
		}
	}
	return rot;
}

template< class Float >
numeric::xyzMatrix<Float>
eigen2xyz( Eigen::Matrix<Float,3,3> const & rot ){
	numeric::xyzMatrix<Float> m;
	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			m(i+1,j+1) = rot(i,j);
		}
	}
	return m;
}

template< class Float >
numeric::xyzTransform<Float>
eigen2xyz( Eigen::Transform<Float,3,Eigen::AffineCompact> const & xin ){
	numeric::xyzTransform<Float> x( eigen2xyz( xin.rotation() ) );
	x.t[0] = xin.translation()[0];
	x.t[1] = xin.translation()[1];
	x.t[2] = xin.translation()[2];
	return x;
}

template< class Float >
void print_eigenxform( Eigen::Transform<Float,3,Eigen::AffineCompact> const & x, std::ostream & out = std::cout ){
	out << x.rotation() << "    " << x.translation().transpose() << std::endl;
}

static void print_header( std::string s, int n=120, int n2=0 ){
	for( int j = 0; j < n2; ++j )
		for( int i = 0; i < n; ++i )
			std::cout << "="; std::cout << std::endl;
	for( int i = 0; i < n/2-(int)s.size()/2-1; ++i )
		std::cout << "=";
	std::cout << " " << s << " ";
	for( int i = 0; i < n/2-(int)s.size()+(int)s.size()/2-1; ++i )
		std::cout << "=";
	std::cout << std::endl;
	for( int j = 0; j < n2; ++j )
		for( int i = 0; i < n; ++i )
			std::cout << "="; std::cout << std::endl;
}

std::string KMGT(double const & x, int const & w=7, int const & d=3);

template<class EigenXform>
float xform_magnitude(
	EigenXform const & x,
	float rg
){
	float err_trans2 = x.translation().squaredNorm();
	float cos_theta = (x.rotation().trace()-1.0)/2.0;

	// sanity check...
	// float cos_theta_slow = cos( Eigen::AngleAxisf( x.rotation() ).angle() );
	// if( fabs( cos_theta_slow-cos_theta ) > 0.001 ){
	// 	std::cout << "ang " << Eigen::AngleAxisf( x.rotation() ).angle()*180.0/M_PI << " " << cos_theta << " " << cos_theta_slow << std::endl;
	// }

	float err_rot = sqrt( std::max( 0.0, 1.0 - cos_theta*cos_theta ) ) * rg;
	if( cos_theta < 0 ) err_rot = rg;
	float err = sqrt( err_trans2 + err_rot*err_rot );
	return err;
}

void pose_to_ala( core::pose::Pose & pose );
void pose_to_gly( core::pose::Pose & pose );
void pose_to_ala( core::pose::Pose & pose, utility::vector1<core::Size> const & res_sel );
void pose_to_val( core::pose::Pose & pose );


std::string
open_for_read_on_path(
	std::vector<std::string> const & path,
	std::string fname,
	utility::io::izstream & in
);

std::string
open_for_write_on_path(
	std::vector<std::string> const & path,
	std::string fname,
	utility::io::ozstream & out,
	bool create_directorys = false
);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// OMG! MOVE ME
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< class Float > Float
sqr( Float const & x ) { return x*x; }




// template< class HBondRay >
// float score_hbond_rays(
// 	HBondRay const & don,
// 	HBondRay const & acc,
// 	float non_directional_fraction = 0.0 // 0.0 - 1.0
// ){
// 	Eigen::Vector3f ho = acc.horb_cen - don.horb_cen;
// 	float ho_dist = ho.norm() - 1.05 ;
// 	ho_dist = ho_dist < 0 ? ho_dist*1.5 : ho_dist; // increase dis pen if too close
// 	float const max_diff = 1.3;
// 	ho_dist = ho_dist >  max_diff ?  max_diff : ho_dist;
// 	ho_dist = ho_dist < -max_diff ? -max_diff : ho_dist;

// 	ho.normalize();
// 	float bho_dot = don.direction.dot(ho);
// 	float hoa_dot = -ho.dot(acc.direction);
	
// 	if( hoa_dot < 0 ) hoa_dot = 0;
// 	if( bho_dot < 0 ) bho_dot = 0;

// 	// if( ho_dist > max_diff ) return 0.0;
// 	// sigmoid -like shape on distance score
// 	float score = -sqr( 1.0 - sqr( ho_dist/max_diff ) );

// 	// score -= bho_dot;
// 	// score -= hoa_dot;
// 	// score /= 3.0;

// 	score *= bho_dot;
// 	score *= hoa_dot;

// 	return score;
// }




/// @brief Returns an map mapping the CA of the atoms given in the subset
/// to themselves
std::map< core::id::AtomID, core::id::AtomID >
residue_subset_to_CA_atom_map( 
    core::select::residue_selector::ResidueSubset const & subset, 
    core::pose::Pose const & pose );


/// @brief Returns the CA rmsd of the residues in the ResidueSubset
/// with or without superposition
core::Real
subset_CA_rmsd(
    core::pose::Pose const & pose1,
    core::pose::Pose const & pose2,
    core::select::residue_selector::ResidueSubset const & subset,
    bool superimpose);


/// Your poses must be identical. Very important
EigenXform
find_xform_from_identical_pose_to_pose( 
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	float align_error = 0.2 );


void
all_by_all_rmsd( 
	std::vector<core::pose::PoseOP> const & poses,
	utility::vector1<utility::vector1<core::Real>> & table );

std::vector<std::vector<std::pair<core::pose::PoseOP, uint64_t>>>
cluster_poses_into_n_bins( 
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n,
	utility::vector1<utility::vector1<core::Real>> & rmsds );

std::vector<core::pose::PoseOP>
cluster_poses_leaving_n( 
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n );

size_t
find_cluster_center( 
	std::vector<uint64_t> indexes,
	utility::vector1<utility::vector1<core::Real>> const & rmsds );

std::vector<core::pose::PoseOP>
cluster_poses_leaving_n_representing_frac(
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n,
	float frac,
	float tol
	);

std::vector<core::pose::PoseOP>
random_selection_poses_leaving_n( 
	std::vector<core::pose::PoseOP> const & poses,
	uint64_t n );


// from riflib/rifdock_subroutines/util.hh


Eigen::Vector3f
pose_center(
    core::pose::Pose const & pose,
    utility::vector1<core::Size> const & useres = utility::vector1<core::Size>()
);



void
get_rg_radius(
    core::pose::Pose const & pose,
    float & rg,
    float & radius,
    utility::vector1<core::Size> const & useres = utility::vector1<core::Size>(),
    bool allatom = false
);


float
angle_between_vectors(
	Eigen::Vector3f const & vec1,
	Eigen::Vector3f const & vec2
);

void 
xform_pose( 
	core::pose::Pose & pose, 
	numeric::xyzTransform<float> s, 
	core::Size sres=1, 
	core::Size eres=0
);

#ifdef USEGRIDSCORE
shared_ptr<protocols::ligand_docking::ga_ligand_dock::GridScorer>
prepare_grid_scorer(
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	std::string const & atype_aas = "ACDEFGHIKLMNPQRSTVWY"
);
#endif


// https://stackoverflow.com/questions/7931358/printing-sizeoft-at-compile-time

// usage: PRINT_SIZE_AS_WARNING<sizeof(int)>()();
template<int N> 
struct PRINT_SIZE_AS_ERROR
{ 
   int operator()() { return (int)std::vector<int>(); } //deliberately causing error
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
}

#endif
