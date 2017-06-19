// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#include <riflib/rif/RifGeneratorUserHotspots.hh>


	#include <ObjexxFCL/format.hh>

	#include <boost/random/mersenne_twister.hpp>
	#include <boost/random/uniform_real.hpp>

	#include <core/id/AtomID.hh>
	#include <core/pose/Pose.hh>
	#include <core/scoring/motif/util.hh>
	#include <core/import_pose/import_pose.hh>

  	#include <numeric/xyzMatrix.hh>
	#include <devel/init.hh>
	#include <riflib/RotamerGenerator.hh>
	#include <riflib/util.hh>
	#include <scheme/numeric/rand_xform.hh>

	#include <scheme/actor/Atom.hh>
	#include <Eigen/SVD>
	#include <Eigen/Core>
	#include <Eigen/Geometry>	

namespace devel {
namespace scheme {
namespace rif {



	void
	RifGeneratorUserHotspots::generate_rif(
		RifAccumulatorP accumulator,
		RifGenParamsP params
	){
		
		typedef ::Eigen::Matrix<float,3,1> Pos;
	
		// some sanity checks
		int const n_hspot_groups = this->opts.hotspot_files.size();
		runtime_assert_msg( n_hspot_groups, "no hotspot group files specified!!" );
		runtime_assert_msg( n_hspot_groups<16, "too many hotspot groups!!" );

		std::cout << "this RIF type doesn't support sat groups!!!" << std::endl;

		std::cout << "RifGeneratorUserHotspots opts:" << std::endl;
		std::cout << "    hotspot_sample_cart_bound:  " << this->opts.hotspot_sample_cart_bound;
		std::cout << "    hotspot_sample_angle_bound: " << this->opts.hotspot_sample_angle_bound;
		std::cout << "    hbond_weight:               " << this->opts.hbond_weight;
		std::cout << "    upweight_multi_hbond:       " << this->opts.upweight_multi_hbond;
		std::cout << "    target_center:              "
			<< this->opts.target_center[0] << " "
			<< this->opts.target_center[1] << " "
			<< this->opts.target_center[2] << std::endl;
		//translation to apply to input hotspot files
		Pos target_vec;
		numeric::xyzVector<double> xyz_tgt_cen( this->opts.target_center[0], this->opts.target_center[1], this->opts.target_center[2] );
		target_vec << opts.target_center[0], opts.target_center[1],opts.target_center[2];


		for( auto s : this->opts.hotspot_files ){
			std::cout << "    hotspot_group:              " << s << std::endl;
		}


		// setup the hacky but fast scorer
		devel::scheme::ScoreRotamerVsTarget<
				VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
			> rot_tgt_scorer;
		{
			std::vector< ::scheme::chemical::HBondRay > target_donors, target_acceptors;
			for( auto ir : params->target_res ){
				::devel::scheme::get_donor_rays   ( *params->target, ir, params->hbopt, target_donors );
				::devel::scheme::get_acceptor_rays( *params->target, ir, params->hbopt, target_acceptors );
			}
			std::cout << "target_donors.size() " << target_donors.size() << " target_acceptors.size() " << target_acceptors.size() << std::endl;
			{
				rot_tgt_scorer.rot_index_p_ = params->rot_index_p;
				rot_tgt_scorer.target_field_by_atype_ = params->field_by_atype;
				rot_tgt_scorer.target_donors_ = target_donors;
				rot_tgt_scorer.target_acceptors_ = target_acceptors;
				rot_tgt_scorer.hbond_weight_ = this->opts.hbond_weight;
				rot_tgt_scorer.upweight_multi_hbond_ = this->opts.upweight_multi_hbond;
				rot_tgt_scorer.upweight_iface_ = 1.0;

			}
		}



    	int const NSAMP = 1000*1000;

    	std::mt19937 rng((unsigned int)time(0) + 934875);
    	float const radius_bound = this->opts.hotspot_sample_cart_bound;
    	float const degrees_bound = this->opts.hotspot_sample_angle_bound;
    	float const radians_bound = degrees_bound * M_PI/180.0;

		// loop over files (one file is one hotspot group)
		for( int i_hotspot_group = 0; i_hotspot_group < this->opts.hotspot_files.size(); ++i_hotspot_group ){

			std::string const & hotspot_file = this->opts.hotspot_files[i_hotspot_group];
			std::cout << hotspot_file << std::endl;
			// read hotspot file into pose
			core::pose::Pose pose;
			core::import_pose::pose_from_file(pose,hotspot_file);

      
			// read in pdb files # i_hotspot_group
			for( int i_hspot_res = 1; i_hspot_res <= pose.size(); ++i_hspot_res ){
				int input_nheavy = pose.residue(i_hspot_res).nheavyatoms();
				EigenXform Xref = ::scheme::chemical::make_stub<EigenXform>(
		            pose.residue(i_hspot_res).xyz( input_nheavy - 2 ) - xyz_tgt_cen,
		            pose.residue(i_hspot_res).xyz( input_nheavy - 1 ) - xyz_tgt_cen,
		            pose.residue(i_hspot_res).xyz( input_nheavy - 0 ) - xyz_tgt_cen
		        );
				
				//get last atom in hotspot residue and iterate over last 3 
      			core::conformation::Atoms::const_iterator iter = pose.residue(i_hspot_res).heavyAtoms_end();
				core::conformation::Atoms::const_iterator end = iter-3;
      	
				//these are the hotspot atoms				
				Pos atom1; Pos atom2; Pos atom3;
				
				//iterate and save xyz into hotspot 
				while(iter >= end){
        			if (iter == end +2){atom1(0,0) = iter->xyz()[0];atom1(1,0) = iter->xyz()[1];atom1(2,0) = iter->xyz()[2];}
        			if (iter == end +1){atom2(0,0) = iter->xyz()[0];atom2(1,0) = iter->xyz()[1];atom2(2,0) = iter->xyz()[2];}
        			if (iter == end +0){atom3(0,0) = iter->xyz()[0];atom3(1,0) = iter->xyz()[1];atom3(2,0) = iter->xyz()[2];}
        			iter --;
      			}
				//translate with the target first
				atom1 = atom1 - target_vec; 
				atom2 = atom2 - target_vec; 
				atom3 = atom3 - target_vec;				

				//calculate centroid of hot_spot res 
				Pos cen_hot = (atom1 + atom2 + atom3)/3;
				
				// for each irot that is the right restype (can be had from rot_intex_p)
				int irot_begin = 0, irot_end = params -> rot_index_p -> size();
				
				for( int irot = irot_begin; irot < irot_end; ++irot ){

					//std::cout << params -> rot_index_p -> resname(irot) << std::endl;
					//std::cout << pose.residue(i_hspot_res).name3() << std::endl;	

					//if (params -> rot_index_p -> resname(irot) == pose.residue(i_hspot_res).name3()){
				  	//std::cout << irot << std::endl;
					::Eigen::Matrix<float,3,3> rif_res; // this is the rif residue last three atoms
						
						// assign rif_res position by column
						int hatoms = params -> rot_index_p -> nheavyatoms(irot);
					std::vector<SchemeAtom> const & rotamer_atoms( params->rot_index_p->atoms(irot) );
					EigenXform Xrotamer = ::scheme::chemical::make_stub<EigenXform>(
                        rotamer_atoms.at( params->rot_index_p->nheavyatoms(irot) - 3 ).position(),
                        rotamer_atoms.at( params->rot_index_p->nheavyatoms(irot) - 2 ).position(),
                        rotamer_atoms.at( params->rot_index_p->nheavyatoms(irot) - 1 ).position()
                    );
					//std::cout << params -> rot_index_p -> resname(irot) << std::endl;
					//std::cout << pose.residue(i_hspot_res).name3() << std::endl;	
					//std::vector<SchemeAtom> const & rotamer_atoms( params->rot_index_p->atoms(irot) );
					if (params -> rot_index_p -> resname(irot) == pose.residue(i_hspot_res).name3()){
				  		

				  		//std::cout << irot << std::endl;
						::Eigen::Matrix<float,3,3> rif_res; // this is the rif residue last three atoms
						
						
						//Eigen::Matrix<float, 3, 1> center_mass;
						// for (int hatom = 0; hatom < hatoms; ++hatom){
						// 		center_mass += rotamer_atoms[hatom].position();
						// }
						// center_mass /= hatoms;
						

						int latoms = params -> rot_index_p -> natoms(irot);
						rif_res << rotamer_atoms[hatoms-1].position()[0],rotamer_atoms[hatoms-2].position()[0],rotamer_atoms[hatoms-3].position()[0],rotamer_atoms[hatoms-1].position()[1],rotamer_atoms[hatoms-2].position()[1],rotamer_atoms[hatoms-3].position()[1],rotamer_atoms[hatoms-1].position()[2],rotamer_atoms[hatoms-2].position()[2],rotamer_atoms[hatoms-3].position()[2];

	     	 			Pos cen_rot = (rif_res.col(0) + rif_res.col(1) + rif_res.col(2))/3;
         				
         				EigenXform x_2_orig = EigenXform::Identity();
         				x_2_orig.translation() = -cen_hot;
         				EigenXform x_2_orig_inverse = x_2_orig.inverse();

						
						//svd superimpose
          				::Eigen::Matrix<float,3,3> cov_mtx;										
          				cov_mtx = (rif_res.col(0) - cen_rot)*(atom1 - cen_hot).transpose() + (rif_res.col(1) - cen_rot)*(atom2 - cen_hot).transpose() + (rif_res.col(2) - cen_rot)*(atom3 - cen_hot).transpose();
            			::Eigen::JacobiSVD<::Eigen::Matrix<float,3,3>> svd(cov_mtx, Eigen::ComputeFullU | Eigen::ComputeFullV);
						::Eigen::Matrix<float,3,3> R_mtx = svd.matrixV() * svd.matrixU().transpose();
						//int R_det = R_mtx.determinant();
						// if ( R_det < 0){
						// 	R_mtx.col(2) = -1*R_mtx.col(2);
						// }
						::Eigen::Matrix<float,3,1> T_mtx = - R_mtx*cen_rot + cen_hot;
        				::Eigen::Matrix<float,3,4> Tran_mtx;					
				  		Tran_mtx.block<3,3>(0,0) = R_mtx;
        				Tran_mtx.block<3,1>(0,3) = T_mtx;

						EigenXform impose;
						impose.matrix() = Tran_mtx;
						EigenXform x_orig_position = EigenXform::Identity();
						
						impose = Xref * Xrotamer.inverse();


						// std::cout << "Xform:" << std::endl;
						// std::cout << Tran_mtx << std::endl;
						// ::Eigen::Matrix<float,3,3> temp;
						// temp.col(0)=Tran_mtx.block<3,1>(0,3);
						// temp.col(1)=Tran_mtx.block<3,1>(0,3);
						// temp.col(2)=Tran_mtx.block<3,1>(0,3);  		
			

					//std::cout << Tran_mtx.block<3,3>(0,0)*rif_res+temp << std::endl;
					//for( auto const & x_perturb : sample_position_deltas ){

					//std::cout << "being parallel block" << std::endl;
					#ifdef USE_OPENMP
					#pragma omp parallel for schedule(dynamic,16)
					#endif

					for(int a = 0; a < NSAMP; ++a){							
						EigenXform x_perturb;
						::scheme::numeric::rand_xform_sphere(rng,x_perturb,radius_bound,radians_bound);

						EigenXform x_position = x_2_orig_inverse * x_perturb  * x_2_orig * impose * x_orig_position;
						//EigenXform x_position = x_perturb * impose * x_orig_position;
						//EigenXform x_position = x_perturb * impose * x_orig_position;
							
						// you can check their "energies" against the target like this, obviously substituting the real rot# and position
						float positioned_rotamer_score = rot_tgt_scorer.score_rotamer_v_target( irot, x_position );
						//std::cout << positioned_rotamer_score << std::endl;
						// add the rotamer to the rif if it's any good
						std::ofstream myfile;

						//if( positioned_rotamer_score > 0) {positioned_rotamer_score = -1;}
						if( positioned_rotamer_score < 0){ // probably want this threshold to be an option or something
							//std::cout << positioned_rotamer_score << std::endl;
							accumulator->insert( x_position, positioned_rotamer_score-4, irot, i_hotspot_group, -1 );
						 // 	std::ostringstream os;
							// os << "hotspot_" << irot << "_" << a << ".pdb";
							// std::string s = os.str();		
							// myfile.open (s, std::fstream::in | std::fstream::out | std::fstream::app);
							
							// std::cout << "MODEL:" << irot << "_" << a << " " << positioned_rotamer_score << std::endl;
							// std::cout << positioned_rotamer_score << std::endl;

							// for( auto a : rotamer_atoms ){
							// 	a.set_position( x_position * a.position() );
							// 	::scheme::actor::write_pdb(myfile, a, params->rot_index_p->chem_index_ );
							// }

							// myfile << positioned_rotamer_score << std::endl;
							//myfile.close();
							// std::cout << "ENDMDL" << std::endl;

							//myfile << "TER" << std::endl;
							// myfile << positioned_rotamer_score << std::endl;
							//myfile.close();
							//std::cout << "ENDMDL" << std::endl;

							
						}
						myfile.close();
					} // end position perturbations

					//myfile.close();
					
				} // check residue type
				} // end loop over rotamers which match hotspot res
				

			} //  end loop over residues in hotspot group

		} // end loop over hotspot groups
		//utility_exit_with_message("done");
		// let the rif builder thing know you're done
		accumulator->checkpoint( std::cout );
		auto rif_ptr = accumulator->rif();
        std::cout << "testing rifbase key iteration" << std::endl;
        int count = 0;
        for( auto key : rif_ptr->key_range() ){
            EigenXform bin_center = rif_ptr->get_bin_center(key);
            // right... the edges of the bins.... this is only *mostly* true
            // runtime_assert( rif_ptr->get_bin_key(bin_center) == key );
            std::cout << "BIN: key " << key << " xform.trans: " << bin_center.translation().transpose() << std::endl;
            auto rotscores = rif_ptr->get_rotamers_for_key(key);
            runtime_assert( rotscores.size() > 0 );
            for( auto rot_score : rotscores ){
                // can't wait for cxx17 structured bindings!!!
                float rotamer_score = rot_score.first;
                int rotamer_number = rot_score.second;
                std::string resn = params->rot_index_p->resname(rotamer_number);
                std::cout << " rotamer " << rotamer_number << " " << resn << ", score " << rotamer_score << std::endl;
            }

        //    if(++count > 10) utility_exit_with_message("aireost");
        }
		
		//utility_exit_with_message("done");

	}


}
}
}


