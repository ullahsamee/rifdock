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

	#include <devel/init.hh>
	#include <riflib/RotamerGenerator.hh>
	#include <riflib/util.hh>

	#include <scheme/actor/Atom.hh>


namespace devel {
namespace scheme {
namespace rif {



	void
	RifGeneratorUserHotspots::generate_rif(
		RifAccumulatorP accumulator,
		RifGenParamsP params
	){
		typedef numeric::xyzVector<core::Real> Vec;


		// some sanity checks
		int const n_hspot_groups = opts.hotspot_files.size();
		runtime_assert_msg( n_hspot_groups, "no hotspot group files specified!!" );
		runtime_assert_msg( n_hspot_groups<16, "too many hotspot groups!!" );
		runtime_assert_msg( accumulator->rif()->has_sat_data_slots(), "This RIF type doesn't support sat groups!!!" );


		// setup the hacky but fast scorer
		devel::scheme::ScoreRotamerVsTarget<
				VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
			> rot_tgt_scorer;
		{
			std::vector< ::scheme::chemical::HBondRay > target_donors, target_acceptors;
			for( auto ir : params->target_res ){
				::devel::scheme::get_donor_rays   ( *params->target, ir, true, target_donors );
				::devel::scheme::get_acceptor_rays( *params->target, ir, true, target_acceptors );
			}
			std::cout << "target_donors.size() " << target_donors.size() << " target_acceptors.size() " << target_acceptors.size() << std::endl;
			{
				rot_tgt_scorer.rot_index_p_ = params->rot_index_p;
				rot_tgt_scorer.target_field_by_atype_ = params->field_by_atype;
				rot_tgt_scorer.target_donors_ = target_donors;
				rot_tgt_scorer.target_acceptors_ = target_acceptors;
				rot_tgt_scorer.hbond_weight_ = opts.hbond_weight;
				rot_tgt_scorer.upweight_multi_hbond_ = opts.upweight_multi_hbond;
				rot_tgt_scorer.upweight_iface_ = 1.0;

			}
		}

		std::vector<EigenXform> sample_position_deltas{ EigenXform::Identity() };
		// fill this in somehow... maybe start with random purterbations... later I can help you do it "right" with a grid of some kind

		// loop over files (one file is one hotspot group)
		for( int i_hotspot_group = 0; i_hotspot_group < opts.hotspot_files.size(); ++i_hotspot_group ){

			std::string const & hotspot_file = opts.hotspot_files[i_hotspot_group];
			core::pose::Pose pose;
			// read hotspot file into pose

			// read in pdb files # i_hotspot_group
			for( int i_hspot_res = 1; i_hspot_res <= pose.n_residue(); ++i_hspot_res ){

				// for each irot that is the right restype (can be had from rot_intex_p)
				int irot_begin = 0, irot_end=0;
				for( int irot = irot_begin; irot < irot_end; ++irot ){

					std::vector<SchemeAtom> const & rotamer_atoms( params->rot_index_p->atoms(irot) );


					// figure out the transform that aligns the rotamer's standard position onto your input hotspot
					// res to align to will be pose.residue(i_hspot_res)
					// this is a dummy
					EigenXform x_orig_position = EigenXform::Identity();

					for( auto const & x_perturb : sample_position_deltas ){

						EigenXform x_position = x_perturb * x_orig_position;

						// you can check their "energies" against the target like this, obviously substituting the real rot# and position
						float positioned_rotamer_score = rot_tgt_scorer.score_rotamer_v_target( irot, x_position );

						// add the rotamer to the rif if it's any good
						if( positioned_rotamer_score < 0 ){ // probably want this threshold to be an option or something
							accumulator->insert( x_position, positioned_rotamer_score, irot, i_hotspot_group, -1 );

							// probably you want to inspect what you're generating...
							// replace cout with a file stream, and you'll get a pdb file
							std::cout << "MODEL" << std::endl;
							for( auto a : rotamer_atoms ){
								a.set_position( x_position * a.position() );
								::scheme::actor::write_pdb( std::cout, a, params->rot_index_p->chem_index_ );
							}
							std::cout << "ENDMDL" << std::endl;

						}

					} // end position perturbations

				} // end loop over rotamers which match hotspot res

			} //  end loop over residues in hotspot group

		} // end loop over hotspot groups

		// let the rif builder thing know you're done
		accumulator->checkpoint( std::cout );
	}


}
}
}

