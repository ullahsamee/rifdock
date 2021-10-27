

#include <basic/options/option_macros.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <riflib/scaffold/nineA_util.hh>
#include <vector>
#include <utility/string_util.hh>

#ifdef GLOBAL_VARIABLES_ARE_BAD
	#ifndef INCLUDED_rif_dock_test_hh_1
	#define INCLUDED_rif_dock_test_hh_1   



OPT_1GRP_KEY(     StringVector , rif_dock, scaffolds )
	OPT_1GRP_KEY(  StringVector, rif_dock, scaffold_res )
	OPT_1GRP_KEY(  StringVector, rif_dock, scaffold_res_fixed )
    OPT_1GRP_KEY(  String      , rif_dock, scaffold_res_pdbinfo_labels )
	OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_res_use_best_guess )
    OPT_1GRP_KEY(  Boolean     , rif_dock, best_guess_mutate_to_val )
	OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_to_ala )
	OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_to_ala_selonly )
	OPT_1GRP_KEY(  Boolean     , rif_dock, replace_orig_scaffold_res )
	OPT_1GRP_KEY(  Boolean     , rif_dock, replace_all_with_ala_1bre )
	OPT_1GRP_KEY(  Boolean     , rif_dock, random_perturb_scaffold )
    OPT_1GRP_KEY(  StringVector, rif_dock, scaffold_clash_contexts )

	OPT_1GRP_KEY(  StringVector, rif_dock, target_bounding_xmaps )
	OPT_1GRP_KEY(  String      , rif_dock, target_pdb )
	OPT_1GRP_KEY(  String      , rif_dock, target_res )
	OPT_1GRP_KEY(  String      , rif_dock, target_rif )
	OPT_1GRP_KEY(  Real        , rif_dock, target_rf_resl )
	OPT_1GRP_KEY(  Integer     , rif_dock, target_rf_oversample )
	OPT_1GRP_KEY(  String      , rif_dock, target_rf_cache )
	OPT_1GRP_KEY(  String      , rif_dock, target_donors )
	OPT_1GRP_KEY(  String      , rif_dock, target_acceptors )
	OPT_1GRP_KEY(  Boolean     , rif_dock, only_load_highest_resl )
    OPT_1GRP_KEY(  Boolean     , rif_dock, dont_load_any_resl )
	OPT_1GRP_KEY(  Boolean     , rif_dock, use_rosetta_grid_energies )
	OPT_1GRP_KEY(  Boolean     , rif_dock, soft_rosetta_grid_energies )

	OPT_1GRP_KEY(  StringVector, rif_dock, data_cache_dir )

	OPT_1GRP_KEY(  Real        , rif_dock, beam_size_M )
    OPT_1GRP_KEY(  Real        , rif_dock, max_beam_multiplier )
    OPT_1GRP_KEY(  Boolean     , rif_dock, multiply_beam_by_seeding_positions )
    OPT_1GRP_KEY(  Boolean     , rif_dock, multiply_beam_by_scaffolds )
	OPT_1GRP_KEY(  Real        , rif_dock, search_diameter )
	OPT_1GRP_KEY(  Real        , rif_dock, hsearch_scale_factor )

	OPT_1GRP_KEY(  Real        , rif_dock, max_rf_bounding_ratio )
	OPT_1GRP_KEY(  Boolean     , rif_dock, make_bounding_plot_data )
	OPT_1GRP_KEY(  Boolean     , rif_dock, align_output_to_scaffold )
	OPT_1GRP_KEY(  Boolean     , rif_dock, output_scaffold_only )
	OPT_1GRP_KEY(  Boolean     , rif_dock, output_full_scaffold_only )
	OPT_1GRP_KEY(  Boolean     , rif_dock, output_full_scaffold )
    OPT_1GRP_KEY(  Boolean     , rif_dock, outputlite )
	OPT_1GRP_KEY(  Boolean     , rif_dock, parallelwrite )
	OPT_1GRP_KEY(  Boolean     , rif_dock, outputsilent )
	OPT_1GRP_KEY(  Integer     , rif_dock, n_pdb_out )
    OPT_1GRP_KEY(  Integer     , rif_dock, n_pdb_out_global )

	OPT_1GRP_KEY(  Real        , rif_dock, rf_resl )
	OPT_1GRP_KEY(  Integer     , rif_dock, rf_oversample )
	OPT_1GRP_KEY(  Boolean     , rif_dock, downscale_atr_by_hierarchy )
	OPT_1GRP_KEY(  Real        , rif_dock, favorable_1body_multiplier )
	OPT_1GRP_KEY(  Real        , rif_dock, favorable_1body_multiplier_cutoff )
	OPT_1GRP_KEY(  Real        , rif_dock, favorable_2body_multiplier )
    OPT_1GRP_KEY(  Real        , rif_dock, rotamer_onebody_inclusion_threshold )
    OPT_1GRP_KEY(  StringVector, rif_dock, rotamer_boltzmann_files )
    OPT_1GRP_KEY(  Boolean     , rif_dock, rotboltz_ignore_missing_rots )

	OPT_1GRP_KEY(  Integer     , rif_dock, rotrf_oversample )
	OPT_1GRP_KEY(  Real        , rif_dock, rotrf_resl )
	OPT_1GRP_KEY(  Real        , rif_dock, rotrf_spread )
	OPT_1GRP_KEY(  Real        , rif_dock, rotrf_scale_atr )
	OPT_1GRP_KEY(  String      , rif_dock, rotrf_cache_dir )

	OPT_1GRP_KEY(  Boolean     , rif_dock, hack_pack )
	OPT_1GRP_KEY(  Boolean     , rif_dock, hack_pack_during_hsearch )
	OPT_1GRP_KEY(  Real        , rif_dock, hack_pack_frac )
	OPT_1GRP_KEY(  Real        , rif_dock, pack_iter_mult )
	OPT_1GRP_KEY(  Integer     , rif_dock, pack_n_iters )
	OPT_1GRP_KEY(  Real       , rif_dock, hackpack_score_cut )
	OPT_1GRP_KEY(  Real        , rif_dock, hbond_weight )
    OPT_1GRP_KEY(  Real        , rif_dock, scaff_bb_hbond_weight )
    OPT_1GRP_KEY(  Boolean     , rif_dock, dump_scaff_bb_hbond_rays )
	OPT_1GRP_KEY(  Real        , rif_dock, upweight_multi_hbond )
	OPT_1GRP_KEY(  Real        , rif_dock, min_hb_quality_for_satisfaction )
	OPT_1GRP_KEY(  Real        , rif_dock, long_hbond_fudge_distance )
	OPT_1GRP_KEY(  Real        , rif_dock, global_score_cut )

	OPT_1GRP_KEY(  Real        , rif_dock, redundancy_filter_mag )
	OPT_1GRP_KEY(  Boolean     , rif_dock, filter_seeding_positions_separately )
	OPT_1GRP_KEY(  Boolean     , rif_dock, filter_scaffolds_separately )

	OPT_1GRP_KEY(  Real        , rif_dock, force_output_if_close_to_input )
	OPT_1GRP_KEY(  Integer     , rif_dock, force_output_if_close_to_input_num )

	OPT_1GRP_KEY(  Real        , rif_dock, upweight_iface )

	OPT_1GRP_KEY(  Boolean     , rif_dock, use_scaffold_bounding_grids )

	OPT_1GRP_KEY(  Boolean     , rif_dock, restrict_to_native_scaffold_res )
	OPT_1GRP_KEY(  Real        , rif_dock, bonus_to_native_scaffold_res )
	OPT_1GRP_KEY(  Boolean     , rif_dock, add_native_scaffold_rots_when_packing )
    OPT_1GRP_KEY(  Boolean     , rif_dock, native_docking )
    OPT_1GRP_KEY(  Real        , rif_dock, ignore_rifres_if_worse_than )

	OPT_1GRP_KEY(  Boolean     , rif_dock, dump_all_rif_rots )
	OPT_1GRP_KEY(  Boolean     , rif_dock, dump_all_rif_rots_into_output )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rif_rots_as_chains )
    OPT_1GRP_KEY(  Boolean     , rif_dock, dump_simple_atoms )
    OPT_1GRP_KEY(  Boolean     , rif_dock, ignore_ala_rifres )

    OPT_1GRP_KEY(  Boolean     , rif_dock, dump_rifgen_hdf5 )
	OPT_1GRP_KEY(  StringVector, rif_dock, dump_rifgen_near_pdb )
	OPT_1GRP_KEY(  Real        , rif_dock, dump_rifgen_near_pdb_dist )
	OPT_1GRP_KEY(  Real        , rif_dock, dump_rifgen_near_pdb_frac )
    OPT_1GRP_KEY(  Boolean     , rif_dock, dump_rifgen_near_pdb_last_atom )
	OPT_1GRP_KEY(  Boolean     , rif_dock, dump_rifgen_text )
    OPT_1GRP_KEY(  IntegerVector, rif_dock, dump_rifgen_for_sat )
    OPT_1GRP_KEY(  Integer     , rif_dock, dump_rifgen_for_sat_models )
    OPT_1GRP_KEY(  String      , rif_dock, dump_rifgen_for_sat_name3 )
    OPT_1GRP_KEY(  Integer     , rif_dock, dump_best_rifgen_rots )
    OPT_1GRP_KEY(  Real        , rif_dock, dump_best_rifgen_rmsd )
	OPT_1GRP_KEY(  String      , rif_dock, score_this_pdb )
	OPT_1GRP_KEY(  String      , rif_dock, dump_pdb_at_bin_center )
    OPT_1GRP_KEY(  Boolean     , rif_dock, test_hackpack )
    OPT_1GRP_KEY(  Boolean     , rif_dock, only_score_input_pos )

	OPT_1GRP_KEY(  String     , rif_dock, dokfile )
	OPT_1GRP_KEY(  String     , rif_dock, outdir )
	OPT_1GRP_KEY(  String     , rif_dock, output_tag )

	OPT_1GRP_KEY(  Boolean    , rif_dock, dont_use_scaffold_loops )
    OPT_1GRP_KEY(  Boolean    , rif_dock, dont_use_scaffold_helices )
    OPT_1GRP_KEY(  Boolean    , rif_dock, dont_use_scaffold_strands )

	OPT_1GRP_KEY(  Boolean    , rif_dock, dump_resfile )
	OPT_1GRP_KEY(  Boolean    , rif_dock, pdb_info_pikaa )
    OPT_1GRP_KEY(  Boolean    , rif_dock, pdb_info_pssm )

	OPT_1GRP_KEY(  Boolean    , rif_dock, cache_scaffold_data )

	OPT_1GRP_KEY(  Real        , rif_dock, tether_to_input_position )

	OPT_1GRP_KEY(  Boolean     , rif_dock, lowres_sterics_cbonly )

	OPT_1GRP_KEY(  Integer     , rif_dock, require_satisfaction )
	OPT_1GRP_KEY(  Integer     , rif_dock, num_hotspots )
	OPT_1GRP_KEY(  Integer     , rif_dock, require_n_rifres )
    OPT_1GRP_KEY(  Integer     , rif_dock, require_hydrophobic_residue_contacts )
    OPT_1GRP_KEY(  Real        , rif_dock, hydrophobic_ddg_cut )
    OPT_1GRP_KEY(  Real        , rif_dock, hydrophobic_ddg_weight )
    OPT_1GRP_KEY(  Real        , rif_dock, one_hydrophobic_better_than )
    OPT_1GRP_KEY(  Real        , rif_dock, two_hydrophobics_better_than )
    OPT_1GRP_KEY(  Real        , rif_dock, three_hydrophobics_better_than )
    OPT_1GRP_KEY(  Boolean     , rif_dock, better_than_must_hbond )
    OPT_1GRP_KEY(  Boolean     , rif_dock, count_all_contacts_as_hydrophobic )
    OPT_1GRP_KEY(  Real        , rif_dock, hydrophobic_ddg_per_atom_cut )
    OPT_1GRP_KEY(  String      , rif_dock, hydrophobic_target_res )
    OPT_1GRP_KEY(  StringVector, rif_dock, ligand_hydrophobic_res_atoms )
    OPT_1GRP_KEY(  Integer     , rif_dock, ligand_require_hydrophobic_residue_contacts )
    OPT_1GRP_KEY(  Real        , rif_dock, ligand_hydrophobic_ddg_weight )
    OPT_1GRP_KEY(  Integer     , rif_dock, num_cation_pi )
    OPT_1GRP_KEY(  Real        , rif_dock, CB_too_close_penalty )
    OPT_1GRP_KEY(  Real        , rif_dock, CB_too_close_dist )
    OPT_1GRP_KEY(  Real        , rif_dock, CB_too_close_resl )
    OPT_1GRP_KEY(  Integer     , rif_dock, CB_too_close_max_target_res_atom_idx )
    OPT_1GRP_KEY(  StringVector, rif_dock, specific_atoms_close_bonus )

	OPT_1GRP_KEY(  Boolean     , rif_dock, use_dl_mix_bb )

	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_fraction )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_then_min_below_thresh )
	OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_at_least )
	OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_at_most )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_min_fraction )
	OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_min_at_least )
    OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_min_at_most )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_fix_target )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_targetbb )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_scaffoldbb )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_allbb )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_cut )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_hard_min )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_score_total )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_score_ddg_only )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_rifres_rifres_weight )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_rifres_scaffold_weight )
	OPT_1GRP_KEY(  String      , rif_dock, rosetta_soft_score )
	OPT_1GRP_KEY(  String      , rif_dock, rosetta_hard_score )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_filter_before )
	OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_filter_n_per_scaffold )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_filter_redundancy_mag )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_filter_even_if_no_score )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_debug_dump_scores )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_score_select_random )
    OPT_1GRP_KEY(  Boolean     , rif_dock, skip_redundancy_filter_before_rosetta )
    OPT_1GRP_KEY(  Boolean     , rif_dock, override_rosetta_pose )

	OPT_1GRP_KEY(  Boolean     , rif_dock, extra_rotamers )
	OPT_1GRP_KEY(  Boolean     , rif_dock, extra_rif_rotamers )
	OPT_1GRP_KEY(  Integer     , rif_dock, always_available_rotamers_level )
    OPT_1GRP_KEY(  Boolean     , rif_dock, packing_use_rif_rotamers )
    OPT_1GRP_KEY(  Boolean     , rif_dock, dump_all_rifdock_rotamers )

    OPT_1GRP_KEY(  Integer     , rif_dock, nfold_symmetry )
    OPT_1GRP_KEY(  RealVector  , rif_dock, symmetry_axis )

    OPT_1GRP_KEY(  Real        , rif_dock, user_rotamer_bonus_constant )
    OPT_1GRP_KEY(  Real        , rif_dock, user_rotamer_bonus_per_chi )


    OPT_1GRP_KEY(  Real        , rif_dock, resl0 )

    OPT_1GRP_KEY(  Integer     , rif_dock, dump_x_frames_per_resl )
    OPT_1GRP_KEY(  Boolean     , rif_dock, dump_only_best_frames )
    OPT_1GRP_KEY(  Integer     , rif_dock, dump_only_best_stride )
    OPT_1GRP_KEY(  String      , rif_dock, dump_prefix )

    OPT_1GRP_KEY(  String      , rif_dock, scaff_search_mode )
    OPT_1GRP_KEY(  String      , rif_dock, nineA_cluster_path )
    OPT_1GRP_KEY(  String      , rif_dock, nineA_baseline_range )

    OPT_1GRP_KEY(  Integer     , rif_dock, low_cut_site )
    OPT_1GRP_KEY(  Integer     , rif_dock, high_cut_site )
    OPT_1GRP_KEY(  Integer     , rif_dock, max_insertion )
    OPT_1GRP_KEY(  Integer     , rif_dock, max_deletion )
    OPT_1GRP_KEY(  Real        , rif_dock, fragment_cluster_tolerance )
    OPT_1GRP_KEY(  Real        , rif_dock, fragment_max_rmsd )
    OPT_1GRP_KEY(  Integer     , rif_dock, max_fragments )
    OPT_1GRP_KEY(  StringVector, rif_dock, morph_rules_files )
    OPT_1GRP_KEY(  String      , rif_dock, morph_silent_file )
    OPT_1GRP_KEY(  String      , rif_dock, morph_silent_archetype )
    OPT_1GRP_KEY(  Real        , rif_dock, morph_silent_max_structures )
    OPT_1GRP_KEY(  Boolean     , rif_dock, morph_silent_random_selection )
    OPT_1GRP_KEY(  Real        , rif_dock, morph_silent_cluster_use_frac )

    OPT_1GRP_KEY(  Boolean     , rif_dock, include_parent )
    OPT_1GRP_KEY(  Boolean     , rif_dock, use_parent_body_energies )

    OPT_1GRP_KEY(  Integer     , rif_dock, dive_resl )
    OPT_1GRP_KEY(  Integer     , rif_dock, pop_resl )
    OPT_1GRP_KEY(  String      , rif_dock, match_this_pdb )
    OPT_1GRP_KEY(  Real        , rif_dock, match_this_rmsd )

    OPT_1GRP_KEY(  String      , rif_dock, rot_spec_fname )
    // constrain file
	OPT_1GRP_KEY(  StringVector, rif_dock, cst_files )

    OPT_1GRP_KEY(  Boolean     , rif_dock, write_seed_to_output )
	OPT_1GRP_KEY(  StringVector, rif_dock, seed_with_these_pdbs )
	OPT_1GRP_KEY(  Boolean     , rif_dock, seed_include_input )

	OPT_1GRP_KEY(  StringVector, rif_dock, seeding_pos )
    OPT_1GRP_KEY(  Real        , rif_dock, patchdock_min_sasa )
    OPT_1GRP_KEY(  Integer     , rif_dock, patchdock_top_ranks )
    OPT_1GRP_KEY(  Boolean     , rif_dock, seeding_by_patchdock )
    OPT_1GRP_KEY(  Boolean     , rif_dock, apply_seeding_xform_after_centering )
    OPT_1GRP_KEY(  String      , rif_dock, xform_pos )
    OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_each_seeding_at_least )
    OPT_1GRP_KEY(  Real        , rif_dock, cluster_score_cut )
    OPT_1GRP_KEY(  Real        , rif_dock, keep_top_clusters_frac )

    OPT_1GRP_KEY(  Boolean     , rif_dock, dump_xform_file )
    OPT_1GRP_KEY(  Real        , rif_dock, dump_override_cart_search_radius )
    OPT_1GRP_KEY(  Real        , rif_dock, dump_override_cart_search_resl )
    OPT_1GRP_KEY(  Real        , rif_dock, dump_override_angle_search_radius )
    OPT_1GRP_KEY(  Real        , rif_dock, dump_override_angle_search_resl )

    OPT_1GRP_KEY(  Real        , rif_dock, unsat_score_scalar )
    OPT_1GRP_KEY(  String      , rif_dock, unsat_helper )
    OPT_1GRP_KEY(  Boolean     , rif_dock, report_common_unsats )
    OPT_1GRP_KEY(  Real        , rif_dock, unsat_score_offset )
    OPT_1GRP_KEY(  Boolean     , rif_dock, unsat_debug )
    OPT_1GRP_KEY(  Boolean     , rif_dock, dump_presatisfied_donors_acceptors )

    OPT_1GRP_KEY(  Real        , rif_dock, burial_target_distance_cut )
    OPT_1GRP_KEY(  Real        , rif_dock, burial_target_neighbor_cut )
    OPT_1GRP_KEY(  Real        , rif_dock, burial_scaffold_distance_cut )
    OPT_1GRP_KEY(  Real        , rif_dock, burial_scaffold_neighbor_cut )
    OPT_1GRP_KEY(  Integer     , rif_dock, require_burial )

    OPT_1GRP_KEY(  Boolean     , rif_dock, force_calculate_sasa )
    OPT_1GRP_KEY(  Real        , rif_dock, sasa_cut )
    OPT_1GRP_KEY(  Real        , rif_dock, score_per_1000_sasa_cut )
    OPT_1GRP_KEY(  String      , rif_dock, skip_sasa_for_res )

    OPT_1GRP_KEY(  String      , rif_dock, buried_list )

    OPT_1GRP_KEY(  StringVector, rif_dock, pdbinfo_requirements )
    OPT_1GRP_KEY(  Integer     , rif_dock, num_pdbinfo_requirements_required )
    OPT_1GRP_KEY(  StringVector, rif_dock, requirement_groups )
    OPT_1GRP_KEY(  IntegerVector, rif_dock, requirements )
    OPT_1GRP_KEY(  String      ,  rif_dock, sat_score_bonus )
    OPT_1GRP_KEY(  String      ,  rif_dock, sat_score_override )

    OPT_1GRP_KEY(  StringVector,  rif_dock, pssm_file )
    OPT_1GRP_KEY(  Real        ,  rif_dock, pssm_weight )
    OPT_1GRP_KEY(  Real        ,  rif_dock, pssm_cutoff )
    OPT_1GRP_KEY(  Boolean     ,  rif_dock, pssm_higher_is_better )
    OPT_1GRP_KEY(  Boolean     ,  rif_dock, pssm_enforce_no_ala )

 

		void register_options() {
			using namespace basic::options;
			using namespace basic::options::OptionKeys;

			NEW_OPT(  rif_dock::scaffolds, "" , utility::vector1<std::string>() );
			NEW_OPT(  rif_dock::scaffold_res, "" , utility::vector1<std::string>() );
			NEW_OPT(  rif_dock::scaffold_res_fixed, "" , utility::vector1<std::string>() );
			NEW_OPT(  rif_dock::scaffold_res_pdbinfo_labels, "Use these comma-separate PDBInfo labels to select designable residues" , "" );
            NEW_OPT(  rif_dock::scaffold_res_use_best_guess, "" , false );
            NEW_OPT(  rif_dock::best_guess_mutate_to_val, "Mutate to polyvaline before making best guess", true );
			NEW_OPT(  rif_dock::scaffold_to_ala, "" , false );
			NEW_OPT(  rif_dock::scaffold_to_ala_selonly, "" , true );
			NEW_OPT(  rif_dock::replace_orig_scaffold_res, "", true );
			NEW_OPT(  rif_dock::replace_all_with_ala_1bre, "" , false );
			NEW_OPT(  rif_dock::random_perturb_scaffold, "" , false );
            NEW_OPT(  rif_dock::scaffold_clash_contexts, "One pdb per scaffold to be loaded as part of the scaffold only for clash checking. Imagine you have homodimers but are docking only one monomer.", utility::vector1<std::string>() );

			NEW_OPT(  rif_dock::target_bounding_xmaps, "" , utility::vector1<std::string>() );
			NEW_OPT(  rif_dock::target_pdb, "" , "" );
			NEW_OPT(  rif_dock::target_res, "" , "" );
			NEW_OPT(  rif_dock::target_rif, "" , "" );
			NEW_OPT(  rif_dock::target_rf_resl, ""       , 0.25 );
			NEW_OPT(  rif_dock::target_rf_oversample, "" , 2 );
			NEW_OPT(  rif_dock::downscale_atr_by_hierarchy, "" , true );
			NEW_OPT(  rif_dock::favorable_1body_multiplier, "Anything with a one-body energy less than favorable_1body_cutoff gets multiplied by this", 1 );
			NEW_OPT(  rif_dock::favorable_1body_multiplier_cutoff, "Anything with a one-body energy less than this gets multiplied by favorable_1body_multiplier", 0 );
			NEW_OPT(  rif_dock::favorable_2body_multiplier, "Anything with a two-body energy less than 0 gets multiplied by this", 1 );
            NEW_OPT(  rif_dock::rotamer_onebody_inclusion_threshold, "Threshold to include residue into 2-body calc. Increase this if 'crazy energy delta'", 8 );
            NEW_OPT(  rif_dock::rotamer_boltzmann_files, "Files that contain rotamer boltzmann penalties for each rifrot at select positions.", utility::vector1<std::string>() );
            NEW_OPT(  rif_dock::rotboltz_ignore_missing_rots, "Ignore mismatches in the number of rotamers. Missing rotamers get score of 0", false );

			NEW_OPT(  rif_dock::target_rf_cache, "" , "NO_CACHE_SPECIFIED_ON_COMMAND_LINE" );
			NEW_OPT(  rif_dock::target_donors, "", "" );
			NEW_OPT(  rif_dock::target_acceptors, "", "" );
			NEW_OPT(  rif_dock::only_load_highest_resl, "Only read in the highest resolution rif", false );
            NEW_OPT(  rif_dock::dont_load_any_resl, "This will certainly crash", false );
			NEW_OPT(  rif_dock::use_rosetta_grid_energies, "Use Frank's grid energies for scoring", false );
			NEW_OPT(  rif_dock::soft_rosetta_grid_energies, "Use soft option for grid energies", false );

			NEW_OPT(  rif_dock::data_cache_dir, "" , utility::vector1<std::string>(1,"./") );
			NEW_OPT(  rif_dock::beam_size_M, "" , 10.000000 );

			NEW_OPT(  rif_dock::max_beam_multiplier, "Maximum beam multiplier", 1 );
			NEW_OPT(  rif_dock::multiply_beam_by_seeding_positions, "Multiply beam size by number of seeding positions", false);
			NEW_OPT(  rif_dock::multiply_beam_by_scaffolds, "Multiply beam size by number of scaffolds", true);
			NEW_OPT(  rif_dock::max_rf_bounding_ratio, "" , 4 );
			NEW_OPT(  rif_dock::make_bounding_plot_data, "" , false );
			NEW_OPT(  rif_dock::align_output_to_scaffold, "" , false );
			NEW_OPT(  rif_dock::output_scaffold_only, "" , false );
			NEW_OPT(  rif_dock::output_full_scaffold_only, "" , false );
			NEW_OPT(  rif_dock::output_full_scaffold, "", false );
            NEW_OPT(  rif_dock::outputlite, "Write the output structures as compressed silent files", false );
	    NEW_OPT(  rif_dock::parallelwrite, "Write the output structures using all available threads", false );
			NEW_OPT(  rif_dock::outputsilent, "", false );
			NEW_OPT(  rif_dock::n_pdb_out, "" , 10 );
            NEW_OPT(  rif_dock::n_pdb_out_global, "Normally n_pdb_out applies to each seeding position, this caps the global", -1);

			NEW_OPT(  rif_dock::rf_resl, ""       , 0.25 );
			NEW_OPT(  rif_dock::rf_oversample, "" , 2 );

			NEW_OPT(  rif_dock::rotrf_oversample, "" , 2 );
			NEW_OPT(  rif_dock::rotrf_resl, "" , 0.3 );
			NEW_OPT(  rif_dock::rotrf_spread, "" , 0.0 );
			NEW_OPT(  rif_dock::rotrf_scale_atr, "" , 1.0 );
			NEW_OPT(  rif_dock::rotrf_cache_dir, "" , "./" );

			NEW_OPT(  rif_dock::hack_pack, "" , true );
			NEW_OPT(  rif_dock::hack_pack_during_hsearch, "hackpack during hsearch", false );
			NEW_OPT(  rif_dock::hack_pack_frac, "" , 0.2 );
			NEW_OPT(  rif_dock::pack_iter_mult, "" , 2.0 );
			NEW_OPT(  rif_dock::pack_n_iters, "" , 1 );
			NEW_OPT(  rif_dock::hackpack_score_cut, "", 0);
			NEW_OPT(  rif_dock::hbond_weight, "" , 2.0 );
            NEW_OPT(  rif_dock::scaff_bb_hbond_weight, "" , 0.0 );
            NEW_OPT(  rif_dock::dump_scaff_bb_hbond_rays, "Dump scaffold backbone hydrogen bond rays", false );
			NEW_OPT(  rif_dock::upweight_multi_hbond, "" , 0.0 );
			NEW_OPT(  rif_dock::min_hb_quality_for_satisfaction, "Minimum fraction of total hbond energy required for satisfaction. Scale -1 to 0", -0.6 );
			NEW_OPT(  rif_dock::long_hbond_fudge_distance, "Any hbond longer than 2A gets moved closer to 2A by this amount for scoring", 0.0 );
			NEW_OPT(  rif_dock::global_score_cut, "" , 0.0 );

			NEW_OPT(  rif_dock::redundancy_filter_mag, "" , 1.0 );
			NEW_OPT(  rif_dock::filter_seeding_positions_separately, "Redundancy filter each seeding position separately", true );
			NEW_OPT(  rif_dock::filter_scaffolds_separately, "Redundancy filter each scaffold separately", true );

			NEW_OPT(  rif_dock::force_output_if_close_to_input, "" , 1.0 );
			NEW_OPT(  rif_dock::force_output_if_close_to_input_num, "" , 0 );

			NEW_OPT(  rif_dock::upweight_iface, "", 1.2 );

			NEW_OPT(  rif_dock::use_scaffold_bounding_grids, "", false );

			NEW_OPT(  rif_dock::search_diameter, "", 150.0 );
			NEW_OPT(  rif_dock::hsearch_scale_factor, "global scaling of rotation/translation search grid", 1.0 );

			NEW_OPT(  rif_dock::restrict_to_native_scaffold_res, "aka structure prediction CHEAT. Depricated. Still allows ALA. Use -native_docking", false );
			NEW_OPT(  rif_dock::bonus_to_native_scaffold_res, "aka favor native CHEAT", -0.3 );
			NEW_OPT(  rif_dock::add_native_scaffold_rots_when_packing, "CHEAT. See -native_docking (which is separate)", false );
            NEW_OPT(  rif_dock::native_docking, "Best way to do docking with only native AAs. Added by bcov 2021.", false );
            NEW_OPT(  rif_dock::ignore_rifres_if_worse_than, "Don't use bad rif residues", 0 );

			NEW_OPT(  rif_dock::dump_all_rif_rots, "", false );
			NEW_OPT(  rif_dock::dump_all_rif_rots_into_output, "dump all rif rots into output", false);
			NEW_OPT(  rif_dock::rif_rots_as_chains, "dump rif rots as chains instead of models, loses resnum if true", false );
            NEW_OPT(  rif_dock::dump_simple_atoms, "Dump simple atom representation of scaffold for debugging.", false );
            NEW_OPT(  rif_dock::ignore_ala_rifres, "If ALA is the result of the rif search. Ignore it.", false );


            NEW_OPT(  rif_dock::dump_rifgen_hdf5, "Dump the rif to an hdf5 file.", false );
			NEW_OPT(  rif_dock::dump_rifgen_near_pdb, "dump rifgen rotamers with same AA type near this single residue", utility::vector1<std::string>());
			NEW_OPT(  rif_dock::dump_rifgen_near_pdb_dist, "", 1.0f );
			NEW_OPT(  rif_dock::dump_rifgen_near_pdb_frac, "", 1.0f );
            NEW_OPT(  rif_dock::dump_rifgen_near_pdb_last_atom, "Use only the last atom to decide if rotamers are close", false );
			NEW_OPT(  rif_dock::dump_rifgen_text, "Dump the rifgen tables within dump_rifgen_near_pdb_dist", false );
            NEW_OPT(  rif_dock::dump_rifgen_for_sat, "Dump a rotamers near a specific sat or multiple sats", utility::vector1<size_t>());
            NEW_OPT(  rif_dock::dump_rifgen_for_sat_models, "Number of rotamers to dump for dump_rifgen_for_sat", 200);
            NEW_OPT(  rif_dock::dump_rifgen_for_sat_name3, "Only dump rotamers with this name3", "");
            NEW_OPT(  rif_dock::dump_best_rifgen_rots, "Dump the best rotamer from the rifgen. This is how many", 0 );
            NEW_OPT(  rif_dock::dump_best_rifgen_rmsd, "The rmsd of the last atom for each aa type used in clustering", 3 );
			NEW_OPT(  rif_dock::score_this_pdb, "Score every residue of this pdb using the rif scoring machinery", "" );
			NEW_OPT(  rif_dock::dump_pdb_at_bin_center, "Dump each residue of this pdb at the rotamer's bin center", "" );
            NEW_OPT(  rif_dock::test_hackpack, "Test the packing objective in the original position too", false );
            NEW_OPT(  rif_dock::only_score_input_pos, "Dont' actually run the protocol, just score the input", false );

			NEW_OPT(  rif_dock::dokfile, "", "default.dok" );
			NEW_OPT(  rif_dock::outdir, "", "./" );
			NEW_OPT(  rif_dock::output_tag, "", "" );

			NEW_OPT(  rif_dock::dont_use_scaffold_loops, "", false );
            NEW_OPT(  rif_dock::dont_use_scaffold_helices, "", false );
            NEW_OPT(  rif_dock::dont_use_scaffold_strands, "", false );

			NEW_OPT(  rif_dock::dump_resfile, "", false );
			NEW_OPT(  rif_dock::pdb_info_pikaa, "", false );
            NEW_OPT(  rif_dock::pdb_info_pssm, "Put a pssm into the pdb info", false );

			NEW_OPT(  rif_dock::cache_scaffold_data, "", false );

			NEW_OPT(  rif_dock::tether_to_input_position, "", -1.0 );

			NEW_OPT(  rif_dock::lowres_sterics_cbonly, "", true );

			NEW_OPT(  rif_dock::require_satisfaction, "", 0 );
			NEW_OPT(  rif_dock::num_hotspots, "Number of hotspots found in Rifdock hotspots. If in doubt, set this to 1000", 0 );
			NEW_OPT(  rif_dock::require_n_rifres, "This sort of works during HackPack", 0 );
            NEW_OPT(  rif_dock::require_hydrophobic_residue_contacts, "How many target res to have at least 0.5 fa_sol, fa_atr, fa_rep with.", 0 );
            NEW_OPT(  rif_dock::hydrophobic_ddg_cut, "Really crappy approximation to hydrophobic ddg", 0 );
            NEW_OPT(  rif_dock::hydrophobic_ddg_weight, "Add the hydrophobic ddg to the rotamer score during hackpack with this weight.", 0 );
            NEW_OPT(  rif_dock::one_hydrophobic_better_than, "Require one rifres to have hydrophobic ddg better than this", 0 );
            NEW_OPT(  rif_dock::two_hydrophobics_better_than, "Require two rifres to have hydrophobic ddg better than this", 0 );
            NEW_OPT(  rif_dock::three_hydrophobics_better_than, "Require three rifres to have hydrophobic ddg better than this", 0 );
            NEW_OPT(  rif_dock::better_than_must_hbond, "*_hydrophobics_better_than must make a hbond to count", 0 );
            NEW_OPT(  rif_dock::count_all_contacts_as_hydrophobic, "This flag should only be used on very polar targets. It counts all contact types toward *_hydrophobics_better_than", 0 );
            NEW_OPT(  rif_dock::hydrophobic_ddg_per_atom_cut, "To be considered for better_than, must have ddg per atom better than this", 0 );
            NEW_OPT(  rif_dock::hydrophobic_target_res, "Comma separated list of residues to consider for hydrophobics. Default is all res", "" );
            NEW_OPT(  rif_dock::ligand_hydrophobic_res_atoms, "Zones of your ligand to be treated like residues for hydrophobic contacts. Format is SEQPOS:ATNAME,SEQPOS:ATNAME. Space separates multiple definitions.", utility::vector1<std::string>() );
            NEW_OPT(  rif_dock::ligand_require_hydrophobic_residue_contacts, "How many of the zones in -ligand_hydrophobic_res_atoms must an output satisfy to count?", 0 );
            NEW_OPT(  rif_dock::ligand_hydrophobic_ddg_weight, "Add the hydrophobic ddg of your ligand zones to the score with this weight during hackpack.", 0 );
            NEW_OPT(  rif_dock::num_cation_pi, "Number of cation pi's in output", 0 );
            NEW_OPT(  rif_dock::CB_too_close_penalty, "Penalty for every scaffold CB closer than -CB_too_close_dist to target.", 0 );
            NEW_OPT(  rif_dock::CB_too_close_dist, "Distance for -CB_too_close_penalty", 6 );
            NEW_OPT(  rif_dock::CB_too_close_resl, "Resolution of CB_too_close grid.", 0.5 );
            NEW_OPT(  rif_dock::CB_too_close_max_target_res_atom_idx, "Atom idx in Rosetta numbering of last heavyatom to use on target. CB is 5. Default is whole sidechain.", 1000 );
            NEW_OPT(  rif_dock::specific_atoms_close_bonus, "Bonus/penalty for specific scaffold atoms near specific target atoms. Use scaffold PDBInfo Labels. LABEL:atom_name."
                                                            " Then, format for space-separated bonuses is: LABEL,bonus,resl,close_dist,max_dist,resnum:atom_name[,resnum:atom_name,...] ."
                                                            " Bonus is linear between max_dist and close_dist then constant closer than close_dist." , utility::vector1<std::string>() );

            



			NEW_OPT(  rif_dock::use_dl_mix_bb, "use phi to decide where d is allow", false );

			NEW_OPT(  rif_dock::rosetta_score_fraction  , "",  0.00 );
			NEW_OPT(  rif_dock::rosetta_score_then_min_below_thresh, "", -9e9 );
			NEW_OPT(  rif_dock::rosetta_score_at_least, "", -1 );
			NEW_OPT(  rif_dock::rosetta_score_at_most, "", 999999999 );
			NEW_OPT(  rif_dock::rosetta_min_fraction  , "",  0.1 );
			NEW_OPT(  rif_dock::rosetta_min_at_least, "", 0 );
            NEW_OPT(  rif_dock::rosetta_min_at_most, "Min at most this many", -1 );
			NEW_OPT(  rif_dock::rosetta_min_targetbb  , "",  false );
			NEW_OPT(  rif_dock::rosetta_min_scaffoldbb  , "",  false );
			NEW_OPT(  rif_dock::rosetta_min_allbb  , "",  false );
			NEW_OPT(  rif_dock::rosetta_min_fix_target, "",  false );
			NEW_OPT(  rif_dock::rosetta_score_cut  , "", -10.0 );
			NEW_OPT(  rif_dock::rosetta_hard_min  , "", false );
			NEW_OPT(  rif_dock::rosetta_score_total  , "", false );
			NEW_OPT(  rif_dock::rosetta_score_ddg_only  , "", false );
			NEW_OPT(  rif_dock::rosetta_score_rifres_rifres_weight, "", 0.75 );
			NEW_OPT(  rif_dock::rosetta_score_rifres_scaffold_weight, "", 0.5 );
			NEW_OPT(  rif_dock::rosetta_soft_score, "", "beta_soft" );
			NEW_OPT(  rif_dock::rosetta_hard_score, "", "beta" );
			NEW_OPT(  rif_dock::rosetta_filter_before, "redundancy filter results before rosetta score", false );
			NEW_OPT(  rif_dock::rosetta_filter_n_per_scaffold, "use with rosetta_filter_before, num to save per scaffold", 300);
			NEW_OPT(  rif_dock::rosetta_filter_redundancy_mag, "use with rosetta_filter_before, redundancy mag on the clustering", 0.5);
			NEW_OPT(  rif_dock::rosetta_filter_even_if_no_score, "Do the filtering for rosetta score and min even if you don't actually score/min", false );
			NEW_OPT(  rif_dock::rosetta_debug_dump_scores, "dump lists of scores around the rosetta score and min", false);
			NEW_OPT(  rif_dock::rosetta_score_select_random, "Select random positions to score rather than best", false);
            NEW_OPT(  rif_dock::skip_redundancy_filter_before_rosetta, "For patchdock method, don't redundancy filter before rosetta score", false );
            NEW_OPT(  rif_dock::override_rosetta_pose, "Override the rosetta score and min output with the usual output", false );

			NEW_OPT(  rif_dock::extra_rotamers, "", true );
			NEW_OPT(  rif_dock::extra_rif_rotamers, "", true );
			NEW_OPT(  rif_dock::always_available_rotamers_level, "", 0 );
	        NEW_OPT(  rif_dock::packing_use_rif_rotamers, "", true );
            NEW_OPT(  rif_dock::dump_all_rifdock_rotamers, "Dump all the rotamers that rifdock reads from the rotamer spec file.", false );

	        NEW_OPT(  rif_dock::nfold_symmetry, "", 1 );
	        NEW_OPT(  rif_dock::symmetry_axis, "", utility::vector1<double>() );

	        NEW_OPT(  rif_dock::user_rotamer_bonus_constant, "", 0 );
			NEW_OPT(  rif_dock::user_rotamer_bonus_per_chi, "", 0 );

			NEW_OPT(  rif_dock::resl0, "", 16 );
			NEW_OPT(  rif_dock::dump_x_frames_per_resl, "Use this to make a movie", 0 );
			NEW_OPT(  rif_dock::dump_only_best_frames, "Only dump the best frames for the movie", false );
			NEW_OPT(  rif_dock::dump_only_best_stride, "When doing dump_only_best_frames, dump every Xth element of the best", 1 );
			NEW_OPT(  rif_dock::dump_prefix, "Convince Brian to make this autocreate the folder", "hsearch" );

			NEW_OPT(  rif_dock::scaff_search_mode, "Which scaffold mode and HSearch do you want? Options: default, morph_dive_pop, nineA_baseline", "default");
			NEW_OPT(  rif_dock::nineA_cluster_path, "Path to cluster database for nineA_baseline.", "" );
			NEW_OPT(  rif_dock::nineA_baseline_range, "format cdindex:low-high (python range style)", "");

			NEW_OPT(  rif_dock::low_cut_site, "The low cut point for fragment insertion, this res and the previous get minimized.", 0 );
			NEW_OPT(  rif_dock::high_cut_site, "The high cut point for fragment insertion, this res and the next get minimized.", 0 );
			NEW_OPT(  rif_dock::max_insertion, "Maximum number of residues to lengthen protein by.", 0 );
			NEW_OPT(  rif_dock::max_deletion, "Maximum number of residues to shorten protein by.", 0 );
			NEW_OPT(  rif_dock::fragment_cluster_tolerance, "RMSD cluster tolerance for fragments.", 0.5 );
			NEW_OPT(  rif_dock::fragment_max_rmsd , "Max RMSD to starting fragment.", 10000 );
			NEW_OPT(  rif_dock::max_fragments, "Maximum number of fragments to find.", 10000000 );
			NEW_OPT(  rif_dock::morph_rules_files, "List of files for each scaffold to specify morph regions", utility::vector1<std::string>() );
			NEW_OPT(  rif_dock::morph_silent_file, "Silent file containing pre-morphed structures. Overrides other options", "" );
			NEW_OPT(  rif_dock::morph_silent_archetype, "PDB to calculate transform difference between input position and silent file", "" );
			NEW_OPT(  rif_dock::morph_silent_max_structures, "Cluster silent file into this many cluster centers", 1000000000 );
			NEW_OPT(  rif_dock::morph_silent_random_selection, "Use random picks instead of clustering to limit silent file", false );
            NEW_OPT(  rif_dock::morph_silent_cluster_use_frac, "Cluster and take the top clusters that make up this frac of total", 1 );

			NEW_OPT(  rif_dock::include_parent, "Include parent fragment in diversified scaffolds.", false );
			NEW_OPT(  rif_dock::use_parent_body_energies, "Don't recalculate 1-/2-body energies for fragment insertions", false );

			NEW_OPT(  rif_dock::dive_resl , "Dive to this depth before diversifying", 5 );
			NEW_OPT(  rif_dock::pop_resl , "Return to this depth after diversifying", 4 );
			NEW_OPT(  rif_dock::match_this_pdb, "Like tether to input position but applied at diversification time.", "" );
			NEW_OPT(  rif_dock::match_this_rmsd, "RMSD for match_this_pdb", 7 );

			NEW_OPT(  rif_dock::rot_spec_fname,"rot_spec_fname","NOT SPECIFIED");
	        // constrain file names
			NEW_OPT(  rif_dock::cst_files, "" , utility::vector1<std::string>() );

            NEW_OPT(  rif_dock::write_seed_to_output, "Write the seeding position index to the scorefile", false );
			NEW_OPT(  rif_dock::seed_with_these_pdbs, "Use these pdbs as seeding positions, use this with tether_to_input_position", utility::vector1<std::string>() );
			NEW_OPT(  rif_dock::seed_include_input, "Include the input scaffold as a seeding position in seed_with_these_pdbs", true );

			NEW_OPT(  rif_dock::seeding_pos, "" , utility::vector1<std::string>() );
            NEW_OPT(  rif_dock::patchdock_min_sasa, "the cutoff sasa value for a valid patchdock output, the default is to use all of them", -1000.0);
            NEW_OPT(  rif_dock::patchdock_top_ranks, "only use the top solutions of patchdock to do refinement, the default is to use all of them", 99999);
            NEW_OPT(  rif_dock::seeding_by_patchdock, "The format of seeding file can be either Rosetta Xform or raw patchdock outputs", true );
            NEW_OPT(  rif_dock::apply_seeding_xform_after_centering, "Apply the seeding position xforms after moving scaffold to center", false );
            NEW_OPT(  rif_dock::xform_pos, "" , "" );
            NEW_OPT(  rif_dock::rosetta_score_each_seeding_at_least, "", -1 );
            NEW_OPT(  rif_dock::cluster_score_cut, "", 0);
            NEW_OPT(  rif_dock::keep_top_clusters_frac, "", 0.5);

            NEW_OPT(  rif_dock::dump_xform_file, "Dump a file to use for -xform_pos. See also dump_override_*", false );
            NEW_OPT(  rif_dock::dump_override_cart_search_radius, "The maximum cartesian distance sampled for the dumped xform_pos.", false );
            NEW_OPT(  rif_dock::dump_override_cart_search_resl, "The cartesian step size sampled for the dumped xform_pos.", false );
            NEW_OPT(  rif_dock::dump_override_angle_search_radius, "The maximum angle sampled for the dumped xform_pos in degrees.", false );
            NEW_OPT(  rif_dock::dump_override_angle_search_resl, "The angular step size sampled for the dumped_xform_pos in degrees.", false );

            NEW_OPT(  rif_dock::unsat_score_scalar, "The buried unsat weights get multiplied by this.", 0 );
            NEW_OPT(  rif_dock::unsat_helper, "Helper file for use with unsats", "" );
            NEW_OPT(  rif_dock::report_common_unsats, "Show probability of burying every unsat across all docks.", false );
            NEW_OPT(  rif_dock::unsat_score_offset, "This gets added to the score of all designs.", 0.0 );
            NEW_OPT(  rif_dock::unsat_debug, "Dump debug info from unsat calculations. Use with -test_hackpack", false );
            NEW_OPT(  rif_dock::dump_presatisfied_donors_acceptors, "Dump the presatisifed donors and acceptors", false );


            NEW_OPT(  rif_dock::burial_target_distance_cut, "Distance cutoff for target burial grid filling.", 4.0 );
            NEW_OPT(  rif_dock::burial_target_neighbor_cut, "Num neighbors to be buried in target burial grid.", 22 );
            NEW_OPT(  rif_dock::burial_scaffold_distance_cut, "Distance cutoff for scaffold burial grid filling.", 5.0 );
            NEW_OPT(  rif_dock::burial_scaffold_neighbor_cut, "Num neighbors to be scaffold in target burial grid.", 17 );
            NEW_OPT(  rif_dock::require_burial, "Require at least this many polar atoms be buried", 0 );

            NEW_OPT(  rif_dock::force_calculate_sasa, "Calculate Sasa even if it won't be used", false );
            NEW_OPT(  rif_dock::sasa_cut, "Anything with a sasa below this value is thrown out", 0 );
            NEW_OPT(  rif_dock::score_per_1000_sasa_cut, "Anything with a score per 1000 sasa units above this value is thrown out", 0 );
            NEW_OPT(  rif_dock::skip_sasa_for_res, "Comma separated list of residues to not include in sasa calculations. (like glycans)", "");

            NEW_OPT(  rif_dock::buried_list, "temp", "" );
            
            NEW_OPT(  rif_dock::pdbinfo_requirements, "Pairs of pdbinfo_label:req1,req2,req3 that specify that a residue with this pdbinfo_label must satisfy these requirements/sats." , utility::vector1<std::string>() );
            NEW_OPT(  rif_dock::num_pdbinfo_requirements_required, "Minimum number of pdbinfo_requirements to satisfy. -1 for all.", -1 );
            NEW_OPT(  rif_dock::requirement_groups, "I want at least 3 of these requirements: 3:1,5,8,10,23. I want less than 2 of these: -2:4,6,7. Space separated. Negative requirements means not this requirement.", utility::vector1<std::string>());
            NEW_OPT(  rif_dock::requirements,        "which rif residue should be in the final output", utility::vector1< int >());
            NEW_OPT(  rif_dock::sat_score_bonus,     "Give bonus to residues that satisfy sat. 0:-2,1:-1.5", "");
            NEW_OPT(  rif_dock::sat_score_override,  "Override score for residues that satisfy sat. 0:-2,1:-1.5", "");

            NEW_OPT(  rif_dock::pssm_file,  "Standard Rosetta formatted PSSM file (21 columns). Either 1 file for all scaffolds or 1 file per scaffold.", utility::vector1<std::string>());
            NEW_OPT(  rif_dock::pssm_weight,  "All values in the PSSM are multiplied by this before being considered scores. Don't forget -pssm_higher_is_better", 1);
            NEW_OPT(  rif_dock::pssm_cutoff,  "If specified, any PSSM value worse than this score denotes a disallowed residue type.", std::numeric_limits<double>::quiet_NaN() );
            NEW_OPT(  rif_dock::pssm_higher_is_better, "By default, PSSM values are interpretted as scores (lower=better). Set this to reverse that (and multiply everything by -1).", false);
            NEW_OPT(  rif_dock::pssm_enforce_no_ala, "If using -pssm_cutoff and ALA is disallowed. Add allowed rotamers during docking. May cause significant slowdown.", false);



		}
	#endif
#endif


#ifndef INCLUDED_rif_dock_test_hh_3
#define INCLUDED_rif_dock_test_hh_3   

struct RifDockOpt
{
	std::vector<std::string> scaffold_fnames;
	std::vector<std::string> scaffold_res_fnames;
	std::vector<std::string> data_cache_path;
	std::vector<std::string> rif_files;

	bool        VERBOSE                              ;
	double      resl0                                ;
	int64_t     DIM                                  ;
	int64_t     DIMPOW2                              ;
	int64_t     beam_size                            ;
    float       max_beam_multiplier                  ;
    bool        multiply_beam_by_seeding_positions   ;
    bool        multiply_beam_by_scaffolds           ;
	bool        replace_all_with_ala_1bre            ;
	bool        lowres_sterics_cbonly                ;
	float       tether_to_input_position_cut         ;
	bool        tether_to_input_position             ;
	float       global_score_cut                     ;
	std::string target_pdb                           ;
	std::string outdir                               ;
	std::string output_tag                           ;
	std::string dokfile_fname                        ;
	bool        dump_all_rif_rots                    ;
	bool        dump_all_rif_rots_into_output        ;
	bool        rif_rots_as_chains                   ;
    bool        dump_simple_atoms                    ;
    bool        ignore_ala_rifres                    ;
    bool        dump_rifgen_hdf5                     ;
	std::vector<std::string> dump_rifgen_near_pdb    ;
	float       dump_rifgen_near_pdb_dist            ;
	float       dump_rifgen_near_pdb_frac            ;
    bool        dump_rifgen_near_pdb_last_atom       ;
	bool        dump_rifgen_text                     ;
    std::vector<size_t> dump_rifgen_for_sat          ;
    int         dump_rifgen_for_sat_models           ;
    std::string dump_rifgen_for_sat_name3            ;
    int         dump_best_rifgen_rots                ;
    float       dump_best_rifgen_rmsd                ;
	std::string score_this_pdb                       ;
	std::string dump_pdb_at_bin_center               ;
    bool        test_hackpack                        ;  
    bool        only_score_input_pos                 ;
	bool        add_native_scaffold_rots_when_packing;
    bool        native_docking                       ;
    float       ignore_rifres_if_worse_than          ;
	bool        restrict_to_native_scaffold_res      ;
	float       bonus_to_native_scaffold_res         ;
	float       hack_pack_frac                       ;
	float       hsearch_scale_factor                 ;
	float       search_diameter                      ;
	bool        use_scaffold_bounding_grids          ;
    utility::vector1<std::string> scaffold_res_pdbinfo_labels;
	bool        scaffold_res_use_best_guess          ;
    bool        best_guess_mutate_to_val             ;
	bool        scaff2ala                            ;
	bool        scaff2alaselonly                     ;
	bool        replace_orig_scaffold_res            ;
	int         require_satisfaction                 ;
	int         num_hotspots                         ;
	int         require_n_rifres                     ;
    int         require_hydrophobic_residue_contacts ;
    float       hydrophobic_ddg_cut                  ;
    float       hydrophobic_ddg_weight               ;
    float       one_hydrophobic_better_than          ;
    float       two_hydrophobics_better_than         ;
    float       three_hydrophobics_better_than       ;
    bool        better_than_must_hbond               ;
    bool        count_all_contacts_as_hydrophobic    ;
    float       hydrophobic_ddg_per_atom_cut         ;
    utility::vector1<int> hydrophobic_target_res     ;
    std::vector<std::string> ligand_hydrophobic_res_atoms;
    int         ligand_require_hydrophobic_residue_contacts;
    float       ligand_hydrophobic_ddg_weight        ;
    int         num_cation_pi                        ;
    float       CB_too_close_penalty                 ;
    float       CB_too_close_dist                    ;
    float       CB_too_close_resl                    ;
    int         CB_too_close_max_target_res_atom_idx ;
    std::vector<std::string> specific_atoms_close_bonus;
	bool 		use_dl_mix_bb						 ;
	float       target_rf_resl                       ;
	bool        align_to_scaffold                    ;
	bool        output_scaffold_only                 ;
	bool        output_full_scaffold_only            ;
	bool        output_full_scaffold                 ;
    bool        outputlite                           ;
    bool 	parallelwrite				 ;	
	bool        outputsilent                         ;
	bool        pdb_info_pikaa                       ;
    bool        pdb_info_pssm                        ;
	bool        dump_resfile                         ;
	std::string target_res_fname                     ;
	int         target_rf_oversample                 ;
	float       max_rf_bounding_ratio                ;
	std::string target_rf_cache                      ;
	std::string target_donors                        ;
	std::string target_acceptors                     ;
	bool        only_load_highest_resl               ;
    bool        dont_load_any_resl                   ;
	bool        use_rosetta_grid_energies            ;
	bool        soft_rosetta_grid_energies           ;
	bool        downscale_atr_by_hierarchy           ;
	float       favorable_1body_multiplier           ;
	float       favorable_1body_multiplier_cutoff    ;
	float       favorable_2body_multiplier           ;
    float       rotamer_onebody_inclusion_threshold  ;
    std::vector<std::string> rotamer_boltzmann_fnames;
    bool        rotboltz_ignore_missing_rots         ;
	bool        random_perturb_scaffold              ;
    std::vector<std::string> scaffold_clash_contexts ;
	bool        dont_use_scaffold_loops              ;
    bool        dont_use_scaffold_helices            ;
    bool        dont_use_scaffold_strands            ;
	bool        cache_scaffold_data                  ;
	float       rf_resl                              ;
	bool        hack_pack                            ;
	bool        hack_pack_during_hsearch             ;
	int         rf_oversample                        ;

	int         rotrf_oversample                     ;
	float       rotrf_resl                           ;
	float       rotrf_spread                         ;
	std::string rotrf_cache_dir                      ;
	float       rotrf_scale_atr                      ;

	float       pack_iter_mult                       ;
	int         pack_n_iters                         ;
	float       hackpack_score_cut                   ;
	float       hbond_weight                         ;
    float       scaff_bb_hbond_weight                ;
    bool        dump_scaff_bb_hbond_rays             ;
	float       upweight_iface                       ;
	float       upweight_multi_hbond                 ;
	float       min_hb_quality_for_satisfaction      ;
	float       long_hbond_fudge_distance            ;
	float       redundancy_filter_mag                ;
	bool        filter_seeding_positions_separately  ;
	bool        filter_scaffolds_separately          ;
	int         force_output_if_close_to_input_num   ;
	float       force_output_if_close_to_input       ;
	int         n_pdb_out                            ;
    int         n_pdb_out_global                     ;
	bool        extra_rotamers                       ;
	bool        extra_rif_rotamers                   ;
	int         always_available_rotamers_level      ;
	int         packing_use_rif_rotamers             ;
    bool        dump_all_rifdock_rotamers            ;

	float       rosetta_score_fraction               ;
	float       rosetta_score_then_min_below_thresh  ;
	float       rosetta_score_at_least               ;
	float       rosetta_score_at_most                ;
	float       rosetta_min_fraction                 ;
	int         rosetta_min_at_least                 ;
    int         rosetta_min_at_most                  ;
	bool        rosetta_min_fix_target               ;
	bool        rosetta_min_targetbb                 ;
	bool        rosetta_min_scaffoldbb               ;
	bool        rosetta_min_allbb                    ;
	float       rosetta_score_cut                    ;
	float       rosetta_hard_min                     ;
	bool        rosetta_score_total                  ;
	bool        rosetta_score_ddg_only               ;
	float       rosetta_score_rifres_rifres_weight   ;
	float       rosetta_score_rifres_scaffold_weight ;
    bool        skip_redundancy_filter_before_rosetta;
    bool        override_rosetta_pose                ;

	bool        rosetta_beta                         ;
	std::string rosetta_soft_score                   ;
	std::string rosetta_hard_score                   ;
	bool        rosetta_filter_before                ;
	int         rosetta_filter_n_per_scaffold        ;
	float       rosetta_filter_redundancy_mag        ;
	bool        rosetta_filter_even_if_no_score      ;
	bool        rosetta_debug_dump_scores            ;
	bool        rosetta_score_select_random                ;

    int         nfold_symmetry                       ;
    std::vector<float> symmetry_axis                 ;

    float       user_rotamer_bonus_constant		     ;
    float       user_rotamer_bonus_per_chi		     ;

    int         dump_x_frames_per_resl				 ;
    bool        dump_only_best_frames                ;
    int         dump_only_best_stride                ;
    std::string dump_prefix                          ;

    std::string scaff_search_mode					 ;
    std::string nineA_cluster_path					 ;
    std::string nineA_baseline_range				 ;

    int         low_cut_site                         ;
    int         high_cut_site                        ;
    int         max_insertion                        ;
    int         max_deletion                         ;
    float       fragment_cluster_tolerance           ;
    float       fragment_max_rmsd                    ;
    int         max_fragments                        ;
    std::vector<std::string> morph_rules_fnames      ;
    std::string morph_silent_file                    ;
    std::string morph_silent_archetype               ;
    int         morph_silent_max_structures          ;
    bool        morph_silent_random_selection        ;
    float       morph_silent_cluster_use_frac        ;

    bool        include_parent                       ;
    bool        use_parent_body_energies             ;

    int         dive_resl                            ;
    int         pop_resl                             ;
    std::string match_this_pdb                       ;
    float       match_this_rmsd                      ;

    std::string rot_spec_fname                       ;
    // constrain file names
	std::vector<std::string> cst_fnames              ;

    bool        write_seed_to_output                  ;
	std::vector<std::string> seed_with_these_pdbs    ;
	bool        seed_include_input                   ;

    std::vector<std::string> seeding_fnames          ;
    std::string xform_fname                          ;
    float       rosetta_score_each_seeding_at_least  ;
    float       cluster_score_cut                    ;
    float       keep_top_clusters_frac               ;
    bool        seeding_by_patchdock                 ;
    bool        apply_seeding_xform_after_centering  ;
    float       patchdock_min_sasa                   ;
    int         patchdock_top_ranks                  ;

    bool        dump_xform_file                      ;
    double      dump_override_cart_search_radius     ;
    double      dump_override_cart_search_resl       ;
    double      dump_override_angle_search_radius    ;
    double      dump_override_angle_search_resl      ;

    float       unsat_score_scalar                   ;
    std::string unsat_helper                         ;
    bool        report_common_unsats                 ;
    float       unsat_score_offset                   ;
    bool        unsat_debug                          ;
    bool        dump_presatisfied_donors_acceptors   ;

    float       burial_target_distance_cut           ;
    float       burial_target_neighbor_cut           ;
    float       burial_scaffold_distance_cut         ;
    float       burial_scaffold_neighbor_cut         ;
    int         require_burial                       ;

    std::string buried_list                          ;
    
    std::vector<std::pair<std::string,std::vector<int>>> pdbinfo_requirements;
    int         num_pdbinfo_requirements_required    ;
    std::vector<std::pair<int,std::vector<int>>> requirement_groups;
    std::vector<int> requirements                    ;
    std::vector<float> sat_score_bonus               ;
    std::vector<bool> sat_score_override             ;

    bool        need_to_calculate_sasa               ;
    float       sasa_cut                             ;
    float       score_per_1000_sasa_cut              ;
    std::set<int> skip_sasa_for_res                  ;


    std::vector<std::string> pssm_file_fnames        ;
    float       pssm_weight                          ;
    float       pssm_cutoff                          ;
    bool        pssm_higher_is_better                ;
    bool        pssm_enforce_no_ala                  ;


    void init_from_cli();


};

#endif

#ifdef GLOBAL_VARIABLES_ARE_BAD

#ifndef INCLUDED_rif_dock_test_hh_4
#define INCLUDED_rif_dock_test_hh_4 

	void RifDockOpt::init_from_cli()
	{
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		runtime_assert( option[rif_dock::target_rif].user() );

		VERBOSE                                = false;
		resl0                                  = option[rif_dock::resl0                              ]();
		DIM                                    = 6;
		DIMPOW2                                = 1<<DIM;
		beam_size                              = int64_t( option[rif_dock::beam_size_M]() * 1000000.0 / DIMPOW2 ) * DIMPOW2;
        max_beam_multiplier                    = option[rif_dock::max_beam_multiplier                ]();
		multiply_beam_by_seeding_positions     = option[rif_dock::multiply_beam_by_seeding_positions ]();
		multiply_beam_by_scaffolds             = option[rif_dock::multiply_beam_by_scaffolds         ]();        
		replace_all_with_ala_1bre              = option[rif_dock::replace_all_with_ala_1bre          ]();

		target_pdb                             = option[rif_dock::target_pdb                         ]();
		lowres_sterics_cbonly                  = option[rif_dock::lowres_sterics_cbonly              ]();
		tether_to_input_position_cut           = option[rif_dock::tether_to_input_position           ]();
		tether_to_input_position               = tether_to_input_position_cut > 0.0;
		global_score_cut                       = option[rif_dock::global_score_cut                   ]();
		outdir                                 = option[rif_dock::outdir                             ]();
		output_tag                             = option[rif_dock::output_tag                         ]();
		dokfile_fname                          = outdir + "/" + option[rif_dock::dokfile             ]();
		dump_all_rif_rots                      = option[rif_dock::dump_all_rif_rots                  ]();
		dump_all_rif_rots_into_output		   = option[rif_dock::dump_all_rif_rots_into_output      ]();
		rif_rots_as_chains                     = option[rif_dock::rif_rots_as_chains                 ]();
        dump_simple_atoms                      = option[rif_dock::dump_simple_atoms                  ]();
        ignore_ala_rifres                      = option[rif_dock::ignore_ala_rifres                  ]();
        dump_rifgen_hdf5                       = option[rif_dock::dump_rifgen_hdf5                   ]();
		dump_rifgen_near_pdb_dist              = option[rif_dock::dump_rifgen_near_pdb_dist          ]();
		dump_rifgen_near_pdb_frac              = option[rif_dock::dump_rifgen_near_pdb_frac          ]();
        dump_rifgen_near_pdb_last_atom         = option[rif_dock::dump_rifgen_near_pdb_last_atom     ]();
		dump_rifgen_text                       = option[rif_dock::dump_rifgen_text                   ]();
        dump_rifgen_for_sat_models             = option[rif_dock::dump_rifgen_for_sat_models         ]();
        dump_rifgen_for_sat_name3              = option[rif_dock::dump_rifgen_for_sat_name3          ]();
        dump_best_rifgen_rots                  = option[rif_dock::dump_best_rifgen_rots              ]();
        dump_best_rifgen_rmsd                  = option[rif_dock::dump_best_rifgen_rmsd              ]();
		score_this_pdb                         = option[rif_dock::score_this_pdb                     ]();
		dump_pdb_at_bin_center                 = option[rif_dock::dump_pdb_at_bin_center             ]();
        test_hackpack                          = option[rif_dock::test_hackpack                      ]();  
        only_score_input_pos                   = option[rif_dock::only_score_input_pos               ]();
		add_native_scaffold_rots_when_packing  = option[rif_dock::add_native_scaffold_rots_when_packing ]();
        native_docking                         = option[rif_dock::native_docking                        ]();
        ignore_rifres_if_worse_than            = option[rif_dock::ignore_rifres_if_worse_than           ]();
		restrict_to_native_scaffold_res        = option[rif_dock::restrict_to_native_scaffold_res       ]();
		bonus_to_native_scaffold_res           = option[rif_dock::bonus_to_native_scaffold_res          ]();
		hack_pack_frac                         = option[rif_dock::hack_pack_frac                        ]();
		hsearch_scale_factor                   = option[rif_dock::hsearch_scale_factor                  ]();
		search_diameter                        = option[rif_dock::search_diameter                       ]();
		use_scaffold_bounding_grids            = option[rif_dock::use_scaffold_bounding_grids           ]();
		scaffold_res_use_best_guess            = option[rif_dock::scaffold_res_use_best_guess           ]();
        best_guess_mutate_to_val               = option[rif_dock::best_guess_mutate_to_val              ]();
		scaff2ala                              = option[rif_dock::scaffold_to_ala                       ]();
		scaff2alaselonly                       = option[rif_dock::scaffold_to_ala_selonly               ]();
		replace_orig_scaffold_res              = option[rif_dock::replace_orig_scaffold_res             ]();
		require_satisfaction                   = option[rif_dock::require_satisfaction                  ]();
		num_hotspots                           = option[rif_dock::num_hotspots                          ]();
		require_n_rifres                       = option[rif_dock::require_n_rifres                      ]();
        require_hydrophobic_residue_contacts   = option[rif_dock::require_hydrophobic_residue_contacts  ]();
        hydrophobic_ddg_cut                    = option[rif_dock::hydrophobic_ddg_cut                   ]();
        hydrophobic_ddg_weight                 = option[rif_dock::hydrophobic_ddg_weight                ]();
        one_hydrophobic_better_than            = option[rif_dock::one_hydrophobic_better_than           ]();
        two_hydrophobics_better_than           = option[rif_dock::two_hydrophobics_better_than          ]();
        three_hydrophobics_better_than         = option[rif_dock::three_hydrophobics_better_than        ]();
        better_than_must_hbond                 = option[rif_dock::better_than_must_hbond                ]();
        count_all_contacts_as_hydrophobic      = option[rif_dock::count_all_contacts_as_hydrophobic     ]();
        hydrophobic_ddg_per_atom_cut           = option[rif_dock::hydrophobic_ddg_per_atom_cut          ]();
        ligand_require_hydrophobic_residue_contacts = option[rif_dock::ligand_require_hydrophobic_residue_contacts]();
        ligand_hydrophobic_ddg_weight          = option[rif_dock::ligand_hydrophobic_ddg_weight         ]();
        num_cation_pi                          = option[rif_dock::num_cation_pi                         ]();
        CB_too_close_penalty                   = option[rif_dock::CB_too_close_penalty                  ]();
        CB_too_close_dist                      = option[rif_dock::CB_too_close_dist                     ]();
        CB_too_close_resl                      = option[rif_dock::CB_too_close_resl                     ]();
        CB_too_close_max_target_res_atom_idx   = option[rif_dock::CB_too_close_max_target_res_atom_idx  ]();
		use_dl_mix_bb						   = option[rif_dock::use_dl_mix_bb							]();
		target_rf_resl                         = option[rif_dock::target_rf_resl                        ]();
		align_to_scaffold                      = option[rif_dock::align_output_to_scaffold              ]();
		output_scaffold_only                   = option[rif_dock::output_scaffold_only                  ]();
		output_full_scaffold_only              = option[rif_dock::output_full_scaffold_only             ]();
		output_full_scaffold                   = option[rif_dock::output_full_scaffold                  ]();
        outputlite                                     = option[rif_dock::outputlite                            ]();
	parallelwrite                                  = option[rif_dock::parallelwrite                         ]();
		outputsilent                           = option[rif_dock::outputsilent                          ]();
		pdb_info_pikaa                         = option[rif_dock::pdb_info_pikaa                        ]();
        pdb_info_pssm                          = option[rif_dock::pdb_info_pssm                         ]();
		dump_resfile                           = option[rif_dock::dump_resfile                          ]();
		target_res_fname                       = option[rif_dock::target_res                            ]();
		target_rf_oversample                   = option[rif_dock::target_rf_oversample                  ]();
		max_rf_bounding_ratio                  = option[rif_dock::max_rf_bounding_ratio                 ]();
		target_rf_cache                        = option[rif_dock::target_rf_cache                       ]();
		target_donors                          = option[rif_dock::target_donors                         ]();
		target_acceptors                       = option[rif_dock::target_acceptors                      ]();		
		only_load_highest_resl                 = option[rif_dock::only_load_highest_resl                ]();
        dont_load_any_resl                     = option[rif_dock::dont_load_any_resl                    ]();
		use_rosetta_grid_energies              = option[rif_dock::use_rosetta_grid_energies             ]();
		soft_rosetta_grid_energies             = option[rif_dock::soft_rosetta_grid_energies            ]();
		downscale_atr_by_hierarchy             = option[rif_dock::downscale_atr_by_hierarchy            ]();
		favorable_1body_multiplier             = option[rif_dock::favorable_1body_multiplier            ]();
		favorable_1body_multiplier_cutoff      = option[rif_dock::favorable_1body_multiplier_cutoff     ]();
		favorable_2body_multiplier             = option[rif_dock::favorable_2body_multiplier            ]();
        rotamer_onebody_inclusion_threshold    = option[rif_dock::rotamer_onebody_inclusion_threshold   ]();
        rotboltz_ignore_missing_rots           = option[rif_dock::rotboltz_ignore_missing_rots          ]();
		random_perturb_scaffold                = option[rif_dock::random_perturb_scaffold               ]();
		dont_use_scaffold_loops                = option[rif_dock::dont_use_scaffold_loops               ]();
        dont_use_scaffold_helices              = option[rif_dock::dont_use_scaffold_helices             ]();
        dont_use_scaffold_strands              = option[rif_dock::dont_use_scaffold_strands             ]();
		cache_scaffold_data                    = option[rif_dock::cache_scaffold_data                   ]();
		rf_resl                                = option[rif_dock::rf_resl                               ]();
		hack_pack                              = option[rif_dock::hack_pack                             ]();
		hack_pack_during_hsearch               = option[rif_dock::hack_pack_during_hsearch              ]();

		rf_oversample                          = option[rif_dock::rf_oversample                         ]();
		redundancy_filter_mag                  = option[rif_dock::redundancy_filter_mag                 ]();
		filter_seeding_positions_separately    = option[rif_dock::filter_seeding_positions_separately   ]();
		filter_scaffolds_separately            = option[rif_dock::filter_scaffolds_separately           ]();
		rotrf_oversample                       = option[rif_dock::rotrf_oversample                      ]();
		rotrf_resl                             = option[rif_dock::rotrf_resl                            ]();
		rotrf_spread                           = option[rif_dock::rotrf_spread                          ]();
		rotrf_cache_dir                        = option[rif_dock::rotrf_cache_dir                       ]();
		rotrf_scale_atr                        = option[rif_dock::rotrf_scale_atr                       ]();
		pack_iter_mult                         = option[rif_dock::pack_iter_mult                        ]();
		pack_n_iters                           = option[rif_dock::pack_n_iters                          ]();
		hackpack_score_cut                     = option[rif_dock::hackpack_score_cut                    ]();
		hbond_weight                           = option[rif_dock::hbond_weight                          ]();
        scaff_bb_hbond_weight                  = option[rif_dock::scaff_bb_hbond_weight                 ]();
        dump_scaff_bb_hbond_rays               = option[rif_dock::dump_scaff_bb_hbond_rays              ]();
		upweight_iface                         = option[rif_dock::upweight_iface                        ]();
		upweight_multi_hbond                   = option[rif_dock::upweight_multi_hbond                  ]();
		min_hb_quality_for_satisfaction        = option[rif_dock::min_hb_quality_for_satisfaction       ]();
		long_hbond_fudge_distance              = option[rif_dock::long_hbond_fudge_distance             ]();
		redundancy_filter_mag                  = option[rif_dock::redundancy_filter_mag                 ]();
		force_output_if_close_to_input_num     = option[rif_dock::force_output_if_close_to_input_num    ]();
		force_output_if_close_to_input         = option[rif_dock::force_output_if_close_to_input        ]();
		n_pdb_out                              = option[rif_dock::n_pdb_out                             ]();
        n_pdb_out_global                       = option[rif_dock::n_pdb_out_global                      ]();
		extra_rotamers                         = option[rif_dock::extra_rotamers                        ]();
		extra_rif_rotamers                     = option[rif_dock::extra_rif_rotamers                    ]();
		always_available_rotamers_level        = option[rif_dock::always_available_rotamers_level       ]();
		packing_use_rif_rotamers               = option[rif_dock::packing_use_rif_rotamers              ]();
        dump_all_rifdock_rotamers              = option[rif_dock::dump_all_rifdock_rotamers             ]();

  		rosetta_score_fraction                 = option[rif_dock::rosetta_score_fraction                ]();
  		rosetta_score_then_min_below_thresh    = option[rif_dock::rosetta_score_then_min_below_thresh   ]();
  		rosetta_score_at_least                 = option[rif_dock::rosetta_score_at_least                ]();
  		rosetta_score_at_most                  = option[rif_dock::rosetta_score_at_most                 ]();
  		rosetta_min_fraction                   = option[rif_dock::rosetta_min_fraction                  ]();
  		rosetta_min_at_least                   = option[rif_dock::rosetta_min_at_least                  ]();
        rosetta_min_at_most                    = option[rif_dock::rosetta_min_at_most                   ]();
  		rosetta_min_fix_target                 = option[rif_dock::rosetta_min_fix_target                ]();
  		rosetta_min_targetbb                   = option[rif_dock::rosetta_min_targetbb                  ]();
  		rosetta_min_scaffoldbb                 = option[rif_dock::rosetta_min_scaffoldbb                ]();
  		rosetta_min_allbb                      = option[rif_dock::rosetta_min_allbb                     ]();
  		rosetta_score_cut                      = option[rif_dock::rosetta_score_cut                     ]();
  		rosetta_hard_min                       = option[rif_dock::rosetta_hard_min                      ]();
  		rosetta_score_total                    = option[rif_dock::rosetta_score_total                   ]();
  		rosetta_score_ddg_only                 = option[rif_dock::rosetta_score_ddg_only                ]();
  		rosetta_score_rifres_rifres_weight     = option[rif_dock::rosetta_score_rifres_rifres_weight    ]();
		rosetta_score_rifres_scaffold_weight   = option[rif_dock::rosetta_score_rifres_scaffold_weight  ]();
        skip_redundancy_filter_before_rosetta  = option[rif_dock::skip_redundancy_filter_before_rosetta ]();
        override_rosetta_pose                  = option[rif_dock::override_rosetta_pose                 ]();
		rosetta_soft_score                     = option[rif_dock::rosetta_soft_score  					]();
		rosetta_hard_score                     = option[rif_dock::rosetta_hard_score 				    ]();
		rosetta_beta                           = option[corrections::beta 								]();
		rosetta_filter_before                  = option[rif_dock::rosetta_filter_before                 ]();
		rosetta_filter_n_per_scaffold          = option[rif_dock::rosetta_filter_n_per_scaffold         ]();
		rosetta_filter_redundancy_mag          = option[rif_dock::rosetta_filter_redundancy_mag         ]();
		rosetta_filter_even_if_no_score        = option[rif_dock::rosetta_filter_even_if_no_score       ]();
		user_rotamer_bonus_constant 		   = option[rif_dock::user_rotamer_bonus_constant 			]();
		user_rotamer_bonus_per_chi 			   = option[rif_dock::user_rotamer_bonus_per_chi 			]();
		rosetta_debug_dump_scores              = option[rif_dock::rosetta_debug_dump_scores             ]();
		rosetta_score_select_random                  = option[rif_dock::rosetta_score_select_random                 ]();

		dump_x_frames_per_resl				   = option[rif_dock::dump_x_frames_per_resl                ]();
		dump_only_best_frames				   = option[rif_dock::dump_only_best_frames                 ]();
		dump_only_best_stride                  = option[rif_dock::dump_only_best_stride                 ]();
		dump_prefix                            = option[rif_dock::dump_prefix                           ]();

		scaff_search_mode					   = option[rif_dock::scaff_search_mode   				    ]();
		nineA_cluster_path					   = option[rif_dock::nineA_cluster_path                    ]();
		nineA_baseline_range				   = option[rif_dock::nineA_baseline_range                  ]();

		low_cut_site                           = option[rif_dock::low_cut_site                          ]();
		high_cut_site                          = option[rif_dock::high_cut_site                         ]();
		max_insertion                          = option[rif_dock::max_insertion                         ]();
		max_deletion                           = option[rif_dock::max_deletion                          ]();
		fragment_cluster_tolerance             = option[rif_dock::fragment_cluster_tolerance            ]();
        fragment_max_rmsd                      = option[rif_dock::fragment_max_rmsd                     ]();
        max_fragments                          = option[rif_dock::max_fragments                         ]();
        morph_silent_file                      = option[rif_dock::morph_silent_file                     ]();
        morph_silent_archetype                 = option[rif_dock::morph_silent_archetype                ]();
        morph_silent_max_structures            = option[rif_dock::morph_silent_max_structures           ]();
        morph_silent_random_selection          = option[rif_dock::morph_silent_random_selection         ]();
        morph_silent_cluster_use_frac          = option[rif_dock::morph_silent_cluster_use_frac         ]();

        include_parent                         = option[rif_dock::include_parent                        ]();
        use_parent_body_energies               = option[rif_dock::use_parent_body_energies              ]();

        dive_resl                              = option[rif_dock::dive_resl                             ]();
        pop_resl                               = option[rif_dock::pop_resl                              ]();
        match_this_pdb                         = option[rif_dock::match_this_pdb                        ]();
        match_this_rmsd                        = option[rif_dock::match_this_rmsd                       ]();

		rot_spec_fname						   = option[rif_dock::rot_spec_fname                        ]();

        write_seed_to_output                    = option[rif_dock::write_seed_to_output                   ]();
		seed_include_input                     = option[rif_dock::seed_include_input                    ]();

		seeding_by_patchdock                    = option[rif_dock::seeding_by_patchdock                 ]();
        apply_seeding_xform_after_centering     = option[rif_dock::apply_seeding_xform_after_centering  ]();
        xform_fname                             = option[rif_dock::xform_pos                            ]();
        rosetta_score_each_seeding_at_least     = option[rif_dock::rosetta_score_each_seeding_at_least  ]();
        cluster_score_cut                       = option[rif_dock::cluster_score_cut                    ]();
        keep_top_clusters_frac                  = option[rif_dock::keep_top_clusters_frac               ]();

        dump_xform_file                         = option[rif_dock::dump_xform_file                      ]();
        dump_override_cart_search_radius        = option[rif_dock::dump_override_cart_search_radius     ]();
        dump_override_cart_search_resl          = option[rif_dock::dump_override_cart_search_resl       ]();
        dump_override_angle_search_radius       = option[rif_dock::dump_override_angle_search_radius    ]();
        dump_override_angle_search_resl         = option[rif_dock::dump_override_angle_search_resl      ]();

        unsat_score_scalar                      = option[rif_dock::unsat_score_scalar                   ]();
        unsat_helper                            = option[rif_dock::unsat_helper                         ]();
        report_common_unsats                    = option[rif_dock::report_common_unsats                 ]();
        unsat_score_offset                      = option[rif_dock::unsat_score_offset                   ]();
        unsat_debug                             = option[rif_dock::unsat_debug                          ]();
        dump_presatisfied_donors_acceptors      = option[rif_dock::dump_presatisfied_donors_acceptors   ]();

        burial_target_distance_cut              = option[rif_dock::burial_target_distance_cut           ]();
        burial_target_neighbor_cut              = option[rif_dock::burial_target_neighbor_cut           ]();
        burial_scaffold_distance_cut            = option[rif_dock::burial_scaffold_distance_cut         ]();
        burial_scaffold_neighbor_cut            = option[rif_dock::burial_scaffold_neighbor_cut         ]();
        require_burial                          = option[rif_dock::require_burial                       ]();

        sasa_cut                                = option[rif_dock::sasa_cut                             ]();
        score_per_1000_sasa_cut                 = option[rif_dock::score_per_1000_sasa_cut              ]();

        buried_list                             = option[rif_dock::buried_list                          ]();
        
        num_pdbinfo_requirements_required       = option[rif_dock::num_pdbinfo_requirements_required    ]();


        pssm_weight                             = option[rif_dock::pssm_weight    ]();
        pssm_cutoff                             = option[rif_dock::pssm_cutoff    ]();
        pssm_higher_is_better                   = option[rif_dock::pssm_higher_is_better    ]();
        pssm_enforce_no_ala                     = option[rif_dock::pssm_enforce_no_ala      ]();



		for( std::string s : option[rif_dock::scaffolds     ]() )     scaffold_fnames.push_back(s);
		for( std::string s : option[rif_dock::scaffold_res  ]() ) scaffold_res_fnames.push_back(s);
		for( std::string s : option[rif_dock::data_cache_dir]() )     data_cache_path.push_back(s);

		for( std::string fn : option[rif_dock::target_bounding_xmaps]() ) rif_files.push_back(fn);
		rif_files.push_back( option[rif_dock::target_rif]() );

		if( scaff2ala && scaff2alaselonly &&  option[rif_dock::scaffold_to_ala_selonly].user() ){
			std::cout << "WARNING: -scaffold_to_ala overrides -scaffold_to_ala_selonly!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		}

		if( rosetta_score_total && rosetta_score_ddg_only ){
			std::cout << "WARNING: rosetta_score_total overrives rosetta_score_ddg_only" << std::endl;
			rosetta_score_ddg_only = false;
		}

		runtime_assert_msg( min_hb_quality_for_satisfaction < 0 && min_hb_quality_for_satisfaction > -1, 
			"-min_hb_quality_for_satisfaction must be between -1 and 0");

        runtime_assert_msg( !( outputsilent && outputlite ), "-outputsilent and -outputlite cannot be chosen together");

        nfold_symmetry = option[rif_dock::nfold_symmetry]();
        symmetry_axis.clear();
        if( option[rif_dock::symmetry_axis]().size() == 3 ){
            symmetry_axis.push_back( option[rif_dock::symmetry_axis]()[1] );
            symmetry_axis.push_back( option[rif_dock::symmetry_axis]()[2] );
            symmetry_axis.push_back( option[rif_dock::symmetry_axis]()[3] );
        } else if( option[rif_dock::symmetry_axis]().size() == 0 ){
            symmetry_axis.push_back(0);
            symmetry_axis.push_back(0);
            symmetry_axis.push_back(1);
        } else {
            std::cout << "bad rif_dock::symmetry_axis option" << std::endl;
            std::exit(-1);
        }


        need_to_calculate_sasa = option[rif_dock::force_calculate_sasa]() || sasa_cut > 0 || score_per_1000_sasa_cut < 0;

        for ( int res : utility::string_split( option[rif_dock::skip_sasa_for_res](), ',', int(0) ) ) {
            skip_sasa_for_res.insert(res);
        }


// Brian

        for ( size_t sat : option[rif_dock::dump_rifgen_for_sat]() ) {
            dump_rifgen_for_sat.push_back(sat);
        }

        if ( option[rif_dock::hydrophobic_target_res]().length() > 0) {
            int t = 0;
            hydrophobic_target_res = utility::string_split(option[rif_dock::hydrophobic_target_res](), ',', t);
        }


        if (option[rif_dock::use_scaffold_bounding_grids]()) {
        	std::cout << "ERROR: use_scaffold_bounding_grids no longer supported. Email bcov@uw.edu" << std::endl;
    		std::exit(-1);
        }


        if (option[rif_dock::nfold_symmetry]() > 1) {
        	std::cout << "ERROR: nfold_symmetry not currently supported. Email bcov@uw.edu" << std::endl;
    		std::exit(-1);
        }

        if ( option[rif_dock::native_docking]() || option[rif_dock::pssm_enforce_no_ala]() ) {
            rotamer_onebody_inclusion_threshold = 9e5; // Must allow everything but can't cross 9e9
        }


        if ( scaff_search_mode == "nineA_baseline" ) {
        	if ( scaffold_fnames.size() > 0 ) {
        		std::cout << "ERROR: can't use -scaffolds with nineA_baseline." << std::endl;
        		std::exit(-1);
        	}

        	std::vector<uint64_t> cdindex_then_clusts = devel::scheme::parse_nineA_baseline_range( nineA_baseline_range );
        	uint64_t num_scaffolds = cdindex_then_clusts.size() - 1;
        	runtime_assert( num_scaffolds > 0 );
        	scaffold_fnames.resize(num_scaffolds);

        }
        for( std::string s : option[rif_dock::scaffold_clash_contexts]() )     scaffold_clash_contexts.push_back(s);

        for( std::string s : option[rif_dock::morph_rules_files ]() ) morph_rules_fnames.push_back(s);

        // constrain file names
		for( std::string s : option[rif_dock::cst_files  ]() ) cst_fnames.push_back(s);

		for( std::string s : option[rif_dock::seed_with_these_pdbs ]() ) seed_with_these_pdbs.push_back(s);

        for( std::string s : option[rif_dock::seeding_pos ]() ) seeding_fnames.push_back(s);

        for( std::string s : option[rif_dock::rotamer_boltzmann_files ]() ) rotamer_boltzmann_fnames.push_back(s);

        for( std::string s : option[rif_dock::pssm_file]() )     pssm_file_fnames.push_back(s);

        for( std::string s : option[rif_dock::ligand_hydrophobic_res_atoms]() ) ligand_hydrophobic_res_atoms.push_back(s);

        for( std::string s : option[rif_dock::specific_atoms_close_bonus]() ) specific_atoms_close_bonus.push_back(s);


        patchdock_min_sasa                      = option[rif_dock::patchdock_min_sasa                  ]();
        patchdock_top_ranks                     = option[rif_dock::patchdock_top_ranks                 ]();
        
        for( int req : option[rif_dock::requirements]() ) requirements.push_back(req);

        for( std::string s : option[rif_dock::dump_rifgen_near_pdb]() ) dump_rifgen_near_pdb.push_back(s);

        if ( ! option[rif_dock::scaffold_res_pdbinfo_labels]().empty() ) {
            scaffold_res_pdbinfo_labels = utility::string_split(option[rif_dock::scaffold_res_pdbinfo_labels](), ',');
        }

        
        for( std::string s : option[rif_dock::pdbinfo_requirements]() ) {
            utility::vector1<std::string> pdbinfo_then_reqs = utility::string_split(s, ':');
            
            if ( pdbinfo_then_reqs.size() != 2 ) {
                std::cout << "ERROR: bad pdbinfo_requirement: " << s << std::endl;
                std::exit(-1);
            }
            utility::vector1<int> req_nos = utility::string_split<int>(pdbinfo_then_reqs[2], ',', int(0));
            std::vector<int> req_nos2;
            for ( int req : req_nos ) req_nos2.push_back( req );
            
            pdbinfo_requirements.push_back(std::pair<std::string,std::vector<int>>( pdbinfo_then_reqs[1], req_nos2 ));
	    }
    
        /////////   sat_score_bonus and sat_score_override   //////////////

        std::string bonus_string = option[rif_dock::sat_score_bonus]();
        std::string override_string = option[rif_dock::sat_score_override]();

        if ( bonus_string.length() > 0 || override_string.length() > 0 ) {

            std::vector<bool> override;
            std::vector<std::string> string;

            if ( bonus_string.length() > 0 ) {
                override.push_back(false);
                string.push_back(bonus_string);
            }
            if ( override_string.length() > 0 ) {
                override.push_back(true);
                string.push_back(override_string);
            }

            std::vector<bool> used;

            for ( int i = 0; i < override.size(); i++ ) {
                try {
                    for ( std::string const & pair : utility::string_split( string[i], ',' ) ) {
                        utility::vector1<std::string> sat_score = utility::string_split( pair, ':' );
                        float score = utility::from_string( sat_score[2], float(0) );
                        int low_sat = 0;
                        int high_sat = 0;
                        if ( sat_score[1].find("-") == std::string::npos ) {
                            low_sat = utility::from_string( sat_score[1], int(0) );
                            high_sat = low_sat;
                        } else {
                            utility::vector1<std::string> sats = utility::string_split( sat_score[1], '-' );
                            low_sat = utility::from_string( sats[1], int(0) );
                            high_sat = utility::from_string( sats[2], int(0) );
                        }

                        for ( int sat = low_sat; sat <= high_sat; sat++ ) {
                            if ( sat >= used.size() ) {
                                used.resize(sat+1);
                                sat_score_bonus.resize(sat+1);
                                sat_score_override.resize(sat+1);
                            }
                            runtime_assert( sat >= 0 );
                            if ( used[sat] ) {
                                utility_exit_with_message("Error, sat repeated twice in bonus/override: " + sat_score[1] );
                            }
                            used[sat] = true;
                            sat_score_bonus[sat] = score;
                            sat_score_override[sat] = override[i];
                        }
                    }
                } catch (...) {
                    utility_exit_with_message("Can't parse bonus/override string: " + string[i]);
                }
            }
        }

        /////////   requirement_groups   //////////////
        for( std::string s : option[rif_dock::requirement_groups]() ) {
            utility::vector1<std::string> num_then_reqs = utility::string_split(s, ':');
            
            if ( num_then_reqs.size() != 2 ) {
                std::cout << "ERROR: bad requirement_group: " << s << std::endl;
                std::exit(-1);
            }
            utility::vector1<int> req_nos = utility::string_split<int>(num_then_reqs[2], ',', int(0));
            std::vector<int> req_nos2;
            for ( int req : req_nos ) req_nos2.push_back( req );

            int num = utility::from_string( num_then_reqs[1], int(0) );
            
            requirement_groups.push_back(std::pair<int,std::vector<int>>( num, req_nos2 ));
        }
    
	}


#endif
#endif


