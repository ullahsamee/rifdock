// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// INC

#define GLOBAL_VARIABLES_ARE_BAD
#include <rif_dock_test.hh>
#undef GLOBAL_VARIABLES_ARE_BAD

	#include <numeric/random/random.hh>

	#include <ObjexxFCL/format.hh>

	#include <boost/foreach.hpp>
	#include <boost/lexical_cast.hpp>
	// #include <boost/random/mersenne_twister.hpp>

	// #include <core/id/AtomID.hh>
	#include <core/import_pose/import_pose.hh>
        #include <core/chemical/ChemicalManager.hh>
        #include <core/chemical/ResidueTypeSet.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pose/util.hh>
	// #include <core/scoring/EnergyGraph.hh>
	// #include <core/scoring/ScoreFunction.hh>
	// #include <core/scoring/ScoreFunctionFactory.hh>
	// #include <core/scoring/hbonds/HBondOptions.hh>
	// #include <core/scoring/methods/EnergyMethodOptions.hh>
	// #include <core/conformation/ResidueFactory.hh>
	// #include <core/kinematics/MoveMap.hh>
	// #include <core/scoring/Energies.hh>
 //    #include <protocols/minimization_packing/MinMover.hh>

	#include <devel/init.hh>
	// #include <riflib/RotamerGenerator.hh>
	#include <riflib/util.hh>

	// #include <numeric/alignment/QCP_Kernel.hh>
	#include <parallel/algorithm>
	#include <exception>
	#include <stdexcept>

	// #include <scheme/actor/Atom.hh>
	// #include <scheme/actor/BackboneActor.hh>
	// #include <scheme/actor/VoxelActor.hh>
	#include <scheme/kinematics/Director.hh>
	// #include <scheme/kinematics/SceneBase.hh>
	// #include <scheme/nest/pmap/OriTransMap.hh>
	// #include <scheme/numeric/rand_xform.hh>
	// // #include <scheme/objective/ObjectiveFunction.hh>
	// #include <scheme/objective/voxel/FieldCache.hh>
	// #include <scheme/objective/voxel/VoxelArray.hh>
	// #include <scheme/objective/hash/XformMap.hh>
	// #include <scheme/objective/storage/RotamerScores.hh>
	// #include <scheme/util/StoragePolicy.hh>
	// #include <scheme/search/HackPack.hh>
	// #include <scheme/objective/integration/SceneObjective.hh>

	#include <riflib/RifFactory.hh>

	#include <utility/file/file_sys_util.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>

	#include <chrono>
	#include <random>


/// Brian
	#include <scheme/objective/hash/XformHash.hh>
	#include <riflib/scaffold/ScaffoldDataCache.hh>
	#include <riflib/scaffold/ScaffoldProviderFactory.hh>
	#include <riflib/BurialManager.hh>
	#include <riflib/UnsatManager.hh>
	#include <riflib/ScoreRotamerVsTarget.hh>


// Task system
	#include <riflib/task/TaskProtocol.hh>
	#include <riflib/rifdock_tasks/HSearchTasks.hh>
	#include <riflib/rifdock_tasks/SetFaModeTasks.hh>
	#include <riflib/rifdock_tasks/HackPackTasks.hh>
	#include <riflib/rifdock_tasks/RosettaScoreAndMinTasks.hh>
	#include <riflib/rifdock_tasks/CompileAndFilterResultsTasks.hh>
	#include <riflib/rifdock_tasks/OutputResultsTasks.hh>
	#include <riflib/rifdock_tasks/UtilTasks.hh>
	#include <riflib/rifdock_tasks/SeedingPositionTasks.hh>
	#include <riflib/rifdock_tasks/MorphTasks.hh>

	#include <riflib/seeding_util.hh>




using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


int main(int argc, char *argv[]) {


	register_options();
	devel::init(argc,argv);


	devel::scheme::print_header( "setup global options" );
	RifDockOpt opt;
	opt.init_from_cli();
	utility::file::create_directory_recursive( opt.outdir );



	#ifdef USE_OPENMP
		omp_lock_t cout_lock, dump_lock;
		omp_init_lock( &cout_lock );
		omp_init_lock( &dump_lock );
	#endif


	using namespace core::scoring;
		using std::cout;
		using std::endl;
		using namespace devel::scheme;
		typedef numeric::xyzVector<core::Real> Vec;
		typedef numeric::xyzMatrix<core::Real> Mat;
		// typedef numeric::xyzTransform<core::Real> Xform;
		using ObjexxFCL::format::F;
		using ObjexxFCL::format::I;
		using devel::scheme::print_header;
		using ::devel::scheme::RotamerIndex;

	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////// static shit
	////////////////////////////////////////////////////////////////////////////////
	typedef ::scheme::util::SimpleArray<3,float> F3;
	typedef ::scheme::util::SimpleArray<3,int> I3;



		::scheme::search::HackPackOpts packopts;
		packopts.pack_n_iters         = opt.pack_n_iters;
		packopts.pack_iter_mult       = opt.pack_iter_mult;
		packopts.hbond_weight         = opt.hbond_weight;
		packopts.upweight_iface       = opt.upweight_iface;
		packopts.upweight_multi_hbond = opt.upweight_multi_hbond;
		packopts.use_extra_rotamers   = opt.extra_rotamers;
		packopts.always_available_rotamers_level = opt.always_available_rotamers_level;
		packopts.packing_use_rif_rotamers = opt.packing_use_rif_rotamers;
		packopts.add_native_scaffold_rots_when_packing = opt.add_native_scaffold_rots_when_packing;
		packopts.rotamer_inclusion_threshold = -0.5;//-0.5
		packopts.rotamer_onebody_inclusion_threshold = 5;//5
		packopts.init_with_best_1be_rots = true;
		packopts.user_rotamer_bonus_constant=opt.user_rotamer_bonus_constant;
		packopts.user_rotamer_bonus_per_chi=opt.user_rotamer_bonus_per_chi;

		std::string const rif_type = get_rif_type_from_file( opt.rif_files.back() );
		BOOST_FOREACH( std::string fn, opt.rif_files ){
			std::string rif_type2 = get_rif_type_from_file( fn );
			runtime_assert_msg( rif_type==rif_type2, "mismatched rif types, expect: " + rif_type + " got: " + rif_type2 + " for " + fn );
		}
		std::cout << "read RIF type: " << rif_type << std::endl;

		cout << "Search Resls: " << opt.resl0;
			std::vector<float> RESLS(1,opt.resl0);
			for( int i = 1; i < opt.rif_files.size(); ++i ){
				RESLS.push_back( RESLS.back()/2.0 );
				cout << " " << RESLS.back();
			}
			cout << endl;

		std::cout << "opt.rosetta_score_fraction: " << opt.rosetta_score_fraction << std::endl;
		std::cout << "opt.rosetta_score_then_min_below_thresh: " << opt.rosetta_score_then_min_below_thresh << std::endl;
		std::cout << "opt.rosetta_score_at_least: " << opt.rosetta_score_at_least << std::endl;
		std::cout << "opt.rosetta_score_at_most: " << opt.rosetta_score_at_most << std::endl;
		std::cout << "opt.rosetta_min_fraction: " << opt.rosetta_min_fraction << std::endl;
		std::cout << "opt.rosetta_min_targetbb: " << opt.rosetta_min_targetbb << std::endl;
		std::cout << "opt.rosetta_min_allbb: " << opt.rosetta_min_allbb << std::endl;
		std::cout << "opt.rosetta_score_cut: " << opt.rosetta_score_cut << std::endl;
		std::cout << "opt.require_satisfaction: " << opt.require_satisfaction << std::endl;

		std::cout << "//////////////////////////// end options /////////////////////////////////" << std::endl;



		// for( int iscaff = 0; iscaff < opt.scaffold_fnames.size(); ++iscaff )
		// {
		// 	std::string scaff_fname = opt.scaffold_fnames.at(iscaff);
		// 	std::cout << scaff_fname << std::endl;
		// 	core::pose::Pose scaffold;
		// 	utility::vector1<core::Size> scaffold_res;
		// 	core::import_pose::pose_from_file(scaffold, scaff_fname);
		// 	scaff_fname = utility::file::file_basename(utility::file_basename(scaff_fname));
		// 	scaffold.dump_pdb(scaff_fname+"_0.pdb");
		// 	scaffold_res = devel::scheme::get_designable_positions_best_guess( scaffold, opt.dont_use_scaffold_loops );
		// 	::devel::scheme::pose_to_ala( scaffold, scaffold_res );
		// 	scaffold.dump_pdb(scaff_fname+"_1.pdb");
		// }
		// utility_exit_with_message("test_scaff sel");


		////////////////////////////// should be no more use of options at this point! ///////////////////////////


		double time_rif=0, time_pck=0, time_ros=0;

		std::mt19937 rng( 0);//std::random_device{}() );


		{
			std::string dokfile_fname_orig = opt.dokfile_fname;
			int i = 2;
			while( utility::file::file_exists(opt.dokfile_fname) ){
				opt.dokfile_fname = dokfile_fname_orig + "." + str(i);
				++i;
			}
			if( i != 2)
				std::cout << "WARNING!" << dokfile_fname_orig << " already exists, using "
			              << opt.dokfile_fname << " instead!" << std::endl;
             else
             	std::cout << "output scores to " << opt.dokfile_fname << std::endl;
		}
		utility::io::ozstream dokout( opt.dokfile_fname );


		devel::scheme::RifFactoryConfig rif_factory_config;
		rif_factory_config.rif_type = rif_type;
		shared_ptr<RifFactory> rif_factory = ::devel::scheme::create_rif_factory( rif_factory_config );


	print_header( "create rotamer index" );
		
		std::cout << "Loading " << opt.rot_spec_fname << "..." << std::endl;
		std::string rot_index_spec_file = opt.rot_spec_fname;

		::scheme::chemical::RotamerIndexSpec rot_index_spec;
		shared_ptr< RotamerIndex > rot_index_p = ::devel::scheme::get_rotamer_index( rot_index_spec_file, opt.use_rosetta_grid_energies, rot_index_spec );
		RotamerIndex & rot_index( *rot_index_p );


		std::cout << "================ RotamerIndex ===================" << std::endl;\
		std::cout << rot_index << std::endl;
		std::cout << "=================================================" << std::endl;


		RotamerRFOpts rotrfopts;
		rotrfopts.oversample     = opt.rotrf_oversample;
		rotrfopts.field_resl     = opt.rotrf_resl;
		rotrfopts.field_spread   = opt.rotrf_spread;
		rotrfopts.data_dir       = opt.rotrf_cache_dir;
		rotrfopts.scale_atr      = opt.rotrf_scale_atr;
		::devel::scheme::RotamerRFTablesManager rotrf_table_manager( rot_index_p, rotrfopts );
		// rotrf_table_manager.preinit_all();
		

		MakeTwobodyOpts make2bopts;
		// hacked by brian             VVVV
		make2bopts.onebody_threshold = 4;
		make2bopts.distance_cut = 15.0;
		make2bopts.hbond_weight = packopts.hbond_weight;
		make2bopts.favorable_2body_multiplier = opt.favorable_2body_multiplier;




	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	print_header( "read and prepare target structure" ); //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	core::pose::Pose target;
	std::vector<SimpleAtom> target_simple_atoms;
	utility::vector1<core::Size> target_res;
	std::vector<HBondRay> target_donors, target_acceptors;
	float rif_radius=0.0, target_redundancy_filter_rg=0.0;
	shared_ptr<BurialManager> burial_manager;
	shared_ptr<UnsatManager> unsat_manager;
	bool donor_acceptors_from_file = false;
	{
		core::import_pose::pose_from_file( target, opt.target_pdb );

		if( opt.use_scaffold_bounding_grids ){
			for( int ir = 1; ir <= target.size(); ++ir ){
				utility::vector1<core::Size> resids(1,ir); // 1-index numbering
				std::vector<SchemeAtom> atoms;
				devel::scheme::get_scheme_atoms( target, resids, atoms );
				int restype = rot_index.chem_index_.resname2num( target.residue(ir).name3() ); // for UNK will be -1
				for( int ia = 0; ia < atoms.size(); ++ia){
					SchemeAtom const & a( atoms[ia] );
					runtime_assert( a.type() > 0 );
					if( a.type() >= 21 ) continue;
					SimpleAtom sa( a.position(), a.type(), restype, ia );
					target_simple_atoms.push_back(sa);
				}
			}
			std::cout << "target_simple_atoms.size() " << target_simple_atoms.size() << std::endl;
		}
		target_res = devel::scheme::get_res( opt.target_res_fname , target, /*nocgp*/false );
		get_rg_radius( target, target_redundancy_filter_rg, rif_radius, target_res, true ); // allatom for target
		rif_radius += 7.0; // hacky guess
		::devel::scheme::HBRayOpts hbopt;


		// hbopt.satisfied_atoms = ::devel::scheme::get_satisfied_atoms(target);

		// utility_exit_with_message("MAKE SURE LKBALL STUFF ISN'T FUCKING UP!!!");

		std::vector<std::pair<int,std::string>> donor_anames;
		std::vector<std::pair<int,std::string>> acceptor_anames;

		BOOST_FOREACH( core::Size ir, target_res ){
			::devel::scheme::get_donor_rays   ( target, ir, hbopt, target_donors, donor_anames );
			::devel::scheme::get_acceptor_rays( target, ir, hbopt, target_acceptors, acceptor_anames );
		}

		if ( opt.target_donors.length() > 0 || opt.target_acceptors.length() > 0 ) {
			donor_acceptors_from_file = true;

			if ( ! ( opt.target_donors.length() > 0 && opt.target_acceptors.length() > 0 ) ) {
				utility_exit_with_message("Must specify both -rifgen:target_donors AND -rifgen:target_acceptors!!!");
			}

			target_donors = load_hbond_rays( opt.target_donors );
			target_acceptors = load_hbond_rays( opt.target_acceptors );
			

		}

		std::cout << "target_donors.size() " << target_donors.size() << " target_acceptors.size() " << target_acceptors.size() << std::endl;

		if ( opt.unsat_orbital_penalty > 0 ) {
			if ( ! donor_acceptors_from_file ) {
				utility_exit_with_message("Must specify both -rifgen:target_donors and -rifgen:target_acceptors to use -unsat_orbital_penalty");
			}



			unsat_manager = make_shared<UnsatManager>( hbond::DefaultUnsatPenalties, rot_index_p, opt.unsat_debug );

			unsat_manager->set_target_donors_acceptors( target, target_donors, target_acceptors, donor_anames, acceptor_anames );
			unsat_manager->find_target_presatisfied( target );

			if (opt.dump_presatisfied_donors_acceptors) {
				unsat_manager->dump_presatisfied();
			}


			BurialOpts burial_opts;
			burial_opts.neighbor_distance_cutoff = opt.unsat_neighbor_cutoff;
			burial_opts.neighbor_count_weights.resize(0);
			for ( int i = 0; i < 20; i++ ) {
				float weight;
				if ( i < opt.unsat_neighbor_cutoff ) {
					weight = 0;
				} else {
					weight = opt.unsat_orbital_penalty;
				}
				burial_opts.neighbor_count_weights.push_back( weight );
			}
			std::cout << "burial weights: " << burial_opts.neighbor_count_weights << std::endl;
			burial_manager = make_shared<BurialManager>( burial_opts, unsat_manager->get_heavy_atom_xyzs() );
			burial_manager->set_target_neighbors( target );


		}
		// {
		// 	{
		// 		std::vector<HBondRay> tmpdon, tmpacc, tmpacclk;
		// 		::devel::scheme::HBRayOpts hbopt;
		// 		hbopt.lkball = false;
		// 		hbopt.withbb = true;
		// 		hbopt.add_acceptor_mid = true;
		// 		hbopt.satisfied_atoms = ::devel::scheme::get_satisfied_atoms(target);
		// 		BOOST_FOREACH( core::Size ir, target_res ){
		// 			::devel::scheme::get_donor_rays   ( target, ir, hbopt, tmpdon );
		// 			::devel::scheme::get_acceptor_rays( target, ir, hbopt, tmpacc );
		// 		}
		// 		hbopt.lkball = true;
		// 		BOOST_FOREACH( core::Size ir, target_res ){
		// 			::devel::scheme::get_acceptor_rays( target, ir, hbopt, tmpacclk );
		// 		}
		// 		target.dump_pdb("target.pdb");
		// 		utility::io::ozstream donout(utility::file_basename(opt.target_pdb)+"_donors.pdb");
		// 		::devel::scheme::dump_hbond_rays( donout, tmpdon, true );
		// 		donout.close();
		// 		utility::io::ozstream accout(utility::file_basename(opt.target_pdb)+"_orb_acceptors.pdb");
		// 		::devel::scheme::dump_hbond_rays( accout, tmpacc, false );
		// 		accout.close();
		// 		utility::io::ozstream accoutlk(utility::file_basename(opt.target_pdb)+"_lkb_acceptors.pdb");
		// 		::devel::scheme::dump_hbond_rays( accoutlk, tmpacclk, false );
		// 		accoutlk.close();
		// 		utility_exit_with_message("testing lkball replace orbs");
		// 	}
		// }
	}

	std::vector<bool> resl_load_map(RESLS.size(), true);	// by default load all resls

	if ( opt.only_load_highest_resl ) {
		for ( int i = 0; i < resl_load_map.size() - 1; i++) {
			resl_load_map[i] = false;
		}
	}

	std::vector< VoxelArrayPtr > target_field_by_atype;
	std::vector< std::vector< VoxelArrayPtr > > target_bounding_by_atype;
	{
		target_bounding_by_atype.resize( RESLS.size() );
		devel::scheme::RosettaFieldOptions rfopts;
		rfopts.field_resl = opt.target_rf_resl;
		rfopts.data_dir = "DUMMY_DATA_DIR_FIXME";
		rfopts.oversample = opt.target_rf_oversample;
		rfopts.block_hbond_sites = false;
		rfopts.max_bounding_ratio = opt.max_rf_bounding_ratio;
		rfopts.fail_if_no_cached_data = true;
		rfopts.repulsive_only_boundary = true;
		rfopts.cache_mismatch_tolerance = 0.01; // this is kinda loose...
		std::string cache_prefix = opt.target_rf_cache;
		devel::scheme::get_rosetta_fields_specified_cache_prefix(
			cache_prefix,
			opt.target_pdb,
			target,
			target_res,
			rfopts,
			target_field_by_atype,
			false
		);


		if( true ){
			// std::cout << "using target bounding grids, generating (or loading) them" << std::endl;
			devel::scheme::RosettaFieldOptions rfopts;
			rfopts.field_resl = opt.target_rf_resl;
			rfopts.data_dir = "DUMMY_DATA_DIR_FIXME";
			rfopts.oversample = opt.target_rf_oversample;
			rfopts.block_hbond_sites = false;
			rfopts.max_bounding_ratio = opt.max_rf_bounding_ratio;
			rfopts.fail_if_no_cached_data = true;
			rfopts.repulsive_only_boundary = true; // default
			devel::scheme::get_rosetta_bounding_fields_from_fba(
				RESLS,
				opt.target_pdb,
				target,
				target_res,
				rfopts,
				target_field_by_atype,
				target_bounding_by_atype,
				false,
				cache_prefix,
				resl_load_map
			);
			runtime_assert( target_bounding_by_atype.size() == RESLS.size() );
			// now scale down the any positive component by 1/RESL if RESL > 1
			if( opt.downscale_atr_by_hierarchy ){
				std::cout << "downscale_atr_by_hierarchy on target bounding steric grids" << std::endl;
				// std::cout << "  zeroing atr component of target bounding steric grids" << std::endl;
				for( int iresl = 0; iresl < RESLS.size(); ++iresl ){
					float correction = 1.0/RESLS[iresl];
					if( correction >= 1.0 ) break;
					BOOST_FOREACH( VoxelArrayPtr vap, target_bounding_by_atype[iresl] ){
						if( vap != nullptr ){
							std::exception_ptr exception = nullptr;
							#ifdef USE_OPENMP
							#pragma omp parallel for schedule(dynamic,64)
							#endif
							for( int k = 0; k < vap->num_elements(); ++k ){
								if( exception ) continue;
								try {
									float & dat = vap->data()[k];
									if( dat < 0 ){
										dat = dat * correction;
										// dat = 0; // for testing w/o attractive sterics
									}
								} catch( std::exception const & ex ) {
									#pragma omp critical
									exception = std::current_exception();
								}
							}
							if( exception ) std::rethrow_exception(exception);
						}
					}
				}
			}
		}
	}


#ifdef USEGRIDSCORE
	shared_ptr<protocols::ligand_docking::ga_ligand_dock::GridScorer> grid_scorer;
	if ( opt.use_rosetta_grid_energies ) {
		print_header( "preparing rosetta energy grids" );
		grid_scorer = prepare_grid_scorer( target, target_res );
	}
#else
	if ( opt.use_rosetta_grid_energies ) {
		utility_exit_with_message( "You must build with -DUSEGRIDSCORE=1 to use -use_rosetta_grid_energies!!!");
	}
#endif


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	print_header( "read in RIFs" ); /////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<shared_ptr<RifBase> > rif_ptrs;
	std::vector<bool> rif_using_rot;
	{
		std::vector<std::string> rif_descriptions( opt.rif_files.size() );
		rif_ptrs.resize( opt.rif_files.size() );
		std::exception_ptr exception = nullptr;
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for( int i_readmap = 0; i_readmap < opt.rif_files.size(); ++i_readmap ){
			if( exception ) continue;
			if ( ! resl_load_map.at(i_readmap) ) continue;
			try {
				std::string const & rif_file = opt.rif_files[i_readmap];
				std::string & rif_dscr = rif_descriptions[i_readmap];
				shared_ptr<RifBase> & rif_ptr = rif_ptrs[i_readmap];
				rif_ptr = rif_factory->create_rif_from_file( rif_file, rif_dscr );
				runtime_assert_msg( rif_ptrs[i_readmap] , "rif creation from file failed! " + rif_file );
				if( opt.VERBOSE ){
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					std::cout << "================= read " << rif_file << "=================" << std::endl
					          << "description:" << std::endl << rif_dscr << std::endl
					          << "load factor: " << rif_ptr->load_factor() << std::endl;
				}
				#ifdef USE_OPENMP
				#pragma omp critical
				#endif
				std::cout << "loaded RIF score for resl " << F(7,3,RESLS[i_readmap])
				          << " raw cart_resl: " << F(7,3,rif_ptr->cart_resl() )
				          << " raw ang_resl: " << F(7,3,rif_ptr->ang_resl() ) << std::endl;

				if (i_readmap == opt.rif_files.size() -1 ) {
					// #pragma omp criticial
					// rif_ptr->super_print( std::cout, rot_index_p );
					// std::ofstream out_file;
					// out_file.open("rif.txt");
					// rif_ptr->super_print( out_file, rot_index_p );
					// out_file.close();

				}
			} catch( std::exception const & ex ) {
				#ifdef USE_OPENMP
				#pragma omp critical
				#endif
				exception = std::current_exception();
			}
		}
		if( exception ) std::rethrow_exception(exception);

		std::cout << "RIF description:" << std::endl << rif_descriptions.back() << std::endl;
		std::cout << "load factor: " << rif_ptrs.back()->load_factor() << std::endl;
		std::cout << "size of value-type: " << rif_ptrs.back()->sizeof_value_type() << std::endl;
		std::cout << "mem_use: " << ::devel::scheme::KMGT( rif_ptrs.back()->mem_use() ) << std::endl;
		std::cout << "===================================================================================" << std::endl;

		rif_using_rot.resize( rot_index_p->size(), false );
		rif_using_rot[ rot_index.ala_rot() ] = true; // always include ala
		rif_ptrs.back()->get_rotamer_ids_in_use( rif_using_rot );
		int Nusingrot = 0;
		for( int i = 0; i < rif_using_rot.size(); ++i ){
			Nusingrot += rif_using_rot[i] ? 1 : 0;
		}
		std::cout << "rif uses: " << Nusingrot << " rotamers " << std::endl;


		if (opt.dump_rifgen_near_pdb.length() > 0) {
			float dump_dist = opt.dump_rifgen_near_pdb_dist;
			float dump_frac = opt.dump_rifgen_near_pdb_frac;
			core::pose::Pose pose = *(core::import_pose::pose_from_file(opt.dump_rifgen_near_pdb));
			core::conformation::Residue const & res = pose.residue(1);

			std::stringstream fname;
			fname << "rifgen_dump_" << opt.dump_rifgen_near_pdb << "_" << boost::str(boost::format("%.2f") % dump_dist ) << ".pdb.gz";

			rif_ptrs.back()->dump_rotamers_near_res( res, fname.str(), dump_dist, dump_frac, rot_index_p );
			if ( opt.dump_rifgen_text ) {
				rif_ptrs.back()->dump_rifgen_text_near_res( res, dump_dist, rot_index_p );
			}

			std::stringstream fname2;
			fname2 << "rifgen_dump_center_" << opt.dump_rifgen_near_pdb << ".pdb.gz";

			rif_ptrs.back()->dump_rotamers_at_bin_center( res, fname2.str(), rot_index_p );

		}

		if ( opt.score_this_pdb.length() > 0) {
			std::cout << "Loading " << opt.score_this_pdb << std::endl;
			core::pose::Pose pose = *(core::import_pose::pose_from_file(opt.score_this_pdb));

		    devel::scheme::ScoreRotamerVsTarget<
			        VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
			    > rot_tgt_scorer;
			    rot_tgt_scorer.rot_index_p_ = rot_index_p;
			    rot_tgt_scorer.target_field_by_atype_ = target_field_by_atype;
			    rot_tgt_scorer.target_donors_ = target_donors;
			    rot_tgt_scorer.target_acceptors_ = target_acceptors;
			    rot_tgt_scorer.hbond_weight_ = packopts.hbond_weight;
			    rot_tgt_scorer.upweight_iface_ = packopts.upweight_iface;
			    rot_tgt_scorer.upweight_multi_hbond_ = packopts.upweight_multi_hbond;
			#ifdef USEGRIDSCORE
			    rot_tgt_scorer.grid_scorer_ = grid_scorer;
			    rot_tgt_scorer.soft_grid_energies_ = opt.soft_rosetta_grid_energies;
			#endif

			std::cout << "Scoring " << opt.score_this_pdb << " against target" << std::endl;


			for ( int ires = 1; ires <= pose.size(); ires++ ) {


				core::conformation::Residue const & res = pose.residue(ires);
				int irot = rot_index_spec.get_matching_rot( res, 5.0f );

				if (irot == -1) {
					std::cout << "Res: " << ires << " No matching rotamer (within 5 deg for all chi)" << std::endl;
					continue;
				} 

				BBActor bb( res );

	            int sat1 = -1, sat2 = -1, hbcount = 0;
	            float score = rot_tgt_scorer.score_rotamer_v_target_sat( irot, bb.position(), sat1, sat2, true, hbcount, 10.0, 4 );

	            std::cout << "Res: " << ires
	                      << " score: " << score
	                      << " irot: " << irot
	                      << " nhbonds: " << hbcount
	                      << " sat1: " << sat1
	                      << " sat2: " << sat2
	                      << std::endl;
	        }

		}
	}


    if( 0 == opt.scaffold_fnames.size() ){
        std::cout << "WARNING: NO SCAFFOLDS!!!!!!" << std::endl;
    }

	for( int iscaff = 0; iscaff < opt.scaffold_fnames.size(); ++iscaff )
	{
		std::string scaff_fname = opt.scaffold_fnames.at(iscaff);
		std::vector<std::string> scaffold_sequence_glob0;				// Scaffold sequence in name3 space
		utility::vector1<core::Size> scaffold_res;//, scaffold_res_all; // Seqposs of residues to design, default whole scaffold
		try {

			ProtocolData pd;

			runtime_assert( rot_index_p );
			std::string scafftag = utility::file_basename( utility::file::file_basename( scaff_fname ) );

			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////   begin scaffold " << scafftag << " " << iscaff << " of " << opt.scaffold_fnames.size() << std::endl;
			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;


			bool needs_scaffold_director = false;

			ScaffoldProviderOP scaffold_provider = get_scaffold_provider(
				iscaff,
				rot_index_p,
				opt,
				make2bopts,
				rotrf_table_manager,
				needs_scaffold_director);

			// General info about a generic scaffold for debugging, cout, and the director
			ScaffoldDataCacheOP test_data_cache = scaffold_provider->get_data_cache_slow( ScaffoldIndex() );
			assert(test_data_cache);
			scaffold_sequence_glob0 = *(test_data_cache->scaffold_sequence_glob0_p);
			scaffold_res = *(test_data_cache->scaffold_res_p);
			float test_scaff_radius = test_data_cache->scaff_radius;
			float test_scaff_redundancy_filter_rg = test_data_cache->scaff_redundancy_filter_rg;
			Eigen::Vector3f test_scaffold_center = test_data_cache->scaffold_center;
			float test_redundancy_filter_rg = std::min( test_scaff_redundancy_filter_rg, target_redundancy_filter_rg );
			std::cout << "using redundancy_filter_rg: ~" << test_redundancy_filter_rg << std::endl;


			shared_ptr<std::vector<EigenXform>> seeding_positions = setup_seeding_positions( opt, pd, scaffold_provider, iscaff );


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "setup scene from scaffold and target" );
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			// SOMETHING WRONG, SCORES OFF BY A LITTLE
			// setup objectives, moved into scaffold loop to guarantee clean slate for each scaff...
			RifSceneObjectiveConfig rso_config;
				rso_config.packopts = &packopts;
				rso_config.rif_ptrs = rif_ptrs;
				rso_config.target_bounding_by_atype = &target_bounding_by_atype;
				rso_config.target_field_by_atype = &target_field_by_atype;
				rso_config.rot_index_p = rot_index_p;
				rso_config.target_donors = &target_donors;
				rso_config.target_acceptors = &target_acceptors;
				rso_config.require_satisfaction = opt.require_satisfaction;
				rso_config.require_n_rifres = opt.require_n_rifres;
#ifdef USEGRIDSCORE
            	rso_config.grid_scorer = grid_scorer;
            	rso_config.soft_grid_energies = opt.soft_rosetta_grid_energies;
#endif
            	rso_config.burial_manager = burial_manager;
            	rso_config.unsat_manager = unsat_manager;

            if ( opt.require_satisfaction > 0 && rif_ptrs.back()->has_sat_data_slots() ) {
            	if ( ! donor_acceptors_from_file && opt.num_hotspots == 0 ) {
            		utility_exit_with_message("New error message to fix an old bug!!! You can fix this error!!!"
            			"\n1. If you are using hotspots, you need to add this flag (and convince Brian/TaYi to fix this)"
            			"\n    -rif_dock:num_hotspots <number of hotspots>"
            			"\n   Feel free to overestimate. 1000 is pretty safe if in doubt."
            			"\n2. Otherwise you need to add these two flags"
            			"\n    -rif_dock:target_donors    <target donors file .pdb.gz>"
            			"\n    -rif_dock:target_acceptors <target acceptors file .pdb.gz>"
            			"\n   These files are already in your rifgen folder. Type ls <rifgen folder> *donor* to find them."
            			);
            	}

            	if ( opt.num_hotspots != 0 ) {
            		rso_config.n_sat_groups = opt.num_hotspots;
            	} else {
            		rso_config.n_sat_groups = target_donors.size() + target_acceptors.size();
            	}


            } else {
            	rso_config.n_sat_groups = 0;
            }
			

				

			ScenePtr scene_prototype;
			std::vector< ObjectivePtr > objectives;
			std::vector< ObjectivePtr > packing_objectives;
			runtime_assert( rif_factory->create_objectives( rso_config, objectives, packing_objectives ) );
			scene_prototype = rif_factory->create_scene();
			runtime_assert_msg( objectives.front()->is_compatible( *scene_prototype ), "objective and scene types not compatible!" );



			ScenePtr scene_minimal( scene_prototype->clone_deep() );
			scene_minimal->add_actor( 0, VoxelActor(target_bounding_by_atype) );


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "setup director based on scaffold and target sizes" ); //////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			shared_ptr<RifDockNestDirector> nest_director;
			DirectorBase director; {
				F3 target_center = pose_center(target);
				float body_radius = std::min( test_scaff_radius, rif_radius );

				double resl0 = opt.resl0;
				double hsearch_scale_factor = opt.hsearch_scale_factor;
				double search_diameter = opt.search_diameter;

				// Ideally one could read these in from the xform file
				if ( opt.xform_fname.length() > 0 ) {
					target_center = F3(0, 0, 0);
	                body_radius = 15.0;
					resl0 = 1;
	                hsearch_scale_factor = 1.2;
	                search_diameter = 4.0;
				}

				double cart_grid = resl0*hsearch_scale_factor/sqrt(3); // 1.5 is a big hack here.... 2 would be more "correct"
				double hackysin = std::min( 1.0, resl0*hsearch_scale_factor/2.0/ body_radius );

				runtime_assert( hackysin > 0.0 );
				double const rot_resl_deg0 = asin( hackysin ) * 180.0 / M_PI;
				int nside = std::ceil( search_diameter / cart_grid );
				std::cout << "search dia.    : " <<  search_diameter << std::endl;
				std::cout << "nside          : " << nside        << std::endl;
				std::cout << "resl0:           " << resl0 << std::endl;
				std::cout << "body_radius:     " << body_radius << std::endl;
				std::cout << "rif_radius:      " << rif_radius << std::endl;
				std::cout << "scaffold_radius: " << test_scaff_radius << std::endl;
				std::cout << "cart_grid:       " << cart_grid  << std::endl;
				std::cout << "rot_resl_deg0:   " << rot_resl_deg0 << std::endl;
				I3 nc( nside, nside, nside );
				F3 lb = target_center + F3( -cart_grid*nside/2.0, -cart_grid*nside/2.0, -cart_grid*nside/2.0 );
				F3 ub = target_center + F3(  cart_grid*nside/2.0,  cart_grid*nside/2.0,  cart_grid*nside/2.0 );
				std::cout << "cart grid ub " << ub << std::endl;
				std::cout << "cart grid lb " << lb << std::endl;
				std::cout << "(ub-lb/nc) = " << ((ub-lb)/nc.template cast<float>()) << std::endl;
				std::cout << "cartcen to corner (cart. covering radius): " << sqrt(3.0)*cart_grid/2.0 << std::endl;
				nest_director = make_shared<RifDockNestDirector>( rot_resl_deg0, lb, ub, nc, 1 );
				std::cout << "NestDirector:" << endl << *nest_director << endl;
				std::cout << "nest size0:    " << nest_director->size(0, RifDockIndex()).nest_index << std::endl;
				std::cout << "size of search space: ~" << float(nest_director->size(0, RifDockIndex()).nest_index)*1024.0*1024.0*1024.0 << " grid points" << std::endl;


				std::vector<DirectorBase> director_list;
				director_list.push_back( nest_director );  // Nest director must come first!!!!
				if ( needs_scaffold_director ) {
					director_list.push_back( make_shared<RifDockScaffoldDirector>(scaffold_provider, 1 ) );
				}
				if ( seeding_positions ) {
					director_list.push_back( make_shared<RifDockSeedingDirector>(seeding_positions, 1, -1 ) );
				}

				director = make_shared<RifDockDirector>(director_list);
			}



			std::vector< ScenePtr > scene_pt( omp_max_threads_1() );
			BOOST_FOREACH( ScenePtr & s, scene_pt ) s = scene_minimal->clone_deep();

			RifDockData rdd {
						iscaff,
						opt,
						RESLS,
						director,
						scene_pt,
						scene_minimal,
						target_simple_atoms,
						target_field_by_atype,
						&target_bounding_by_atype,
						&target_donors,
 						&target_acceptors,
 						target_redundancy_filter_rg,
 						target,
 						rot_index_p,
 						rotrf_table_manager,
 						objectives,
 						packing_objectives,
 						packopts,
 						rif_ptrs,
 						rso_config,
 						rif_factory,
 						nest_director->nest(),
    					#ifdef USE_OPENMP
 							dump_lock,
 						#endif
 						dokout,
 						scaffold_provider
#ifdef USEGRIDSCORE
    				,   grid_scorer
#endif
			};



			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "perform test with scaffold in original position" ); //////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			global_set_fa_mode( false, rdd );
    		test_data_cache->setup_onebody_tables( rot_index_p, opt);
			cout << std::endl;
			cout << "scores for scaffold in original position: " << std::endl;
			{
				EigenXform x(EigenXform::Identity());
				x.translation() = test_scaffold_center;
				director->set_scene( RifDockIndex(), 0, *scene_minimal);
				scene_minimal->set_position(1,x);
				for(int i = 0; i < RESLS.size(); ++i){
					if ( ! resl_load_map.at(i) ) continue;
					std::vector<float> sc = objectives[i]->scores(*scene_minimal);
					cout << "input bounding score " << i << " " << F(7,3,RESLS[i]) << " "
					     << F( 7, 3, sc[0]+sc[1] ) << " "
					     << F( 7, 3, sc[0]       ) << " "
					     << F( 7, 3, sc[1]       ) << endl;

				}
				if ( opt.test_hackpack ) {
					scaffold_provider->setup_twobody_tables( ScaffoldIndex() );
					if ( opt.unsat_orbital_penalty > 0 ) {
						scaffold_provider->setup_twobody_tables_per_thread( ScaffoldIndex() );
					}

					std::vector<float> sc = packing_objectives.back()->scores(*scene_minimal);
					cout << "input bounding score pack " << F(7,3,RESLS.back()) << " "
										     << F( 7, 3, sc[0]+sc[1] ) << " "
										     << F( 7, 3, sc[0]       ) << " "
										     << F( 7, 3, sc[1]       ) << endl;

				}


			}
			// std::cout << "scores for scaffold in original position: " << std::endl;
   //          {

   //  			test_data_cache->setup_twobody_tables( rot_index_p, opt, make2bopts, rotrf_table_manager);
   //              // EigenXform x(EigenXform::Identity());
   //              // x.translation() = test_scaffold_center;
   //              director->set_scene( RifDockIndex(4361221, 269, ScaffoldIndex()), 0, *scene_minimal);
   //              // scene_minimal->set_position(1,x);
   //              for(int i = 5; i < RESLS.size(); ++i){
   //                  std::vector<float> sc = packing_objectives.back()->scores(*scene_minimal);
   //                  std::cout << "input bounding score " << i << " " << F(7,3,RESLS[i]) << " "
   //                       << F( 7, 3, sc[0]+sc[1] ) << " "
   //                       << F( 7, 3, sc[0]       ) << " "
   //                       << F( 7, 3, sc[1]       ) << std::endl;

   //              }

   //                  devel::scheme::ScoreRotamerVsTarget<
   //      VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
   //  > rot_tgt_scorer;
   //  rot_tgt_scorer.rot_index_p_ = rot_index_p;
   //  rot_tgt_scorer.target_field_by_atype_ = target_field_by_atype;
   //  rot_tgt_scorer.target_donors_ = target_donors;
   //  rot_tgt_scorer.target_acceptors_ = target_acceptors;
   //  rot_tgt_scorer.hbond_weight_ = packopts.hbond_weight;
   //  rot_tgt_scorer.upweight_iface_ = packopts.upweight_iface;
   //  rot_tgt_scorer.upweight_multi_hbond_ = packopts.upweight_multi_hbond;
			// 	BBActor bb = scene_minimal->template get_actor<BBActor>(1,6);
			// 	float const recalc_rot_v_tgt = rot_tgt_scorer.score_rotamer_v_target( 277, bb.position(), 10.0, 4 );
			// 	std::cout << recalc_rot_v_tgt << std::endl;

   //          }


			int final_resl = rdd.RESLS.size() - 1;

			std::vector<shared_ptr<Task>> task_list;


			if ( opt.xform_fname.length() > 0) {
				create_rifine_task( task_list, rdd );
			} else {
				if ( opt.scaff_search_mode == "morph_dive_pop" ) {
					create_dive_pop_hsearch_task( task_list, rdd); 
				} else {


					task_list.push_back(make_shared<DiversifyBySeedingPositionsTask>()); // this is a no-op if there are no seeding positions
					task_list.push_back(make_shared<DiversifyByNestTask>( 0 ));

					task_list.push_back(make_shared<HSearchInit>( ));
					for ( int i = 0; i <= final_resl; i++ ) {
						task_list.push_back(make_shared<HSearchScoreAtReslTask>( i, i, opt.tether_to_input_position_cut ));

						if (opt.hack_pack_during_hsearch) {
							task_list.push_back(make_shared<SortByScoreTask>( ));
							task_list.push_back(make_shared<FilterForHackPackTask>( 1, rdd.packopts.pack_n_iters, rdd.packopts.pack_iter_mult ));
							task_list.push_back(make_shared<HackPackTask>( i, i, opt.global_score_cut )); 
						}

						task_list.push_back(make_shared<HSearchFilterSortTask>( i, opt.beam_size / opt.DIMPOW2, opt.global_score_cut, i < final_resl ));

						if (opt.dump_x_frames_per_resl > 0) {
							task_list.push_back(make_shared<DumpHSearchFramesTask>( i, i, opt.dump_x_frames_per_resl, opt.dump_only_best_frames, opt.dump_only_best_stride, 
								                                                    opt.dump_prefix + "_" + test_data_cache->scafftag + boost::str(boost::format("_resl%i")%i) ));
						}
						if ( i < final_resl ) {
							task_list.push_back(make_shared<HSearchScaleToReslTask>( i, i+1, opt.DIMPOW2, opt.global_score_cut )); 
						} 
					}
					task_list.push_back(make_shared<HSearchFinishTask>( opt.global_score_cut )); 
				}


				task_list.push_back(make_shared<SetFaModeTask>( true ));

				if ( opt.hack_pack ) {
					task_list.push_back(make_shared<FilterForHackPackTask>( opt.hack_pack_frac, rdd.packopts.pack_n_iters, rdd.packopts.pack_iter_mult ));
					task_list.push_back(make_shared<HackPackTask>(  final_resl, final_resl, opt.global_score_cut )); 
				}

				bool do_rosetta_score = opt.rosetta_score_fraction > 0 || opt.rosetta_score_then_min_below_thresh > -9e8 || opt.rosetta_score_at_least > 0;
				     do_rosetta_score = do_rosetta_score && opt.hack_pack;
				bool do_rosetta_min   = rdd.opt.rosetta_min_fraction > 0.0 && do_rosetta_score;

				if ( do_rosetta_score ) {
					if (opt.rosetta_filter_before) {
						task_list.push_back(make_shared<CompileAndFilterResultsTask>( final_resl, final_resl, opt.rosetta_filter_n_per_scaffold, opt.rosetta_filter_redundancy_mag, 
																				      0, 0, opt.filter_seeding_positions_separately, opt.filter_scaffolds_separately )); 
					} 
					else {
						task_list.push_back(make_shared<FilterForRosettaScoreTask>( opt.rosetta_score_fraction,  opt.rosetta_score_then_min_below_thresh, opt.rosetta_score_at_least, 
							                                                        opt.rosetta_score_at_most, opt.rosetta_score_select_random )); 
					}

					if (opt.rosetta_debug_dump_scores) task_list.push_back(make_shared<DumpScoresTask>( "hackpack_scores.dat")); 
					if (opt.rosetta_debug_dump_scores) task_list.push_back(make_shared<DumpRotScoresTask>( "hackpack_rot_scores.dat", false, final_resl)); 

					task_list.push_back(make_shared<RosettaScoreTask>( final_resl, opt.rosetta_score_cut, do_rosetta_min, !do_rosetta_min)); 

					if (opt.rosetta_debug_dump_scores) task_list.push_back(make_shared<DumpScoresTask>( "rosetta_scores.dat")); 
				}

				if ( do_rosetta_min ) {
					task_list.push_back(make_shared<FilterForRosettaMinTask>( opt.rosetta_min_fraction, opt.rosetta_min_at_least ));
					task_list.push_back(make_shared<RosettaMinTask>( final_resl, opt.rosetta_score_cut, true )); 

					if (opt.rosetta_debug_dump_scores) task_list.push_back(make_shared<DumpScoresTask>( "rosetta_min_scores.dat")); 
				}
				
				task_list.push_back(make_shared<CompileAndFilterResultsTask>( final_resl, final_resl, opt.n_pdb_out, opt.redundancy_filter_mag, opt.force_output_if_close_to_input_num, 
					                                                          opt.force_output_if_close_to_input, opt.filter_seeding_positions_separately, 
					                                                          opt.filter_scaffolds_separately ));
				task_list.push_back(make_shared<OutputResultsTask>( final_resl, final_resl));
			}


			TaskProtocol protocol( task_list );


			shared_ptr<std::vector<SearchPoint>> starting_point = make_shared<std::vector<SearchPoint>>( );
			starting_point->push_back(SearchPoint(RifDockIndex()));

			ThreePointVectors input;
			input.search_points = starting_point;
			ThreePointVectors results = protocol.run( input, rdd, pd );

			time_rif += pd.time_rif;
			time_pck += pd.time_pck;
			time_ros += pd.time_ros;





		} catch( std::exception const & ex ) {
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "error (below) on scaffold " << scaff_fname << " (will continue with others, if any)" << std::endl;
			std::cout << ex.what() << std::endl;
			std::cout << "scene residue numering (may help debug):" << std::endl;
			for( int i = 1; i <= scaffold_res.size(); ++i ){
				std::cout << "scene res numbering: " << i-1 << " " << scaffold_sequence_glob0.at(scaffold_res[i]-1) << " pose number: " << scaffold_res[i] << std::endl;
			}
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		} catch ( ... ) {
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "unknown error on scaffold " << scaff_fname << ", will continue with others, if any." << std::endl;
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		}


	} // end scaffold loop


	dokout.close();






	#ifdef USE_OPENMP
		omp_destroy_lock( &cout_lock );
		omp_destroy_lock( &dump_lock );
	#endif

	std::cout << "rif_dock_test_DONE" << std::endl;

	return 0;
 }
