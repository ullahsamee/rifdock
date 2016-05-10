// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// INC
	#include <numeric/random/random.hh>

	#include <ObjexxFCL/format.hh>

	#include <basic/options/option_macros.hh>

	#include <boost/foreach.hpp>
	#include <boost/lexical_cast.hpp>
	// #include <boost/random/mersenne_twister.hpp>

	#include <core/id/AtomID.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pose/util.hh>
	#include <core/scoring/ScoreFunction.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	#include <core/scoring/hbonds/HBondOptions.hh>
	#include <core/scoring/methods/EnergyMethodOptions.hh>
	#include <core/conformation/ResidueFactory.hh>
	#include <protocols/simple_moves/MinMover.hh>
	#include <core/kinematics/MoveMap.hh>
	#include <core/scoring/Energies.hh>

	#include <devel/init.hh>
	#include <riflib/RotamerGenerator.hh>
	#include <riflib/rosetta_field.hh>
	#include <riflib/util.hh>
	#include <riflib/rotamer_energy_tables.hh>

	// #include <numeric/alignment/QCP_Kernel.hh>
	#include <parallel/algorithm>
	#include <exception>
	#include <stdexcept>

	#include <scheme/actor/Atom.hh>
	#include <scheme/actor/BackboneActor.hh>
	#include <scheme/actor/VoxelActor.hh>
	#include <scheme/kinematics/Director.hh>
	#include <scheme/kinematics/SceneBase.hh>
	#include <scheme/nest/pmap/OriTransMap.hh>
	// #include <scheme/objective/ObjectiveFunction.hh>
	#include <scheme/objective/voxel/FieldCache.hh>
	// #include <scheme/objective/voxel/VoxelArray.hh>
	// #include <scheme/objective/hash/XformMap.hh>
	// #include <scheme/objective/storage/RotamerScores.hh>
	#include <scheme/util/StoragePolicy.hh>
	#include <scheme/search/HackPack.hh>
	#include <scheme/objective/integration/SceneObjective.hh>

	#include <riflib/RifFactory.hh>

	#include <utility/file/file_sys_util.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>

	#include <chrono>

OPT_1GRP_KEY(     StringVector , rif_dock, scaffolds )
	OPT_1GRP_KEY(  StringVector, rif_dock, scaffold_res )
	OPT_1GRP_KEY(  StringVector, rif_dock, scaffold_res_fixed )
	OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_to_ala )
	OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_to_ala_selonly )

	OPT_1GRP_KEY(  StringVector, rif_dock, target_bounding_xmaps )
	OPT_1GRP_KEY(  String      , rif_dock, target_pdb )
	OPT_1GRP_KEY(  String      , rif_dock, target_res )
	OPT_1GRP_KEY(  String      , rif_dock, target_rif )
	OPT_1GRP_KEY(  Real        , rif_dock, target_rf_resl )
	OPT_1GRP_KEY(  Integer     , rif_dock, target_rf_oversample )
	OPT_1GRP_KEY(  String      , rif_dock, target_rf_cache )

	OPT_1GRP_KEY(  StringVector, rif_dock, data_cache_dir )

	OPT_1GRP_KEY(  Real        , rif_dock, beam_size_M )
	OPT_1GRP_KEY(  Real        , rif_dock, search_diameter )
	OPT_1GRP_KEY(  Real        , rif_dock, hsearch_scale_factor )

	OPT_1GRP_KEY(  Real        , rif_dock, max_rf_bounding_ratio )
	OPT_1GRP_KEY(  Boolean     , rif_dock, make_bounding_plot_data )
	OPT_1GRP_KEY(  Boolean     , rif_dock, align_output_to_scaffold )
	OPT_1GRP_KEY(  Integer     , rif_dock, n_pdb_out )

	OPT_1GRP_KEY(  Real        , rif_dock, rf_resl )
	OPT_1GRP_KEY(  Integer     , rif_dock, rf_oversample )
	OPT_1GRP_KEY(  Boolean     , rif_dock, downscale_atr_by_hierarchy )

	OPT_1GRP_KEY(  Integer     , rif_dock, rotrf_oversample )
	OPT_1GRP_KEY(  Real        , rif_dock, rotrf_resl )
	OPT_1GRP_KEY(  Real        , rif_dock, rotrf_spread )
	OPT_1GRP_KEY(  Real        , rif_dock, rotrf_scale_atr )
	OPT_1GRP_KEY(  String      , rif_dock, rotrf_cache_dir )

	OPT_1GRP_KEY(  Boolean     , rif_dock, hack_pack )
	OPT_1GRP_KEY(  Real        , rif_dock, hack_pack_frac )
	OPT_1GRP_KEY(  Real        , rif_dock, pack_iter_mult )
	OPT_1GRP_KEY(  Real        , rif_dock, hbond_weight )
	OPT_1GRP_KEY(  Real        , rif_dock, upweight_multi_hbond )
	OPT_1GRP_KEY(  Real        , rif_dock, global_score_cut )

	OPT_1GRP_KEY(  Integer     , rif_dock, n_result_limit )
	OPT_1GRP_KEY(  Real        , rif_dock, redundancy_filter_mag )

	OPT_1GRP_KEY(  Real        , rif_dock, force_output_if_close_to_input )
	OPT_1GRP_KEY(  Integer     , rif_dock, force_output_if_close_to_input_num )

	OPT_1GRP_KEY(  Real        , rif_dock, upweight_iface )

	OPT_1GRP_KEY(  Boolean     , rif_dock, use_scaffold_bounding_grids )

	OPT_1GRP_KEY(  Boolean     , rif_dock, restrict_to_native_scaffold_res )
	OPT_1GRP_KEY(  Real        , rif_dock, bonus_to_native_scaffold_res )
	OPT_1GRP_KEY(  Boolean     , rif_dock, add_native_scaffold_rots_when_packing )

	OPT_1GRP_KEY(  Boolean     , rif_dock, dump_all_rif_rots )

	OPT_1GRP_KEY(  String     , rif_dock, dokfile )
	OPT_1GRP_KEY(  String     , rif_dock, outdir )
	OPT_1GRP_KEY(  String     , rif_dock, output_tag )

	OPT_1GRP_KEY(  Boolean    , rif_dock, dont_use_scaffold_loops )

	OPT_1GRP_KEY(  Boolean    , rif_dock, full_scaffold_output )
	OPT_1GRP_KEY(  Boolean    , rif_dock, dump_resfile )
	OPT_1GRP_KEY(  Boolean    , rif_dock, pdb_info_pikaa )

	OPT_1GRP_KEY(  Boolean    , rif_dock, cache_scaffold_data )

	OPT_1GRP_KEY(  Real        , rif_dock, tether_to_input_position )

	OPT_1GRP_KEY(  Boolean     , rif_dock, lowres_sterics_cbonly )

	OPT_1GRP_KEY(  Integer     , rif_dock, require_satisfaction )

	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_fraction )		
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_then_min_below_thresh )
	OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_at_least )
	OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_at_most )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_min_fraction )	
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_cut )	



	void register_options() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		NEW_OPT(  rif_dock::scaffolds, "" , utility::vector1<std::string>() );
		NEW_OPT(  rif_dock::scaffold_res, "" , utility::vector1<std::string>() );
		NEW_OPT(  rif_dock::scaffold_res_fixed, "" , utility::vector1<std::string>() );
		NEW_OPT(  rif_dock::scaffold_to_ala, "" , false );
		NEW_OPT(  rif_dock::scaffold_to_ala_selonly, "" , true );

		NEW_OPT(  rif_dock::target_bounding_xmaps, "" , utility::vector1<std::string>() );
		NEW_OPT(  rif_dock::target_pdb, "" , "" );
		NEW_OPT(  rif_dock::target_res, "" , "" );
		NEW_OPT(  rif_dock::target_rif, "" , "" );
		NEW_OPT(  rif_dock::target_rf_resl, ""       , 0.25 );
		NEW_OPT(  rif_dock::target_rf_oversample, "" , 2 );
		NEW_OPT(  rif_dock::downscale_atr_by_hierarchy, "" , true );

		NEW_OPT(  rif_dock::target_rf_cache, "" , "NO_CACHE_SPECIFIED_ON_COMMAND_LINE" );

		NEW_OPT(  rif_dock::data_cache_dir, "" , utility::vector1<std::string>(1,"./") );
		NEW_OPT(  rif_dock::beam_size_M, "" , 10.000000 );
		NEW_OPT(  rif_dock::max_rf_bounding_ratio, "" , 4 );
		NEW_OPT(  rif_dock::make_bounding_plot_data, "" , false );
		NEW_OPT(  rif_dock::align_output_to_scaffold, "" , false );
		NEW_OPT(  rif_dock::n_pdb_out, "" , 10 );

		NEW_OPT(  rif_dock::rf_resl, ""       , 0.25 );
		NEW_OPT(  rif_dock::rf_oversample, "" , 2 );

		NEW_OPT(  rif_dock::rotrf_oversample, "" , 2 );
		NEW_OPT(  rif_dock::rotrf_resl, "" , 0.3 );
		NEW_OPT(  rif_dock::rotrf_spread, "" , 0.0 );
		NEW_OPT(  rif_dock::rotrf_scale_atr, "" , 1.0 );
		NEW_OPT(  rif_dock::rotrf_cache_dir, "" , "./" );

		NEW_OPT(  rif_dock::hack_pack, "" , true );
		NEW_OPT(  rif_dock::hack_pack_frac, "" , 0.2 );
		NEW_OPT(  rif_dock::pack_iter_mult, "" , 2.0 );
		NEW_OPT(  rif_dock::hbond_weight, "" , 2.0 );
		NEW_OPT(  rif_dock::upweight_multi_hbond, "" , 0.0 );
		NEW_OPT(  rif_dock::global_score_cut, "" , 0.0 );

		NEW_OPT(  rif_dock::n_result_limit, "" , 2000000000 );

		NEW_OPT(  rif_dock::redundancy_filter_mag, "" , 1.0 );

		NEW_OPT(  rif_dock::force_output_if_close_to_input, "" , 1.0 );
		NEW_OPT(  rif_dock::force_output_if_close_to_input_num, "" , 0 );

		NEW_OPT(  rif_dock::upweight_iface, "", 1.2 );

		NEW_OPT(  rif_dock::use_scaffold_bounding_grids, "", false );

		NEW_OPT(  rif_dock::search_diameter, "", 150.0 );
		NEW_OPT(  rif_dock::hsearch_scale_factor, "global scaling of rotation/translation search grid", 1.0 );

		NEW_OPT(  rif_dock::restrict_to_native_scaffold_res, "aka structure prediction CHEAT", false );
		NEW_OPT(  rif_dock::bonus_to_native_scaffold_res, "aka favor native CHEAT", -0.3 );
		NEW_OPT(  rif_dock::add_native_scaffold_rots_when_packing, "CHEAT", true );

		NEW_OPT(  rif_dock::dump_all_rif_rots, "", false );

		NEW_OPT(  rif_dock::dokfile, "", "default.dok" );
		NEW_OPT(  rif_dock::outdir, "", "./" );
		NEW_OPT(  rif_dock::output_tag, "", "" );

		NEW_OPT(  rif_dock::dont_use_scaffold_loops, "", false );

		NEW_OPT(  rif_dock::full_scaffold_output, "", false );
		NEW_OPT(  rif_dock::dump_resfile, "", false );
		NEW_OPT(  rif_dock::pdb_info_pikaa, "", false );

		NEW_OPT(  rif_dock::cache_scaffold_data, "", false );

		NEW_OPT(  rif_dock::tether_to_input_position, "", -1.0 );

		NEW_OPT(  rif_dock::lowres_sterics_cbonly, "", true );

		NEW_OPT(  rif_dock::require_satisfaction, "", 0 );

		NEW_OPT(  rif_dock::rosetta_score_fraction  , "",  0.00 );		
		NEW_OPT(  rif_dock::rosetta_score_then_min_below_thresh, "", -9e9 );
		NEW_OPT(  rif_dock::rosetta_score_at_least, "", -1 );
		NEW_OPT(  rif_dock::rosetta_score_at_most, "", 999999999 );
		NEW_OPT(  rif_dock::rosetta_min_fraction  , "",  0.1 );
		NEW_OPT(  rif_dock::rosetta_score_cut  , "", -10.0 );




	}

using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;

// ::scheme::util::SimpleArray<3,float>
Eigen::Vector3f
pose_center(
	core::pose::Pose const & pose,
	utility::vector1<core::Size> const & useres = utility::vector1<core::Size>()
){
	typedef numeric::xyzVector<core::Real> Vec;
	Vec cen(0,0,0);
	int count = 0;
	for( int ir = 1; ir <= pose.n_residue(); ++ir ) {
		if( useres.size()==0 || std::find(useres.begin(),useres.end(),ir)!=useres.end() ){
			for( int ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ){
				cen += pose.xyz(core::id::AtomID(ia,ir));
				++count;
			}
		// } else {
			// std::cout << "pose_center skip " << ir << std::endl;
		}
	}
	cen /= double(count);
	// ::scheme::util::SimpleArray<3,float> center;
	Eigen::Vector3f center;
	center[0] = cen[0];
	center[1] = cen[1];
	center[2] = cen[2];
	return center;
}


void
get_rg_radius(
	core::pose::Pose const & pose,
	float & rg,
	float & radius,
	utility::vector1<core::Size> const & useres = utility::vector1<core::Size>(),
	bool allatom = false
){
	Eigen::Vector3f centmp = pose_center( pose, useres );
	numeric::xyzVector<double> cen;
	float maxdis = -9e9, avgdis2 = 0.0;
	for( int i = 0; i < 3; ++i ) cen[i] = centmp[i];
	for( int iri = 1; iri <= useres.size(); ++iri ){
		int ir = useres[iri];
		if( allatom ){
			for( int ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ){
				numeric::xyzVector<double> coord = pose.residue(ir).xyz(ia);
				avgdis2 += cen.distance_squared( coord );
				maxdis = std::max( maxdis, (float)cen.distance( coord ) );
			}
		} else {
			numeric::xyzVector<double> coord;
			if(      pose.residue(ir).has("CB") ) coord = pose.residue(ir).xyz("CB");
			else if( pose.residue(ir).has("CA") ) coord = pose.residue(ir).xyz("CA");
			else                                  coord = pose.residue(ir).nbr_atom_xyz();
			avgdis2 += cen.distance_squared( coord );
			maxdis = std::max( maxdis, (float)cen.distance( coord ) );
		}
	}
	avgdis2 /= useres.size();
	rg = sqrt( avgdis2 );
	radius = maxdis;
}



#pragma pack (push, 4) // allows size to be 12 rather than 16
struct SearchPoint {
	float score;
	uint64_t index;
	SearchPoint() : score(9e9), index(0) {}
	SearchPoint(uint64_t i) : score(9e9), index(i) {}
	bool operator < (SearchPoint const & o) const {
		return score < o.score;
	}
};
#pragma pack (pop)

struct SearchPointWithRots {
	float score;
	uint32_t prepack_rank;
	uint64_t index;
	shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
	core::pose::PoseOP pose_ = nullptr;
	SearchPointWithRots() : score(9e9), prepack_rank(0), index(0), rotamers_(nullptr) {}
	SearchPointWithRots(uint64_t i, uint32_t orank) : score(9e9), prepack_rank(orank), index(i), rotamers_(nullptr) {}
	// ~SearchPointWithRots() { delete rotamers_; }
	void checkinit() { if( rotamers_==nullptr ) rotamers_ = make_shared< std::vector< std::pair<intRot,intRot> > > ();	}
	std::vector< std::pair<intRot,intRot> > & rotamers() { checkinit(); return *rotamers_; }
	std::vector< std::pair<intRot,intRot> > const & rotamers() const { runtime_assert(rotamers_!=nullptr); return *rotamers_; }
	size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
	bool operator < (SearchPointWithRots const & o) const {
		return score < o.score;
	}
	friend void swap(SearchPointWithRots & a, SearchPointWithRots & b){
        std::swap( a.score, b.score );
        std::swap( a.prepack_rank, b.prepack_rank );
        std::swap( a.index, b.index );
        std::swap( a.rotamers_, b.rotamers_ );
        std::swap( a.pose_, b.pose_ );
    }
};


void xform_pose( core::pose::Pose & pose, numeric::xyzTransform<float> s, core::Size sres=1, core::Size eres=0 ) {
  if(eres==0) eres = pose.n_residue();
  for(core::Size ir = sres; ir <= eres; ++ir) {
    for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s*pose.xyz(aid) );
    }
  }
}


struct RifDockResult {
	float dist0, packscore, nopackscore, rifscore, stericscore;
	uint64_t isamp, scene_index;
	uint32_t prepack_rank;
	float cluster_score;
	bool operator< ( RifDockResult const & o ) const { return packscore < o.packscore; }
	shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
	core::pose::PoseOP pose_ = nullptr;
	size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
	std::vector< std::pair<intRot,intRot> > const & rotamers() const { assert(rotamers_!=nullptr); return *rotamers_; }
};

// how can I fix this??? make the whole prototype into a class maybe???
// what does it do?
// 	set and rescore scene with nopackscore, record more score detail
//  compute dist0
//	select results with some redundancy filtering
template<
	class EigenXform,
	// class Scene,
	class ScenePtr,
	class ObjectivePtr
>
void
awful_compile_output_helper(
	int64_t isamp,
	int resl,
	std::vector< SearchPointWithRots > const & packed_results,
	std::vector< ScenePtr > & scene_pt,
	shared_ptr< ::scheme::kinematics::Director<EigenXform> > director,
	float redundancy_filter_rg,
	float redundancy_filter_mag,
	Eigen::Vector3f scaffold_center,
	std::vector< std::vector< RifDockResult > > & allresults_pt,
	             std::vector< RifDockResult >   & selected_results,
	std::vector< std::pair< EigenXform, int64_t > > & selected_xforms,
	int n_pdb_out,
	#ifdef USE_OPENMP
		omp_lock_t & dump_lock,
	#endif
	ObjectivePtr objective,
	int & nclose,
	int nclosemax,
	float nclosethresh
){
	SearchPointWithRots const & sp = packed_results[isamp];
	if( sp.score >= 0.0f ) return;
	ScenePtr scene_minimal( scene_pt[omp_get_thread_num()] );
	director->set_scene( sp.index, resl, *scene_minimal );
	std::vector<float> sc = objective->scores(*scene_minimal);
	float const nopackscore = sc[0]+sc[1]; //result.sum();
	float const rifscore = sc[0]; //result.template get<MyScoreBBActorRIF>();
	float const stericscore = sc[1]; //result.template get<MyClashScore>();
	float dist0; {
		EigenXform x = scene_minimal->position(1);
		x.translation() -= scaffold_center;
		dist0 = ::devel::scheme::xform_magnitude( x, redundancy_filter_rg );
	}

	RifDockResult r; // float dist0, packscore, nopackscore, rifscore, stericscore;
	r.isamp = isamp;
	r.prepack_rank = sp.prepack_rank;
	r.scene_index = sp.index;
	r.packscore = sp.score;
	r.nopackscore = nopackscore;
	r.rifscore = rifscore;
	r.stericscore = stericscore;
	r.dist0 = dist0;
	r.cluster_score = 0.0;
	r.pose_ = sp.pose_;
	allresults_pt.at( omp_get_thread_num() ).push_back( r ); // recorded w/o rotamers here

	bool force_selected = ( dist0 < nclosethresh && ++nclose < nclosemax ); // not thread-safe... is this important?

	if( selected_results.size() < n_pdb_out || force_selected ){

		EigenXform xposition1 = scene_minimal->position(1);
		EigenXform xposition1inv = xposition1.inverse();

		float mindiff_candidate = 9e9;
		int64_t i_closest_result;
		typedef std::pair< EigenXform, int64_t > XRpair;
		BOOST_FOREACH( XRpair const & xrp, selected_xforms ){
			EigenXform const & xsel = xrp.first;
			EigenXform const xdiff = xposition1inv * xsel;
			float diff = devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg );
			if( diff < mindiff_candidate ){
				mindiff_candidate = diff;
				i_closest_result = xrp.second;
			}
			// todo: also compare AA composition of rotamers
		}

		if( mindiff_candidate < redundancy_filter_mag ){ // redundant result
			selected_results[i_closest_result].cluster_score += 1.0; //sp.score==0.0 ? nopackscore : sp.score;
		}

		if( mindiff_candidate > redundancy_filter_mag || force_selected ){

			#ifdef USE_OPENMP
			omp_set_lock( &dump_lock );
			#endif
			{
				// std::cout << "checking again to add selected " << selected_xforms.size() << " " << omp_get_thread_num() << std::endl;
				float mindiff_actual = 9e9;
				BOOST_FOREACH( XRpair const & xrp, selected_xforms ){
					EigenXform const & xsel = xrp.first;
					EigenXform const xdiff = xposition1inv * xsel;
					mindiff_actual = std::min( mindiff_actual, devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg ) );
				}
				if( mindiff_actual > redundancy_filter_mag || force_selected ){
					if( redundancy_filter_mag > 0.0001 )
						selected_xforms.push_back( std::make_pair(xposition1,(int64_t)selected_results.size()) );
					r.rotamers_ = sp.rotamers_;
					selected_results.push_back( r ); // recorded with rotamers here
				} else {
					// std::cout << " second check failed" << std::endl;
				}
			}
			#ifdef USE_OPENMP
			omp_unset_lock( &dump_lock );
			#endif

		} // end if( mindiff > redundancy_filter_mag ){

	} // end 	if( selected_xforms.size() < n_pdb_out || force_selected )

}

int main(int argc, char *argv[]) {

	#ifdef USE_OPENMP
		omp_lock_t cout_lock, dump_lock;
		omp_init_lock( &cout_lock );
		omp_init_lock( &dump_lock );
	#endif

	register_options();
	devel::init(argc,argv);

	using basic::options::option;
		using namespace basic::options::OptionKeys;
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


		typedef ::scheme::nest::NEST< 6,
									  EigenXform,
									  ::scheme::nest::pmap::OriTransMap,
									  ::scheme::util::StoreNothing, // do not store a transform in the Nest
									  uint64_t,
									  float,
									  false // do not inherit from NestBase
									 > NestOriTrans6D;

		typedef ::scheme::kinematics::NestDirector< NestOriTrans6D > DirectorOriTrans6D;
		typedef shared_ptr< ::scheme::kinematics::Director<EigenXform> > DirectorBase;

		typedef ::scheme::actor::BackboneActor<EigenXform> BBActor;
		typedef ::scheme::actor::VoxelActor<EigenXform,float> VoxelActor;

		typedef ::scheme::actor::SimpleAtom< Eigen::Vector3f > SimpleAtom;


	print_header( "setup global options" );
	runtime_assert( option[rif_dock::target_rif].user() );
		std::vector<std::string> data_cache_path;
		bool lowres_sterics_cbonly = option[rif_dock::lowres_sterics_cbonly]();
		float tether_to_input_position_cut = option[rif_dock::tether_to_input_position]();
		bool tether_to_input_position = tether_to_input_position_cut > 0.0;
		float global_score_cut = option[rif_dock::global_score_cut]();
		BOOST_FOREACH( std::string dir, option[rif_dock::data_cache_dir]() ) data_cache_path.push_back( dir );
		std::string const outdir = option[rif_dock::outdir]();
		std::string const output_tag = option[rif_dock::output_tag]();
		utility::file::create_directory_recursive( outdir );
		std::string dokfile_fname = outdir + "/" + option[rif_dock::dokfile]();
		bool const dump_all_rif_rots = option[ rif_dock::dump_all_rif_rots ]();
		bool const add_native_scaffold_rots_when_packing = option[ rif_dock::add_native_scaffold_rots_when_packing ]();
		bool const restrict_to_native_scaffold_res = option[ rif_dock::restrict_to_native_scaffold_res ]();
		float const bonus_to_native_scaffold_res = option[ rif_dock::bonus_to_native_scaffold_res ]();
		float const hack_pack_frac = option[ rif_dock::hack_pack_frac ]();
		float const hsearch_scale_factor = option[ rif_dock::hsearch_scale_factor ]();
		float const search_diameter = option[ rif_dock::search_diameter ]();
		bool const use_scaffold_bounding_grids = option[rif_dock::use_scaffold_bounding_grids]();
		double const resl0 = 16.0;
		float const upweight_iface = option[rif_dock::upweight_iface]();
		std::vector<std::string> rif_files;
		BOOST_FOREACH( std::string fn, option[rif_dock::target_bounding_xmaps]() ){
			rif_files.push_back( fn );
		}
		rif_files.push_back( option[rif_dock::target_rif]() );
		std::string rif_type = get_rif_type_from_file( rif_files.back() );
		BOOST_FOREACH( std::string fn, rif_files ){
			std::string rif_type2 = get_rif_type_from_file( fn );
			runtime_assert_msg( rif_type==rif_type2, "mismatched rif types, expect: " + rif_type + " got: " + rif_type2 + " for " + fn );
		}
		std::cout << "read RIF type: " << rif_type << std::endl;
		int const require_satisfaction = option[rif_dock::require_satisfaction]();

		cout << "Search Resls: " << resl0;
		std::vector<float> RESLS(1,resl0);
		for( int i = 1; i <= option[rif_dock::target_bounding_xmaps]().size(); ++i ){
			RESLS.push_back( RESLS.back()/2.0 );
			cout << " " << RESLS.back();
		}
		cout << endl;
		int64_t const DIM = 6;
		int64_t const DIMPOW2 = 1<<DIM;
		int64_t const beam_size = int64_t( option[rif_dock::beam_size_M]() * 1000000.0 / DIMPOW2 ) * DIMPOW2;

		bool const VERBOSE = false;

		bool const scaff2ala        = option[rif_dock::scaffold_to_ala]();
		bool const scaff2alaselonly = option[rif_dock::scaffold_to_ala_selonly]();
		if( scaff2ala && scaff2alaselonly &&  option[rif_dock::scaffold_to_ala_selonly].user() ){
			std::cout << "WARNING: -scaffold_to_ala overrides -scaffold_to_ala_selonly!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		}
		bool const replace_all_with_ala_1bre = true; // todo: decide how to handle this

		::scheme::search::HackPackOpts hackpackopts;
		hackpackopts.pack_iter_mult = option[rif_dock::pack_iter_mult]();
		hackpackopts.hbond_weight   = option[rif_dock::hbond_weight]();
		hackpackopts.upweight_iface = upweight_iface;
		hackpackopts.upweight_multi_hbond = option[rif_dock::upweight_multi_hbond]();

		float const target_rf_resl = option[rif_dock::target_rf_resl]()<=0.0 ? RESLS.back()/2.0 : option[rif_dock::target_rf_resl]();

		bool align_to_scaffold = option[rif_dock::align_output_to_scaffold]();
		float rosetta_score_fraction = option[rif_dock::rosetta_score_fraction]();				
		float rosetta_score_then_min_below_thresh = option[rif_dock::rosetta_score_then_min_below_thresh]();
		float rosetta_score_at_least = option[rif_dock::rosetta_score_at_least]();
		float rosetta_score_at_most  = option[rif_dock::rosetta_score_at_most]();		
		float rosetta_min_fraction = option[rif_dock::rosetta_min_fraction]();		
		float rosetta_score_cut = option[rif_dock::rosetta_score_cut]();		
		std::cout << "rosetta_score_fraction: " << rosetta_score_fraction << std::endl;		
		std::cout << "rosetta_score_then_min_below_thresh: " << rosetta_score_then_min_below_thresh << std::endl;
		std::cout << "rosetta_score_at_least: " << rosetta_score_at_least << std::endl;
		std::cout << "rosetta_score_at_most: " << rosetta_score_at_most << std::endl;		
		std::cout << "rosetta_min_fraction: " << rosetta_min_fraction << std::endl;
		std::cout << "rosetta_score_cut: " << rosetta_score_cut << std::endl;
		double time_rif=0, time_pck=0, time_ros=0;

		bool pdb_info_pikaa = option[rif_dock::pdb_info_pikaa]();

		std::cout << "//////////////////////////// end options /////////////////////////////////" << std::endl;

		////////////////////////////// should be no more use of options at this point! ///////////////////////////





		utility::io::ozstream dokout( dokfile_fname );



		devel::scheme::RifFactoryConfig rif_factory_config;
		rif_factory_config.rif_type = rif_type;
		shared_ptr<RifFactory> rif_factory = ::devel::scheme::create_rif_factory( rif_factory_config );



	print_header( "create rotamer index" );
		shared_ptr< RotamerIndex > rot_index_p = make_shared< RotamerIndex >();
		RotamerIndex & rot_index( *rot_index_p );
		::devel::scheme::get_rotamer_index( rot_index );
		std::cout << "================ RotamerIndex ===================" << std::endl;
		std::cout << rot_index << std::endl;
		// {
		// 	utility::io::ozstream out("rot_index.pdb");
		// 	rot_index.dump_pdb( out );
		// 	utility_exit_with_message("ortsdn");
		// }
		std::cout << "=================================================" << std::endl;

		RotamerRFOpts rotrfopts;
		rotrfopts.oversample     = option[rif_dock::rotrf_oversample]();
		rotrfopts.field_resl     = option[rif_dock::rotrf_resl]();
		rotrfopts.field_spread   = option[rif_dock::rotrf_spread]();
		rotrfopts.data_dir       = option[rif_dock::rotrf_cache_dir]();
		rotrfopts.scale_atr      = option[rif_dock::rotrf_scale_atr]();
		::devel::scheme::RotamerRFTablesManager rotrftables( rot_index_p, rotrfopts );
		// rotrftables.preinit_all();




	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	print_header( "read and prepare target structure" ); //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	core::pose::Pose target;
	std::string target_fname = option[rif_dock::target_pdb]();
	std::vector<SimpleAtom> target_simple_atoms;
	utility::vector1<core::Size> target_res;
	std::vector<HBondRay> target_donors, target_acceptors;
	float rif_radius=0.0, target_redundancy_filter_rg=0.0;
	{
		core::import_pose::pose_from_file( target, target_fname );

		if( use_scaffold_bounding_grids ){
			for( int ir = 1; ir <= target.n_residue(); ++ir ){
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
		std::string target_res_fname = "";
		if( option[rif_dock::target_res].user() ){
			target_res_fname = option[rif_dock::target_res]();
		}
		target_res = devel::scheme::get_res( target_res_fname , target );
		get_rg_radius( target, target_redundancy_filter_rg, rif_radius, target_res, true ); // allatom for target
		rif_radius += 7.0; // hacky guess
		BOOST_FOREACH( core::Size ir, target_res ){
			::devel::scheme::get_donor_rays   ( target, ir, true, target_donors );
			::devel::scheme::get_acceptor_rays( target, ir, true, target_acceptors );
		}
		std::cout << "target_donors.size() " << target_donors.size() << " target_acceptors.size() " << target_acceptors.size() << std::endl;
	}
	std::vector< VoxelArrayPtr > target_field_by_atype;
	std::vector< std::vector< VoxelArrayPtr > > target_bounding_by_atype;
	{
		target_bounding_by_atype.resize( RESLS.size() );
		devel::scheme::RosettaFieldOptions rfopts;
		rfopts.field_resl = target_rf_resl;
		rfopts.data_dir = "DUMMY_DATA_DIR_FIXME";
		rfopts.oversample = option[rif_dock::target_rf_oversample]();
		rfopts.block_hbond_sites = false;
		rfopts.max_bounding_ratio = option[rif_dock::max_rf_bounding_ratio]();
		rfopts.fail_if_no_cached_data = true;
		rfopts.repulsive_only_boundary = true;
		rfopts.cache_mismatch_tolerance = 0.01; // this is kinda loose...
		std::string cache_prefix = option[rif_dock::target_rf_cache]();
		devel::scheme::get_rosetta_fields_specified_cache_prefix(
			cache_prefix,
			target_fname,
			target,
			target_res,
			rfopts,
			target_field_by_atype,
			false
		);


		// if( use_scaffold_bounding_grids ){
		// 	std::cout << "not using target steric grids" << std::endl;
		// } else {
		if( true ){
			// std::cout << "using target bounding grids, generating (or loading) them" << std::endl;
			devel::scheme::RosettaFieldOptions rfopts;
			rfopts.field_resl = target_rf_resl;
			rfopts.data_dir = "DUMMY_DATA_DIR_FIXME";
			rfopts.oversample = option[rif_dock::target_rf_oversample]();
			rfopts.block_hbond_sites = false;
			rfopts.max_bounding_ratio = option[rif_dock::max_rf_bounding_ratio]();
			rfopts.fail_if_no_cached_data = true;
			rfopts.repulsive_only_boundary = true; // default
			devel::scheme::get_rosetta_bounding_fields_from_fba(
				RESLS,
				target_fname,
				target,
				target_res,
				rfopts,
				target_field_by_atype,
				target_bounding_by_atype,
				false,
				cache_prefix
			);
			runtime_assert( target_bounding_by_atype.size() == RESLS.size() );
			// now scale down the any positive component by 1/RESL if RESL > 1
			if( option[rif_dock::downscale_atr_by_hierarchy]() ){
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

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	print_header( "read in RIFs" ); /////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<shared_ptr<RifBase> > rif_ptrs;
	std::vector<bool> rif_using_rot;
	{
		std::vector<std::string> rif_descriptions( rif_files.size() );
		rif_ptrs.resize( rif_files.size() );
		std::exception_ptr exception = nullptr;
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for( int i_readmap = 0; i_readmap < rif_files.size(); ++i_readmap ){
			if( exception ) continue;
			try {
				std::string const & rif_file = rif_files[i_readmap];
				std::string & rif_dscr = rif_descriptions[i_readmap];
				shared_ptr<RifBase> & rif_ptr = rif_ptrs[i_readmap];
				rif_ptr = rif_factory->create_rif_from_file( rif_file, rif_dscr );
				runtime_assert_msg( rif_ptrs[i_readmap] , "rif creation from file failed! " + rif_file );
				if( VERBOSE ){
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
	}



	for( int iscaff = 1; iscaff <= option[rif_dock::scaffolds]().size(); ++iscaff )
	{
		std::string scaff_fname = option[rif_dock::scaffolds]().at(iscaff);
		std::vector<std::string> scaffold_sequence_glob0;
		utility::vector1<core::Size> scaffold_res;//, scaffold_res_all;
		try {

			runtime_assert( rot_index_p );
			std::string scafftag = utility::file_basename( utility::file::file_basename( scaff_fname ) );

			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////   begin scaffold " << scafftag << " " << iscaff << " of " << option[rif_dock::scaffolds]().size() << std::endl;
			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;

			core::pose::Pose scaffold, scaffold_centered, scaffold_full_centered, both_pose, both_full_pose;
			float scaff_radius = 0.0;
			float redundancy_filter_rg = 0.0;

			std::vector<int> scaffres_g2l, scaffres_l2g;
			std::vector<bool> scaffuseres;
			Eigen::Vector3f scaffold_center;
			std::vector<Vec> scaffca;
			std::vector<std::vector<float> > scaffold_onebody_glob0;
			std::vector<std::vector<float> > local_onebody;
			std::vector< std::pair<int,int> > local_rotamers;
			typedef ::scheme::objective::storage::TwoBodyTable<float> TBT;
			shared_ptr<TBT> scaffold_twobody = make_shared<TBT>( scaffold.n_residue(), rot_index.size()  );
			shared_ptr<TBT> local_twobody;
			{
				core::import_pose::pose_from_file( scaffold, scaff_fname );
				scaffold_full_centered = scaffold;

				for( int ir = 1; ir <= scaffold.n_residue(); ++ir ){
					scaffold_sequence_glob0.push_back( scaffold.residue(ir).name3() );
				}

				std::string scaff_res_fname = "";
				if( option[rif_dock::scaffold_res].user() ){
					if( option[rif_dock::scaffold_res]().size() == option[rif_dock::scaffolds]().size() ){
						scaff_res_fname = option[rif_dock::scaffold_res]().at(iscaff);
					} else if( option[rif_dock::scaffold_res]().size() == 1 ){
						scaff_res_fname = option[rif_dock::scaffold_res]().at(1);
					} else {
						utility_exit_with_message( "-scaffold_res list not same length as -scaffolds list" );
					}
					scaffold_res = devel::scheme::get_res( scaff_res_fname , scaffold );
				} else {
					scaffold_res = devel::scheme::get_res_by_sasa( scaffold, option[rif_dock::dont_use_scaffold_loops]() );
				}
				if     ( scaff2ala )        ::devel::scheme::pose_to_ala( scaffold );
				else if( scaff2alaselonly ) ::devel::scheme::pose_to_ala( scaffold, scaffold_res );

				float scaff_redundancy_filter_rg=0;
				get_rg_radius( scaffold, scaff_redundancy_filter_rg, scaff_radius, scaffold_res, false ); // not allatom for scaff
				redundancy_filter_rg = std::min( scaff_redundancy_filter_rg, target_redundancy_filter_rg );
				std::cout << "scaffold selected region rg: " << scaff_redundancy_filter_rg << ", radius: " << scaff_radius << std::endl;
				std::cout << "using redundancy_filter_rg: " << redundancy_filter_rg << std::endl;

				int count = 0;
				scaffres_g2l.resize(scaffold.n_residue(),-1);
				scaffuseres .resize(scaffold.n_residue(),false);
				for( auto ir : scaffold_res ){
					scaffres_g2l[ir-1] = count++;
					scaffres_l2g.push_back(ir-1);
					scaffuseres[ir-1] = true;
				}

				std::cout << "scaffold: " << scaff_fname << " nres: " << scaffold.n_residue() << " using_res: " << scaffold_res.size() << std::endl;
				scaffold_center = pose_center(scaffold,scaffold_res);
				scaffold_centered = scaffold;
				for( int ir = 1; ir <= scaffold.n_residue(); ++ir ){
					Vec tmp( scaffold_center[0], scaffold_center[1], scaffold_center[2] );
					for( int ia = 1; ia <= scaffold.residue_type(ir).natoms(); ++ia ){
						core::id::AtomID aid(ia,ir);
						scaffold_centered.set_xyz( aid, scaffold.xyz(aid) - tmp );
					}
					for( int ia = 1; ia <= scaffold_full_centered.residue_type(ir).natoms(); ++ia ){
						core::id::AtomID aid(ia,ir);
						scaffold_full_centered.set_xyz( aid, scaffold_full_centered.xyz(aid) - tmp );
					}
				}

				both_pose      = scaffold_centered;
				both_full_pose = scaffold_full_centered;
				core::pose::append_pose_to_pose( both_pose, target );
				core::pose::append_pose_to_pose( both_full_pose, target );
				runtime_assert( both_pose.n_residue() == scaffold.n_residue() + target.n_residue() );
				runtime_assert( both_pose.n_residue() == both_full_pose.n_residue() );



				for( int ir = 1; ir <= scaffold.n_residue(); ++ir ){
					scaffca.push_back( scaffold.residue(ir).xyz("CA") );
					// scaffold_res_all.push_back(ir);
				}
				std::string scaff_tag = utility::file_basename( scaff_fname );
				std::string cachefile = "__1BE_" + scaff_tag + (replace_all_with_ala_1bre?"_ALLALA":"") + ".bin.gz";
				if( ! option[rif_dock::cache_scaffold_data]() ) cachefile = "";
				get_onebody_rotamer_energies( scaffold, rot_index, scaffold_onebody_glob0, data_cache_path, cachefile, replace_all_with_ala_1bre );

				if( restrict_to_native_scaffold_res ){
					std::cout << "KILLING NON-NATIVE ROTAMERS ON SCAFFOLD!!!" << std::endl;
					for( int ir = 0; ir < scaffold_onebody_glob0.size(); ++ir ){
						for( int irot = 0; irot < rot_index.size(); ++irot ){
							if( rot_index.resname(irot) != scaffold_sequence_glob0.at(ir) && rot_index.resname(irot) != "ALA" ){
								scaffold_onebody_glob0[ir][irot] = 9e9;
							}
						}
					}
				}
				if( bonus_to_native_scaffold_res != 0 ){
					std::cout << "adding to native scaffold res 1BE " << bonus_to_native_scaffold_res << std::endl;
					for( int ir = 0; ir < scaffold_onebody_glob0.size(); ++ir ){
						for( int irot = 0; irot < rot_index.size(); ++irot ){
							if( rot_index.resname(irot) == scaffold_sequence_glob0.at(ir) ){
								scaffold_onebody_glob0[ir][irot] += bonus_to_native_scaffold_res;
							}
						}
					}
				}

				for( int i = 0; i < scaffres_l2g.size(); ++i ){
					local_onebody.push_back( scaffold_onebody_glob0.at( scaffres_l2g.at(i) ) );
				}
				for( int i = 0; i < scaffres_g2l.size(); ++i ){
					if( scaffres_g2l[i] < 0 ){
						BOOST_FOREACH( float & f, scaffold_onebody_glob0[i] ) f = 9e9;
					}
				}
				std::string cachefile2b = "__2BE_" + scaff_tag + "_GLOBAL" + ".bin.gz";
				if( ! option[rif_dock::cache_scaffold_data]() ) cachefile2b = "";
				MakeTwobodyOpts make2bopts;
				make2bopts.onebody_threshold = 2.0;
				make2bopts.distance_cut = 15.0;
				make2bopts.hbond_weight = hackpackopts.hbond_weight;
				std::string dscrtmp;
				get_twobody_tables( data_cache_path, cachefile2b, dscrtmp, scaffold, rot_index, scaffold_onebody_glob0, rotrftables, make2bopts, *scaffold_twobody );
				std::cout << "twobody memuse: " << (float)scaffold_twobody->twobody_mem_use()/1000.0/1000.0 << "M" << std::endl;
				for( int i = 0; i < scaffold_onebody_glob0.size(); ++i ){
					runtime_assert( scaffold_onebody_glob0[i].size() == rot_index.size() );
					for( int j = 0; j < scaffold_onebody_glob0[i].size(); ++j ){
						scaffold_onebody_glob0[i][j] = rif_using_rot[j] ? scaffold_onebody_glob0[i][j] : 9e9;
					}
				}
				for( int i = 0; i < local_onebody.size(); ++i ){
					runtime_assert( local_onebody[i].size() == rot_index.size() );
					for( int j = 0; j < local_onebody[i].size(); ++j ){
						local_onebody[i][j] = rif_using_rot[j] ? local_onebody[i][j] : 9e9;
					}
				}

				local_rotamers.clear();
				for( int i = 0; i < local_onebody.size(); ++i ){
					int iresglobal = scaffres_l2g.at(i);
					std::string name3 = scaffold_sequence_glob0.at(iresglobal);
					std::pair<int,int> ib = rot_index.index_bounds( name3 );
					// std::cout << "local_rotamers " << i << " " << iresglobal << " " << name3 << " " << ib.first << " " << ib.second << std::endl;
					local_rotamers.push_back( ib );
				}

				local_twobody = scaffold_twobody->create_subtable( scaffuseres, scaffold_onebody_glob0, 2.0 );
				std::cout << "filt_2b memuse: " << (float)local_twobody->twobody_mem_use()/1000.0/1000.0 << "M" << std::endl;
				std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! must fix issue with non-global 2B table calculation, seems to use scaffold_res when it shouldn't" << endl;
				std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

				// todo: prune twobody table???
			}

			// SOMETHING WRONG, SCORES OFF BY A LITTLE
			// setup objectives, moved into scaffold loop to guarantee clean slate for each scaff...
			RifSceneObjectiveConfig rso_config;
				rso_config.hackpackopts = &hackpackopts;
				rso_config.rif_ptrs = rif_ptrs;
				rso_config.target_bounding_by_atype = &target_bounding_by_atype;
				rso_config.target_field_by_atype = &target_field_by_atype;
				rso_config.local_onebody = &local_onebody;
				rso_config.local_rotamers = &local_rotamers;
				rso_config.local_twobody = local_twobody;
				rso_config.rot_index_p = rot_index_p;
				rso_config.target_donors = &target_donors;
				rso_config.target_acceptors = &target_acceptors;
				rso_config.add_native_scaffold_rots_when_packing = add_native_scaffold_rots_when_packing;
				rso_config.n_sat_groups = target_donors.size() + target_acceptors.size();
				rso_config.require_satisfaction = require_satisfaction;

			ScenePtr scene_prototype;
			std::vector< ObjectivePtr > objectives;
			ObjectivePtr packing_objective;
			runtime_assert( rif_factory->create_objectives( rso_config, objectives, packing_objective ) );
			scene_prototype = rif_factory->create_scene();
			runtime_assert_msg( objectives.front()->is_compatible( *scene_prototype ), "objective and scene types not compatible!" );





			print_header( "setup 3D rosetta_field grids for scaffold" );
			std::vector< VoxelArrayPtr > scaffold_field_by_atype;
			std::vector< std::vector< VoxelArrayPtr > > scaffold_bounding_by_atype;
			std::vector< SimpleAtom > scaffold_simple_atoms, scaffold_simple_atoms_all;
			if( use_scaffold_bounding_grids ){
				scaffold_bounding_by_atype.resize( RESLS.size() );
				float const rf_resl = option[rif_dock::rf_resl]()==0.0 ? RESLS.back()/2.0 : option[rif_dock::rf_resl]();
				devel::scheme::RosettaFieldOptions rfopts;
				rfopts.field_resl = rf_resl;
				rfopts.data_dir = "DUMMY_DATA_DIR_FIXME";
				rfopts.oversample = option[rif_dock::rf_oversample]();
				rfopts.block_hbond_sites = false;
				rfopts.max_bounding_ratio = option[rif_dock::max_rf_bounding_ratio]();
				rfopts.repulsive_only_boundary = true; // default
				devel::scheme::get_rosetta_bounding_fields(
					RESLS,
					scaff_fname+"_CEN"+(scaff2ala?"_ALLALA":""),
					scaffold_centered,
					scaffold_res,
					rfopts,
					scaffold_field_by_atype,
					scaffold_bounding_by_atype,
					false
				);
				runtime_assert( scaffold_bounding_by_atype.size() == RESLS.size() );
				// now scale down the any positive component by 1/RESL if RESL > 1
				if( option[rif_dock::downscale_atr_by_hierarchy]() ){
					std::cout << "downscale_atr_by_hierarchy on scaffold bounding steric grids" << std::endl;
					// std::cout << "zeroing atr component of scaffold steric grids" << std::endl;
					for( int iresl = 0; iresl < RESLS.size(); ++iresl ){
						float correction = 1.0/RESLS[iresl];
						if( correction >= 1.0 ) break;
						BOOST_FOREACH( VoxelArrayPtr vap, scaffold_bounding_by_atype[iresl] ){
							if( vap != nullptr ){
								#ifdef USE_OPENMP
								#pragma omp parallel for schedule(dynamic,64)
								#endif
								for( int k = 0; k < vap->num_elements(); ++k ){
									float & dat = vap->data()[k];
									if( dat < 0 ){
										dat = dat * correction;
									}
								}
							}
						}
					}
				}
			} else {
				std::cout << "not using scaffold bounding grids" << std::endl;
				for( int ir = 1; ir <= scaffold_centered.n_residue(); ++ir ){
					utility::vector1<core::Size> resids(1,ir); // 1-index numbering
					{
						std::vector<SchemeAtom> scaff_res_atoms;
						if( !lowres_sterics_cbonly && std::find( scaffold_res.begin(), scaffold_res.end(), ir ) != scaffold_res.end() ){
							devel::scheme::get_scheme_atoms( scaffold_centered, resids, scaff_res_atoms, true );
						} else { // is not selected residue
							devel::scheme::get_scheme_atoms_cbonly( scaffold_centered, resids, scaff_res_atoms );
						}
						int restype = rot_index.chem_index_.resname2num( scaffold_centered.residue(ir).name3() ); // for UNK will be -1
						for( int ia = 0; ia < scaff_res_atoms.size(); ++ia){
							SchemeAtom const & a( scaff_res_atoms[ia] );
							runtime_assert( a.type() > 0 );
							if( a.type() >= 21 ) continue;
							SimpleAtom sa( a.position(), a.type(), restype, ia );
							scaffold_simple_atoms.push_back(sa);
						}
					}
					{
						std::vector<SchemeAtom> all_scaff_res_atoms;
						devel::scheme::get_scheme_atoms( scaffold_centered, resids, all_scaff_res_atoms, false );
						int restype = rot_index.chem_index_.resname2num( scaffold_centered.residue(ir).name3() ); // for UNK will be -1
						for( int ia = 0; ia < all_scaff_res_atoms.size(); ++ia){
							SchemeAtom const & a( all_scaff_res_atoms[ia] );
							runtime_assert( a.type() > 0 );
							if( a.type() >= 21 ) continue;
							SimpleAtom sa( a.position(), a.type(), restype, ia );
							scaffold_simple_atoms_all.push_back(sa);
						}
					}
				}
				std::cout << "scaffold_simple_atoms " << scaffold_simple_atoms.size() << std::endl;

			}

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "setup scene from scaffold and target" );
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			ScenePtr scene_minimal( scene_prototype->clone_deep() );
			ScenePtr scene_full( scene_prototype->clone_deep() );
			{
				for( int ir = 1; ir <= scaffold.n_residue(); ++ir ){
					Vec N  = scaffold_centered.residue(ir).xyz("N" );
					Vec CA = scaffold_centered.residue(ir).xyz("CA");
					Vec C  = scaffold_centered.residue(ir).xyz("C" );

					// todo map res indices, must also edit onebody_energies
					BBActor bbactor( N, CA, C, '-', '-', scaffres_g2l[ir-1] );
					runtime_assert( bbactor.index_ == scaffres_g2l[ir-1] );


					scene_full->add_actor(1,bbactor);
					if( std::find(scaffold_res.begin(),scaffold_res.end(),ir)!=scaffold_res.end() ){
						scene_minimal->add_actor(1,bbactor);
					}
				}

				if( use_scaffold_bounding_grids ){
					BOOST_FOREACH( SimpleAtom const & sa, target_simple_atoms )	scene_minimal->add_actor( 0, sa );
					runtime_assert( scene_minimal->template num_actors<SimpleAtom>(0) == target_simple_atoms.size() );
					scene_minimal->add_actor( 1, VoxelActor(scaffold_bounding_by_atype) );
				} else {
					BOOST_FOREACH( SimpleAtom const & sa, scaffold_simple_atoms ) scene_minimal->add_actor( 1, sa );
					runtime_assert( scene_minimal->template num_actors<SimpleAtom>(1) == scaffold_simple_atoms.size() );
					scene_minimal->add_actor( 0, VoxelActor(target_bounding_by_atype) );
				}


			}
			cout << "scores for scaffold in original position: " << std::endl;
			{
				EigenXform x(EigenXform::Identity());
				x.translation() = scaffold_center;
				scene_minimal->set_position(1,x);
				for(int i = 0; i < RESLS.size(); ++i){
					std::vector<float> sc = objectives[i]->scores(*scene_minimal);
					cout << "input bounding score " << i << " " << F(7,3,RESLS[i]) << " "
					     << F( 7, 3, sc[0]+sc[1] ) << " "
					     << F( 7, 3, sc[0]       ) << " "
					     << F( 7, 3, sc[1]       ) << endl;
				}

			}

			// utility_exit_with_message("FOO");

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "setup director based on scaffold and target sizes" ); //////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			DirectorBase director; {
				F3 target_center = pose_center(target);
				float const body_radius = std::min( scaff_radius, rif_radius );
				double const cart_grid = resl0*hsearch_scale_factor/sqrt(3); // 1.5 is a big hack here.... 2 would be more "correct"
				double const hackysin = std::min( 1.0, resl0*hsearch_scale_factor/2.0/ body_radius );
				runtime_assert( hackysin > 0.0 );
				double const rot_resl_deg0 = asin( hackysin ) * 180.0 / M_PI;
				int nside = std::ceil( search_diameter / cart_grid );
				std::cout << "search dia.    : " <<  search_diameter << std::endl;
				std::cout << "nside          : " << nside        << std::endl;
				std::cout << "resl0:           " << resl0 << std::endl;
				std::cout << "body_radius:     " << body_radius << std::endl;
				std::cout << "rif_radius:      " << rif_radius << std::endl;
				std::cout << "scaffold_radius: " << scaff_radius << std::endl;
				std::cout << "cart_grid:       " << cart_grid  << std::endl;
				std::cout << "rot_resl_deg0:   " << rot_resl_deg0 << std::endl;
				I3 nc( nside, nside, nside );
				F3 lb = target_center + F3( -cart_grid*nside/2.0, -cart_grid*nside/2.0, -cart_grid*nside/2.0 );
				F3 ub = target_center + F3(  cart_grid*nside/2.0,  cart_grid*nside/2.0,  cart_grid*nside/2.0 );
				std::cout << "cart grid ub " << ub << std::endl;
				std::cout << "cart grid lb " << lb << std::endl;
				std::cout << "(ub-lb/nc) = " << ((ub-lb)/nc.template cast<float>()) << std::endl;
				std::cout << "cartcen to corner (cart. covering radius): " << sqrt(3.0)*cart_grid/2.0 << std::endl;
				shared_ptr<DirectorOriTrans6D> director_concrete = make_shared<DirectorOriTrans6D>( rot_resl_deg0, lb, ub, nc, 1 );
				std::cout << "Director:" << endl << *director_concrete << endl;
				director = director_concrete;
				std::cout << "nest size0:    " << director->size(0) << std::endl;
				std::cout << "size of search space: ~" << float(director->size(0))*1024.0*1024.0*1024.0 << " grid points" << std::endl;
			}

			std::vector< SearchPointWithRots > packed_results;
			std::vector< ScenePtr > scene_pt( omp_max_threads_1() );
			int64_t non0_space_size = 0;
			int64_t npack = 0;
			int64_t total_search_effort = 0;
			{
		        std::chrono::time_point<std::chrono::high_resolution_clock> start_rif = std::chrono::high_resolution_clock::now();

				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				print_header( "perform hierarchical search" ); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				std::vector< std::vector< SearchPoint > > samples( RESLS.size() );
				{
					samples[0].resize( director->size(0) );
					for( uint64_t i = 0; i < director->size(0); ++i ) samples[0][i] = SearchPoint( i );
					BOOST_FOREACH( ScenePtr & s, scene_pt ) s = scene_minimal->clone_shallow();
					for( int iresl = 0; iresl < RESLS.size(); ++iresl )
					{
						cout << "HSearsh stage " << iresl+1 << " resl " << F(5,2,RESLS[iresl]) << " begin threaded sampling, " << KMGT(samples[iresl].size()) << " samples: ";
						int64_t const out_interval = samples[iresl].size()/50;
						std::exception_ptr exception = nullptr;
					    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
				        start = std::chrono::high_resolution_clock::now();
				        total_search_effort += samples[iresl].size();

						#ifdef USE_OPENMP
						#pragma omp parallel for schedule(dynamic,64)
						#endif
						for( int64_t i = 0; i < samples[iresl].size(); ++i ){
							if( exception ) continue;
							try {
								if( i%out_interval==0 ){ cout << '*'; cout.flush();	}
								uint64_t const isamp = samples[iresl][i].index;
								ScenePtr tscene( scene_pt[omp_get_thread_num()] );
								director->set_scene( isamp, iresl, *tscene );

								if( tether_to_input_position ){
									EigenXform x = tscene->position(1);
									x.translation() -= scaffold_center;
									float xmag =  xform_magnitude( x, redundancy_filter_rg );
									if( xmag > tether_to_input_position_cut + RESLS[iresl] ){
										samples[iresl][i].score = 9e9;
										continue;
									} else {
										// std::cout << "inbounds " << iresl << " " << xform_magnitude( tscene->position(1), redundancy_filter_rg ) << std::endl;
									}
								}
								samples[iresl][i].score = objectives[iresl]->score( *tscene );
							} catch( std::exception const & ex ) {
								#ifdef USE_OPENMP
								#pragma omp critical
								#endif
								exception = std::current_exception();
							}
						}
						if( exception ) std::rethrow_exception(exception);
				      	end = std::chrono::high_resolution_clock::now();
						std::chrono::duration<double> elapsed_seconds_rif = end-start;
						float rate = (double)samples[iresl].size()/ elapsed_seconds_rif.count()/omp_max_threads();
						cout << endl;// << "done threaded sampling, partitioning data..." << endl;

						SearchPoint max_pt, min_pt;
						int64_t len = samples[iresl].size();
						if( samples[iresl].size() > beam_size/DIMPOW2 ){
							__gnu_parallel::nth_element( samples[iresl].begin(), samples[iresl].begin()+beam_size/DIMPOW2, samples[iresl].end() );
							len = beam_size/DIMPOW2;
							min_pt = *__gnu_parallel::min_element( samples[iresl].begin(), samples[iresl].begin()+len );
							max_pt = *(samples[iresl].begin()+beam_size/DIMPOW2);
						} else {
							min_pt = *__gnu_parallel::min_element( samples[iresl].begin(), samples[iresl].end() );
							max_pt = *__gnu_parallel::max_element( samples[iresl].begin(), samples[iresl].end() );
						}

						cout << "HSearsh stage " << iresl+1 << " complete, resl. " << F(7,3,RESLS[iresl]) << ", "
							  << " " << KMGT(samples[iresl].size()) << ", promote: " << F(9,6,min_pt.score) << " to "
							  << F(9,6, std::min(global_score_cut,max_pt.score)) << " rate " << rate << "/s/t " << std::endl;

						if( iresl+1 == samples.size() ) break;

						for( int64_t i = 0; i < len; ++i ){
							uint64_t isamp0 = samples[iresl][i].index;
							if( samples[iresl][i].score >= global_score_cut ) continue;
							if( iresl == 0 ) ++non0_space_size;
							for( uint64_t j = 0; j < DIMPOW2; ++j ){
								uint64_t isamp = isamp0 * DIMPOW2 + j;
								samples[iresl+1].push_back( SearchPoint(isamp) );
							}
						}
						runtime_assert_msg( samples[iresl+1].size() , "search fail, no valid samples!" );
						samples[iresl].clear();

					}
					std::cout << "full sort of final samples" << std::endl;
					__gnu_parallel::sort( samples.back().begin(), samples.back().end() );
				}
				std::chrono::duration<double> elapsed_seconds_rif = std::chrono::high_resolution_clock::now()-start_rif;
				time_rif += elapsed_seconds_rif.count();

				std::cout << "total non-0 space size was approx " << float(non0_space_size)*1024.0*1024.0*1024.0 << " grid points" << std::endl;
				std::cout << "total search effort " << KMGT(total_search_effort) << std::endl;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////         HACK PACK           /////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		        std::chrono::time_point<std::chrono::high_resolution_clock> start_pack = std::chrono::high_resolution_clock::now();

				if( option[rif_dock::hack_pack]() ){

					if( use_scaffold_bounding_grids ){
						if( 0 == scene_minimal->template num_actors<SimpleAtom>(0) ){
							BOOST_FOREACH( SimpleAtom const & sa, target_simple_atoms )	scene_minimal->add_actor( 0, sa );
							runtime_assert( scene_minimal->template num_actors<SimpleAtom>(0) == target_simple_atoms.size() );
						}
					} else {
						// for final stage, use all scaffold atoms, not just CB ones
						runtime_assert( scene_minimal->clear_actors<SimpleAtom>( 1 ) );
						runtime_assert( scene_minimal->template num_actors<SimpleAtom>(1) == 0 );
						BOOST_FOREACH( SimpleAtom const & sa, scaffold_simple_atoms_all ) scene_minimal->add_actor( 1, sa );
						runtime_assert( scene_minimal->template num_actors<SimpleAtom>(1) == scaffold_simple_atoms_all.size() );
						// these should be shallow copies in scene_pt
						// so editing scene_minimal will change all conformations
						runtime_assert( scene_pt.front()->template num_actors<SimpleAtom>(1) == scaffold_simple_atoms_all.size() );
					}

					// if( scene_minimal->template num_actors<SimpleAtom>(0) == 0 ){
					// 	for(int i = 0; i < 10; ++i) std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!! hackpack add sterics back !!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
					// 	BOOST_FOREACH( SimpleAtom const & sa, target_simple_atoms )	scene_minimal->add_actor( 0, sa );
					// }
					// runtime_assert( scene_minimal->template num_actors<SimpleAtom>(0) == target_simple_atoms.size() );

				    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
			        start = std::chrono::high_resolution_clock::now();

					size_t n_packsamp = 0;
					for( n_packsamp; n_packsamp < samples.back().size(); ++n_packsamp ){
						if( samples.back()[n_packsamp].score > 0 ) break;
					}
					int const config = RESLS.size()-1;
					npack = std::min( n_packsamp, (size_t)(total_search_effort * ( hack_pack_frac / hackpackopts.pack_iter_mult) ) );
					packed_results.resize( npack );
					print_header( "hack-packing top " + KMGT(npack) );
					std::cout << "packing options: " << hackpackopts << std::endl;
					std::cout << "packing w/rif rofts ";
					int64_t const out_interval = std::max<int64_t>(1,npack/100);
					std::exception_ptr exception = nullptr;
					#ifdef USE_OPENMP
					#pragma omp parallel for schedule(dynamic,64)
					#endif
					for( int ipack = 0; ipack < npack; ++ipack ){
						if( exception ) continue;
						try {
							if( ipack%out_interval==0 ){ cout << '*'; cout.flush();	}
							uint64_t const isamp = samples.back()[ipack].index;
							if( samples.back()[ipack].score > global_score_cut ) continue;
							packed_results[ ipack ].index = isamp;
							packed_results[ ipack ].prepack_rank = ipack;
							ScenePtr tscene( scene_pt[omp_get_thread_num()] );
							director->set_scene( isamp, RESLS.size()-1, *tscene );
							packed_results[ ipack ].score = packing_objective->score_with_rotamers( *tscene, packed_results[ ipack ].rotamers() );
						} catch( std::exception const & ex ) {
							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							exception = std::current_exception();
						}
					}
					if( exception ) std::rethrow_exception(exception);
			        end = std::chrono::high_resolution_clock::now();

					std::cout << std::endl;
					std::cout << "full sort of packed samples" << std::endl;
					__gnu_parallel::sort( packed_results.begin(), packed_results.end() );

					std::chrono::duration<double> elapsed_seconds_pack = end-start;
					std::cout << "packing rate: " << (double)npack/elapsed_seconds_pack.count()                   << " iface packs per second" << std::endl;
					std::cout << "packing rate: " << (double)npack/elapsed_seconds_pack.count()/omp_max_threads() << " iface packs per second per thread" << std::endl;

				} else {
					packed_results.resize( samples.back().size() );
					#ifdef USE_OPENMP
					#pragma omp parallel for schedule(dynamic,1024)
					#endif
					for( int ipack = 0; ipack < packed_results.size(); ++ipack ){
						packed_results[ipack].score = 0;//samples.back()[ipack].score;
						packed_results[ipack].index = samples.back()[ipack].index;
					}
				}

				std::chrono::duration<double> elapsed_seconds_pack = std::chrono::high_resolution_clock::now()-start_pack;
				time_pck += elapsed_seconds_pack.count();
			}






			bool const do_rosetta_score = rosetta_score_fraction > 0 || rosetta_score_then_min_below_thresh > -9e8 || rosetta_score_at_least > 0;

			if( do_rosetta_score && option[rif_dock::hack_pack]() ){

				std::chrono::time_point<std::chrono::high_resolution_clock> start_rosetta = std::chrono::high_resolution_clock::now();

				int n_score_calculations = 0;

				for( int do_min = 0; do_min < 2; ++do_min ){

					std::chrono::duration<double> time_copy, time_score, time_min;

					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					if( do_min ) print_header( "rosetta min and score" ); ////////////////////////////////////////////////////////////////////
					else         print_header( "rosetta score" ); ////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					runtime_assert_msg( target_res.size() == 1, "rosetta_minim_below_thresh is intended for use with small molecules" );

					std::vector<protocols::simple_moves::MinMoverOP> minmover_pt(omp_max_threads());
					std::vector<core::scoring::ScoreFunctionOP> scorefunc_pt(omp_max_threads());
					std::vector<core::pose::Pose> both_full_per_thread(omp_max_threads());
					std::vector<core::pose::Pose> both_per_thread     (omp_max_threads());
					std::vector<core::pose::Pose> target_pt           (omp_max_threads());
					std::vector<core::pose::Pose> work_pose_pt        (omp_max_threads());					
					for( int i = 0; i < omp_max_threads(); ++i){
						both_full_per_thread[i] = both_full_pose;
						both_per_thread[i] = both_pose;
						scorefunc_pt[i] = core::scoring::get_score_function();
						scorefunc_pt[i]->set_etable( "FA_STANDARD_SOFT" );
						if( !do_min ){
							scorefunc_pt[i]->set_weight( core::scoring::fa_dun, scorefunc_pt[i]->get_weight(core::scoring::fa_dun)*0.8 );
							scorefunc_pt[i]->set_weight( core::scoring::fa_rep, scorefunc_pt[i]->get_weight(core::scoring::fa_rep)*0.8 );
						} else {
							scorefunc_pt[i]->set_weight( core::scoring::fa_dun, scorefunc_pt[i]->get_weight(core::scoring::fa_dun)*1.2 );							
							scorefunc_pt[i]->set_weight( core::scoring::fa_rep, scorefunc_pt[i]->get_weight(core::scoring::fa_rep)*1.2 );
						}
						scorefunc_pt[i]->set_weight( core::scoring::hbond_sc, scorefunc_pt[i]->get_weight(core::scoring::hbond_sc)*1.3 );
						// scorefunc_pt[i]->set_weight( core::scoring::fa_rep, scorefunc_pt[i]->get_weight(core::scoring::fa_rep)*0.67 );
						// scorefunc_pt[i]->set_weight( core::scoring::fa_dun, scorefunc_pt[i]->get_weight(core::scoring::fa_dun)*0.67 );
						// core::scoring::methods::EnergyMethodOptions opts = scorefunc_pt[i]->energy_method_options();
						// core::scoring::hbonds::HBondOptions hopts = opts.hbond_options();
						// hopts.use_hb_env_dep( false );
						// opts.hbond_options( hopts );
						// scorefunc_pt[i]->set_energy_method_options( opts );
						core::kinematics::MoveMapOP movemap = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );
						movemap->set_chi(true);
						movemap->set_jump(true);
						movemap->set_bb(false);
						minmover_pt[i] = protocols::simple_moves::MinMoverOP(
							new protocols::simple_moves::MinMover( movemap, scorefunc_pt[i], "dfpmin_armijo_nonmonotone", 0.001, true ) );
					}

				    std::chrono::time_point<std::chrono::high_resolution_clock> startall = std::chrono::high_resolution_clock::now();

					size_t n_scormin = 0;
					if( do_min ){
						// min take ~10x score time, so do on 1/10th of the scored
						n_scormin = n_score_calculations * rosetta_min_fraction;
					} else {
						// for scoring, use user cut
						n_scormin = rosetta_score_fraction/40.0 * total_search_effort;
						if( rosetta_score_then_min_below_thresh > -9e8 ){
							for( n_scormin=0; n_scormin < packed_results.size(); ++n_scormin ){
								if( packed_results[n_scormin].score > rosetta_score_then_min_below_thresh )
									break;
							}
						}
						n_scormin = std::min<int>( std::max<int>( n_scormin, rosetta_score_at_least ), rosetta_score_at_most );
						n_scormin = std::min<int>( n_scormin, packed_results.size() );
						n_score_calculations = n_scormin;
					}
					packed_results.resize(n_scormin);
					int64_t const out_interval = std::max<int64_t>(1,n_scormin/50);
					if( do_min) std::cout << "rosetta min on "   << KMGT(n_scormin) << ": ";
					else        std::cout << "rosetta score on " << KMGT(n_scormin) << ": ";					
					#ifdef USE_OPENMP
					#pragma omp parallel for schedule(dynamic,1)
					#endif
					for( int imin = 0; imin < n_scormin; ++imin )
					{
						if( imin%out_interval==0 ){ cout << '*'; cout.flush();	}

						int const ithread = omp_get_thread_num();

						director->set_scene( packed_results[imin].index, RESLS.size()-1, *scene_pt[ithread] );
						EigenXform xposition1 = scene_pt[ithread]->position(1);
						EigenXform xalignout = EigenXform::Identity();
						if( align_to_scaffold ) xalignout = xposition1.inverse();

					    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
				        start = std::chrono::high_resolution_clock::now();

						// core::pose::PoseOP pose_to_min_ptr = core::pose::PoseOP(new core::pose::Pose);
						core::pose::Pose & pose_to_min( work_pose_pt[ithread] );

						if( option[rif_dock::full_scaffold_output]() ) pose_to_min = both_full_per_thread[ithread];
						else                                           pose_to_min = both_per_thread[ithread];
						xform_pose( pose_to_min, eigen2xyz(xalignout)            , pose_to_min.n_residue(), pose_to_min.n_residue()   );
						xform_pose( pose_to_min, eigen2xyz(xalignout*xposition1) , 1                      , pose_to_min.n_residue()-1 );


						// place the rotamers
						core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");

						for( int ipr = 0; ipr < packed_results[imin].numrots(); ++ipr ){
							int ires = scaffres_l2g.at( packed_results[imin].rotamers().at(ipr).first );
							int irot =                  packed_results[imin].rotamers().at(ipr).second;
							core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rot_index.resname(irot)) );
							pose_to_min.replace_residue( ires+1, *newrsd, true );
							for( int ichi = 0; ichi < rot_index.nchi(irot); ++ichi ){
								pose_to_min.set_chi( ichi+1, ires+1, rot_index.chi( irot, ichi ) );
							}
						}
						end = std::chrono::high_resolution_clock::now();
						std::chrono::duration<double> elapsed_seconds_copypose = end-start;
						#pragma omp critical
						time_copy += elapsed_seconds_copypose;
						// pose_to_min.dump_pdb("test_pre.pdb");

						if( do_min ){
							// std::cout << "MIN!" << std::endl;
					        start = std::chrono::high_resolution_clock::now();
							minmover_pt[ithread]->apply( pose_to_min );
					        end = std::chrono::high_resolution_clock::now();
							std::chrono::duration<double> elapsed_seconds_min = end-start;
							#pragma omp critical
							time_min += elapsed_seconds_min;
						} else {
							// std::cout << "SCORE!" << std::endl;							
					        start = std::chrono::high_resolution_clock::now();
							scorefunc_pt[ithread]->score( pose_to_min );
					        end = std::chrono::high_resolution_clock::now();
							std::chrono::duration<double> elapsed_seconds_score = end-start;
							#pragma omp critical
							time_score += elapsed_seconds_score;
						}

						float total_lj_neg  = scorefunc_pt[ithread]->get_weight(core::scoring::fa_atr) * pose_to_min.energies().total_energies()[core::scoring::fa_atr];
						      total_lj_neg += scorefunc_pt[ithread]->get_weight(core::scoring::fa_rep) * pose_to_min.energies().total_energies()[core::scoring::fa_rep];
						      total_lj_neg = std::max(0.0f,total_lj_neg);
						packed_results[imin].score  = 1.00*total_lj_neg;
						packed_results[imin].score += 1.00*pose_to_min.energies().total_energies()[core::scoring::hbond_sc]; // last res is ligand
						packed_results[imin].score += 1.00*pose_to_min.energies().residue_total_energy(pose_to_min.n_residue());
						for( int ipr = 0; ipr < packed_results[imin].numrots(); ++ipr ){
							int ires = scaffres_l2g.at( packed_results[imin].rotamers().at(ipr).first );
							packed_results[imin].score += 0.5*pose_to_min.energies().residue_total_energy(ires);
						}


						if( do_min && packed_results[imin].score < rosetta_score_cut ){
							packed_results[imin].pose_ = core::pose::PoseOP( new core::pose::Pose(pose_to_min) );
						}

						// #pragma omp critical
						// pose_to_min.energies().show(cout,pose_to_min.n_residue());
						// pose_to_min.dump_pdb("test_post.pdb");
						// #pragma omp critical
						// std::cout << imin << " prescore: " << e_pre_min << " minscore: " << pose_to_min.energies().total_energy() << std::endl;
						// utility_exit_with_message("testing");
					}
					cout << endl;
					__gnu_parallel::sort( packed_results.begin(), packed_results.end() );
					{
						size_t n_scormin = 0;
						for( n_scormin; n_scormin < packed_results.size(); ++n_scormin ){
							if( do_min && packed_results[n_scormin].score > rosetta_score_cut ) break;
						}
						packed_results.resize(n_scormin);
					}

				    std::chrono::time_point<std::chrono::high_resolution_clock> stopall = std::chrono::high_resolution_clock::now();
					std::chrono::duration<double> elap_sec = stopall - startall;

					if( do_min ){
						std::cout << "total pose copy time: " << time_copy.count() << "s total min time: " << time_min.count() << "s" << " walltime total " << elap_sec.count() << std::endl;
						std::cout << "min rate: "   << (double)n_scormin / elap_sec.count()                     << " sc/rb minimizations per second" << std::endl;
						std::cout << "min rate: "   << (double)n_scormin / elap_sec.count() / omp_max_threads() << " sc/rb minimizations per second per thread" << std::endl;
					} else {
						std::cout << "total pose copy time: " << time_copy.count() << "s total score time: " << time_score.count() << "s" << " walltime total " << elap_sec.count() << std::endl;						
						std::cout << "score rate: " << (double)n_scormin / elap_sec.count()                     << " rosetta scores per second" << std::endl;
						std::cout << "score rate: " << (double)n_scormin / elap_sec.count() / omp_max_threads() << " rosetta scores per second per thread" << std::endl;						
					}
				}

				std::chrono::duration<double> elapsed_seconds_rosetta = std::chrono::high_resolution_clock::now()-start_rosetta;
				time_ros += elapsed_seconds_rosetta.count();

			}






			/*
			 this needs to get fixed
			 the job of this code:
			 (1) do redundancy filtering
			 (2) build selected_results, allresults from packed_results
			 could probably split these up?
			       packed_results is
			       		struct SearchPointWithRots {
							float score;
							uint32_t prepack_rank;
							uint64_t index;
							shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
							core::pose::PoseOP pose_ = nullptr;
					allresults is
						struct RifDockResult {
							float dist0, packscore, nopackscore, rifscore, stericscore;
							uint64_t isamp, scene_index;
							uint32_t prepack_rank;
							float cluster_score;
							bool operator< ( RifDockResult const & o ) const { return packscore < o.packscore; }
							shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
							core::pose::PoseOP pose_ = nullptr;
							size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
							std::vector< std::pair<intRot,intRot> > const & rotamers() const { assert(rotamers_!=nullptr); return *rotamers_; }


			*/
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "compile and filter results" ); ///////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::vector< RifDockResult > selected_results, allresults;
			{
				int64_t Nout = packed_results.size(); //samples.back().size(); //std::min(samples.back().size(),(size_t)1000000);
				Nout = std::min( (int64_t)option[rif_dock::n_result_limit](), Nout );

				std::vector< std::vector< RifDockResult > > allresults_pt( omp_max_threads() );
				std::vector< std::pair< EigenXform, int64_t > > selected_xforms;
				selected_xforms.reserve(65536); // init big to reduce liklihood of resizes
				float redundancy_filter_mag = option[ rif_dock::redundancy_filter_mag ]();
				int nclose = 0;
				int nclosemax      = option[rif_dock::force_output_if_close_to_input_num]();
				float nclosethresh = option[rif_dock::force_output_if_close_to_input]();
				int n_pdb_out = option[rif_dock::n_pdb_out]();
				std::cout << "redundancy_filter_mag " << redundancy_filter_mag << "A \"rmsd\"" << std::endl;
				int64_t Nout_singlethread = std::min( (int64_t)10000, Nout );

				std::cout << "going throuth 10K results (1 thread): ";
				int64_t out_interval = 10000/81;
				for( int64_t isamp = 0; isamp < Nout_singlethread; ++isamp ){
					if( isamp%out_interval==0 ){ cout << '*'; cout.flush();	}
					awful_compile_output_helper< EigenXform, ScenePtr, ObjectivePtr >(
						isamp, RESLS.size()-1, packed_results, scene_pt, director,
						redundancy_filter_rg, redundancy_filter_mag, scaffold_center,
						allresults_pt, selected_results, selected_xforms, n_pdb_out,
						#ifdef USE_OPENMP
							dump_lock,
						#endif
						objectives.back(), nclose, nclosemax, nclosethresh
					);
				}
				std::cout << std::endl;

				std::cout << "going throuth all results (threaded): ";
				out_interval = Nout / 82;
				std::exception_ptr exception = nullptr;
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,8)
				#endif
				for( int64_t isamp = Nout_singlethread; isamp < Nout; ++isamp ){
					if( exception ) continue;
					try{
						if( isamp%out_interval==0 ){ cout << '*'; cout.flush();	}
						awful_compile_output_helper< EigenXform, ScenePtr, ObjectivePtr >(
							isamp, RESLS.size()-1, packed_results, scene_pt, director,
							redundancy_filter_rg, redundancy_filter_mag, scaffold_center,
							allresults_pt, selected_results, selected_xforms, n_pdb_out,
							#ifdef USE_OPENMP
								dump_lock,
							#endif
							objectives.back(), nclose, nclosemax, nclosethresh
						);
					} catch(...) {
						#pragma omp critical
						exception = std::current_exception();
					}
				}
				if( exception ) std::rethrow_exception(exception);
				std::cout << std::endl;

				std::cout << "sort compiled results" << std::endl;
				BOOST_FOREACH( std::vector<RifDockResult> const & rs, allresults_pt ){
					BOOST_FOREACH( RifDockResult const & r, rs ){
						allresults.push_back( r );
					}
				}
				__gnu_parallel::sort( allresults.begin(), allresults.end() );


			} // end result compilation loop


			std::cout << "allresults.size(): " << allresults.size() << " selected_results.size(): " << selected_results.size() << std::endl;

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "timing info" ); //////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			std::cout<<"total RIF     time: "<<KMGT(time_rif)<<" fraction: "<<time_rif/(time_rif+time_pck+time_ros)<<std::endl;
			std::cout<<"total Pack    time: "<<KMGT(time_pck)<<" fraction: "<<time_pck/(time_rif+time_pck+time_ros)<<std::endl;
			std::cout<<"total Rosetta time: "<<KMGT(time_ros)<<" fraction: "<<time_ros/(time_rif+time_pck+time_ros)<<std::endl;			

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "output results" ); //////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			BOOST_FOREACH( RifDockResult const & r, allresults ){
				; // nothing with all results ATM
			}



			if( align_to_scaffold )	std::cout << "ALIGN TO SCAFFOLD" << std::endl;
			else                    std::cout << "ALIGN TO TARGET"   << std::endl;
			for( int i_selected_result = 0; i_selected_result < selected_results.size(); ++i_selected_result ){
				RifDockResult const & selected_result = selected_results.at( i_selected_result );

				std::string pdboutfile = outdir + "/" + scafftag + "_" + devel::scheme::str(i_selected_result,9)+".pdb";
				if( output_tag.size() ){
					pdboutfile = outdir + "/" + scafftag+"_" + output_tag + "_" + devel::scheme::str(i_selected_result,9)+".pdb";
				}

				std::ostringstream oss;
		        oss << "rif score: " << I(4,i_selected_result)
		            << " rank "       << I(9,selected_result.isamp)
		            << " dist0:    "  << F(7,4,selected_result.dist0)
		            << " packscore: " << F(7,3,selected_result.packscore)
		            << " score: "     << F(7,3,selected_result.nopackscore)
		            << " rif: "       << F(7,3,selected_result.rifscore)
		            << " steric: "    << F(7,3,selected_result.stericscore)
		            << " cluster: "   << I(7,selected_result.cluster_score)
		            << " rifrank: "   << I(7,selected_result.prepack_rank) << " " << F(7,5,(float)selected_result.prepack_rank/(float)npack)
		            << " " << pdboutfile
		            << std::endl;
		        std::cout << oss.str();
		        dokout << oss.str(); dokout.flush();

				 // crappy pdb io
		        {

					director->set_scene( selected_result.scene_index, RESLS.size()-1, *scene_full    );
					director->set_scene( selected_result.scene_index, RESLS.size()-1, *scene_minimal );

					EigenXform xposition1 = scene_full->position(1);
					EigenXform xalignout = EigenXform::Identity();
					if( align_to_scaffold ){
						xalignout = xposition1.inverse();
					}

					std::ostringstream packout, allout;
					std::map< int, std::string > pikaa;
					for( int i_actor = 0; i_actor < scene_minimal->template num_actors<BBActor>(1); ++i_actor ){
						BBActor bba = scene_minimal->template get_actor<BBActor>(1,i_actor);
						int const ires = scaffres_l2g.at( bba.index_ );

						// if( dump_all_rif_rots )
						{
							std::vector< std::pair< float, int > > rotscores;
							rif_ptrs.back()->get_rotamers_for_xform( bba.position(), rotscores );
							typedef std::pair<float,int> PairFI;
							BOOST_FOREACH( PairFI const & p, rotscores ){
								int const irot = p.second;
								float const sc = p.first + scaffold_onebody_glob0.at( ires ).at( irot );
								if( sc < 0 ){
									allout << "MODEL" << endl;
									BOOST_FOREACH( SchemeAtom a, rot_index.rotamers_.at( irot ).atoms_ ){
										a.set_position( xalignout * bba.position() * a.position() ); // is copy
										a.nonconst_data().resnum = ires;
										::scheme::actor::write_pdb( allout, a, nullptr );
									}
									allout << "ENDMDL" << endl;
									char oneletter = rot_index_p->oneletter(irot);
									if( std::find( pikaa[ires+1].begin(), pikaa[ires+1].end(), oneletter ) == pikaa[ires+1].end() ){
										pikaa[ires+1] += oneletter;
									}
								}

							}
						}

						int packed_rot = -1;
						for( int ipr = 0; ipr < selected_result.numrots(); ++ipr ){
							// std::cout << "checking rots " << sp.rotamers()[ipr].first << " " << scaffres_g2l[ires] << std::endl;
							if( selected_result.rotamers().at(ipr).first == scaffres_g2l.at( ires ) ){
								packed_rot = selected_result.rotamers().at(ipr).second;
							}
						}
						if( packed_rot >= 0 ){
							// packout << "MODEL" << endl;
							BOOST_FOREACH( SchemeAtom a, rot_index.rotamers_.at( packed_rot ).atoms_ ){
								a.set_position( xalignout * bba.position() * a.position() ); // is copy
								a.nonconst_data().resnum = ires;
								::scheme::actor::write_pdb( packout, a, nullptr );
							}
							packout << "TER" << endl;
						}

					}

					core::pose::Pose pose_from_rif;
					if( option[rif_dock::full_scaffold_output]() ) pose_from_rif = both_full_pose;
					else                                           pose_from_rif = both_pose;
					xform_pose( pose_from_rif, eigen2xyz(xalignout)           , scaffold.n_residue()+1, pose_from_rif.n_residue() );
					xform_pose( pose_from_rif, eigen2xyz(xalignout*xposition1),                      1,     scaffold.n_residue() );

					// place the rotamers
					core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
					std::ostringstream resfile, expdb;
					resfile << "ALLAA" << std::endl;
					resfile << "start" << std::endl;
					expdb << "rif_residues ";

					for( int ipr = 0; ipr < selected_result.numrots(); ++ipr ){
						int ires = scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
						int irot =                  selected_result.rotamers().at(ipr).second;
						core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rot_index.resname(irot)) );
						pose_from_rif.replace_residue( ires+1, *newrsd, true );
						resfile << ires+1 << " A NATRO" << std::endl;
						expdb << ires+1 << (ipr+1<selected_result.numrots()?",":""); // skip comma on last one
						for( int ichi = 0; ichi < rot_index.nchi(irot); ++ichi ){
							pose_from_rif.set_chi( ichi+1, ires+1, rot_index.chi( irot, ichi ) );
						}
					}

					core::pose::Pose & pose_to_dump( *(selected_result.pose_ ? selected_result.pose_.get() : &pose_from_rif) );
					utility::io::ozstream out1( pdboutfile );
					// scene_full->set_position( 1, xalignout * xposition1 );
					// write_pdb( out1, dynamic_cast<Scene&>(*scene_full), rot_index.chem_index_ );
					if( pdb_info_pikaa ){
						for( auto p : pikaa ){
							std::sort( p.second.begin(), p.second.end() );
							pose_to_dump.pdb_info()->add_reslabel(p.first, "PIKAA" );
							pose_to_dump.pdb_info()->add_reslabel(p.first, p.second );
						}
					} else {
						for( auto p : pikaa ){
							pose_to_dump.pdb_info()->add_reslabel(p.first, "RIFRES" );
						}
					}
					// if( selected_result.pose_ ){
					// 	for( auto p : pikaa ){
					// 		std::cout << "residue " << p.first << " " << selected_result.pose_->residue(p.first).name() << " fa_rep: "
					// 		          << selected_result.pose_->energies().residue_total_energies(p.first)[core::scoring::fa_rep] << std::endl;
					// 	}
					// }

					out1 << expdb.str() << std::endl;
					pose_to_dump.dump_pdb(out1);

					out1.close();

					if( option[rif_dock::dump_resfile]() ){
						utility::io::ozstream out1res( outdir + "/" + scafftag+"_"+devel::scheme::str(i_selected_result,9)+".resfile");
						out1res << resfile.str();
						out1res.close();
					}

					if( dump_all_rif_rots ){
						// utility_exit_with_message("this is not currently implemented, ask Will");
						utility::io::ozstream out2( outdir + "/" + scafftag+"_allrifrots_"+devel::scheme::str(i_selected_result,9)+".pdb");
						out2 << allout.str();
						out2.close();
					}

					// utility::io::ozstream out4(scafftag+"_pack_rot_"+devel::scheme::str(i_selected_result,9)+".pdb");
					// out4 << packout.str();
					// out4.close();

				} // end crappy pdb io

			}

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
