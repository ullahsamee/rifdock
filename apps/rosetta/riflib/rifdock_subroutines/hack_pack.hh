// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_hack_pack_hh
#define INCLUDED_riflib_rifdock_subroutines_hack_pack_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/rifdock_subroutines/meta.hh>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;

template<class DirectorBase, class ScaffoldProvider >
struct HackPackData {
    RifDockOpt & opt;
    std::vector<float> & RESLS;
    DirectorBase & director;
    int64_t & total_search_effort;
    std::vector< devel::scheme::ScenePtr > & scene_pt;
    devel::scheme::ScenePtr & scene_minimal;
    std::vector<devel::scheme::SimpleAtom> & target_simple_atoms;
    // std::vector<devel::scheme::SimpleAtom > & scaffold_simple_atoms_all;
    int64_t & npack;
    ::scheme::search::HackPackOpts & packopts;
    devel::scheme::ObjectivePtr & packing_objective;
    shared_ptr< std::vector< _SearchPointWithRots<DirectorBase> > > & hsearch_results_p;
    shared_ptr<ScaffoldProvider> scaffold_provider;
};



template<class DirectorBase, class ScaffoldProvider >
void
hack_pack(
    shared_ptr<std::vector< _SearchPointWithRots<DirectorBase> > > & packed_results_p,
    HackPackData<DirectorBase, ScaffoldProvider> & d) {


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

    typedef ::scheme::util::SimpleArray<3,float> F3;
    typedef ::scheme::util::SimpleArray<3,int> I3;

    typedef _SearchPointWithRots<DirectorBase> SearchPointWithRots;
    typedef shared_ptr<ScaffoldProvider> ScaffoldProviderOP;

    typedef _DirectorBigIndex<DirectorBase> DirectorIndex;
    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;


    if( d.opt.hack_pack ){

        std::vector< SearchPointWithRots > & hsearch_results = *d.hsearch_results_p;
        packed_results_p = make_shared<std::vector< SearchPointWithRots >>();
        std::vector< SearchPointWithRots > & packed_results = *packed_results_p;

        // if( d.opt.use_scaffold_bounding_grids ){
        //     if( 0 == d.scene_minimal->template num_actors<SimpleAtom>(0) ){
        //         // Brian - If this gets uncommented, it's probably wrong. This was written assuming
        //         //  that the scene_pts are all clones of each other
        //         BOOST_FOREACH( SimpleAtom const & sa, d.target_simple_atoms ) d.scene_minimal->add_actor( 0, sa );
        //         runtime_assert( d.scene_minimal->template num_actors<SimpleAtom>(0) == d.target_simple_atoms.size() );
        //     }
        // } else {
            // for final stage, use all scaffold atoms, not just CB ones

        d.scaffold_provider->set_fa_mode(true);


        /// delete this later


        shared_ptr<ParametricScene> scene_minimal_typed( std::dynamic_pointer_cast<ParametricScene>(d.scene_minimal));
        scene_minimal_typed->replace_body(1, d.scaffold_provider->get_scaffold(scaffold_index_default_value( ScaffoldIndex())));

        for ( ScenePtr scene : d.scene_pt ) {
            shared_ptr<ParametricScene> scene_typed( std::dynamic_pointer_cast<ParametricScene>(scene));
            scene_typed->replace_body(1, d.scaffold_provider->get_scaffold(scaffold_index_default_value( ScaffoldIndex())));
        }


        // delete ///////////////////

            // runtime_assert( d.scene_minimal->template clear_actors<SimpleAtom>( 1 ) );
            // runtime_assert( d.scene_minimal->template num_actors<SimpleAtom>(1) == 0 );
            // BOOST_FOREACH( SimpleAtom const & sa, d.scaffold_simple_atoms_all ) d.scene_minimal->add_actor( 1, sa );
            // runtime_assert( d.scene_minimal->template num_actors<SimpleAtom>(1) == d.scaffold_simple_atoms_all.size() );
            // these should be shallow copies in d.scene_pt
            // so editing d.scene_minimal will change all conformations
            // runtime_assert( d.scene_pt.front()->template num_actors<SimpleAtom>(1) == d.scaffold_simple_atoms_all.size() );
        // }






        std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
        start = std::chrono::high_resolution_clock::now();

        size_t n_packsamp = 0;
        for( n_packsamp; n_packsamp < hsearch_results.size(); ++n_packsamp ){
            if( hsearch_results[n_packsamp].score > 0 ) break;
        }
        int const config = d.RESLS.size()-1;
        d.npack = std::min( n_packsamp, (size_t)(d.total_search_effort *
            ( d.opt.hack_pack_frac / (d.packopts.pack_n_iters*d.packopts.pack_iter_mult)) ) );

        packed_results.resize( d.npack );
        print_header( "hack-packing top " + KMGT(d.npack) );
        std::cout << "packing options: " << d.packopts << std::endl;
        std::cout << "packing w/rif rofts ";
        int64_t const out_interval = std::max<int64_t>(1,d.npack/100);
        std::exception_ptr exception = nullptr;
        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic,64)
        #endif
        for( int ipack = 0; ipack < d.npack; ++ipack ){
            if( exception ) continue;
            try {
                if( ipack%out_interval==0 ){ cout << '*'; cout.flush(); }
                DirectorIndex const isamp = hsearch_results[ipack].index;
                if( hsearch_results[ipack].score > d.opt.global_score_cut ) continue;
                packed_results[ ipack ].index = isamp;
                packed_results[ ipack ].prepack_rank = ipack;
                ScenePtr tscene( d.scene_pt[omp_get_thread_num()] );
                d.director->set_scene( isamp, d.RESLS.size()-1, *tscene );
                packed_results[ ipack ].score = d.packing_objective->score_with_rotamers( *tscene, packed_results[ ipack ].rotamers() );
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
        std::cout << "packing rate: " << (double)d.npack/elapsed_seconds_pack.count()                   << " iface packs per second" << std::endl;
        std::cout << "packing rate: " << (double)d.npack/elapsed_seconds_pack.count()/omp_max_threads() << " iface packs per second per thread" << std::endl;



    } else {
        packed_results_p = d.hsearch_results_p;

    }












}










#endif