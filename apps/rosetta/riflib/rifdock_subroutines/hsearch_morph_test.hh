// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_hsearch_morph_test_hh
#define INCLUDED_riflib_rifdock_subroutines_hsearch_morph_test_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/rifdock_subroutines/meta.hh>
#include <riflib/rifdock_subroutines/hsearch_original.hh>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


// template<__DirectorBase>
// using HsearchFunctionType = typedef


template<class DirectorBase, class ScaffoldProvider >
bool
hsearch_morph_test(
    shared_ptr<std::vector< _SearchPointWithRots<DirectorBase> > > & hsearch_results_p,
    HsearchData<DirectorBase, ScaffoldProvider > & d) {


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

    typedef _DirectorBigIndex<DirectorBase> DirectorIndex;
    typedef tmplSearchPoint<DirectorIndex> SearchPoint;

    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;

//////////////////////////////////////////////////////////////////////////////////////


            // these are in rosetta numbers
            // int rmsd_resid1 = 107;
            // int rmsd_resid2 = 54;
            // int rmsd_resid3 = 13;


            // utility::vector1<core::Size> rmsd_resids { rmsd_resid1, rmsd_resid2, rmsd_resid3 };
            // std::vector<SchemeAtom> rmsd_atoms;
            // devel::scheme::get_scheme_atoms_cbonly( scaffold_unmodified_from_file, rmsd_resids, rmsd_atoms);
            
            // SchemeAtom rmsd_cb1 = rmsd_atoms[0];
            // SchemeAtom rmsd_cb2 = rmsd_atoms[1];
            // SchemeAtom rmsd_cb3 = rmsd_atoms[2];

            // std::cout << "Checking that we have the right CBetas for rmsd" << std::endl;
            // std::cout << I(5, rmsd_resid1) << " " << rmsd_cb1.position().transpose() << std::endl;
            // std::cout << I(5, rmsd_resid2) << " " << rmsd_cb2.position().transpose() << std::endl;
            // std::cout << I(5, rmsd_resid3) << " " << rmsd_cb3.position().transpose() << std::endl;


//////////////////////////////////////////////////////////////////////////////////////

using ::scheme::scaffold::BOGUS_INDEX;
using ::scheme::scaffold::TreeIndex;
using ::scheme::scaffold::TreeLimits;

    shared_ptr<MorphingScaffoldProvider> morph_provider = std::dynamic_pointer_cast<MorphingScaffoldProvider>(d.scaffold_provider);

    morph_provider->test_make_children( TreeIndex(0, 0) );

    TreeLimits limits = morph_provider->get_scaffold_index_limits();
    uint64_t num_scaffolds = limits[1];

    for ( uint64_t scaffno = 0; scaffno < num_scaffolds; scaffno++ ) {
        TreeIndex ti(1, scaffno);
        ScaffoldDataCacheOP sdc = morph_provider->get_data_cache_slow(ti);
        sdc->setup_onebody_tables( d.rot_index_p, d.opt);
    }



    std::vector< std::vector< SearchPoint > > samples( d.RESLS.size() );

    bool search_failed = false;
    {
        samples[0].resize( ::scheme::kinematics::bigindex_nest_part(d.director->size(0))*num_scaffolds );
        uint64_t index_count = 0;
        for ( uint64_t scaffno = 0; scaffno < num_scaffolds; scaffno++ ) {
            for( uint64_t i = 0; i < ::scheme::kinematics::bigindex_nest_part(d.director->size(0)); ++i ) {
                samples[0][index_count++] = SearchPoint( DirectorIndex( i, TreeIndex(1, scaffno)) );
            }
        }
        BOOST_FOREACH( ScenePtr & s, d.scene_pt ) s = d.scene_minimal->clone_specific_deep(std::vector<uint64_t> {1});
        for( int iresl = 0; iresl < d.RESLS.size(); ++iresl )
        {
            cout << "HSearsh stage " << iresl+1 << " resl " << F(5,2,d.RESLS[iresl]) << " begin threaded sampling, " << KMGT(samples[iresl].size()) << " samples: ";
            int64_t const out_interval = samples[iresl].size()/50;
            std::exception_ptr exception = nullptr;
            std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
            start = std::chrono::high_resolution_clock::now();
            d.total_search_effort += samples[iresl].size();

            std::vector< double > rmsds( samples[iresl].size(), 1000 );

            bool answer_exists = false;

            #ifdef USE_OPENMP
            #pragma omp parallel for schedule(dynamic,64)
            #endif
            for( int64_t i = 0; i < samples[iresl].size(); ++i ){
                if( exception ) continue;
                try {
                    if( i%out_interval==0 ){ cout << '*'; cout.flush(); }
                    DirectorIndex const isamp = samples[iresl][i].index;

                    ScenePtr tscene( d.scene_pt[omp_get_thread_num()] );
                    d.director->set_scene( isamp, iresl, *tscene );




                    // the real rif score!!!!!!
                    samples[iresl][i].score = d.objectives[iresl]->score( *tscene );// + tot_sym_score;





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
            if( samples[iresl].size() > d.opt.beam_size/d.opt.DIMPOW2 ){
                __gnu_parallel::nth_element( samples[iresl].begin(), samples[iresl].begin()+d.opt.beam_size/d.opt.DIMPOW2, samples[iresl].end() );
                len = d.opt.beam_size/d.opt.DIMPOW2;
                min_pt = *__gnu_parallel::min_element( samples[iresl].begin(), samples[iresl].begin()+len );
                max_pt = *(samples[iresl].begin()+d.opt.beam_size/d.opt.DIMPOW2);
            } else {
                min_pt = *__gnu_parallel::min_element( samples[iresl].begin(), samples[iresl].end() );
                max_pt = *__gnu_parallel::max_element( samples[iresl].begin(), samples[iresl].end() );
            }

            cout << "HSearsh stage " << iresl+1 << " complete, resl. " << F(7,3,d.RESLS[iresl]) << ", "
                  << " " << KMGT(samples[iresl].size()) << ", promote: " << F(9,6,min_pt.score) << " to "
                  << F(9,6, std::min(d.opt.global_score_cut,max_pt.score)) << " rate " << KMGT(rate) << "/s/t " << std::endl;

            // cout << "Answer: " << ( answer_exists ? "exists" : "doesn't exist" ) << std::endl;

            if( iresl+1 == samples.size() ) break;

            for( int64_t i = 0; i < len; ++i ){
                uint64_t isamp0 = ::scheme::kinematics::bigindex_nest_part(samples[iresl][i].index);
                if( samples[iresl][i].score >= d.opt.global_score_cut ) continue;
                if( iresl == 0 ) ++d.non0_space_size;
                for( uint64_t j = 0; j < d.opt.DIMPOW2; ++j ){
                    uint64_t isamp = isamp0 * d.opt.DIMPOW2 + j;
                    samples[iresl+1].push_back( SearchPoint(DirectorIndex(isamp, ::scheme::kinematics::bigindex_scaffold_index(samples[iresl][i].index))) );
                }
            }
            if( 0 == samples[iresl+1].size() ){
                search_failed = true;
                std::cout << "search fail, no valid samples!" << std::endl;
                break;
            }
            samples[iresl].clear();

        }
        if( search_failed ) return false;
        std::cout << "full sort of final samples" << std::endl;
        __gnu_parallel::sort( samples.back().begin(), samples.back().end() );
    }
    if( search_failed ) return false;

    std::cout << "total non-0 space size was approx " << float(d.non0_space_size)*1024.0*1024.0*1024.0 << " grid points" << std::endl;
    std::cout << "total search effort " << KMGT(d.total_search_effort) << std::endl;


    hsearch_results_p = make_shared<std::vector< SearchPointWithRots >>();
    std::vector< SearchPointWithRots > & hsearch_results = *hsearch_results_p;


    hsearch_results.resize( samples.back().size() );
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1024)
    #endif
    for( int ipack = 0; ipack < hsearch_results.size(); ++ipack ){
        hsearch_results[ipack].score = samples.back()[ipack].score;
        hsearch_results[ipack].index = samples.back()[ipack].index;
    }



    return true;



}




#endif