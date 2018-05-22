// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_task_types_hh
#define INCLUDED_riflib_task_types_hh


#include <riflib/types.hh>
#include <rif_dock_test.hh>


#include <core/pose/Pose.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rotamer_energy_tables.hh>
#include <scheme/search/HackPack.hh>
#include <riflib/RifBase.hh>
#include <riflib/RifFactory.hh>

#include <utility/io/ozstream.hh>

#ifdef USEGRIDSCORE
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#endif

#include <chrono>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;

namespace devel {
namespace scheme {




template<class _DirectorBigIndex>
struct tmplSearchPointWithRots;

template<class _DirectorBigIndex>
struct tmplRifDockResult;


#pragma pack (push, 4) // allows size to be 12 rather than 16
template<class _DirectorBigIndex>
struct tmplSearchPoint {
    typedef _DirectorBigIndex DirectorBigIndex;
    typedef tmplSearchPoint<DirectorBigIndex> This;
    float score;
    uint64_t sasa;
    DirectorBigIndex index;
    tmplSearchPoint() : score(9e9) {}
    tmplSearchPoint(DirectorBigIndex i) : score(9e9), index(i) {}
    bool operator < (This const & o) const {
        return score < o.score;
    }
    This& operator=( tmplSearchPointWithRots<DirectorBigIndex> const & ot ) {
        score = ot.score;
        index = ot.index;
        sasa = ot.sasa;
        return *this;
    }
    This& operator=( tmplRifDockResult<DirectorBigIndex> const & ot ) {
        score = ot.score;
        index = ot.index;
        sasa = ot.sasa;
        return *this;
    }
};
#pragma pack (pop)


template<class _DirectorBigIndex>
struct tmplSearchPointWithRots {
    typedef _DirectorBigIndex DirectorBigIndex;
    typedef tmplSearchPointWithRots<DirectorBigIndex> This;
    float score;
    uint64_t sasa;
    uint32_t prepack_rank;
    DirectorBigIndex index;
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
    core::pose::PoseOP pose_ = nullptr;
    tmplSearchPointWithRots() : score(9e9), prepack_rank(0), rotamers_(nullptr) {}
    tmplSearchPointWithRots(DirectorBigIndex i, uint32_t orank) : score(9e9), prepack_rank(orank), index(i), rotamers_(nullptr) {}
    // ~SearchPointWithRots() { delete rotamers_; }
    void checkinit() { if( rotamers_==nullptr ) rotamers_ = make_shared< std::vector< std::pair<intRot,intRot> > > ();  }
    std::vector< std::pair<intRot,intRot> > & rotamers() { checkinit(); return *rotamers_; }
    std::vector< std::pair<intRot,intRot> > const & rotamers() const { runtime_assert(rotamers_!=nullptr); return *rotamers_; }
    size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
    bool operator < (This const & o) const {
        return score < o.score;
    }
    This& operator=( tmplSearchPoint<DirectorBigIndex> const & ot ) {
        score = ot.score;
        index = ot.index;
        sasa = ot.sasa;
        return *this;
    }
    This& operator=( tmplRifDockResult<DirectorBigIndex> const & ot ) {
        score = ot.score;
        prepack_rank = ot.prepack_rank;
        index = ot.index;
        rotamers_ = ot.rotamers_;
        pose_ = ot.pose_;
        sasa = ot.sasa;
        return *this;
    }
};

template<class _DirectorBigIndex>
struct tmplRifDockResult {
    typedef _DirectorBigIndex DirectorBigIndex;
    typedef tmplRifDockResult<DirectorBigIndex> This;
    float dist0, nopackscore, rifscore, stericscore;
    float score; // formerly packscore
    uint64_t sasa;
    uint64_t isamp;
    DirectorBigIndex index; // formerly index
    uint32_t prepack_rank;
    float cluster_score;
    bool operator< ( This const & o ) const { return score < o.score; }
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
    core::pose::PoseOP pose_ = nullptr;
    size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
    std::vector< std::pair<intRot,intRot> > const & rotamers() const { assert(rotamers_!=nullptr); return *rotamers_; }

    This& operator=( tmplSearchPoint<DirectorBigIndex> const & ot ) {
        score = ot.score;
        index = ot.index;
        sasa = ot.sasa;
        return *this;
    }
    This& operator=( tmplSearchPointWithRots<DirectorBigIndex> const & ot ) {
        score = ot.score;
        prepack_rank = ot.prepack_rank;
        index = ot.index;
        rotamers_ = ot.rotamers_;
        pose_ = ot.pose_;
        sasa = ot.sasa;
        return *this;
    }
};

struct SasaComparator
{

    template <typename AnyPoint>
    inline bool operator()(AnyPoint const & lhs, AnyPoint const & rhs)
    {
      return lhs.sasa > rhs.sasa;
    }
};


struct ScorePer1000SasaComparator
{

    template <typename AnyPoint>
    inline bool operator()(AnyPoint const & lhs, AnyPoint const & rhs)
    {
      return lhs.score / lhs.sasa < rhs.score / rhs.sasa;
    }
};



// Convenience templates for the above templated containers

template <class __Director>
using _SearchPointWithRots = tmplSearchPointWithRots<_DirectorBigIndex<__Director>>;
typedef _SearchPointWithRots<DirectorBase> SearchPointWithRots;

template <class __Director>
using _RifDockResult = tmplRifDockResult<_DirectorBigIndex<__Director>>;
typedef _RifDockResult<DirectorBase> RifDockResult;

template <class __Director>
using _SearchPoint = tmplSearchPoint<_DirectorBigIndex<__Director>>;
typedef _SearchPoint<DirectorBase> SearchPoint;






struct RifDockData {
    int iscaff;
    RifDockOpt & opt;
    std::vector<float> & RESLS;
    DirectorBase & director;
    std::vector< ScenePtr > & scene_pt;
    ScenePtr & scene_minimal;
    std::vector<SimpleAtom> & target_simple_atoms;
    std::vector< VoxelArrayPtr > & target_field_by_atype;
    std::vector< std::vector< VoxelArrayPtr > > const * target_bounding_by_atype;
    std::vector< ::scheme::chemical::HBondRay > * target_donors;
    std::vector< ::scheme::chemical::HBondRay > * target_acceptors;
    RifScoreRotamerVsTarget const & rot_tgt_scorer;
    float & target_redundancy_filter_rg;
    core::pose::Pose & target;
    shared_ptr< RotamerIndex > & rot_index_p;
    RotamerRFTablesManager & rotrf_table_manager;
    std::vector< ObjectivePtr > & objectives;
    std::vector< ObjectivePtr > & packing_objectives;
    ::scheme::search::HackPackOpts & packopts;
    std::vector<shared_ptr<RifBase> > & rif_ptrs;
    RifSceneObjectiveConfig & rso_config;
    shared_ptr<RifFactory> & rif_factory;
    NestOriTrans6D const & nest;
    #ifdef USE_OPENMP
        omp_lock_t & dump_lock;
    #endif

    utility::io::ozstream & dokout;
    ScaffoldProviderOP scaffold_provider;
    shared_ptr<BurialManager> burial_manager;
    shared_ptr<UnsatManager> unsat_manager;

#ifdef USEGRIDSCORE
    shared_ptr<protocols::ligand_docking::ga_ligand_dock::GridScorer> grid_scorer;
#endif
};

struct ProtocolData {

// Data related to the original RifDock as written by Will
    int64_t non0_space_size;
    int64_t total_search_effort;
    int64_t npack;

    double time_rif;
    double time_pck; 
    double time_ros;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_rif;

    double hsearch_rate; // instantaneous rate of hsearch

// for constraints
    std::vector< ScaffoldIndex > unique_scaffolds;

// for hsearch
    double beam_multiplier;

// for seeding positions
    std::vector<std::string> seeding_tags;



    ProtocolData() :
    non0_space_size(0),
    total_search_effort(0),
    npack(0),
    time_rif(0),
    time_pck(0),
    time_ros(0),
    hsearch_rate(0),
    beam_multiplier(1)



    {}



};





}}



#endif