// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_util_hh
#define INCLUDED_riflib_rifdock_subroutines_util_hh


#include <riflib/types.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

#include <numeric/xyzTransform.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <scheme/kinematics/Director.hh>

#include <rif_dock_test.hh>
#include <riflib/rotamer_energy_tables.hh>
#include <riflib/RifFactory.hh>



using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


inline
Eigen::Vector3f
pose_center(
    core::pose::Pose const & pose,
    utility::vector1<core::Size> const & useres = utility::vector1<core::Size>()
){
    typedef numeric::xyzVector<core::Real> Vec;
    Vec cen(0,0,0);
    int count = 0;
    for( int ir = 1; ir <= pose.size(); ++ir ) {
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

inline
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


inline
void xform_pose( core::pose::Pose & pose, numeric::xyzTransform<float> s, core::Size sres=1, core::Size eres=0 ) {
  if(eres==0) eres = pose.size();
  for(core::Size ir = sres; ir <= eres; ++ir) {
    for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s*pose.xyz(aid) );
    }
  }
}




template<class _DirectorIndex>
struct tmplRifDockResult {
    typedef _DirectorIndex DirectorIndex;
    typedef tmplRifDockResult<DirectorIndex> This;
    float dist0, packscore, nopackscore, rifscore, stericscore;
    uint64_t isamp;
    DirectorIndex scene_index;
    uint32_t prepack_rank;
    float cluster_score;
    bool operator< ( This const & o ) const { return packscore < o.packscore; }
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
    core::pose::PoseOP pose_ = nullptr;
    size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
    std::vector< std::pair<intRot,intRot> > const & rotamers() const { assert(rotamers_!=nullptr); return *rotamers_; }
};




#pragma pack (push, 4) // allows size to be 12 rather than 16
template<class _DirectorIndex>
struct tmplSearchPoint {
    typedef _DirectorIndex DirectorIndex;
    typedef tmplSearchPoint<DirectorIndex> This;
    float score;
    DirectorIndex index;
    tmplSearchPoint() : score(9e9) {
        index = ::scheme::kinematics::director_index_default_value(index);
    }
    tmplSearchPoint(DirectorIndex i) : score(9e9), index(i) {}
    bool operator < (This const & o) const {
        return score < o.score;
    }
};
#pragma pack (pop)



template<class _DirectorIndex>
struct tmplSearchPointWithRots {
    typedef _DirectorIndex DirectorIndex;
    typedef tmplSearchPointWithRots<DirectorIndex> This;
    float score;
    uint32_t prepack_rank;
    DirectorIndex index;
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
    core::pose::PoseOP pose_ = nullptr;
    tmplSearchPointWithRots() : score(9e9), prepack_rank(0), rotamers_(nullptr) {
        index = ::scheme::kinematics::director_index_default_value(index);
    }
    tmplSearchPointWithRots(DirectorIndex i, uint32_t orank) : score(9e9), prepack_rank(orank), index(i), rotamers_(nullptr) {}
    // ~SearchPointWithRots() { delete rotamers_; }
    void checkinit() { if( rotamers_==nullptr ) rotamers_ = make_shared< std::vector< std::pair<intRot,intRot> > > ();  }
    std::vector< std::pair<intRot,intRot> > & rotamers() { checkinit(); return *rotamers_; }
    std::vector< std::pair<intRot,intRot> > const & rotamers() const { runtime_assert(rotamers_!=nullptr); return *rotamers_; }
    size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
    bool operator < (This const & o) const {
        return score < o.score;
    }
    friend void swap(This & a, This & b){
        std::swap( a.score, b.score );
        std::swap( a.prepack_rank, b.prepack_rank );
        std::swap( a.index, b.index );
        std::swap( a.rotamers_, b.rotamers_ );
        std::swap( a.pose_, b.pose_ );
    }
};


// Convenience templates for the above templated containers

template <class __Director>
using _SearchPointWithRots = tmplSearchPointWithRots<_DirectorBigIndex<__Director>>;

template <class __Director>
using _RifDockResult = tmplRifDockResult<_DirectorBigIndex<__Director>>;

template <class __Director>
using _SearchPoint = tmplSearchPoint<_DirectorBigIndex<__Director>>;




template<class DirectorBase, class ScaffoldProvider >
struct RifDockData {
    RifDockOpt & opt;
    std::vector<float> & RESLS;
    DirectorBase & director;
    std::vector< devel::scheme::ScenePtr > & scene_pt;
    devel::scheme::ScenePtr & scene_minimal;
    std::vector<devel::scheme::SimpleAtom> & target_simple_atoms;
    std::vector< devel::scheme::VoxelArrayPtr > & target_field_by_atype;
    std::vector< std::vector< devel::scheme::VoxelArrayPtr > > const * target_bounding_by_atype;
    std::vector< ::scheme::chemical::HBondRay > * target_donors;
    std::vector< ::scheme::chemical::HBondRay > * target_acceptors;
    float & target_redundancy_filter_rg;
    core::pose::Pose & target;
    shared_ptr< devel::scheme::RotamerIndex > & rot_index_p;
    ::devel::scheme::RotamerRFTablesManager & rotrf_table_manager;
    std::vector< devel::scheme::ObjectivePtr > & objectives;
    devel::scheme::ObjectivePtr & packing_objective;
    ::scheme::search::HackPackOpts & packopts;
    std::vector<shared_ptr<devel::scheme::RifBase> > & rif_ptrs;

    shared_ptr<ScaffoldProvider> scaffold_provider;
};











#endif