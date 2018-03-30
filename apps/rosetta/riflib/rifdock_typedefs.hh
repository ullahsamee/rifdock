// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rifdock_typedefs_hh
#define INCLUDED_riflib_rifdock_typedefs_hh

#include <scheme/types.hh>

#include <scheme/actor/BackboneActor.hh>
#include <scheme/actor/VoxelActor.hh>
#include <scheme/actor/Atom.hh>
#include <scheme/kinematics/Scene.hh>
#include <scheme/scaffold/ScaffoldProviderBase.hh>
#include <scheme/nest/NEST.hh>
#include <scheme/nest/pmap/OriTransMap.hh>
#include <scheme/kinematics/Director.hh>
#include <scheme/objective/integration/SceneObjective.hh>
#include <scheme/chemical/RotamerIndex.hh>


#include <boost/mpl/vector.hpp>

#include <unordered_map>

namespace devel {
namespace scheme {


struct RIFAnchor {
    RIFAnchor() {}
};

struct ScaffoldDataCache;


// Typedefs that should never change

typedef ::scheme::actor::BackboneActor<EigenXform> BBActor;

typedef ::scheme::actor::VoxelActor<EigenXform,float> VoxelActor;

typedef ::scheme::actor::SimpleAtom< Eigen::Vector3f > SimpleAtom;


// Typedefs related to the Hierarchical Search Scene

typedef ::scheme::actor::Score_Voxel_vs_Atom<
        VoxelActor,
        SimpleAtom,
        false
    > MyClashScore;

typedef boost::mpl::vector<
            BBActor,
            SimpleAtom,
            VoxelActor,
            RIFAnchor
        > ParametricSceneContainers;

typedef ::scheme::kinematics::impl::Conformation<
    ParametricSceneContainers,
    ScaffoldDataCache > ParametricSceneConformation; 

typedef shared_ptr<ParametricSceneConformation> ParametricSceneConformationOP;
typedef shared_ptr<ParametricSceneConformation const > ParametricSceneConformationCOP;


typedef ::scheme::kinematics::Scene<
        ParametricSceneConformation,
        EigenXform
    > ParametricScene;
    

// Typedefs related to the Hierarchical Search ScaffoldProvider

typedef ::scheme::scaffold::TreeScaffoldProvider<ParametricSceneConformation> ScaffoldProvider;
typedef shared_ptr<ScaffoldProvider> ScaffoldProviderOP;

typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;




// If you add something to the index, you must follow these rules
// 1. Keep your item lightweight
// 2. Add your item to the hash functions
// 3. RifDockIndex() must provide a valid index

struct RifDockIndex {
    uint64_t nest_index;
    uint16_t seeding_index;
    ScaffoldIndex scaffold_index;

    RifDockIndex() :
        nest_index(0),
        seeding_index(0),
        scaffold_index(ScaffoldIndex()) {}
        

    RifDockIndex( 
        uint64_t nest_index_in,
        uint64_t seeding_index_in,
        ScaffoldIndex scaffold_index_in
        ) : 
        nest_index( nest_index_in ),
        seeding_index( seeding_index_in ),
        scaffold_index( scaffold_index_in ) {}

    bool operator==(RifDockIndex const & o) {
      return (
        nest_index == o.nest_index &&
        seeding_index == o.seeding_index &&
        scaffold_index == o.scaffold_index
        );
    }

};

inline
std::ostream & operator << ( std::ostream & out, RifDockIndex const & rdi ){
    out << "RifDockIndex: " << rdi.nest_index << " " << rdi.seeding_index << " " << rdi.scaffold_index;
    return out;
}


// Typedefs related to the Hierarchical Search Director

typedef ::scheme::nest::NEST< 6,
            devel::scheme::EigenXform,
            ::scheme::nest::pmap::OriTransMap,
            ::scheme::util::StoreNothing, // do not store a transform in the Nest
            uint64_t,
            float,
            false // do not inherit from NestBase
            > NestOriTrans6D;


// typedef ::scheme::kinematics::NestDirector< NestOriTrans6D > DirectorOriTrans6D;
typedef ::scheme::kinematics::NestDirector< NestOriTrans6D, RifDockIndex> RifDockNestDirector;

typedef ::scheme::kinematics::ScaffoldDirector< EigenXform, ScaffoldProvider, RifDockIndex > RifDockScaffoldDirector;

typedef ::scheme::kinematics::SeedingDirector< EigenXform, std::vector<EigenXform>, RifDockIndex > RifDockSeedingDirector;

typedef ::scheme::kinematics::CompositeDirector< EigenXform, RifDockIndex > RifDockDirector;

typedef shared_ptr<::scheme::kinematics::Director<EigenXform, RifDockIndex>> DirectorBase;


typedef shared_ptr< ::scheme::kinematics::SceneBase<EigenXform,uint64_t> > ScenePtr;
typedef shared_ptr< ::scheme::objective::integration::SceneOjbective<EigenXform,uint64_t> > ObjectivePtr;


struct SelectiveRifDockIndexHasher {
    SelectiveRifDockIndexHasher( 
        bool treat_nests_differently,
        bool treat_seeds_differently,
        bool treat_scaffolds_differently
        ) :
        treat_nests_differently_(treat_nests_differently),
        treat_seeds_differently_(treat_seeds_differently),
        treat_scaffolds_differently_(treat_scaffolds_differently)
        {}

    size_t operator() (devel::scheme::RifDockIndex const & rdi) const {

        using boost::hash;
        using boost::hash_combine;

        std::size_t the_hash = 0;

        if ( treat_nests_differently_ ) {
            boost::hash<int> hasher;
            hash_combine(the_hash, hasher(rdi.nest_index));
        }
        if ( treat_seeds_differently_ ) {
            boost::hash<int> seeding_index_hasher;
            hash_combine(the_hash, seeding_index_hasher(rdi.seeding_index));
        }
        if ( treat_scaffolds_differently_ ) {      
            std::hash<devel::scheme::ScaffoldIndex> scaffold_index_hasher;
            hash_combine(the_hash, scaffold_index_hasher(rdi.scaffold_index));
        }

        return the_hash;
    }
private:
    bool treat_nests_differently_;
    bool treat_seeds_differently_;
    bool treat_scaffolds_differently_;
};

struct SelectiveRifDockIndexEquater {
    SelectiveRifDockIndexEquater( 
        bool treat_nests_differently,
        bool treat_seeds_differently,
        bool treat_scaffolds_differently
        ) :
        treat_nests_differently_(treat_nests_differently),
        treat_seeds_differently_(treat_seeds_differently),
        treat_scaffolds_differently_(treat_scaffolds_differently)
        {}

    bool operator() (devel::scheme::RifDockIndex const & rdi1, devel::scheme::RifDockIndex const & rdi2) const {

        if ( treat_nests_differently_ ) {
            if ( rdi1.nest_index != rdi2.nest_index ) return false;
        }
        if ( treat_seeds_differently_ ) {
            if ( rdi1.seeding_index != rdi2.seeding_index ) return false;
        }
        if ( treat_scaffolds_differently_ ) {
            if ( ! ( rdi1.scaffold_index == rdi2.scaffold_index ) ) return false;
        }

        return true;
    }
private:
    bool treat_nests_differently_;
    bool treat_seeds_differently_;
    bool treat_scaffolds_differently_;
};


template< class __AnyPoint >
using _AnyPointVectorsMap = std::unordered_map< RifDockIndex, std::vector<__AnyPoint>, SelectiveRifDockIndexHasher, SelectiveRifDockIndexEquater >;



typedef ::scheme::actor::Atom<
    ::Eigen::Vector3f
> SchemeAtom;

typedef ::scheme::chemical::ChemicalIndex<
    ::scheme::chemical::AtomData
> ChemicalIndex;

struct RosettaRotamerGenerator;

typedef ::scheme::chemical::RotamerIndex<
    SchemeAtom,
    RosettaRotamerGenerator,
    EigenXform
> RotamerIndex;



}
}

// namespace std {

//     template <>
//     struct hash<devel::scheme::RifDockIndex>
//     {
//         std::size_t operator()(const devel::scheme::RifDockIndex& rdi) const {
//             using std::size_t;
//             using boost::hash;
//             using boost::hash_combine;

//             std::size_t the_hash = 0;

//             boost::hash<int> hasher;
//             hash_combine(the_hash, hasher(rdi.nest_index));
//             boost::hash<int> seeding_index_hasher;
//             hash_combine(the_hash, seeding_index_hasher(rdi.seeding_index));
//             std::hash<devel::scheme::ScaffoldIndex> scaffold_index_hasher;
//             hash_combine(the_hash, scaffold_index_hasher(rdi.scaffold_index));

//             return the_hash;
//         }
//     };

// }


#endif
