// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


// This file should only be included in .cc files!!!!!

#ifndef INCLUDED_riflib_rifdock_subroutines_sketchy_include_hh
#define INCLUDED_riflib_rifdock_subroutines_sketchy_include_hh



#include <rif_dock_test.hh>
    #include <numeric/random/random.hh>

    #include <ObjexxFCL/format.hh>

    #include <boost/foreach.hpp>
    #include <boost/lexical_cast.hpp>
    // #include <boost/random/mersenne_twister.hpp>

    #include <core/id/AtomID.hh>
    #include <core/import_pose/import_pose.hh>
    #include <core/pose/Pose.hh>
    #include <core/pose/PDBInfo.hh>
    #include <core/pose/util.hh>
    #include <core/scoring/EnergyGraph.hh>
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
    #include <scheme/numeric/rand_xform.hh>
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
    #include <random>


    #include <scheme/objective/hash/XformHash.hh>


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




typedef int32_t intRot;





using ::scheme::make_shared;
using ::scheme::shared_ptr;













#endif