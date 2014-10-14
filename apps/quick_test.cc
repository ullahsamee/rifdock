#include "gtest_util.hh"

// include desired tests here:

// #include "util/hash.gtest.cc"
// #include <util/StoragePolicy.gtest.cc>
// #include <util/template_loop.gtest.cc> // 3.5 vs 2.2 wo Eigen
// #include <util/SimpleArray.gtest.cc>

// #include "numeric/euler_angles.gtest.cc"

// #include <nest/NEST.gtest.cc> // 4.479 vs. 3.0 wo Eigen
// #include <nest/NEST_neighbor.gtest.cc> // 4.567 // 3.26 wo Eigen
// #include <nest/MultiNest.gtest.cc>

// #include <nest/maps/parameter_maps.gtest.cc> // 7.546 // 4.0 wo Eigen
// #include <nest/maps/parameter_maps_test_nbrcell.gtest.cc> // 4.472 // 3.1 wo eigen
// #include <nest/maps/SphereDodec.gtest.cc>
// #include <nest/maps/SphereQuad.gtest.cc>
// #include <nest/maps/HecatonicosachoronMap.gtest.cc>
// #include <nest/maps/QuaternionMap.gtest.cc>
// #include <nest/maps/EulerAnglesMap.gtest.cc>
// #include <nest/maps/TetracontoctachoronMap.gtest.cc>

// #include "numeric/geom_4d.gtest.cc"
// #include "numeric/bcc_lattice.gtest.cc"
// #include "numeric/bcc_lattice_orientation.gtest.cc"

// #include <util/meta/util.gtest.cc>
// #include <util/container/ContainerInteractions.gtest.cc>
// #include <util/meta/InstanceMap.gtest.cc>
// #include <util/meta/InstanceMap_container.gtest.cc>
// #include <util/meta/InstanceMap_numeric.gtest.cc>

// #include <objective/ObjectiveFunction.gtest.cc>
// #include "objective/voxel/VoxelArray.gtest.cc"
// #include "objective/voxel/FieldCache.gtest.cc"

#include <kinematics/Scene.gtest.cc>
// #include <kinematics/Scene_test_eigen.gtest.cc>
// #include <kinematics/SceneIterator.gtest.cc>

// #include <io/dump_pdb_atom.gtest.cc>

// #include <kinematics/Scene_test_objective.gtest.cc>

// #include <actor/BBStub.gtest.cc>

// #include "rosetta/score/AnalyticEvaluation.gtest.cc"
// #include "rosetta/score/RosettaField.gtest.cc"

// #include "chemical/ligand_factory.gtest.cc"

// #include "objective/methods/hbond_5dof.gtest.cc"

int main(int argc, char **argv)
{
	std::vector<std::string> args;
	for(int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
	std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
	init_gtest_tests(args);
	return run_gtest_tests();
}
