#include "gtest_util.hh"

// include desired tests here:

// #include "scheme/util/hash.gtest.cc"
// #include "scheme/util/StoragePolicy.gtest.cc"
// #include "scheme/util/template_loop.gtest.cc" // 3.5 vs 2.2 wo Eigen
// #include "scheme/util/SimpleArray.gtest.cc"

// #include "scheme/numeric/euler_angles.gtest.cc"

// #include "scheme/dock/fftdock.gtest.cc"

// #include "scheme/nest/NEST.gtest.cc" // 4.479 vs. 3.0 wo Eigen
// #include "scheme/nest/NEST_neighbor.gtest.cc" // 4.567 // 3.26 wo Eigen
// #include "scheme/nest/MultiNest.gtest.cc"

// #include "scheme/nest/maps/parameter_maps.gtest.cc" // 7.546 // 4.0 wo Eigen
// #include "scheme/nest/maps/parameter_maps_test_nbrcell.gtest.cc" // 4.472 // 3.1 wo eigen
// #include "scheme/nest/maps/SphereDodec.gtest.cc"
// #include "scheme/nest/maps/SphereQuad.gtest.cc"
// #include "scheme/nest/maps/HecatonicosachoronMap.gtest.cc"
// #include "scheme/nest/maps/QuaternionMap.gtest.cc"
// #include "scheme/nest/maps/EulerAnglesMap.gtest.cc"
// #include "scheme/nest/maps/TetracontoctachoronMap.gtest.cc"
#include "scheme/nest/maps/Rotation1DMap.gtest.cc"

// #include "scheme/numeric/geom_4d.gtest.cc"
// #include "scheme/numeric/bcc_lattice.gtest.cc"
// #include "scheme/numeric/bcc_lattice_orientation.gtest.cc"

// #include "scheme/util/meta/util.gtest.cc"
// #include "scheme/util/container/ContainerInteractions.gtest.cc"
// #include "scheme/util/meta/InstanceMap.gtest.cc"
// #include "scheme/util/meta/InstanceMap_container.gtest.cc"
// #include "scheme/util/meta/InstanceMap_numeric.gtest.cc"

// #include "scheme/objective/ObjectiveFunction.gtest.cc"
// #include "scheme/objective/voxel/VoxelArray.gtest.cc"
// #include "scheme/objective/voxel/FieldCache.gtest.cc"

// #include "scheme/kinematics/Scene.gtest.cc"
// #include "scheme/kinematics/Scene_test_eigen.gtest.cc"
// #include "scheme/kinematics/SceneIterator.gtest.cc"

// #include "scheme/io/dump_pdb_atom.gtest.cc"

// #include "scheme/kinematics/Scene_test_objective.gtest.cc"

// #include "scheme/actor/BBStub.gtest.cc"

// #include "scheme/rosetta/score/AnalyticEvaluation.gtest.cc"
// #include "scheme/rosetta/score/RosettaField.gtest.cc"

// #include "scheme/chemical/ligand_factory.gtest.cc"

// #include "scheme/objective/methods/hbond_5dof.gtest.cc"

int main(int argc, char **argv)
{
	std::vector<std::string> args;
	for(int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
	std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
	init_gtest_tests(args);
	return run_gtest_tests();
}
