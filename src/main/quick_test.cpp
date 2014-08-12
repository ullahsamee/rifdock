#include <main/gtest_util.hh>

// include desired tests here:

// #include <util/StoragePolicy_test.cc>
// #include <util/template_loop_test.cc> // 3.5 vs 2.2 wo Eigen
// #include <util/SimpleArray_test.cc>

// #include <nest/NEST_test.cc> // 4.479 vs. 3.0 wo Eigen
// #include <nest/NEST_neighbor_test.cc> // 4.567 // 3.26 wo Eigen
// #include <nest/MultiNest_test.cc>

// #include <nest/maps/parameter_maps_test.cc> // 7.546 // 4.0 wo Eigen
// #include <nest/maps/parameter_maps_test_nbrcell.cc> // 4.472 // 3.1 wo eigen
// #include <nest/maps/SphereDodec_test.cc>
// #include <nest/maps/SphereQuad_test.cc>

// #include <util/meta/util_test.cc>
// #include <util/container/ContainerInteractions_test.cc>
// #include <util/meta/InstanceMap_test.cc>
// #include <util/meta/InstanceMap_container_test.cc>
// #include <util/meta/InstanceMap_numeric_test.cc>

// #include <objective/ObjectiveFunction_test.cc>
// #include "objective/voxel/VoxelArray_test.cc"
#include "objective/voxel/FieldCache_test.cc"
// #include "objective/rosetta/AnalyticEvaluation_test.cc"
// #include "objective/rosetta/RosettaField_test.cc"

// #include <kinematics/Scene_test.cc>
// #include <kinematics/Scene_test_eigen.cc>
// #include <kinematics/SceneIterator_test.cc>

// #include <kinematics/Scene_test_objective.cc>

// #include <actor/ActorBBStub_test.cc>

int main(int argc, char **argv)
{
	std::vector<std::string> args;
	for(int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
	std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
	init_gtest_tests(args);
	return run_gtest_tests();
}
