#include <main/gtest_util.hh>

// must include all tests here:
#include <nest/NEST_test.hh>
#include <util/dilated_int_test.hh>

int main(int argc, char **argv)
{
	std::vector<std::string> args;
	for(int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
	std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
	init_gtest_tests(args);
	return run_gtest_tests();
}
