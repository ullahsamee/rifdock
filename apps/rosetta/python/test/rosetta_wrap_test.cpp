#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/annotated_sequence.hh>

#include <iostream>
#include <vector>

void make_a_pdb(std::string s){
	std::cout << "make pose" << std::endl;
	core::pose::Pose pose;
	std::cout << "make_pose_from_sequence" << std::endl;	
	core::pose::make_pose_from_sequence(pose,s,core::chemical::FA_STANDARD);
	std::cout << "dumping pdb" << std::endl;	
	pose.dump_pdb("test.pdb");
	std::cout << "done" << std::endl;	
}

namespace py = pybind11;

PYBIND11_PLUGIN(example) {
    py::module m("example", "pybind11 example plugin");

    m.def("make_a_pdb", &make_a_pdb, "test make a pdb from sequence");

    return m.ptr();
}
