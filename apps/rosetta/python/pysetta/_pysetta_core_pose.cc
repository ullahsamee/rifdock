#include <pybind11/pybind11.h>

#include <core/pose/Pose.hh>

using namespace ::core::pose;

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_PLUGIN(_pysetta_core_pose) {
    py::module m("_pysetta_core_pose", "rosetta core::pose");

	py::class_< Pose, std::shared_ptr<Pose> >(m, "Pose")
		.def( py::init<>() )
		.def( "dump_pdb", (bool (Pose::*)( std::string const & file_name, std::string const & tag) const)
			 &Pose::dump_pdb, "dump Pose to pdb file", py::arg("file_name"), py::arg("tag") = "1" )
	;



    return m.ptr();
}

