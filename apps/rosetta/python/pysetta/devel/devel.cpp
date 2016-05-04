#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <devel/init.hh>

namespace py = pybind11;

void rosetta_init( std::vector<std::string> args ) {
	std::cout << "initializing rosetta with args:" << std::endl;
	char** cargs = new char*[args.size()];
	int i = 0;
	for( auto const & arg : args ){
		std::cout << arg << std::endl;
		cargs[i] = new char[arg.size()+1];
		for( int j = 0; j < arg.size(); ++j ){
			cargs[i][j] = arg[j];
		}
		cargs[i][arg.size()] = 0;
		++i;
	}
	::devel::init( args.size(), cargs ); 
}

PYBIND11_PLUGIN(_devel) {
    py::module m("_devel", "rosetta devel namespace");
    m.def("init", &rosetta_init, "init rosetta from 'cli args'");
    return m.ptr();
}

