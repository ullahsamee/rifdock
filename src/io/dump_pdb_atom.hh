#ifndef INCLUDED_io_dump_pdb_atom_HH
#define INCLUDED_io_dump_pdb_atom_HH

#include <iostream>
#include <cassert>

namespace scheme { namespace io {

inline void dump_pdb_atom(
	std::ostream & out,
	bool nothet,
	size_t iatom,
	std::string atom,
	std::string res,
	char chain,
	size_t ires,
	double x, double y, double z,
	double o, double b,
	std::string elem
){
	assert( atom.size()<5);
	assert( res.size()<4);
	assert( x<10000 && x > -1000 );
	assert( y<10000 && y > -1000 );
	assert( z<10000 && z > -1000 );
	// cout << "ATOM   1604  C   GLU A 220       5.010  12.933   1.553  1.00 41.10           C" << endl;
	char buf[128];
	snprintf(buf,128,"%s%5lu %4s %3s %c%4lu    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
		nothet?"ATOM  ":"HETATM",
		iatom,
		atom.c_str(),
		res.c_str(),
		chain,
		ires,
		x,y,z
		,o,b,elem.c_str()
	);
	out << buf;
}

template<class XYZ> void dump_pdb_atom(std::ostream & out,XYZ const & xyz){
	dump_pdb_atom(out,true,0,"ATOM","RES",' ',0,xyz[0],xyz[1],xyz[2],1,0,"");
}
template<class XYZ> void dump_pdb_atom(std::ostream & out,std::string aname,XYZ const & xyz){
	dump_pdb_atom(out,true,0,aname,"RES",' ',0,xyz[0],xyz[1],xyz[2],1,0,"");
}

}
}

#endif
