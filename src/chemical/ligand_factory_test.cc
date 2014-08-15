#include <gtest/gtest.h>

#include "actor/Atom.hh"
#include "ligand_factory.hh"
#include "util/SimpleArray.hh"
#include <iterator>     // std::back_inserter
#include <fstream>

#include <boost/foreach.hpp>

namespace scheme { namespace chemical { namespace test {

using std::cout;
using std::endl;

TEST(ligand_factory,make_atom_pdbline){
	typedef util::SimpleArray<3,double> Position;
	typedef actor::Atom<Position> Atom;
	LigandFactory<Atom> f;

	std::string l;
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( false,  f.make_atom_pdbline(l).data().ishet );
	l="HETATM    7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( true ,  f.make_atom_pdbline(l).data().ishet );

	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( 7,  f.make_atom_pdbline(l).data().atomnum );
	l="ATOM  999999 C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( 999999,  f.make_atom_pdbline(l).data().atomnum );
	l="ATOM  9999999C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( 999999,  f.make_atom_pdbline(l).data().atomnum );

	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( "C3", f.make_atom_pdbline(l).data().atomname );
	l="ATOM      7 ATOMABTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( "ATOM",  f.make_atom_pdbline(l).data().atomname );

	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( "BTN",  f.make_atom_pdbline(l).data().resname );
	l="ATOM      7 ATOMABTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( "ABTN",  f.make_atom_pdbline(l).data().resname );

	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( 'X',  f.make_atom_pdbline(l).data().chain );
	l="ATOM      7  C3  BTN     1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( ' ',  f.make_atom_pdbline(l).data().chain );

	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( 1,  f.make_atom_pdbline(l).data().resnum );
	l="ATOM      7  C3  BTN X9999      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( 9999,  f.make_atom_pdbline(l).data().resnum );
	l="ATOM      7  C3  BTN X9999      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( 9999,  f.make_atom_pdbline(l).data().resnum );
	l="ATOM      7  C3  BTN X999999    -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( 999999,  f.make_atom_pdbline(l).data().resnum );

	l="ATOM      7  C3  BTN X   1  -99999.999 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( -99999.999,  f.make_atom_pdbline(l).position()[0] );
	l="ATOM      7  C3  BTN X   1      -0.470-999.999   5.377  1.00 20.00           C"; ASSERT_EQ(   -999.999,  f.make_atom_pdbline(l).position()[1] );
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087-999.999  1.00 20.00           C"; ASSERT_EQ(   -999.999,  f.make_atom_pdbline(l).position()[2] );

	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377111.00 20.00           C"; ASSERT_EQ( 111,  f.make_atom_pdbline(l).data().occ );
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00120.00           C"; ASSERT_EQ( 120,  f.make_atom_pdbline(l).data().bfac );

	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( "C",  f.make_atom_pdbline(l).data().elem );
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00  AAAAAAAAAA"; ASSERT_EQ( "AAAAAAAAAA",  f.make_atom_pdbline(l).data().elem );


	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( l, io::dump_pdb_atom(f.make_atom_pdbline(l)) );
	l="ATOM      7  C   BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( l, io::dump_pdb_atom(f.make_atom_pdbline(l)) );
	l="ATOM      7  C3A BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( l, io::dump_pdb_atom(f.make_atom_pdbline(l)) );
	l="ATOM      7 AC3A BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C"; ASSERT_EQ( l, io::dump_pdb_atom(f.make_atom_pdbline(l)) );


	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00"; ASSERT_EQ( false,  f.make_atom_pdbline(l).data().ishet );
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00"; ASSERT_EQ( 7,  f.make_atom_pdbline(l).data().atomnum );
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00"; ASSERT_EQ( "C3", f.make_atom_pdbline(l).data().atomname );
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00"; ASSERT_EQ( "BTN",  f.make_atom_pdbline(l).data().resname );
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00"; ASSERT_EQ( 'X',  f.make_atom_pdbline(l).data().chain );
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00"; ASSERT_EQ( 1,  f.make_atom_pdbline(l).data().resnum );
	l="ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00"; ASSERT_EQ( "",  f.make_atom_pdbline(l).data().elem );

}

TEST(ligand_factory,make_biotin){
	typedef util::SimpleArray<3,double> Position;
	typedef actor::Atom<Position> Atom;
	LigandFactory<Atom> f;

	std::vector<Atom> btn;
	f.make_biotin_minimal(std::back_inserter(btn));
	// BOOST_FOREACH(Atom a,btn) cout << io::dump_pdb_atom(a) << endl;
	ASSERT_EQ( io::dump_pdb_atom(btn[0]), "ATOM      1  N1  BTN X   1       0.696 -12.422   3.375  1.00 20.00           N" );
	ASSERT_EQ( io::dump_pdb_atom(btn[1]), "ATOM      2  S1  BTN X   1       0.576  -9.666   5.336  1.00 20.00           S" );
	ASSERT_EQ( io::dump_pdb_atom(btn[2]), "ATOM      3  C1  BTN X   1      -0.523 -10.824   6.189  1.00 20.00           C" );
	ASSERT_EQ( io::dump_pdb_atom(btn[3]), "ATOM      4  N2  BTN X   1      -1.324 -12.123   4.201  1.00 20.00           N" );
	ASSERT_EQ( io::dump_pdb_atom(btn[4]), "ATOM      5  C2  BTN X   1      -0.608 -12.327   3.072  1.00 20.00           C" );
	ASSERT_EQ( io::dump_pdb_atom(btn[5]), "ATOM      6  O1  BTN X   1      -1.125 -12.422   1.933  1.00 20.00           O" );
	ASSERT_EQ( io::dump_pdb_atom(btn[6]), "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C" );
	ASSERT_EQ( io::dump_pdb_atom(btn[7]), "ATOM      8  C4  BTN X   1       0.953 -12.267   4.780  1.00 20.00           C" );
	ASSERT_EQ( io::dump_pdb_atom(btn[8]), "ATOM      9  C5  BTN X   1       1.765 -11.040   5.134  1.00 20.00           C" );
	ASSERT_EQ( io::dump_pdb_atom(btn[9]), "ATOM     10  C6  BTN X   1      -1.836 -10.395   6.850  1.00 20.00           C" );

}

TEST(ligand_factory,make_aas){
	typedef util::SimpleArray<3,double> Position;
	typedef actor::Atom<Position> Atom;
	LigandFactory<Atom> f;

	std::vector<Atom> PHE;
	f.make_atoms( std::back_inserter(PHE), "PHE", false );
	ASSERT_EQ( PHE.size(), 7 );
	BOOST_FOREACH(Atom a,PHE) ASSERT_GT( a.type(), 0 );
	// std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,PHE) io::dump_pdb_atom(out,a); out.close();

	std::vector<Atom> PHE_H;
	f.make_atoms( std::back_inserter(PHE_H), "PHE", true );
	ASSERT_EQ( PHE_H.size(), 12 );
	BOOST_FOREACH(Atom a,PHE_H) ASSERT_GT( a.type(), 0 );
	// std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,PHE_H) io::dump_pdb_atom(out,a); out.close();

	std::vector<Atom> TYR;
	f.make_atoms( std::back_inserter(TYR), "TYR", false );
	ASSERT_EQ( TYR.size(), 8 );
	BOOST_FOREACH(Atom a,TYR) ASSERT_GT( a.type(), 0 );
	// std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,TYR) io::dump_pdb_atom(out,a); out.close();

	std::vector<Atom> TYR_H;
	f.make_atoms( std::back_inserter(TYR_H), "TYR", true );
	// std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,TYR_H) io::dump_pdb_atom(out,a); out.close();
	// BOOST_FOREACH(Atom a,TYR_H) cout << a << endl;
	ASSERT_EQ( TYR_H.size(), 12 );
	BOOST_FOREACH(Atom a,TYR_H) ASSERT_GT( a.type(), 0 );

	std::vector<Atom> ASP_H;
	f.make_atoms( std::back_inserter(ASP_H), "ASP", true );
	// std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,ASP_H) io::dump_pdb_atom(out,a); out.close();
	// BOOST_FOREACH(Atom a,ASP_H) cout << a << endl;
	ASSERT_EQ( ASP_H.size(), 4 );
	BOOST_FOREACH(Atom a,ASP_H) ASSERT_GT( a.type(), 0 );

	std::vector<Atom> GLU_H;
	f.make_atoms( std::back_inserter(GLU_H), "GLU", true );
	// std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,GLU_H) io::dump_pdb_atom(out,a); out.close();
	// BOOST_FOREACH(Atom a,GLU_H) cout << a << endl;
	ASSERT_EQ( GLU_H.size(), 4 );
	BOOST_FOREACH(Atom a,GLU_H) ASSERT_GT( a.type(), 0 );

	std::vector<Atom> ASN_H;
	f.make_atoms( std::back_inserter(ASN_H), "ASN", true );
	// std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,ASN_H) io::dump_pdb_atom(out,a); out.close();
	// BOOST_FOREACH(Atom a,ASN_H) cout << a << endl;
	ASSERT_EQ( ASN_H.size(), 6 );
	BOOST_FOREACH(Atom a,ASN_H) ASSERT_GT( a.type(), 0 );

	std::vector<Atom> TRP_H;
	f.make_atoms( std::back_inserter(TRP_H), "TRP", true );
	// std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,TRP_H) io::dump_pdb_atom(out,a); out.close();
	// BOOST_FOREACH(Atom a,TRP_H) cout << a << endl;
	ASSERT_EQ( TRP_H.size(), 16 );
	BOOST_FOREACH(Atom a,TRP_H) ASSERT_GT( a.type(), 0 );

	std::vector<Atom> TRP;
	f.make_atoms( std::back_inserter(TRP), "TRP", false );
	// std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,TRP) io::dump_pdb_atom(out,a); out.close();
	// BOOST_FOREACH(Atom a,TRP) cout << a << endl;
	ASSERT_EQ( TRP.size(), 10 );
	BOOST_FOREACH(Atom a,TRP) ASSERT_GT( a.type(), 0 );

}

}
}
}
