#include <gtest/gtest.h>

#include "actor/Atom.hh"
#include "chemical/ligand_factory.hh"
#include "numeric/euler_angles.hh"

#include <Eigen/Geometry>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>


#include <fstream>

namespace scheme {

using std::cout;
using std::endl;

template<class T>
std::string str(T const & t){ return boost::lexical_cast<std::string>(t); }

template<class M,class T>
std::vector<T> operator*(M const & m, std::vector<T> const & v){
	std::vector<T> r(v);
	BOOST_FOREACH(T & t,r) t = T(t,m);
	return r;
}

template<class T>
void dump_pdb(std::vector<T> const & v,std::string fname){
	std::ofstream out(fname.c_str());
	BOOST_FOREACH(T const & t,v) io::dump_pdb_atom(out,t);
	out.close();
}

TEST(hbond_5dof,test){
	using namespace Eigen;
	typedef Transform<double,3,AffineCompact> Xform;
	typedef Vector3d Position;
	typedef actor::Atom<Position> Atom;
	chemical::LigandFactory<Atom> f;

	std::vector<Atom> gly,trp,trp1;
	f.make_atoms( std::back_inserter(trp), "TRP", true );
	f.make_atoms( std::back_inserter(gly), "GLY", true );
	gly = AngleAxisd(0.5*M_PI,Vector3d::UnitY()) * gly;

	boost::random::mt19937 rng((unsigned int)time(0));
	boost::normal_distribution<> gauss;
	boost::uniform_real<> uniform;


	Eigen::Quaterniond q( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
	q.normalize();
	Xform frame1( q.matrix() ); frame1.translation() = Vector3d(7,2,3);
	// frame1.makeAffine();

	dump_pdb(frame1*gly,"test1.pdb");
	double e0 = uniform(rng)*2*M_PI;
	double e2 = uniform(rng)*M_PI;
	for(double e1 = 0; e1 < 2*M_PI; e1+=M_PI/8){

		Matrix3d r; numeric::from_euler_angles(Vector3d(e0,e1,e2),r);
		Xform x(r);
		x.translation() = Vector3d(1,2,3);
		x = x.inverse();
		x = frame1*x;

		trp1 = x*trp;
		dump_pdb(trp1,"test"+str(e1)+".pdb");

		Matrix3d m = ( x.inverse()*frame1 ).linear();
		Vector3d euler; numeric::euler_angles( m, euler );
		printf("%7.4f %7.4f %7.4f\n",euler[0],euler[1],euler[2]);

	}





}

}
