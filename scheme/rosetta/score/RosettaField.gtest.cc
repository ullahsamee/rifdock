#include <gtest/gtest.h>

#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "scheme/rosetta/score/RosettaField.hh"
#include "scheme/objective/voxel/FieldCache.hh"
#include "scheme/actor/Atom.hh"
// #include <Eigen/core>

namespace scheme { namespace rosetta { namespace score { namespace test {

using std::cout;
using std::endl;


// 1  CNH2 
// 2  COO  
// 3  CH1  
// 4  CH2  
// 5  CH3  
// 6  aroC 
// 7  Ntrp 
// 8  Nhis 
// 9  NH2O 
// 10 Nlys 
// 11 Narg 
// 12 Npro 
// 13 OH   
// 14 ONH2 
// 15 OOC  
// 16 Oaro 
// 17 S    
// 18 Nbb  
// 19 CAbb 
// 20 CObb 
// 21 OCbb 
// 22 Phos 
// 23 Pbb  
// 24 Hpol 
// 25 Hapo 
// 26 Haro 

// ATOM      1  N1  BTN X   1       0.696 -12.422   3.375  1.00 20.00           N  
// ATOM      2  S1  BTN X   1       0.576  -9.666   5.336  1.00 20.00           S  
// ATOM      3  C1  BTN X   1      -0.523 -10.824   6.189  1.00 20.00           C  
// ATOM      4  N2  BTN X   1      -1.324 -12.123   4.201  1.00 20.00           N  
// ATOM      5  C2  BTN X   1      -0.608 -12.327   3.072  1.00 20.00           C  
// ATOM      6  O1  BTN X   1      -1.125 -12.422   1.933  1.00 20.00           O  
// ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00           C  
// ATOM      8  C4  BTN X   1       0.953 -12.267   4.780  1.00 20.00           C  
// ATOM      9  C5  BTN X   1       1.765 -11.040   5.134  1.00 20.00           C  
// ATOM     10  C6  BTN X   1      -1.836 -10.395   6.850  1.00 20.00           C  

// ATOM  C6  CH2   X   -0.04
// ATOM  C1  CH1   X   0.03
// ATOM  S1  S     X   -0.15
// ATOM  C5  CH2   X   0.01
// ATOM  C2  CH1   X   0.05
// ATOM  N1  Ntrp  X   -0.29
// ATOM  C4  aroC  X   0.30
// ATOM  O1  OH    X   -0.25
// ATOM  N2  Ntrp  X   -0.29
// ATOM  C3  CH1   X   0.06


TEST(RosettaField,test_btn){
	typedef util::SimpleArray<3,float> F3;
	typedef actor::Atom<F3> Atom;
	std::vector<Atom> atoms;
	atoms.push_back( Atom( F3( 0.696,-12.422,3.375), 7  ));
	atoms.push_back( Atom( F3( 0.576, -9.666,5.336), 17 ));
	atoms.push_back( Atom( F3(-0.523,-10.824,6.189), 3  ));
	atoms.push_back( Atom( F3(-1.324,-12.123,4.201), 7  ));
	atoms.push_back( Atom( F3(-0.608,-12.327,3.072), 3  ));
	atoms.push_back( Atom( F3(-1.125,-12.422,1.933), 13 ));
	atoms.push_back( Atom( F3(-0.470,-12.087,5.377), 3  ));
	atoms.push_back( Atom( F3( 0.953,-12.267,4.780), 6  ));
	atoms.push_back( Atom( F3( 1.765,-11.040,5.134), 4  ));
	atoms.push_back( Atom( F3(-1.836,-10.395,6.850), 4  ));

	RosettaField<Atom> rf(atoms);
	F3 lb(9e9,9e9,9e9),ub(-9e9,-9e9,-9e9);
	BOOST_FOREACH(Atom const & a,rf.atoms_){ 
		lb = lb.min(a.position());
		ub = ub.max(a.position());
	}
	// cout << "LB " << lb << endl;
	// cout << "UB " << ub << endl;

	for(int atype = 1; atype <= 25; ++atype){
		RosettaFieldAtype<Atom> rfa(rf,1);
		objective::voxel::FieldCache3D<float> rc(rfa,lb-6.0f,ub+6.0f,1.0);
		objective::voxel::BoundingFieldCache3D<float> brc(rc,2.0,1.0);

		boost::random::mt19937 rng((unsigned int)time(0));
		boost::uniform_real<> uniform;
		for(int i = 0; i < 10000; ++i){
			F3 idx = F3( uniform(rng), uniform(rng), uniform(rng) ) * (rc.ub_-rc.lb_) + rc.lb_;
			ASSERT_LE( brc[idx], rc[idx] );
		}

	}

	int atype = 1;
	RosettaFieldAtype<Atom> rfa(rf,1);
	objective::voxel::FieldCache3D<float> rc(rfa,lb-6.0f,ub+6.0f,0.25);
	objective::voxel::BoundingFieldCache3D<float> brc(rc,1.0,0.25);

		// size_t nbz=0, naz=0;
		// for(size_t i = 0; i < rc.num_elements(); ++i){
		// 	if(rc.data()[i]>0) ++naz;
		// 	if(rc.data()[i]<0) ++nbz;
		// }
		// cout << rc.num_elements() << " " << (float)naz/rc.num_elements() << " " << (float)nbz/rc.num_elements() << endl;

	std::ofstream out("/tmp/btn1.dat");
	cout << rc.shape()[0] << " " << rc.shape()[1] << " " << rc.shape()[2] << endl;
	for(size_t i = 0; i < rc.num_elements(); ++i) out << rc.data()[i] << endl;
	out.close();

	std::ofstream out2("/tmp/btn2.dat");
	cout << brc.shape()[0] << " " << brc.shape()[1] << " " << brc.shape()[2] << endl;
	for(size_t i = 0; i < brc.num_elements(); ++i) out2 << brc.data()[i] << endl;
	out2.close();


}


}}}}

