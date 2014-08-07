#include <gtest/gtest.h>

#include "kinematics/Scene_io.hh"
#include "actor/ActorConcept_io.hh"
#include "objective/ObjectiveFunction.hh"
#include "objective/ObjectiveVisitor.hh"
#include "numeric/X1dim.hh"

#include <io/dump_pdb_atom.hh>

#include <boost/foreach.hpp>



#include <stdint.h>
#include <fstream>

#include <Eigen/Geometry>

namespace scheme { namespace kinematics { namespace test_eigen {

using std::cout;
using std::endl;
using boost::tie;

using Eigen::AngleAxis;
typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;
typedef Eigen::Vector3d Vec;
}}} namespace Eigen {
scheme::kinematics::test_eigen::Xform inverse(scheme::kinematics::test_eigen::Xform const & x){ return x.inverse(); }
std::ostream & operator<<(std::ostream & out,scheme::kinematics::test_eigen::Xform const & x){
	return out << "t("<<x.translation().transpose()<<")";
}

} namespace scheme { namespace kinematics { namespace test_eigen {


TEST(Scene_eigen,eigen_transform){
	// Xform x(Xform::Identity());
	AngleAxis<double> aa(0,Vec(1,0,0));
	Xform x(aa);
	//cout << x << endl;
	Eigen::Vector3d v(0,0,0);
	x*v;
	Eigen::RowVector3d rv(0,0,0);	
	// x*rv; // not allowed
}

typedef size_t Index;
typedef std::pair<size_t,size_t> Index2;
typedef std::pair<Index2,Index2> Index4;

struct Xactor {
	typedef Xform Position;
	Xform position_;
	int data_;
	Xactor() : position_(Xform::Identity()) {}
	Xactor(Position const & p, int d) : position_(p),data_(d) {}
	Xactor(Xactor const & a,Position const & moveby){ position_ = moveby*a.position(); data_ = a.data_; }
	void set_position( Position const & pos	){ position_ = pos; }
	Position const & position() const { return position_; }
	// bool operator==(Xactor const & o) const { return o.position_.isApprox(position_) && o.data_==data_; }
};
std::ostream & operator<<(std::ostream & out,Xactor const & x){
	return out << "Xactor( " << x.position() << " " << x.data_ << " )";
}

void dump_pdb(std::ostream & out,Xactor const & a){
	cout << a << endl;
	io::dump_pdb_atom(out, "CA", a.position() * Vec(0,0,0) );
	io::dump_pdb_atom(out, "O" , a.position() * Vec(1.5,0,0) );	
	io::dump_pdb_atom(out, "NI", a.position() * Vec(0,1.5,0) );
	io::dump_pdb_atom(out, "N" , a.position() * Vec(0,0,1.5) );		
}

////////////// test scores ////////////////////////

	struct ScoreX {
		typedef double Result;
		typedef Xactor Interaction;
		static std::string name(){ return "ScoreADI"; }
		template<class Config>
		Result operator()(Interaction const & a, Config const& ) const {
			return a.data_;
		}
	};
	std::ostream & operator<<(std::ostream & out,ScoreX const& s){ return out << s.name(); }

	double distance( Xform const & a, Xform const & b ){
		return (a.translation()-b.translation()).norm();
	}

	struct ScoreXX {
		static size_t ncalls;
		typedef double Result;
		typedef Xactor Actor1;
		typedef Xactor Actor2;
		typedef std::pair<Xactor,Xactor> Interaction;
		static std::string name(){ return "ScoreADIADI"; }
		template<class Config>
		Result operator()(Actor1 const & a1, Actor2 const & a2, Config const& ) const {
			++ncalls;
			// cout << a.first << " " << a.second << endl;
			return distance(a1.position(),a2.position());
		}
		template<class Config>
		Result operator()(Interaction const & i, Config const& c) const {
			return this->template operator()<Config>(i.first,i.second,c);
		}
	};
	size_t ScoreXX::ncalls = 0;
	std::ostream & operator<<(std::ostream & out,ScoreXX const& s){ return out << s.name(); }


	struct Config {};

/////////////////// tests //////////////////////////

TEST(Scene_eigen,basic){
	typedef	objective::ObjectiveFunction<
		m::vector<
			ScoreX,
			ScoreXX
		>,
		Config
	> ObjFun;
	typedef ObjFun::Results Results;
	ObjFun score;

	typedef m::vector< Xactor > Actors;
	typedef Conformation<Actors> Conf;
	typedef Scene<Conf,Xform,size_t> Scene;

	Scene scene(3);
	ASSERT_EQ( score(scene).sum(), 0 );

	scene.mutable_conformation_asym(0).add_actor( Xactor(Xform::Identity(),7) );
	scene.set_position(0,Xform::Identity());


	std::ofstream out("test.pdb");
	dump_pdb(out,scene);
	out.close();

}

// TEST(Scene_eigen,symmetry){
// 	typedef	objective::ObjectiveFunction<
// 		m::vector<
// 			ScoreADI,
// 			ScoreADC,
// 			ScoreADIADI,
// 			ScoreADCADI
// 		>,
// 		Config
// 	> ObjFun;
// 	typedef ObjFun::Results Results;

// 	ObjFun score;

// 	typedef m::vector< ADI, ADC > Actors;
// 	typedef Conformation<Actors> Conformation;
// 	typedef Scene<Conformation,X1dim,size_t> Scene;

// 	Scene scene(2);
// 	scene.add_symframe(10);

// 	ASSERT_EQ( score(scene).sum(), 0 );

// 	scene.mutable_conformation_asym(0).add_actor( ADI(0,1) );
// 	scene.mutable_conformation_asym(1).add_actor( ADI(1,2) );
// 	// check that symmetric interactiosn are downweighted by 0.5: 40/2 = 20
// 	ASSERT_EQ( score(scene), Results(3,0,1+20,0) ); 

// 	scene.mutable_conformation_asym(1).add_actor( ADI(0,1) );
// 	scene.mutable_conformation_asym(1).add_actor( ADI(1,2) );
// 	// check that symmetric interactiosn are downweighted by 0.5: 160/2 = 80
// 	ASSERT_EQ( score(scene), Results(6,0,2+80,0) );

// 	// scene.mutable_conformation_asym(1).add_actor( ADC(0,'1') );
// 	// scene.mutable_conformation_asym(1).add_actor( ADC(0,'2') );
// 	// ASSERT_EQ( score(scene), Results(20,3,16,80) );

// }







// TEST(Scene_eigen,DISABLED_performance){
// // TEST(Scene_eigen,performance){
// 	// TODO: speed up SceneIter iteration 
// 	//       iteration seems to take about 100 cycles per score call overhead
// 	//       much of this is probably all the conditions for symmetry checks
// 	//       could template out these and have both sym and asym scenes?
// 	// FIX with visitation pattern, seems at least 10x faster

// 	std::cout << "This test performs 301.934M score calls, should take about a second when compiled with optimizations." << std::endl;

// 	typedef	objective::ObjectiveFunction<
// 		m::vector<
// 			ScoreADI,
// 			ScoreADC,
// 			ScoreADIADI,
// 			ScoreADCADI
// 		>,
// 		Config
// 	> ObjFun;
// 	typedef ObjFun::Results Results;

// 	ObjFun score;

// 	typedef m::vector< ADI, ADC > Actors;
// 	typedef Conformation<Actors> Conformation;
// 	typedef Scene<Conformation,X1dim,uint32_t> Scene;

// 	ScoreADIADI obj;
// 	Config c;
// 	objective::ObjectiveVisitor<ScoreADIADI,Config> visitor(obj,c);

// 	Scene scene; {
// 		Scene::Index const NBOD = 10;
// 		Scene::Index const NSYM = 20;
// 		Scene::Index const NACT = 400;
// 		for(Scene::Index i = 0; i < NBOD; ++i) scene.add_body();		
// 		for(Scene::Index i = 0; i < NSYM-1; ++i) scene.add_symframe(i+1);
// 		for(Scene::Index i = 0; i < NBOD; ++i){
// 			for(Scene::Index j = 0; j < NACT; ++j){
// 				scene.mutable_conformation_asym(i).add_actor( ADI(j,i) );
// 			}
// 		}
// 	}

// 	cout << score(scene).get<ScoreADIADI>() << " " << (double)ScoreADIADI::ncalls/1000000.0 << "M" << endl;
// 	ScoreADIADI::ncalls = 0;
// 	return;


// 	if(false)
// 	{
// 			// typedef ADI Actor1;
// 			// typedef ADI Actor2;
// 			// typedef Scene::Position Position;
// 			// Scene::Index const NBOD = scene.bodies_.size();
// 			// Scene::Index const NSYM = scene.symframes_.size()+1;
// 			// for(Scene::Index i1 = 0; i1 < NBOD*NSYM; ++i1){
// 			// 	Conformation const & c1 = scene.conformation(i1);
// 			// 	Position     const & p1 =     scene.position(i1);
// 			// 	Scene::Index const NACT1 = c1.get<Actor1>().size();
// 			// 	for(Scene::Index i2 = 0; i2 < NBOD*NSYM; ++i2){
// 			// 		if( i1 >= NBOD && i2 >= NBOD ) continue;
// 			// 		if( i2 <= i1 ) continue;
// 			// 		Conformation const & c2 = scene.conformation(i2);
// 			// 		Position     const & p2 =     scene.position(i2);
// 			// 		Scene::Index const NACT2 = c2.get<Actor2>().size();					
// 			// 		for(Scene::Index j1 = 0; j1 < NACT1; ++j1){
// 			// 			Actor1 a1( c1.get<Actor1>()[j1], p1 );
// 			// 			for(Scene::Index j2 = 0; j2 < NACT2; ++j2){
// 			// 				Actor2 a2( c2.get<Actor2>()[j2], p2 );
// 			// 				visitor( a1, a2, i1<NBOD&&i2<NBOD?1.0:0.5 );
// 			// 			}
// 			// 		}
// 			// 	}
// 			// }
// 			// cout << visitor.result_ << " " << (double)ScoreADIADI::ncalls/1000000.0 << endl;
// 			// ScoreADIADI::ncalls = 0;

// 	}

// 	if(false)
// 		performance_test_helper(scene,visitor);
// 	ScoreADIADI::ncalls = 0;


// 	scene.visit(visitor);
// 	cout << visitor.result_ << " " << (double)ScoreADIADI::ncalls/1000000.0 << "M" << endl;
// 	ScoreADIADI::ncalls = 0;

// }

}
}
}
