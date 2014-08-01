#include <gtest/gtest.h>

#include "kinematics/Scene_io.hh"
#include "actor/ActorConcept_io.hh"
#include "objective/ObjectiveFunction.hh"
#include "numeric/X1dim.hh"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/fusion/include/io.hpp>


namespace scheme { namespace kinematics { namespace test {

using std::cout;
using std::endl;
using boost::tie;
using numeric::X1dim;

typedef actor::ActorConcept<X1dim,int> ADI;
typedef actor::ActorConcept<X1dim,char> ADC;	
typedef size_t Index;
typedef std::pair<size_t,size_t> Index2;
typedef std::pair<Index2,Index2> Index4;

////////////// test scores ////////////////////////

	struct ScoreADI {
		typedef double Result;
		typedef ADI Interaction;
		static std::string name(){ return "ScoreADI"; }
		template<class Config>
		Result operator()(Interaction const & a, Config const& ) const {
			return a.data_;
		}
	};
	std::ostream & operator<<(std::ostream & out,ScoreADI const& s){ return out << s.name(); }

	struct ScoreADIADI {
		typedef double Result;
		typedef std::pair<ADI,ADI> Interaction;
		static std::string name(){ return "ScoreADIADI"; }
		template<class Config>
		Result operator()(Interaction const & a, Config const& ) const {
			// cout << a.first << " " << a.second << endl;
			return distance(a.first.position(),a.second.position());
		}
	};
	std::ostream & operator<<(std::ostream & out,ScoreADIADI const& s){ return out << s.name(); }

	struct ScoreADC {
		typedef double Result;
		typedef ADC Interaction;
		static std::string name(){ return "ScoreADC"; }
		template<class Config>
		Result operator()(Interaction const & a, Config const& ) const {
			switch(a.data_){
				case '1': return 1; break;
				case '2': return 2; break;
				case '3': return 3; break;
				case '4': return 4; break;
			};
			return 0;
		}
	};
	std::ostream & operator<<(std::ostream & out,ScoreADC const& s){ return out << s.name(); }

	struct ScoreADCADI {
		typedef double Result;
		typedef std::pair<ADC,ADI> Interaction;
		static std::string name(){ return "ScoreADCADI"; }
		template<class Config>
		Result operator()(Interaction const & a, Config const& ) const {
			return 10*distance(a.first.position(),a.second.position());
		}
	};
	std::ostream & operator<<(std::ostream & out,ScoreADCADI const& s){ return out << s.name(); }

	struct Config {};

/////////////////// tests //////////////////////////

TEST(SceneObjective,basic){
	typedef	objective::ObjectiveFunction<
		m::vector<
			ScoreADI,
			ScoreADC,
			ScoreADIADI,
			ScoreADCADI
		>,
		Config
	> ObjFun;
	typedef ObjFun::Results Results;
	ObjFun score;

	typedef m::vector< ADI, ADC > Actors;
	typedef Conformation<Actors> Conformation;
	typedef Scene<Conformation,X1dim,size_t> Scene;

	Scene scene(3);
	ASSERT_EQ( score(scene).sum(), 0 );

	scene.mutable_conformation_asym(0).add_actor( ADI(0,1) );
	scene.mutable_conformation_asym(0).add_actor( ADI(0,2) );
	scene.mutable_conformation_asym(0).add_actor( ADI(0,3) );
	scene.mutable_conformation_asym(0).add_actor( ADI(0,4) );
	ASSERT_EQ( score(scene), Results(10,0,0,0) );

	scene.mutable_conformation_asym(1).add_actor( ADI(1,1) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,2) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,3) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,4) );
	ASSERT_EQ( score(scene), Results(20,0,16,0) );

	scene.mutable_conformation_asym(0).add_actor( ADC(0,'1') );
	scene.mutable_conformation_asym(0).add_actor( ADC(0,'2') );
	ASSERT_EQ( score(scene), Results(20,3,16,80) );

	scene.mutable_conformation_asym(2).add_actor( ADI(2,100) );
	ASSERT_EQ( score(scene), Results(120,3,28,120) );

	scene.mutable_conformation_asym(2).add_actor( ADI(2,1000) );
	ASSERT_EQ( score(scene), Results(1120,3,40,160) );

}

TEST(SceneObjective,symmetry){
	typedef	objective::ObjectiveFunction<
		m::vector<
			ScoreADI,
			ScoreADC,
			ScoreADIADI,
			ScoreADCADI
		>,
		Config
	> ObjFun;
	typedef ObjFun::Results Results;

	ObjFun score;

	typedef m::vector< ADI, ADC > Actors;
	typedef Conformation<Actors> Conformation;
	typedef Scene<Conformation,X1dim,size_t> Scene;

	Scene scene(2);
	scene.add_symframe(10);

	ASSERT_EQ( score(scene).sum(), 0 );

	scene.mutable_conformation_asym(0).add_actor( ADI(0,1) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,2) );
	// check that symmetric interactiosn are downweighted by 0.5: 40/2 = 20
	ASSERT_EQ( score(scene), Results(3,0,1+20,0) ); 

	scene.mutable_conformation_asym(1).add_actor( ADI(0,1) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,2) );
	// check that symmetric interactiosn are downweighted by 0.5: 160/2 = 80
	ASSERT_EQ( score(scene), Results(6,0,2+80,0) );

	// scene.mutable_conformation_asym(1).add_actor( ADC(0,'1') );
	// scene.mutable_conformation_asym(1).add_actor( ADC(0,'2') );
	// ASSERT_EQ( score(scene), Results(20,3,16,80) );

}


}
}
}
