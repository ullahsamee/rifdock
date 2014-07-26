#include <gtest/gtest.h>
#include <kinematics/Scene.hh>
#include <boost/make_shared.hpp>

namespace scheme {
namespace kinematics {

using std::cout;
using std::endl;

struct X1dim {
	double val_;
	X1dim() : val_(0) {}
	X1dim(double d) : val_(d) {}
	bool operator==(X1dim const & o) const { return o.val_==val_; }
};
X1dim operator*(X1dim a, X1dim b){ return X1dim(a.val_+b.val_); }
std::ostream & operator<<(std::ostream & out, X1dim const & x){ return out << x.val_; }


TEST(Scene,bodies_hold_actors)
{
	typedef ActorConcept<X1dim,int> ADI;
	typedef ActorConcept<X1dim,char> ADC;	
	typedef m::vector< ADI, ADC > Actors;
	typedef m::vector< ADI, std::pair<ADI,ADC> > Interactions;

	typedef ConformationConcept<Actors,std::vector<m::_1>,X1dim> Conformation;
	shared_ptr<Conformation> conf = make_shared<Conformation>();
	conf->add_actor( ADI(0.0,0) );
	conf->add_actor( ADI(1.0,0) );
	conf->add_actor( ADC(2.0,'0') );
	shared_ptr<Conformation> conf2 = make_shared<Conformation>();
	conf2->add_actor( ADI(0.0,1) );
	conf2->add_actor( ADI(1.0,1) );
	conf2->add_actor( ADC(2.0,'1') );

	// typedef Scene<Conformation,Interactions> Scene;
	// Scene scene;
	// scene.add_body(conf);
	// scene.add_body(conf2);

	// ASSERT_EQ( scene.body(0).get_actor<ADI>(0), ADI(0.0,0)    );
	// ASSERT_EQ( scene.body(0).get_actor<ADI>(1), ADI(1.0,0)    );
	// ASSERT_EQ( scene.body(0).get_actor<ADC>(0), ADC(2.0,'0')  );
	// ASSERT_EQ( scene.body(1).get_actor<ADI>(0), ADI(0.0,1)    );
	// ASSERT_EQ( scene.body(1).get_actor<ADI>(1), ADI(1.0,1)    );
	// ASSERT_EQ( scene.body(1).get_actor<ADC>(0), ADC(2.0,'1')  );

	// ASSERT_THROW( scene.body(2), std::out_of_range );
	// ASSERT_THROW( scene.body(1).get_actor<ADI>(2), std::out_of_range );
	// // body.set_position(1.0);
	// // cout << body.get_actor<ADI>(0) << endl;	

	// cout << scene << endl;

}

}
}
