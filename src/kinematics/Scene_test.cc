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
};
X1dim operator*(X1dim a, X1dim b){ return X1dim(a.val_+b.val_); }
std::ostream & operator<<(std::ostream & out, X1dim const & x){ return out << x.val_; }


TEST(Scene,basic_test){
	typedef ActorConcept<X1dim,int> ADI;
	typedef ActorConcept<X1dim,char> ADC;	
	typedef m::vector< ADI, ADC > Actors;
	typedef Scene<Actors> Scene;
	Scene scene;
	cout << scene << endl;

	boost::shared_ptr<Conformation<Actors> > conf = boost::make_shared<Conformation<Actors> >();
	conf->get<ADI>().push_back( ADI(0.0,0) );
	Body<Actors,X1dim> body(conf);

	cout << ADI(0,0) << endl;
	cout << body.get_positioned_actor<ADI>(0) << endl;
	body.set_position(1.0);
	cout << body.get_positioned_actor<ADI>(0) << endl;	
}

}
}
