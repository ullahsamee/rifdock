#include <gtest/gtest.h>
#include <kinematics/Scene.hh>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include "boost/tuple/tuple.hpp"




namespace scheme { namespace kinematics {

using std::cout;
using std::endl;
using boost::tie;

struct X1dim {
	double val_;
	X1dim() : val_(0) {}
	X1dim(double d) : val_(d) {}
	bool operator==(X1dim const & o) const { return o.val_==val_; }
};
X1dim operator*(X1dim a, X1dim b){ return X1dim(a.val_+b.val_); }
std::ostream & operator<<(std::ostream & out, X1dim const & x){ return out << x.val_; }

	namespace TEST_BOUNDS {
		typedef ActorConcept<X1dim,int> ADI;
		typedef ActorConcept<X1dim,char> ADC;	
		typedef m::vector< ADI, ADC > Actors;

		typedef Conformation<Actors> Conformation;
		typedef Scene<Conformation,X1dim,size_t> Scene;
		typedef size_t Index;
		typedef std::pair<size_t,size_t> Index2;
		typedef std::pair<Index2,Index2> Index4;

		typedef std::pair<ADI,ADC> Interaction;
		typedef SceneIter<Scene,Interaction > Iter2;
		typedef get_placeholder_type<Interaction,Index>::type PH;

		TEST(Scene,interactions_empty_onebody){
			Scene scene;
			SceneIter<Scene,ADI> beg,end;
			tie(beg,end) = scene.get_interactions<ADI>();
			cout << *beg << " / " << *end << endl;
			++beg;
			cout << *beg << " / " << *end << endl;
			++beg;
			cout << *beg << " / " << *end << endl;
			++beg;
			cout << *beg << " / " << *end << endl;			
			++beg;

		}

	}

// TEST(Scene,bodies_hold_actors)
// {
// 	using std::make_pair;

// 	typedef ActorConcept<X1dim,int> ADI;
// 	typedef ActorConcept<X1dim,char> ADC;	
// 	typedef m::vector< ADI, ADC > Actors;
// 	typedef m::vector< ADI, std::pair<ADI,ADC> > Interactions;

// 	typedef Conformation<Actors> Conformation;
// 	typedef Scene<Conformation,X1dim,size_t> Scene;
// 	typedef size_t Index;
// 	typedef std::pair<size_t,size_t> Index2;
// 	typedef std::pair<Index2,Index2> Index4;

// 	typedef std::pair<ADI,ADC> Interaction;
// 	typedef SceneIter<Scene,Interaction > Iter2;
// 	typedef get_placeholder_type<Interaction,Index>::type PH;

// 	shared_ptr<Conformation> conf = make_shared<Conformation>();
// 	conf->add_actor( ADI(0.0,0) );
// 	// conf->add_actor( ADI(0.0,1) );
// 	// conf->add_actor( ADI(0.0,2) );
// 	// conf->add_actor( ADC(0.0,'0') );
// 	shared_ptr<Conformation> conf2 = make_shared<Conformation>();
// 	conf2->add_actor( ADI(1.0,0) );
// 	// conf2->add_actor( ADI(1.0,1) );
// 	conf2->add_actor( ADC(1.0,'0') );
// 	// conf2->add_actor( ADC(1.0,'1') );

// 	Scene scene;
// 	scene.add_body(conf);
// 	scene.add_body(conf2);

// 	SceneIter<Scene,ADI> beg,end;
// 	tie(beg,end) = scene.get_interactions<ADI>();
// 	cout << *beg++ << " / " << *end << endl;
// 	cout << *beg++ << " / " << *end << endl;
// 	cout << *beg++ << " / " << *end << endl;
// 	cout << *beg++ << " / " << *end << endl;			
// 	// BOOST_FOREACH( Index2 ip, iters1 ) cout << ip << endl;

// 	// {
// 	// 	std::pair<Iter2,Iter2> iters2 = scene.get_interactions<Interaction>();
// 	// 	BOOST_FOREACH( Index4 ip, iters2 ){
// 	// 		Interaction::first_type a1;
// 	// 		Interaction::second_type a2;
// 	// 		boost::tie(a1,a2) = scene.get_interaction<Interaction>(ip);
// 	// 		cout << ip.first.first << " " << ip.second.first << " " << a1 << "    ";
// 	// 		cout << ip.first.second << " " << ip.second.second << " " << a2 << endl;
// 	// 	}
// 	// }

// 	// {
// 	// 	SceneIter<Scene,Interaction > beg,end;
// 	// 	boost::tie(beg,end) = scene.get_interactions<Interaction>();
// 	// 	ASSERT_EQ( *beg, PH(make_pair(make_pair(0,0),make_pair(0,0))) ); ++beg;
// 	// }

// 	// ASSERT_EQ( scene.body(0).get_actor<ADI>(0), ADI(0.0,0)    );
// 	// ASSERT_EQ( scene.body(0).get_actor<ADI>(1), ADI(1.0,0)    );
// 	// ASSERT_EQ( scene.body(0).get_actor<ADC>(0), ADC(2.0,'0')  );
// 	// ASSERT_EQ( scene.body(1).get_actor<ADI>(0), ADI(0.0,1)    );
// 	// ASSERT_EQ( scene.body(1).get_actor<ADI>(1), ADI(1.0,1)    );
// 	// ASSERT_EQ( scene.body(1).get_actor<ADC>(0), ADC(2.0,'1')  );

// 	// ASSERT_THROW( scene.body(2), std::out_of_range );
// 	// ASSERT_THROW( scene.body(1).get_actor<ADI>(2), std::out_of_range );
// 	// // body.set_position(1.0);
// 	// // cout << body.get_actor<ADI>(0) << endl;	

// 	// cout << scene << endl;

// }

}
}
