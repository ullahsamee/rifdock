#include <gtest/gtest.h>
#include <objective/ObjectiveFunction.hh>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/size.hpp>
#include <boost/static_assert.hpp>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/io.hpp>


namespace scheme {
namespace objective {

namespace mpl = boost::mpl;
namespace bf = boost::fusion;
using std::cout;
using std::endl;

struct ScoreInt {
	typedef double Result;
	typedef int Petal1;
	typedef int Petals;
	static std::string name(){ return "ScoreInt"; }
	template<class Config>
	void operator()(int a, Result & result, Config const& c) const {
		result += a*c.scale;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreInt const& si){ return out << si.name(); }

struct ScoreInt2 {
	typedef double Result;
	typedef int Petal1;
	typedef int Petals;
	static std::string name(){ return "ScoreInt2"; }	
	template<class Config>
	void operator()(int a, Result & result, Config const& c) const {
		result += 2*a*c.scale;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreInt2 const& si){ return out << si.name(); }

struct ScoreDouble {
	typedef double Result;
	typedef double Petal1;
	typedef double Petals;
	static std::string name(){ return "ScoreDouble"; }	
	template<class Config>
	void operator()(double a, Result & result, Config const& c) const {
		result += a*c.scale;
	}
};
std::ostream & operator<<(std::ostream & out,ScoreDouble const& si){ return out << si.name(); }

struct ScoreIntDouble {
	typedef double Result;
	typedef int Petal1;
	typedef double Petal2;
	typedef std::pair<int,double> Petals;
	static std::string name(){ return "ScoreIntDouble"; }	
	template<class Config>
	void operator()(int i,double a, Result & result, Config const& c) const {
		result += i*a*c.scale;
	}
	template<class Config>
	void operator()(Petals const & p, Result & result, Config const& c) const {
		this->operator()(p.first,p.second,result,c);
	}
};
std::ostream & operator<<(std::ostream & out,ScoreIntDouble const& si){ return out << si.name(); }

struct ConfigTest {
	double scale;
	ConfigTest():scale(1){}
};

TEST(ObjectiveFunction,basic_tests){
	typedef	ObjectiveFunction<
		mpl::list<
			ScoreInt,
			ScoreInt2,
			ScoreDouble,
			ScoreIntDouble
		>,
		ConfigTest
	> ObjFun;
	ObjFun score;
	typedef ObjFun::Results Results;
	typedef util::meta::InstanceMap< mpl::vector<int,double,std::pair<int,double> >, std::vector<mpl::_1> > PetalSource;
	PetalSource petal_source;
	ASSERT_EQ( Results(0,0,0,0), score(petal_source) );
	petal_source.get<int>().push_back(1);
	ASSERT_EQ( Results(1,2,0,0), score(petal_source) );
	petal_source.get<double>().push_back(1.234);
	petal_source.get<double>().push_back(2);	
	petal_source.get<double>().push_back(-2);
	ASSERT_EQ( Results(1,2,1.234,0), score(petal_source) );
	petal_source.get<std::pair<int,double> >().push_back(std::make_pair(1,1.5));
	ASSERT_EQ( Results(1,2,1.234,1.5), score(petal_source) );

	score.config_.scale = 2.0;
	ASSERT_EQ( Results(2,4,2.468,3.0), score(petal_source) );

	// cout << score << endl;
}


TEST(ObjectiveFunction,test_results){
	typedef	ObjectiveFunction< mpl::vector< ScoreDouble, ScoreInt >, ConfigTest	> ObjFun;
	ObjFun score;
	typedef util::meta::InstanceMap< mpl::vector<int,double>, std::vector<mpl::_1> > PetalSource;
	PetalSource petal_source;
	petal_source.get<int>().push_back(1);
	petal_source.get<int>().push_back(7);	
	petal_source.get<double>().push_back(1.2345);	
	ObjFun::Results weights(2.0);
	ObjFun::Results results = score(petal_source);
	float tot = results;
	ASSERT_EQ( 9.2345f, tot );
	ASSERT_EQ( 9.2345, results );		
	ASSERT_EQ( 18.469, (double)(results*weights) );
	ASSERT_EQ( 18.469, (double)(results+results) );
	ASSERT_EQ( 27.7035, (double)(results*weights + results) );
	ASSERT_EQ( 9.2345, (double)(results*weights - results) );
	ASSERT_EQ( 9.2345, (double)(results*weights/weights) );
}

}
}


