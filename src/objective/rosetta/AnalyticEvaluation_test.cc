#include <gtest/gtest.h>

#include "objective/rosetta/AnalyticEvaluation.hh"

#include "objective/rosetta/EtableParamsOnePair_init.hh"


namespace objective {
namespace rosetta {

using std::cout;
using std::endl;

TEST(RosettaEtable,init_params){
	std::vector<EtableParamsOnePair> params;
	init_EtableParamsOnePair(params);
	cout << params.size() << endl;

	double sol,atr,rep;
	double dis,dis2,inv_dis2;
	
	dis = 3.3; dis2 = dis*dis; inv_dis2 = 1.0/dis2;
	lk_evaluation( params[0], 3.0, 1/9.0, sol ); 
	lj_evaluation( params[0], dis, dis2, inv_dis2, atr, rep);
	cout << atr << " " << rep << " " << sol << endl;


}


}
}

