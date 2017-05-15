#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
using namespace std;
#include <iomanip>
#include "lattice.h"
#include "berry_tests.h"

bool IsOdd (int i) {
    return ((i%2)==1);
}
int main(){
//    double theta=0.5*M_PI, alpha=1.0;
//    int num_core=1;

//	ifstream infile("params");
//	int Ne, invNu, nWarmups, nBins, nMeas, nSteps, seed;
//	infile>>Ne>>invNu;
//	infile>>nWarmups>>nMeas>>nSteps>>nBins;
//	infile>>seed;
//	LATTICE_PARAMS params(Ne);
//	infile>>params.testing;

//	for(int i=0;i<2;i++){
//		params.w_delta=0.;
//		LATTICE ll(params);
//		coul_energy(ll, nWarmups, nMeas, nSteps, nBins, "out");
//		params.w_delta+=0.5;
//	}

//	ParticleHoleSym2();
//	Explicit();
//    double theta=0.333333333333*M_PI, alpha=1.0;
    
//    CFL_berry_phases_parallel("params_e", "ne69", num_core, "loop1", theta, alpha);
    
//    pairamplitude_MC2("paramsne", false, num_core, vector<int>{1,3,5,7,9}, vector<double>{4,4,4,4,4});
//    pairamplitude_MC("paramsne", false, num_core, vector<int>{1,3,5,7,9});
    
    
//    data test;
//    two_holes("params_lau", "", 5000, test);
//    vector<data> datas;
    
    
//    //the laughlin-hole part
//    string inputfile="params_lau5";
//    double leg=0.3;
//    
//    //data test;
//    //two_holes(inputfile ,"", 5000, test);
//    
//    vector<data> datas;
//    for (int i=5; i<6; i++) {
//        double steplength=0.05;
//        int ncore=5;
//        cout<<"stplength="<<steplength<<" length="<<leg<<"*"<<0.1*i<<" (i="<<i<<")."<<endl;
//        laughlinberryphase(inputfile ,"", vector<double>{leg, 0.1*i}, steplength, 0, 0, ncore);
//    }
//    //the laughlin-hole part ends.
    
    
    single_run();
    
//    LatticeSumHoplist("paramsne");
//
//    pairamplitude_ExplicitLatticeSum3();
    
//    pairamplitudeold("paramsne", false, 1, true, false);
}