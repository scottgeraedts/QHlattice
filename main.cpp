#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
using namespace std;
#include <iomanip>
#include "lattice.h"
#include "berry_tests.h"


bool IsOdd (int i) {
    return ((i%2)==1);
}
int main(){
//    parallel_energy(2, "params");
    
    //pair-amplitude.
//    pairamplitude_MC("params", false, 2, vector<int>{1,3,15,19});
//    parallel_ce_pa(2, vector<int>{1,3,15,19}, "params");
    cout<<"pairamplitude2"<<endl;
    pairamplitude_ExplicitLatticeSum2();
//    cout<<"pairamplitude3"<<endl;
//    pairamplitude_ExplicitLatticeSum3();
//    
    
}
