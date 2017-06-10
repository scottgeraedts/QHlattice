#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
using namespace std;
#include <iomanip>
#include "lattice.h"
#include "berry_tests.h"


bool IsOdd (int i) {
    return ((i%2)==1);
}
int main(){
    parallel_ce_pa(2, vector<int>{0,1,2,3,4,5,6,7}, "params");
    
//    cout<<"pairamplitude2"<<endl;
//    pairamplitude_ExplicitLatticeSum2(5);
//    cout<<"pairamplitude3"<<endl;
//    pairamplitude_ExplicitLatticeSum3(5);
    
}
