#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
using namespace std;
#include <iomanip>
#include "lattice.h"
#include "berry_tests.h"
#include "time.h"


bool IsOdd (int i) {
    return ((i%2)==1);
}
int main(){
    void Ed_Coulomb();
    
    Ed_Coulomb();
}
void Ed_Coulomb(){
    vector<NQ> NQCE(3);
    NQCE[0].N=0; NQCE[1].N=1; NQCE[2].N=2;
    NQCE[0].Q=4.5; NQCE[1].Q=4.5; NQCE[2].Q=4.;
    
    int n=0;
    vector<NQ> NQPA(n);
//    for (int i=0; i<n; i++) {
//        NQPA[i].N=2*i+1;
//        NQPA[i].Q=4.;
//        if (NQPA[i].N>8) {
//            NQPA[i].Q=4.25;
//        }
//    }

    int ncore=20; double shift=0.25;
    for (int i=-1; i<=5; i++) {
        parallel_ce_pa(ncore, NQCE, NQPA, shift, "params_cfl", -1);
    }
}
/*
void test_exact(double shift){
    cout<<"*****shift="<<shift<<"*****"<<endl;
//    testlatticepp(shift);
    
    int n=9;
    vector<NQ> N_Q(n);
    for (int i=0; i<n; i++) {
        N_Q[i].N=i;
        N_Q[i].Q=-1.;
    }
    
    vector<double> pa;
    pa=
    pairamplitude_ExplicitLatticeSum3(3,shift,N_Q);
    for (int i=0; i<pa.size(); i++) {
        cout<<"pa["<<i<<"]="<<setprecision(10)<<pa[i]<<endl;
    }
    cout<<endl;
//    pa=
//    pairamplitude_ExplicitLatticeSum3(3,shift,N_Q);
//    for (int i=0; i<pa.size(); i++) {
//        cout<<"pa["<<i<<"]="<<setprecision(10)<<pa[i]<<endl;
//    }
//    cout<<endl;
}
*/
