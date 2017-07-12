#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
using namespace std;
#include <iomanip>
#include "lattice.h"
#include "berry_tests.h"


bool IsOdd (int i) {
    return ((i%2)==1);
}
int main(){
    void Ed_Coulomb();
    void test_exact(double);
    
//    test_exact(0.25);
    Ed_Coulomb();
//    test_exact(0.00);
//    test_exact(0.25);
//    test_exact(0.30);
//    test_exact(0.50);
}
void Ed_Coulomb(){
    vector<NQ> NQCE(2);
    NQCE[0].N=0; NQCE[1].N=2;
    NQCE[0].Q=4.75; NQCE[1].Q=4.75;
    
    int n=4;
    vector<NQ> NQPA(n);
    for (int i=0; i<n; i++) {
        NQPA[i].N=2*i+1;
        NQPA[i].Q=4.;
        if (NQPA[i].N>8) {
            NQPA[i].Q=4.25;
        }
    }
    
    parallel_ce_pa(10, NQCE, NQPA, true, 0.25, "params_cfl");
}
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
