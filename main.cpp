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
    
    Ed_Coulomb();
//    test_exact(0.00);
//    test_exact(0.25);
}
void Ed_Coulomb(){
    vector<NQ> NQCE(2);
    NQCE[0].N=0; NQCE[1].N=2;
    NQCE[0].Q=4.75; NQCE[1].Q=4.75;
    
    int n=8;
    vector<NQ> NQPA(n);
    for (int i=0; i<n; i++) {
        NQPA[i].N=i;
        NQPA[i].Q=5.;
        if (i==6) NQPA[i].Q=5.25;
        if (i==7) NQPA[i].Q=5.50;
    }
    
    parallel_ce_pa(10, NQCE, NQPA, true, 0.25, "paramscfl");
}
void test_exact(double shift){
    cout<<"*****shift="<<shift<<"*****"<<endl;
    testlatticepp(shift);
    
    vector<NQ> N_Q(8);
    for (int i=0; i<8; i++) {
        N_Q[i].N=i;
        N_Q[i].Q=-1.;
    }
    
    vector<double> pa;
    pa=
    pairamplitude_ExplicitLatticeSum2(3,shift,N_Q);
    for (int i=0; i<pa.size(); i++) {
        cout<<"pa["<<i<<"]="<<setprecision(10)<<pa[i]<<endl;
    }
    cout<<endl;
    pa=
    pairamplitude_ExplicitLatticeSum3(3,shift,N_Q);
    for (int i=0; i<pa.size(); i++) {
        cout<<"pa["<<i<<"]="<<setprecision(10)<<pa[i]<<endl;
    }
    cout<<endl;
}
