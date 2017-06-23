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
    Ed_Coulomb();
}
void Ed_Coulomb(){
    parallel_ce_pa(2, vector<int>{}, "params_cfl");
}
void test_exact(){
    cout<<"pairamplitude2"<<endl;
    pairamplitude_ExplicitLatticeSum2(3);
    int m1=2, m2=7, NPhi=5, n=0;
    onebody(m1, m2, NPhi, n);
    cout<<"pairamplitude3"<<endl;
    pairamplitude_ExplicitLatticeSum3(5);
}
void test_int_bc(){
    vector<int> PP = vector<int>{0,1,2,3,4,5,6,7};
    vector<double> value (PP.size(), 0.);
    vector<double> tmp;
    
    int invNu=3, Ne=2, NPhi=Ne*invNu;
    int N=10;
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            double x,y;
            
            if (i==0&&j==0&&N==0) {
                x=0.; y=0.;
            }
            else {
                x=i/(1.*N*NPhi);
                y=j/(1.*N*NPhi);
            }
            
            tmp=pairamplitude_ExplicitLatticeSum2(invNu, x, y, PP);
            for (int k=0; k<value.size(); k++) {
                value[k]+=tmp[k];
            }
        }
    }
    
    double x=0., y=0;
    tmp=pairamplitude_ExplicitLatticeSum2(invNu, x, y, PP);
    for (int k=0; k<value.size(); k++) {
        value[k]+=tmp[k];
    }
    
    for (int i=0; i<value.size(); i++) {
        if (N==0)
            cout<<"pairamplitude["<<i<<"] = "<<setprecision(15)<<value[i]<<endl;
        else
            cout<<"pairamplitude["<<i<<"] = "<<setprecision(15)<<value[i]/(1.*N*N)<<endl;
    }
}

