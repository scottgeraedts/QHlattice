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
    
    clock_t t1,t2,t3,t4;
    
//    t1=clock();
//    parallel_ce_pa(2, NQCE, NQPA, false, 0., "params_cfl00");//Ne=37, 50
//    t2=clock();
//    parallel_ce_pa(2, NQCE, NQPA, false, 0., "params_cfl10");//Ne=69, 50
//    t3=clock();
//    parallel_ce_pa(2, NQCE, NQPA, false, 0., "params_cfl20");//Ne=97, 50
//    t4=clock();
//    cout<<"\nTIME="<<((float)t2-(float)t1)/CLOCKS_PER_SEC<<" "<<((float)t3-(float)t2)/CLOCKS_PER_SEC<<" "<<((float)t4-(float)t3)/CLOCKS_PER_SEC<<endl;
//    
//    t1=clock();
//    parallel_ce_pa(2, NQCE, NQPA, false, 0., "params_cfl01");//Ne=37, 200
//    t2=clock();
//    parallel_ce_pa(2, NQCE, NQPA, false, 0., "params_cfl11");//Ne=69, 200
//    t3=clock();
//    parallel_ce_pa(2, NQCE, NQPA, false, 0., "params_cfl21");//Ne=97, 200
//    t4=clock();
//    cout<<"\nTIME="<<((float)t2-(float)t1)/CLOCKS_PER_SEC<<" "<<((float)t3-(float)t2)/CLOCKS_PER_SEC<<" "<<((float)t4-(float)t3)/CLOCKS_PER_SEC<<endl;
    

    t1=clock();
    pomeranchuk_instability(2, NQCE, "params_cfl", vector<double>{1,1,1});
    t2=clock();
    cout<<"\nTIME="<<((float)t2-(float)t1)/CLOCKS_PER_SEC<<endl;
    
//    NQCE[0].Q=10.; NQCE[1].Q=10.; NQCE[2].Q=10.;
//    t1=clock();
//    pomeranchuk_instability(2, NQCE, "params_cfll", vector<double>{1,1});
//    t2=clock();
//    cout<<"\nTIME="<<((float)t2-(float)t1)/CLOCKS_PER_SEC<<endl;
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
