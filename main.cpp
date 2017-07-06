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
    void testonebody();
    
//    Ed_Coulomb();
    test_exact(0.00);
    test_exact(0.25);
}
void Ed_Coulomb(){
    parallel_ce_pa(10, vector<int>{0,2} ,vector<int>{1,3}, true, 0., "params_cfl");
}
void test_exact(double shift){
//    for (double shift=0.; shift<=0.5; shift+=0.25) {
//        cout<<"***** shift="<<shift<<" *****"<<endl;
//        testlatticepp(shift);
//        cout<<"**********"<<endl<<endl;
//    }
    
    cout<<"*****shift="<<shift<<"*****"<<endl;
    testlatticepp(shift);
    
//    TODO:generalized bc works for 2 particle, but not for 3 particle. don't know y.
    
//    vector<double> pa;
//    pa=
//    pairamplitude_ExplicitLatticeSum2(3,shift,shift,vector<int>{0,1,2,3,4,5,6,7});
//    for (int i=0; i<pa.size(); i++) {
//        cout<<"pa["<<i<<"]="<<setprecision(10)<<pa[i]<<endl;
//    }
//    cout<<endl;
//    pa=
//    pairamplitude_ExplicitLatticeSum3(3,shift,shift,vector<int>{0,1,2,3,4,5,6,7});
//    for (int i=0; i<pa.size(); i++) {
//        cout<<"pa["<<i<<"]="<<setprecision(10)<<pa[i]<<endl;
//    }
//    cout<<endl;
    
}
/*
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
void testonebody(){
    ifstream inf("para2");
    int Nphi, m1, m2, m3, m4, m5, m6, m7, m8, Kx, Ky, round, lat_scale;
    double Kq;
    complex<double> tmp=0.;
    
    inf>>Nphi;
    inf>>m1>>m2>>m3>>m4;
    inf>>Kx>>Ky>>Kq;
    inf>>round;
    inf>>lat_scale;
    
    double shift;
    
    cout<<"test one body"<<endl;
    shift=0.;
    cout<<"*****shift="<<shift<<"*****"<<endl;
    onebody(m1, m2, Nphi, shift, lat_scale);
    shift=0.25;
    cout<<"*****shift="<<shift<<"*****"<<endl;
    onebody(m1, m2, Nphi, shift, lat_scale);
    shift=0.5;
    cout<<"*****shift="<<shift<<"*****"<<endl;
    onebody(m1, m2, Nphi, shift, lat_scale);
    
}
*/
