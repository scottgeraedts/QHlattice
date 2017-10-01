#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
using namespace std;
#include <iomanip>
#include "lattice.h"
#include "berry_tests.h"
#include "ctime"

bool IsOdd (int i) {
    return ((i%2)==1);
}
int main(){
//    void Ed_Coulomb();
//    Ed_Coulomb();
    
    void Ed_Coulomb2();
//    Ed_Coulomb2();
    
    void sq3();
//    sq3();
    
    void testins();
//    testins();
    
//    time_t t1, t2;
//    t1 = time(NULL);
//    t2 = time(NULL);
//    cout<<"time eclips="<<(double)t2-(double)t1<<endl;    
    
    void quarterbp_e();
    void quarterbp_h();
    quarterbp_e();
    quarterbp_h();
}
void quarterbp_e(){
    int ncore=1;
    double theta=0.5*M_PI, alpha=1.0;
    string outputfile="4bp/quater_berry";
    CFL_berry_phases_parallel("params12e", outputfile, ncore, "Scott1", theta, alpha);
//    CFL_berry_phases_parallel("params12e", outputfile, ncore, "Scott2", theta, alpha);
//    CFL_berry_phases_parallel("params12e", outputfile, ncore, "Scott3", theta, alpha);
//    CFL_berry_phases_parallel("params12e", outputfile, ncore, "Scott4", theta, alpha);
//    CFL_berry_phases_parallel("params12e", outputfile, ncore, "Scott5", theta, alpha);
}
void quarterbp_h(){
    int ncore=1;
    double theta=0.5*M_PI, alpha=1.0;
    string outputfile="4bp/quater_berry";
    CFL_berry_phases_parallel("params12h", outputfile, ncore, "loop1_new", theta, alpha);
//    CFL_berry_phases_parallel("params12h", outputfile, ncore, "loop2_new", theta, alpha);
//    CFL_berry_phases_parallel("params12h", outputfile, ncore, "loop3_new", theta, alpha);
}
void Ed_Coulomb(){
    vector<double> screen;
    for (int i=0; i<100; i++) {
        screen.push_back(0.05*i);
    }
    vector<NQ> NQCE(3);
    vector<NQ> NQPA(0);
    
    NQCE[0].N=0; NQCE[0].Q=4.50; NQCE[0].screen=screen;
    NQCE[1].N=1; NQCE[1].Q=4.25; NQCE[1].screen=screen;
    NQCE[2].N=2; NQCE[2].Q=4.00; NQCE[2].screen=screen;
    
    int ncore=10; double shift=0.2;
    for (int i=7; i>=0; i--) {
        int num=0;
        parallel_ce_pa(ncore, NQCE, NQPA, shift, "params_cfl", i, num*ncore);
    }
    
}
void Ed_Coulomb2(){
    vector<double> screen;
    for (int i=0; i<100; i++) {
        screen.push_back(0.05*i);
    }
    vector<NQ> NQCE(3);
    vector<NQ> NQPA(0);
    NQCE[0].N=0; NQCE[0].Q=4.50; NQCE[0].screen=screen;
    NQCE[1].N=1; NQCE[1].Q=4.25; NQCE[1].screen=screen;
    NQCE[2].N=2; NQCE[2].Q=4.00; NQCE[2].screen=screen;
    
    int ncore=10; double shift=0.2;
    int N=8;
    
    for (int i=-1; i<N; i++) {
        int num=0;
        //parallel_ce_pa(ncore, NQCE, NQPA, shift, "params_cfl", i, num*ncore);
    }
    
    vector<double> a(N,1); a.push_back(N);
    int num=26;
    pomeranchuk_instability(ncore, NQCE, "params_cfl", a, num);
    
//    int num=0;
//    parallel_ce_pa(ncore, NQCE, NQPA, shift, "params_cfl", -1, num*ncore);
//    parallel_ce_pa(ncore, NQCE, NQPA, shift, "params_cfl", 0, num*ncore);
//    parallel_ce_pa(ncore, NQCE, NQPA, shift, "params_cfl", 1, num*ncore);
}
void sq3(){
    int num=10;
    structurefactor("params_squ", num);
//    structurefactor("params_hex", num);
//    structurefactor("params_a5a3", num);
//    structurefactor("params_a5a15", num);
//    structurefactor("params_a6a1", num);
}
void testins(){
    vector<double> screen;
    screen.push_back(0.05*50);
    
    vector<NQ> NQCE(3);
    vector<NQ> NQPA(0);
    NQCE[0].N=0; NQCE[0].Q=4.50; NQCE[0].screen=screen;
    NQCE[1].N=1; NQCE[1].Q=4.25; NQCE[1].screen=screen;
    NQCE[2].N=2; NQCE[2].Q=4.00; NQCE[2].screen=screen;
    
    int ncore=10; double shift=0.2;
    int N=1;
    
    vector<double> a(N,1);
    for (int i=1; i<=8; i++) {
        double tmp=0.25*i;
        vector<double> b=a; b.push_back(tmp);
        pomeranchuk_instability(ncore, NQCE, "params_cfltestins", b, i*ncore);
    }
}
