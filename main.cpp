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
    vector<double> screen;
    for (int i=0; i<10; i++) {
        screen.push_back(0.5*i);
    }
    vector<NQ> NQCE(5);
    vector<NQ> NQPA(0);
    
    NQCE[0].N=0; NQCE[0].Q=10.0; NQCE[0].screen=screen;
    NQCE[1].N=1; NQCE[1].Q=5.50; NQCE[1].screen=screen;
    NQCE[2].N=2; NQCE[2].Q=5.00; NQCE[2].screen=screen;
    NQCE[3].N=2; NQCE[3].Q=4.25; NQCE[3].screen=vector<double>{0.,0.1};
    NQCE[4].N=2; NQCE[4].Q=4.00; NQCE[4].screen=vector<double>{0.,0.1};

    int ncore=10; double shift=0.25;
    for (int i=-1; i<=5; i++) {
        parallel_ce_pa(ncore, NQCE, NQPA, shift, "params_cfl", i, 0);
    }
    for (int i=-1; i<=5; i++) {
        //parallel_ce_pa(ncore, NQCE, NQPA, shift, "params_cfl2", i);
    }
}
