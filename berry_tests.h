#ifndef BERRY_TESTS_H
#define BERRY_TESTS_H

#include "utils.h"
#include "lattice.h"
#include <complex>
#include <Eigen/Core>

struct data{
    int num;
    double position[2];
    double amp[5];
    double ang[5];
    double energy;
    double ang_trace;
    double det;
    double dfromnorm;
};
struct oneconfig{
    vector<vector<int> > deletelist, addlist;
    int num;
};
struct hop{
    vector<int> list;//this is momentum list. interaction = v_{1234} c1^+ c2^+ c3 c4.
    vector<complex<double>> ele;//this is a list for hopping amplitude, for different alpha.
};
/*
//Some Utilities
void hermitianize(Eigen::MatrixXcd &x);
double phasemod(complex<double> in);

//Some Checkings, including ortho, auto-correlation
void check_orthogonality(string type);
void single_run(string filename="params", bool trace=false);
void structurefactor(string, int num_core=1);
*/

//Energetics
void parallel_ce_pa(int ncore, vector<int> PP, bool bo_shift=false, double shift=0.25, string filename="params");//energy & p.a. calculator for high LL.
/*
void parallel_energy(int ncore, string filename="params");//energy calculator for high LL.
void coul_energy(LATTICE& lattice, int nWarmup, int nMeas, int nSteps, int nBins, string filename);
void plot_CFL_coule_vsdbar(int grid, int Ne, int nMeas, int nBins);
void CFL_ne5_energy_var(int nMeas, int nBins, int num_core=1);
void findstate();

//Laughlin-Hole Berry Phase
void laughlinberryphase(string input_name, string output_name, vector<double> length, double steplength, int change_nMeas, int change_Ne, int num_core, double theta=0.5*M_PI, double alpha=1.0);
void laughlin_bp_single_state(int gs, vector<double> length, double steplength, vector<data> &datas);
void two_holes(string input_name, string str, int nmeasurement, data& test);
void test_error(int ne, double loop, double steplength, int nMea, int ncore, string test, int num_core);

//CFL Berry Phase
void CFL_berry_phases_parallel(string params_name, string output_name, int num_core, string kind, double theta=0.5*M_PI, double alpha=1.0);
void get_dlist(string holes, int tempNe, string kind, vector< vector<int> > &extra_ds, vector< vector<int> > &remove_ds, vector< vector<int> > &old_ds);

//Particle-Hole Symmetry
void ParticleHoleSym();
void ParticleHoleSymBackwards();
void ParticleHoleSym2();
void Explicit();

*/

//pairamplitude (testings)
inline double Laguerrel(int N, double x);
complex<double> interaction(int m1, int m3, int m4, int No, vector<double> vpseu, string type="pa");//this gives 'ED interaction matrix' for square torus.
complex<double> latticepp(LATTICE ll, int m1, int m2, int m3, int m4, string type);//this gives 'sigma function lattice sum interaction matrix'.
void testlatticepp(double shift=0.);//This function tested that the matrix element V1234 got from 'sigma function lattice sum' is correct.

//pairamplitude
void LatticeSumHoplist(string);//generalize hoplist from lattice summation.
void pairamplitude_MC(string filename, bool trace, int num_core, vector<int> PP);
void pairamplitude_MC2(string filename, bool trace, int num_core, vector<int> PP, vector<double> Q);
//2 and 3 particle explicit lattice sum for pair-amplitude.
void pairamplitude_ExplicitLatticeSum2(int invNu);
void pairamplitude_ExplicitLatticeSum3(int invNu);
vector<double> pairamplitude_ExplicitLatticeSum2(int invNu, double shift1, double shift2, vector<int> PP);

void onebody(int m1, int m2, int NPhi, double shift, int lat_scale);

//old functions
void pairamplitudeold(string filename, bool trace, int num_core, bool pseu, bool mc);

#endif
