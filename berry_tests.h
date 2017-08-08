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

//Some Utilities
void hermitianize(Eigen::MatrixXcd &x);
double phasemod(complex<double> in);
inline double Laguerrel(int N, double x);

//Some Checkings, including ortho, auto-correlation
void single_run(string filename="params", bool trace=false);
void structurefactor(string, int num_core=1);

//Energetics
void parallel_ce_pa(int ncore, vector<NQ> CE, vector<NQ> PP, double shift=0.25, string filename="params", int ind=-1);//energy & p.a. calculator for high LL.

//Pomeranchuk Instability.
void print_d(vector<vector<int>>);
vector<vector<vector<int>>> output_dset(int);
void pomeranchuk_instability(int ncore, vector<NQ> CE, string filename, vector<double> a);

//Laughlin-Hole Berry Phase
void laughlinberryphase(string input_name, string output_name, vector<double> length, double steplength, int change_nMeas, int change_Ne, int num_core, double theta=0.5*M_PI, double alpha=1.0);
void laughlin_bp_single_state(int gs, vector<double> length, double steplength, vector<data> &datas);
void two_holes(string input_name, string str, int nmeasurement, data& test);

//CFL Berry Phase
void CFL_berry_phases_parallel(string params_name, string output_name, int num_core, string kind, double theta=0.5*M_PI, double alpha=1.0);
void get_dlist(string holes, int tempNe, string kind, vector< vector<int> > &extra_ds, vector< vector<int> > &remove_ds, vector< vector<int> > &old_ds);

//Particle-Hole Symmetry
void ParticleHoleSym();
void ParticleHoleSymBackwards();
void ParticleHoleSym2();
void Explicit();

//Explicit lattice summation
complex<double> interaction(int m1, int m3, int m4, int No, vector<double> vpseu, string type="pa");
//this gives 'ED interaction matrix' for square torus.
complex<double> latticepp(LATTICE ll, int m1, int m2, int m3, int m4, string type);
//this gives 'sigma function lattice sum interaction matrix'.
void testlatticepp(double shift=0.);
vector<double> pairamplitude_ExplicitLatticeSum2(int invNu, double shift, vector<NQ> PP);
vector<double> pairamplitude_ExplicitLatticeSum3(int invNu, double shift, vector<NQ> PP);

#endif
