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

//Some Utilities
void hermitianize(Eigen::MatrixXcd &x);
double phasemod(complex<double> in);

//Some Checkings, including ortho, auto-correlation
void check_orthogonality(string type);
void single_run(string filename="params");
void structurefactor(string);
void structurefactor(string, int);

//Energetics
void coul_energy(LATTICE& lattice, int nWarmup, int nMeas, int nSteps, int nBins, string filename);
void plot_CFL_coule_vsdbar(int grid, int Ne, int nMeas, int nBins);
void CFL_ne5_energy_var(int nMeas, int nBins, int num_core=1);
void findstate();

//Laughlin-Hole Berry Phase
void laughlinberryphase(string params_name, vector<double> length, double steplength, int change_nMeas, int change_Ne, int num_core, double theta=0.5*M_PI, double alpha=1.0);
void laughlin_bp_single_state(int gs, vector<double> length, double steplength, vector<data> &datas);
void two_holes(string str, int nmeasurement, data& test);
void test_error(int ne, double loop, double steplength, int nMea, int ncore, string test, int num_core);

//CFL Berry Phase
void CFL_berry_phases_parallel(string params_name, string output_name, int num_core, string kind, double theta=0.5*M_PI, double alpha=1.0);

//Particle-Hole Symmetry
void ParticleHoleSym();
void ParticleHoleSym2();
void testIQHwf();
complex<double> landauwf(int Nphi, int n, vector<double> latticeshift, vector<int> z, double theta=0.5*M_PI, double alpha=1.0);

//...
void GetCoefficient(vector<int> input);

#endif
