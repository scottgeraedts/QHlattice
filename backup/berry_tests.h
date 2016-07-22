#ifndef BERRY_TESTS_H
#define BERRY_TESTS_H

#include "utils.h"
#include "lattice.h"

struct data{
    int num;
    double position[2];
    double amp[3];
    double ang[3];
    double energy;
    double ang_trace;
    double det;
    double dfromnorm;
};
void two_holes(string str, int nmeasurement, data& test);
complex<double> chop(complex<double> input);
Eigen::MatrixXcd chop(Eigen::MatrixXcd mat);
void hermitianize(Eigen::MatrixXcd &x);
void testeigen();
void test_largesize();
void single_run();
void plot_CFL_coule_vsdbar(int grid);
void coul_energy_CFL_dbar(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins, double* dbar_parameter);
void coul_energy_laughlin(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins);
double phasemod(complex<double> in);
int supermod(int k, int n);


void laughlin_bp_single_state(int gs, vector<double> length, double steplength, vector<data> &datas);
void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne, int num_core);
void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne);
void test_error(int ne, double loop, int nMea, int ncore, string filename, string test);

#endif
