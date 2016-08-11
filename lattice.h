#ifndef LATTICE_H
#define LATTICE_H

#include <stdio.h>
#include "utils.h"
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include "MersenneTwister.h"
#include <numeric>
#include <deque>
//#include "weierstrass.h"

using namespace std;

//this number is multiplied by every term in the determinants to ensure that they don't overflow
//const double in_determinant_rescaling=0.2;
const double in_determinant_rescaling=1.;

extern"C"{
	void z_function_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z);
	void z_function_with_modular_transform_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z, int *sl2z);
	complex<double> lattice_z_(int *NPhi, int *x, int *y, complex<double> *l1, complex<double> *l2, int * use_table);
	double v_coulomb_(int *NPhi, int *m, int *n, complex<double> *l1, complex <double> *l2); //coulomb2_m.f90
	double new_v_coulomb_(int *NPhi, int *m, int *n, complex<double> *l1, complex <double> *l2); //coulomb2_m.f90
	void setup_z_function_table_();
	void set_l_(int *NPhi, complex<double> *l1, complex <double> *l2);
	void setup_laughlin_state_(int *Ne, int *invNu, int *sl2z, int *k); //wl_tools.f90
	void get_laughlin_cm_(int *x, complex<double> *wf);
}

class LATTICE{
    friend class berry_phase;
public:
	LATTICE();
//    LATTICE(int, int, bool, string, int seed, int gs=0);
    LATTICE(int Ne_t, int invNu_t, bool testing_t=false, string type_t="CFL", double theta=0.5*M_PI, double alpha=1.0, int seed=0, int gs_t=0);
    LATTICE(int Ne_t, int invNu_t, bool testing_t, string type_t, int seed, int gs_t) : LATTICE(Ne_t, invNu_t, testing_t, type_t, 0.5*M_PI, 1.0, seed, gs_t){} ;
    // dbar = (dbar_parameter[0]/NPhi*L1, dbar_parameter[1]/NPhi*L2);
	~LATTICE();

	int Ne;
	double running_weight;//running_weight is a global variable. need reset in every run.
    int tries,accepts;
//    Eigen::MatrixXcd det_M;
    
    vector<double> hole;
    bool fermions,holes_set;

	//stepping functions
    vector<double> dbar_parameter;
	bool testing; // output 'running_weight' and 'get_weight()', useful in debug.
	void step(int);// step(int Nsteps); Nsetps = total MC steps. tries:steps, accepts:updated steps.
	double get_weight(const vector< vector<int> > &zs);  
	complex<double> get_wf(const vector< vector<int> > &zs);
//    void make_CFL_COM(complex<double>& value);
    void make_CFL_det(Eigen::MatrixXcd& newMatrix, vector<int> newloc, int electron, complex<double>& value);

	//utility functions
	complex<double> modded_lattice_z(int x, int y);
	void print_ds();//.

	//initialization related functions
	void make_fermi_surface(double* center_frac, int N);
	void reset();
	vector <vector<int> > get_ds();
	void change_dbar_parameter(double dbarx, double dbary);
	vector<double> get_dbar_parameter();
	void set_ds(vector< vector<int> > ds);//sets a custom set of composite fermion dipole moments
	void set_hole(vector<double> temphole);
	vector<double> get_hole();

	//measurement related functions
	double coulomb_energy();
	double threebody();
	void update_structure_factors();
	void print_structure_factors(int nMeas);
	vector <vector<int> > get_locs();
	complex<double> formfactor(int qx, int qy);
	complex<double> rhoq(int qx ,int qy, const vector< vector<int> > &zs);
	
private:
	void sum_locs(int []);
	void setup_coulomb();
//	void setup_weierstrass();
	int simple_update();// returns '1' if updated, '0' if not updated.
	vector<int> random_move(const vector<int> &oldsite);
	int p(int); int m(int);
	void cold_start();
	void hotter_start();
	void det_helper(const vector<int> &z1, const vector<int> &z2, const vector<int> &d, vector<int> &z);
    double det_helper(int z1, int z2, int d, double dbar_parameter);
//    complex<double> jies_weierstrass(double x, double y);
    void check_sanity();

    complex<double> L1,L2;
	int NPhi, invNu;
	string type;
    int gs;
    double theta, alpha;
	vector <vector <double> > coulomb_table,sq2;
	vector <vector <complex<double> > > sq;
	vector <vector <vector <vector< complex<double> > > > > sq3;
	vector <vector<int> > sx,sx2;
	vector <vector <complex<double> > > shifted_ztable;
	Eigen::MatrixXcd oldMatrix;
	complex<double> oldDeterminant;
	Eigen::FullPivLU<Eigen::MatrixXcd> detSolver;	
	MTRand ran;
	vector< vector<int> > locs;//an integer defined on an NPhi lattice
	vector< vector<double> > ws, ws0;
	int one,zero;
	vector< complex<double> > omega;
    vector<vector<int> > ds;//an integer defined on an Ne lattice
    vector<int> dsum;//an integer defined on an NPhi lattice
};

#endif
