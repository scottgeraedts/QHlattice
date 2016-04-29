#ifndef LATTICE_H
#define LATTICE_H

#include "utils.h"
#include <Eigen/LU>
#include "MersenneTwister.h"
#include <numeric>
#include <deque>
#include "weierstrass.h"

using namespace std;

extern"C"{
	void z_function_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z);
	void z_function_with_modular_transform_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z, int *sl2z);
	complex<double> lattice_z_(int *NPhi, int *x, int *y, complex<double> *l1, complex<double> *l2, int * use_table);
	double v_coulomb_(int *NPhi, int *m, int *n, complex<double> *l1, complex <double> *l2); //coulomb2_m.f90
	void setup_z_function_table_();
	void set_l_(int *NPhi, complex<double> *l1, complex <double> *l2);
	void setup_laughlin_state_(int *Ne, int *invNu, int *sl2z, int *k); //wl_tools.f90
	void get_laughlin_cm_(int *x, complex<double> *wf);
}

class LATTICE{
public:
	LATTICE(int, int, int); // LATTICE(int Nphi_t, int invNu_t, int seed);
    // dbar = (dbar_parameter[0]/NPhi*L1, dbar_parameter[1]/NPhi*L2);
	~LATTICE();

	void step(int);// step(int Nsteps); Nsetps = total MC steps. tries:steps, accepts:updated steps.
	double get_weight(const vector< vector<int> > &zs);  
	complex<double> get_wf(const vector< vector<int> > &zs);
	void make_fermi_surface(double* center_frac, int N);
	complex<double> modded_lattice_z(int x, int y);
	void print_ds();//.
	double coulomb_energy();
	double threebody();
	void update_structure_factors();
	void print_structure_factors(int nMeas);
	void reset();
	void change_dbar_parameter(double dbarx, double dbary);
	void set_ds(vector< vector<int> > ds);//sets a custom set of composite fermion dipole moments
	vector <vector<int> > get_locs();

	vector <vector<int> > get_ds();
	
	int Ne;
	bool testing; // output 'running_weight' and 'get_weight()', useful in debug.
	double running_weight;
    int tries,accepts;
    double dbar_parameter[2];
    
    void make_CFL_COM(complex<double>& value);
    void make_CFL_det(Eigen::MatrixXcd& newMatrix, vector<int> newloc, int electron, complex<double>& value);
    
private:
	void sum_locs(int []);
	void setup_coulomb();
	void setup_weierstrass();
	int simple_update();// returns '1' if updated, '0' if not updated.
	vector<int> random_move(const vector<int> &oldsite);
	int p(int); int m(int);
	void cold_start();
	void det_helper(const vector<int> &z1, const vector<int> &z2, const vector<int> &d, vector<int> &z);
    double det_helper(int z1, int z2, int d, double dbar_parameter);
    complex<double> jies_weierstrass(double x, double y);
    weierstrass weiers;

	int NPhi, invNu;
	complex<double> L1,L2;
	vector <vector <double> > coulomb_table,sq2;
	vector <vector <complex<double> > > sq;
	vector <vector <vector <vector< complex<double> > > > > sq3;
	vector < vector<int> > sx,sx2;
	vector <vector <int> > ds;
	vector <vector <complex<double> > > shifted_ztable;
	vector<int> dsum;
	Eigen::MatrixXcd oldMatrix;
	complex<double> oldDeterminant;
	Eigen::FullPivLU<Eigen::MatrixXcd> detSolver;	
	MTRand ran;
	bool fermions;
	vector< vector<int> > locs;
	vector< vector<double> > ws;
	string type;
	int one,zero;
	int * sl2z;
	vector< complex<double> > omega;
};
#endif
