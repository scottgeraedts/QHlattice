#ifndef LATTICE_H
#define LATTICE_H

#include <stdio.h>
#include "utils.h"
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include "MersenneTwister.h"
#include <numeric>
#include <deque>
#include <unordered_set>
#include <cmath>
#include <iomanip>

using namespace std;

//this number is multiplied by every term in the determinants to ensure that they don't overflow
//const double in_determinant_rescaling=0.2;

double get_in_det_rescaling(int Ne, int invNu);

extern"C"{
	void z_function_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z);
	void z_function_with_modular_transform_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z, int *sl2z);
	complex<double> lattice_z_(int *NPhi, int *x, int *y, complex<double> *l1, complex<double> *l2, int * use_table);
//	double v_coulomb_(int *NPhi, int *m, int *n, complex<double> *l1, complex <double> *l2); //coulomb2_m.f90
	double new_v_coulomb_(int *NPhi, int *m, int *n, complex<double> *l1, complex <double> *l2); //coulomb2_m.f90
	void setup_z_function_table_();
	void set_l_(int *NPhi, complex<double> *l1, complex <double> *l2);
	void setup_laughlin_state_(int *Ne, int *invNu, int *sl2z, int *k); //wl_tools.f90
	void get_laughlin_cm_(int *x, complex<double> *wf);
	void jacobi_theta_(int *n, complex<double> *z, complex<double> *tau, complex<double> *theta, int *sum);
}

class LATTICE_PARAMS{
public:
	int Ne, invNu, seed, gs;
	double theta, alpha, rescale;
	complex<double> w_delta, dbar_delta;
	bool testing, trace;
	string type;
	LATTICE_PARAMS( int Ne_t){
		Ne=Ne_t;
		invNu=2;
		seed=0;
		gs=0;
		theta=M_PI*0.5;
		alpha=1.0;
		testing=false;
		trace=false;
		type="CFL";
		w_delta=0;
		dbar_delta=0;
		rescale=1;
	}
};
	
class LATTICE{
public:
	LATTICE();
	LATTICE(LATTICE_PARAMS params);
//    LATTICE(int Ne_t, int invNu_t, bool testing_t=false, string type_t="CFL", int seed=0, int gs_t=0, double theta=0.5*M_PI, double alpha=1.0);
    LATTICE(int Ne_t, int invNu_t, bool testing_t=false, string type_t="CFL", int seed=0, int gs_t=0, double theta=0.5*M_PI, double alpha=1.0, bool trace_t=false);
	~LATTICE();

	int Ne, NPhi;
	double running_weight;//running_weight is a global variable. need reset in every run.
    int tries,accepts;
    complex<double> getL(int dir);
    
    vector<double> hole;
    bool fermions,holes_set;
    
    double in_determinant_rescaling;
	//stepping functions
    vector<double> dbar_parameter;
	bool testing;
	void step(int);// step(int Nsteps); Nsetps = total MC steps. tries:steps, accepts:updated steps.
	double get_weight(const vector< vector<int> > &zs);  
	complex<double> get_wf(const vector< vector<int> > &zs);
    void make_CFL_det(Eigen::MatrixXcd& newMatrix, vector<int> newloc, int electron, complex<double>& value, const vector< vector<int> > &zs);

	//utility functions
	complex<double> modded_lattice_z(int x, int y);
	void print_ds();
    void print_ws();
    void print_landautable();
    void print_laguerreltable();
    vector<vector<vector<vector<double>>>> laguerretable;
    vector<vector<vector<complex<double>>>> landautable;
    
    void setup_laguerre_con();
    void setup_laguerre_lat();
    void setup_landautable();

	//initialization related functions
	void make_fermi_surface(double* center_frac, int N);
	void reset();
	void reset(const vector< vector<int> > &zs);
	vector <vector<int> > get_ds();
	void change_dbar_parameter(double dbarx, double dbary);
	vector<double> get_dbar_parameter();
    
    //initialization or reset related functions
	void set_ds(vector< vector<int> > ds);
    void set_ws(vector<vector<double>> ws);
    void set_zeros(vector<double> zeros0);//for "FilledLL" state only, set zeros[0][0][0], zeros[0][0][1];
	void set_hole(vector<double> temphole);
    //*********************
    //To Avoid Bugs, 'set_ws' must be followed by 'set_ds', 'change_dbar_parameter' must following 'set_ds'.
    
	vector<double> get_hole();

	//measurement related functions
	double coulomb_energy();
	double threebody();
	void update_structure_factors();
	void print_structure_factors(int nMeas, string filename="");
//    double pairamplitude(int n, double alpha);
    double pairamplitude(int n, int a);
    
	vector <vector<int> > get_locs();
	complex<double> formfactor(int qx, int qy);
//    double formfactor(double qx, double qy);
	complex<double> rhoq(int qx ,int qy, const vector< vector<int> > &zs);
    vector<int> dsum;//an integer defined on an NPhi lattice
    
    //Filled Landau Level Wavefunction (not IQH wf)
    complex<double> FilledLL(vector<vector<int>> zs);
	//another version of filled Landau level wavefunction, this one is normalized
    complex<double> FilledLL2(vector<vector<int>> zs);
	//product of 2 CFL wavefunctions
    complex<double> doubled_CFL(const vector<vector<int>> &zs);
    
    //stuff for lattice_wrapper
    complex<double> update_weight(const vector<vector<int>> &zs, int electron, vector<int> newz);
    void update();
	static vector< vector<int> > hot_start(int NPhi_t, int Ne_t, MTRand &ran);
	static vector<int> random_move(const vector<int> &oldsite, int NPhi_t, MTRand &ran_t);
    
    int trace;
	
private:
	void init(int seed);
    double get_in_det_rescaling(int Ne, int invNu);
	void sum_locs(int []);
	void setup_coulomb();
    
	int simple_update();// returns '1' if updated, '0' if not updated.
	int p(int); int m(int);
	void cold_start();
	void det_helper(const vector<int> &z1, const vector<int> &z2, const vector<int> &d, vector<int> &z);
    double det_helper(int z1, int z2, int d, double dbar_parameter);
    void check_sanity();

    complex<double> L1,L2;
	int invNu;
	complex<double> w_delta, dbar_delta; //an offset to the COM zeros, and to the zeros in the determinant
	string type;
    int gs;
    double theta, alpha;
    vector <double> pair_amp;
	vector <vector <double> > coulomb_table, sq2, sq2_mqy;
	vector <vector <complex<double> > > sq, sq_mqy;//'minus qy', qy<=0.
	vector <vector <vector <vector< complex<double> > > > > sq3;
    vector <vector <double> > SMAq;
	vector <vector<int> > sx,sx2;
	vector <vector <complex<double> > > shifted_ztable;
	Eigen::MatrixXcd oldMatrix, newMatrix;
	complex<double> oldDeterminant, newDeterminant;
	Eigen::FullPivLU<Eigen::MatrixXcd> detSolver;
//    Eigen::FullPivLU<Eigen::MatrixXcd> detSolver_FLL;
	MTRand ran;
	vector<vector<int>> locs;//an integer defined on an NPhi lattice
	vector<vector<double>> ws, ws0;
    vector<vector<vector<double>>> zeros;//a parameter for 'FilledLL' state.
	int one,zero;
	vector< complex<double> > omega;
//    int sqgrid=1;//the q in structure factor can be finer than Nphi*Nphi lattice.
//    vector< complex<double> > omegasq;
    
    vector<vector<int> > ds;//an integer defined on an Ne lattice
    
};

class wf_info{
public:
	wf_info();
	wf_info(bool conj, bool denom, int start, int end, int sign);
	LATTICE wf;
	bool conj,denom;
	vector< vector<int> > make_zs(const vector< vector<int> > &inzs);
	int start, end;
	int sign;
};

class LATTICE_WRAPPER{
public:
	LATTICE_WRAPPER(int Ne, vector<vector<wf_info>> &wfs_t, int seed, bool testing);
	int step(int);
	int step_fromwf(int);
	complex<double> get_wf(), get_wf(const vector<vector<int> > &tempzs);
	complex<double> get_wf(int i, int j), get_wf(int i, int j, const vector<vector<int>> &tempzs);
	void reset();
	void acceptance_rate();
	vector< vector<int> > get_zs();
	vector< vector<int> > hot_start();
	vector<vector<wf_info>> wfs;
private:
	
	bool testing;
	int Ne; //assume Ne=NPhi
	vector< vector<int> > zs;
	vector<int> random_move(vector<int>);
	int tries,accepts;
	double running_weight;
	MTRand ran;
	vector<vector<complex<double>>> oldweight;
};

#endif
