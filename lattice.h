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
struct NQ{
    int N;
    double Q;
    vector<double> screen;
};
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
    LATTICE(int Ne_t, int invNu_t, bool testing_t=false, string type_t="CFL", int seed=0, int gs_t=0, double theta=0.5*M_PI, double alpha=1.0, double shift_x=0., double shift_y=0., bool trace_t=false, bool correlatedsampling=false);
	~LATTICE();

	int Ne, NPhi;
	double running_weight;//running_weight is a global variable. need reset in every run.
    double c_running_weight;//the running_weight for correlated sampling.
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
    complex<double> get_laughlinwf(vector<vector<double> > z);
    void make_CFL_det(Eigen::MatrixXcd& newMatrix, vector<int> newloc, int electron, complex<double>& value, const vector< vector<int> > &zs);
    double updateweight(const int &electron, const vector<int> &newloc);

	//utility functions
	void print_ds();
    void print_ws();
    
    //boundary conditions
    vector<double> get_shift();
    vector<double> get_bc();//if CFL related,bc=sum_of_w-sum_of_d=wsum-dsum/NPhi.

	//initialization related functions
	void make_fermi_surface(double* center_frac, int N);
	void reset();
	void reset(const vector< vector<int> > &zs);
	vector <vector<int> > get_ds();
	void change_dbar_parameter(double dbarx, double dbary);
	vector<double> get_dbar_parameter();
    int trace;
    
    //initialization or reset related functions
	void set_ds(vector< vector<int> > ds);
    void set_zeros(vector<double> zeros0);//for "FilledLL" state only, set zeros[0][0][0], zeros[0][0][1];
	void set_hole(vector<double> temphole);
	vector<double> get_hole();

	//measurement related functions
	double threebody();
    
	vector <vector<int> > get_locs();
    void change_locs(const vector<vector<int>> &locs_t);
    void update_matdet(bool=true);
    void init_matdet();
    
	complex<double> formfactor(int qx, int qy);
	complex<double> rhoq(int qx ,int qy, const vector< vector<int> > &zs);
    vector<int> dsum;//an integer defined on an NPhi lattice

	//product of 2 CFL wavefunctions
    complex<double> doubled_CFL(const vector<vector<int>> &zs);
    
    //stuff for lattice_wrapper
    complex<double> update_weight(const vector<vector<int>> &zs, int electron, vector<int> newz);
    void update();
	static vector< vector<int> > hot_start(int NPhi_t, int Ne_t, MTRand &ran);
	static vector<int> random_move(const vector<int> &oldsite, int NPhi_t, MTRand &ran_t);
    
    //stuff for correlatedsampling
    void setup_nonsamplestates(vector<vector<vector<int>>>, vector<double>);
    vector<double> get_runweis();//return runing_weight of CFL excited states.
    vector<double> get_ratio();
    
    //structure factor
    void update_structure_factors();
    void print_structure_factors(int nMeas, string filename="");
    
    //Landau table
    void setup_landautable();
    vector<vector<vector<complex<double>>>> landautable;
    
    //high LL Coulomb and Pair-amplitude things.
    //these two functions set up table used by coulomb/pair-amplitude. type = "ce" or "pa".
    void setup_table(NQ, vector<vector<vector<double>>>&, vector<vector<vector<double>>>&, string typee="ce");
    void setup_tables(vector<NQ>, string typee="ce");
    
    void setup_coulomb0();
    vector<vector<double>> coulomb_table;
    vector<vector<vector<vector<double>>>> coulomb_tableHLL;//[LLN][Screen][x][y]
    vector<vector<vector<vector<double>>>> PA_table;
    vector<double> coulomb_energy(int, string typee="ce");
    void shortrange(int ind, vector<double>&, vector<double>&, string);
    double get_correlated_weight();
	
private:
    bool samplestate;
	void init(int seed);
    double get_in_det_rescaling(int Ne, int invNu);
    double shiftx, shifty;//should = b.c.
    //TODO:const initialization.
	void sum_locs(int []);
    
	int simple_update();// returns '1' if updated, '0' if not updated.
    
	complex<double> modded_lattice_z(int x, int y);
	complex<double> modded_lattice_z_nodbar(int x, int y);
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
	vector <vector <double> > sq2, sq2_mqy;
	vector <vector <complex<double> > > sq, sq_mqy;//'minus qy', qy<=0.
	vector <vector <vector <vector< complex<double> > > > > sq3;
	vector <vector<int> > sx,sx2;
//    void shift_ws();//boundary conditions.

	vector <vector <complex<double> > > shifted_ztable;
	int omega_dbar_range;
	vector< vector<complex<double> > > omega_dbar;
    
    vector<vector<complex<double>>> Ftable;
    vector<vector<vector<vector<double>>>> coulomb_qtable;
    vector<vector<vector<vector<double>>>> PA_qtable;
    vector<double> CE_cutoff;
    vector<double> PA_cutoff;

	Eigen::MatrixXcd oldMatrix, newMatrix;
	complex<double> oldDeterminant, newDeterminant;
	Eigen::FullPivLU<Eigen::MatrixXcd> detSolver;
	MTRand ran;
	vector<vector<int>> locs;//integers defined on an NPhi lattice
	vector<vector<double>> ws, ws0;
    vector<double> wsum;
    vector<vector<vector<double>>> zeros;//a parameter for 'FilledLL' state.
	int one,zero;
	vector< complex<double> > omega;
    
    vector<vector<int> > ds;//an integer defined on an Ne lattice
    
    //stuff for correlated sampling (private).
    bool correlatedsampling;
    vector<LATTICE> nonsamplestates;
    vector<double> ratio;
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
