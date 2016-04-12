#ifndef LATTICE_H
#define LATTICE_H

#include "utils.h"
#include "MersenneTwister.h"
#include <numeric>
#include <deque>

using namespace std;

extern"C"{
	void z_function_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z);
	complex<double> lattice_z_(int *NPhi, int *x, int *y, complex<double> *l1, complex<double> *l2, int * use_table);
	double v_coulomb_(int *NPhi, int *m, int *n, complex<double> *l1, complex <double> *l2);
	void setup_z_function_table_();
	void set_l_(int *NPhi, complex<double> *l1, complex <double> *l2);
	void setup_laughlin_state_(int *Ne, int *invNu, int *sl2z, int *k);
	void get_laughlin_cm_(int *x, complex<double> *wf);
}

class LATTICE{
public:
	LATTICE(int, int, int);
	~LATTICE();
	void step(int);
	double get_weight();
	void sum_locs(int []);
	int simple_update();
	vector<int> random_move(const vector<int> &oldsite);
	double coulomb_energy();
	void setup_coulomb();
	void setup_weierstrass();
	void update_structure_factors();
	void print_structure_factors(int nMeas);
	
	int Ne;
	bool testing;
	double running_weight;
	int tries,accepts;
	
private:
	int p(int); int m(int);

	int NPhi, invNu;
	complex<double> L1,L2;
	vector <vector <double> > coulomb_table,sq2;
	vector <vector <complex<double> > > sq;
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
