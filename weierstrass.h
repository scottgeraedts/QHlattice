//
//  weierstrass.h
//  weierstrass
//
//  Created by Jie Wang on 3/31/16.
//  Copyright © 2016 Jie Wang. All rights reserved.
//
#ifndef weierstrass_h
#define weierstrass_h

#include <complex>
#include <vector>

using namespace std;

extern"C"{
    void z_function_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z);//weierstrass_gaussian @ x*l1+y*l2.
    
    complex<double> lattice_z_(int *NPhi, int *x, int *y, complex<double> *l1, complex<double> *l2, int * use_table);//weierstrass_gaussian @ x*l1/No+y*l2/No.
    
    void z_function_with_modular_transform_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z, int *sl2z);
    
    //z_function_ and lattice_z_ are the same when z=(m*l1+n*l2)/N0, m,n are integers.
    
}


class weierstrass{
private:
    complex<double> tau;
    complex<double> ii;
    vector<complex<double> > zeros;
    double pre;
    
public:
    weierstrass(double, double, int, vector<complex<double> >);
    weierstrass(complex<double>, complex<double>, vector<complex<double> >);
    weierstrass();
    void run_weierstrass();
    complex<double> weta1,weta2,Gbar;
    int A;
    complex<double> omega1, omega2;
    int Nphi;
    double alpha, theta;
    double pi;
    
    void geteta_1();
    void geteta_2();
    complex<double> wsigma(complex<double> z);
    complex<double> wsigma_gaussian(complex<double> z);
};

#endif /* weierstrass_h */
