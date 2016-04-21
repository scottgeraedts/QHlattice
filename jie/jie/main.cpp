//
//  main.cpp
//  jie
//
//  Created by Jie Wang on 4/14/16.
//  Copyright © 2016 Jie Wang. All rights reserved.
//

#include <iostream>
#include<fstream>
#include <vector>
#include <complex>
using namespace std;

int main(int argc, const char * argv[]) {
    complex<double> l1,l2; l1=1; l2 = complex<double>(0,1);
    int Ne=1000;
    //vector<vector<complex<double> > > ds;
    vector<complex<double> > ds;
    void make_fermi_surface(complex<double>, complex<double>, int, double*, vector<complex<double> >&);
    
    //complex<double> metric[2*2]; metric[0]=1.; metric[1]=0.; metric[2]=0.; metric[3]=1.;
    
    double center_frac[2]; center_frac[0]=0.5; center_frac[1]=0.5;
    make_fermi_surface(l1, l2, Ne, center_frac, ds);
    
    
    ofstream myfile;
    myfile.open ("example.txt");
    for (int i=0; i<ds.size(); i++) {
        myfile << real(ds[i])<<"   "<<imag(ds[i]);
        myfile<<endl;
    }
    
    myfile.close();
    
    
    
    
    return 0;
}

/*
make_fermi_surface(l1,l2,Ne,center_frac,ds);
 l1, l2 are primitive lattice. Ne is #e.
 center_frac = (x0,y0). x0*l1+y0*l2 is the center position of fermi surface.
 ds contains d s .
*/
void make_fermi_surface(complex<double> l1, complex<double>l2, int Ne, double* center_frac, vector<complex<double> >& ds){
    vector<vector<int> > d_list;
    double x0,y0; x0=center_frac[0]; y0=center_frac[1];
    //initial sub-lattice: L_{mn}/Ne where d s lives on.
    vector<int> d; d=vector<int>(2,0);
    d_list.push_back(d);
    for (int i=-Ne; i<Ne; i++) {
        for (int j=-Ne; j<Ne; j++) {
            if (i==0 && j==0) continue;
            vector<int> d;
            d.push_back(i); d.push_back(j);
            d_list.push_back(d);
            d.clear();
        }
    }
    for (int i=0; i<d_list.size(); i++) {
        d_list[i][0]+=(int)(x0*Ne); d_list[i][1]+=(int)(y0*Ne);
    }
    
    //fill electrons one by one, fill those close to center first.
    int x,y; double min,min_;
    vector< vector<int> >::iterator it;
    
    for (int k=0; k<Ne; k++) {
        it=d_list.begin();
        
        x=d_list[0][0]; y=d_list[0][1]; min = norm(l1*((double)x/(1.*Ne)-x0)+l2*((double)y/(1.*Ne)-y0));
        for (int i=1; i<d_list.size(); i++) {
            x=d_list[i][0]; y=d_list[i][1]; min_ = norm(l1*((double)x/(1.*Ne)-x0)+l2*((double)y/(1.*Ne)-y0));
            if (min_<min) {
                it=d_list.begin()+i;
                min=min_;
            }
        }

        vector<complex<double> > tmp;
        ds.push_back(1.*(*it)[0]*l1/(1.*Ne)+1.*(*it)[1]*l2/(1.*Ne));
        it=d_list.erase(it);
    }
    
}

