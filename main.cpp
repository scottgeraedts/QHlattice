using namespace std;
#include "lattice.h"
#include "berry_tests.h"
#include <iomanip>
bool IsOdd (int i) {
    return ((i%2)==1);
}
int main(){
    int num_core=2;
//    CFL_berry_phases_parallel("params_ne8", "ne8bp", num_core, "loop1");//params_name, output_name, num_core, kind.
//    plot_CFL_coule_vsdbar(2, 21, 5000, 10);//grid, Ne, nMeas, nBins.
//    CFL_ne5_energy_var(5000, 100);//nMeas, nBins, num_core.
//    CFL_ne5_energy_var(5000, 100, 6);//nMeas, nBins.
//    CFL_ne5_energy_var(50, 10, 6);//nMeas, nBins.
    
    
//    
//    int Ne1=3, Ne2=4, invNu=2;
//    vector<vector<int>> z(Ne1+Ne2, vector<int>(2));
//    for (unsigned i=0; i<Ne1+Ne2; i++) {
//        z[i][0]=i; z[i][1]=i;
//    }
//    LATTICE cfl1(Ne1, invNu, false, "CFL", 0, 0), cfl2(Ne2, invNu, false, "CFL", 0, 0);
//    complex<double> ret=asymwf(cfl1, cfl2, z);
//    cout<<ret<<endl;
    
    ParticleHoleSym();
    
}
void phase_variance(){
    ofstream outfile("M");
    int Ne=8,invNu=2,nWarmup=5000,nMeas=100,nSteps=20,nBins=1000,seed=0;
    outfile<<"Ne=8, invNu=2, nWarmup="<<nWarmup<<", nMeas="<<nMeas<<", nSteps="<<nSteps<<", nBins="<<nBins<<endl;
    bool testing=false;
    int Nphi=Ne*invNu;
    //initialize MC object
    
    //this instance of LATTICE is only to set up the circular fermi surface of tempNe electrons
    LATTICE templl(Ne+1, invNu, testing, "CFL", seed, 0);
    vector<vector<int> > old_ds=templl.get_ds(), extra_ds;
    //old_dbar is the center of the circular fermi surface
    vector<double> old_dbar=templl.get_dbar_parameter();
    
    extra_ds.push_back(vector<int>{1,-1});
    extra_ds.push_back(vector<int>{1,1});
    
    vector<vector<int>> ds0=old_ds, ds1=old_ds;
    ds0.erase(remove(ds0.begin(), ds0.end(), extra_ds[0]), ds0.end());//missing {1,-1}.
    ds1.erase(remove(ds1.begin(), ds1.end(), extra_ds[1]), ds1.end());//missing {1,1}.
    
    vector<LATTICE> ll(invNu), pp(invNu);
    for (int i=0; i<invNu; i++) {
        ll[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
        ll[i].set_ds(ds0);
        pp[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
        pp[i].set_ds(ds1);
    }
    
    //many body K.
    for (int i=0; i<invNu; i++) {
        for (int j=0; j<2; j++) {
            if (ll[i].dsum[j]%invNu!=0 || pp[i].dsum[j]%invNu!=0) {
                cout<<"dsum mod invNu != 0, somewhere wrong!"<<endl;
                exit(0);
            }
        }
    }
    int dKx=ll[0].dsum[0]/invNu-pp[0].dsum[0]/invNu, dKy=ll[0].dsum[1]/invNu-pp[0].dsum[1]/invNu;
    //alphabar = K1 L1/Nphi + K2 L2/Nphi = d/invNu. (K1,K2) are many body momentums.
    //So for d = d1 L1/Ne + d2 L2/Ne => K1=sum d1, K2=sum d2.
    //dKx/y is divided by invNu, because in lattice.cpp dsum is defined on L/Nphi lattice.
    
    //Monte Carlo Part & Output.
    for (unsigned nbin=0; nbin<nBins; nbin++) {
        //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
        vector<Eigen::MatrixXcd> overlaps(4, Eigen::MatrixXcd::Zero(invNu, invNu));
        
        for (int i=0; i<invNu; i++) {
            ll[i].reset();
            pp[i].reset();
            ll[i].step(nWarmup);
        }
        
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<invNu; i++)
                ll[i].step(nSteps);
            
            for (int i=0; i<invNu; i++) {
                for (int j=0; j<invNu; j++) {
                    complex<double> temp;
                    temp=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[0](i,j)+=temp*ll[i].rhoq(dKx,dKy,ll[i].get_locs());
                    overlaps[1](i,j)+=norm(temp);
                    temp=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[2](i,j)+=temp;
                    overlaps[3](i,j)+=norm(temp);
                }
            }
        }
        
        for (int l=0; l<4; l++) overlaps[l]/=(1.*nMeas);
        overlaps[0]=overlaps[0].array()/overlaps[1].array().sqrt();
        overlaps[2]=overlaps[2].array()/overlaps[3].array().sqrt();
        hermitianize(overlaps[2]);
        
        //berry matrix.
        Eigen::MatrixXcd berrymatrix=overlaps[2].inverse() * overlaps[0];
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix);
        outfile<<nbin<<" "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[1])<<endl;
        outfile<<nbin<<" "<<berrymatrix(0,0).real()<<" "<<berrymatrix(0,0).imag()<<" "<<berrymatrix(1,1).real()<<" "<<berrymatrix(1,1).imag()<<" "<<berrymatrix(0,1).real()<<" "<<berrymatrix(0,1).imag()<<" "<<berrymatrix(1,0).real()<<" "<<berrymatrix(1,0).imag()<<endl;
        outfile<<endl;
    }
    outfile.close();
}
void onestep(int ne, string output_name){
    ofstream outfile(output_name.c_str());
    int Ne=ne,invNu=2,nWarmup=5000,nMeas=500,nSteps=20,nBins=1,seed=0,Nphi=Ne*invNu;
    outfile<<"Ne="<<Ne<<", invNu=2, nWarmup="<<nWarmup<<", nMeas="<<nMeas<<", nSteps="<<nSteps<<", nBins="<<nBins<<endl;
    bool testing=false;
    //initialize MC object
    
    vector<vector<vector<int>>> counter(Ne, vector<vector<int>>(Nphi, vector<int>(Nphi)));
    
    //set ds.
    vector<vector<vector<int>>> ds(2);
    if (Ne==2) {
        ds[0].push_back(vector<int>{0, 0}); ds[0].push_back(vector<int>{1, 0});
        ds[1].push_back(vector<int>{0, 0}); ds[1].push_back(vector<int>{0, 1});
    }
    if (Ne==3) {
        ds[0].push_back(vector<int>{ 1, 0});
        ds[0].push_back(vector<int>{-1, 0});
        ds[0].push_back(vector<int>{ 0, 1});
        ds[1].push_back(vector<int>{ 1, 0});
        ds[1].push_back(vector<int>{ 0, 1});
        ds[1].push_back(vector<int>{ 0, 0});
    }
    if (Ne==4) {
        ds[0].push_back(vector<int>{0, 0});
        ds[0].push_back(vector<int>{0, 1});
        ds[0].push_back(vector<int>{0, -1});
        ds[1].push_back(vector<int>{0, 0});
        ds[1].push_back(vector<int>{0, 1});
        ds[1].push_back(vector<int>{0, -1});
        
        ds[0].push_back(vector<int>{-1, 0});
        ds[1].push_back(vector<int>{-1, 1});
    }
    if (Ne==8) {
        LATTICE templl(9, 2, false, "CFL", 0, 0);
        vector<vector<int> > old_ds=templl.get_ds(), extra_ds;
        templl.print_ds();
        
        ds[0]=old_ds; ds[1]=old_ds;
        
        extra_ds.push_back(vector<int>{1,1});
        extra_ds.push_back(vector<int>{1,-1});
        extra_ds.push_back(vector<int>{-1,-1});
        extra_ds.push_back(vector<int>{-1,1});
        
        ds[0].erase(remove(ds[0].begin(),ds[0].end(),extra_ds[0]),ds[0].end());
        ds[1].erase(remove(ds[1].begin(),ds[1].end(),extra_ds[1]),ds[1].end());
        
        //        ds[0].erase(remove(ds[0].begin(),ds[0].end(),extra_ds[1]),ds[0].end());
        //        ds[1].erase(remove(ds[1].begin(),ds[1].end(),extra_ds[2]),ds[1].end());
        //
        //        ds[0].erase(remove(ds[0].begin(),ds[0].end(),extra_ds[2]),ds[0].end());
        //        ds[1].erase(remove(ds[1].begin(),ds[1].end(),extra_ds[3]),ds[1].end());
        //
        //        ds[0].erase(remove(ds[0].begin(),ds[0].end(),extra_ds[3]),ds[0].end());
        //        ds[1].erase(remove(ds[1].begin(),ds[1].end(),extra_ds[0]),ds[1].end());
    }
    if (Ne==20) {
        LATTICE templl(21, 2, false, "CFL", 0, 0);
        vector<vector<int> > old_ds=templl.get_ds();
        
        ds[0]=old_ds; ds[1]=old_ds;
        ds[0].erase(remove(ds[0].begin(),ds[0].end(),vector<int>{0,-1}),ds[0].end());
        ds[1].erase(remove(ds[1].begin(),ds[1].end(),vector<int>{0,-2}),ds[1].end());
    }
    
    //set cfl wfs.
    vector<LATTICE> ll(invNu), pp(invNu);
    for (int i=0; i<invNu; i++) {
        ll[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
        ll[i].set_ds(ds[1]);
        pp[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
        pp[i].set_ds(ds[0]);
    }
    //many body K.
    for (int i=0; i<invNu; i++) {
        for (int j=0; j<2; j++) {
            if (ll[i].dsum[j]%invNu!=0 || pp[i].dsum[j]%invNu!=0) {
                cout<<"dsum mod invNu != 0, somewhere wrong!"<<endl;
                exit(0);
            }
        }
    }
    int dKx=ll[0].dsum[0]/invNu-pp[0].dsum[0]/invNu, dKy=ll[0].dsum[1]/invNu-pp[0].dsum[1]/invNu;
    //alphabar = K1 L1/Nphi + K2 L2/Nphi = d/invNu. (K1,K2) are many body momentums.
    //So for d = d1 L1/Ne + d2 L2/Ne => K1=sum d1, K2=sum d2.
    //dKx/y is divided by invNu, because in lattice.cpp dsum is defined on L/Nphi lattice.
    
    //Monte Carlo Part & Output.
    for (unsigned nbin=0; nbin<nBins; nbin++) {
        //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
        vector<Eigen::MatrixXcd> overlaps(4, Eigen::MatrixXcd::Zero(invNu, invNu));
        vector<Eigen::MatrixXcd> overlaps_density(4, Eigen::MatrixXcd::Zero(invNu, invNu));
        
        for (int i=0; i<invNu; i++) {
            ll[i].reset();
            pp[i].reset();
            ll[i].step(nWarmup);
        }
        
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<invNu; i++)
                ll[i].step(nSteps);
            
            vector<vector<int>> locs=ll[0].get_locs();
            for (int i=0; i<locs.size(); i++) {
                counter[i][locs[i][0]][locs[i][1]]++;
            }
            
            for (int i=0; i<invNu; i++) {
                for (int j=0; j<invNu; j++) {
                    complex<double> temp;
                    //overlap.
                    temp=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[0](i,j)+=temp;//excluding density operator.
                    overlaps[1](i,j)+=norm(temp);
                    temp=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[2](i,j)+=temp;
                    overlaps[3](i,j)+=norm(temp);
                    //overlap_density.
                    temp=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps_density[0](i,j)+=temp*ll[i].rhoq(dKx,dKy,ll[i].get_locs());//including density operator.
                    overlaps_density[1](i,j)+=norm(temp);
                    temp=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps_density[2](i,j)+=temp;
                    overlaps_density[3](i,j)+=norm(temp);
                }
            }
        }
        
        for (int l=0; l<4; l++) overlaps[l]/=(1.*nMeas);
        overlaps[0]=overlaps[0].array()/overlaps[1].array().sqrt();// <ll|pp>
        overlaps[2]=overlaps[2].array()/overlaps[3].array().sqrt();// <ll|ll>
        hermitianize(overlaps[2]);
        
        for (int l=0; l<4; l++) overlaps_density[l]/=(1.*nMeas);
        overlaps_density[0]=overlaps_density[0].array()/overlaps_density[1].array().sqrt();// <ll|rhoq|pp>
        overlaps_density[2]=overlaps_density[2].array()/overlaps_density[3].array().sqrt();// <ll|ll>
        hermitianize(overlaps_density[2]);
        
        //        //berry matrix.
        //        Eigen::MatrixXcd berrymatrix=overlaps[2].inverse() * overlaps[0];
        //        Eigen::MatrixXcd berrymatrix_density=overlaps_density[2].inverse() * overlaps_density[0];
        ////        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix);
        ////        outfile<<nbin<<" "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[1])<<endl;
        //        outfile<<nbin<<" "<<berrymatrix(0,0).real()<<" "<<berrymatrix(0,0).imag()<<" "<<berrymatrix(1,1).real()<<" "<<berrymatrix(1,1).imag()<<" "<<berrymatrix(0,1).real()<<" "<<berrymatrix(0,1).imag()<<" "<<berrymatrix(1,0).real()<<" "<<berrymatrix(1,0).imag()<<endl;
        //        outfile<<nbin<<" "<<berrymatrix_density(0,0).real()<<" "<<berrymatrix_density(0,0).imag()<<" "<<berrymatrix_density(1,1).real()<<" "<<berrymatrix_density(1,1).imag()<<" "<<berrymatrix_density(0,1).real()<<" "<<berrymatrix_density(0,1).imag()<<" "<<berrymatrix_density(1,0).real()<<" "<<berrymatrix_density(1,0).imag()<<endl;
        //        outfile<<endl;
        
        //outputs: <ll|pp> matrix, <ll|rhoq|pp>, <ll|ll>.
        outfile<<nbin<<" ";
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                outfile<<overlaps[0](i,j).real()<<" "<<overlaps[0](i,j).imag()<<" ";
            }
        }
        outfile<<endl;
        outfile<<nbin<<" ";
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                outfile<<overlaps_density[0](i,j).real()<<" "<<overlaps_density[0](i,j).imag()<<" ";
            }
        }
        outfile<<endl;
        outfile<<nbin<<" ";
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                outfile<<overlaps[2](i,j).real()<<" "<<overlaps[2](i,j).imag()<<" ";
            }
        }
        outfile<<endl;
    }
    
    outfile.close();
}