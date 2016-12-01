#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
using namespace std;
#include <iomanip>
#include "lattice.h"
#include "berry_tests.h"

bool IsOdd (int i) {
    return ((i%2)==1);
}
int main(){
//    double theta=0.5*M_PI, alpha=1.0;
    int num_core=2;
//    CFL_berry_phases_parallel("params", "ne8bp", num_core, "fullloop", theta, alpha);//params_name, output_name, num_core, kind.
//	ParticleHoleSym2();

//    double theta=0.333333333333*M_PI, alpha=1.0;

//    vector<double> length;
//    double steplength=0.1;
    
//    structurefactor("params_lau20_1_1");
//    structurefactor("params_sq69");
//    single_run();
//    structurefactor("params_sq10", 2);
    
//    structurefactor("params_sq49", num_core);
//    structurefactor("params_sq36", num_core);
//    structurefactor("params_sq25", num_core);
//    structurefactor("params_sq16", num_core);
//    structurefactor("params_sq9", num_core);
//    structurefactor("params_sq64", num_core);
    
//    structurefactor("params_sq6new", num_core);
//    structurefactor("params_sq6asp12test", num_core);
//    structurefactor("params_sq6asp15test", num_core);
//    structurefactor("params_sq6theta1test", num_core);
//    structurefactor("params_sq6asp2test", num_core);
    
//    structurefactor("params_sq4", num_core);
    
//    single_run("paramsne20test");
//    single_run("paramsne15test");
    
//    structurefactor("paramsne2test", num_core);
//    structurefactor("paramsne3test", num_core);
//    structurefactor("paramsne4test", num_core);
//    structurefactor("paramsne5test", num_core);
    
//    structurefactor("paramsne4theta1", num_core);
//    structurefactor("paramsne4theta2", num_core);
//    structurefactor("paramsne4theta3", num_core);
//    structurefactor("paramsne4theta4", num_core);
//    structurefactor("paramsne4theta5", num_core);
//    structurefactor("paramsne4theta6", num_core);
//    structurefactor("paramsne4theta7", num_core);
//    structurefactor("paramsne4theta8", num_core);
    
//    structurefactor("paramsne4squ", num_core);
//    structurefactor("paramsne4hex", num_core);
    
//    structurefactor("params_sq6asp95new", num_core);
//    structurefactor("params_sq6hexasp95new", num_core);
    
//    pairamplitude("paramsne4", true, num_core);
//    pairamplitude("paramsne6", true, num_core);//...
//    pairamplitude("paramsne10", true, num_core);
//    pairamplitude("paramsne20", true, num_core);
    
    pairamplitude("paramsne2", true, num_core);//...
    
}
void outformfactor(){
    int Ne=8, invNu=2, Nphi=Ne*invNu, seed=0;
    LATTICE lat(Ne, invNu, false, "CFL", seed, 0);
    cout<<lat.formfactor(0,0)<<endl;
    ofstream outfile("formfactor");
    for (int i=0; i<Nphi+1; i++) {
        outfile<<abs(lat.formfactor(0,i))<<endl;
    }
    outfile.close();
}
void test_Scott_finding8(){
    int Ne=8, invNu=2;
    int nMeas=5000, nWarmup=1000000, nSteps=20, seed=0;
    vector<vector<int>> ds_tmp(Ne+1, vector<int>(2));
    vector<vector<vector<int>>> ds(3);
    
    for (int i=-1; i<2; i++) {
        for (int j=-1; j<2; j++) {
            ds_tmp[(i+1)*3+j+1][0]=i; ds_tmp[(i+1)*3+j+1][1]=j;
        }
    }
    for (int i=0; i<3; i++) ds[i]=ds_tmp;
    ds[0].erase(remove(ds[0].begin(),ds[0].end(),vector<int>{1,1}),ds[0].end());
    ds[0].erase(remove(ds[0].begin(),ds[0].end(),vector<int>{1,-1}),ds[0].end());
    ds[0].push_back(vector<int>{2, 1});
    ds[1].erase(remove(ds[1].begin(),ds[1].end(),vector<int>{1,0}),ds[1].end());
    ds[1].erase(remove(ds[1].begin(),ds[1].end(),vector<int>{1,-1}),ds[1].end());
    ds[1].push_back(vector<int>{2,0});
    ds[2].erase(remove(ds[2].begin(),ds[2].end(),vector<int>{0,-1}),ds[2].end());
    
    //    for (int i=0; i<3; i++) {
    //        cout<<"i="<<i<<endl;
    //        for (int j=0; j<ds[i].size(); j++) {
    //            cout<<ds[i][j][0]<<" "<<ds[i][j][1]<<endl;
    //        }
    //        cout<<endl;
    //    }
    
    vector<LATTICE> psi(3);
    for (int i=0; i<3; i++) {
        psi[i]=LATTICE(Ne, invNu, false, "CFL", seed, 0); psi[i].set_ds(ds[i]);
    }
    
    Eigen::Matrix3cd overlap, overlap2;
    for (int m=0; m<3; m++) {
        psi[m].reset();
        psi[m].step(nWarmup);
        for (int nmea=0; nmea<nMeas; nmea++) {
            psi[m].step(nSteps);
            for (int n=0; n<=m; n++) {
                complex<double> tmp=psi[n].get_wf(psi[m].get_locs())/psi[m].get_wf(psi[m].get_locs());
                overlap(n,m)+=tmp;
                overlap2(n,m)+=norm(tmp);
            }
        }
    }
    overlap.array()/=(1.*nMeas); overlap2.array()/=(1.*nMeas); overlap.array()/=overlap2.array().sqrt();
    
    for (int m=0; m<3; m++) {
        for (int n=m; n<3; n++) {
            overlap(n,m)=conj(overlap(m,n));
        }
    }
    Eigen::ComplexEigenSolver<Eigen::Matrix3cd> es(overlap);
    cout<<"overlap=\n"<<overlap<<endl;
    cout<<"det="<<overlap.determinant()<<endl;
    cout<<endl<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<abs(es.eigenvalues()[2])<<endl;
}
void test_Scott_finding4(){
    int Ne=4, invNu=2;
    int nMeas=5000, nWarmup=1000000, nSteps=20, seed=0;
    vector<vector<int>> ds_tmp(4, vector<int>(2));
    for (int i=0; i<2; i++) for (int j=0; j<2; j++) {ds_tmp[i*2+j][0]=i; ds_tmp[i*2+j][1]=j;}
    
    vector<vector<vector<int>>> ds(2, ds_tmp);

    ds[1].erase(remove(ds[1].begin(),ds[1].end(),vector<int>{1,1}),ds[1].end());
    ds[1].erase(remove(ds[1].begin(),ds[1].end(),vector<int>{0,0}),ds[1].end());
    ds[1].push_back(vector<int>{2,1}); ds[1].push_back(vector<int>{-1,0});

//    for (int i=0; i<2; i++) {
//        cout<<"i="<<i<<endl;
//        for (int j=0; j<ds[i].size(); j++) {
//            cout<<ds[i][j][0]<<" "<<ds[i][j][1]<<endl;
//        }
//        cout<<endl;
//    }
    
    vector<LATTICE> psi(2);
    for (int i=0; i<2; i++) {
        psi[i]=LATTICE(Ne, invNu, false, "CFL", seed, 0); psi[i].set_ds(ds[i]);
    }

    Eigen::Matrix2cd overlap, overlap2;
    for (int m=0; m<2; m++) {
        psi[m].reset();
        psi[m].step(nWarmup);
        for (int nmea=0; nmea<nMeas; nmea++) {
            psi[m].step(nSteps);
            for (int n=0; n<2; n++) {
                complex<double> tmp=psi[n].get_wf(psi[m].get_locs())/psi[m].get_wf(psi[m].get_locs());
                overlap(n,m)+=tmp;
                overlap2(n,m)+=norm(tmp);
            }
        }
    }
    overlap.array()/=(1.*nMeas); overlap2.array()/=(1.*nMeas); overlap.array()/=overlap2.array().sqrt();
    
    complex<double> tmp=(overlap(0,1)+conj(overlap(1,0)))/2.;
    overlap(0,1)=tmp; overlap(1,0)=conj(tmp);
    overlap(0,0)=1.; overlap(1,1)=1.;
    
    cout<<"overlap=\n"<<overlap<<endl;
//    cout<<"det="<<overlap.determinant()<<endl;
}
void phase_variance(){
    ofstream outfile("M");
    int Ne=8,invNu=2,nWarmup=5000,nMeas=100,nSteps=20,nBins=1000,seed=0;
    outfile<<"Ne=8, invNu=2, nWarmup="<<nWarmup<<", nMeas="<<nMeas<<", nSteps="<<nSteps<<", nBins="<<nBins<<endl;
    bool testing=false;
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
    for (int nbin=0; nbin<nBins; nbin++) {
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
    for (int nbin=0; nbin<nBins; nbin++) {
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
            for (int i=0; i<(signed)locs.size(); i++) {
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
