#include "lattice.h"
//#include "berry_phase.h"

complex<double> chop(complex<double> input){
    double rea=real(input), ima=imag(input);
    if (rea<pow(10,-5)) {
        rea=0.;
    }
    if (ima<pow(10,-5)) {
        ima=0.;
    }
    return rea+complex<double>(0,1)*ima;
}

Eigen::MatrixXcd chop(Eigen::MatrixXcd mat){
    int m=mat.rows(), n=mat.cols();
    Eigen::MatrixXcd ret = Eigen::MatrixXcd(m,n);
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            ret(i,j)=chop(mat(i,j));
        }
    }
    return ret;
}
void hermitianize(Eigen::MatrixXcd &x){
    x(0,1)=0.5*(x(0,1)+conj(x(1,0)));
    x(1,0)=conj(x(0,1));
    x(0,2)=0.5*(x(0,2)+conj(x(2,0)));
    x(2,0)=conj(x(0,2));
    x(1,2)=0.5*(x(1,2)+conj(x(2,1)));
    x(2,1)=conj(x(1,2));
    
}
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

int main(){
    void error_loop(); void error_Ne(); void error_nMear1(); void error_nMear2(); void error_nMear3(); void error_nMear4(); void error_normality(); void error_steplength();

    void error_ne_sl_para(); void error_steplength_ll_ne10_para();
//    error_nMear1();
//    error_nMear2();
//    error_nMear3();
//    error_nMear4();
//    error_normality();
//    error_steplength();
//    error_loop();
//    error_Ne();
    
    void error_nMear_slne4();
//    error_nMear_slne4();
    
    
    void error_ne_sl();
    void error_ne_ll();
    void error_loop();
//    error_ne_sl();
//    error_ne_ll();
//    error_loop();
    
    void error_ne_ll_para(); void error_ne_ll_para2();
//    error_ne_ll_para();
//    error_ne_ll_para2();
    
    void error_steplength_ll_ne10_para2();
//    error_steplength_ll_ne10_para2();
    
    
    void error_nMear_para_slne4();
    void error_nMear_para_llne4();
//    error_nMear_para_slne4();
//    error_nMear_para_llne4();
    
    void error_normality_ll_para(); void error_normality_sl_para();
//    error_normality_sl_para();
//    error_normality_ll_para();
    
    void error_normality_ll_(); void error_normality_sl_();
//    error_normality_sl_();
//    error_normality_ll_();
    
    void error_steplength_sl_ne10_para1(); void error_steplength_sl_ne10_para2();
//    error_steplength_sl_ne10_para1();
//    error_steplength_sl_ne10_para2();
    
    void error_loop_ne6(); void error_loop_ne8(); void error_loop_ne10();
//    error_loop_ne8();
//    error_loop_ne6();
    error_loop_ne10();
    
    void error_steplength_sl_ne10_para3(); void error_steplength_sl_ne10_3();
//    error_steplength_sl_ne10_para3();
//    error_steplength_sl_ne10_3();
}

void test_parallel(int ncore){
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne);
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne, int num_core);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    double steplength=0.01;
    
    double btime, etime;
    
    btime=omp_get_wtime();
    laughlinberryphase(length, steplength, datas, 0, 0);
    etime=omp_get_wtime();
    cout<<"non-parallel time cost = "<<etime-btime<<endl;
    
    cout<<"\n-----parallel programming, ncore="<<ncore<<"-----"<<endl;
    btime=omp_get_wtime();
    laughlinberryphase_parallel(length, steplength, datas, 0, 0, ncore);
    etime=omp_get_wtime();
    cout<<"parallel time cost = "<<etime-btime<<endl;
}

void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne){
    int Ne,Ne_t,invNu,nWarmup,nMeas,nMeas_t,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile("params");
    infile>>Ne_t>>invNu;
    infile>>nWarmup>>nMeas_t>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    if (change_nMeas==0) nMeas=nMeas_t;
    else nMeas=change_nMeas;
    if (change_Ne==0) Ne=Ne_t;
    else Ne=change_Ne;
    
    
    vector<vector<double> > holes; vector<int> Grid(2);
    for (int i=0; i<2; i++) {Grid[i]=(int)(length[i]/steplength);}
    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=steplength*i; a[1]=0.;                  holes.push_back(a);}
    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=steplength*i; a[0]=length[0];           holes.push_back(a);}
    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=length[0]-steplength*i; a[1]=length[1]; holes.push_back(a);}
    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=length[1]-steplength*i; a[0]=0.;        holes.push_back(a);}
    int nds=holes.size();
    vector<vector<double> > holes2(nds, vector<double>(2,0));
    int supermod(int k, int n);
    for(int i=0;i<nds;i++) holes2[supermod(i-1,nds)]=holes[i];//(holes[b],holes2[b]) = (holes[b],holes[b+1]).
    
    vector<LATTICE> ll(3), pp(3);
    for (int i=0; i<3; i++) {ll[i]=LATTICE(Ne, invNu, testing, type, seed, i); pp[i]=LATTICE(Ne, invNu, testing, type, seed, i);}
    
    vector<vector<Eigen::MatrixXcd > > overlaps;
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
    for (int b=0; b<nds; b++) {
        vector<Eigen::MatrixXcd> aa;
        Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(3,3);
        for (int i=0; i<4; i++) aa.push_back(a);
        overlaps.push_back(aa);
    }
    
    for(int b=0; b<nds; b++) {
        for (int i=0; i<3; i++) {
            ll[i].set_hole(holes[b]);
            pp[i].set_hole(holes2[b]);
            ll[i].reset(); ll[i].step(nWarmup);
            pp[i].reset();
        }
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<3; i++) ll[i].step(nSteps);
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    vector<complex<double> > temp(2);
                    temp[0]=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[b][0](i,j)+=temp[0];
                    overlaps[b][1](i,j)+=norm(temp[0]);
                    temp[1]=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[b][2](i,j)+=temp[1];
                    overlaps[b][3](i,j)+=norm(temp[1]);
                }
            }
        }
        for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {overlaps[b][0](i,j)/=sqrt(abs(overlaps[b][1](i,j))); overlaps[b][2](i,j)/=sqrt(abs(overlaps[b][3](i,j)));}
        
        hermitianize(overlaps[b][2]);
    }
    
    vector<Eigen::MatrixXcd> berrymatrix_step(nds);
    for (int b=0; b<nds; b++) berrymatrix_step[b] = overlaps[b][2].inverse() * overlaps[b][0];
    
    Eigen::Matrix3cd berrymatrix_integral = Eigen::Matrix3cd::Identity(3,3);
    vector<double> phases(3,0.);
    datas.clear();
    for (int b=0; b<nds; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step[b]);
        data tmp;
        berrymatrix_integral *= berrymatrix_step[b];
        for (int i=0; i<3; i++) {
            phases[i]+=arg(es.eigenvalues()[i]);
            tmp.num = b; tmp.amp[i] = abs(es.eigenvalues()[i]); tmp.ang[i] = arg(es.eigenvalues()[i]);
        }
        
        //insert dfromnorm.
        double normeigenvalue=0., normmatrix=0.;
        for (int i=0; i<3; i++) {
            normeigenvalue+=sqrt(norm(es.eigenvalues()[i]));
        }
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                normmatrix+=sqrt(norm(berrymatrix_step[b](i,j)));
            }
        }
        tmp.dfromnorm=normmatrix-normeigenvalue;
        
        datas.push_back(tmp);
    }
    
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_integral);
    datas[0].ang_trace = arg(berrymatrix_integral.trace());
    datas[0].det = arg(berrymatrix_integral.determinant());
    cout<<"\n\n Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<endl;
    cout<<"phase sum = "<<phases[0]<<" "<<phases[1]<<" "<<phases[2]<<"\nphase average = "<<(phases[0]+phases[1]+phases[2])/3<<endl;
    cout<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
    cout<<"amp(berrymatrix_integral.eigenvalue) = "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[0])<<endl;
    cout<<"arg(berrymatrix_integral.eigenvalue) = "<<arg(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[2])<<endl;
    cout<<"sum arg(berrymatrix_integral.eigenvalue) = "<< ( arg(es.eigenvalues()[0])+arg(es.eigenvalues()[1])+arg(es.eigenvalues()[2]) )/3 <<endl;
    cout<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
    cout<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
    cout<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl;
}

void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne, int num_core){
    int Ne,Ne_t,invNu,nWarmup,nMeas,nMeas_t,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile("params");
    infile>>Ne_t>>invNu;
    infile>>nWarmup>>nMeas_t>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    if (change_nMeas==0) nMeas=nMeas_t;
    else nMeas=change_nMeas;
    if (change_Ne==0) Ne=Ne_t;
    else Ne=change_Ne;
    
    vector<vector<double> > holes; vector<int> Grid(2);
    for (int i=0; i<2; i++) {Grid[i]=(int)(length[i]/steplength);}
    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=steplength*i; a[1]=0.;                  holes.push_back(a);}
    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=steplength*i; a[0]=length[0];           holes.push_back(a);}
    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=length[0]-steplength*i; a[1]=length[1]; holes.push_back(a);}
    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=length[1]-steplength*i; a[0]=0.;        holes.push_back(a);}
    int nds=holes.size();
    vector<vector<double> > holes2(nds, vector<double>(2,0));
    int supermod(int k, int n);
    for(int i=0;i<nds;i++) holes2[supermod(i-1,nds)]=holes[i];//(holes[b],holes2[b]) = (holes[b],holes[b+1]).
    
    //     vector<LATTICE> ll(3), pp(3);
    //     for (int i=0; i<3; i++) {ll[i]=LATTICE(Ne, invNu, testing, type, seed, i); pp[i]=LATTICE(Ne, invNu, testing, type, seed, i);}
    omp_set_num_threads(num_core);
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(3)), pp(num_core, vector<LATTICE>(3));//do this to avoid wrong memory access since openmp share memory.
    for (int k=0; k<num_core; k++) for (int i=0; i<3; i++) {ll[k][i]=LATTICE(Ne, invNu, testing, type, k, i); pp[k][i]=LATTICE(Ne, invNu, testing, type, k, i);}
    //seed = k.
    
    vector<vector<Eigen::MatrixXcd > > overlaps;
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
    for (int b=0; b<nds; b++) {
        vector<Eigen::MatrixXcd> aa;
        Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(3,3);
        for (int i=0; i<4; i++) aa.push_back(a);
        overlaps.push_back(aa);
    }
    
    
    for(int b=0; b<nds; b++) {
        //        for (int i=0; i<3; i++) {
        //            ll[i].set_hole(holes[b]);
        //            pp[i].set_hole(holes2[b]);
        //            ll[i].reset(); ll[i].step(nWarmup);
        //            pp[i].reset();
        //        }
        for (int k=0; k<num_core; k++) {
            for (int i=0; i<3; i++) {
                ll[k][i].set_hole(holes[b]);
                pp[k][i].set_hole(holes2[b]);
                ll[k][i].reset(); ll[k][i].step(nWarmup);
                pp[k][i].reset();
            }
        }
        
        vector<vector<Eigen::MatrixXcd> > overlaps_tmp(num_core, vector<Eigen::MatrixXcd>(4));
        for (int i=0; i<num_core; i++) for (int j=0; j<4; j++) overlaps_tmp[i][j].setZero(3,3);
        
        
        //         ofstream bbout("test_openmp");
#pragma omp parallel for
        for (int k=0; k<nMeas; k++) {
            //             bbout<<omp_get_thread_num()<<endl;
            int coren = omp_get_thread_num();
            for (int i=0; i<3; i++) ll[coren][i].step(nSteps);
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    vector<complex<double> > temp(2);
                    temp[0]=pp[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                    overlaps_tmp[coren][0](i,j)+=temp[0];
                    overlaps_tmp[coren][1](i,j)+=norm(temp[0]);
                    temp[1]=ll[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                    overlaps_tmp[coren][2](i,j)+=temp[1];
                    overlaps_tmp[coren][3](i,j)+=norm(temp[1]);
                }
            }
        }
        
        //         for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
        for (int l=0; l<4; l++) for (int i=0; i<num_core; i++) overlaps[b][l]+=overlaps_tmp[i][l]/(1.*nMeas);
        
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {overlaps[b][0](i,j)/=sqrt(abs(overlaps[b][1](i,j))); overlaps[b][2](i,j)/=sqrt(abs(overlaps[b][3](i,j)));}
        hermitianize(overlaps[b][2]);
    }
    
    
    vector<Eigen::MatrixXcd> berrymatrix_step(nds);
    for (int b=0; b<nds; b++) berrymatrix_step[b] = overlaps[b][2].inverse() * overlaps[b][0];
    
    Eigen::Matrix3cd berrymatrix_integral = Eigen::Matrix3cd::Identity(3,3);
    vector<double> phases(3,0.);
    datas.clear();
    for (int b=0; b<nds; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step[b]);
        data tmp;
        berrymatrix_integral *= berrymatrix_step[b];
        for (int i=0; i<3; i++) {
            phases[i]+=arg(es.eigenvalues()[i]);
            tmp.num = b; tmp.amp[i] = abs(es.eigenvalues()[i]); tmp.ang[i] = arg(es.eigenvalues()[i]);
        }
        
        //insert dfromnorm.
        double normeigenvalue=0., normmatrix=0.;
        for (int i=0; i<3; i++) {
            normeigenvalue+=sqrt(norm(es.eigenvalues()[i]));
        }
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                normmatrix+=sqrt(norm(berrymatrix_step[b](i,j)));
            }
        }
        tmp.dfromnorm=normmatrix-normeigenvalue;
        
        datas.push_back(tmp);
    }
    
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_integral);
    datas[0].ang_trace = arg(berrymatrix_integral.trace());
    datas[0].det = arg(berrymatrix_integral.determinant());
    cout<<"\n\n Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<endl;
    cout<<"phase sum = "<<phases[0]<<" "<<phases[1]<<" "<<phases[2]<<"\nphase average = "<<(phases[0]+phases[1]+phases[2])/3<<endl;
    cout<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
    cout<<"amp(berrymatrix_integral.eigenvalue) = "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[0])<<endl;
    cout<<"arg(berrymatrix_integral.eigenvalue) = "<<arg(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[2])<<endl;
    cout<<"sum arg(berrymatrix_integral.eigenvalue) = "<< ( arg(es.eigenvalues()[0])+arg(es.eigenvalues()[1])+arg(es.eigenvalues()[2]) )/3 <<endl;
    cout<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
    cout<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
    cout<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl;
    
}


double phasemod(complex<double> in){
    double out=arg(in);
    if(out<0) return out+2*M_PI;
    else if (out>2*M_PI) return out-2*M_PI;
    else return out;
    return out;
}

void error_nMear_slne4(){
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    double steplength=0.01;
    ofstream bout("test_error_nMeas_slne4_new");
    
    for (int nMeas=2000; nMeas<10000; nMeas+=2000){
        laughlinberryphase(length, steplength, datas, nMeas, 0);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<nMeas<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}

void error_nMear_llne4(){
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.5; length[1]=0.5;
    double steplength=0.01;
    ofstream bout("test_error_nMeas_llne4");
    
    for (int nMeas=20000; nMeas<100000; nMeas+=20000){
        laughlinberryphase(length, steplength, datas, nMeas, 0);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<nMeas<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}


void error_nMear_para_slne4_para(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    double steplength=0.01;
    ofstream bout("test_error_nMeas_para_slne4");
    
    for (int nMeas=10000; nMeas<200000; nMeas+=20000){
        laughlinberryphase_parallel(length, steplength, datas, nMeas, 0, 10);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<nMeas<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_nMear_para_slne4(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    double steplength=0.01;
    ofstream bout("test_error_nMeas_slne4_para");
    
    for (int nMeas=10000; nMeas<200000; nMeas+=20000){
        laughlinberryphase_parallel(length, steplength, datas, nMeas, 0, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<nMeas<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_nMear_para_llne4(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.5; length[1]=0.5;
    double steplength=0.05;
    ofstream bout("test_error_nMeas_llne4_para");
    
    for (int nMeas=10000; nMeas<200000; nMeas+=20000){
        laughlinberryphase_parallel(length, steplength, datas, nMeas, 0, 10);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<nMeas<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_ne_sl_para(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    double steplength=0.01;
    ofstream bout("test_error_ne_sl_para");
    
    for (int ne=2; ne<22; ne=ne+2) {
        laughlinberryphase_parallel(length, steplength, datas, 0, ne, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<ne<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_ne_sl(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    double steplength=0.01;
    ofstream bout("test_error_ne_sl");
    
    for (int ne=2; ne<22; ne=ne+2) {
        laughlinberryphase(length, steplength, datas, 0, ne);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<ne<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_ne_ll(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.5; length[1]=0.5;
    double steplength=0.05;
    ofstream bout("test_error_ne_sl");
    
    for (int ne=2; ne<22; ne=ne+2) {
        laughlinberryphase(length, steplength, datas, 0, ne);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<ne<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_ne_ll_para(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.5; length[1]=0.5;
    double steplength=0.01;
    ofstream bout("test_error_ne_ll_para");
    
    for (int ne=12; ne<22; ne=ne+2) {
        laughlinberryphase_parallel(length, steplength, datas, 0, ne, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<ne<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_ne_ll_para2(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.5; length[1]=0.5;
    double steplength=0.01;
    ofstream bout("test_error_ne_ll_para2");
    
    for (int ne=14; ne<22; ne=ne+2) {
        laughlinberryphase_parallel(length, steplength, datas, 0, ne, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<ne<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_steplength_ll_ne10_para(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.5; length[1]=0.5;
    
    ofstream bout("test_error_steplength_ll_ne10_para2");
    for (double steplength=0.17; steplength<0.25; steplength+=0.02) {
        laughlinberryphase_parallel(length, steplength, datas, 0, 0, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}

void error_steplength_ll_ne10_para2(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.5; length[1]=0.5;
    
    ofstream bout("test_error_steplength_ll_ne10_para_new2");
    for (double steplength=0.003; steplength<0.01; steplength+=0.002) {
        laughlinberryphase_parallel(length, steplength, datas, 0, 0, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}

void error_steplength_sl_ne10_para1(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    
    ofstream bout("test_error_steplength_sl_ne10_para1_new");
    for (double steplength=0.0001; steplength<0.001; steplength+=0.0002) {
        laughlinberryphase_parallel(length, steplength, datas, 0, 0, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_steplength_sl_ne10_para2(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    
    ofstream bout("test_error_steplength_sl_ne10_para2_new");
    for (double steplength=0.001; steplength<0.05; steplength+=0.002) {
        laughlinberryphase_parallel(length, steplength, datas, 0, 0, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}

void error_steplength_sl_ne10_para3(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    
    ofstream bout("test_error_steplength_sl_ne10_para3");
    for (double steplength=0.0023; steplength<0.05; steplength+=0.002) {
        laughlinberryphase_parallel(length, steplength, datas, 0, 0, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}
void error_steplength_sl_ne10_3(){
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int, int);
    vector<double> length(2); vector<data> datas;
    length[0]=0.05; length[1]=0.05;
    
    ofstream bout("test_error_steplength_sl_ne10_3");
    for (double steplength=0.0023; steplength<0.05; steplength+=0.002) {
        laughlinberryphase(length, steplength, datas, 0, 0);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}

void error_loop(){
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int, int);
    vector<double> length(2); vector<data> datas;
    double steplength=0.01;
    ofstream bout("test_error_loop_ne4");
    for (double x=0.05; x<0.8; x+=0.05) {
        length[0]=x; length[1]=x;
        laughlinberryphase(length, steplength, datas, 0, 0);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<x<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}

void error_loop_ne6(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    double steplength=0.01;
    ofstream bout("test_error_loop_ne6_para");
    for (double x=0.05; x<0.8; x+=0.05) {
        length[0]=x; length[1]=x;
        laughlinberryphase_parallel(length, steplength, datas, 0, 6, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<x<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}

void error_loop_ne8(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    double steplength=0.01;
    ofstream bout("test_error_loop_ne8_para");
    for (double x=0.05; x<0.8; x+=0.05) {
        length[0]=x; length[1]=x;
        laughlinberryphase_parallel(length, steplength, datas, 0, 8, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<x<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}

void error_loop_ne10(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    double steplength=0.01;
    ofstream bout("test_error_loop_ne10_para");
    for (double x=0.05; x<0.8; x+=0.05) {
        length[0]=x; length[1]=x;
        laughlinberryphase_parallel(length, steplength, datas, 0, 10, 5);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<x<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
}

void error_normality_ll_para(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    //    double steplength=0.01;
    length[0]=0.5; length[1]=0.5;
    
    ofstream bout("test_error_normality_ll_para");
    for (double steplength=0.01; steplength<0.25; steplength+=0.01) {
        laughlinberryphase_parallel(length, steplength, datas, 0, 0, 5);
        bout<<steplength<<endl;
        for (int i=0; i<datas.size(); i++) {
            bout<<datas[i].num<<" "<<datas[i].dfromnorm<<endl;
        }
        bout<<endl;
    }
}

void error_normality_sl_para(){
    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
    vector<double> length(2); vector<data> datas;
    //    double steplength=0.01;
    length[0]=0.05; length[1]=0.05;
    
    ofstream bout("test_error_normality_sl_para");
    for (double steplength=0.005; steplength<0.05; steplength+=0.005) {
        laughlinberryphase_parallel(length, steplength, datas, 0, 0, 5);
        bout<<steplength<<endl;
        for (int i=0; i<datas.size(); i++) {
            bout<<datas[i].num<<" "<<datas[i].dfromnorm<<endl;
        }
        bout<<endl;
    }
}

void error_normality_ll_(){
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int, int);
    vector<double> length(2); vector<data> datas;
    //    double steplength=0.01;
    length[0]=0.5; length[1]=0.5;
    
    ofstream bout("test_error_normality_ll");
    for (double steplength=0.01; steplength<0.25; steplength+=0.01) {
        laughlinberryphase(length, steplength, datas, 0, 0);
        bout<<steplength<<endl;
        for (int i=0; i<datas.size(); i++) {
            bout<<datas[i].num<<" "<<datas[i].dfromnorm<<endl;
        }
        bout<<endl;
    }
}

void error_normality_sl(){
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int, int);
    vector<double> length(2); vector<data> datas;
    //    double steplength=0.01;
    length[0]=0.05; length[1]=0.05;
    
    ofstream bout("test_error_normality_sl");
    for (double steplength=0.005; steplength<0.05; steplength+=0.005) {
        laughlinberryphase(length, steplength, datas, 0, 0);
        bout<<steplength<<endl;
        for (int i=0; i<datas.size(); i++) {
            bout<<datas[i].num<<" "<<datas[i].dfromnorm<<endl;
        }
        bout<<endl;
    }
}