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
struct data{
    int num;
    double position[2];
    double amp[3];
    double ang[3];
};
int main(){
//    berry_phase bp(20);
//    bp.two_full_braiding();

//    berry_phase bp(16);
//    bp.two_full_braiding();
    
//	void single_run();
//	single_run();

//    void two_holes_scott();
//    two_holes_scott();
    string str;
    void two_holes(string str, int nmeasurement, data&);//str="test" or "".
    
//    ofstream bout("test_may20");
//    for (int i=500; i<10000; i=i+50) {
//        data test;
//        two_holes("test", i, test);
//        bout<<i<<" "<<test.amp[0]<<" "<<test.amp[1]<<" "<<test.amp[2]<<" "<<test.ang[0]<<" "<<test.ang[1]<<" "<<test.ang[2]<<endl;
//    }
    
//    data test;
//    two_holes("", 0, test);
    
    
    void laughlin_bp_single_state(int gs, vector<double> length, double steplength, vector<data> &datas);//gs = 0, 1, 2, labeling ground state.
    vector<double> length(2); double steplength = 0.01;
//    vector<vector<data> > datas; for (int i=0; i<3; i++) {vector<data> tmp; datas.push_back(tmp);}
    vector<data> datas;
    vector<vector<data> > datass(3);
    
    length[0]=0.5; length[1]=0.8;
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas);
    laughlinberryphase(length, steplength, datas);
    ofstream bout("output");
    for (int i=0; i<datas.size(); i++) {
        for (int gs=0; gs<3; gs++) {
            bout<<datas[i].num<<" "<<datas[i].position[0]<<" "<<datas[i].position[1]<<" "<<gs<<" "<<datas[i].amp[gs]<<" "<<datas[i].ang[gs]<<endl;
        }
    }
    
//    //use model w.f. rather than berry matrix.
//    int gs=0;
//    length[0]=0.5; length[1]=0.8;
//    for (int gs=0; gs<3; gs++) {
//        laughlin_bp_single_state(gs, length, steplength, datass[gs]);
//        double angle=0.;
//        for (int i=0; i<datass[gs].size(); i++) {
//            angle+=datass[gs][i].ang[gs];
//        }
//        cout<<" gs = "<<gs<<" total phase = "<<angle<<endl;
//    }
    
//    //model w.f.
//    int gs=0;
//    length[0]=0.5; length[1]=0.5;
//    ofstream bout("teststeplength");
//    for (double steplength=0.01; steplength<0.5; steplength=steplength+0.02) {
//        laughlin_bp_single_state(gs, length, steplength, datas);
//        double minamp = datas[0].amp[gs], maxamp = datas[0].amp[gs], ang=0.;
//        for (int i=0; i<datas.size(); i++) {
//            ang+=datas[i].ang[gs];
//            if (datas[i].amp[gs]>maxamp) {maxamp=datas[i].amp[gs];}
//            if (datas[i].amp[gs]<minamp) {minamp=datas[i].amp[gs];}
//        }
//        bout<<steplength<<" "<<maxamp<<" "<<minamp<<" "<<ang<<endl;
//    }
//    
//    ofstream bout("berry_laughlin_single_state");
//    for (int i=0; i<datas.size(); i++) {
//        bout<<datas[i].position[0]<<" "<<datas[i].position[1]<<" "<<datas[i].amp[gs]<<" "<<datas[i].ang[gs]<<endl;
//    }
    

    void testeigen();
//    testeigen();
    
    void test_largesize();
//    test_largesize();
}

void testeigen(){
    Eigen::Matrix2cd MM;
    complex<double> ii = complex<double> (0,1);
    MM<<1,2.+ii,3.-ii,4;
    cout<<"MM=\n"<<MM<<endl;
    cout<<"MM.adjoint =\n "<<MM.adjoint()<<endl;
//    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(MM);
//    cout<<"\n eigenvectors = \n"<<es.eigenvectors()<<endl;//eigenvectors are columns.
//    cout<<"eigenvalues = \n"<<es.eigenvalues()<<endl;
//    cout<<endl;
//    cout<<"eigenvectors . eigenvectors.adjoint = \n"<<es.eigenvectors()*es.eigenvectors().adjoint()<<endl;
//    cout<<"eigenvectors.adjoint . eigenvectors = \n"<<es.eigenvectors().adjoint()*es.eigenvectors()<<endl;
//    cout<<endl;
//    cout<<"eigenvectors.adjoit() * MM * eigenvectors = \n"<<es.eigenvectors().adjoint()*MM*es.eigenvectors()<<endl;
//    
//    Eigen::Matrix2f mat;
//    mat<<1,4,5,9;
//    cout<<"\nmat=\n"<<mat;
//    mat=mat*mat+mat;
////    cout<<"\nsqrt(mat)=\n"<<mat;
//    
//    Eigen::Vector2cd V; V(0)=1.; V(1)=ii;
//    cout<<"(1,1).squared norm = "<<V.squaredNorm()<<endl;
    
//    Eigen::Matrix2cd M2;
//    M2.col(0)=V;
//    cout<<"\nM2=\n"<<M2<<endl;
    
//    int supermod(int k, int n);
//    cout<<"supermod = "<<supermod(5,4)<<endl;
}
void test_largesize(){
    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    vector<double> hole0(2), hole1(2);
    hole0[0]=0.; hole0[1]=0.; hole1[0]=0.; hole1[1]=0.01;
    
    ofstream bout("test_large_size");
    vector<LATTICE> ll(invNu),ll2(invNu);
    Eigen::MatrixXcd overlaps = Eigen::MatrixXcd::Zero(invNu, invNu);
    
    for (int ne=15; ne<23; ne++) {
        for(int gs=0;gs<invNu;gs++){
            ll[gs]=LATTICE(ne,invNu,testing,type,seed,gs);
            ll[gs].set_hole(hole0);
            ll2[gs]=LATTICE(ne,invNu,testing,type,seed,gs);
        }
        complex<double> berry=0.;
        for(int gs1=0;gs1<invNu;gs1++){
            ll[gs1].reset();
            ll[gs1].step(nWarmup);
            complex<double> btemp;
            for(int i=0;i<nMeas;i++){
                ll[gs1].step(nSteps);
                //            energy+=ll[gs1].coulomb_energy();
                for(int gs2=0;gs2<invNu;gs2++){
                    ll2[gs2].set_hole(hole1);
                    ll2[gs2].reset();
                    berry=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
                    overlaps(gs1,gs2)+=berry;
                }
            }
        }
        
        overlaps/=(1.*nMeas);
        cout<<"overlapmatrix = \n"<<overlaps<<endl;
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps);
        bout<<ne<<" "<<abs(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[1])<<" "<<abs(es.eigenvalues()[2])<<" "<<arg(es.eigenvalues()[2])<<endl;
    }
}

void laughlin_bp_single_state(int gs, vector<double> length, double steplength, vector<data> &datas){
    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object

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
    
    LATTICE ll(Ne, invNu, testing, type, seed, gs), pp(Ne, invNu, testing, type, seed, gs);
    vector<vector<complex<double> > > overlaps;
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|psi(xb)|psi(xb+1)|^2>, overlaps[b][2](i,j)=overlaps[b][0](i,j)/sqrt{overlaps[b][1](i,j)}.
    for (int b=0; b<nds; b++) {
        vector<complex<double> > aa;
        aa.push_back(0.); aa.push_back(0.); aa.push_back(0.);
        overlaps.push_back(aa);
    }
    
    for(int b=0;b<nds;b++){
        ll.set_hole(holes[b]); pp.set_hole(holes2[b]);
        ll.reset(); ll.step(nWarmup);
        for(int k=0;k<nMeas;k++){
            ll.step(nSteps);
            complex<double> temp=pp.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
            overlaps[b][0]+=temp; overlaps[b][1]+=norm(temp);
        }
        for (int l=0; l<3; l++) {
            overlaps[b][l]/=(1.*nMeas);
        }
        overlaps[b][2] = overlaps[b][0]/sqrt(overlaps[b][1]);
    }
    
    datas.clear();
    for (int b=0; b<nds; b++) {
        data tmp;
        tmp.position[0]=holes[b][0]; tmp.position[1]=holes[b][1];
        tmp.amp[gs]=abs(overlaps[b][2]); tmp.ang[gs]=arg(overlaps[b][2]);
        datas.push_back(tmp);
    }
}

void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas){
    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
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
    
    vector<LATTICE> ll, pp;
    for (int i=0; i<3; i++) {
        LATTICE a(Ne, invNu, testing, type, seed, i);
        ll.push_back(a); pp.push_back(a);//ll=psi(xb), pp=psi(xb+1);
    }
    
    vector<vector<Eigen::MatrixXcd > > overlaps;
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2](i,j)=overlaps[b][0](i,j)/sqrt{overlaps[b][1](i,j)}, overlaps[b][3]=<psi(xb)|psi(xb)>.
    for (int b=0; b<nds; b++) {
        vector<Eigen::MatrixXcd> aa;
        Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(3,3);
        aa.push_back(a); aa.push_back(a); aa.push_back(a); aa.push_back(a);
        overlaps.push_back(aa);
    }
    
    for(int b=0; b<nds; b++) {
        for (int i=0; i<3; i++) {
            ll[i].set_hole(holes[b]);
            pp[i].set_hole(holes2[b]);
            ll[i].reset(); ll[i].step(nWarmup);
        }
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<3; i++) {
                ll[i].step(nSteps);
            }
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    vector<complex<double> > temp(2);
                    temp[0]=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[b][0](i,j)+=temp[0];
                    overlaps[b][1](i,j)+=norm(temp[0]);
                    if (i>=j) {
                        temp[1]=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                        overlaps[b][3](i,j)+=temp[1];
                    }
                }
            }
        }
        for (int i=0; i<3; i++) for (int j=i+1; j<3; j++) overlaps[b][3](i,j) = conj(overlaps[b][3](j,i));
        for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
//        cout<<"\n b = "<<b<<"overlaps[b][3] = \n"<<overlaps[b][3]<<endl;
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) overlaps[b][2](i,j)=overlaps[b][0](i,j)/sqrt(abs(overlaps[b][1](i,j)));
    }

    //comment this paragraph if want orthogonal.
    vector<Eigen::Matrix3cd> alphas;
    for (int b=0; b<nds; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b][3]);
        Eigen::Matrix3cd alpha;
        for (int i=0; i<3; i++) {
            Eigen::VectorXcd V = es.eigenvectors().col(i);
            complex<double> temp = es.eigenvalues()[i]*V.squaredNorm();
            V/=sqrt(temp);
            alpha.col(i)=V;
        }
        alphas.push_back(alpha);
    }
    
//    ofstream output("berrymatrixstep");
//    output<<"Ne = "<<Ne<<" length[0] = "<<length[0]<<" length[1] = "<<length[1]<<endl;
    Eigen::Matrix3cd berrymatrix_integral = Eigen::Matrix3cd::Identity(3,3);
    datas.clear();
    for (int b=0; b<nds; b++) {
        int supermod(int k, int n);
        Eigen::MatrixXcd berrymatrix_step = alphas[b].adjoint() * overlaps[b][2] * alphas[supermod(b+1,nds)];//for orthogonal.
//        Eigen::MatrixXcd berrymatrix_step = overlaps[b][2];//for non-orthogonal.
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step);
        
        data tmp; tmp.position[0] = holes[b][0]; tmp.position[1] = holes[b][1]; tmp.num=b;
        for (int gs=0; gs<3; gs++) {
            tmp.amp[gs] = abs(es.eigenvalues()[gs]); tmp.ang[gs] = arg(es.eigenvalues()[gs]);
        }
        datas.push_back(tmp);
        berrymatrix_integral *= berrymatrix_step;
//        output<<"b = "<<b<<"berrymatrix_step = \n"<<berrymatrix_step<<endl;
    }
    cout<<"\ntrace of total berry matrix = "<<berrymatrix_integral.trace()<<" ang(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
//    cout<<"\nberrymatrix_integral = \n"<<berrymatrix_integral<<endl;
    
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es_integral(berrymatrix_integral);
//    cout<<abs(es_integral.eigenvalues()[0])<<" "<<arg(es_integral.eigenvalues()[0])<<" "<<abs(es_integral.eigenvalues()[1])<<" "<<arg(es_integral.eigenvalues()[1])<<" "<<abs(es_integral.eigenvalues()[2])<<" "<<arg(es_integral.eigenvalues()[2])<<endl;
    double totalberryphase=0.; for (int i=0; i<3; i++) totalberryphase+=arg(es_integral.eigenvalues()[i])/3.;
//    cout<<"total berry phase = "<<totalberryphase<<endl;
    
    
    Eigen::MatrixXcd test;
    for (int i=0; i<10; i++) {
//        test = alphas[i].adjoint() * overlaps[i][3] * alphas[i];
//        cout<<" \n i = "<<i<<" test = \n"<<chop(test)<<endl;//expect to print out identity matrix.
//        cout<<"\n i = "<<i<<" alpha.adjoint * alpha = \n"<<alphas[i].adjoint()*alphas[i]<<endl;
    }
}

void two_holes(string str, int nmeasurement, data& test){
    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    vector<vector<double> > holes;
    double x=0.;
    vector<double> a(2);
    while (x<0.05) {
        a[0]=x; a[1]=0.;
        holes.push_back(a);
        x+=0.005;
    }
    int nds=holes.size();
    
    if (str=="test") {
        nMeas=nmeasurement;
        nds=1;
    }
    
    ofstream bout("twoholelaughlinnew");
    vector<LATTICE> ll(invNu),ll2(invNu);//ll is psi(x), ll2 is psi(x).
    vector<vector<Eigen::MatrixXcd> > overlaps;
    //overlaps[b][0]=<psi(0)|psi(xb)>, overlaps[b][1]=<|<psi(0)|psi(xb)>|^2>, overlaps[b][2](i,j)=overlaps[b][0](i,j)/sqrt{overlaps[b][1](i,j)}.
    for (int b=0; b<nds; b++) {
        vector<Eigen::MatrixXcd> aa;
        Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(3,3);
        aa.push_back(a); aa.push_back(a); aa.push_back(a);
        overlaps.push_back(aa);
    }
    
    for(int gs=0;gs<invNu;gs++){
        ll[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
        ll[gs].set_hole(holes[0]);
        ll2[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
    }
    
    for(int gs1=0;gs1<invNu;gs1++){
        ll[gs1].reset();
        ll[gs1].step(nWarmup);
        for(int i=0;i<nMeas;i++){
            ll[gs1].step(nSteps);
            for(int gs2=0;gs2<invNu;gs2++){
                for(int b=0;b<nds;b++){
                    ll2[gs2].set_hole(holes[b]);
                    ll2[gs2].reset();
                    complex<double> temp=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
                    overlaps[b][0](gs1,gs2)+=temp;
                    overlaps[b][1](gs1,gs2)+=norm(temp);
                }
            }
        }
        cout<<energy[gs1]/(1.*nMeas)<<endl;
    }
//    cout<<"\noverlaps[10][0]=\n"<<overlaps[10][0]<<endl;
//    cout<<"\noverlaps[10][1]=\n"<<overlaps[10][1]<<endl;
    
    /*
    for (int b=0; b<nds; b++) {
        for (int gs1=0;gs1<invNu;gs1++) {
            ll[gs1].reset();
            ll[gs1].step(nWarmup);
            for (int gs2=0; gs2<invNu; gs2++) {
                ll2[gs2].set_hole(holes[b]);
                ll2[gs2].reset();
                for (int i=0;i<nMeas;i++) {
                    ll[gs1].step(nSteps);
                    complex<double> temp=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
                    overlaps[b][0](gs1,gs2)+=temp;
                    overlaps[b][1](gs1,gs2)+=norm(temp);
                }
            }
        }
    }
    cout<<"\noverlaps[10][0]=\n"<<overlaps[10][0]<<endl;
    cout<<"\noverlaps[10][1]=\n"<<overlaps[10][1]<<endl;
     */
    
    Eigen::Matrix3cd berryloop = Eigen::Matrix3cd::Identity(3,3);
    for (int b=0; b<nds; b++) {
        for (int l=0; l<3; l++) overlaps[b][l]/=(1.*nMeas);
        for (int gs1=0; gs1<invNu; gs1++) {
            for (int gs2=0; gs2<invNu; gs2++) {
                double denorm = abs(overlaps[b][1](gs1,gs2));
                overlaps[b][2](gs1,gs2) = overlaps[b][0](gs1,gs2)/sqrt(denorm);
            }
        }
        berryloop*=overlaps[b][2];
    }

//    cout<<"\nberryloop=\n"<<berryloop<<endl;
//    cout<<"determinant=\n"<<arg(berryloop.determinant())<<endl;
    
    for (int b=0; b<nds; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b][2]);
        bout<<holes[b][0]<<" "<<holes[b][1]<<" "<<abs(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[1])<<" "<<abs(es.eigenvalues()[2])<<" "<<arg(es.eigenvalues()[2])<<endl;
//        for (int i=0; i<3; i++) {
//            test.amp[i]=abs(es.eigenvalues()[i]);
//            test.ang[i]=arg(es.eigenvalues()[i]);
//        }
    }
    
}

void two_holes_scott(){
    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    vector<vector<double> > holes;
    double x=0.;
    vector<double> a(2);
    while (x<0.01) {
        a[0]=x; a[1]=0.;
        holes.push_back(a);
        x+=0.02;
    }
    int nds=holes.size();
    
    
    ofstream bout("twoholelaughlin");
    vector<LATTICE> ll(invNu),ll2(invNu);
    vector<Eigen::MatrixXcd> overlaps12(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
    vector<Eigen::MatrixXcd> overlaps12_2(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
    Eigen::MatrixXcd overlaps11=Eigen::MatrixXcd::Zero(invNu, invNu);
    Eigen::MatrixXcd overlaps11_2=Eigen::MatrixXcd::Zero(invNu, invNu);
	
    vector<Eigen::MatrixXcd> overlaps22(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
    vector<Eigen::MatrixXcd> overlaps22_2(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
    
	for(int gs=0;gs<invNu;gs++){
		ll[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
		ll[gs].set_hole(holes[0]);
		ll2[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
    }
    
    complex<double> berry;
	for(int gs1=0;gs1<invNu;gs1++){    
	    ll[gs1].reset();
	    ll[gs1].step(nWarmup);
        for(int i=0;i<nMeas;i++){
            ll[gs1].step(nSteps);
            //	        energy+=ll[gs1].coulomb_energy();
            for(int gs2=0;gs2<invNu;gs2++){
                berry=ll[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
                overlaps11(gs1,gs2)+=berry;
                overlaps11_2(gs1,gs2)+=norm(berry);
                for(int b=0;b<nds;b++){
					ll2[gs2].set_hole(holes[b]);
                    ll2[gs2].reset();
                    berry=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
                    overlaps12[b](gs1,gs2)+=berry;
                    overlaps12_2[b](gs1,gs2)+=norm(berry);
                }
            }
        }
	}
    cout<<"energies:"<<endl;
    vector<double> energy(invNu,0);
    for (int gs1=0; gs1<invNu; gs1++) {
        for (int b=0; b<nds; b++) {
            
            ll2[gs1].set_hole(holes[b]);
            ll2[gs1].reset();
            ll2[gs1].step(nWarmup);
            for (int i=0; i<nMeas; i++) {
                ll2[gs1].step(nSteps);
                energy[gs1]+=ll2[gs1].coulomb_energy();
                for (int gs2=0; gs2<invNu; gs2++) {
                	berry=ll2[gs2].get_wf(ll2[gs1].get_locs())/ll2[gs1].get_wf(ll2[gs1].get_locs());
                    overlaps22[b](gs1,gs2)+=berry;
                    overlaps22_2[b](gs1,gs2)+=norm(berry);
                }
                
            }
        }
        cout<<energy[gs1]/(1.*nMeas)<<endl;
    }
        
    overlaps11/=(1.*nMeas);
    overlaps11_2/=(1.*nMeas);
	overlaps11=overlaps11.cwiseQuotient(overlaps11_2.cwiseSqrt());
	hermitianize(overlaps11);

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es11(overlaps11);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es,es22;
    Eigen::MatrixXcd diaged_berry;
    for (int b=0; b<nds; b++) {
        overlaps12[b]/=(1.*nMeas); 
        overlaps22[b]/=(1.*nMeas);
        overlaps12_2[b]/=(1.*nMeas);
        overlaps22_2[b]/=(1.*nMeas);
//		cout<<overlaps12[b]<<endl<<endl;
//		cout<<overlaps12_2[b]<<endl<<endl;
		
        overlaps12[b]=overlaps12[b].cwiseQuotient(overlaps12_2[b].cwiseSqrt());
        overlaps22[b]=overlaps22[b].cwiseQuotient(overlaps22_2[b].cwiseSqrt());
		hermitianize(overlaps12[b]);
		hermitianize(overlaps22[b]);
		cout<<overlaps12[b]<<endl<<endl;
		cout<<overlaps11<<endl<<endl;
		cout<<overlaps12[b]<<endl<<endl;
//		cout<<"------------------------------------"<<endl;
		
		es22.compute(overlaps22[b]);
		
//        
//        Eigen::MatrixXcd esll_diag = Eigen::MatrixXcd::Zero(invNu, invNu);
//        Eigen::MatrixXcd esll2_diag = Eigen::MatrixXcd::Zero(invNu, invNu);
//        esll_diag=chop(esll.eigenvectors().adjoint() * overlapsll * esll.eigenvectors());
//        esll2_diag=chop(esll2.eigenvectors().adjoint() * overlapsll2[b] * esll2.eigenvectors());
//        
//        cout<<"\nb="<<b<<endl;
//        cout<<"\noverlapsll = \n"<<overlapsll<<endl;
//        cout<<"\noverlapsll2["<<b<<"] = \n"<<overlapsll2[b]<<endl;
//        cout<<"esll.eigenvector.adjoint() * esll * esll.eigenvector = \n"<<esll_diag<<endl;
//        cout<<"esll2.eigenvector.adjoint() * esll2 * esll2.eigenvector = \n"<<esll2_diag<<endl;
//        cout<<endl;
//        
		Eigen::MatrixXcd eigvecs=es22.eigenvectors();
		cout<<eigvecs<<endl;
        diaged_berry = es22.eigenvectors().adjoint() * overlaps12[b] *  es22.eigenvectors();
//        cout<<"\noverlaps[b]="<<overlaps[b]<<endl;
        cout<<"diaged_berry="<<endl<<diaged_berry<<endl;
		es.compute(overlaps12[b]);
        bout<<holes[b][0]<<" "<<holes[b][1]<<" "<<abs(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[1])<<" "<<abs(es.eigenvalues()[2])<<" "<<arg(es.eigenvalues()[2])<<endl;
    }
    

}
void single_run(){
	int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
	bool testing;
	string type;
	ifstream infile("params");
	infile>>Ne>>invNu;
	infile>>nWarmup>>nMeas>>nSteps>>nBins;
	infile>>seed;
	infile>>testing;
	infile>>type;
	//initialize MC object
    
    int gs=0;
	LATTICE ll(Ne,invNu, testing, type, seed, gs);
	vector<double> hole(2); hole[0]=0.1; hole[1]=0.2;
	ll.set_hole(hole);
	ofstream outfile("out"),eout("energy");
    
    ll.print_ds();
    for(int s=0;s<nBins;s++){
        
        ll.reset();
        ll.step(nWarmup);
        double E=0,E2=0;
        //double P=0,P2=0,three=0;
        double e; //,p;
        complex<double> berry_phase(0,0);
        deque<double> e_tracker, p_tracker;
        int Ntrack=10;
        vector<double> autocorr_e(Ntrack,0), autocorr_p(Ntrack,0);
        for(int i=0;i<nMeas;i++){
            ll.step(nSteps);
            e=ll.coulomb_energy();
            E+=e;
            E2+=e*e;
            //			p=log(ll.running_weight);
            //			P+=p;
            eout<<e<<endl;
            //			//autocorrelations
            //			e_tracker.push_front(e);
            //			p_tracker.push_front(p);
            //			if(i>=Ntrack){
            //				for(int j=0;j<Ntrack;j++){
            //					autocorr_e[j]+=e_tracker[j]*e;
            //					autocorr_p[j]+=p_tracker[j]*p;
            //				}
            //				e_tracker.pop_back();
            //				p_tracker.pop_back();
            //			}
            
            ll.update_structure_factors();
        }
        outfile<<E/(1.*nMeas*ll.Ne)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))/(1.*ll.Ne)<<" "<<real(berry_phase)/(1.*nMeas)<<" "<<imag(berry_phase)/(1.*nMeas)<<endl;
        cout<<"acceptance rate: "<<(1.*ll.accepts)/(1.*ll.tries)<<endl;
        
        //		ofstream auto_out("auto");
        //		for(int j=0;j<Ntrack;j++){
        //			auto_out<<j+1<<" ";
        //			auto_out<<autocorr_e[j]/(1.*(nMeas-Ntrack))<<" "<<pow(E/(1.*nMeas),2)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))<<" ";
        //			auto_out<<autocorr_p[j]/(1.*(nMeas-Ntrack))<<" "<<pow(P/(1.*nMeas),2)<<" "<<(P2/(1.*nMeas)-pow(P/(1.*nMeas),2))<<" ";
        //			auto_out<<endl;
        
        //		}
    }
    outfile<<endl;
    ll.print_structure_factors(nMeas*nBins);
    eout.close();
    outfile.close();
}

void coul_energy_CFL_dbar(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins, double* dbar_parameter){
    double sumE=0;
    edbar.change_dbar_parameter(dbar_parameter[0],dbar_parameter[1]);
    for (int s=0; s<nBins; s++) {
        edbar.reset();
        edbar.step(nWarmup);
        double E=0, E2=0, e;
        for(int i=0;i<nMeas;i++){
            edbar.step(nSteps);
            e=edbar.coulomb_energy();
            E+=e;
            E2+=e*e;
        }
        sumE+=E/(1.*nMeas*edbar.Ne);
    }
    ave_E=sumE/(1.*nBins);
}

void coul_energy_laughlin(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins){
    double sumE=0;
    for (int s=0; s<nBins; s++) {
        edbar.reset();
        edbar.step(nWarmup);
        double E=0, E2=0, e;
        for(int i=0;i<nMeas;i++){
            edbar.step(nSteps);
            e=edbar.coulomb_energy();
            E+=e;
            E2+=e*e;
        }
        sumE+=E/(1.*nMeas*edbar.Ne);
        cout<<E/(1.*nMeas*edbar.Ne)<<endl;
    }
    ave_E=sumE/(1.*nBins);
}

double phasemod(complex<double> in){
    double out=arg(in);
        if(out<0) return out+2*M_PI;
        else if (out>2*M_PI) return out-2*M_PI;
        else return out;
    return out;
}
