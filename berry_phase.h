#ifndef BERRY_PHASE_H
#define BERRY_PHASE_H


//class berry_phase{
//    friend class LATTICE;
//private:
//    int NPhi, invNu, nWarmup, nMeas, nSteps, nBins, seed;
//    int Ne;
//    bool testing;
//    string type;
//    
//    vector <vector <int> > tempds;
//    
//public:
//    berry_phase(int);
//    int nds;
//    vector<int> dpoint;
//    vector<vector<int> > dcenter, d2, d3;
//    double energy;
//    complex<double> berry2,berry3;
//    
//    double phasemod(complex<double>);
//    
//    void two_full_braiding();
//    void two_half_braiding();
//    void three_full_braiding();
//    void laughlinberryphase();
//};
//berry_phase::berry_phase(int nds_t){
//    ifstream infile("params");
//    infile>>Ne>>invNu;
//    infile>>nWarmup>>nMeas>>nSteps>>nBins;
//    infile>>seed;
//    infile>>testing;
//    infile>>type;
//    
//    if (type=="CFL") {
//        NPhi=Ne*invNu;
//    }
//    
//    if(NPhi!=78){
//        cout<<"the berry phase calculator only works for 78 flux quanta right now!"<<endl;
//        exit(0);
//    }
//    
////    if(NPhi!=46){
////        cout<<"the berry phase calculator only works for 46 flux quanta right now!"<<endl;
////        exit(0);
////    }
//    
//    dpoint=vector<int>(2);
//    nds=nds_t;
//    
//    dcenter=vector <vector<int> >(nds,vector<int>(2));
//    d2=vector <vector<int> >(nds,vector<int>(2));
//    d3=vector <vector<int> >(nds,vector<int>(2));
//    
//    vector<int>a(2);
//    a[0]=4; a[1]=0;
////    a[]={4,0};
////    dcenter={4,0};
//    dcenter[0]=a;
////    dcenter.push_back(a);
//    
//    a[0]=4; a[1]=1;
//    dcenter[1]=a;
//    a[0]=3; a[1]=2;
//    dcenter[2]=a;
//    a[0]=2; a[1]=3;
//    dcenter[3]=a;
//    a[0]=1; a[1]=4;
//    dcenter[4]=a;
//    a[0]=0; a[1]=4;
//    dcenter[5]=a;
//    a[0]=-1; a[1]=4;
//    dcenter[6]=a;
//    a[0]=-2; a[1]=3;
//    dcenter[7]=a;
//    a[0]=-3; a[1]=2;
//    dcenter[8]=a;
//    a[0]=-4; a[1]=1;
//    dcenter[9]=a;
//    a[0]=-4; a[1]=0;
//    dcenter[10]=a;
//    a[0]=-4; a[1]=-1;
//    dcenter[11]=a;
//    a[0]=-3; a[1]=-2;
//    dcenter[12]=a;
//    a[0]=-2; a[1]=-3;
//    dcenter[13]=a;
//    a[0]=-1; a[1]=-4;
//    dcenter[14]=a;
//    a[0]=0; a[1]=-4;
//    dcenter[15]=a;
//    a[0]=1; a[1]=-4;
//    dcenter[16]=a;
//    a[0]=2; a[1]=-3;
//    dcenter[17]=a;
//    a[0]=3; a[1]=-2;
//    dcenter[18]=a;
//    a[0]=4; a[1]=-1;
//    dcenter[19]=a;
//    
//    
////    dcenter[0]={4,0};
////    dcenter[1]={4,1};
////    cout<<"dcenter[1]={4,1}"<<dcenter[1][0]<<" "<<dcenter[1][1]<<endl;
////    dcenter[2]={3,2};
////    dcenter[3]={2,3};
////    dcenter[4]={1,4};
////    dcenter[5]={0,4};
////    dcenter[6]={-1,4};
////    dcenter[7]={-2,3};
////    dcenter[8]={-3,2};
////    dcenter[9]={-4,1};
////    dcenter[10]={-4,0};
////    dcenter[11]={-4,-1};
////    dcenter[12]={-3,-2};
////    dcenter[13]={-2,-3};
////    dcenter[14]={-1,-4};
////    dcenter[15]={0,-4};
////    dcenter[16]={1,-4};
////    dcenter[17]={2,-3};
////    dcenter[18]={3,-2};
////    dcenter[19]={4,-1};
//
////    dcenter[0]={3,0};
////    dcenter[1]={3,1};
////    dcenter[2]={2,2};
////    dcenter[3]={1,3};
////    dcenter[4]={0,3};
////    dcenter[5]={-1,3};
////    dcenter[6]={-2,2};
////    dcenter[7]={-3,1};
////    dcenter[8]={-3,0};
////    dcenter[9]={-3,-1};
////    dcenter[10]={-2,-2};
////    dcenter[11]={-1,-3};
////    dcenter[12]={0,-3};
////    dcenter[13]={1,-3};
////    dcenter[14]={2,-2};
////    dcenter[15]={3,-1};
//    
//    int supermod(int k, int n);
//    for(int i=0;i<nds;i++){
//        d2[supermod(i-1,nds)]=dcenter[i];
//        d3[supermod(i+1,nds)]=dcenter[i];
//    }
//}
//double berry_phase::phasemod(complex<double> in){
//    double out=arg(in);
////    if(out<0) return out+2*M_PI;
////    else if (out>2*M_PI) return out-2*M_PI;
////    else return out;
//    return out;
//}
//void berry_phase::two_full_braiding(){
//    cout<<"start working"<<endl;
//    double tot_berry_phase2=0.;
//    double tot_berry_phase3=0.;
//    
//    ofstream bout("berry");
//    LATTICE ll(Ne,invNu, testing, type, seed);
//    LATTICE ll2(Ne,invNu, testing, type, seed);
////    LATTICE ll3(Ne,invNu, testing, type, seed);
//    
//    for(int b=0;b<nds;b++){
//        //stuff for initializing with a custom set of ds (for the Berry phase calculation)
//        double center_frac[2]={0.,0.}; //might need to change this for different N, but it works for N=22
//        ll.make_fermi_surface(center_frac, ll.Ne-2);
//        tempds=ll.get_ds();//ds are defined on Ne lattice.
//        ll.print_ds();
//        dpoint[0]=dcenter[b][0]; dpoint[1]=dcenter[b][1];
//        tempds.push_back(dpoint);
//        dpoint[0]=-dcenter[b][0]; dpoint[1]=-dcenter[b][1];
//        tempds.push_back(dpoint);
//        ll.set_ds(tempds);
//        
//        ll2.make_fermi_surface(center_frac,ll.Ne-2);
//        tempds=ll2.get_ds();
//        dpoint[0]=d2[b][0]; dpoint[1]=d2[b][1];
//        tempds.push_back(dpoint);
//        dpoint[0]=-d2[b][0]; dpoint[1]=-d2[b][1];
//        tempds.push_back(dpoint);
//        ll2.set_ds(tempds);
//        
////        ll3.make_fermi_surface(center_frac,ll.Ne-2);
////        tempds=ll3.get_ds();
////        dpoint[0]=d3[b][0]; dpoint[1]=d3[b][1];
////        tempds.push_back(dpoint);
////        dpoint[0]=-d3[b][0]; dpoint[1]=-d3[b][1];
////        tempds.push_back(dpoint);
////        ll3.set_ds(tempds);
//        
//        energy=0;
//        berry2=0; berry3=0;
//        ll.step(nWarmup);
//        complex<double> btemp;
//        for(int i=0;i<nMeas;i++){
//            ll.step(nSteps);
//            energy+=ll.coulomb_energy();
//            berry2+=ll2.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
//            btemp=ll2.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
//			cout<<"output:"<<ll2.get_wf(ll.get_locs())<<" "<<ll.get_wf(ll.get_locs())<<" "<<abs(btemp)<<" "<<arg(btemp)<<endl;
////            berry3+=ll3.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
//        }
//        bout<<dcenter[b][0]<<" "<<dcenter[b][1]<<" "<<d2[b][0]<<" "<<d2[b][1]<<"   "<<abs(berry2)/(1.*nMeas)<<"   "<<phasemod(berry2)<<"   "<<energy/(1.*nMeas*ll.Ne)<<endl;
////        bout<<dcenter[b][0]<<" "<<dcenter[b][1]<<" "<<d3[b][0]<<" "<<d3[b][1]<<"   "<<abs(berry3)/(1.*nMeas)<<"   "<<phasemod(berry3)<<"   "<<energy/(1.*nMeas*ll.Ne)<<endl;
////        tot_berry_phase2+=phasemod(berry2);
////        tot_berry_phase3+=phasemod(berry3);
//
//    }
//    
////    cout<<"total berry phase2 = "<<tot_berry_phase2<<" total berry phase3 = "<<tot_berry_phase3<<endl;
//}


//----------laughlin berry phase----------

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
    for(int i=0;i<x.rows();i++){
        for(int j=i+1;j<x.rows();j++){
            x(i,j)=0.5*(x(i,j)+conj(x(j,i)));
            x(j,i)=conj(x(i,j));
        }
    }
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


//void laughlin_bp_single_state(int gs, vector<double> length, double steplength, vector<data> &datas){
//    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
//    bool testing;
//    string type;
//    ifstream infile("params");
//    infile>>Ne>>invNu;
//    infile>>nWarmup>>nMeas>>nSteps>>nBins;
//    infile>>seed;
//    infile>>testing;
//    infile>>type;
//    //initialize MC object
//    
//    vector<vector<double> > holes; vector<int> Grid(2);
//    for (int i=0; i<2; i++) {Grid[i]=(int)(length[i]/steplength);}
//    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=steplength*i; a[1]=0.;                  holes.push_back(a);}
//    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=steplength*i; a[0]=length[0];           holes.push_back(a);}
//    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=length[0]-steplength*i; a[1]=length[1]; holes.push_back(a);}
//    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=length[1]-steplength*i; a[0]=0.;        holes.push_back(a);}
//    int nds=holes.size();
//    vector<vector<double> > holes2(nds, vector<double>(2,0));
//    int supermod(int k, int n);
//    for(int i=0;i<nds;i++) holes2[supermod(i-1,nds)]=holes[i];//(holes[b],holes2[b]) = (holes[b],holes[b+1]).
//    
//    LATTICE ll(Ne, invNu, testing, type, seed, gs), pp(Ne, invNu, testing, type, seed, gs);
//    vector<vector<complex<double> > > overlaps;
//    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|psi(xb)|psi(xb+1)|^2>, overlaps[b][2](i,j)=overlaps[b][0](i,j)/sqrt{overlaps[b][1](i,j)}.
//    for (int b=0; b<nds; b++) {
//        vector<complex<double> > aa;
//        aa.push_back(0.); aa.push_back(0.); aa.push_back(0.);
//        overlaps.push_back(aa);
//    }
//    
//    for(int b=0;b<nds;b++){
//        ll.set_hole(holes[b]); pp.set_hole(holes2[b]);
//        ll.reset(); ll.step(nWarmup);
//        for(int k=0;k<nMeas;k++){
//            ll.step(nSteps);
//            complex<double> temp=pp.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
//            overlaps[b][0]+=temp; overlaps[b][1]+=norm(temp);
//        }
//        for (int l=0; l<3; l++) {
//            overlaps[b][l]/=(1.*nMeas);
//        }
//        overlaps[b][2] = overlaps[b][0]/sqrt(overlaps[b][1]);
//        // 	    cout<<holes[b][0]<<" "<<holes[b][1]<<" "<<overlaps[b][2]<<endl;
//    }
//    datas.clear();
//    for (int b=0; b<nds; b++) {
//        data tmp;
//        tmp.position[0]=holes[b][0]; tmp.position[1]=holes[b][1];
//        tmp.amp[gs]=abs(overlaps[b][2]); tmp.ang[gs]=arg(overlaps[b][2]);
//        datas.push_back(tmp);
//    }
//    
//    double phase=0.;
//    for (int b=0; b<nds; b++) {
//        phase+=datas[b].ang[gs];
//    }
//    cout<<"\n\nphase = "<<phase<<endl;
//    cout<<endl;
//    
//}
//
//void two_holes(string str, int nmeasurement, data& test){
//    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
//    bool testing;
//    string type;
//    ifstream infile("params");
//    infile>>Ne>>invNu;
//    infile>>nWarmup>>nMeas>>nSteps>>nBins;
//    infile>>seed;
//    infile>>testing;
//    infile>>type;
//    //initialize MC object
//    
//    vector<vector<double> > holes;
//    double x=0.;
//    vector<double> a(2);
//    while (x<0.05) {
//        a[0]=x; a[1]=0.;
//        holes.push_back(a);
//        x+=0.005;
//    }
//    int nds=holes.size();
//    
//    if (str=="test") {
//        nMeas=nmeasurement;
//        nds=1;
//    }
//    
//    ofstream bout("twoholelaughlinnew");
//    vector<LATTICE> ll(invNu),ll2(invNu);//ll is psi(x), ll2 is psi(x).
//    vector<vector<Eigen::MatrixXcd> > overlaps;
//    //overlaps[b][0]=<psi(0)|psi(xb)>, overlaps[b][1]=<|<psi(0)|psi(xb)>|^2>.
//    for (int b=0; b<nds; b++) {
//        vector<Eigen::MatrixXcd> aa;
//        Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(3,3);
//        aa.push_back(a); aa.push_back(a);
//        overlaps.push_back(aa);
//    }
//    
//    for(int gs=0;gs<invNu;gs++){
//        ll[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
//        ll[gs].set_hole(holes[0]);
//        ll2[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
//    }
//    
//    for(int gs1=0;gs1<invNu;gs1++){
//        ll[gs1].reset();
//        ll[gs1].step(nWarmup);
//        for(int i=0;i<nMeas;i++){
//            ll[gs1].step(nSteps);
//            for(int gs2=0;gs2<invNu;gs2++){
//                for(int b=0;b<nds;b++){
//                    ll2[gs2].set_hole(holes[b]);
//                    ll2[gs2].reset();
//                    complex<double> temp=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
//                    overlaps[b][0](gs1,gs2)+=temp;
//                    overlaps[b][1](gs1,gs2)+=norm(temp);
//                }
//            }
//        }
//    }
//    
//    /*
//     for (int b=0; b<nds; b++) {
//     for (int gs1=0;gs1<invNu;gs1++) {
//     ll[gs1].reset();
//     ll[gs1].step(nWarmup);
//     for (int gs2=0; gs2<invNu; gs2++) {
//     ll2[gs2].set_hole(holes[b]);
//     ll2[gs2].reset();
//     for (int i=0;i<nMeas;i++) {
//     ll[gs1].step(nSteps);
//     complex<double> temp=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
//     overlaps[b][0](gs1,gs2)+=temp;
//     overlaps[b][1](gs1,gs2)+=norm(temp);
//     }
//     }
//     }
//     }
//     cout<<"\noverlaps[10][0]=\n"<<overlaps[10][0]<<endl;
//     cout<<"\noverlaps[10][1]=\n"<<overlaps[10][1]<<endl;
//     */
//    
//    Eigen::Matrix3cd berryloop = Eigen::Matrix3cd::Identity(3,3);
//    for (int b=0; b<nds; b++) {
//        for (int l=0; l<2; l++) overlaps[b][l]/=(1.*nMeas);
//        for (int gs1=0; gs1<invNu; gs1++) {
//            for (int gs2=0; gs2<invNu; gs2++) {
//                overlaps[b][0](gs1,gs2)/=sqrt(abs(overlaps[b][1](gs1,gs2)));
//            }
//        }
//        berryloop*=overlaps[b][0];
//    }
//    
//    //    cout<<"\nberryloop=\n"<<berryloop<<endl;
//    //    cout<<"determinant=\n"<<arg(berryloop.determinant())<<endl;
//    
//    for (int b=0; b<nds; b++) {
//        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b][2]);
//        bout<<holes[b][0]<<" "<<holes[b][1]<<" "<<abs(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[1])<<" "<<abs(es.eigenvalues()[2])<<" "<<arg(es.eigenvalues()[2])<<endl;
//        //        for (int i=0; i<3; i++) {
//        //            test.amp[i]=abs(es.eigenvalues()[i]);
//        //            test.ang[i]=arg(es.eigenvalues()[i]);
//        //        }
//    }
//    
//}
//
//void two_holes_scott(){
//    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
//    bool testing;
//    string type;
//    ifstream infile("params");
//    infile>>Ne>>invNu;
//    infile>>nWarmup>>nMeas>>nSteps>>nBins;
//    infile>>seed;
//    infile>>testing;
//    infile>>type;
//    //initialize MC object
//    
//    vector<vector<double> > holes;
//    double x=0.;
//    vector<double> a(2);
//    while (x<0.01) {
//        a[0]=x; a[1]=0.;
//        holes.push_back(a);
//        x+=0.02;
//    }
//    int nds=holes.size();
//    
//    
//    ofstream bout("twoholelaughlin");
//    vector<LATTICE> ll(invNu),ll2(invNu);
//    vector<Eigen::MatrixXcd> overlaps12(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
//    vector<Eigen::MatrixXcd> overlaps12_2(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
//    Eigen::MatrixXcd overlaps11=Eigen::MatrixXcd::Zero(invNu, invNu);
//    Eigen::MatrixXcd overlaps11_2=Eigen::MatrixXcd::Zero(invNu, invNu);
//    
//    vector<Eigen::MatrixXcd> overlaps22(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
//    vector<Eigen::MatrixXcd> overlaps22_2(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
//    
//    for(int gs=0;gs<invNu;gs++){
//        ll[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
//        ll[gs].set_hole(holes[0]);
//        ll2[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
//    }
//    
//    complex<double> berry;
//    for(int gs1=0;gs1<invNu;gs1++){
//        ll[gs1].reset();
//        ll[gs1].step(nWarmup);
//        for(int i=0;i<nMeas;i++){
//            ll[gs1].step(nSteps);
//            //	        energy+=ll[gs1].coulomb_energy();
//            for(int gs2=0;gs2<invNu;gs2++){
//                berry=ll[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
//                overlaps11(gs1,gs2)+=berry;
//                overlaps11_2(gs1,gs2)+=norm(berry);
//                for(int b=0;b<nds;b++){
//                    ll2[gs2].set_hole(holes[b]);
//                    ll2[gs2].reset();
//                    berry=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
//                    overlaps12[b](gs1,gs2)+=berry;
//                    overlaps12_2[b](gs1,gs2)+=norm(berry);
//                }
//            }
//        }
//    }
//    cout<<"energies:"<<endl;
//    vector<double> energy(invNu,0);
//    for (int gs1=0; gs1<invNu; gs1++) {
//        for (int b=0; b<nds; b++) {
//            
//            ll2[gs1].set_hole(holes[b]);
//            ll2[gs1].reset();
//            ll2[gs1].step(nWarmup);
//            for (int i=0; i<nMeas; i++) {
//                ll2[gs1].step(nSteps);
//                energy[gs1]+=ll2[gs1].coulomb_energy();
//                for (int gs2=0; gs2<invNu; gs2++) {
//                    berry=ll2[gs2].get_wf(ll2[gs1].get_locs())/ll2[gs1].get_wf(ll2[gs1].get_locs());
//                    overlaps22[b](gs1,gs2)+=berry;
//                    overlaps22_2[b](gs1,gs2)+=norm(berry);
//                }
//                
//            }
//        }
//        cout<<energy[gs1]/(1.*nMeas)<<endl;
//    }
//    
//    overlaps11/=(1.*nMeas);
//    overlaps11_2/=(1.*nMeas);
//    overlaps11=overlaps11.cwiseQuotient(overlaps11_2.cwiseSqrt());
//    //	hermitianize(overlaps11);
//    
//    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es11(overlaps11);
//    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es,es22;
//    Eigen::MatrixXcd diaged_berry;
//    for (int b=0; b<nds; b++) {
//        overlaps12[b]/=(1.*nMeas);
//        overlaps22[b]/=(1.*nMeas);
//        overlaps12_2[b]/=(1.*nMeas);
//        overlaps22_2[b]/=(1.*nMeas);
//        //		cout<<overlaps12[b]<<endl<<endl;
//        //		cout<<overlaps12_2[b]<<endl<<endl;
//        
//        overlaps12[b]=overlaps12[b].cwiseQuotient(overlaps12_2[b].cwiseSqrt());
//        overlaps22[b]=overlaps22[b].cwiseQuotient(overlaps22_2[b].cwiseSqrt());
//        //		hermitianize(overlaps12[b]);
//        //		hermitianize(overlaps22[b]);
//        cout<<overlaps12[b]<<endl<<endl;
//        cout<<overlaps11<<endl<<endl;
//        cout<<overlaps12[b]<<endl<<endl;
//        //		cout<<"------------------------------------"<<endl;
//        
//        es22.compute(overlaps22[b]);
//        
//        //
//        //        Eigen::MatrixXcd esll_diag = Eigen::MatrixXcd::Zero(invNu, invNu);
//        //        Eigen::MatrixXcd esll2_diag = Eigen::MatrixXcd::Zero(invNu, invNu);
//        //        esll_diag=chop(esll.eigenvectors().adjoint() * overlapsll * esll.eigenvectors());
//        //        esll2_diag=chop(esll2.eigenvectors().adjoint() * overlapsll2[b] * esll2.eigenvectors());
//        //
//        //        cout<<"\nb="<<b<<endl;
//        //        cout<<"\noverlapsll = \n"<<overlapsll<<endl;
//        //        cout<<"\noverlapsll2["<<b<<"] = \n"<<overlapsll2[b]<<endl;
//        //        cout<<"esll.eigenvector.adjoint() * esll * esll.eigenvector = \n"<<esll_diag<<endl;
//        //        cout<<"esll2.eigenvector.adjoint() * esll2 * esll2.eigenvector = \n"<<esll2_diag<<endl;
//        //        cout<<endl;
//        //
//        Eigen::MatrixXcd eigvecs=es22.eigenvectors();
//        cout<<eigvecs<<endl;
//        diaged_berry = es22.eigenvectors().adjoint() * overlaps12[b] *  es22.eigenvectors();
//        //        cout<<"\noverlaps[b]="<<overlaps[b]<<endl;
//        cout<<"diaged_berry="<<endl<<diaged_berry<<endl;
//        es.compute(overlaps12[b]);
//        bout<<holes[b][0]<<" "<<holes[b][1]<<" "<<abs(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[1])<<" "<<abs(es.eigenvalues()[2])<<" "<<arg(es.eigenvalues()[2])<<endl;
//    }
//}

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
    LATTICE ll(Ne, invNu, testing, type, seed, gs);
    ofstream outfile("out");
    ofstream eout("energy");
    ll.print_ds();
    for(int s=0;s<nBins;s++){
        ll.change_dbar_parameter(s*0.1,s*0.1);
        ll.reset();
        ll.step(nWarmup);
        double E=0,E2=0;
        double P=0,P2=0,three=0;
        double e,p;
        complex<double> berry_phase(0,0);
        deque<double> e_tracker, p_tracker;
        int Ntrack=10;
        vector<double> autocorr_e(Ntrack,0), autocorr_p(Ntrack,0);
        for(int i=0;i<nMeas;i++){
            ll.step(nSteps);
            e=ll.coulomb_energy();
            E+=e;
            E2+=e*e;
            //            p=log(ll.running_weight); //this is a bug since changes are made in lattice.cpp, caus running_weight is logarithm of psi square.
            p=ll.running_weight;
            P+=p;
            P2+=p*p;
            eout<<e<<endl;
            
            //autocorrelations
            e_tracker.push_front(e);
            p_tracker.push_front(p);
            if(i>=Ntrack){
                for(int j=0;j<Ntrack;j++){
                    autocorr_e[j]+=e_tracker[j]*e;
                    autocorr_p[j]+=p_tracker[j]*p;
                }
                e_tracker.pop_back();
                p_tracker.pop_back();
            }
            
            //structure factor
            //            ll.update_structure_factors();
        }
        outfile<<E/(1.*nMeas*ll.Ne)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))/(1.*ll.Ne)<<" "<<real(berry_phase)/(1.*nMeas)<<" "<<imag(berry_phase)/(1.*nMeas)<<endl;
        cout<<"acceptance rate: "<<(1.*ll.accepts)/(1.*ll.tries)<<endl;
        
        ofstream auto_out("auto");
        for(int j=0;j<Ntrack;j++){
            auto_out<<j+1<<" ";
            auto_out<<autocorr_e[j]/(1.*(nMeas-Ntrack))<<" "<<pow(E/(1.*nMeas),2)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))<<" ";
            auto_out<<autocorr_p[j]/(1.*(nMeas-Ntrack))<<" "<<pow(P/(1.*nMeas),2)<<" "<<(P2/(1.*nMeas)-pow(P/(1.*nMeas),2))<<" ";
            auto_out<<endl;
        }
        
    }
    
    outfile<<endl;
    //    ll.print_structure_factors(nMeas*nBins);
    eout.close();
    outfile.close();
    
}
void single_run_jie(){
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
    LATTICE ll(Ne, invNu, testing, type, seed, gs);
    ofstream outfile("out");
    ofstream eout("energy");
    ll.print_ds();
    
    ofstream auto_out("auto_jie");
    for (int step=1; step<=10; step+=1) {
        nSteps=step;
//        ll.change_dbar_parameter(s*0.1,s*0.1);
        ll.reset();
        ll.step(nWarmup);
        double E=0,E2=0;
        double P=0,P2=0,three=0;
        double e,p;
        complex<double> berry_phase(0,0);
        deque<double> e_tracker, p_tracker;
        int Ntrack=10;
        vector<double> autocorr_e(Ntrack,0), autocorr_p(Ntrack,0);
        
        for(int i=0;i<nMeas;i++){
            ll.step(nSteps);
            e=ll.coulomb_energy();
            E+=e;
            E2+=e*e;
            //            p=log(ll.running_weight); //this is a bug since changes are made in lattice.cpp, caus running_weight is logarithm of psi square.
            p=ll.running_weight;
            P+=p;
            P2+=p*p;
            eout<<e<<endl;
            
            //autocorrelations
            e_tracker.push_front(e);
            p_tracker.push_front(p);
            if(i>=Ntrack){
                for(int j=0;j<Ntrack;j++){
                    autocorr_e[j]+=e_tracker[j]*e;
                    autocorr_p[j]+=p_tracker[j]*p;
                }
                e_tracker.pop_back();
                p_tracker.pop_back();
            }
            
            //structure factor
            //            ll.update_structure_factors();
        }
        
        //output auto
        auto_out<<step<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))<<" "<<(P2/(1.*nMeas)-pow(P/(1.*nMeas),2))<<endl;
    }
    
    outfile<<endl;
    //    ll.print_structure_factors(nMeas*nBins);
    eout.close();
    outfile.close();
    
}

//void plot_CFL_coule_vsdbar(int grid){
//    ofstream bout("CFL_coule_vsdbar");
//    void coul_energy_CFL_dbar(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins, double* dbar_parameter);
//    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
//    bool testing;
//    string type;
//    ifstream infile("params");
//    infile>>Ne>>invNu;
//    infile>>nWarmup>>nMeas>>nSteps>>nBins;
//    infile>>seed;
//    infile>>testing;
//    infile>>type;
//    //initialize MC object
//    
//    int Nphi=Ne*invNu; int gs=0;
//    LATTICE edbar(Ne, invNu, testing, type, seed, gs);
//    
//    for (int i=0; i<grid*Nphi; i++) {
//        for (int j=0; j<grid*Nphi; j++) {
//            double ave_e, dbar_parameter[2];
//            dbar_parameter[0]=1.*i/(1.*Nphi*grid); dbar_parameter[1]=1.*j/(1.*Nphi*grid);
//            coul_energy_CFL_dbar(edbar, ave_e, nWarmup, nMeas, nSteps, nBins, dbar_parameter);
//            bout<<dbar_parameter[0]<<" "<<dbar_parameter[1]<<" "<<ave_e<<endl;
//        }
//    }
//    
//}
//
//void coul_energy_CFL_dbar(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins, double* dbar_parameter){
//    double sumE=0;
//    edbar.change_dbar_parameter(dbar_parameter[0],dbar_parameter[1]);
//    for (int s=0; s<nBins; s++) {
//        edbar.reset();
//        edbar.step(nWarmup);
//        double E=0, E2=0, e;
//        for(int i=0;i<nMeas;i++){
//            edbar.step(nSteps);
//            e=edbar.coulomb_energy();
//            E+=e;
//            E2+=e*e;
//        }
//        sumE+=E/(1.*nMeas*edbar.Ne);
//    }
//    ave_E=sumE/(1.*nBins);
//}
//
//void coul_energy_laughlin(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins){
//    double sumE=0;
//    for (int s=0; s<nBins; s++) {
//        edbar.reset();
//        edbar.step(nWarmup);
//        double E=0, E2=0, e;
//        for(int i=0;i<nMeas;i++){
//            edbar.step(nSteps);
//            e=edbar.coulomb_energy();
//            E+=e;
//            E2+=e*e;
//        }
//        sumE+=E/(1.*nMeas*edbar.Ne);
//        //        cout<<E/(1.*nMeas*edbar.Ne)<<endl;
//    }
//    ave_E=sumE/(1.*nBins);
//}

        
double phasemod(complex<double> in){
    double out=arg(in);
    if(out<0) return out+2*M_PI;
    else if (out>2*M_PI) return out-2*M_PI;
    else return out;
    return out;
}

//---------- parallel programming code ------//
void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne, int num_core){
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
    
    //if change_nMeas==0/ change_Ne==0, nMeas/ Ne is set according to 'params'.
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
    
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
    vector<vector<Eigen::MatrixXcd > > overlaps(nds, vector<Eigen::MatrixXcd>(4, Eigen::MatrixXcd::Zero(3, 3)));
    
    omp_set_num_threads(num_core);
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(3)), pp(num_core, vector<LATTICE>(3));//do this to avoid wrong memory access since openmp share memory.
    for (int k=0; k<num_core; k++) for (int i=0; i<3; i++) {ll[k][i]=LATTICE(Ne, invNu, testing, type, seed, i); pp[k][i]=LATTICE(Ne, invNu, testing, type, seed, i);}
    
    //parallel programming begin.
#pragma omp parallel for
    for(int b=0; b<nds; b++) {
        int coren = omp_get_thread_num();
        for (int i=0; i<3; i++) {
            ll[coren][i].set_hole(holes[b]);
            pp[coren][i].set_hole(holes2[b]);
            ll[coren][i].reset(); ll[coren][i].step(nWarmup);
            pp[coren][i].reset();
        }
        
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<3; i++) ll[coren][i].step(nSteps);
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    vector<complex<double> > temp(2);
                    temp[0]=pp[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                    overlaps[b][0](i,j)+=temp[0];
                    overlaps[b][1](i,j)+=norm(temp[0]);
                    temp[1]=ll[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                    overlaps[b][2](i,j)+=temp[1];
                    overlaps[b][3](i,j)+=norm(temp[1]);
                }
            }
        }
        for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
        overlaps[b][0]=overlaps[b][0].array()/overlaps[b][1].array().sqrt();
        overlaps[b][2]=overlaps[b][2].array()/overlaps[b][3].array().sqrt();
        
        hermitianize(overlaps[b][2]);
    }
    //parallel programming end.
    
    
    vector<Eigen::MatrixXcd> berrymatrix_step(nds);
    for (int b=0; b<nds; b++) berrymatrix_step[b] = overlaps[b][2].inverse() * overlaps[b][0];
    
    Eigen::MatrixXcd berrymatrix_integral = Eigen::MatrixXcd::Identity(invNu, invNu);
    vector<double> phases(invNu, 0.);
    datas.clear();
    for (int b=0; b<nds; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step[b]);
        data tmp;
        berrymatrix_integral *= berrymatrix_step[b];
        for (int i=0; i<invNu; i++) {
            phases[i]+=arg(es.eigenvalues()[i]);
            tmp.num = b; tmp.amp[i] = abs(es.eigenvalues()[i]); tmp.ang[i] = arg(es.eigenvalues()[i]);
        }
        
        // dfromnorm. calculates deviation from normality.
        double normeigenvalue=0., normmatrix=0.;
        for (int i=0; i<invNu; i++) {
            normeigenvalue+=sqrt(norm(es.eigenvalues()[i]));
        }
        for (int i=0; i<invNu; i++) {
            for (int j=0; j<invNu; j++) {
                normmatrix+=sqrt(norm(berrymatrix_step[b](i,j)));
            }
        }
        tmp.dfromnorm=normmatrix-normeigenvalue;
        
        datas.push_back(tmp);
    }
    
    double avephase;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_integral);
    datas[0].ang_trace = arg(berrymatrix_integral.trace());
    datas[0].det = arg(berrymatrix_integral.determinant());
    cout<<"\n\n Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<endl;
    cout<<"phase sum = "; for (int i=0; i<invNu; i++) {cout<<phases[i]<<" "; avephase+=phases[i]/(1.*invNu);} cout<<"\nphase average = "<<avephase<<endl;
    cout<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
    cout<<"amp(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<abs(es.eigenvalues()[i])<<" "; cout<<endl;
    cout<<"arg(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<arg(es.eigenvalues()[i])<<" ";cout<<endl;
    avephase=0.; for (int i=0; i<invNu; i++) avephase+=arg(es.eigenvalues()[i])/(1.*invNu); cout<<"sum arg(berrymatrix_integral.eigenvalue) = "<<avephase<<endl;
    cout<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
    cout<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
    cout<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl;
    
}

void test_error(int ne, double loop, double steplength, int nMea, int ncore, string test, int num_core){
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne, int num_core);
    vector<double> length(2); vector<data> datas;
    length[0]=loop; length[1]=loop;
    
    ofstream bout("test");
    if (test=="steplength") {
        for (double steplength=0.001; steplength<loop; steplength+=0.002) {
            laughlinberryphase(length, steplength, datas, nMea, ne, num_core);
            vector<double> phase(3);
            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
            bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
        }
    }
    if (test=="ne") {
        for (int ne=2; ne<22; ne=ne+2) {
            laughlinberryphase(length, steplength, datas, nMea, ne, num_core);
            vector<double> phase(3);
            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
            bout<<ne<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
        }
    }
    if (test=="nMea") {
        for (int nMeas=10; nMeas<20; nMeas+=2){
            laughlinberryphase(length, steplength, datas, nMea, ne, num_core);
            vector<double> phase(3);
            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
            bout<<nMeas<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
        }
    }
    if (test=="loop") {
        for (double x=0.05; x<0.8; x+=0.05) {
            length[0]=x; length[1]=x;
            laughlinberryphase(length, steplength, datas, nMea, ne, num_core);
            vector<double> phase(3);
            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
            bout<<x<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
        }
    }
    if (test=="normality") {
        for (double steplength=0.01; steplength<0.25; steplength+=0.01) {
            laughlinberryphase(length, steplength, datas, 0, 0, num_core);
            
            bout<<steplength<<endl;
            for (int i=0; i<datas.size(); i++) {
                bout<<datas[i].num<<" "<<datas[i].dfromnorm<<endl;
            }
            bout<<endl;
        }
    }
    
}

#endif