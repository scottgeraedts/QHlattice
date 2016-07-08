//#include "berry_tests.h"

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
	ofstream outfile("out"),eout("energy");
    ll.print_ds();
    for(int s=0;s<nBins;s++){
        
		ll.change_dbar_parameter(s*0.1,s*0.1);
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

void plot_CFL_coule_vsdbar(int grid){
    ofstream bout("CFL_coule_vsdbar");
    void coul_energy_CFL_dbar(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins, double* dbar_parameter);
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
    
    int Nphi=Ne*invNu; int gs=0;
    LATTICE edbar(Ne, invNu, testing, type, seed, gs);
    
    for (int i=0; i<grid*Nphi; i++) {
        for (int j=0; j<grid*Nphi; j++) {
            double ave_e, dbar_parameter[2];
            dbar_parameter[0]=1.*i/(1.*Nphi*grid); dbar_parameter[1]=1.*j/(1.*Nphi*grid);
            coul_energy_CFL_dbar(edbar, ave_e, nWarmup, nMeas, nSteps, nBins, dbar_parameter);
            bout<<dbar_parameter[0]<<" "<<dbar_parameter[1]<<" "<<ave_e<<endl;
        }
    }
    
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
//        cout<<E/(1.*nMeas*edbar.Ne)<<endl;
    }
    ave_E=sumE/(1.*nBins);
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
    //overlaps[b][0]=<psi(0)|psi(xb)>, overlaps[b][1]=<|<psi(0)|psi(xb)>|^2>.
    for (int b=0; b<nds; b++) {
        vector<Eigen::MatrixXcd> aa;
        Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(3,3);
        aa.push_back(a); aa.push_back(a);
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
    }
    
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
        for (int l=0; l<2; l++) overlaps[b][l]/=(1.*nMeas);
        for (int gs1=0; gs1<invNu; gs1++) {
            for (int gs2=0; gs2<invNu; gs2++) {
                overlaps[b][0](gs1,gs2)/=sqrt(abs(overlaps[b][1](gs1,gs2)));
            }
        }
        berryloop*=overlaps[b][0];
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

double phasemod(complex<double> in){
    double out=arg(in);
        if(out<0) return out+2*M_PI;
        else if (out>2*M_PI) return out-2*M_PI;
        else return out;
    return out;
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
        // 	    cout<<holes[b][0]<<" "<<holes[b][1]<<" "<<overlaps[b][2]<<endl;
    }
    datas.clear();
    for (int b=0; b<nds; b++) {
        data tmp;
        tmp.position[0]=holes[b][0]; tmp.position[1]=holes[b][1];
        tmp.amp[gs]=abs(overlaps[b][2]); tmp.ang[gs]=arg(overlaps[b][2]);
        datas.push_back(tmp);
    }
    
    double phase=0.;
    for (int b=0; b<nds; b++) {
        phase+=datas[b].ang[gs];
    }
    cout<<"\n\nphase = "<<phase<<endl;
    cout<<endl;
    
}
