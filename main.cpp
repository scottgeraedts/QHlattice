#include "lattice.h"
#include "berry_phase.h"

int main(){
//    berry_phase bp(20);
//    bp.two_full_braiding();

//    berry_phase bp(16);
//    bp.two_full_braiding();
    
//	void single_run();
//	single_run();

    void two_holes(string str);
    string str;
    str = "test_hermitian";
//    str = "";
    two_holes(str);
    
//    void laughlinberryphase();
//    laughlinberryphase();
    
//    void testeigen();
//    testeigen();
    
//    void test_largesize();
//    test_largesize();
    

}

void testeigen(){
    Eigen::Matrix2cd MM;
    complex<double> ii = complex<double> (0,1);
    MM<<1,2.+ii,2.-ii,4;
    cout<<"MM=\n"<<MM<<endl;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(MM);
    cout<<"\n eigenvectors = \n"<<es.eigenvectors()<<endl;//eigenvectors are columns.
    cout<<"eigenvalues = \n"<<es.eigenvalues()<<endl;
    cout<<endl;
    cout<<"eigenvectors . eigenvectors.adjoint = \n"<<es.eigenvectors()*es.eigenvectors().adjoint()<<endl;
    cout<<"eigenvectors.adjoint . eigenvectors = \n"<<es.eigenvectors().adjoint()*es.eigenvectors()<<endl;
    cout<<endl;
    cout<<"eigenvectors.adjoit() * MM * eigenvectors = \n"<<es.eigenvectors().adjoint()*MM*es.eigenvectors()<<endl;
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
void laughlinberryphase(){
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
    
    int Grid=20;
    vector<vector<double> > holes;
    for (int i=0; i<Grid; i++) {vector<double> a(2); a[0]=0.5/(1.*Grid)*i;  a[1]=0.;     holes.push_back(a);}
    for (int i=0; i<Grid; i++) {vector<double> a(2); a[0]=0.5; a[1]=0.5/(1.*Grid)*i;     holes.push_back(a);}
    for (int i=0; i<Grid; i++) {vector<double> a(2); a[0]=0.5-0.5/(1.*Grid)*i; a[1]=0.5; holes.push_back(a);}
    for (int i=0; i<Grid; i++) {vector<double> a(2); a[0]=0.; a[1]=0.5-0.5/(1.*Grid)*i;  holes.push_back(a);}
//    cout<<"holes.size()="<<holes.size()<<endl;
//    for (int i=0; i<holes.size(); i++) {cout<<holes[i][0]<<" "<<holes[i][1]<<endl;}
    int nds=holes.size();
    
    vector<vector<double> > holes2(nds, vector<double>(2,0));
    
    int supermod(int k, int n);
    for(int i=0;i<nds;i++) holes2[supermod(i-1,nds)]=holes[i];
    
    //ll_0, ll_1, ll_2 are 3 gs whose hole given by holes[b];
    //pp_0, pp_1, pp_2 are 3 gs whose hole given by holes2[b];
    ofstream bout("berry_laughlin");
    
    vector<LATTICE> ll(3), pp(3);
    for (int i=0; i<3; i++) {
        ll[i]=LATTICE(Ne, invNu, testing, type, seed, i);
        pp[i]=LATTICE(Ne, invNu, testing, type, seed, i);
    }
    
    for(int b=0;b<nds;b++){
        Eigen::MatrixXcd berrymatrix = Eigen::MatrixXcd::Zero(3,3);
        for (int i=0; i<3; i++) {
            ll[i].set_hole(holes[b]);
            pp[i].set_hole(holes2[b]);
            ll[i].reset(); ll[i].step(nWarmup);
        }
        for(int k=0;k<nMeas;k++){
            for (int i=0; i<3; i++) ll[i].step(nSteps);
            //berrymatrix(m,n) := conj(ll_m) . pp_n = |ll_m|^2 . pp_n/ll_m;
            //generalization of: berry2+=ll2.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) berrymatrix(i,j)+=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
            }
        }
        berrymatrix/=(1.*nMeas);
        double phasemod(complex<double> in);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix);
//        bout<<holes[b][0]<<" "<<holes[b][1]<<" "<<holes2[b][0]<<" "<<holes2[b][1]<<"   "<<sqrt(norm(berrymatrix.determinant()))<<endl;
        cout<<holes[b][0]<<" "<<holes[b][1]<<" "<<holes2[b][0]<<" "<<holes2[b][1]<<"   "<<abs(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[1])<<" "<<abs(es.eigenvalues()[2])<<" "<<arg(es.eigenvalues()[2])<<endl;
    }
    
}
void two_holes(string str){
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
    while (x<0.2) {
        a[0]=x; a[1]=0.;
        holes.push_back(a);
        x+=0.02;
    }
    int nds=holes.size();
    
    
    ofstream bout("twoholelaughlin");
    vector<LATTICE> ll(invNu),ll2(invNu);
    vector<Eigen::MatrixXcd> overlaps(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
//    vector<Eigen::MatrixXcd> overlapsll(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
    Eigen::MatrixXcd overlapsll=Eigen::MatrixXcd::Zero(invNu, invNu);
    vector<Eigen::MatrixXcd> overlapsll2(nds,Eigen::MatrixXcd::Zero(invNu,invNu));
    
	for(int gs=0;gs<invNu;gs++){
		ll[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
		ll[gs].set_hole(holes[0]);
		ll2[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
    }
    
//    double energy=0.;
    complex<double> berry;
    complex<double> berryll, berryll2;
	for(int gs1=0;gs1<invNu;gs1++){    
	    ll[gs1].reset();
	    ll[gs1].step(nWarmup);
        for(int i=0;i<nMeas;i++){
            ll[gs1].step(nSteps);
            //	        energy+=ll[gs1].coulomb_energy();
            for(int gs2=0;gs2<invNu;gs2++){
                if (str=="test_hermitian") {
                    overlapsll(gs1,gs2)+=ll[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
                }
                else if(str=="" && gs2>=gs1){
                    overlapsll(gs1,gs2)+=ll[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
                }
                for(int b=0;b<nds;b++){
                    ll2[gs2].set_hole(holes[b]);
                    ll2[gs2].reset();
                    complex<double> w0,w2; w0=holes[0][0]*ll[gs1].L1+holes[0][1]*ll[gs1].L2; w2=holes[b][0]*ll[gs1].L1+holes[b][1]*ll[gs1].L2;
                    berry=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs())*exp(-1./6.*(norm(w2)-norm(w0)));//is that correct?
                    overlaps[b](gs1,gs2)+=berry;
                }
            }
        }
	}
    
    cout<<overlapsll/(1.*nMeas)<<endl;
    
    
    for (int gs1=0; gs1<invNu; gs1++) {
        for (int b=0; b<nds; b++) {
            ll2[gs1].set_hole(holes[b]);
            ll2[gs1].reset();
            ll2[gs1].step(nWarmup);
            for (int i=0; i<nMeas; i++) {
                ll2[gs1].step(nSteps);
                
                for (int gs2=0; gs2<invNu; gs2++) {
                    if (str=="test_hermitian") {
                        berryll2=ll2[gs2].get_wf(ll2[gs1].get_locs())/ll2[gs1].get_wf(ll2[gs1].get_locs());
                        overlapsll2[b](gs1,gs2)+=berryll2;
                    }
                    else if(str=="" && gs2>=gs1){
                        berryll2=ll2[gs2].get_wf(ll2[gs1].get_locs())/ll2[gs1].get_wf(ll2[gs1].get_locs());
                        overlapsll2[b](gs1,gs2)+=berryll2;
                    }
                }
                
            }
        }
    }
    
    if (str!="test_hermitian") {
        for (int gs1=0; gs1<invNu; gs1++) {
            for (int gs2=gs1-1; gs2>=0; gs2--) {
                for (int b=0; b<nds; b++) {
                    overlapsll(gs1,gs2) = conj(overlapsll(gs2,gs1));
                    overlapsll2[b](gs1,gs2) = conj(overlapsll2[b](gs2,gs1));
                }
            }
        }
    }
    
    overlapsll/=(1.*nMeas);
    for (int b=0; b<nds; b++) {
        overlaps[b]/=(1.*nMeas); overlapsll2[b]/=(1.*nMeas);
//        cout<<overlaps[b]<<endl;
//        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b]);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> esll(overlapsll);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> esll2(overlapsll2[b]);
        
        Eigen::MatrixXcd esll_diag = Eigen::MatrixXcd::Zero(invNu, invNu);
        Eigen::MatrixXcd esll2_diag = Eigen::MatrixXcd::Zero(invNu, invNu);
        
//        cout<<endl;
//        bout<<"\noverlapsll2["<<b<<"] = \n"<<overlapsll2[b]<<endl;
//        bout<<"\noverlapsll = \n"<<overlapsll<<endl;
//        cout<<"esll.eigenvector.adjoint() * esll * esll.eigenvector = \n"<<esll.eigenvectors().adjoint() * overlapsll * esll.eigenvectors()<<endl;
//        cout<<"esll2.eigenvector.adjoint() * esll2 * esll2.eigenvector = \n"<<esll2.eigenvectors().adjoint() * overlapsll2[b] * esll2.eigenvectors()<<endl;
//        cout<<endl;
        
        Eigen::MatrixXcd diaged_berry = Eigen::MatrixXcd::Zero(invNu, invNu);
        diaged_berry = esll.eigenvectors().adjoint() * overlaps[b] * esll2.eigenvectors();
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(diaged_berry);
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
    double ave_E=0.;
    
    ll.print_ds();
    for(int s=0;s<nBins;s++){
        
        ll.reset();
        ll.step(nWarmup);
        double E=0,E2=0,P=0,P2=0,three=0;
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
    //    if(out<0) return out+2*M_PI;
    //    else if (out>2*M_PI) return out-2*M_PI;
    //    else return out;
    return out;
}
