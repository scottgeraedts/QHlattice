#include "lattice.h"
#include "berry_phase.h"

int main(){
//    berry_phase bp(20);
//    bp.two_full_braiding();

//    berry_phase bp(16);
//    bp.two_full_braiding();
    
//	void single_run();
//	single_run();

    void laughlinberryphase();
    laughlinberryphase();
    
//    
//    int NPhi,invNu,nWarmup,nMeas,nSteps,nBins,seed;
//    bool testing;
//    string type;
//    ifstream infile("params");
//    infile>>NPhi>>invNu;
//    infile>>nWarmup>>nMeas>>nSteps>>nBins;
//    infile>>seed;
//    infile>>testing;
//    infile>>type;
//    //initialize MC object
//    
//    LATTICE ll(NPhi,invNu, testing, type, seed);
    
 /*   
	int NPhi,invNu,nWarmup,nMeas,nSteps,nBins,seed;
	bool testing;
	string type;
	ifstream infile("params");
	infile>>NPhi>>invNu; 
	infile>>nWarmup>>nMeas>>nSteps>>nBins;
	infile>>seed;
	infile>>testing;
	infile>>type;
	//initialize MC object

    ofstream out_co_dbar("co_dbar");
    
    double Ne=1.*NPhi/(1.*invNu);
    double dbar_parameter[2] = {0., 1.};
    dbar_parameter[0]=1; dbar_parameter[1]=0;
	LATTICE ll(NPhi,invNu, testing, type, seed);

     //test coulomb energy for laughlim m=3 state. For 4 particles, the energy should be -0.4141710479.
    double eval=0.;
    void coul_energy_laughlin(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins);
    coul_energy_laughlin(ll, eval, nWarmup, nMeas, nSteps, nBins);
    cout<<"laughlin state coulomb energy is"<<eval<<endl;
   */  
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
    
    vector<vector<double> > holes;
    double x=0., y=0.;
    while (x<0.5) {
        vector<double> a(2);
        a[0]=x; a[1]=0.;
        holes.push_back(a);
        x+=0.05;
    }
    y=0.05;
    while (y<0.5) {
        vector<double> a(2);
        a[0]=0.5; a[1]=y;
        holes.push_back(a);
        y+=0.05;
    }
    x=0.45;
    while (x>0.) {
        vector<double> a(2);
        a[0]=x; a[1]=0.5;
        holes.push_back(a);
        x-=0.05;
    }
    y=0.45;
    while (y>0.) {
        vector<double> a(2);
        a[0]=0.; a[1]=y;
        holes.push_back(a);
        y-=0.05;
    }
    holes.pop_back();
//    cout<<"holes.size()="<<holes.size()<<endl;
//    for (int i=0; i<holes.size(); i++) {cout<<holes[i][0]<<" "<<holes[i][1]<<endl;}
    int nds=holes.size();
    
    
    int gs=0;
//    vector<vector<double> > holes(nds, vector<double>(2,0));
    
//    for (int i=0; i<nds; i++) {
//        holes[i][0]=0.1*i+0.1;
//        holes[i][1]=0.1*i+0.1;
//    }
    
    /*
     inialization holes here.
     */
    
    vector<vector<double> > holes2(nds, vector<double>(2,0));
    vector<vector<double> > holes3(nds, vector<double>(2,0));
    
    int supermod(int k, int n);
    for(int i=0;i<nds;i++){
        holes2[supermod(i-1,nds)]=holes[i];
        holes3[supermod(i+1,nds)]=holes[i];
    }
        
    ofstream bout("berry_laughlin");
    LATTICE ll (Ne, invNu, testing, type, seed, gs);
    LATTICE ll2(Ne, invNu, testing, type, seed, gs);
    LATTICE ll3(Ne, invNu, testing, type, seed, gs);
    
    double energy=0.;
    complex<double> berry2=0., berry3=0.;
    for(int b=0;b<nds;b++){
        ll.set_hole(holes[b]);
        ll2.set_hole(holes2[b]);
        ll3.set_hole(holes3[b]);
        
        double phasemod(complex<double> in);
        ll.reset();
        ll.step(nWarmup);
        complex<double> btemp;
        for(int i=0;i<nMeas;i++){
            ll.step(nSteps);
            energy+=ll.coulomb_energy();
            berry2+=ll2.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
//            cout<<"berry2 advance abs= "<<abs(ll2.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs()))<<endl;
//            btemp=ll2.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
//            cout<<"output:"<<ll2.get_wf(ll.get_locs())<<" "<<ll.get_wf(ll.get_locs())<<" "<<abs(btemp)<<" "<<arg(btemp)<<endl;
            berry3+=ll3.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
        }
        
//        cout<<"ll.get_wf="<<ll.get_wf(ll.get_locs())<<endl;
//        cout<<"ll2.get_wf="<<ll2.get_wf(ll.get_locs())<<endl;
//        cout<<"ll3.get_wf="<<ll3.get_wf(ll.get_locs())<<endl;
//        cout<<"ll2.Ne="<<ll2.Ne<<endl;
//        cout<<"print out holes2"<<endl;
//        cout<<ll2.hole[0]<<" "<<ll2.hole[1]<<endl;
//        cout<<"finish printing"<<endl;
//        cout<<"print ll2. holes_set = "<<ll2.holes_set<<endl;
        
        bout<<holes[b][0]<<" "<<holes[b][1]<<" "<<holes2[b][0]<<" "<<holes2[b][1]<<"   "<<abs(berry2)/(1.*nMeas)<<"   "<<phasemod(berry2)<<"   "<<energy/(1.*nMeas*ll.Ne)<<endl;
        bout<<holes[b][0]<<" "<<holes[b][1]<<" "<<holes3[b][0]<<" "<<holes3[b][1]<<"   "<<abs(berry3)/(1.*nMeas)<<"   "<<phasemod(berry3)<<"   "<<energy/(1.*nMeas*ll.Ne)<<endl;
        //        tot_berry_phase2+=phasemod(berry2);
        //        tot_berry_phase3+=phasemod(berry3);
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