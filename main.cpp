#include "lattice.h"
#include "berry_phase.h"

int main(){
	void single_run();
	single_run();
//    berry_phase bp(16);
//    bp.two_full_braiding();
    
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

void single_run(){
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

	LATTICE ll(NPhi,invNu, testing, type, seed);
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
