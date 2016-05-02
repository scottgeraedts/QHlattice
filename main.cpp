#include "lattice.h"
void single_run(LATTICE &ll, int nWarmup, int nMeas, int nSteps, int nBins){
	ofstream outfile("out"),eout("energy");
    double ave_E=0.;

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

double phasemod(complex<double> in){ 
	double out=arg(in);
	if(out<0) return out+2*M_PI;
	else if (out>2*M_PI) return out-2*M_PI;
	else return out;
}
void berry_phase(){
	int NPhi,invNu,nWarmup,nMeas,nSteps,nBins,seed;
	bool testing;
	string type;
	ifstream infile("params");
	infile>>NPhi>>invNu; 
	infile>>nWarmup>>nMeas>>nSteps>>nBins;
	infile>>seed;
	infile>>testing;
	infile>>type;

	LATTICE ll(NPhi,invNu, testing, type, seed);
	LATTICE ll2(NPhi,invNu, testing, type, seed);
	LATTICE ll3(NPhi,invNu, testing, type, seed);
	
	//these are the coordinates of the "extra" electrons that we are taking around the Fermi surface
	//for each MC run, we can compute two Berry phases, the results of moving electrons either clockwise or counterclockwise around the fermi surface
	vector<vector<int> > dcenter, d2, d3; 
	
	//here I'm initializing the ds for a system with 22 electrons
	if(NPhi!=46){
		cout<<"the berry phase calculator only works for 44 flux quanta right now!"<<endl;
		exit(0);
	}
	int nds=16;
	dcenter=vector <vector<int> >(nds,vector<int>(2));
	d2=vector <vector<int> >(nds,vector<int>(2));
	d3=vector <vector<int> >(nds,vector<int>(2));
	dcenter[0]={3,0};
	dcenter[1]={3,1};
	dcenter[2]={2,2};
	dcenter[3]={1,3};
	dcenter[4]={0,3};
	dcenter[5]={-1,3};
	dcenter[6]={-2,2};
	dcenter[7]={-3,1};
	dcenter[8]={-3,0};
	dcenter[9]={-3,-1};
	dcenter[10]={-2,-2};
	dcenter[11]={-1,-3};
	dcenter[12]={0,-3};
	dcenter[13]={1,-3};
	dcenter[14]={2,-2};
	dcenter[15]={3,-1};
	for(int i=0;i<nds;i++){
		d2[supermod(i-1,nds)]=dcenter[i];
		d3[supermod(i+1,nds)]=dcenter[i];
	}
	
	double energy;
	complex<double> berry2,berry3;
	vector <vector <int> > tempds;
	vector<int> dpoint(2);
	ofstream bout("berry");
	for(int b=0;b<nds;b++){
		//stuff for initializing with a custom set of ds (for the Berry phase calculation)
		double center_frac[2]={0.,0.}; //might need to change this for different N, but it works for N=22
		ll.make_fermi_surface(center_frac, ll.Ne-2);
		tempds=ll.get_ds();
		dpoint[0]=dcenter[b][0]; dpoint[1]=dcenter[b][1];
		tempds.push_back(dpoint);
		dpoint[0]=-dcenter[b][0]; dpoint[1]=-dcenter[b][1];
		tempds.push_back(dpoint);
		ll.set_ds(tempds);
	
		ll2.make_fermi_surface(center_frac,ll.Ne-2);
		tempds=ll2.get_ds();
		dpoint[0]=d2[b][0]; dpoint[1]=d2[b][1];
		tempds.push_back(dpoint);
		dpoint[0]=-d2[b][0]; dpoint[1]=-d2[b][1];
		tempds.push_back(dpoint);
		ll2.set_ds(tempds);

		ll3.make_fermi_surface(center_frac,ll.Ne-2);
		tempds=ll3.get_ds();
		dpoint[0]=d3[b][0]; dpoint[1]=d3[b][1];
		tempds.push_back(dpoint);
		dpoint[0]=-d3[b][0]; dpoint[1]=-d3[b][1];
		tempds.push_back(dpoint);
		ll3.set_ds(tempds);
		
		energy=0;
		berry2=0; berry3=0;
		ll.step(nWarmup);
		for(int i=0;i<nMeas;i++){
			ll.step(nSteps);
			energy+=ll.coulomb_energy();
			berry2+=ll2.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
			berry3+=ll3.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
		}
		bout<<dcenter[b][0]<<" "<<dcenter[b][1]<<" "<<d2[b][0]<<" "<<d2[b][1]<<" "<<abs(berry2)/(1.*nMeas)<<" "<<phasemod(berry2)<<" "<<energy/(1.*nMeas*ll.Ne)<<endl;
		bout<<dcenter[b][0]<<" "<<dcenter[b][1]<<" "<<d3[b][0]<<" "<<d3[b][1]<<" "<<abs(berry3)/(1.*nMeas)<<" "<<phasemod(berry3)<<" "<<energy/(1.*nMeas*ll.Ne)<<endl;
	}
	
	
		
}
int main(){
	berry_phase();
//	int NPhi,invNu,nWarmup,nMeas,nSteps,nBins,seed;
//	bool testing;
//	string type;
//	ifstream infile("params");
//	infile>>NPhi>>invNu; 
//	infile>>nWarmup>>nMeas>>nSteps>>nBins;
//	infile>>seed;
//	infile>>testing;
//	infile>>type;
//	//initialize MC object
//	
//    ofstream out_co_dbar("co_dbar");
//    
//    double Ne=1.*NPhi/(1.*invNu);
//    double dbar_parameter[2] = {0., 1.};
//    dbar_parameter[0]=1; dbar_parameter[1]=0;
//	LATTICE ll(NPhi,invNu, testing, type, seed);

////    ofstream eout_dbar("energy_dbar");
//    
//    
////    void coul_energy_dbar(LATTICE& edbar, double&, int nWarmup, int nMeas, int nSteps, int nBins, double* dbar_parameter );
//    
////    for (int i=0; i<5; i++) {
////        for (int j=0; j<5; j++) {
//////            cout<<"\ndbar_parameter = ("<<i<<", "<<j<<")/NPhi"<<endl;
////			dbar_parameter[0]=0.2*i;
////			dbar_parameter[1]=0.2*j;
////            coul_energy_dbar(ll,ave_E,nWarmup,nMeas,nSteps,nBins,dbar_parameter);
//////            cout<<"coulomb energy = "<<ave_E<<endl;
////            out_co_dbar<<i<<"   "<<j<<"   "<<ave_E<<endl;
////        }
////    }
//    
//	single_run(ll,nWarmup,nMeas,nSteps,nBins);    
}

void coul_energy_dbar(LATTICE& edbar, double& ave_E, int nWarmup, int nMeas, int nSteps, int nBins, double* dbar_parameter ){
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
//        cout<<E/(1.*nMeas*edbar.Ne)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))/(1.*edbar.Ne)<<endl;
        sumE+=E/(1.*nMeas*edbar.Ne);
    }
    ave_E=sumE/(1.*nBins);
//    outfile<<E/(1.*nMeas*ll.Ne)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))/(1.*ll.Ne)<<endl;
//    cout<<"acceptance rate: "<<(1.*ll.accepts)/(1.*ll.tries)<<endl;
}
