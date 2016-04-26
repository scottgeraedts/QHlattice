#include "lattice.h"
int main(){
	int NPhi,invNu,nWarmup,nSteps,nMeas,nBins,seed;
	ifstream infile("params");
	infile>>NPhi>>invNu;
	infile>>nWarmup>>nMeas>>nSteps>>nBins;
	infile>>seed;
	//initialize MC object
	
	ofstream outfile("out"),eout("energy");
	LATTICE ll(NPhi,invNu, seed);
	//take steps
//	ll.reset();
//	vector< vector<int> > tlocs=ll.get_locs();
//	cout<<ll.get_weight(tlocs)<<endl;
//	for(int i=0;i<ll.Ne;i++) tlocs[i][0]++;
//	cout<<ll.get_weight(tlocs)<<endl;
//	for(int i=0;i<ll.Ne;i++) tlocs[i][0]++;
//	cout<<ll.get_weight(tlocs)<<endl;
//	cout<<endl;
//	tlocs=ll.get_locs();
//	cout<<ll.get_weight(tlocs)<<endl;
//	for(int i=0;i<ll.Ne;i++) tlocs[i][1]++;
//	cout<<ll.get_weight(tlocs)<<endl;
//	for(int i=0;i<ll.Ne;i++) tlocs[i][1]++;   
//	cout<<ll.get_weight(tlocs)<<endl;
	
	for(int s=0;s<nBins;s++){
    
		ll.reset();		
		ll.step(nWarmup);
		double E=0,E2=0,P=0,P2=0,three=0;
		double e,p;
		deque<double> e_tracker, p_tracker;
		int Ntrack=10;
		vector<double> autocorr_e(Ntrack,0), autocorr_p(Ntrack,0);
		for(int i=0;i<nMeas;i++){
			cout<<ll.get_weight(ll.get_locs())<<" "<<norm(ll.get_wf(ll.get_locs(),ll.get_ds()))*exp(-ll.Ne*ll.Ne)<<endl;
			ll.step(nSteps);
			e=ll.coulomb_energy();
			E+=e;
			E2+=e*e;
			p=log(ll.running_weight);
			P+=p;
			P2+=p*p;
			three+=ll.threebody();
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
			
			ll.update_structure_factors();
		}
		outfile<<E/(1.*nMeas*ll.Ne)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))/(1.*ll.Ne)<<" "<<three/(1.*nMeas)<<endl;
		cout<<"acceptance rate: "<<(1.*ll.accepts)/(1.*ll.tries)<<endl;
	//	cout<<"almost done"<<endl;
	
		ofstream auto_out("auto");
		for(int j=0;j<Ntrack;j++){
			auto_out<<j+1<<" ";
			auto_out<<autocorr_e[j]/(1.*(nMeas-Ntrack))<<" "<<pow(E/(1.*nMeas),2)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))<<" ";
			auto_out<<autocorr_p[j]/(1.*(nMeas-Ntrack))<<" "<<pow(P/(1.*nMeas),2)<<" "<<(P2/(1.*nMeas)-pow(P/(1.*nMeas),2))<<" ";
			auto_out<<endl;
		}
	}
	cout<<"start print"<<endl;
	ll.print_structure_factors(nMeas*nBins);
	eout.close();
	outfile.close();
}
