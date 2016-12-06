#include "lattice.h"

LATTICE_WRAPPER::LATTICE_WRAPPER(int Ne_t, vector<vector<wf_info>> &wfs_t, int seed, bool testing_t){

	Ne=Ne_t;
	wfs=wfs_t;
	ran.seed(seed);
	zs=LATTICE::hot_start(Ne,Ne, ran);
	testing=testing_t;
	
	//set up running weights
	oldweight=vector<vector<complex<double>>>(wfs.size());
	for(int i=0; i<(signed)wfs.size(); i++){
		oldweight[i]=vector<complex<double>>(wfs[i].size());
		for(int j=0; j<(signed)wfs[i].size(); j++){
			oldweight[i][j]=get_wf(i,j,zs);
		}
	}
}

int LATTICE_WRAPPER::step(int nSteps){
	int electron;
	bool print=false;
	vector<int> oldz, newz, tempz(2);
	vector< vector<int> >::iterator itf;
	vector< vector<int> > tempzs;	
	complex<double> prod,prob,oldprob,oldprod,temp;	
	double normprob;
	vector<vector<complex<double>>> newweight=oldweight;

	for(int i=0;i<nSteps*Ne;i++){
		electron=ran.randInt(Ne-1);
		newz=LATTICE::random_move(zs[electron], Ne, ran);
//		itf=find(zs.begin(),zs.end(),newz);
//		if(itf!=zs.end()) continue;

		prob=0;	oldprob=0;
		for(auto it=wfs.begin();it!=wfs.end(); ++it){
			prod=1; oldprod=1;
			for(auto it2=it->begin(); it2!=it->end(); ++it2){
				if(electron>=it2->start and electron<it2->end){
					if(norm(oldweight[it-wfs.begin()][it2-it->begin()]) > 1e-14){
						tempz[0]=newz[0]*it2->sign;
						tempz[1]=newz[1]*it2->sign;
						temp=it2->wf.update_weight( it2->make_zs(zs), electron-it2->start, tempz)*oldweight[it-wfs.begin()][it2-it->begin()];
					}else{
						tempzs=zs;
						tempzs[electron]=newz;
						temp=get_wf(it-wfs.begin(), it2-it->begin(), tempzs);
					}
					if(print) cout<<"weight for "<<it-wfs.begin()<<" "<<it2-it->begin()<<" "<<temp<<endl;
				}else{
					temp=oldweight[it-wfs.begin()][it2-it->begin()];
				}
				
				newweight[it-wfs.begin()][it2-it->begin()]=temp;
				if( it2->denom){
					prod/=temp;
					oldprod/=oldweight[it-wfs.begin()][it2-it->begin()];
				}else{
					prod*=temp;
					oldprod*=oldweight[it-wfs.begin()][it2-it->begin()];
				}
			}
			prob+=prod;
			oldprob+=oldprod;
		}
		normprob=norm(prob/oldprob);	

		//cout<<electron<<" "<<newz[0]<<" "<<newz[1]<<" "<<normprob<<endl;

		bool update=false;
		if(normprob>=1) update=true;
		else if(ran.rand()<normprob) update=true;
	
		if(update){
			zs[electron]=newz;
			running_weight+=log(normprob);

			if(testing) cout<<"accepted: "<<running_weight<<" "<<log(norm(get_wf()))<<endl;
			for(auto it=wfs.begin();it!=wfs.end(); ++it){
				for(auto it2=it->begin(); it2!=it->end(); ++it2){			 
					if(electron>=it2->start and electron<it2->end) it2->wf.update();
					oldweight[it-wfs.begin()][it2-it->begin()]=newweight[it-wfs.begin()][it2-it->begin()];
				}
			}
			
			accepts++;
		}
		tries++;
			
	}
	return 1;
}

int LATTICE_WRAPPER::step_fromwf(int nSteps){
	int electron;
	vector<int> oldz, newz, tempz(2);
	vector< vector<int> >::iterator itf;
	double oldweight,newweight;	

	for(int i=0;i<nSteps*Ne;i++){
		electron=ran.randInt(Ne-1);
		newz=LATTICE::random_move(zs[electron], Ne, ran);
		oldz=zs[electron];
//		itf=find(zs.begin(),zs.end(),newz);
//		if(itf!=zs.end()) continue;

		oldweight=norm(get_wf());
		zs[electron]=newz;
		newweight=norm(get_wf());
					
		bool update=false;
		
		//cout<<electron<<" "<<newz[0]<<" "<<newz[1]<<" "<<newweight/oldweight<<endl;
		if(newweight/oldweight>=1) update=true;
		else if(ran.rand()<newweight/oldweight) update=true;
	
		if(update){
			running_weight+=log((newweight))-log((oldweight));

			if(testing) cout<<"accepted: "<<running_weight<<" "<<log(norm(get_wf()))<<endl;
			
			accepts++;
		}else{
			zs[electron]=oldz;
		}
		tries++;
	}
	return 1;	
}
void LATTICE_WRAPPER::reset(){
	for(auto it=wfs.begin();it!=wfs.end(); ++it)
		for(auto it2=it->begin(); it2!=it->end(); ++it2)
			it2->wf.reset( it2->make_zs(zs) );
	running_weight=log(norm(get_wf() ) );
	
}	

complex<double> LATTICE_WRAPPER::get_wf(){
	return get_wf(zs);
}
complex<double> LATTICE_WRAPPER::get_wf(const vector<vector<int>> &tempzs){
	complex<double> out=0,temp, prod, bigprod;
	for(auto it=wfs.begin();it!=wfs.end(); ++it){
		prod=1;
		for(auto it2=it->begin(); it2!=it->end(); ++it2){
			temp=it2->wf.get_wf( it2->make_zs(tempzs) );
			if( it2->conj ) temp=conj(temp);
			if( it2->denom ) prod/=temp;
			else prod*=temp;
		}
		bigprod*=prod;
		out+=abs(prod);
	} 
	return 0.01*out+bigprod;
}	

complex<double> LATTICE_WRAPPER::get_wf(int i, int j){
	return get_wf(i,j,zs);
}
complex<double> LATTICE_WRAPPER::get_wf(int i, int j, const vector<vector<int>> &tempzs){
	return wfs[i][j].wf.get_wf(wfs[i][j].make_zs(tempzs));
}

vector< vector<int> > LATTICE_WRAPPER::get_zs(){ return zs; }
void LATTICE_WRAPPER::acceptance_rate(){
	cout<<"acceptance rate: "<<(1.*accepts)/(1.*tries)<<endl;
}

wf_info::wf_info(){ }
wf_info::wf_info(bool conj_t, bool denom_t, int start_t, int end_t, int sign_t):conj(conj_t),denom(denom_t),start(start_t),end(end_t),sign(sign_t){

}
vector< vector<int> > wf_info::make_zs(const vector< vector<int> > &inzs){
	vector< vector<int> > outzs(end-start, vector<int>(2));
	copy(inzs.begin()+start, inzs.begin()+end, outzs.begin());
	for(auto it=outzs.begin(); it!=outzs.end(); ++it){
		(*it)[0]*=sign;
		(*it)[1]*=sign;
	}
	return outzs;
}
