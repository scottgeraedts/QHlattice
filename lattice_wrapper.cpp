#include "lattice.h"

LATTICE_WRAPPER::LATTICE_WRAPPER(int Ne_t, vector<wf_info> &wfs_t, int seed, bool testing_t){

	Ne=Ne_t;
	wfs=wfs_t;
	ran.seed(seed);
	zs=LATTICE::hot_start(Ne,Ne, ran);
	testing=testing_t;
}

int LATTICE_WRAPPER::step(int nSteps){
	int electron;
	vector<int> oldz, newz, tempz(2);
	vector< vector<int> >::iterator itf;
	vector< vector<int> > tempzs;	
	double prob=0,temp;	
	int count=0;

	for(int i=0;i<nSteps*Ne;i++){
		electron=ran.randInt(Ne-1);
		newz=LATTICE::random_move(zs[electron], Ne, ran);
		itf=find(zs.begin(),zs.end(),newz);
		if(itf!=zs.end()) continue;

		prob=0;			
		for(auto it=wfs.begin();it!=wfs.end(); ++it){
			if(electron>=(*it).start and electron<(*it).end){
				tempz[0]=newz[0]*(*it).sign;
				tempz[1]=newz[1]*(*it).sign;
				temp=(*it).wf.update_weight( (*it).make_zs(zs), electron-(*it).start, tempz);
				if( (*it).denom) prob-=temp;
				else prob+=temp;
			}
		}
			
		bool update=false;
		if(prob>0) update=true;
		else if(ran.rand()<exp(prob)) update=true;
	
		if(update){
			zs[electron]=newz;
			running_weight+=prob;

			if(testing) cout<<"accepted: "<<running_weight<<" "<<log(norm(get_wf()))<<endl;
			for(auto it=wfs.begin();it!=wfs.end(); ++it) 
				if(electron>=(*it).start and electron<(*it).end) (*it).wf.update();
			
			count++;
		}
			
	}
	return count;
}

void LATTICE_WRAPPER::reset(){
	for(auto it=wfs.begin();it!=wfs.end(); ++it)
		(*it).wf.reset( (*it).make_zs(zs) );
	running_weight=log(norm(get_wf() ) );
	
}	

complex<double> LATTICE_WRAPPER::get_wf(){
	complex<double> out=1,temp;
	for(auto it=wfs.begin();it!=wfs.end(); ++it){
		temp=(*it).wf.get_wf( (*it).make_zs(zs) );
		if( (*it).conj ) temp=conj(temp);
		if( (*it).denom ) out/=temp;
		else out*=temp;
	} 
	return out;
}	


vector< vector<int> > LATTICE_WRAPPER::get_zs(){ return zs; }

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
