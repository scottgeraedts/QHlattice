#include "lattice.h"

int supermod(int x, int NPhi){	return (x%NPhi+NPhi)%NPhi; }
	
LATTICE::LATTICE(int NPhi_t, int invNu_t, int seed):NPhi(NPhi_t),invNu(invNu_t){
	//various parameters from input file
	L1=sqrt(2*M_PI*NPhi)/sqrt(2.);//these are only used for calls to Duncan's functions, if you use them in other places there will be problems due to 
	L2=complex<double> (0,real(L1));//the different definitions of magnetic length
	ran.seed(seed);
	if(NPhi%invNu) cout<<"NPhi not divisible by filling!"<<endl;
	Ne=NPhi/invNu;
	fermions=true;
	testing=false;
	type="laughlin";

//	cout<<NPhi<<" "<<invNu<<" "<<Ne<<" "<<L1<<" "<<L2<<endl;
	one=1; zero=0; //useful for fortran calls
	locs=vector< vector<int> >(Ne, vector<int>(2,0));//initalize locations of all electrons
	for(int i=0;i<Ne;i++){
		locs[i][0]=i;
		locs[i][1]=i;
	}

	ws=vector< vector<double> > (invNu, vector<double>(2,0) );
	for( int i=0;i<invNu;i++) ws[i][0]=( (i+0.5)/(1.*invNu)-0.5);
	
	//********calls to duncan's functions
	set_l_(&NPhi, &L1, &L2);
	setup_z_function_table_();
	sl2z=new int[4];
	sl2z[0]=1; sl2z[1]=0; sl2z[2]=0; sl2z[3]=1;
	setup_laughlin_state_(&Ne,&invNu,sl2z,&zero);
	running_weight=get_weight();
//	cout<<"starting weight "<<running_weight<<endl;

	//*****some counters
	tries=0; accepts=0;
	setup_coulomb();
	omega=vector <complex<double> >(NPhi);
	for(int i=0;i<NPhi;i++) omega[i]=polar(1.,2*M_PI*i/(1.*NPhi));
	sq=vector<vector<complex<double> > > (NPhi, vector<complex<double> >(NPhi,0));
	sq2=vector<vector<double> > (NPhi, vector<double>(NPhi,0));
}
void LATTICE::step(int Nsteps){
	
	for(int i=0;i<Nsteps*Ne;i++){
		tries++;
		accepts+=simple_update();
	}
	
}
int LATTICE::simple_update(){

	//find a new position for a randomly chosen electron
	int electron=ran.randInt(Ne-1);
	vector<int> newloc=random_move(locs[electron]);
	vector< vector<int> >::iterator it=find(locs.begin(),locs.end(),newloc);
	if(it!=locs.end()) return 0;

//	cout<<"start"<<endl;
	//figure out the probability difference from the vandermode part
	double prob=1;
	complex<double> temp;
	int xi,yi;
	for(int i=0;i<Ne;i++){
		//divide off old part
		if(i==electron) continue;
		xi=(locs[i][0]-locs[electron][0]);
		yi=(locs[i][1]-locs[electron][1]);
		if(i>electron){
			xi=-xi; yi=-yi;
		}		
		temp=lattice_z_(&NPhi,&xi,&yi,&L1,&L2,&one);
//		z_function_(&x,&y,&L1,&L2,&one,&NPhi,&temp);
		prob/=norm( pow(temp,invNu) );

		//multiply new part
		xi=(locs[i][0]-newloc[0]);
		yi=(locs[i][1]-newloc[1]);
		if(i>electron){
			xi=-xi; yi=-yi;
		}		
		temp=lattice_z_(&NPhi,&xi,&yi,&L1,&L2,&one);
//		z_function_(&x,&y,&L1,&L2,&one,&NPhi,&temp);
		prob*=norm( pow(temp,invNu) );
	}
//	cout<<"vandermond "<<prob<<endl;

	//figure out the probability difference from the COM part
	int oldCOM[2], newCOM[2];
	sum_locs(oldCOM);
	newCOM[0]=oldCOM[0]-locs[electron][0]+newloc[0];
	newCOM[1]=oldCOM[1]-locs[electron][1]+newloc[1];

	get_laughlin_cm_(oldCOM,&temp);
	prob/=norm(temp);	
	get_laughlin_cm_(newCOM,&temp);
	prob*=norm(temp);	
//	double x,y;
//	for( int i=0;i<invNu;i++){
//		x=oldCOM[0]/(1.*NPhi)-ws[i][0];
//		y=oldCOM[1]/(1.*NPhi)-ws[i][1];
//		z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
//		
//		prob/=norm(temp);

//		x=newCOM[0]/(1.*NPhi)-ws[i][0];
//		y=newCOM[1]/(1.*NPhi)-ws[i][1];
//		z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
//		prob*=norm(temp);
//	}

	//update or not
	bool update=false;
	if(prob>1) update=true;
	else if( ran.rand()<prob) update=true;
	
	if(update){
		locs[electron]=newloc;
		running_weight*=prob;
		if(testing) cout<<running_weight<<" "<<get_weight()<<endl;
//		for(int i=0;i<Ne;i++) cout<<locs[i][0]<<" "<<locs[i][1]<<endl;
		return 1;
	}
	else return 0;
}
//chooses a new site, given an old site
vector<int> LATTICE::random_move( const vector<int> &in){
	vector<int>newloc(2);

	//the stupidest way to do this
//	newloc[0]=ran.randInt(NPhi-1);
//	newloc[1]=ran.randInt(NPhi-1);

	double hoplength=2;
	int n=pow(2*hoplength+1,2)-1;
	vector<int> newx(n),newy(n);
	vector<double> newprob(n);
	int xi,yi;
	for(int i=0;i<=2*hoplength;i++){
		xi=i-hoplength+in[0]; 
		for(int j=0;j<=2*hoplength;j++){
			yi=j-hoplength+in[1];
			if(xi==in[0] && yi==in[1]) continue;
			newx[i]=supermod(xi,NPhi);
			newy[i]=supermod(yi,NPhi);
//			newprob.push_back( exp(-0.5*( pow(in[0]-i,2)+pow(in[1]-j,2) )/pow(hoplength,2) ) );
			newprob[i]=1.;
		}
	}
	double r=ran.rand();
	r*=accumulate(newprob.begin(),newprob.end(),0.);	
//	cout<<accumulate(newprob.begin(),newprob.end(),0.)<<endl;	
	
//	for(int i=0;i<(signed) newx.size();i++) newprob.push_back(1/(1.*newx.size()));
//	for(int i=0;i<(signed) newx.size();i++) cout<<newx[i]<<" "<<newy[i]<<" "<<in[0]<<" "<<in[1]<<" "<<newprob[i]<<endl;
	double runningprob=0;
	for(int i=0;i<(signed) newx.size();i++){
		runningprob+=newprob[i];
//		cout<<r<<" "<<runningprob<<" "<<newx[i]<<" "<<newy[i]<<" "<<in[0]<<" "<<in[1]<<endl;
		if(r<runningprob){
			newloc[0]=newx[i];
			newloc[1]=newy[i];
			return newloc;
		}
	}
	cout<<"never found a move!"<<endl;
	return newloc;
}	

int LATTICE::p(int site){
	if(site==NPhi-1) return 0;
	else return site+1;
}
int LATTICE::m(int site){
	if(site==0) return NPhi-1;
	else return site-1;
}
double LATTICE::get_weight(){
	double out=1,x,y;
	complex<double> temp;
	//vandermonde piece
	for( int i=0;i<Ne;i++){
		for( int j=i+1;j<Ne;j++){
			x=(locs[i][0]-locs[j][0])/(1.*NPhi);
			y=(locs[i][1]-locs[j][1])/(1.*NPhi);
			z_function_(&x,&y,&L1,&L2,&one,&NPhi,&temp);
			if(type=="laughlin") out*=norm( pow(temp,invNu) );
		}
	}
	int COM[2];
	sum_locs(COM);
	for( int i=0;i<invNu;i++){
		x=COM[0]/(1.*NPhi)-ws[i][0];
		y=COM[1]/(1.*NPhi)-ws[i][1];
		z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
		out*=norm(temp);
	}
	return out*exp(-Ne*Ne);
} 
void LATTICE::sum_locs(int out[]){
	out[0]=0; out[1]=0;
	for( int i=0;i<Ne;i++){
		out[0]+=locs[i][0];
		out[1]+=locs[i][1];
	}
	out[0]=out[0]%NPhi;
	out[1]=out[1]%NPhi;
}

///***********MEASUREMENT FUNCTIONS *******///////////////
//right now I'm using Duncan's function for this
//double LATTICE::two_body_coulomb(int dx, int dy){
//	double out=0;
//	if(dx>NPhi/2) dx-=NPhi;
//	if(dx<-NPhi/2) dx+=NPhi;	
//	if(dy>NPhi/2) dy-=NPhi;
//	if(dy<-NPhi/2) dy+=NPhi;	

//	if(dx==0 && dy==0) out=-	
//	for(int mx=0;mx>-1;mx++){
//		for(int my=0; my>-1;my++){

//		}	
//	}
//}
void LATTICE::setup_coulomb(){
	coulomb_table=vector<vector<double> >(NPhi,vector<double>(NPhi,0));
//	coulomb_table[0][0]=v_coulomb_(&NPhi,&NPhi,&NPhi,&L1,&L2);
	for(int m=0;m<NPhi;m++){
		for(int n=0;n<NPhi;n++){
			coulomb_table[m][n]=v_coulomb_(&NPhi,&m,&n,&L1,&L2);
		}
	}
}
double LATTICE::coulomb_energy(){
	double out=0.5*Ne*coulomb_table[0][0];
	int m,n;
	for(int i=0;i<Ne;i++){
		for(int j=0;j<i;j++){
			m=((locs[i][0]-locs[j][0])%NPhi+NPhi)%NPhi;
			n=((locs[i][1]-locs[j][1])%NPhi+NPhi)%NPhi;
//			cout<<m<<" "<<n<<endl;
			out+=coulomb_table[m][n];
		}
	}
	return out;
}
void LATTICE::update_structure_factors(){
	double kappa=2.*M_PI/(1.*NPhi);
	complex<double>temp;
	for(int qx=0;qx<NPhi;qx++){
		for(int qy=0; qy<NPhi; qy++){
			temp=0;
			for(int i=0;i<Ne;i++){
				temp+=omega[supermod(qx*locs[i][0]+qy*locs[i][1],NPhi)];
				//polar(1.,qx*kappa*locs[i][0]+qy*kappa*locs[i][1]);
			}
			sq[qx][qy]+=temp;
			sq2[qx][qy]+=norm(temp);
		}
	}
}
void LATTICE::print_structure_factors(int nMeas){
	ofstream sqout("sq");
	ofstream sqout2("sq2");
//	double kappa=2.*M_PI/(1.*NPhi);
	for(int qx=0;qx<NPhi;qx++){
		for(int qy=0; qy<NPhi; qy++){
			sqout2<<sq2[qx][qy]/(1.*nMeas)<<" ";
			sqout<<abs(sq[qx][qy]/(1.*nMeas))<<" ";
		}
		sqout<<endl;
		sqout2<<endl;
	}
	sqout.close();
	sqout2.close();
}
LATTICE::~LATTICE(){ delete []sl2z; }
