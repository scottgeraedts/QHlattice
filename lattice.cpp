#include "lattice.h"

int supermod(int k, int n){	return ((k %= n) < 0) ? k+n : k; }

LATTICE::LATTICE(){
	Ne=0;
}
LATTICE::LATTICE(int Ne_t, int invNu_t, bool testing_t=false, string type_t="CFL", int seed=0, int gs):Ne(Ne_t),invNu(invNu_t),type(type_t){
	//various parameters from input file
	testing=testing_t;
	NPhi=Ne*invNu;
	if(type=="laughlin-hole"){
		NPhi++;
	}
	L1=sqrt(2*M_PI*NPhi)/sqrt(2.);//these are only used for calls to Duncan's functions, if you use them in other places there will be problems due to 
	L2=complex<double> (0,real(L1));//the different definitions of magnetic length
	ran.seed(seed);
	fermions=true;
    /*
	//setting up jie's function, might never use this
    vector<complex<double> > zeros; for(int i=0; i<NPhi; i++) zeros.push_back(complex<double>(0,0));
	weiers=weierstrass(.5*L1, .5*L2, zeros);
     */
    
	one=1; zero=0; //useful for fortran calls

	//****initialize z's, w's, d's
	locs=vector< vector<int> >(Ne, vector<int>(2,0));//initalize locations of all electrons

	//setting the ws. Note that the sum of these is ALWAYS zero, adding things like composite fermion momenta or holes doesn't change this.
	//in the y direction these take the values gs*L/invNu, where gs in (0,invNu-1) is an integer which labels the ground state
	ws=vector< vector<double> > (invNu, vector<double>(2,0) );
	for( int i=0;i<invNu;i++){
        ws[i][0]=( (i+0.5)/(1.*invNu)-0.5);
        ws[i][1]=gs/(1.*invNu);
//        ws[i][1]=( (i+0.5)/(1.*invNu)-0.5);
//        ws[i][0]=gs/(1.*invNu);
	}

	double center_frac[2]={0.,0.};
	if(type=="CFL"){
		if(Ne%2==0){ center_frac[0]=0.5/(1.*Ne); center_frac[1]=0.5/(1.*Ne);}
		make_fermi_surface(center_frac, Ne);
				
		print_ds();
		dsum=vector<int>(2,0);
		for(int i=0;i<Ne;i++){
			dsum[0]+=ds[i][0]*invNu; dsum[1]+=ds[i][1]*invNu;
            //dsum=NPhi * Ne dbar, i.e. it is the point on the lattice of electrons where the TOTAL d lives
            //'ds' is defined on L/Ne lattice, 'dsum' in this way is defined on L/Nphi lattice.
		}
		change_dbar_parameter(dsum[0]/(1.*Ne),dsum[1]/(1.*Ne));

		//the average d should also be on a lattice point, so dividing dsum by Ne should yield an integer
//		if(dsum[0]%Ne || dsum[1]%Ne) cout<<"Warning! The average of the ds is not on a lattice point! "<<dsum[0]<<" "<<dsum[1]<<endl;
//		cout<<"dsum: "<<dsum[0]<<" "<<dsum[1]<<endl;
	}
	holes_set=false;
	
	//********calls to duncan's functions
	set_l_(&NPhi, &L1, &L2);
	setup_z_function_table_();
	int *sl2z=new int[4];
	sl2z[0]=1; sl2z[1]=0; sl2z[2]=0; sl2z[3]=1;
	if(type!="laughlin-hole") setup_laughlin_state_(&Ne,&invNu,sl2z,&zero);
//	cout<<"starting weight "<<running_weight<<endl;

	cout<<"starting coulomb setup"<<endl;
	setup_coulomb();
	cout<<"done coulomb setup"<<endl;
	omega=vector <complex<double> >(2*NPhi);
	for(int i=0;i<2*NPhi;i++) omega[i]=polar(1.,M_PI*i/(1.*NPhi)); //note that spacings of these is pi/N, not 2pi/N

	//*****some counters
	sq=vector<vector<complex<double> > > (NPhi, vector<complex<double> >(NPhi,0));
	sq2=vector<vector<double> > (NPhi, vector<double>(NPhi,0));
	sq3=vector <vector< vector< vector <complex<double> > > > >(NPhi, vector <vector <vector< complex<double> > > >(NPhi, vector <vector <complex<double> > >(NPhi, vector<complex<double> >(NPhi,0))));
	delete [] sl2z;
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
	double prob=0;
	complex<double> temp;

	//***************hole part
	double dx,dy; //TODO: this calls z_function every time, should speed that up
	if(type=="laughlin-hole"){
		dx=(locs[electron][0]/(1.*NPhi)-hole[0]);	
		dy=(locs[electron][1]/(1.*NPhi)-hole[1]);	
		z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
		prob-=log(norm(temp));
		dx=(newloc[0]/(1.*NPhi)-hole[0]);	
		dy=(newloc[1]/(1.*NPhi)-hole[1]);	
		z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
		prob+=log(norm(temp));
		
	}		
		
	//***************vandermode part
	int vandermonde_exponent=invNu;
	if(type=="CFL") vandermonde_exponent-=2;
	int xi,yi;
	
	if(vandermonde_exponent!=0){
		for(int i=0;i<Ne;i++){
			//divide off old part
			if(i==electron) continue;
			xi=(locs[i][0]-locs[electron][0]);
			yi=(locs[i][1]-locs[electron][1]);
			if(i>electron){
				xi=-xi; yi=-yi;
			}		
			temp=lattice_z_(&NPhi,&xi,&yi,&L1,&L2,&one);
			prob-=log(norm( pow(temp,vandermonde_exponent) ));

			//multiply new part
			xi=(locs[i][0]-newloc[0]);
			yi=(locs[i][1]-newloc[1]);
			if(i>electron){
				xi=-xi; yi=-yi;
			}		
			temp=lattice_z_(&NPhi,&xi,&yi,&L1,&L2,&one);
			prob+=log(norm( pow(temp,vandermonde_exponent) ));
		}
	}

	//***************COM PART
	//figure out the probability difference from the COM part


	int oldCOM[2], newCOM[2];
	sum_locs(oldCOM);// 'locs' is defined on L/Nphi lattice. So does COM.
	if(type=="CFL"){
		oldCOM[0]-=dsum[0]/invNu;
		oldCOM[1]-=dsum[1]/invNu;
	    //The reason for dividing invNu is: sum_{i=1}^{invNu}w_i = sum_{j=1}^{Ne}d_j, so dsum/invNu Actually means average of w_i. The w_i in this code sums to 0.
	}
	newCOM[0]=oldCOM[0]-locs[electron][0]+newloc[0];
	newCOM[1]=oldCOM[1]-locs[electron][1]+newloc[1];

	if(type=="laughlin-hole"){
		double dx,dy;
		for( int i=0;i<invNu;i++){
			dx=oldCOM[0]/(1.*NPhi)-ws[i][0]+hole[0]/(1.*invNu);
			dy=oldCOM[1]/(1.*NPhi)-ws[i][1]+hole[1]/(1.*invNu);
			z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
			prob-=log(norm(temp));
			dx=newCOM[0]/(1.*NPhi)-ws[i][0]+hole[0]/(1.*invNu);
			dy=newCOM[1]/(1.*NPhi)-ws[i][1]+hole[1]/(1.*invNu);
			z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
			prob+=log(norm(temp));
		}
	}else{

		get_laughlin_cm_(oldCOM,&temp);
		prob-=log(norm(temp));	
		get_laughlin_cm_(newCOM,&temp);
		prob+=log(norm(temp));	
	}
	///***********determinant part
	complex<double> newDeterminant;
	Eigen::MatrixXcd newMatrix=oldMatrix;
    if(type=="CFL"){
		
		make_CFL_det(newMatrix, newloc, electron, newDeterminant);
		prob+=log(norm(newDeterminant/oldDeterminant));
	}
	  
    //*******************update or not
	bool update=false;
	if(prob>0) update=true;
	else if( ran.rand()<exp(prob)) update=true;
	
	if(update){
		locs[electron]=newloc;
		running_weight+=prob;
//		cout<<prob<<endl;
		if(testing) cout<<running_weight<<" "<<get_weight(locs)<<" "<<log(norm(get_wf(locs)))<<endl;
//		cout<<"new:"<<endl<<newMatrix<<endl;
		oldDeterminant=newDeterminant;
        oldMatrix=newMatrix;
//        for(int i=0;i<Ne;i++) cout<<locs[i][0]<<" "<<locs[i][1]<<endl;
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

	int hoplength=2;
//	if(Ne<10) hoplength=2;
	int n=pow(2*hoplength+1,2)-1;
	vector<int> newx(n),newy(n);
	vector<double> newprob(n);
	int xi,yi,counter=0;
	for(int i=0;i<=2*hoplength;i++){
		xi=i-hoplength+in[0]; 
		for(int j=0;j<=2*hoplength;j++){
			yi=j-hoplength+in[1];
			if(xi==in[0] && yi==in[1]) continue;
			newx[counter]=supermod(xi,NPhi);
			newy[counter]=supermod(yi,NPhi);
//			newprob.push_back( exp(-0.5*( pow(in[0]-i,2)+pow(in[1]-j,2) )/pow(hoplength,2) ) );
			newprob[counter]=1.;
			counter++;
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
double LATTICE::get_weight(const vector< vector<int> > &zs){
	double out=0,x,y;
	complex<double> temp,temp2;
	//hole piece
	if(type=="laughlin-hole"){
		for(int i=0;i<Ne;i++){
			x=(locs[i][0]/(1.*NPhi)-hole[0]);	
			y=(locs[i][1]/(1.*NPhi)-hole[1]);	
			z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
			out+=log(norm(temp));
		}		
	}		
		

	//vandermonde piece
	int vandermonde_exponent=invNu;
	if(type=="CFL") vandermonde_exponent-=2;
	complex<double> z;
	for( int i=0;i<Ne;i++){
		for( int j=i+1;j<Ne;j++){
			x=(zs[i][0]-zs[j][0])/(1.*NPhi);
			y=(zs[i][1]-zs[j][1])/(1.*NPhi);
			z_function_(&x,&y,&L1,&L2,&one,&NPhi,&temp);
//			temp=jies_weierstrass(x,y);
			out+=log(norm( pow(temp,vandermonde_exponent) ));
		}
	}
	
	//com part
	int COM[2]={0,0};
	for( int i=0;i<Ne;i++){
		COM[0]+=zs[i][0];
		COM[1]+=zs[i][1];
	}

	if(type=="CFL"){
		COM[0]-=dsum[0]/invNu;//this is correct since dsum lives on the same lattice as the zs
		COM[1]-=dsum[1]/invNu;
	}
	for( int i=0;i<invNu;i++){
		x=COM[0]/(1.*NPhi)-ws[i][0];
		y=COM[1]/(1.*NPhi)-ws[i][1];
		if(type=="laughlin-hole"){
			x+=hole[0]/(1.*invNu);
			y+=hole[1]/(1.*invNu);
		}
		z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
		out+=log(norm(temp));
	}
	//determinant part
	double oldDivisor;
	if(type=="CFL"){
		complex<double> product;
		Eigen::MatrixXcd M(Ne,Ne);
		for(int i=0;i<Ne;i++){
			for(int j=0;j<Ne;j++){
				product=1;
				for(int k=0;k<Ne;k++){
					if(k==i) continue;
//                    x=det_helper(zs[i][0],zs[k][0],ds[j][0],dsum[0])/(1.*NPhi);
//                    y=det_helper(zs[i][1],zs[k][1],ds[j][1],dsum[1])/(1.*NPhi);
                    x=det_helper(zs[i][0],zs[k][0],ds[j][0],dbar_parameter[0])/(1.*NPhi);
                    y=det_helper(zs[i][1],zs[k][1],ds[j][1],dbar_parameter[1])/(1.*NPhi);
					z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
//					cout<<temp<<" "<<x<<" "<<y<<endl;
					product*=temp;
				}
				//this part is only valid on a square torus!
				M(i,j)=product*pow(in_determinant_rescaling,Ne-1)*polar(1., 2*M_PI*NPhi*(zs[i][1]*ds[j][0] - zs[i][0]*ds[j][1])/(2.*invNu*NPhi*Ne) );
			}
		}
		detSolver.compute(M);
		temp=detSolver.determinant(); 
		//ameliorate large determinants by shrinking before taking their norm
		oldDivisor=abs(real(temp))+abs(imag(temp));
		out+=log(norm(temp/oldDivisor))+2*log(oldDivisor);
//		temp2=temp*exp(out);
//		out=log(real(temp2*conj(temp)));
	}		
	return out;
} 
//given both a set of positions and a set of ds, computes the wavefunction (NOT the norm of the wavefunction)
//also different compared to get_weight: this function uses the precomputed lattice_z data
complex<double> LATTICE::get_wf(const vector< vector<int> > &zs){
	complex<double> out=1,temp;
	int ix,iy;

	double x,y;
	//hole part
	if(type=="laughlin-hole"){
		for(int i=0;i<Ne;i++){
			x=(zs[i][0]/(1.*NPhi)-hole[0]);	
			y=(zs[i][1]/(1.*NPhi)-hole[1]);	
			z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
			out*=temp;
		}		
	}		
    
	//vandermonde piece
	int vandermonde_exponent=invNu;
	if(type=="CFL") vandermonde_exponent-=2;
	if(vandermonde_exponent!=0){
		for( int i=0;i<Ne;i++){
			for( int j=i+1;j<Ne;j++){
				ix=(zs[i][0]-zs[j][0]);
				iy=(zs[i][1]-zs[j][1]);
				out*=pow(lattice_z_(&NPhi,&ix,&iy,&L1,&L2,&one),vandermonde_exponent);
			}
		}
	}
	
	//COM piece
	int COM[2]={0,0};
	for( int i=0;i<Ne;i++){
		COM[0]+=zs[i][0];
		COM[1]+=zs[i][1];
	}
	if(type=="laughlin-hole"){
		double dx,dy;
		for( int i=0;i<invNu;i++){
			dx=COM[0]/(1.*NPhi)-ws[i][0]+hole[0]/(1.*invNu);
			dy=COM[1]/(1.*NPhi)-ws[i][1]+hole[1]/(1.*invNu);
			z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
			out*=temp;
		}	
	}
    else{
		if(type=="CFL"){
			COM[0]-=dsum[0]/invNu;
			COM[1]-=dsum[1]/invNu;
		}
		get_laughlin_cm_(COM,&temp);
		out*=temp;
	}

	if(type=="CFL"){
		complex<double> product;
		vector<int> z(2);
		Eigen::MatrixXcd M(Ne,Ne);
		for(int i=0;i<Ne;i++){
			for(int j=0;j<Ne;j++){
				product=1;
				for(int k=0;k<Ne;k++){
					if(k==i) continue;
					det_helper(zs[i],zs[k],ds[j],z);
					temp=modded_lattice_z(z[0],z[1]);
					product*=temp;
				}
				M(i,j)=product*polar(pow(in_determinant_rescaling,Ne-1), 2*M_PI*NPhi*(zs[i][1]*ds[j][0] - zs[i][0]*ds[j][1])/(2.*invNu*NPhi*Ne) );
			}
		}
		detSolver.compute(M);
		out=out*detSolver.determinant();
	}		
	return conj(out);
} 

void LATTICE::sum_locs(int out[]){
	out[0]=0; out[1]=0;
	for( int i=0;i<Ne;i++){
		out[0]+=locs[i][0];
		out[1]+=locs[i][1];
	}
//	out[0]=out[0]%NPhi;
//	out[1]=out[1]%NPhi;
}

vector< vector<int> > LATTICE::get_locs(){ return locs; }


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
	complex<double> eL1=L1*sqrt(2.);
	complex<double> eL2=L2*sqrt(2.);
	for(int m=0;m<NPhi;m++){
		for(int n=0;n<NPhi;n++){
			coulomb_table[m][n]=v_coulomb_(&NPhi,&m,&n,&eL1,&eL2);
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
/*
make_fermi_surface(l1,l2,Ne,center_frac,ds);
 l1, l2 are primitive lattice. N is #e.
 N might differ from Ne because if you want to calculate Berry phases, you don't want to place all the electrons
 center_frac = (x0,y0). x0*l1+y0*l2 is the center position of fermi surface.
 ds contains d s, defined on L/Ne lattice.
*/
void LATTICE::make_fermi_surface(double* center_frac, int N){
	ds.clear();
    vector<vector<int> > d_list;
    double x0,y0; x0=center_frac[0]; y0=center_frac[1];
    //initial sub-lattice: L_{mn}/Ne where d s lives on.
    vector<int> d=vector<int>(2,0);
    d_list.push_back(d);
    for (int i=-Ne; i<Ne; i++) {
        for (int j=-Ne; j<Ne; j++) {
            if (i==0 && j==0) continue;
            vector<int> d;
            d.push_back(i); d.push_back(j);
            d_list.push_back(d);
            d.clear();
        }
    }
//uncomment this if the lattice should be shifted    
    for (int i=0; i<(signed)d_list.size(); i++) {
        d_list[i][0]+=(int)(x0*Ne); d_list[i][1]+=(int)(y0*Ne);
    }
    //fill electrons one by one, fill those close to center first.
    int x,y; double min,min_;
    vector< vector<int> >::iterator it;
    
    for (int k=0; k<N; k++) {
        it=d_list.begin();
        
        x=d_list[0][0]; y=d_list[0][1]; min = norm(L1*((double)x/(1.*Ne)-x0)+L2*((double)y/(1.*Ne)-y0));
        for (int i=1; i<(signed)d_list.size(); i++) {
            x=d_list[i][0]; y=d_list[i][1]; min_ = norm(L1*((double)x/(1.*Ne)-x0)+L2*((double)y/(1.*Ne)-y0));
            if (min_<min) {
                it=d_list.begin()+i;
                min=min_;
            }
        }
        vector<complex<double> > tmp;
//        ds.push_back(1.*(*it)[0]*l1/(1.*Ne)+1.*(*it)[1]*l2/(1.*Ne));
		ds.push_back(*it);
        it=d_list.erase(it);
    }
    
}
void LATTICE::print_ds(){
	ofstream dout("ds");
	for (unsigned int i = 0; i < ds.size(); i += 1)
		dout<<ds[i][0]<<" "<<ds[i][1]<<endl;
	dout.close();
}
double LATTICE::threebody(){
	int has_right,has_top,has_diag;
	int three_counter=0, two_counter=0;
	double out=0;
	for(int i=0;i<Ne;i++){
		has_right=0;
		has_top=0;
		has_diag=0;
		for(int j=0;j<Ne;j++){
			if(locs[j][0]==p(p(locs[i][0])) && locs[j][1]==locs[i][1]) has_right=1;
			if(locs[j][0]==locs[i][0] && locs[j][1]==p(p(locs[i][1]))) has_top=1;
			if(locs[j][0]==p(p(locs[i][0])) && locs[j][1]==m(m(locs[i][1])) ) has_diag=1;
		}
		if(has_right && has_top) three_counter++;
		two_counter+=has_right+has_top+has_diag;
	}
	out=(three_counter-0.5*(two_counter)+0.5*Ne/(1.*NPhi))/(1.*NPhi);
	return out;
}

void LATTICE::update_structure_factors(){
	vector< vector<complex<double> > >temp(NPhi,vector<complex<double> >(NPhi,0));
	for(int qx=0;qx<NPhi;qx++){
		for(int qy=0; qy<NPhi; qy++){
			for(int i=0;i<Ne;i++){
				temp[qx][qy]+=omega[supermod(2*(qx*locs[i][0]+qy*locs[i][1]),2*NPhi)];
				//polar(1.,qx*kappa*locs[i][0]+qy*kappa*locs[i][1]);
			}
			sq[qx][qy]+=temp[qx][qy];
			sq2[qx][qy]+=norm(temp[qx][qy]);
		}
	}
//	int qx3,qy3;
//	for(int qx1=0;qx1<NPhi;qx1++){
//		for(int qx2=0;qx2<NPhi;qx2++){
//			for(int qy1=0;qy1<NPhi;qy1++){
//				for(int qy2=0;qy2<NPhi;qy2++){
//					qx3=supermod(-qx1-qx2,NPhi);
//					qy3=supermod(-qy1-qy2,NPhi);
//					sq3[qx1][qx2][qy1][qy2]+=temp[qx1][qy1]*temp[qx2][qy2]*temp[qx3][qy3];
//				}
//			}
//		}
//	}

}
			
void LATTICE::print_structure_factors(int nMeas){
	ofstream sqout("sq"), sqout2("sq2"),sqout3("sq3");
	for(int qx=0;qx<NPhi;qx++){
		for(int qy=0; qy<NPhi; qy++){
			sqout2<<sq2[qx][qy]/(1.*nMeas)<<" ";
			sqout<<abs(sq[qx][qy]/(1.*nMeas))<<" ";
		}
		sqout<<endl;
		sqout2<<endl;
	}
	
//	for(int qx1=0;qx1<NPhi;qx1++){
//		for(int qy1=0;qy1<NPhi;qy1++){
//			for(int qx2=0;qx2<NPhi;qx2++){
//				for(int qy2=0;qy2<NPhi;qy2++)
//					sqout3<<real(sq3[qx1][qx2][qy1][qy2])/(1.*nMeas)<<" ";
//				sqout3<<endl;
//			}
//		}
//	}
	sqout3.close();
	sqout.close();
	sqout2.close();
}
void LATTICE::reset(){
	tries=0; accepts=0;
	hotter_start();
	//**** setting up the initial determinant matrix
	running_weight=0;
	int site=0;
	int initial_state_counter=0;
	bool repeat=true;
	
    while(repeat){
   		repeat=false;
		
		if(type=="CFL"){
			complex<double> temp,product;
			double x,y;
			oldMatrix=Eigen::MatrixXcd(Ne,Ne);
			for(int i=0;i<Ne;i++){
				for(int j=0;j<Ne;j++){
					product=1;
					for(int k=0;k<Ne;k++){
						if(k==i) continue;
                        x=det_helper(locs[i][0],locs[k][0],ds[j][0],dbar_parameter[0])/(1.*NPhi);
                        y=det_helper(locs[i][1],locs[k][1],ds[j][1],dbar_parameter[1])/(1.*NPhi);
						z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
						product*=temp;
					}
					oldMatrix(i,j)=product*polar(pow(in_determinant_rescaling,Ne-1), M_PI*(locs[i][1]*ds[j][0] - locs[i][0]*ds[j][1])/(1.*NPhi) );
				}
			}
			detSolver.compute(oldMatrix);
			oldDeterminant=detSolver.determinant();
//			cout<<in_determinant_rescaling<<" "<<oldDeterminant<<endl;
		}
		running_weight=get_weight(locs);
		cout<<running_weight<<" "<<get_wf(locs)<<endl;
		if(get_wf(locs)==complex<double>(0,0)){
			//for some sizes the configuration specified by cold_start has zero weight
			//if that happens fiddle around until you find a better configuration, thats why theres a while loop
			cout<<"warning! needed to fiddle with the initial configuration "<<running_weight<<endl;
			for(int i=0;i<Ne;i++) cout<<locs[i][0]<<" "<<locs[i][1]<<endl;
			locs[site][0]=locs[p(site)][0];
			site++;
			if(site==Ne) cout<<"couldn't easily find a good configuration!"<<endl;
			repeat=true;
		}
		initial_state_counter++;
		if(initial_state_counter>1000){
			cout<<"couldn't find a good starting configuration"<<endl;
			exit(0);
		}
	}
    check_sanity();
}
//checks a few different things to make sure that they make sense
void LATTICE::check_sanity(){
	if(type!="CFL" && type!="laughlin" && type!= "laughlin-hole"){
		cout<<"type not recognized: "<<type<<endl;
		exit(0);
	}
	if(type=="CFL" && shifted_ztable.size()==0){
		cout<<"CFL but haven't set up the d's"<<endl;
		exit(0);
	}
	if(type=="laughlin-hole" && holes_set==false){
		cout<<"laughlin-hole but haven't set the hole position"<<endl;
		exit(0);
	}
	
}
//changes the dbar (which ordinarily is the sum of d's) to whatever we want
void LATTICE::change_dbar_parameter(double dbarx, double dbary){
	if(type!="CFL"){
		cout<<"changing dbar, but not a CFL"<<endl;
		//exit(0);
	}
	dbar_parameter=vector<double>(2);
	dbar_parameter[0]=dbarx;
	dbar_parameter[1]=dbary;
	//cout<<"dbar "<<dbar_parameter[0]<<" "<<dbar_parameter[1]<<endl;

	complex<double> temp;
	double x,y;
	shifted_ztable=vector< vector< complex<double> > > (NPhi,vector< complex<double> >(NPhi,0));

	for(int ix=0;ix<NPhi;ix++){
		x=(ix+dbar_parameter[0])/(1.*NPhi);
		for(int iy=0;iy<NPhi;iy++){
			y=(iy+dbar_parameter[1])/(1.*NPhi);
			z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
			shifted_ztable[ix][iy]=temp;
		}
	}
}
void LATTICE::set_ds(vector< vector<int> > tds){
	ds=tds;
	dsum=vector<int>(2,0);
	for(int i=0;i<Ne;i++){
		dsum[0]+=ds[i][0]*invNu; dsum[1]+=ds[i][1]*invNu;
        //dsum=NPhi * Ne dbar, i.e. it is the point on the lattice of electrons where the TOTAL d lives
        //'ds' is defined on L/Ne lattice, 'dsum' in this way is defined on L/Nphi lattice.
	}
	//print_ds();
	change_dbar_parameter(dsum[0]/(1.*Ne),dsum[1]/(1.*Ne));	
	reset();
}
vector< vector<int> > LATTICE::get_ds(){ return ds; }
vector<double> LATTICE::get_dbar_parameter(){return dbar_parameter;}

void LATTICE::set_hole(vector<double> temphole){
	hole=temphole;
	holes_set=true;
}
vector<double> LATTICE::get_hole(){ return hole; }

inline void LATTICE::det_helper(const vector<int> &z1, const vector<int> &z2, const vector<int> &d, vector<int> &z){
	z[0]=z1[0]-z2[0]-d[0]*invNu;
	z[1]=z1[1]-z2[1]-d[1]*invNu;
}
inline double LATTICE::det_helper(int z1, int z2, int d, double dbarp){ return z1-z2-d*invNu+dbarp;}

complex<double> LATTICE::jies_weierstrass(double x, double y){
	complex<double> z(x,y);
	complex<double> out=weiers.wsigma(z)*exp(-0.5*pow(z,2)*weiers.Gbar/(1.))*exp(-0.5*z*conj(z)/(1.*NPhi));
	return out;
}
//puts all the electrons in a diagonal line, a pretty low weight starting position but works for most cases
void LATTICE::cold_start(){
	for(int i=0;i<Ne;i++){
		locs[i][0]=i*invNu;
		locs[i][1]=i*invNu;
	}
}
//try to figure out a starting configuration that will have a higher weight
void LATTICE::hotter_start(){
	bool repeat=true;
	for(int i=0;i<Ne;i++){
		while(repeat){
			locs[i][0]=ran.randInt(NPhi-1);
			locs[i][1]=ran.randInt(NPhi-1);
			repeat=false;
			for(int j=0;j<i;j++)
				if(locs[i][0]==locs[j][0] && locs[i][1]==locs[j][1]) repeat=true;
		}
		repeat=true;
	}
}

//call's duncan's lattice_z function, if the arguments x or y are outside of the range (0,NPhi) it shifts them into that range
//and multiplies by the appropriate phase
//only works for a square torus
//on a square torus, the phase is always +- 1? 
complex<double> LATTICE::modded_lattice_z(int x, int y){// why need this function? because want to use z_function table, which is given in first BZ.
	int modx=supermod(x,NPhi);
	int mody=supermod(y,NPhi);
//	complex<double> out=lattice_z_(&NPhi,&modx,&mody,&L1,&L2,&one);
	complex<double> out=shifted_ztable[modx][mody];
	int j=(modx-x)/NPhi, k=(mody-y)/NPhi;
//	out*=omega[supermod(-(mody+dbar_parameter[1])*j+(modx+dbar_parameter[0])*k,2*NPhi)];
	out*=polar( 1., (-(mody+dbar_parameter[1])*j+(modx+dbar_parameter[0])*k)*M_PI/(1.*NPhi));
//	out*=polar(1.,-M_PI/(1.*NPhi)*(mody*j-modx*k));
//	cout<<polar(1.,-M_PI/(1.*NPhi)*(y*j-x*k))<<endl;
	if(j%2 || k%2) return -out;
	else return out;
}

void LATTICE::make_CFL_det(Eigen::MatrixXcd& newMatrix, vector<int> newloc, int electron, complex<double>& new_det){
    complex<double> product;
    vector<int> z(2);
    complex<double> temp;
    
    for(int i=0;i<Ne;i++){
        for(int j=0;j<Ne;j++){
            if(i==electron){
                product=1.;
                for(int k=0;k<Ne;k++){
                    if(k==i) continue;
                    det_helper(newloc,locs[k],ds[j],z);
					temp=modded_lattice_z(z[0],z[1]);
                    product*=temp;
                }
                newMatrix(i,j)=pow(in_determinant_rescaling,Ne-1)*product*omega[supermod(newloc[1]*ds[j][0] - newloc[0]*ds[j][1],2*NPhi)];
            }
            else if(i!=electron){//all other elements just need to be updated by the ratio of a sigma function
                if(newMatrix(i,j)==0.){//if zero, need to recompute the whole thing
                    product=1.;
                    for(int k=0;k<Ne;k++){
                        if(k==i) continue;
                        if(k==electron){
                            det_helper(locs[i],newloc,ds[j],z);
                        }
                        else{
							det_helper(locs[i],locs[k],ds[j],z);
                        }
                        temp=modded_lattice_z(z[0],z[1]);
                        product*=temp;
                    }
                    newMatrix(i,j)=pow(in_determinant_rescaling,Ne-1)*product*omega[supermod(locs[i][1]*ds[j][0] - locs[i][0]*ds[j][1],2*NPhi) ];
                }
                else if(newMatrix(i,j)!=0.){
                    det_helper(locs[i],locs[electron],ds[j],z);
                    temp=modded_lattice_z(z[0],z[1]);
                    newMatrix(i,j)/=temp;
					det_helper(locs[i],newloc,ds[j],z);
                    temp=modded_lattice_z(z[0],z[1]);
                    newMatrix(i,j)*=temp;
                }
            }
        }
    }
    
    detSolver.compute(newMatrix);
    new_det=detSolver.determinant();
    //		cout<<electron<<endl;
    //		cout<<oldMatrix<<endl<<endl<<newMatrix<<endl;
}

LATTICE::~LATTICE(){  }
