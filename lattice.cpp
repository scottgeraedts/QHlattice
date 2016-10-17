#include "lattice.h"
int supermod(int k, int n){	return ((k %= n) < 0) ? k+n : k; }

LATTICE::LATTICE() {
    Ne=0;
}
LATTICE::LATTICE(int Ne_t, int invNu_t, bool testing_t, string type_t, int seed, int gs_t, double theta_t, double alpha_t):Ne(Ne_t),invNu(invNu_t),testing(testing_t),type(type_t),gs(gs_t),theta(theta_t),alpha(alpha_t){
	//various parameters from input file
	NPhi=Ne*invNu;
	if(type=="laughlin-hole") NPhi++;
    if (type=="FilledLL" && invNu!=1) {
        cout<<"FilledLL state, however invNu!=1."<<endl;
        exit(0);
    }
    
    //set in_determinant_rescaling.
    //It has been set for more Ne for invNu=2, some Ne for invNu=4.
    if (type=="CFL")
        in_determinant_rescaling=get_in_det_rescaling(Ne, invNu);
    
    //make l1, l2 through theta, alpha.
    L1=sqrt(1.*M_PI*NPhi/sin(theta))*alpha;
    L2=sqrt(1.*M_PI*NPhi/sin(theta))/alpha*polar(1.,theta);

	ran.seed(seed);
	fermions=true;
    
	one=1; zero=0; //useful for fortran calls

	//****initialize z's, w's, d's
	locs=vector< vector<int> >(Ne, vector<int>(2,0)); //initalize locations of all electrons

	double center_frac[2]={0., 0.};
    dbar_parameter=vector<double>(2);
    
	if(type=="CFL"){
		if(Ne%2==0){center_frac[0]=0.5/(1.*Ne); center_frac[1]=0.5/(1.*Ne);}
		make_fermi_surface(center_frac, Ne);
        print_ds();
	}
    else if (type=="laughlin" || type=="laughlin-hole" || type=="FilledLL") {
        ds.clear();
    }
    else {cout<<"recognized type."<<endl; exit(0);}
    
    //setting the ws. Note that the sum of these is ALWAYS zero, adding things like composite fermion momenta or holes doesn't change this.
    //in the y direction these take the values gs*L/invNu, where gs in (0,invNu-1) is an integer which labels the ground state
    ws0=vector< vector<double> > (invNu, vector<double>(2,0) );
    for( int i=0;i<invNu;i++){
        ws0[i][0]=((i+0.5)/(1.*invNu)-0.5);
        ws0[i][1]=gs/(1.*invNu);
    }
    
    if (type=="CFL") set_ds(ds);//set ds, and reset ws.
    else if (type=="FilledLL") set_zeros(vector<double>{0., 0.});
    else ws=ws0;
    
	holes_set=false;
    //*********************
    //To Avoid Bugs, 'set_ws' must be followed by 'set_ds', 'change_dbar_parameter' must following 'set_ds'.

	//********calls to duncan's functions
	set_l_(&NPhi, &L1, &L2);
	setup_z_function_table_();
	int *sl2z=new int[4];
	sl2z[0]=1; sl2z[1]=0; sl2z[2]=0; sl2z[3]=1;
	if(type!="laughlin-hole") setup_laughlin_state_(&Ne,&invNu,sl2z,&gs);
//	cout<<"starting weight "<<running_weight<<endl;
    delete [] sl2z;

	//*****some counters
	setup_coulomb();
    omega=vector<complex<double>>(2*NPhi);
	for(int i=0;i<2*NPhi;i++) omega[i]=polar(1.,M_PI*i/(1.*NPhi));

    sq=vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi, 0));
    sq2=vector<vector<double>>(NPhi, vector<double>(NPhi, 0));
    sq_mqy=vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi, 0));
    sq2_mqy=vector<vector<double>>(NPhi, vector<double>(NPhi, 0));
    SMAq=vector<vector<double>>(NPhi, vector<double>(NPhi, 0.));
//	sq3=vector<vector<vector<vector<complex<double>>>>>(NPhi, vector<vector<vector<complex<double>>>>(NPhi, vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi,0))));
    //comment out sq3, otherwise memory exceed easilly for large invNu state.
}
double LATTICE::get_in_det_rescaling(int Ne, int invNu){
    double rescaling=1.;
    if (invNu==1) {
        rescaling=1.0;
    }
    else if (invNu==2) {
        if (Ne<15) rescaling=1.;
        else if (Ne>=15 && Ne<=38) rescaling=0.28;
        else {rescaling=0.28; cout<<"Please set in_determinant_rescaling if doing berry phase."<<endl;}
    }
    else if (invNu==4) {
        if (Ne<15) rescaling=0.2;
        else if (Ne>=15 && Ne<=19) rescaling=0.12;
        else if (Ne==20) rescaling=0.105;
        else if (Ne==21) rescaling=0.095;
        else {cout<<"Please set in_determinant_rescaling."<<endl; exit(0);}
    }
    else {cout<<"Please set in_determinant_rescaling."<<endl; exit(0);}
    return rescaling;
}
void LATTICE::step(int Nsteps){
	for(int i=0;i<Nsteps*Ne;i++){
		tries++;
		accepts+=simple_update();
//        cout<<"acctpes="<<accepts<<endl;
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
    //move the declaration of 'newDeterminant' and 'newMatrix' from ''***Determinant Part'' to Here.
    complex<double> newDeterminant;
    Eigen::MatrixXcd newMatrix=oldMatrix;

    if (type=="FilledLL") {
        //Here I just use 'newDeterminant' to store the FilledLL wf value, though it has nothing to do with determinant.
        for (int n=0; n<NPhi; n++) {
            newMatrix(electron, n)=1.;
            for (int i=0; i<NPhi; i++) {
                newMatrix(electron, n)*=modded_lattice_z(newloc[0]-i,newloc[1]-n)*polar(1., (zeros[n][i][0]*newloc[1]-zeros[n][i][1]*newloc[0])*M_PI/NPhi);
            }
        }
        detSolver.compute(newMatrix);
        newDeterminant=detSolver.determinant();
        prob+=log(norm(newDeterminant/oldDeterminant));
    }
    else {
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
        int oldCOM[2], newCOM[2];
        sum_locs(oldCOM);
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
            double dx,dy;
            for( int i=0;i<invNu;i++){
                dx=oldCOM[0]/(1.*NPhi)-ws[i][0];
                dy=oldCOM[1]/(1.*NPhi)-ws[i][1];
                z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                prob-=log(norm(temp));
                dx=newCOM[0]/(1.*NPhi)-ws[i][0];
                dy=newCOM[1]/(1.*NPhi)-ws[i][1];
                z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                prob+=log(norm(temp));
            }
            //		get_laughlin_cm_(oldCOM,&temp);
            //		prob-=log(norm(temp));	
            //		get_laughlin_cm_(newCOM,&temp);
            //		prob+=log(norm(temp));	
        }
        ///***********determinant part
//        complex<double> newDeterminant;
//        Eigen::MatrixXcd newMatrix=oldMatrix;
        if(type=="CFL"){
            make_CFL_det(newMatrix, newloc, electron, newDeterminant);
            prob+=log(norm(newDeterminant/oldDeterminant));
        }
        //******phase. no, since calculating weight.
    }
	  
    //*******************update or not
	bool update=false;
	if(prob>0) update=true;
	else if(ran.rand()<exp(prob)) update=true;
	
	if(update){
		locs[electron]=newloc;
		running_weight+=prob;
//		cout<<prob<<endl;
		if(testing) cout<<running_weight<<" "<<get_weight(locs)<<" "<<log(norm(get_wf(locs)))<<endl;
        oldDeterminant=newDeterminant;
        oldMatrix=newMatrix;
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
    
    if (type=="FilledLL") {
        out=log(norm(FilledLL(zs)));
    }
    else {
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
                out+=log(norm( pow(temp,vandermonde_exponent) ));
            }
        }
        
        //com part
        int COM[2]={0,0};
        for( int i=0;i<Ne;i++){
            COM[0]+=zs[i][0];
            COM[1]+=zs[i][1];
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
                    M(i,j)=product*polar(pow(in_determinant_rescaling,Ne-1), 2*M_PI*NPhi*(zs[i][1]*ds[j][0] - zs[i][0]*ds[j][1])/(2.*invNu*NPhi*Ne) );
                }
            }
            detSolver.compute(M);
            temp=detSolver.determinant();
            //ameliorate large determinants by shrinking before taking their norm
            oldDivisor=abs(real(temp))+abs(imag(temp));
            out+=log(norm(temp/oldDivisor))+2*log(oldDivisor);
        }
        
        //phase previously missing.
        vector<double> wsum(2);
        for (int i=0; i<invNu; i++) {
            wsum[0]+=ws[i][0];
            wsum[1]+=ws[i][1];
        }
        complex<double> w_comp = wsum[0]*L1+wsum[1]*L2, zcom_comp = 1.*COM[0]/(1.*NPhi)*L1+1.*COM[1]/(1.*NPhi)*L2, tmp;
        if (type=="laughlin" || type=="laughlin-hole") {
            tmp=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
        }
        else if(type=="CFL") {
            complex<double> dsum_comp = 1.*dsum[0]/NPhi*L1 + 1.*dsum[1]/NPhi*L2;
            tmp=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
        }
        out+=log(norm(tmp));
    }
    
	return out;
} 
//given both a set of positions and a set of ds, computes the wavefunction (NOT the norm of the wavefunction)
complex<double> LATTICE::get_wf(const vector< vector<int> > &zs){
    if (zs.size()!=Ne) {
        cout<<"cannot get_wf because zs.size()!=Ne"<<endl;
        exit(0);
    }
	complex<double> out=1,temp;
	int ix,iy;
	double x,y;
    
    if (type=="FilledLL") {
        out=FilledLL(zs);
    }
    else {
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
        //    cout<<" vandermonde piece = "<<out<<endl;
        complex<double> assist=out;
        
        //COM piece
        int COM[2]={0,0};
        for( int i=0;i<Ne;i++){
            COM[0]+=zs[i][0];
            COM[1]+=zs[i][1];
        }
        double dx,dy;
        if(type=="laughlin-hole"){
            for( int i=0;i<invNu;i++){
                dx=COM[0]/(1.*NPhi)-ws[i][0]+hole[0]/(1.*invNu);
                dy=COM[1]/(1.*NPhi)-ws[i][1]+hole[1]/(1.*invNu);
                z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                out*=temp;
            }
        }
        else{
            for( int i=0;i<invNu;i++){
                dx=COM[0]/(1.*NPhi)-ws[i][0];
                dy=COM[1]/(1.*NPhi)-ws[i][1];
                z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                out*=temp;
            }
        }
        //    cout<<" com piece = "<<out/assist<<endl; assist=out;
        
        //Determinant piece
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
                    M(i,j)=product*polar(pow(in_determinant_rescaling, Ne-1), 2*M_PI*NPhi*(zs[i][1]*ds[j][0] - zs[i][0]*ds[j][1])/(2.*invNu*NPhi*Ne) );
                }
            }
            detSolver.compute(M);
            out=out*detSolver.determinant();
            //        cout<<"det piece = "<<detSolver.determinant()<<endl;
        }
        
        //phase previously missing.
        vector<double> wsum(2);
        for (int i=0; i<invNu; i++) {
            wsum[0]+=ws[i][0];
            wsum[1]+=ws[i][1];
        }
        complex<double> w_comp = wsum[0]*L1+wsum[1]*L2, zcom_comp = 1.*COM[0]/NPhi*L1+1.*COM[1]/NPhi*L2;
        if (type=="laughlin" || type=="laughlin-hole") {
            out*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
        }
        else if(type=="CFL") {
            complex<double> dsum_comp = 1.*dsum[0]/NPhi*L1 + 1.*dsum[1]/NPhi*L2;
            out*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
        }
    }
    return out;
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

vector< vector<int> > LATTICE::get_locs(){return locs;}

vector< vector<int> > LATTICE::get_ds(){return ds;}//'ds' are private, so call 'get_ds' function.

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
			coulomb_table[m][n]=new_v_coulomb_(&NPhi,&m,&n,&eL1,&eL2);
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
 x0*l1+y0*l2 is the center position of fermi sea, therefore x0 y0 are fractional number.
 ds are defined on L/Ne lattice.
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
		ds.push_back(*it);
        it=d_list.erase(it);
    }
//    cout<<"ds.size="<<ds.size()<<endl;
}
void LATTICE::print_ds(){
	ofstream dout("ds");
	for (unsigned int i = 0; i < ds.size(); i += 1)
		dout<<ds[i][0]<<" "<<ds[i][1]<<endl;
	dout.close();
}
void LATTICE::print_ws(){
    ofstream wsout("ws");
    for (unsigned i=0; i<ws.size(); i++) {
        wsout<<ws[i][0]<<" "<<ws[i][1]<<endl;
    }
    wsout.close();
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
	vector<vector<complex<double>>>temp(NPhi,vector<complex<double>>(NPhi,0.));
    vector<vector<complex<double>>>temp2(NPhi,vector<complex<double>>(NPhi,0.));
	for(int qx=0;qx<NPhi;qx++){
		for(int qy=0; qy<NPhi; qy++){
			for(int i=0;i<Ne;i++){
                temp[qx][qy]+=omega[supermod(2*(qx*locs[i][1]-qy*locs[i][0]),2*NPhi)];
                //Using this definition, 'rho(q)' is taken as e^{iq\times z} where q has spatial index.
                temp2[qx][qy]+=omega[supermod(2*(qx*locs[i][1]+qy*locs[i][0]),2*NPhi)];
                //temp2 restore qy<=0 components.
			}
            sq[qx][qy]+=temp[qx][qy];
            sq2[qx][qy]+=norm(temp[qx][qy]);
            sq_mqy[qx][qy]+=temp2[qx][qy];
            sq2_mqy[qx][qy]+=norm(temp2[qx][qy]);
            SMAq[qx][qy]+=coulomb_energy()*norm(temp[qx][qy]);
		}
	}
//    cout<<sq[0][0]<<endl;
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
void LATTICE::print_structure_factors(int nMeas, string filename){
	ofstream sqout(type+"sq/sq"+filename), sqout2(type+"sq/sq2"+filename), sqout3(type+"sq/sq3"+filename);
//    ofstream sqoutfl(type+"sq/sqfl"+filename), sqout2fl(type+"sq/sq2fl"+filename);//'fl'=first line.
	for(int qx=0;qx<NPhi;qx++){
		for(int qy=0; qy<NPhi; qy++){
			sqout2<<sq2[qx][qy]/(1.*nMeas)<<" ";
			sqout<<abs(sq[qx][qy]/(1.*nMeas))<<" ";
		}
		sqout<<endl;
		sqout2<<endl;
	}
    ofstream sqout_mqy(type+"sq/sq_mqy"+filename), sqout2_mqy(type+"sq/sq2_mqy"+filename), sqout3_mqy(type+"sq/sq3_mqy"+filename);
    for(int qx=0;qx<NPhi;qx++){
        for(int qy=0; qy<NPhi; qy++){
            sqout2_mqy<<sq2_mqy[qx][qy]/(1.*nMeas)<<" ";
            sqout_mqy<<abs(sq_mqy[qx][qy]/(1.*nMeas))<<" ";
        }
        sqout_mqy<<endl;
        sqout2_mqy<<endl;
    }
    ofstream smaout(type+"sq/sma"+filename);
    for (int qx=0; qx<NPhi; qx++) {
        for (int qy=0; qy<NPhi; qy++) {
            smaout<<SMAq[qx][qy]/(1.*nMeas)<<" ";
        }
        smaout<<endl;
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
//    sqoutfl.close();
//    sqout2fl.close();
	sqout3.close();
    sqout3_mqy.close();
	sqout.close();
	sqout2.close();
    sqout_mqy.close();
    sqout2_mqy.close();
    smaout.close();
}

complex<double> LATTICE::formfactor(int qx, int qy){
    complex<double> out=0;
    complex<double> temp;
    double kappa=2*M_PI/real(L1);
    for(int px=-9;px<10;px++){
        temp=exp(-0.25*pow(kappa*(px*NPhi-qx),2));
        for(int py=-9;py<10;py++){
            out+=temp*exp(-0.25*pow(kappa*(py*NPhi-qy),2));
        }
    }
    return out;
}//need to change if want to do non-square case.

complex<double> LATTICE::rhoq(int qx, int qy, const vector< vector<int> > &zs){
	complex<double> out=0;
	for(int i=0;i<Ne;i++){
//		out+=omega[supermod((2*qx*locs[i][0]+2*qy*locs[i][1]), 2*NPhi)];//need to be checked.
        out+=omega[supermod((2*qx*locs[i][1]-2*qy*locs[i][0]), 2*NPhi)];
	}
//	return out/(formfactor(qx,qy)*(1.*NPhi));
    return out;
}
void LATTICE::reset(){
	tries=0; accepts=0;
	cold_start();
    
	//**** setting up the initial determinant matrix
	running_weight=0;
	int site=0;
	int initial_state_counter=0;
	
    while(running_weight==0){
	//for some sizes the configuration specified by cold_start has zero weight
	//if that happens fiddle around until you find a better configuration, thats why theres a while loop
		locs[site][0]=locs[p(site)][0];
		site++;
		if(site==Ne) cout<<"couldn't easily find a good configuration!"<<endl;
		
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
        if (type=="FilledLL") {
            oldMatrix=Eigen::MatrixXcd(NPhi, NPhi);
            oldDeterminant=get_wf(locs);//it is not determinant in this case, just a name trick.
            for (int m=0; m<NPhi; m++) {
                for (int n=0; n<NPhi; n++) {
                    oldMatrix(m,n)=1.;
                    for (int i=0; i<NPhi; i++) {
                        oldMatrix(m,n)*=modded_lattice_z(locs[m][0]-i, locs[m][1]-n)*polar(1., 1.*(i*locs[m][1]-n*locs[m][0])*M_PI/NPhi/NPhi);
                        //zeros[gs][i][0]=i; zeros[gs][i][1]=gs;
                    }
                }
            }
        }
        
		initial_state_counter++;
		if(initial_state_counter>1000){
			cout<<"couldn't find a good starting configuration"<<endl;
			exit(0);
		}
	}
    
    sq=vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi,0));
    sq2=vector<vector<double>>(NPhi, vector<double>(NPhi,0));
    sq_mqy=vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi,0));
    sq2_mqy=vector<vector<double>>(NPhi, vector<double>(NPhi,0));
    SMAq=vector<vector<double>>(NPhi, vector<double>(NPhi,0));
    
    check_sanity();
}
//checks a few different things to make sure that they make sense
void LATTICE::check_sanity(){
	if(type!="CFL" && type!="laughlin" && type!= "laughlin-hole" && type!="FilledLL"){
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
    if(type=="FilledLL" && (zeros.size()!=NPhi || Ne!=NPhi)) {
        cout<<"Filled LL zeros size or Ne/NPhi wrong."<<endl;
        exit(0);
    }
}
//changes the dbar (which ordinarily is the sum of d's) to whatever we want
void LATTICE::change_dbar_parameter(double dbarx, double dbary){
	if(type!="CFL" && type!="FilledLL"){
		cout<<"changing dbar, but not CFL nor FilledLL."<<endl;
		//exit(0);
	}
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
    if (tds.size()!=Ne) {
        cout<<"cannot set ds, ds.size() wrong."<<endl;
        exit(0);
    }
    ds=tds; dsum=vector<int>(2,0);
	for(int i=0;i<Ne;i++){
        dsum[0]+=ds[i][0]*invNu; dsum[1]+=ds[i][1]*invNu; //'ds' is on L/Ne lattice, 'dsum' is on L/Nphi lattice.
	}
	//print_ds();
	change_dbar_parameter(dsum[0]/(1.*Ne),dsum[1]/(1.*Ne));
    
    //reset ws, according to ds.
    ws.clear(); ws=ws0;
    if (type=="CFL") {
        for (int i=0; i<ws.size(); i++) {
            ws[i][0]+=dsum[0]/(1.*invNu*NPhi); ws[i][1]+=dsum[1]/(1.*invNu*NPhi);// if 'd = wsum'.
        }
    }
}
void LATTICE::set_ws(vector<vector<double>> ws_t){
    if (ws_t.size()!=invNu) {
        cout<<"cannot set ws, ws.size() wrong."<<endl;
        exit(0);
    }
    for (unsigned i=0; i<invNu; i++) {
        if (ws_t[i].size()!=2) {
            cout<<"cannot set ws, ws[i].size() wrong."<<endl;
            exit(0);
        }
    }
    ws0=ws_t; ws=ws0;
    if (type=="CFL") {
        for (int i=0; i<ws.size(); i++) {
            ws[i][0]+=dsum[0]/(1.*invNu*NPhi); ws[i][1]+=dsum[1]/(1.*invNu*NPhi);// if 'd = wsum'.
        }
    }
}
void LATTICE::set_zeros(vector<double> zeros0){
    if (zeros0.size()!=2) {
        cout<<"cannot set zeros, since its dim != 2."<<endl;
        exit(0);
    }
    zeros=vector<vector<vector<double>>> (NPhi, vector<vector<double>>(NPhi, vector<double>(2)));
    for (int gs=0; gs<NPhi; gs++) {
        for( int i=0;i<NPhi;i++){
            zeros[gs][i][0]=i/(1.*NPhi)+zeros0[0];
            zeros[gs][i][1]=gs/(1.*NPhi)+zeros0[1];
        }
    }
    change_dbar_parameter(-zeros[0][0][0]*NPhi, -zeros[0][0][1]*NPhi);//the aim to do this is to initialize 'shifted_ztable', and 'modded_lattice_z'.
}
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

void LATTICE::cold_start(){
	for(int i=0;i<Ne;i++){
		locs[i][0]=i;
		locs[i][1]=i;
	}
}

//call's duncan's lattice_z function, if the arguments x or y are outside of the range (0,NPhi) it shifts them into that range
//and multiplies by the appropriate phase
//only works for a square torus
complex<double> LATTICE::modded_lattice_z(int x, int y){
	int modx=supermod(x,NPhi);
	int mody=supermod(y,NPhi);
//	complex<double> out=lattice_z_(&NPhi,&modx,&mody,&L1,&L2,&one);
	complex<double> out=shifted_ztable[modx][mody];
	int j=(modx-x)/NPhi, k=(mody-y)/NPhi;
//	out*=omega[supermod(-(mody+dbar_parameter[1])*j+(modx+dbar_parameter[0])*k,2*NPhi)];
	out*=polar(1., (-(mody+dbar_parameter[1])*j+(modx+dbar_parameter[0])*k)*M_PI/(1.*NPhi));
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
//    cout<<electron<<endl;
//    cout<<oldMatrix<<endl<<endl<<newMatrix<<endl;
}

vector<double> LATTICE::get_dbar_parameter(){
    return dbar_parameter;
}

complex<double> LATTICE::FilledLL(vector<vector<int>> z){
    Eigen::MatrixXcd landauwf(NPhi, NPhi);
    for (int m=0; m<NPhi; m++) {
        for (int n=0; n<NPhi; n++) {
            landauwf(m,n)=1.;
            for (int i=0; i<zeros[n].size(); i++) {
                complex<double> temp=modded_lattice_z(z[m][0]-i, z[m][1]-n);
                landauwf(m,n)*=temp*polar(1., (zeros[n][i][0]*z[m][1]-zeros[n][i][1]*z[m][0])*M_PI/NPhi);
            }
        }
    }
//    cout<<"matrix="<<endl<<landauwf<<endl;
    detSolver.compute(landauwf);
    return detSolver.determinant();
}

LATTICE::~LATTICE(){  }
