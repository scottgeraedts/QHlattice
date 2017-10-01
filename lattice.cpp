#include "lattice.h"
//int supermod(int k, int n){	return ((k %= n) < 0) ? k+n : k; }
inline double laguerre(int n, double x){
    double result=0;
    for (int k=0; k<=n; k++) {
        if (k%2==0) result+=pow(x,k)*comb(n,k)/factorial(k);
        else result-=pow(x,k)*comb(n,k)/factorial(k);
    }
    return result;
}
LATTICE::LATTICE():shiftx(0.),shifty(0.){
    Ne=0;
}
LATTICE::LATTICE(LATTICE_PARAMS params):shiftx(0.),shifty(0.){
	Ne=params.Ne;
	invNu=params.invNu;
	testing=params.testing;
	type=params.type;
	gs=params.gs;
	theta=params.theta;
	alpha=params.alpha;
	trace=params.trace;
	w_delta=params.w_delta;
	dbar_delta=params.dbar_delta;
	testing=params.testing;
	in_determinant_rescaling=params.rescale;

	init(params.seed);
}
LATTICE::LATTICE(int Ne_t, int invNu_t, bool testing_t, string type_t, int seed, int gs_t, double theta_t, double alpha_t, double shift_x, double shift_y, bool trace_t, bool correlatedsampling_t):Ne(Ne_t),invNu(invNu_t),testing(testing_t),type(type_t),gs(gs_t),theta(theta_t),alpha(alpha_t),shiftx(shift_x),shifty(shift_y),trace(trace_t),correlatedsampling(correlatedsampling_t){

 	init(seed);
    in_determinant_rescaling=get_in_det_rescaling(Ne, invNu);
}
void LATTICE::init(int seed){

	if(trace)
		cout<<"WARNING! Trace=true is no longer supported, use lattice_wrapper instead"<<endl;
    if (correlatedsampling && (type!="CFL")) {
        cout<<"correlatedsampling is enabled for CFL only"<<endl;
        exit(0);
    }

	//various parameters from input file
 	NPhi=Ne*invNu;
	if(type=="laughlin-hole") NPhi++;
    if (type=="FilledLL" && invNu!=1) {
        cout<<"FilledLL state, however invNu!=1."<<endl;
        exit(0);
    }

    //make l1, l2 through theta, alpha.
    L1=sqrt(1.*M_PI*NPhi/sin(theta))*alpha;
    L2=sqrt(1.*M_PI*NPhi/sin(theta))/alpha*polar(1.,theta);
 	
	//********calls to duncan's functions
	set_l_(&NPhi, &L1, &L2);
    setup_z_function_table_();
	int *sl2z=new int[4];
	sl2z[0]=1; sl2z[1]=0; sl2z[2]=0; sl2z[3]=1;
	if(type!="laughlin-hole") setup_laughlin_state_(&Ne,&invNu,sl2z,&gs);
//	cout<<"starting weight "<<running_weight<<endl;
    delete [] sl2z;
    
	ran.seed(seed);
	fermions=true;
    
	one=1; zero=0; //useful for fortran calls

	//****initialize z's, w's, d's
	if(type=="doubledCFL") locs=vector< vector<int> >(2*Ne, vector<int>(2,0)); //initalize locations of all electrons
	else locs=vector< vector<int> >(Ne, vector<int>(2,0)); //initalize locations of all electrons

	double center_frac[2]={0., 0.};
    dbar_parameter=vector<double>(2,0.);
    
	if(type=="CFL" or type=="doubledCFL"){
		if(Ne%2==0){center_frac[0]=0.5/(1.*Ne); center_frac[1]=0.5/(1.*Ne);}
		make_fermi_surface(center_frac, Ne);
        print_ds();
	}
    else if (type=="laughlin" || type=="laughlin-hole" || type=="FilledLL" ) {
        ds.clear();
    }
    else {cout<<"recognized type."<<endl; exit(0);}
    
    //setting the ws. Note that the sum of these is ALWAYS zero, adding things like composite fermion momenta or holes doesn't change this.
    //in the y direction these take the values gs*L/invNu, where gs in (0,invNu-1) is an integer which labels the ground state
    ws0=vector< vector<double> > (invNu, vector<double>(2,0) );
    for( int i=0;i<invNu;i++){
        ws0[i][0]=((i+0.5)/(1.*invNu)-0.5)+real(w_delta)*lil_sign(gs)*lil_sign(i);
        ws0[i][1]=gs/(1.*invNu)+imag(w_delta)*lil_sign(gs)*lil_sign(i);
    }
    
    //change the boundary condition. For CFL state, bc is set by set_ds().
    ws=ws0;
    for (int i=0; i<invNu; i++) {
        ws[i][0]+=shiftx/(1.*invNu);
        ws[i][1]+=shifty/(1.*invNu);
    }
    
    //alternatively, in the CFL case we can access different ground states by moving the ds (to do this uncomment the next line)
   // if(type=="CFL") for(int i=0; i<Ne; i++) ds[i][1]+=gs;
    if (type=="CFL" or type=="doubledCFL")
        set_ds(ds);//set ds, and reset ws.
    else if (type=="FilledLL")
        set_zeros(vector<double>{0., 0.});
    
    wsum=vector<double>(2,0.);
    for (int i=0; i<invNu; i++) {
        wsum[0]+=ws[i][0];
        wsum[1]+=ws[i][1];
    }
    
	holes_set=false;
    //*********************
    //To Avoid Bugs, 'set_ws' must be followed by 'set_ds', 'change_dbar_parameter' must following 'set_ds'.
    
    omega=vector<complex<double>>(2*NPhi);
	for(int i=0;i<2*NPhi;i++) omega[i]=polar(1.,M_PI*i/(1.*NPhi));

    setup_coulomb0();
    sq=vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi, 0));
    sq2=vector<vector<double>>(NPhi, vector<double>(NPhi, 0));
    sq_mqy=vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi, 0));
    sq2_mqy=vector<vector<double>>(NPhi, vector<double>(NPhi, 0));
    sq3=vector<vector<vector<vector<complex<double>>>>> (NPhi, vector<vector<vector<complex<double>>>>(NPhi, vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi,0.))));
}
vector<double> LATTICE::get_shift(){return vector<double>{shiftx, shifty};}
vector<double> LATTICE::get_bc(){if (type=="CFL" or type=="doubledCFL") return vector<double>{wsum[0]-dsum[0]/(1.*NPhi), wsum[1]-dsum[1]/(1.*NPhi)}; else return wsum;}
double LATTICE::get_in_det_rescaling(int Ne, int invNu){
    double rescaling=1.;
    if (type=="CFL") {
        if (invNu==1) {
            rescaling=1.0;
        }
        else if (invNu==2) {
            if (Ne<40) rescaling=0.25;
            else if (Ne>=40 && Ne<=50) rescaling=0.2;
            else if (Ne<90) rescaling=0.18;
            else if (Ne==97) rescaling=0.15;
            else {rescaling=0.15; cout<<"Please set in_determinant_rescaling if doing berry phase."<<endl;}
        }
        else if (invNu==4) {
            if (Ne<15) rescaling=0.1;
            else if (Ne>=15 && Ne<=21) rescaling=0.08;
            else {cout<<"Please set in_determinant_rescaling."<<endl; exit(0);}
        }
    }
    else if (type=="laughlin"||type=="laughlin-hole") {
        if (invNu==3) {
            if (Ne<35) {
                rescaling=0.87;
            }
            else if (Ne==50) {
                rescaling=0.9;
            }
            else {cout<<"Please set in_determinant_rescaling."<<endl; exit(0);}
        }
        if (invNu==5) {
            if (Ne<=20) {
                rescaling=0.7;
            }
            else if (Ne==50) {
                rescaling=0.812;
            }
            else {cout<<"Please set in_determinant_rescaling."<<endl; exit(0);}
        }
    }
    return rescaling;
}
void LATTICE::step(int Nsteps){
	for(int i=0;i<Nsteps*Ne;i++){
		tries++;
		accepts+=simple_update();
	}
}
int LATTICE::simple_update(){
	//find a new position for a randomly chosen electron
    int electron;
	if(type=="doubledCFL") electron=ran.randInt(2*Ne-1);
	else electron=ran.randInt(Ne-1);
	vector<int> newloc=random_move(locs[electron], NPhi, ran);
	vector<vector<int>>::iterator it=find(locs.begin(),locs.end(),newloc);
	if(it!=locs.end()) return 0;
    bool update=false;
    
    //newMatrix=oldMatrix, newDet=oldDet; they're global.
    update_matdet(false);
    if (correlatedsampling) for (int i=0; i<nonsamplestates.size(); i++) nonsamplestates[i].update_matdet(false);
    
    double prob=updateweight(electron, newloc);
    //this function changes newMatrix. whether oldMatrix is changed or not dependent on if update.
    
    if (testing && !correlatedsampling)
        cout<<running_weight<<" "<<get_weight(locs)<<" "<<log(norm(get_wf(locs)))<<endl<<endl;
    else if (testing && correlatedsampling) {
        cout<<c_running_weight<<" "<<get_correlated_weight()<<endl;
        cout<<"gs ="<<running_weight<<" "<<get_weight(locs)<<" "<<log(norm(get_wf(locs)))<<endl;
        for (int i=0; i<nonsamplestates.size(); i++) {
            cout<<"i="<<i<<" "<<nonsamplestates[i].running_weight<<" "<<nonsamplestates[i].get_weight(locs)<<" "<<log(norm(nonsamplestates[i].get_wf(locs)))<<" ";
            
            vector<vector<int>> tmp=nonsamplestates[i].get_locs();
            for (int j=0; j<tmp.size(); j++) {
                //cout<<tmp[j][0]<<" "<<tmp[j][1]<<" ";
            }
            cout<<endl;
        }
    }
    
    if (!correlatedsampling) {
        //*******************update or not(prob)
        if(prob>0) update=true;
        else if(ran.rand()<exp(prob)) update=true;
        
        if(update){
            locs[electron]=newloc;
            running_weight+=prob;
            update_matdet();
        }
        //else if(type=="doubledCFL") locs[electron]=oldloc;
    }
    else if (correlatedsampling) {
        double prob_c=0.;//correlated sampling probablity;
        
        int Nn=nonsamplestates.size();
        vector<double> dif_wt(Nn, 0.), ins_wt(Nn, 0.);
        //dif_wt[i]=running_weight[i]-running_weight; ins_wt[i]=updateweight[i]. i labels nonsampled sates.
        
        for (int i=0; i<Nn; i++) {
            dif_wt[i]=nonsamplestates[i].running_weight-running_weight;
            ins_wt[i]=nonsamplestates[i].updateweight(electron, newloc);
        }
        //till this line, gs and excited states' newMatrix/newDeterminant got updated.
        
        double den=ratio[nonsamplestates.size()], num=ratio[nonsamplestates.size()]*exp(prob);
        for (int i=0; i<Nn; i++) {
            den+=ratio[i]*exp(dif_wt[i]);
            num+=ratio[i]*exp(ins_wt[i]+dif_wt[i]);//TODO:correlated-sampling
        }
        prob_c=log(num/den);
        
        //*******************update or not (prob_c)
        if (prob_c>0) update=true;
        else if(ran.rand()<exp(prob_c)) update=true;
        
        if (update) {
            running_weight+=prob;
            c_running_weight+=prob_c;
            
            locs[electron]=newloc;
            update_matdet();//update Matrix/Det for ground state.
            for (int i=0; i<nonsamplestates.size(); i++) {
                nonsamplestates[i].change_locs(locs);
                nonsamplestates[i].running_weight+=ins_wt[i];
                //cout<<"nonsamplestates[i].running_weight="<<nonsamplestates[i].running_weight<<endl;
                nonsamplestates[i].update_matdet();
            }
        }
    }
    
    return update? 1:0;
}
double LATTICE::updateweight(const int &electron, const vector<int> &newloc){
    double out=0.;
    complex<double> temp;
    
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
        out+=log(norm(newDeterminant/oldDeterminant));
    }else if(type=="doubledCFL"){
        //oldloc=locs[electron];
        vector<vector<int>> newlocs=locs;
        newlocs[electron]=newloc;
        newDeterminant=get_wf(newlocs);
        
        //locs[electron]=newloc;
        //newDeterminant=get_wf(locs);
        out+=log(norm(newDeterminant/oldDeterminant));
    }else {
        //***************hole part
        double dx,dy;
        if(type=="laughlin-hole"){
            dx=(locs[electron][0]/(1.*NPhi)-hole[0]);
            dy=(locs[electron][1]/(1.*NPhi)-hole[1]);
            z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
            out-=log(norm(temp));
            dx=(newloc[0]/(1.*NPhi)-hole[0]);
            dy=(newloc[1]/(1.*NPhi)-hole[1]);
            z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
            out+=log(norm(temp));
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
                out-=log(norm( pow(temp,vandermonde_exponent) ));
                
                //multiply new part
                xi=(locs[i][0]-newloc[0]);
                yi=(locs[i][1]-newloc[1]);
                if(i>electron){
                    xi=-xi; yi=-yi;
                }
                temp=lattice_z_(&NPhi,&xi,&yi,&L1,&L2,&one);
                out+=log(norm( pow(temp,vandermonde_exponent) ));
            }
        }
        
        //***************COM PART
        vector<double> wsum(2); for (int i=0; i<invNu; i++) for (int j=0; j<2; j++) wsum[j]+=ws[i][j];
        complex<double> w_comp = wsum[0]*L1+wsum[1]*L2, dsum_comp;
        if (type=="CFL") {
            dsum_comp = 1.*dsum[0]/NPhi*L1 + 1.*dsum[1]/NPhi*L2;
        }
        
        int oldCOM[2], newCOM[2];
        sum_locs(oldCOM);
        newCOM[0]=oldCOM[0]-locs[electron][0]+newloc[0];
        newCOM[1]=oldCOM[1]-locs[electron][1]+newloc[1];
        
        if (!trace) {
            if(type=="laughlin-hole"){
                double dx,dy;
                for( int i=0;i<invNu;i++){
                    dx=oldCOM[0]/(1.*NPhi)-ws[i][0]+hole[0]/(1.*invNu);
                    dy=oldCOM[1]/(1.*NPhi)-ws[i][1]+hole[1]/(1.*invNu);
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    out-=log(norm(temp));
                    dx=newCOM[0]/(1.*NPhi)-ws[i][0]+hole[0]/(1.*invNu);
                    dy=newCOM[1]/(1.*NPhi)-ws[i][1]+hole[1]/(1.*invNu);
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    out+=log(norm(temp));
                }
            }else{
                double dx,dy;
                for( int i=0;i<invNu;i++){
                    dx=oldCOM[0]/(1.*NPhi)-ws[i][0];
                    dy=oldCOM[1]/(1.*NPhi)-ws[i][1];
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    out-=log(norm(temp));
                    dx=newCOM[0]/(1.*NPhi)-ws[i][0];
                    dy=newCOM[1]/(1.*NPhi)-ws[i][1];
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    //while (temp<1e-15) {TODO: what if wf(COM) = 0.}
                    if (abs(temp)<1e-15) {
//                        cout<<"COM=0"<<endl;
//                        exit(0);
                        return -10000;
                    }
                    out+=log(norm(temp));
                }
                //get_laughlin_cm_(oldCOM,&temp);
                //prob-=log(norm(temp));
                //get_laughlin_cm_(newCOM,&temp);
                //prob+=log(norm(temp));
            }
        }
        else {
            if (type=="laughlin-hole") {
                cout<<" 'trace' is not set up for laughlin-hole"<<endl;
                exit(0);
            }
            complex<double> w_comp0=w_comp;
            complex<double> sum=0, product;
            for (int k=0; k<invNu; k++) {
                product=1.;
                for( int i=0;i<invNu;i++){
                    dx=oldCOM[0]/(1.*NPhi)-ws[i][0];
                    dy=oldCOM[1]/(1.*NPhi)-ws[i][1]-k/(1.*invNu);
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    product*=temp;
                }
                w_comp=w_comp0+1.*k*L2;
                
                complex<double> zcom_comp = 1.*oldCOM[0]/NPhi*L1+1.*oldCOM[1]/NPhi*L2;
                if (type=="laughlin") product*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
                else if (type=="CFL") product*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
                sum+=product;
            }
            out-=log(norm(sum));
           	sum=0.;
            for (int k=0; k<invNu; k++) {
                product=1.;
                for( int i=0;i<invNu;i++){
                    dx=newCOM[0]/(1.*NPhi)-ws[i][0];
                    dy=newCOM[1]/(1.*NPhi)-ws[i][1]-k/(1.*invNu);
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    product*=temp;
                }
                w_comp=w_comp0+1.*k*L2;
                
                complex<double> zcom_comp = 1.*newCOM[0]/NPhi*L1+1.*newCOM[1]/NPhi*L2;
                if (type=="laughlin") product*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
                else if (type=="CFL") product*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
                sum+=product;
            }
            out+=log(norm(sum));
        }
        
        ///***********determinant part
        if(type=="CFL"){
            make_CFL_det(newMatrix, newloc, electron, newDeterminant, locs);
            out+=log(norm(newDeterminant/oldDeterminant));
        }
    }
    return out;
}
//very similar to simple update (above), but it uses a provided set of zs, instead of this objects version
//and it is also provided with an position to update
complex<double> LATTICE::update_weight(const vector< vector<int> > &zs, int electron, vector<int> newloc){
	double out=1.;
	complex<double> temp;

    //***************hole part
    double dx,dy;
    if(type=="laughlin-hole"){
        dx=(zs[electron][0]/(1.*NPhi)-hole[0]);
        dy=(zs[electron][1]/(1.*NPhi)-hole[1]);
        z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
		out/=abs(temp);
        dx=(newloc[0]/(1.*NPhi)-hole[0]);
        dy=(newloc[1]/(1.*NPhi)-hole[1]);
        z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
		out*=abs(temp);
    }
       
    //***************vandermode part
    int vandermonde_exponent=invNu;
    if(type=="CFL") vandermonde_exponent-=2;
    int xi,yi;
    
    if(vandermonde_exponent!=0){
        for(int i=0;i<Ne;i++){
            //divide off old part
            if(i==electron) continue;
            xi=(zs[i][0]-zs[electron][0]);
            yi=(zs[i][1]-zs[electron][1]);
            if(i>electron){
                xi=-xi; yi=-yi;
            }
            temp=lattice_z_(&NPhi,&xi,&yi,&L1,&L2,&one);
            out/=abs(pow(temp,vandermonde_exponent));
            
            //multiply new part
            xi=(zs[i][0]-newloc[0]);
            yi=(zs[i][1]-newloc[1]);
            if(i>electron){
                xi=-xi; yi=-yi;
            }
            temp=lattice_z_(&NPhi,&xi,&yi,&L1,&L2,&one);
            out*=abs(pow(temp,vandermonde_exponent));
        }
    }
    //***************COM PART
    int oldCOM[2]={0,0}, newCOM[2];
    for(auto it=zs.begin(); it!=zs.end(); ++it){
    	oldCOM[0]+=(*it)[0];
    	oldCOM[1]+=(*it)[1];
    }
    newCOM[0]=oldCOM[0]-zs[electron][0]+newloc[0];
    newCOM[1]=oldCOM[1]-zs[electron][1]+newloc[1];
    
    if(!trace){
		if(type=="laughlin-hole"){
		    double dx,dy;
		    for( int i=0;i<invNu;i++){
		        dx=oldCOM[0]/(1.*NPhi)-ws[i][0]+hole[0]/(1.*invNu);
		        dy=oldCOM[1]/(1.*NPhi)-ws[i][1]+hole[1]/(1.*invNu);
		        z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
		        out/=abs(temp);
		        dx=newCOM[0]/(1.*NPhi)-ws[i][0]+hole[0]/(1.*invNu);
		        dy=newCOM[1]/(1.*NPhi)-ws[i][1]+hole[1]/(1.*invNu);
		        z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
				out*=abs(temp);
		    }
		}else{
		    for( int i=0;i<invNu;i++){
		        xi=oldCOM[0]-ws[i][0]*NPhi;
		        yi=oldCOM[1]-ws[i][1]*NPhi;
				temp=modded_lattice_z_nodbar(xi, yi);
		        out/=abs(temp);
		        xi=newCOM[0]-ws[i][0]*NPhi;
		        yi=newCOM[1]-ws[i][1]*NPhi;
				temp=modded_lattice_z_nodbar(xi, yi);
				out*=abs(temp);
		    }
		
			//the jie phase
//		    vector<double> wsum(2); for (int i=0; i<invNu; i++) for (int j=0; j<2; j++) wsum[j]+=ws[i][j];
//		    complex<double> w_comp = wsum[0]*L1+wsum[1]*L2, dsum_comp;
//		    complex<double> zcom_comp=1.*(newCOM[0]-oldCOM[0])/NPhi*L1+1.*(newCOM[1]-oldCOM[1])/NPhi*L2;
//		    if (type=="CFL") dsum_comp = 1.*dsum[0]/NPhi*L1 + 1.*dsum[1]/NPhi*L2;
//			//SLOW LINE
//            if (type=="laughlin") out*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
//            else if (type=="CFL") out*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
    	}
	}else {
        if (type=="laughlin-hole") {
            cout<<" 'trace' is not set up for laughlin-hole"<<endl;
            exit(0);
        }
        vector<double> wsum(2); for (int i=0; i<invNu; i++) for (int j=0; j<2; j++) wsum[j]+=ws[i][j];
        complex<double> w_comp = wsum[0]*L1+wsum[1]*L2, dsum_comp;
        complex<double> w_comp0=w_comp;
        complex<double> sum=0, product;
        for (int k=0; k<invNu; k++) {
        	product=1.;
            for( int i=0;i<invNu;i++){
                dx=oldCOM[0]/(1.*NPhi)-ws[i][0];
                dy=oldCOM[1]/(1.*NPhi)-ws[i][1]-k/(1.*invNu);
                z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                product*=temp;
            }
            w_comp=w_comp0+1.*k*L2;
            complex<double> zcom_comp = 1.*oldCOM[0]/NPhi*L1+1.*oldCOM[1]/NPhi*L2;
            if (type=="laughlin") product*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
            else if (type=="CFL") product*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
            sum+=product;
        }
        out/=abs(sum);
       	sum=0.;
        for (int k=0; k<invNu; k++) {
        	product=1.;
            for( int i=0;i<invNu;i++){
                dx=newCOM[0]/(1.*NPhi)-ws[i][0];
                dy=newCOM[1]/(1.*NPhi)-ws[i][1]-k/(1.*invNu);
                z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                product*=temp;
            }
            w_comp=w_comp0+1.*k*L2;
            
            complex<double> zcom_comp = 1.*newCOM[0]/NPhi*L1+1.*newCOM[1]/NPhi*L2;
            if (type=="laughlin") product*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
            else if (type=="CFL") product*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
            sum+=product;
        }
		out*=abs(sum);
    }

    if(type=="CFL"){
    	newMatrix=oldMatrix;
        make_CFL_det(newMatrix, newloc, electron, newDeterminant, zs);
        out*=abs(newDeterminant)/abs(oldDeterminant);
    }
	return out;	  
}

//chooses a new site, given an old site
vector<int> LATTICE::random_move( const vector<int> &in, int NPhi_t, MTRand &ran_t){
	vector<int>newloc(2);

	//the stupidest way to do this
//	int rint=ran_t.randInt( NPhi_t*NPhi_t-1);
//	newloc[0]=rint/NPhi_t;
//	newloc[1]=rint%NPhi_t;
//	return newloc;

	int hoplength=1;
	int n=pow(2*hoplength+1,2)-1;
	vector<int> newx(n),newy(n);
	vector<double> newprob(n);
	int xi,yi,counter=0;
	for(int i=0;i<=2*hoplength;i++){
		xi=i-hoplength+in[0]; 
		for(int j=0;j<=2*hoplength;j++){
			yi=j-hoplength+in[1];
			if(xi==in[0] && yi==in[1]) continue;
			newx[counter]=supermod(xi,NPhi_t);
			newy[counter]=supermod(yi,NPhi_t);
//			newprob.push_back( exp(-0.5*( pow(in[0]-i,2)+pow(in[1]-j,2) )/pow(hoplength,2) ) );
			newprob[counter]=1.;
			counter++;
		}
	}
	double r=ran_t.rand();
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
    
    if (type=="doubledCFL"){
    	out=log(norm(doubled_CFL(zs)));
    }else {
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
                temp=pow(temp,vandermonde_exponent);
                
                //......
                if (type=="laughlin"||type=="laughlin-hole") {
                    //temp*=pow(in_determinant_rescaling, vandermonde_exponent*(Ne-1));
                    temp*=pow(in_determinant_rescaling, Ne-1);
                }
                out+=log(norm( temp ));
            }
        }
        
        //com part
        int COM[2]={0,0};
        for( int i=0;i<Ne;i++){
            COM[0]+=zs[i][0];
            COM[1]+=zs[i][1];
        }
        vector<double> wsum(2); for (int i=0; i<invNu; i++) for (int j=0; j<2; j++) wsum[j]+=ws[i][j];
        complex<double> w_comp = wsum[0]*L1+wsum[1]*L2, zcom_comp = 1.*COM[0]/NPhi*L1+1.*COM[1]/NPhi*L2, dsum_comp;
        if (type=="CFL") dsum_comp = 1.*dsum[0]/NPhi*L1 + 1.*dsum[1]/NPhi*L2;
        
        if (!trace) {
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
        }
        else {
            if (type=="laughlin-hole") {
                cout<<" trace is not set up for laughlin-hole"<<endl;
                exit(0);
            }
            complex<double> w_comp0=w_comp;
            complex<double> sum=0., product;
            for (int k=0; k<invNu; k++) {
            	product=1.;
                for( int i=0;i<invNu;i++){
                    x=COM[0]/(1.*NPhi)-ws[i][0];
                    y=COM[1]/(1.*NPhi)-ws[i][1]-k/(1.*invNu);
                    z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
                    w_comp=w_comp0+1.*k*L2;
                    product*=temp;
                }
                if (type=="laughlin") product*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
                else if (type=="CFL") product*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
                sum+=product;
            }
            out+=log(norm(sum));
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
                        x=det_helper(zs[i][0],zs[k][0],ds[j][0],dbar_parameter[0])/(1.*NPhi);
                        y=det_helper(zs[i][1],zs[k][1],ds[j][1],dbar_parameter[1])/(1.*NPhi);
                        z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
                        product*=temp;
                    }
                    //this part is only valid on a square torus!
                    M(i,j)=product*polar(pow(in_determinant_rescaling,Ne-1), 2*M_PI*NPhi*(zs[i][1]*ds[j][0] - zs[i][0]*ds[j][1])/(2.*invNu*NPhi*Ne));
                }
            }
//            M/=pow(10,19);
//            cout<<"Matrix=\n"<<M<<endl;
            detSolver.compute(M);
            temp=detSolver.determinant();
//            cout<<"det="<<temp<<endl;exit(0);
            //ameliorate large determinants by shrinking before taking their norm
            oldDivisor=abs(real(temp))+abs(imag(temp));
            out+=log(norm(temp/oldDivisor))+2*log(oldDivisor);
        }
        
    }
	return out;
}
//given both a set of positions and a set of ds, computes the wavefunction (NOT the norm of the wavefunction)
complex<double> LATTICE::get_wf(const vector< vector<int> > &zs){
    if (type!="doubledCFL" and (signed)zs.size()!=Ne) {
        cout<<"cannot get_wf because zs.size()!=Ne "<<zs.size()<<" "<<Ne<<endl;
        exit(0);
    }
    complex<double> out=1,temp;
    int ix,iy;
    double x,y;
    
    if (type=="doubledCFL") out=doubled_CFL(zs);
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
                    out*=pow(lattice_z_(&NPhi,&ix,&iy,&L1,&L2,&one), vandermonde_exponent);
                    
                    if (type=="laughlin"||type=="laughlin-hole") {
                        out*=pow(in_determinant_rescaling, Ne-1);
                    }
                }
            }
        }
        
        //COM piece (together with phase)
        int COM[2]={0,0};
        for( int i=0;i<Ne;i++){
            COM[0]+=zs[i][0];
            COM[1]+=zs[i][1];
        }
        double dx,dy;
        vector<double> wsum(2); for (int i=0; i<invNu; i++) for (int j=0; j<2; j++) wsum[j]+=ws[i][j];
        complex<double> w_comp = wsum[0]*L1+wsum[1]*L2, zcom_comp = 1.*COM[0]/NPhi*L1+1.*COM[1]/NPhi*L2, dsum_comp;
        if (type=="CFL") {
            dsum_comp = 1.*dsum[0]/NPhi*L1 + 1.*dsum[1]/NPhi*L2;
        }
        
        if (!trace) {
            if(type=="laughlin-hole"){
                for( int i=0;i<invNu;i++){
                    dx=COM[0]/(1.*NPhi)-ws[i][0]+hole[0]/(1.*invNu);
                    dy=COM[1]/(1.*NPhi)-ws[i][1]+hole[1]/(1.*invNu);
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    out*=temp;
                }
                out*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
            }
            else{
                for( int i=0;i<invNu;i++){
                    dx=COM[0]/(1.*NPhi)-ws[i][0];
                    dy=COM[1]/(1.*NPhi)-ws[i][1];
                    //SLOW LINE
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    out*=temp;
                }
                if (type=="laughlin") out*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
                else if (type=="CFL") out*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
            }
        }
        else {
            if (type=="laughlin-hole") {
                cout<<" 'trace' is not set up for laughlin-hole"<<endl;
                exit(0);
            }
            complex<double> w_comp0=w_comp;
            complex<double> sum=0., product;
            for (int k=0; k<invNu; k++) {
                product=1.;
                for( int i=0;i<invNu;i++){
                    dx=COM[0]/(1.*NPhi)-ws[i][0];
                    dy=COM[1]/(1.*NPhi)-ws[i][1]-k/(1.*invNu);
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    w_comp=w_comp0+1.*k*L2;
                    product*=temp;
                }
                if (type=="laughlin") product*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
                else if (type=="CFL") product*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
                if(k==1) sum+=(1.*trace)*product;
                else sum+=product;
            }
            out*=sum;
        }
        //  cout<<" com piece = "<<out<<endl;
        
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
                    //SLOW LINE
                    M(i,j)=product*polar(pow(in_determinant_rescaling, Ne-1), 2*M_PI*NPhi*(zs[i][1]*ds[j][0] - zs[i][0]*ds[j][1])/(2.*invNu*NPhi*Ne) );
                }
            }
            
            detSolver.compute(M);
            
            complex<double> tmp=detSolver.determinant();
            out*=tmp;
        }
        
    }
    return out;
}
complex<double> LATTICE::get_laughlinwf(vector<vector<double> > z){
    complex<double> out=1.,temp;
    double x,y;
    
    if (z.size()!=Ne) {
        cout<<"z.size()!=Ne"<<endl;
        exit(0);
    }
    
    //vandermonde piece
    int vandermonde_exponent=invNu;
    for( int i=0;i<Ne;i++){
        for( int j=i+1;j<Ne;j++){
            x=z[j][0]-z[i][0];
            y=z[j][1]-z[i][1];
            z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
            out*=pow(temp, vandermonde_exponent);
            out*=pow(in_determinant_rescaling, Ne-1);
        }
    }
    
    //COM piece (together with phase)
    double COM[2]={0.,0.};
    for( int i=0;i<Ne;i++){
        COM[0]+=z[i][0];
        COM[1]+=z[i][1];
    }
    double dx,dy;
    vector<double> wsum(2); for (int i=0; i<invNu; i++) for (int j=0; j<2; j++) wsum[j]+=ws[i][j];
    complex<double> w_comp = wsum[0]*L1+wsum[1]*L2, zcom_comp = 1.*COM[0]*L1+1.*COM[1]*L2;
    for( int i=0;i<invNu;i++){
        dx=COM[0]-ws[i][0];
        dy=COM[1]-ws[i][1];
        z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
        out*=temp;
    }
    
    out*=exp(1./(2.*NPhi)*( conj(w_comp)*zcom_comp - (w_comp)*conj(zcom_comp) ));
    
    return out;
}
void LATTICE::sum_locs(int out[]){
	out[0]=0; out[1]=0;
	for( int i=0;i<Ne;i++){
		out[0]+=locs[i][0];
		out[1]+=locs[i][1];
	}
}
vector< vector<int> > LATTICE::get_locs(){return locs;}
vector< vector<int> > LATTICE::get_ds(){return ds;}

//Coulomb Energy/ Pair-Amplitude Calculator
void LATTICE::setup_coulomb0(){
	coulomb_table=vector<vector<double> >(NPhi,vector<double>(NPhi,0));
	//coulomb_table[0][0]=v_coulomb_(&NPhi,&NPhi,&NPhi,&L1,&L2);
    complex<double> eL1=L1;
    complex<double> eL2=L2;
	for(int m=0;m<NPhi;m++){
		for(int n=0;n<NPhi;n++){
			coulomb_table[m][n]=new_v_coulomb_(&NPhi,&m,&n,&eL1,&eL2);
		}
	}
}
//void LATTICE::setup_table(int n, double cutoff, vector<vector<double>>& qtable, vector<vector<double>>& output, string typee, double screen){
//    if (typee!="ce" && typee!="pa") {
//        cout<<"only coulomb energy and pair-amplitude are calculable."<<endl;
//        exit(0);
//    }
//    
//    //boundary conditions.
//    complex<double> ph1, ph2;
//    ph1=polar(1., +2.*M_PI*get_bc()[1]);
//    ph2=polar(1., -2.*M_PI*get_bc()[0]);
//    //    cout<<"***** wsum *****"<<endl;
//    //    cout<<wsum[0]<<" "<<wsum[1]<<endl;
//    //    cout<<"**********"<<endl<<endl;
//    
//    qtable=vector<vector<double>>(NPhi,vector<double>(NPhi,0.));
//    Ftable=vector<vector<complex<double>>>(NPhi,vector<complex<double>>(NPhi,0.));
//    
//    int round=5;
//    for (int qx=0; qx<round*NPhi; qx++) {
//        for (int qy=0; qy<round*NPhi; qy++) {
//            int Qx=qx, Qy=qy;
//            if (2*Qx>round*NPhi) Qx-=round*NPhi;
//            if (2*Qy>round*NPhi) Qy-=round*NPhi;
//            double q2=2.*norm(Qx/(1.*NPhi)*L1+Qy/(1.*NPhi)*L2);
//            
//            int phy_qx=qx%NPhi, phy_qy=qy%NPhi;
//            if (2*phy_qx>NPhi) phy_qx-=NPhi;
//            if (2*phy_qy>NPhi) phy_qy-=NPhi;
//            int ind_m = (Qx-phy_qx)/NPhi;
//            int ind_n = (Qy-phy_qy)/NPhi;
//            
//            double sign = pow(-1., phy_qx*ind_n-phy_qy*ind_m)*pow(-1., NPhi*(ind_m*ind_n+ind_m+ind_n) );
//            Ftable[qx%NPhi][qy%NPhi] += sign * exp(-0.25*q2) * pow(ph1, ind_m) * pow(ph2, ind_n);
//            
//            if (typee=="ce")
//                qtable[qx%NPhi][qy%NPhi]+=1.*exp(-screen*sqrt(q2))/sqrt(q2)*pow(laguerre(n,0.5*q2),2)*exp(-0.5*q2)/(1.*NPhi);
//            else if (typee=="pa")
//                qtable[qx%NPhi][qy%NPhi]+=2.*laguerre(n, q2)*exp(-0.5*q2)/(1.*NPhi);
//            //NPhi comes from (2pi)/(2pi*Nphi), where numerator is from def of interaction (?), denormator is from sumq.
//        }
//    }
//    if (typee=="ce")
//        qtable[0][0]=0.;//Coulomb energy exclude (0,0) point. Replaced by Madelung energy.
//    
//    output=vector<vector<double>>(NPhi, vector<double>(NPhi,0.));
//    for (int qx=0; qx<NPhi; qx++)
//        for (int qy=0; qy<NPhi; qy++)
//            for (int i=0; i<NPhi; i++)
//                for (int j=0; j<NPhi; j++) {
//                    
//                    int Qx=qx,Qy=qy;
//                    if (2*Qx>NPhi) Qx-=NPhi;
//                    if (2*Qy>NPhi) Qy-=NPhi;
//                    complex<double> z=Qx/(1.*NPhi)*L1+Qy/(1.*NPhi)*L2;
//                    double x=sqrt(2.)*abs(z);
//                    
//                    double epsilon=1e-30;
//                    if (norm(Ftable[qx][qy])<epsilon) {
//                        continue;//exclude points that [f0]_N=0.
//                    }
//                    else if (cutoff<0. || x<cutoff) {
//                        //if cutoff is less than zero, no BZ cutoff.
//                        output[i][j]+=qtable[qx][qy]/norm(Ftable[qx][qy]) * cos( (2.*M_PI)/(1.*NPhi)*(qx*j-qy*i) );
//                    }
//                }
//    
//    //if calculating the LLL Coulomb energy, use Duncan's table;
//    if (typee=="ce" && n==0) {
//        //output=coulomb_table;
//        //output[0][0]=0.;
//    }
//}
void LATTICE::setup_table(NQ NQCE, vector<vector<vector<double>>>& qtable, vector<vector<vector<double>>>& output, string typee){
    if (typee!="ce" && typee!="pa") {
        cout<<"only coulomb energy and pair-amplitude are calculable."<<endl;
        exit(0);
    }
    
    int n=NQCE.N;
    double cutoff=NQCE.Q;
    vector<double> screen=NQCE.screen;
    
    //boundary conditions.
    complex<double> ph1, ph2;
    ph1=polar(1., +2.*M_PI*get_bc()[1]);
    ph2=polar(1., -2.*M_PI*get_bc()[0]);
    //cout<<"***** wsum *****"<<endl;
    //cout<<wsum[0]<<" "<<wsum[1]<<endl;
    //cout<<"**********"<<endl<<endl;
    
    qtable=vector<vector<vector<double>>>(screen.size(),vector<vector<double>>(NPhi,vector<double>(NPhi,0.)));
    Ftable=vector<vector<complex<double>>>(NPhi,vector<complex<double>>(NPhi,0.));
    
    int round=5;
    for (int qx=0; qx<round*NPhi; qx++) {
        for (int qy=0; qy<round*NPhi; qy++) {
            int Qx=qx, Qy=qy;
            if (2*Qx>round*NPhi) Qx-=round*NPhi;
            if (2*Qy>round*NPhi) Qy-=round*NPhi;
            double q2=2.*norm(Qx/(1.*NPhi)*L1+Qy/(1.*NPhi)*L2);
            
            int phy_qx=qx%NPhi, phy_qy=qy%NPhi;
            if (2*phy_qx>NPhi) phy_qx-=NPhi;
            if (2*phy_qy>NPhi) phy_qy-=NPhi;
            int ind_m = (Qx-phy_qx)/NPhi;
            int ind_n = (Qy-phy_qy)/NPhi;
            
            double sign = pow(-1., phy_qx*ind_n-phy_qy*ind_m)*pow(-1., NPhi*(ind_m*ind_n+ind_m+ind_n) );
            Ftable[qx%NPhi][qy%NPhi] += sign * exp(-0.25*q2) * pow(ph1, ind_m) * pow(ph2, ind_n);
            
            for (int s=0; s<screen.size(); s++) {
                if (typee=="ce")
                    qtable[s][qx%NPhi][qy%NPhi]+=1.*exp(-screen[s]*sqrt(q2))/sqrt(q2)*pow(laguerre(n,0.5*q2),2)*exp(-0.5*q2)/(1.*NPhi);
                else if (typee=="pa")
                    qtable[s][qx%NPhi][qy%NPhi]+=2.*laguerre(n, q2)*exp(-0.5*q2)/(1.*NPhi);
                //NPhi comes from (2pi)/(2pi*Nphi), where numerator is from def of interaction (?), denormator is from sumq.
            }

        }
    }
    if (typee=="ce") {
        for (int s=0; s<screen.size(); s++)
            qtable[s][0][0]=0.;//Coulomb energy exclude (0,0) point. Replaced by Madelung energy.
    }
    
    output=vector<vector<vector<double>>>(screen.size(),vector<vector<double>>(NPhi,vector<double>(NPhi,0.)));
    for (int qx=0; qx<NPhi; qx++)
        for (int qy=0; qy<NPhi; qy++)
            for (int i=0; i<NPhi; i++)
                for (int j=0; j<NPhi; j++) {
                    
                    int Qx=qx,Qy=qy;
                    if (2*Qx>NPhi) Qx-=NPhi;
                    if (2*Qy>NPhi) Qy-=NPhi;
                    complex<double> z=Qx/(1.*NPhi)*L1+Qy/(1.*NPhi)*L2;
                    double x=sqrt(2.)*abs(z);
                    
                    double epsilon=1e-30;
                    if (norm(Ftable[qx][qy])<epsilon) {
                        continue;//exclude points that [f0]_N=0.
                    }
                    else if (cutoff<0. || x<cutoff) {
                        for (int s=0; s<screen.size(); s++) {
                            //if cutoff is less than zero, no BZ cutoff.
                            output[s][i][j]+=qtable[s][qx][qy]/norm(Ftable[qx][qy]) * cos( (2.*M_PI)/(1.*NPhi)*(qx*j-qy*i) );
                        }
                    }
                }
    
    //if calculating the LLL Coulomb energy, use Duncan's table;
    if (typee=="ce" && n==0) {
        //output=coulomb_table;
        //output[0][0]=0.;
    }
}
void LATTICE::setup_tables(vector<NQ> N_Q, string typee) {
    
    if (typee=="ce") {
        coulomb_tableHLL.clear();//[NQ_ind][Screen_ind][3][4]
        coulomb_qtable.clear();
        CE_cutoff.clear();
        for (int i=0; i<N_Q.size(); i++) CE_cutoff.push_back(N_Q[i].Q);
    }
    else if (typee=="pa") {
        PA_table.clear();
        PA_qtable.clear();
        PA_cutoff.clear();
        for (int i=0; i<N_Q.size(); i++) PA_cutoff.push_back(N_Q[i].Q);
    }
    
    for (int i=0; i<N_Q.size(); i++) {
        vector<vector<vector<double>>> qtable, table;
        setup_table(N_Q[i], qtable, table, "ce");
        if (typee=="ce") {
            coulomb_qtable.push_back(qtable);
            coulomb_tableHLL.push_back(table);
        }
        else if (typee=="pa"){
            setup_table(N_Q[i], qtable, table, "pa");
            PA_qtable.push_back(qtable);
            PA_table.push_back(table);
        }
    }
}
vector<double> LATTICE::coulomb_energy(int ind, string typee){
    vector<double>out(coulomb_tableHLL[ind].size());
    //double out=0.;
    int m,n;
    for(int i=0;i<Ne;i++){
        for(int j=0;j<i;j++){
            m=supermod(locs[i][0]-locs[j][0], NPhi);
            n=supermod(locs[i][1]-locs[j][1], NPhi);
            
            for (int s=0; s<coulomb_tableHLL[ind].size(); s++) {
                if (typee=="ce")
                    out[s]+=coulomb_tableHLL[ind][s][m][n];
                else if (typee=="pa")
                    out[s]+=PA_table[ind][s][m][n];
            }
        }
    }
    if (typee=="ce") {
        for (int s=0; s<coulomb_tableHLL[ind].size(); s++) {
            out[s]+=0.5*Ne*coulomb_table[0][0];
        }
    }
    return out;
}
Eigen::VectorXd LATTICE::coulomb_energy_eigen(int ind, string typee){
    vector<double>out=coulomb_energy(ind, typee);
    Eigen::VectorXd output=Eigen::VectorXd::Zero(out.size());
    for (int i=0; i<out.size(); i++) output(i)=out[i];
    return output;
}
void LATTICE::shortrange(int ind, vector<double>& value, vector<double>& error, string typee){
    value=vector<double>(coulomb_tableHLL[ind].size());
    error=vector<double>(coulomb_tableHLL[ind].size());
    
    double Smax, nu=1./(1.*invNu), Sinf=nu*(1.-nu);
    if (type=="CFL") Smax=0.8;
    else Smax=0.5;
    double bound=abs(Smax-Sinf);
    if (Sinf<bound) bound=Sinf;
    
    for (int qx=0; qx<NPhi; qx++)
        for (int qy=0; qy<NPhi; qy++) {
            
            int Qx=qx, Qy=qy;
            if (2*Qx>NPhi) Qx-=NPhi;
            if (2*Qy>NPhi) Qy-=NPhi;
            
            complex<double> z=Qx/(1.*NPhi)*L1+Qy/(1.*NPhi)*L2;
            double x=sqrt(2.)*abs(z);
            
            for (int s=0; s<coulomb_tableHLL[ind].size(); s++) {
                if (typee=="ce" && x>CE_cutoff[ind]) {
                    value[s]+=-0.5*NPhi/pow(1.*invNu,2)*coulomb_qtable[ind][s][qx][qy];
                    error[s]+= 0.5*NPhi*abs(coulomb_qtable[ind][s][qx][qy])*bound;//NPhi comes from structure factor. From exp(ir)->S(q), there is a factor of NPhi. 
                }
                else if (typee=="pa" && x>PA_cutoff[ind]){
                    value[s]+=-0.5*NPhi/pow(1.*invNu,2)*PA_qtable[ind][s][qx][qy];
                    error[s]+= 0.5*NPhi*abs(PA_qtable[ind][s][qx][qy])*bound;
                }
            }
        }
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
    //vector<vector<complex<double>>>temp2(NPhi,vector<complex<double>>(NPhi,0.));
    for(int qx=0;qx<NPhi;qx++){
        for(int qy=0; qy<NPhi; qy++){
            for(int i=0;i<Ne;i++){
                temp[qx][qy]+=omega[supermod(2*(qx*locs[i][1]-qy*locs[i][0]),2*NPhi)];
                //Using this definition, 'rho(q)' is taken as e^{iq\times z} where q has spatial index.
                //temp2[qx][qy]+=omega[supermod(2*(qx*locs[i][1]+qy*locs[i][0]),2*NPhi)];
                //temp2 restore qy<=0 components.
            }
            sq[qx][qy]+=temp[qx][qy];
            sq2[qx][qy]+=norm(temp[qx][qy]);
            //sq_mqy[qx][qy]+=temp2[qx][qy];
            //sq2_mqy[qx][qy]+=norm(temp2[qx][qy]);
        }
    }
    
    int qx3,qy3;
    for(int qx1=0;qx1<NPhi;qx1++){
        for(int qx2=0;qx2<NPhi;qx2++){
            for(int qy1=0;qy1<NPhi;qy1++){
                for(int qy2=0;qy2<NPhi;qy2++){
                    qx3=supermod(-qx1-qx2,NPhi);
                    qy3=supermod(-qy1-qy2,NPhi);
                    sq3[qx1][qx2][qy1][qy2]+=temp[qx1][qy1]*temp[qx2][qy2]*temp[qx3][qy3];
                }
            }
        }
    }
}
void LATTICE::print_structure_factors(int nMeas, string filename){
    ofstream sqout("SQ/sq"+filename), sqout2("SQ/sq2"+filename), sqout3("SQ/sq3"+filename);
    for (int qx=0; qx<NPhi; qx++) {
        for (int qy=0; qy<NPhi; qy++) {
            sqout2<<sq2[qx][qy]/(1.*nMeas)<<" ";
            sqout<<abs(sq[qx][qy])/(1.*nMeas)<<" ";
        }
        sqout<<endl;
        sqout2<<endl;
    }
    for(int qx1=0; qx1<NPhi; qx1++) {
        for(int qy1=0; qy1<NPhi; qy1++) {
            for(int qx2=0; qx2<NPhi; qx2++) {
                for(int qy2=0; qy2<NPhi; qy2++) {
                    sqout3<<real(sq3[qx1][qx2][qy1][qy2])/(1.*nMeas)<<" ";
                }
                sqout3<<endl;
            }
        }
    }
    sqout.close();
    sqout2.close();
    sqout3.close();
    
//    ofstream sqout(type+"sq/sq"+filename), sqout2(type+"sq/sq2"+filename), sqout3(type+"sq/sq3"+filename);
//    //    ofstream sqoutfl(type+"sq/sqfl"+filename), sqout2fl(type+"sq/sq2fl"+filename);//'fl'=first line.
//    for(int qx=0;qx<NPhi;qx++){
//        for(int qy=0; qy<NPhi; qy++){
//            sqout2<<sq2[qx][qy]/(1.*nMeas)<<" ";
//            sqout<<abs(sq[qx][qy]/(1.*nMeas))<<" ";
//        }
//        sqout<<endl;
//        sqout2<<endl;
//    }
//    ofstream sqout_mqy(type+"sq/sq_mqy"+filename), sqout2_mqy(type+"sq/sq2_mqy"+filename), sqout3_mqy(type+"sq/sq3_mqy"+filename);
//    for(int qx=0;qx<NPhi;qx++){
//        for(int qy=0; qy<NPhi; qy++){
//            sqout2_mqy<<sq2_mqy[qx][qy]/(1.*nMeas)<<" ";
//            sqout_mqy<<abs(sq_mqy[qx][qy]/(1.*nMeas))<<" ";
//        }
//        sqout_mqy<<endl;
//        sqout2_mqy<<endl;
//    }
    
    
//    ofstream smaout(type+"sq/sma"+filename);
//    for (int qx=0; qx<NPhi; qx++) {
//        for (int qy=0; qy<NPhi; qy++) {
//            smaout<<SMAq[qx][qy]/(1.*nMeas)<<" ";
//        }
//        smaout<<endl;
//    }
    
//    	for(int qx1=0;qx1<NPhi;qx1++){
//    		for(int qy1=0;qy1<NPhi;qy1++){
//    			for(int qx2=0;qx2<NPhi;qx2++){
//    				for(int qy2=0;qy2<NPhi;qy2++)
//    					sqout3<<real(sq3[qx1][qx2][qy1][qy2])/(1.*nMeas)<<" ";
//    				sqout3<<endl;
//    			}
//    		}
//    	}
    //    sqoutfl.close();
    //    sqout2fl.close();
//    sqout3.close();
//    sqout3_mqy.close();
//    sqout.close();
//    sqout2.close();
//    sqout_mqy.close();
//    sqout2_mqy.close();
//    smaout.close();
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
    
    if (type=="laughlin" || type=="laughlin-hole") {
        locs=hot_start(locs.size(), NPhi, ran);
        cold_start();//TODO: when calculating ne=3, invnu=5, it seems that we have to start from cold_start(). otherwise, could't easily found good configuration.
    }
    else {
        locs=hot_start(locs.size(), NPhi, ran);
        //cold_start();
    }
    
	//**** setting up the initial determinant matrix
	running_weight=-1e11;
    c_running_weight=-1e11;
	int site=0;
	int initial_state_counter=0;
    
    while( (!correlatedsampling && running_weight<-1e10) || (correlatedsampling && c_running_weight<-1e10) ){
	//for some sizes the configuration specified by cold_start has zero weight
	//if that happens fiddle around until you find a better configuration, thats why theres a while loop
        
        if (site==Ne-1) locs[site][0]=locs[0][0];
        else locs[site][0]=locs[site+1][0];
        site++;

        if(site==Ne) {
            cout<<"couldn't easily find a good configuration!"<<endl;
            exit(0);
        }
		
        //initialize the Matrix/ Determinant for the ground state.
        if(type=="CFL") init_matdet();
        
        //initialize running_weight for the ground state. If correlated-sampling=true, excited both c_running_weight and running_weight for the excited states are initialized too (in function get_correlated_weight).
        //get_correlated_weight does a lot for excited states: synthocize locs, initialize new/old Matrix/Det, and the running_weight. Returns c_running_weight.
        running_weight=get_weight(locs);
        if (correlatedsampling) {
            for (int i=0; i<nonsamplestates.size(); i++) {
                nonsamplestates[i].change_locs(locs);
                nonsamplestates[i].init_matdet();
                nonsamplestates[i].running_weight=nonsamplestates[i].get_weight(locs);
            }
            c_running_weight=get_correlated_weight();
        }
        
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
        }else if(type=="doubledCFL"){
        	oldDeterminant=get_wf(locs);
        }
        
		initial_state_counter++;
		if(initial_state_counter>5000){
			cout<<"couldn't find a good starting configuration"<<endl;
			exit(0);
		}
	}

	//cout<<"starting weight"<<running_weight<<endl;
    if (testing && !correlatedsampling)
        cout<<"starting weight "<<running_weight<<" "<<get_weight(locs)<<" "<<log(norm(get_wf(locs)))<<endl;
    else if (testing && correlatedsampling) {
        cout<<c_running_weight<<" "<<get_correlated_weight()<<endl;
        cout<<"gs starting weight "<<running_weight<<" "<<get_weight(locs)<<" "<<log(norm(get_wf(locs)))<<endl;
        cout<<"es starting weight ";
        for (int i=0; i<nonsamplestates.size(); i++) {
            cout<<nonsamplestates[i].running_weight<<" "<<nonsamplestates[i].get_weight(locs)<<" "<<log(norm(nonsamplestates[i].get_wf(locs)))<<endl;
        }
    }

    //sq=vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi,0));
    //sq2=vector<vector<double>>(NPhi, vector<double>(NPhi,0));
    //sq_mqy=vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi,0));
    //sq2_mqy=vector<vector<double>>(NPhi, vector<double>(NPhi,0));
    check_sanity();
}
void LATTICE::init_matdet(){
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
    newMatrix=oldMatrix;
    oldDeterminant=detSolver.determinant();
    newDeterminant=oldDeterminant;
}
//same as above, only sets up the starting determinants
void LATTICE::reset(const vector< vector<int> > &zs){
	if(type=="CFL"){
		complex<double> temp,product;
		double x,y;
		oldMatrix=Eigen::MatrixXcd(Ne,Ne);
		for(int i=0;i<Ne;i++){
			for(int j=0;j<Ne;j++){
				product=1;
				for(int k=0;k<Ne;k++){
					if(k==i) continue;
		            x=det_helper(zs[i][0],zs[k][0],ds[j][0],dbar_parameter[0])/(1.*NPhi);
		            y=det_helper(zs[i][1],zs[k][1],ds[j][1],dbar_parameter[1])/(1.*NPhi);
					z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
					product*=temp;
				}
				oldMatrix(i,j)=product*polar(pow(in_determinant_rescaling,Ne-1), M_PI*(zs[i][1]*ds[j][0] - zs[i][0]*ds[j][1])/(1.*NPhi) );
			}
		}
		detSolver.compute(oldMatrix);
		newMatrix=oldMatrix;
		oldDeterminant=detSolver.determinant();
		newDeterminant=oldDeterminant;
	}
    check_sanity();
}
//checks a few different things to make sure that they make sense
void LATTICE::check_sanity(){
	if(type!="CFL" && type!="laughlin" && type!= "laughlin-hole" && type!="FilledLL" &&   type!="doubledCFL"){
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
    if(type=="FilledLL" && (zeros.size()!=(unsigned)NPhi || Ne!=NPhi)) {
        cout<<"Filled LL zeros size or Ne/NPhi wrong."<<endl;
        exit(0);
    }
    double pre=1e-15;
    double tmpx=fmod(get_shift()[0]-get_bc()[0],1.);
    double tmpy=fmod(get_shift()[1]-get_bc()[1],1.);
    if ( abs(tmpx)>pre or abs(tmpy)>pre ) {
        cout<<"shift!=bc."<<endl;
        cout<<"shift="<<get_shift()[0]<<" "<<get_shift()[1]<<" bc="<<get_bc()[0]<<" "<<get_bc()[1]<<endl;
        cout<<get_shift()[0]-get_bc()[0]<<" "<<get_shift()[1]-get_bc()[1]<<endl;
        exit(0);
    }
    
    //check boundary conditions for correlated sampling.
    if (correlatedsampling) {
        bool tmp=true;
        for (int i=0; i<nonsamplestates.size(); i++) {
            if ( abs(nonsamplestates[i].get_bc()[0]-shiftx)>1e-15 or abs(nonsamplestates[i].get_shift()[0]-shiftx)>1e-15 ) tmp=false;
            if ( abs(nonsamplestates[i].get_bc()[1]-shifty)>1e-15 or abs(nonsamplestates[i].get_shift()[1]-shifty)>1e-15 ) tmp=false;
        }
        if (!tmp) {
            cout<<"Excited states' shift != B.C."<<endl;
            cout<<"shift="<<shiftx<<" "<<shifty<<endl;
            for (int i=0; i<nonsamplestates.size(); i++) {
                cout<<nonsamplestates[i].get_shift()[0]<<" "<<nonsamplestates[i].get_shift()[1]<<" "<<nonsamplestates[i].get_bc()[0]<<" "<<nonsamplestates[i].get_bc()[1]<<endl;
            }
            exit(0);
        }
        if (ratio.size()!=nonsamplestates.size()+1) {
            cout<<"ratio.size() wrong."<<endl;
            exit(0);
        }
    }
}
//changes the dbar (which ordinarily is the sum of d's) to whatever we want
void LATTICE::change_dbar_parameter(double dbarx, double dbary){
	if(type!="CFL" && type!="FilledLL" && type!="doubledCFL"){
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
	omega_dbar_range=2;
	omega_dbar=vector<vector<complex<double>>>(2,vector<complex<double>>(2*omega_dbar_range+1));
	for(int i=-omega_dbar_range; i<=omega_dbar_range; i++)
		for(int j=0; j<2; j++) omega_dbar[j][i+omega_dbar_range]=polar(1.,M_PI*dbar_parameter[j]*i/(1.*NPhi));
	
}
void LATTICE::set_ds(vector< vector<int> > tds){
    if (tds.size()!=(unsigned)Ne) {
        cout<<"cannot set ds, ds.size() wrong."<<endl;
        exit(0);
    }
    ds=tds; dsum=vector<int>(2,0);
	for(int i=0;i<Ne;i++){
        dsum[0]+=ds[i][0]*invNu; dsum[1]+=ds[i][1]*invNu; //'ds' is on L/Ne lattice, 'dsum' is on L/Nphi lattice.
	}
	print_ds();
	change_dbar_parameter(dsum[0]/(1.*Ne)+real(dbar_delta),dsum[1]/(1.*Ne)+imag(dbar_delta));
    
    //reset ws, according to ds.
    ws.clear(); ws=ws0;
    if (type=="CFL" or type=="doubledCFL") {
        for (unsigned int i=0; i<ws.size(); i++) {
            ws[i][0]+=dsum[0]/(1.*invNu*NPhi); ws[i][1]+=dsum[1]/(1.*invNu*NPhi);// if 'd = wsum'.
        }
    }
    
    //shift ws.
    for (int i=0; i<invNu; i++) {
        ws[i][0]+=shiftx/(1.*invNu);
        ws[i][1]+=shifty/(1.*invNu);
    }
    wsum=vector<double>(2,0.);
    for (int i=0; i<invNu; i++) {
        wsum[0]+=ws[i][0];
        wsum[1]+=ws[i][1];
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

complex<double> LATTICE::getL(int dir){
	if(dir==1) return L1;
	else return L2;
}
void LATTICE::cold_start(){
	for(int i=0;i<(signed)locs.size();i++){
		locs[i][0]=i;
		locs[i][1]=i;
	}
}
vector< vector<int> > LATTICE::hot_start(int Ne_t, int NPhi_t, MTRand &ran_t){
	unordered_set<int> locs_table;
	bool found;
	int x,y;
	vector< vector<int> > new_locs(Ne_t, vector<int>(2));
	for(int i=0;i<(signed)new_locs.size();i++){
		found=true;
		while(found){
			x=ran_t.randInt(NPhi_t-1);
			y=ran_t.randInt(NPhi_t-1);
			if(locs_table.count(NPhi_t*x+y)==0) found=false;
		}
		locs_table.insert(NPhi_t*x+y);
		new_locs[i][0]=x;
		new_locs[i][1]=y;
		//cout<<x<<" "<<y<<endl;
	}
	return new_locs;
}

//call's duncan's lattice_z function, if the arguments x or y are outside of the range (0,NPhi) it shifts them into that range
//and multiplies by the appropriate phase
//TODO:only works for a square torus
complex<double> LATTICE::modded_lattice_z(int x, int y){
	int modx=supermod(x,NPhi);
	int mody=supermod(y,NPhi);
//	complex<double> out=lattice_z_(&NPhi,&modx,&mody,&L1,&L2,&one);
	complex<double> out=shifted_ztable[modx][mody];
	int j=(modx-x)/NPhi, k=(mody-y)/NPhi;
	
	//the reason this isn't tabulated is because dbar_parameter doesn't have to live on a lattice
	//but could make a new table, it would make the code slightly faster
//	out*=omega_dbar[supermod(-(mody)*j+(modx)*k,2*NPhi)];
	//SLOW LINE
	out*=omega[supermod(-(mody)*j+(modx)*k,2*NPhi)]*omega_dbar[0][k+omega_dbar_range]*omega_dbar[1][-j+omega_dbar_range];
//	out*=polar(1., (-(mody+dbar_parameter[1])*j+(modx+dbar_parameter[0])*k)*M_PI/(1.*NPhi));
	if(j%2 || k%2) return -out;
	else return out;
}

//similar to above but for the COM part
complex<double> LATTICE::modded_lattice_z_nodbar(int x, int y){
	int modx=supermod(x,NPhi);
	int mody=supermod(y,NPhi);
	complex<double> out=lattice_z_(&NPhi,&modx,&mody,&L1,&L2,&one);
//	complex<double> out=shifted_ztable[modx][mody];
	int j=(modx-x)/NPhi, k=(mody-y)/NPhi;
	
	//the reason this isn't tabulated is because dbar_parameter doesn't have to live on a lattice
	//but could make a new table, it would make the code slightly faster
	out*=-omega[supermod(-(mody)*j+(modx)*k,2*NPhi)];

	if(j%2 || k%2) return -out;
	else return out;
}
void LATTICE::make_CFL_det(Eigen::MatrixXcd& new_Matrix, vector<int> newloc, int electron, complex<double>& new_det, const vector<vector<int>>& zs){
    complex<double> product, temp;
    vector<int> z(2);
    
    for(int i=0;i<Ne;i++){
        for(int j=0;j<Ne;j++){
            if(i==electron){
                product=1.;
                for(int k=0;k<Ne;k++){
                    if(k==i) continue;
                    det_helper(newloc,zs[k],ds[j],z);
                    product*=modded_lattice_z(z[0],z[1]);
                }
                new_Matrix(i,j)=pow(in_determinant_rescaling,Ne-1)*product*omega[supermod(newloc[1]*ds[j][0] - newloc[0]*ds[j][1],2*NPhi)];
            }
            else if(i!=electron){//all other elements just need to be updated by the ratio of a sigma function
                det_helper(zs[i],zs[electron],ds[j],z);
                temp=modded_lattice_z(z[0],z[1]);
                if( temp==0. || new_Matrix(i,j)==0. || Ne<=2){
                    product=1.;
                    for(int k=0;k<Ne;k++){
                        if(k==i) continue;
                        if(k==electron)
                            det_helper(zs[i],newloc,ds[j],z);
                        else
                            det_helper(zs[i],zs[k],ds[j],z);
                        product*=modded_lattice_z(z[0],z[1]);
                    }
                    new_Matrix(i,j)=pow(in_determinant_rescaling,Ne-1)*product*omega[supermod(zs[i][1]*ds[j][0] - zs[i][0]*ds[j][1],2*NPhi) ];
                }
                else{
                    new_Matrix(i,j)/=temp;
					det_helper(zs[i],newloc,ds[j],z);
                    new_Matrix(i,j)*=modded_lattice_z(z[0],z[1]);
                }
            }
        }
    }
    detSolver.compute(new_Matrix);
    new_det=detSolver.determinant();
}
vector<double> LATTICE::get_dbar_parameter(){
    return dbar_parameter;
}
void LATTICE::setup_landautable(){
    landautable=vector<vector<vector<complex<double>>>>(NPhi, vector<vector<complex<double>>>(NPhi, vector<complex<double>>(NPhi, 1.)));
    //landautable[m][x][y], m landau state, (x, y) lattice point.
    
    vector<vector<vector<double>>> ZEROS(NPhi, vector<vector<double>>(NPhi, vector<double>(2, 0.)));
    //ZEROS[m][ind][i]: m-state index, ind-zero index, i-x or y.
    for (int ind_state=0; ind_state<NPhi; ind_state++) {
        for (int ind_zero=0; ind_zero<NPhi; ind_zero++) {
            
            ZEROS[ind_state][ind_zero][0]=(ind_zero+0.5)/(1.*NPhi)-0.5 + get_bc()[0]/(1.*NPhi);
            ZEROS[ind_state][ind_zero][1]=ind_state/(1.*NPhi) + get_bc()[1]/(1.*NPhi);
            
        }
    }
    
    vector<complex<double>> AVEZERO(NPhi, 0.);
    for (int ind_state=0; ind_state<NPhi; ind_state++) {
        vector<double> tmp(2, 0.);
        for (int ind_zero=0; ind_zero<NPhi; ind_zero++) {
            tmp[0]+=ZEROS[ind_state][ind_zero][0];
            tmp[1]+=ZEROS[ind_state][ind_zero][1];
        }
        AVEZERO[ind_state] = (tmp[0]*L1+tmp[1]*L2)/(1.*NPhi);
    }
 
    for (int ind_state=0; ind_state<NPhi; ind_state++)
        for (int x=0; x<NPhi; x++)
            for (int y=0; y<NPhi; y++) {
                
                for (int ind=0; ind<NPhi; ind++) {
                    double dx=x/(1.*NPhi)-ZEROS[ind_state][ind][0], dy=y/(1.*NPhi)-ZEROS[ind_state][ind][1];
                    complex<double> temp;
                    z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
                    landautable[ind_state][x][y]*=temp;
                }
                
                complex<double> z = x/(1.*NPhi)*L1+y/(1.*NPhi)*L2;
                complex<double> temp = conj(AVEZERO[ind_state])*z - AVEZERO[ind_state]*conj(z);
                landautable[ind_state][x][y]*=exp(0.5*temp);
                
            }
    
}
//given both a set of positions and a set of ds, computes the wavefunction (NOT the norm of the wavefunction)
complex<double> LATTICE::doubled_CFL(const vector< vector<int> > &bothzs){
    if (bothzs.size()!=(unsigned)Ne*2) {
        cout<<"cannot get_wf because zs.size()!=2Ne"<<endl;
        exit(0);
    }
	complex<double> out=1,temp;
    
	vector<vector<vector<int> > > onezs(2,bothzs);
	
	vector<vector<int> > zs;
	onezs[0].resize(Ne);
	onezs[1].erase(onezs[1].begin(),onezs[1].begin()+Ne);
	for(int i=0;i<Ne;i++){
		onezs[1][i][0]*=-1; //=supermod(-onezs[1][i][0],NPhi);
		onezs[1][i][1]*=-1; //=supermod(-onezs[1][i][1],NPhi);
	}
	for(int copy=0;copy<2;copy++){
		zs=onezs[copy];
//		for(int i=0;i<Ne;i++) cout<<"("<<zs[i][0]<<","<<zs[i][1]<<") ";
//		cout<<endl;
        //COM piece
        int COM[2]={0,0};
        for( int i=0;i<Ne;i++){
            COM[0]+=zs[i][0];
            COM[1]+=zs[i][1];
        }
        double dx,dy;
        for( int i=0;i<invNu;i++){
            dx=COM[0]/(1.*NPhi)-ws[i][0];
            dy=COM[1]/(1.*NPhi)-ws[i][1]-copy/(1.*invNu);
            z_function_(&dx,&dy,&L1,&L2,&zero,&NPhi,&temp);
            out*=temp;
        }
//        cout<<" com piece = "<<out<<endl; 
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
        
        //phase previously missing.
        vector<double> wsum(2);
        for (int i=0; i<invNu; i++) {
            wsum[0]+=ws[i][0];
            wsum[1]+=ws[i][1]+copy/(1.*invNu);
        }
        complex<double> w_comp = wsum[0]*L1+wsum[1]*L2, zcom_comp = 1.*COM[0]/NPhi*L1+1.*COM[1]/NPhi*L2;
        complex<double> dsum_comp = 1.*dsum[0]/NPhi*L1 + 1.*dsum[1]/NPhi*L2;
        out*=exp(1./(2.*NPhi)*( conj(w_comp - dsum_comp)*zcom_comp - (w_comp - dsum_comp)*conj(zcom_comp) ));
    }
    return out;
}
void LATTICE::update(){
    if(type=="CFL"){
		oldMatrix=newMatrix;
		oldDeterminant=newDeterminant;
	}
}
//for Correlated Sampling
vector<double> LATTICE::get_ratio(){return ratio;}
void LATTICE::change_locs(const vector<vector<int>> &locs_t){locs=locs_t;}
void LATTICE::update_matdet(bool bo){
    if (bo) {
        oldMatrix=newMatrix;
        oldDeterminant=newDeterminant;
    }
    else {
        newMatrix=oldMatrix;
        newDeterminant=oldDeterminant;
    }
}
void LATTICE::setup_nonsamplestates(vector<vector<vector<int>>> ds_t, vector<double> ratio_t){
    nonsamplestates=vector<LATTICE>(ds_t.size());
    for (int i=0; i<nonsamplestates.size(); i++) {
        bool tr=false, cor=false;//trace, correlatedsampling (can be true, does not affect result).
        int seed=0;//seed is not important.
        nonsamplestates[i]=LATTICE(Ne, invNu, 0, "CFL", seed, gs, theta, alpha, shiftx, shifty, tr, cor);
        if (ds_t[i].size()!=Ne) {cout<<"dset size wrong"<<endl; exit(0);}
        nonsamplestates[i].set_ds(ds_t[i]);
    }
    ratio=vector<double>(nonsamplestates.size()+1, 0.);
    ratio=ratio_t;//ratio is weight for the correlated sampling function.
}
double LATTICE::get_correlated_weight(){
    double out=ratio[nonsamplestates.size()]*exp(running_weight);
    double temp=1.;
    for (int i=0; i<nonsamplestates.size(); i++) {
        //out+=ratio[i]*exp(nonsamplestates[i].running_weight);
        
        temp+=ratio[i]/ratio[nonsamplestates.size()]*exp(nonsamplestates[i].running_weight-running_weight);
        //TODO:correlated-sampling
    }
    return log(ratio[nonsamplestates.size()])+running_weight+log(temp);
}
vector<double> LATTICE::get_runweis(){
    vector<double> out(nonsamplestates.size(), 0.);
    for (int i=0; i<nonsamplestates.size(); i++) {
        out[i]=nonsamplestates[i].running_weight;
    }
    return out;
}
vector<double> LATTICE::get_w(){
    //returns un-normalized w_i, i=0,1,...,Nn.
    int Nn=nonsamplestates.size();
    vector<double> out(Nn+1, 0.);
    
    vector<double> weit=get_runweis(); weit.push_back(running_weight);
    double den=0.;
    for (int i=0; i<Nn+1; i++) den+=ratio[i]*exp(weit[i]-weit[Nn]);
    for (int i=0; i<Nn+1; i++) out[i]=exp(weit[i]-weit[Nn])/den;
    return out;
}
LATTICE::~LATTICE(){
    if(testing)
        cout<<"acceptance rate: "<<accepts/(1.*tries)<<endl;
}
