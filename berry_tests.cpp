#include "berry_tests.h"

void hermitianize(Eigen::MatrixXcd &x){
	for(int i=0;i<x.rows();i++){
		for(int j=i+1;j<x.rows();j++){
			x(i,j)=0.5*(x(i,j)+conj(x(j,i)));
			x(j,i)=conj(x(i,j));
		}
	}
}
double phasemod(complex<double> in){
    double out=arg(in);
    if(out<0) return out+2*M_PI;
    else if (out>2*M_PI) return out-2*M_PI;
    else return out;
    return out;
}
inline double Laguerrel(int N, double x){
    vector<double> v;
    v.push_back(1); v.push_back(1-x);
    for (int i=1; i<N; i++) {
        double temp=((2.0*i+1-x)*v[i]-i*v[i-1])/(1.0+i);
        v.push_back(temp);
    }
    return v[N];
}


void single_run(string filename, bool trace){
    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    double theta_t, alpha_t;
    string type;
    
    ifstream infile(filename);
    infile>>Ne>>invNu>>theta_t>>alpha_t;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    cout<<"Ne="<<Ne<<" invNu="<<invNu<<" nMeas="<<nMeas<<" nSteps="<<nSteps<<endl;
    
    int gs=0;
    int NPhi=Ne*invNu;
    if (type=="laughlin-hole") NPhi++;

    double theta=theta_t*M_PI, alpha=alpha_t;
    
//    LATTICE ds_generator(Ne, invNu, testing, type, seed, gs, theta, alpha);
//    vector< vector<int> > old_ds=ds_generator.get_ds();
//    vector<vector<int>> old_ds;
//    if (Ne==21) {
//        old_ds.clear();
//        for (int i=-3; i<4; i++) {
//            for (int j=-1; j<2; j++) {
//                old_ds.push_back(vector<int>{i,j});
//            }
//        }
//    }
//    if (Ne==17) {
//        old_ds.clear();
//        for (int i=-2; i<3; i++) {
//            for (int j=0; j<2; j++) {
//                old_ds.push_back(vector<int>{i,j});
//            }
//        }
//        for (int i=-1; i<2; i++) {
//            old_ds.push_back(vector<int>{i,-1});
//            old_ds.push_back(vector<int>{i,2});
//        }
//        old_ds.push_back(vector<int>{2,-1});
//    }
    
    LATTICE ll(Ne, invNu, testing, type, seed, gs, theta, alpha, trace);
    ll.setup_coulomb0();
    ll.setup_tables(vector<int>{0,2}, "ce");
//    ll.set_ds(old_ds);
//    ll.print_ds();
    
    ofstream outfile("out"), eout("energy");

    for(int s=0;s<nBins;s++){
        //ll.change_dbar_parameter(s*0.1,s*0.1);
        ll.reset();
        
        ll.step(nWarmup);
        double E=0,E2=0;
        double EE=0, EE2=0;
        double P=0;
        double e,ee,p;
        complex<double> berry_phase(0,0);
        deque<double> e_tracker, p_tracker;
        int Ntrack=50;
        vector<double> autocorr_e(Ntrack,0), autocorr_p(Ntrack,0);
        
        for(int i=0;i<nMeas;i++){
            ll.step(nSteps);
            e=ll.coulomb_energy(0, "ce");
            ee=ll.coulomb_energy(1, "ce");
            E+=e;EE+=ee;
            E2+=e*e;EE2+=ee*ee;
            p=ll.running_weight;
            P+=p;
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
            //ll.update_structure_factors();
        }
        //ll.print_structure_factors(nMeas, "_"+to_string((long long int)s));
        
        outfile<<E/(1.*nMeas*ll.Ne)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))/(1.*ll.Ne)<<" "<<real(berry_phase)/(1.*nMeas)<<" "<<imag(berry_phase)/(1.*nMeas)<<endl;
        cout<<"acceptance rate: "<<(1.*ll.accepts)/(1.*ll.tries)<<endl;
        //while doing experiment on standard error, i found we should use the follows as error. (ed result for 4/12 is -0.414171)
        cout<<"n=0 Landau Level, coulomb1"<<endl;
        cout<<"E="<<E/(1.*nMeas*ll.Ne)<<" var="<<sqrt(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))/sqrt(1.*nMeas)/(1.*ll.Ne)<<endl;
        cout<<"n=0 Landau Level, coulomb2"<<endl;
        cout<<"E="<<EE/(1.*nMeas*ll.Ne)<<" var="<<sqrt(EE2/(1.*nMeas)-pow(EE/(1.*nMeas),2))/sqrt(1.*nMeas)/(1.*ll.Ne)<<endl;
        
        ofstream auto_out("auto");
        for(int j=0;j<Ntrack;j++){
            auto_out<<j+1<<" ";
            auto_out<<autocorr_e[j]/(1.*(nMeas-Ntrack))<<" "<<pow(E/(1.*nMeas),2)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))<<" ";
            //auto_out<<autocorr_p[j]/(1.*(nMeas-Ntrack))<<" "<<pow(P/(1.*nMeas),2)<<" "<<(P2/(1.*nMeas)-pow(P/(1.*nMeas),2))<<" ";
            auto_out<<endl;
        }
    }
    outfile<<endl;
    eout.close();
    outfile.close();
}
void parallel_ce_pa(int ncore, vector<int> CE, vector<int> PP, bool bo_shift, double shift, string filename){
    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    double theta_t, alpha_t;
    string type;
    
    ifstream infile(filename);
    infile>>Ne>>invNu>>theta_t>>alpha_t;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    cout<<"Ne="<<Ne<<" invNu="<<invNu<<" nMeas="<<nMeas<<" nSteps="<<nSteps<<" nBins="<<nBins<<" ncore="<<ncore<<endl;
    
    int gs=0;
    int NPhi=Ne*invNu;
    if (type=="laughlin-hole") NPhi++;
    
    double theta=theta_t*M_PI, alpha=alpha_t;
    
    vector<vector<int>> ds;
    if (Ne==8) {
        ds=vector<vector<int>> (8,vector<int>(2,0));
        for (int i=0; i<8; i++) {
            ds[i][0]=i/3-1;
            ds[i][1]=i%3-1;
        }
    }
    if (Ne==16) {
        ds.clear();
        for (int i=-2; i<=2; i++) {
            vector<int> temp1(2,0), temp2(2,0);
            temp1[0]=i; temp2[0]=i;
            temp1[1]=0; temp2[1]=-1;
            ds.push_back(temp1);
            ds.push_back(temp2);
        }
        for (int i=-1; i<=1; i++) {
            vector<int> temp1(2,0), temp2(2,0);
            temp1[0]=i; temp2[0]=i;
            temp1[1]=1; temp2[1]=-2;
            ds.push_back(temp1);
            ds.push_back(temp2);
        }
    }
//    for (int i=0; i<ds.size(); i++) {
//        cout<<ds[i][0]<<" "<<ds[i][1]<<endl;
//    }
//    exit(0);
    
    vector<LATTICE> ll(ncore);
    for (int i=0; i<ncore; i++) {
        seed=i;
        ll[i]=LATTICE(Ne, invNu, testing, type, seed, gs, theta, alpha, false);
        
        ll[i].set_lat_scalex(1);//force lattice sum on Nphi*Nphi lattice.
        ll[i].set_lat_scaleq(1);//force lattice sum on Nphi*Nphi lattice.
        
        if (!bo_shift) {
            if (NPhi%2==0)
                shift=0.25;
            else
                shift=0.;
        }
        
        ll[i].shift_ws(shift);
        
        if (Ne==8 && type=="CFL") ll[i].set_ds(ds);
        if (Ne==16&& type=="CFL") ll[i].set_ds(ds);
        
        
        ll[i].setup_tables(CE, "ce");
        ll[i].setup_tables(PP, "pa");
        
    }
    
    int Coul_type=CE.size();
    vector<vector<double>> E(Coul_type, vector<double>(nBins,0.)), EE(Coul_type, vector<double>(nBins,0.));
    
    int PA_type=PP.size();
    vector<vector<double>> PA(PA_type, vector<double>(nBins,0.)), PAPA(PA_type, vector<double>(nBins,0.));

    omp_set_num_threads(ncore);
#pragma omp parallel for
    for(int s=0;s<nBins;s++){
        int coren=omp_get_thread_num();
        ll[coren].reset();
        ll[coren].step(nWarmup);
        
        for(int i=0;i<nMeas;i++){
            ll[coren].step(nSteps);
            
            //Coulomb Energy
            for (int c=0; c<Coul_type; c++) {
                double e=ll[coren].coulomb_energy(c,"ce");
                E[c][s]+=e;
                EE[c][s]+=e*e;
            }
            
            //Pair Amplitude
            for (int p=0; p<PA_type; p++) {
                double pa=ll[coren].coulomb_energy(p,"pa");
                PA[p][s]+=pa;
                PAPA[p][s]+=pa*pa;
            }

        }
    }
    
    vector<double> Etotal(Coul_type,0.), EEtotal(Coul_type,0.);
    vector<double> PAtotal(PA_type,0.), PAPAtotal(PA_type,0.);
    
    for (int s=0; s<nBins; s++) {
        for (int i=0; i<Coul_type; i++) {
            Etotal[i]+=E[i][s];
            EEtotal[i]+=EE[i][s];
        }
        for (int p=0; p<PA_type; p++) {
            PAtotal[p]+=PA[p][s];
            PAPAtotal[p]+=PAPA[p][s];
        }
    }
    
    //while doing experiment on standard error, i found we should use the follows as error. (ed result for 4/12 is -0.414171)
    ofstream outfile("out_"+filename);
    ofstream outpa("out_pa_"+filename);
    outfile<<"Ne="<<Ne<<" invNu="<<invNu<<" nMeas="<<nMeas<<" nBins="<<nBins<<endl;
    outpa<<"Ne="<<Ne<<" invNu="<<invNu<<" nMeas="<<nMeas<<" nBins="<<nBins<<endl;
    outfile<<"shift="<<ll[0].get_shift()<<endl<<endl;
    outpa<<"shift="<<ll[0].get_shift()<<endl<<endl;
    
    
    nMeas*=nBins;
    
    for (int c=0; c<Coul_type; c++) {
        double short_value, short_error;
        ll[0].shortrange(c, short_value, short_error, "ce");
        short_value/=(1.*Ne);
        short_error/=(1.*Ne);
        double MCerror=sqrt(EEtotal[c]/(1.*nMeas)-pow(Etotal[c]/(1.*nMeas),2))/sqrt(1.*nMeas)/(1.*Ne);
        
        outfile<<"n="<<CE[c]<<" COULOMB-ENERGY."<<endl;
        outfile<<"cutoff="<<ll[0].CE_cutoff[c]<<endl;
        outfile<<"truncation value  ="<<setprecision(10)<<short_value<<endl;
        outfile<<"truncation error  = "<<setprecision(10)<<short_error<<endl;
        outfile<<"MonteCarlo error  = "<<setprecision(10)<<MCerror<<endl;
        
        outfile<<"E="<<setprecision(10)<<Etotal[c]/(1.*nMeas*Ne)+short_value<<" var="<<sqrt( MCerror*MCerror + short_error*short_error)<<endl<<endl;
    }
    for (int p=0; p<PA_type; p++) {
        double short_value, short_error;
        ll[0].shortrange(p, short_value, short_error, "pa");
        short_value/=(1.*Ne);
        short_error/=(1.*Ne);
        double MCerror=sqrt(PAPAtotal[p]/(1.*nMeas)-pow(PAtotal[p]/(1.*nMeas),2))/sqrt(1.*nMeas)/(1.*Ne);
        
        outpa<<"n="<<PP[p]<<" PAIR-AMPLITUDE."<<endl;
        outpa<<"cutoff="<<ll[0].PA_cutoff[p]<<endl;
        outpa<<"truncation value  ="<<setprecision(10)<<short_value<<endl;
        outpa<<"truncation error  = "<<setprecision(10)<<short_error<<endl;
        outpa<<"MonteCarlo error  = "<<setprecision(10)<<MCerror<<endl;
        
        outpa<<"PA="<<setprecision(10)<<PAtotal[p]/(1.*nMeas*Ne)+short_value<<" var="<<sqrt( MCerror*MCerror + short_error*short_error)<<endl<<endl;
    }
    outfile.close();
    outpa.close();
}
void structurefactor(string intputfilename, int num_core){//fielname='params_sq_...'.
    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    double theta_t, theta, alpha;
    string type;
    ifstream infile(intputfilename);
    infile>>Ne>>invNu>>theta_t>>alpha;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    theta=theta_t*M_PI;
    
    omp_set_num_threads(num_core);
    
    vector<LATTICE> ll(num_core);
    for (int i=0; i<num_core; i++) {
        ll[i]=LATTICE(Ne, invNu, testing, type, i, 0, theta, alpha);//gs=0, seed=i.
    }
    if (type=="CFL") {
        vector<vector<int>> ds(Ne, vector<int>(2));
        if (Ne==64 || Ne==49 || Ne==36 || Ne==25 || Ne==16 || Ne==9 || Ne==4) {
            int dis=(int)(sqrt(Ne));
            for (int y=0; y<dis; y++) {for (int x=0; x<dis; x++) {ds[x+y*dis][0]=x-dis/2; ds[x+y*dis][1]=y-dis/2;}}
            for (int i=0; i<num_core; i++) {
                ll[i].set_ds(ds);
            }
            ll[0].print_ds();
        }
    }
    
#pragma omp parallel for
    for(int s=0;s<nBins;s++){
        int coren=omp_get_thread_num();
//        printf("coren=%d\n",coren);
//        ll[coren].change_dbar_parameter(s*0.1,s*0.1);
        ll[coren].reset();
        ll[coren].step(nWarmup);

        for(int i=0;i<nMeas;i++){
            ll[coren].step(nSteps);
            ll[coren].update_structure_factors();
        }
        ll[coren].print_structure_factors(nMeas, intputfilename+"_"+to_string((long long int)(s)));
    }
}
//laughlin-hole berry phase.
void laughlinberryphase(string input_name, string output_name, vector<double> length, double steplength, int change_nMeas, int change_Ne, int num_core, double theta, double alpha){
    //input file is 'params_name'.
    //For some Ne, kind specifies certain type of berry phase loops.
    vector<data> datas;
    //output will be three files.
    //The first contains phase, amp, energy for each step each bin, and final berry matrix.
    //The second is same as the first, but in a format friend to Mathematica.
    //The third is matrix element in each step each bin.
    
//    ofstream outfile(output_name.c_str());
//    string outfile2name=output_name+"_Mathmatica";//output data in format easy for Mathematica to deal with.
//    ofstream outfile2(outfile2name.c_str());
//    string outfile3name=output_name+"_Matrix";//output matrix element in each step.
//    ofstream outfile3(outfile3name.c_str());
    
    int Ne,Ne_t,invNu,nWarmup,nMeas,nMeas_t,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile(input_name);
    infile>>Ne_t>>invNu;
    infile>>nWarmup>>nMeas_t>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    //if change_nMeas==0/ change_Ne==0, nMeas/ Ne is set according to 'params'.
    if (change_nMeas==0) nMeas=nMeas_t;
    else nMeas=change_nMeas;
    if (change_Ne==0) Ne=Ne_t;
    else Ne=change_Ne;
    
    vector<vector<double> > holes; vector<int> Grid(2);
    for (int i=0; i<2; i++) {Grid[i]=(int)(length[i]/steplength);}
    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=steplength*i; a[1]=0.;                  holes.push_back(a);}
    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=steplength*i; a[0]=length[0];           holes.push_back(a);}
    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=length[0]-steplength*i; a[1]=length[1]; holes.push_back(a);}
    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=length[1]-steplength*i; a[0]=0.;        holes.push_back(a);}
    int nds=holes.size();
    vector<vector<double> > holes2(nds, vector<double>(2,0));
    int supermod(int k, int n);
    for(int i=0;i<nds;i++) holes2[supermod(i-1,nds)]=holes[i];//(holes[b],holes2[b]) = (holes[b],holes[b+1]).
    
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
    vector<vector<Eigen::MatrixXcd > > overlaps(nds, vector<Eigen::MatrixXcd>(4, Eigen::MatrixXcd::Zero(invNu, invNu)));
    
    omp_set_num_threads(num_core);
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(invNu)), pp(num_core, vector<LATTICE>(invNu));//do this to avoid wrong memory access since openmp share memory.
    for (int k=0; k<num_core; k++) for (int i=0; i<invNu; i++) {ll[k][i]=LATTICE(Ne, invNu, testing, type, seed, i, theta, alpha); pp[k][i]=LATTICE(Ne, invNu, testing, type, seed, i, theta, alpha);}
    
    vector<double> trace_phase(nBins,0.);
    
    for (int nbin=0; nbin<nBins; nbin++) {
        //parallel programming begin.
#pragma omp parallel for
        for(int b=0; b<nds; b++) {
            
            int coren = omp_get_thread_num();
            for (int i=0; i<invNu; i++) {
                ll[coren][i].set_hole(holes[b]);
                pp[coren][i].set_hole(holes2[b]);
                ll[coren][i].reset(); ll[coren][i].step(nWarmup);
                pp[coren][i].reset();
            }
            
            for (int k=0; k<nMeas; k++) {
                for (int i=0; i<invNu; i++) ll[coren][i].step(nSteps);
                for (int i=0; i<invNu; i++) {
                    for (int j=0; j<invNu; j++) {
                        vector<complex<double> > temp(2);
                        temp[0]=pp[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                        overlaps[b][0](i,j)+=temp[0];
                        overlaps[b][1](i,j)+=norm(temp[0]);
                        temp[1]=ll[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                        overlaps[b][2](i,j)+=temp[1];
                        overlaps[b][3](i,j)+=norm(temp[1]);
                    }
                }
            }
            for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
            overlaps[b][0]=overlaps[b][0].array()/overlaps[b][1].array().sqrt();
            overlaps[b][2]=overlaps[b][2].array()/overlaps[b][3].array().sqrt();
            
            hermitianize(overlaps[b][2]);
        }
        //parallel programming end.

        
        vector<Eigen::MatrixXcd> berrymatrix_step(nds);
        for (int b=0; b<nds; b++) berrymatrix_step[b] = overlaps[b][2].inverse() * overlaps[b][0];
        
        Eigen::MatrixXcd berrymatrix_integral = Eigen::MatrixXcd::Identity(invNu, invNu);
        vector<double> phases(invNu, 0.);
        datas.clear();
        for (int b=0; b<nds; b++) {
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step[b]);
            data tmp;
            berrymatrix_integral *= berrymatrix_step[b];
            for (int i=0; i<invNu; i++) {
                phases[i]+=arg(es.eigenvalues()[i]);
                tmp.num = b; tmp.amp[i] = abs(es.eigenvalues()[i]); tmp.ang[i] = arg(es.eigenvalues()[i]);
            }
            
            // dfromnorm. calculates deviation from normality.
            double normeigenvalue=0., normmatrix=0.;
            for (int i=0; i<invNu; i++) {
                normeigenvalue+=sqrt(norm(es.eigenvalues()[i]));
            }
            for (int i=0; i<invNu; i++) {
                for (int j=0; j<invNu; j++) {
                    normmatrix+=sqrt(norm(berrymatrix_step[b](i,j)));
                }
            }
            tmp.dfromnorm=normmatrix-normeigenvalue;
            datas.push_back(tmp);
        }
        
        //----------
        double avephase=0.;
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_integral);
        
//        //write into outfile.
//        outfile<<"----------\nnBin="<<nbin<<" theta, alpha="<<theta/M_PI<<"pi, "<<alpha<<endl;
//        outfile<<"Ne="<<Ne<<" nMea="<<nMeas<<" length="<<length[0]<<" "<<length[1]<<" ncore="<<num_core<<endl;
//        for (int b=0; b<nds; b++) {
//            outfile<<b<<" phases = "; for (int i=0; i<invNu; i++) outfile<<datas[b].ang[i]<<" ";
//            outfile<<endl;
//            outfile<<b<<" amplitude = "; for (int i=0; i<invNu; i++) outfile<<datas[b].amp[i]<<" ";
//            outfile<<endl;
//            //                outfile<<"\n          energy = "<<energy[b]/(1.*nMeas*Ne)<<endl;
//        }
//        outfile<<"phase sum = "; for (int i=0; i<invNu; i++) {outfile<<phases[i]<<" "; avephase+=phases[i]/(1.*invNu);} outfile<<"\nphase average = "<<avephase<<endl;
//        
//        datas[0].ang_trace = arg(berrymatrix_integral.trace());
//        //    datas[0].det = arg(berrymatrix_integral.determinant());
//        outfile<<endl;
//        outfile<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
//        outfile<<"amp(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) outfile<<abs(es.eigenvalues()[i])<<" "; outfile<<endl;
//        outfile<<"arg(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) outfile<<arg(es.eigenvalues()[i])<<" "; outfile<<endl;
//        avephase=0.; for (int i=0; i<invNu; i++) avephase+=arg(es.eigenvalues()[i])/(1.*invNu); outfile<<"ave arg(berrymatrix_integral.eigenvalue) = "<<avephase<<endl;
//        outfile<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
//        outfile<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
//        outfile<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl<<endl;
//        
//        //write into outfile2. Same data, just for Mathematica convenience.
//        outfile2<<"nBin="<<nbin<<", Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<" theta, alpha="<<theta/M_PI<<"pi, "<<alpha<<endl;
//        for (int b=0; b<nds; b++) {
//            for (int i=0; i<invNu; i++) outfile2<<datas[b].ang[i]<<" ";//output phases in each step.
//            outfile2<<endl;
//            for (int i=0; i<invNu; i++) outfile2<<datas[b].amp[i]<<" ";//output amplitude in each step.
//            outfile2<<endl;
//            //                outfile2<<energy[b]/(1.*nMeas*Ne)<<endl;
//        }
//        
//        //write into outfile3. Matrix Element in each step.
//        outfile3<<"nBin="<<nbin<<", Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<" theta, alpha="<<theta/M_PI<<"pi, "<<alpha<<endl;
//        for (int b=0; b<nds; b++) {
//            outfile3<<b<<" "<<real(berrymatrix_step[b](0,0))<<" "<<imag(berrymatrix_step[b](0,0))<<" "<<real(berrymatrix_step[b](1,1))<<" "<<imag(berrymatrix_step[b](1,1))<<" "<<real(berrymatrix_step[b](1,0))<<" "<<imag(berrymatrix_step[b](1,0))<<" "<<real(berrymatrix_step[b](0,1))<<" "<<imag(berrymatrix_step[b](0,1))<<endl;
//        }

        //output (not ofstream)
//    datas[0].ang_trace = arg(berrymatrix_integral.trace());
//    datas[0].det = arg(berrymatrix_integral.determinant());
//    cout<<"\n\n Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<endl;
//    cout<<"phase sum = "; for (int i=0; i<invNu; i++) {cout<<phases[i]<<" "; avephase+=phases[i]/(1.*invNu);} cout<<"\nphase average = "<<avephase<<endl;
//    cout<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
//    cout<<"amp(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<abs(es.eigenvalues()[i])<<" "; cout<<endl;
//    cout<<"arg(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<arg(es.eigenvalues()[i])<<" ";cout<<endl;
//    avephase=0.; for (int i=0; i<invNu; i++) avephase+=arg(es.eigenvalues()[i])/(1.*invNu); cout<<"sum arg(berrymatrix_integral.eigenvalue) = "<<avephase<<endl;
//    cout<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
//    cout<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
//    cout<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl;
        
        if (arg(berrymatrix_integral.trace())<0)
            trace_phase[nbin]=arg(berrymatrix_integral.trace())+2.*M_PI;
        else
            trace_phase[nbin]=arg(berrymatrix_integral.trace());
        
    }
    
    cout<<"Ne="<<Ne<<" invNu="<<invNu<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<endl;
    
    double phase_mean=0., phase2_mean=0.;
    for (int nbin=0; nbin<nBins; nbin++) {
        phase_mean+=trace_phase[nbin]/(1.*nBins);
        phase2_mean+=trace_phase[nbin]*trace_phase[nbin]/(1.*nBins);
    }
    cout<<"phase mean:\nphase std:\n"<<phase_mean<<"\n"<<sqrt( phase2_mean-pow(phase_mean,2) )/sqrt(1.*nBins)<<endl;
    
//    outfile.close();
//    outfile2.close();
//    outfile3.close();
}
void laughlin_bp_single_state(int gs, vector<double> length, double steplength, vector<data> &datas){
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
    
    vector<vector<double> > holes; vector<int> Grid(2);
    for (int i=0; i<2; i++) {Grid[i]=(int)(length[i]/steplength);}
    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=steplength*i; a[1]=0.;                  holes.push_back(a);}
    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=steplength*i; a[0]=length[0];           holes.push_back(a);}
    for (int i=0; i<Grid[0]; i++) {vector<double> a(2); a[0]=length[0]-steplength*i; a[1]=length[1]; holes.push_back(a);}
    for (int i=0; i<Grid[1]; i++) {vector<double> a(2); a[1]=length[1]-steplength*i; a[0]=0.;        holes.push_back(a);}
    int nds=holes.size();
    vector<vector<double> > holes2(nds, vector<double>(2,0));
    int supermod(int k, int n);
    for(int i=0;i<nds;i++) holes2[supermod(i-1,nds)]=holes[i];//(holes[b],holes2[b]) = (holes[b],holes[b+1]).
    
    LATTICE ll(Ne, invNu, testing, type, seed, gs), pp(Ne, invNu, testing, type, seed, gs);
    vector<vector<complex<double> > > overlaps;
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|psi(xb)|psi(xb+1)|^2>, overlaps[b][2](i,j)=overlaps[b][0](i,j)/sqrt{overlaps[b][1](i,j)}.
    for (int b=0; b<nds; b++) {
        vector<complex<double> > aa;
        aa.push_back(0.); aa.push_back(0.); aa.push_back(0.);
        overlaps.push_back(aa);
    }
    
    for(int b=0;b<nds;b++){
        ll.set_hole(holes[b]); pp.set_hole(holes2[b]);
        ll.reset(); ll.step(nWarmup);
        for(int k=0;k<nMeas;k++){
            ll.step(nSteps);
            complex<double> temp=pp.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
            overlaps[b][0]+=temp; overlaps[b][1]+=norm(temp);
        }
        for (int l=0; l<3; l++) {
            overlaps[b][l]/=(1.*nMeas);
        }
        overlaps[b][2] = overlaps[b][0]/sqrt(overlaps[b][1]);
        // 	    cout<<holes[b][0]<<" "<<holes[b][1]<<" "<<overlaps[b][2]<<endl;
    }
    datas.clear();
    for (int b=0; b<nds; b++) {
        data tmp;
        tmp.position[0]=holes[b][0]; tmp.position[1]=holes[b][1];
        tmp.amp[gs]=abs(overlaps[b][2]); tmp.ang[gs]=arg(overlaps[b][2]);
        datas.push_back(tmp);
    }
    
    double phase=0.;
    for (int b=0; b<nds; b++) {
        phase+=datas[b].ang[gs];
    }
    cout<<"\n\nphase = "<<phase<<endl;
    cout<<endl;
    
}
void two_holes(string input_name, string str, int nmeasurement, data& test){
    int Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile(input_name);
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    vector<vector<double> > holes;
    double x=0.;
    vector<double> a(2);
    while (x<0.1) {
        a[0]=x; a[1]=0.;
        holes.push_back(a);
        x+=0.005;
    }
    int nds=holes.size();
    
    if (str=="test") {
        nMeas=nmeasurement;
        nds=1;
    }
    
    ofstream bout("twoholelaughlinnew"+to_string((long int)invNu));
    vector<LATTICE> ll(invNu),ll2(invNu);//ll is psi(x), ll2 is psi(x).
    vector<vector<Eigen::MatrixXcd> > overlaps;
    //overlaps[b][0]=<psi(0)|psi(xb)>, overlaps[b][1]=<|<psi(0)|psi(xb)>|^2>.
    for(int b=0; b<nds; b++){
        vector<Eigen::MatrixXcd> aa;
        Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(3,3);
        aa.push_back(a); aa.push_back(a);
        overlaps.push_back(aa);
    }
    
    for(int gs=0;gs<invNu;gs++){
        ll[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
        ll[gs].set_hole(holes[0]);
        ll2[gs]=LATTICE(Ne,invNu,testing,type,seed,gs);
    }
    
//    for(int gs1=0;gs1<invNu;gs1++){
//        ll[gs1].reset();
//        ll[gs1].step(nWarmup);
//        for(int i=0;i<nMeas;i++){
//            ll[gs1].step(nSteps);
//            for(int gs2=0;gs2<invNu;gs2++){
//                for(int b=0;b<nds;b++){
//                    ll2[gs2].set_hole(holes[b]);
//                    ll2[gs2].reset();
//                    complex<double> temp=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
//                    overlaps[b][0](gs1,gs2)+=temp;
//                    overlaps[b][1](gs1,gs2)+=norm(temp);
//                }
//            }
//        }
//    }
    
    for (int round=0; round<nBins; round++) {
        int gs1=0, gs2=0;
        
        for (int b=0; b<nds; b++) {
            overlaps[b][0](gs1,gs2)=0.;
            overlaps[b][1](gs1,gs2)=0.;
        }
        
        ll[gs1].reset();
        ll[gs1].step(nWarmup);
        for(int i=0;i<nMeas;i++){
            ll[gs1].step(nSteps);
            for(int b=0;b<nds;b++){
                ll2[gs2].set_hole(holes[b]);
                ll2[gs2].reset();
                complex<double> temp=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
                overlaps[b][0](gs1,gs2)+=temp;
                overlaps[b][1](gs1,gs2)+=norm(temp);
            }
        }
        for (int b=0; b<nds; b++) {
            overlaps[b][0](gs1,gs2)/=(1.*nMeas);
            overlaps[b][1](gs1,gs2)/=(1.*nMeas);
            overlaps[b][0](gs1,gs2)/=sqrt(abs(overlaps[b][1](gs1,gs2)));
            bout<<abs(overlaps[b][0](0,0))<<" "<<arg(overlaps[b][0](0,0))<<endl;
        }
//        cout<<"round="<<round<<endl;
    }
    
////    Eigen::Matrix3cd berryloop = Eigen::Matrix3cd::Identity(3,3);
//    for (int b=0; b<nds; b++) {
//        for (int l=0; l<2; l++) overlaps[b][l]/=(1.*nMeas);
//        for (int gs1=0; gs1<invNu; gs1++) {
//            for (int gs2=0; gs2<invNu; gs2++) {
//                overlaps[b][0](gs1,gs2)/=sqrt(abs(overlaps[b][1](gs1,gs2)));
//            }
//        }
////        berryloop*=overlaps[b][0];
//        
//        //output the amplitude and phase of each step.
//        bout<<abs(overlaps[b][0](0,0))<<" "<<arg(overlaps[b][0](0,0))<<endl;
//    }
//    bout.close();
    
    
//        cout<<"\nberryloop=\n"<<berryloop<<endl;
//        cout<<"determinant=\n"<<arg(berryloop.determinant())<<endl;
//    
//    for (int b=0; b<nds; b++) {
//        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b][2]);
//        bout<<holes[b][0]<<" "<<holes[b][1]<<" "<<abs(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[1])<<" "<<abs(es.eigenvalues()[2])<<" "<<arg(es.eigenvalues()[2])<<endl;
//        //        for (int i=0; i<3; i++) {
//        //            test.amp[i]=abs(es.eigenvalues()[i]);
//        //            test.ang[i]=arg(es.eigenvalues()[i]);
//        //        }
//    }
    
}
//CFL Berry Phase Calculator, Core Part.
void CFL_berry_phases_parallel(string params_name, string output_name, int num_core, string kind, double theta, double alpha){
    //input file is 'params_name'.
    //For some Ne, kind specifies certain type of berry phase loops.
    vector<data> datas;
    //output will be three files.
    //The first contains phase, amp, energy for each step each bin, and final berry matrix.
    //The second is same as the first, but in a format friend to Mathematica.
    //The third is matrix element in each step each bin.
    output_name+="_"+kind;
    ofstream outfile(output_name.c_str());
    string outfile2name=output_name+"_Mathmatica";//output data in format easy for Mathematica to deal with.
    ofstream outfile2(outfile2name.c_str());
    string outfile3name=output_name+"_Matrix";//output matrix element in each step.
    ofstream outfile3(outfile3name.c_str());
    
    int supermod(int k, int n);
    int tempNe,Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    string type;
    
    //get inputs from the params file
    ifstream infile(params_name.c_str());
    //number of electrons, and inverse of filling fraction
    infile>>tempNe>>invNu;
    //monte carlo parameters
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    //random seed
    infile>>seed;
    //if 1, will calculate the energy every step and print it.
    //This is helpful for debugging but slows down the code so should be set to 0 if the code is working
    //*IMPORTANT: every time you try a larger system size than anyone has tried before (the largest I have done is 57)
    //*you should rerun with testing=1, and make sure all the columns in the output are the same.
    //This is because increasing the size may cause floating point overflow errors
    //*
    infile>>testing;
    //  a string which chooses mode to run the code in. right now there are 4 choices:
    //  twod: put two electrons outside the circular fermi surface, move them both around
    //  oned: move one electron outside the fermi surface
    //  mtwod: removes two electrons from the circular fermi surface, i.e. adds two "holes", moves these holes around
    //  oned: removes one electron from the fermi surface
    infile>>type;
//	infile>>kind;    //......
    //tempNe is the number of electrons in the circular part of the Fermi surface
    //Ne is the total number of electrons including the extra electrons/holes
    if(type=="twod") Ne=tempNe+2;
    else if(type=="oned") Ne=tempNe+1;
    else if(type=="moned") Ne=tempNe-1;
    else if(type=="mtwod") Ne=tempNe-2;
    else if(type=="holes_and_particles") Ne=tempNe;
    else{
        cout<<"unrecognized type"<<endl;
        exit(0);
    }
    //holes==true if we are removing electrons, false otherwise
    string holes="particles";
    if(type[0]=='m') holes="holes";
    else if(type=="holes_and_particles") holes="holes_and_particles";
    
    //this instance of LATTICE is only to set up the circular fermi surface of tempNe electrons
    LATTICE templl(tempNe, invNu, testing, "CFL", seed, 0);
    //old_ds is the maximal circular fermi surface with tempNe electrons
    vector<vector<int> > old_ds=templl.get_ds(), extra_ds, remove_ds;
    //old_dbar is the center of the circular fermi surface
    vector<double> old_dbar=templl.get_dbar_parameter();
    templl.print_ds(); //exit(0);

	get_dlist(holes, tempNe, kind, extra_ds, remove_ds, old_ds);	
	
    templl.set_ds(old_ds);
//    templl.print_ds(); exit(0);
//    exit(0);
    
    int nds=extra_ds.size();
    int dsteps=nds; //nds
    
    //num_core threads to parallel the code.
    //ll is the object we will do monte carlo on, pp is the object with the electrons (or holes) shifted by one space
    omp_set_num_threads(num_core);
    
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(invNu)), pp(num_core, vector<LATTICE>(invNu));
    //do this to avoid accessing wrong memory since openmp shares memory.
    for (int k=0; k<num_core; k++)
        for (int i=0; i<invNu; i++) {
            ll[k][i]=LATTICE(Ne, invNu, testing, "CFL", seed, i, theta, alpha);
            pp[k][i]=LATTICE(Ne, invNu, testing, "CFL", seed, i, theta, alpha);
        }
    
    for (int nbin=0; nbin<nBins; nbin++) {
        //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
        vector<vector<Eigen::MatrixXcd > > overlaps(nds, vector<Eigen::MatrixXcd>(4, Eigen::MatrixXcd::Zero(invNu,invNu) ) );
        
        //parallel programming begin.
        vector<double> energy(dsteps);
#pragma omp parallel for
        for(int b=0; b<dsteps; b++) {
            int coren = omp_get_thread_num();
            //            int dKx,dKy;
            vector<vector<int> > new_ds_ll, new_ds_pp;
            new_ds_ll=old_ds;
            new_ds_pp=old_ds;
            
            //depending on the mode, this adds one or two electons just outside the Fermi surface
            if(holes=="particles"){
                new_ds_ll.push_back(extra_ds[b]);
                new_ds_pp.push_back(extra_ds[supermod(b+1,nds)]);
                //                dKx=2*(extra_ds[supermod(b+1,nds)][0]-extra_ds[b][0]);
                //                dKy=2*(extra_ds[supermod(b+1,nds)][1]-extra_ds[b][1]);
                if(type=="twod"){
                    //                    dKx=0; dKy=0;
                    new_ds_ll.push_back(extra_ds[supermod(b+nds/2,nds)]);
                    new_ds_pp.push_back(extra_ds[supermod(b+nds/2+1,nds)]);
                }
                //this removes one or two electrons from the list of ds, if we are doing holes
            }
            if(holes=="holes"){
                new_ds_ll.erase(remove(new_ds_ll.begin(),new_ds_ll.end(),extra_ds[b]),new_ds_ll.end());
                new_ds_pp.erase(remove(new_ds_pp.begin(),new_ds_pp.end(),extra_ds[supermod(b+1,nds)]),new_ds_pp.end());
                //                dKx=-2*(extra_ds[supermod(b+1,nds)][0]-extra_ds[b][0]);
                //                dKy=-2*(extra_ds[supermod(b+1,nds)][1]-extra_ds[b][1]);
                if(type=="mtwod"){
                    //                    dKx=0; dKy=0;
                    new_ds_ll.erase(remove(new_ds_ll.begin(),new_ds_ll.end(),extra_ds[supermod(b+nds/2,nds)]),new_ds_ll.end());
                    new_ds_pp.erase(remove(new_ds_pp.begin(),new_ds_pp.end(),extra_ds[supermod(b+nds/2+1,nds)]),new_ds_pp.end());
                }
            }
            if(holes=="holes_and_particles"){
                new_ds_ll.push_back(extra_ds[b]);
                new_ds_ll.erase(remove(new_ds_ll.begin(),new_ds_ll.end(),remove_ds[b]),new_ds_ll.end());
			}            	
            
            for (int i=0; i<invNu; i++) {
                pp[coren][i].set_ds(new_ds_pp);
                ll[coren][i].set_ds(new_ds_ll);
                ll[coren][i].reset(); pp[coren][i].reset();
                ll[coren][i].step(nWarmup);
            }
            //        cout<<"warmed up"<<endl;
            //        pp[0].print_ds();
            //        energy=0;
            energy[b]=0.;
            
            //many body K.
            for (int i=0; i<invNu; i++) {
                for (int j=0; j<2; j++) {
                    if (ll[coren][i].dsum[j]%invNu!=0 || pp[coren][i].dsum[j]%invNu!=0) {
                        cout<<"dsum mod invNu != 0, somewhere wrong!"<<endl;
                        exit(0);
                    }
                }
            }
            int dKx=ll[coren][0].dsum[0]/invNu-pp[coren][0].dsum[0]/invNu, dKy=ll[coren][0].dsum[1]/invNu-pp[coren][0].dsum[1]/invNu;
            //alphabar = K1 L1/Nphi + K2 L2/Nphi = d/invNu. (K1,K2) are many body momentums.
            //So for d = d1 L1/Ne + d2 L2/Ne => K1=sum d1, K2=sum d2.
            //dKx/y is divided by invNu, because in lattice.cpp dsum is defined on L/Nphi lattice.
            
            for (int k=0; k<nMeas; k++) {
                for (int i=0; i<invNu; i++)
                    ll[coren][i].step(nSteps);
                
                energy[b]+=ll[coren][0].coulomb_energy(0);
                
                for (int i=0; i<invNu; i++) {
                    for (int j=0; j<invNu; j++) {
                        vector<vector<int>> locs=ll[coren][i].get_locs();
                        int Ne=ll[coren][i].Ne, NPhi=ll[coren][i].NPhi;
                        
                        complex<double> density_matrix=0.;
//                        for (int gs=0; gs<invNu; gs++) density_matrix+=ll[coren][i].rhoq(dKx,dKy+gs*Ne,locs);//**********this is where i want to make change.
//                        density_matrix/=ll[coren][i].formfactor(dKx,dKy);//**********for the time being, I get rid of form factor.
                        
//                        if (i==j) density_matrix=ll[coren][i].rhoq(dKx,dKy,locs);
//                        else density_matrix=ll[coren][i].rhoq(dKx,dKy+Ne,locs);
                        
                        if (i==j) density_matrix=ll[coren][i].rhoq(dKx,dKy,locs);
                        else density_matrix=0.;
                        
//                        density_matrix=ll[coren][i].rhoq(dKx,dKy,locs);
                        
                        complex<double> temp=pp[coren][j].get_wf(locs)/ll[coren][i].get_wf(locs);
                        overlaps[b][0](i,j)+=temp*density_matrix;// <ll|rhoq|pp>
                        //if(i==j) cout<<norm(temp*density_matrix)<<" "<<sqrt(norm(temp))<<endl;
                        overlaps[b][1](i,j)+=norm(temp);
                        //                        temp=ll[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                        //                        overlaps[b][2](i,j)+=temp;// <ll|ll>
                        //                        overlaps[b][3](i,j)+=norm(temp);
                    }
                }
            }

            for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
            overlaps[b][0]=overlaps[b][0].array()/overlaps[b][1].array().sqrt();
            //Unlike laughlin-hole states, CFL states are orthogonal, we actually do not need overlaps[b][2].
            
            //            overlaps[b][2]=overlaps[b][2].array()/overlaps[b][3].array().sqrt();//
            //            hermitianize(overlaps[b][2]);
            //            cout<<"energy: "<<energy/(1.*nMeas*Ne)<<endl;
        }
        //parallel programming end.
        
        vector<Eigen::MatrixXcd> berrymatrix_step(dsteps);
        //        for (int b=0; b<dsteps; b++) berrymatrix_step[b] = overlaps[b][2].inverse() * overlaps[b][0];
        for (int b=0; b<dsteps; b++) berrymatrix_step[b] = overlaps[b][0];//CFL states are orthogonal.
        
        Eigen::MatrixXcd berrymatrix_integral = Eigen::MatrixXcd::Identity(invNu, invNu);
        vector<double> phases(invNu, 0.);
        datas.clear();//clear datas.
        
        stringstream filename;
        filename<<"ampout"<<Ne;
        ofstream ampout(filename.str().c_str(),ios::app);
        complex<double> complex_d1,complex_d2;
        double offset;
        
        for (int b=0; b<dsteps; b++) {
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step[b]);
            data tmp;
            berrymatrix_integral *= berrymatrix_step[b];
            for (int i=0; i<invNu; i++) {
                phases[i]+=arg(es.eigenvalues()[i]);
                tmp.num = b; tmp.amp[i] = abs(es.eigenvalues()[i]); tmp.ang[i] = arg(es.eigenvalues()[i]);
            }
            
            //write amplitude and phase of each step
            if(tempNe%2) offset=0;
            else offset=0.5;
            if(holes=="holes" or holes=="particles"){
		        complex_d1=(1.*extra_ds[b][0]-offset)/ll[0][0].getL(1)+(1.*extra_ds[b][1]-offset)/ll[0][0].getL(2);
		        complex_d2=(1.*extra_ds[(b+1)%dsteps][0]-offset)/ll[0][0].getL(1)+(1.*extra_ds[(b+1)%dsteps][1]-offset)/ll[0][0].getL(2);
			}else if(holes=="holes_and_particles"){
		        complex_d1=(1.*extra_ds[b][0]-offset)/ll[0][0].getL(1)+(1.*extra_ds[b][1]-offset)/ll[0][0].getL(2);
		        complex_d2=(1.*remove_ds[b][0]-offset)/ll[0][0].getL(1)+(1.*remove_ds[b][1]-offset)/ll[0][0].getL(2);
			}else{ cout<<"invalid holes"<<endl; }			
            ampout<<abs(complex_d1)<<" "<<abs(complex_d2)<<" "<<arg(complex_d1)<<" "<<arg(complex_d2)<<" "<<abs(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[0])<<endl;
            //            // dfromnorm. calculates deviation from normality.
            //            double normeigenvalue=0., normmatrix=0.;
            //            for (int i=0; i<invNu; i++) {normeigenvalue+=sqrt(norm(es.eigenvalues()[i]));}
            //            for (int i=0; i<invNu; i++) {
            //                for (int j=0; j<invNu; j++) {normmatrix+=sqrt(norm(berrymatrix_step[b](i,j)));}
            //            }
            //            tmp.dfromnorm=normmatrix-normeigenvalue;
            datas.push_back(tmp);
        }
        ampout.close();
        double avephase=0.;
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_integral);
        
        //write into outfile.
        outfile<<"----------\nnBin="<<nbin<<" theta, alpha="<<theta/M_PI<<"pi, "<<alpha<<endl;
        outfile<<"Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<endl;
        for (int b=0; b<dsteps; b++) {
            outfile<<"step = ("<<extra_ds[b][0]<<", "<<extra_ds[b][1]<<"),\n          phases = "; for (int i=0; i<invNu; i++) outfile<<datas[b].ang[i]<<" ";
            outfile<<"\n          amplitude = "; for (int i=0; i<invNu; i++) outfile<<datas[b].amp[i]<<" ";
            outfile<<"\n          energy = "<<energy[b]/(1.*nMeas*Ne)<<endl;
        }
        avephase=0.;
        outfile<<"phase sum = "; for (int i=0; i<invNu; i++) {outfile<<phases[i]<<" "; avephase+=phases[i]/(1.*invNu);} outfile<<"\nphase average = "<<avephase<<endl;
        
        datas[0].ang_trace = arg(berrymatrix_integral.trace());
        datas[0].det = arg(berrymatrix_integral.determinant());
        outfile<<endl;
        outfile<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
        outfile<<"amp(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) outfile<<abs(es.eigenvalues()[i])<<" "; outfile<<endl;
        outfile<<"arg(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) outfile<<arg(es.eigenvalues()[i])<<" "; outfile<<endl;
        avephase=0.; for (int i=0; i<invNu; i++) avephase+=arg(es.eigenvalues()[i])/(1.*invNu); outfile<<"ave arg(berrymatrix_integral.eigenvalue) = "<<avephase<<endl;
        outfile<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
        outfile<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
        outfile<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl<<endl;
        
        //write into outfile2. Same data, just for Mathematica convenience.
        outfile2<<"nBin="<<nbin<<", Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<" theta, alpha="<<theta/M_PI<<"pi, "<<alpha<<endl;
        for (int b=0; b<dsteps; b++) {
            for (int i=0; i<invNu; i++) outfile2<<datas[b].ang[i]<<" ";//output phases in each step.
            outfile2<<endl;
            for (int i=0; i<invNu; i++) outfile2<<datas[b].amp[i]<<" ";//output amplitude in each step.
            outfile2<<endl;
            outfile2<<energy[b]/(1.*nMeas*Ne)<<endl;
        }
        
        //write into outfile3. Matrix Element in each step.
        outfile3<<"nBin="<<nbin<<", Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<" theta, alpha="<<theta/M_PI<<"pi, "<<alpha<<endl;
        for (int b=0; b<dsteps; b++) {
            outfile3<<b<<" "<<real(berrymatrix_step[b](0,0))<<" "<<imag(berrymatrix_step[b](0,0))<<" "<<real(berrymatrix_step[b](1,1))<<" "<<imag(berrymatrix_step[b](1,1))<<" "<<real(berrymatrix_step[b](1,0))<<" "<<imag(berrymatrix_step[b](1,0))<<" "<<real(berrymatrix_step[b](0,1))<<" "<<imag(berrymatrix_step[b](0,1))<<endl;
        }
    }
    outfile.close();
    outfile2.close();
    outfile3.close();
}

//Particle Hole Symmetry (Ne9, maximal symmetric ds).
void ParticleHoleSym(){
    int Ne, invNu, seed, nMeas, nWarmup, nSteps, nBins; bool testing; string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    Ne=8; invNu=2;
    int Ne1=Ne/2, Ne2=Ne-Ne1;
    
    //cfl1 is the wavefunction that we will project into filled landau level.
    //We will see if overlap with cfl2 after projection is close to 1 or not.
    vector<LATTICE> cfl1(invNu), cfl2(invNu);
    for (int gs=0; gs<invNu; gs++) {
        cfl1[gs]=LATTICE(Ne1, invNu, testing, "CFL", seed, gs);
        cfl2[gs]=LATTICE(Ne2, invNu, testing, "CFL", seed, gs);
    }
    LATTICE doubledcfl(Ne1, invNu, testing, "doubledCFL", seed, 0);
    LATTICE FLL(Ne, 1, testing, "laughlin", seed, 0);//Filled LL Wavefunction.
//    cout<<"testing = "<<testing<<endl;
//    cfl1[1].print_ws();
    
    //monte carlo.
    complex<double> testdouble=0, testdouble_denom=0;
    for (int nbin=0; nbin<nBins; nbin++) {
        vector<Eigen::MatrixXcd> overlaps(2, Eigen::MatrixXcd::Zero(invNu, invNu));
        //overlaps[0][m][n]=<FLL|cfl1[m]*cfl2[n]>, overlaps[1][m][n]=<|FLL|cfl1[m]*cfl2[n]|^2>.
        
        FLL.reset();
        FLL.step(nWarmup);
        for (int nmea=0; nmea<nMeas; nmea++) {
            FLL.step(nSteps);
            vector<vector<int>> z=FLL.get_locs(), z1=z, z2=z;
            z1.resize(Ne1);
            z2.erase(z2.begin(), z2.begin()+Ne1);

			//inversion symmetry
            for(int j=0;j<(signed)z2.size();j++){
            	z2[j][0]=-z2[j][0];
            	z2[j][1]=-z2[j][1];
            }
            
            for (int m=0; m<invNu; m++) {
                for (int n=0; n<invNu; n++) {
                    complex<double> tmp=cfl1[m].get_wf(z1)*cfl2[n].get_wf(z2)/FLL.get_wf(z);
                    overlaps[0](m,n)+=tmp;
                    overlaps[1](m,n)+=norm(tmp);
                    if(m==0 and n==1){
		                tmp=doubledcfl.get_wf(z)/FLL.get_wf(z);
		                testdouble+=tmp;
		                testdouble_denom+=norm(tmp);
					}
                	//cout<<norm(tmp)<<endl;
                	if(m==0 and n==1) cout<<doubledcfl.get_wf(z)<<" "<<cfl1[m].get_wf(z1)*cfl2[n].get_wf(z2)<<endl;
                }
            }
            
        }
        
        for (int l=0; l<2; l++) overlaps[l]/=(1.*nMeas);
        overlaps[0]=overlaps[0].array()/overlaps[1].array().sqrt();
        testdouble/=(1.*nMeas);
        testdouble_denom/=(1.*nMeas);
        
        cout<<"nbin="<<nbin<<endl;
        for (int m=0; m<invNu; m++) {
            for (int n=0; n<invNu; n++) {
                cout<<"m="<<m<<" ,n="<<n<<" ,overlap="<<abs(overlaps[0](m,n))<<endl;
                cout<<"m="<<m<<" ,n="<<n<<" ,1-|overlap|="<<1-sqrt(norm(overlaps[0](m,n)))*sqrt(comb(Ne,Ne/2))<<endl;
                cout<<abs(testdouble/sqrt(testdouble_denom))<<endl;
                cout<<endl;
            }
        }
        cout<<endl;
    }
}
//Particle Hole Symmetry (Ne9, maximal symmetric ds).
//same as above but uses the product of 2 cfls as the weight
void ParticleHoleSymBackwards(){
    int Ne, invNu, seed, nMeas, nWarmup, nSteps, nBins; bool testing; string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    int Ne1=Ne/2;

    LATTICE cfl=LATTICE(Ne1, invNu, testing, "CFL", seed, 0);
    LATTICE FLL(Ne, 1, testing, "laughlin", seed, 0);//Filled LL Wavefunction.
    LATTICE doubledcfl(Ne1,invNu,testing,"doubledCFL",seed,0);
//    cout<<"testing = "<<testing<<endl;
//    cfl1[1].print_ws();
    
    //monte carlo.
    for (int nbin=0; nbin<nBins; nbin++) {
        vector<complex<double> > overlaps(2, 0);
        //overlaps[0][m][n]=<FLL|cfl1[m]*cfl2[n]>, overlaps[1][m][n]=<|FLL|cfl1[m]*cfl2[n]|^2>.
        
        doubledcfl.reset();
        doubledcfl.step(nWarmup);
        for (int nmea=0; nmea<nMeas; nmea++) {
            doubledcfl.step(nSteps);
            vector<vector<int>> z=doubledcfl.get_locs(), z1=z, z2=z;
//	   		for(int i=0;i<Ne;i++) cout<<"("<<z[i][0]<<","<<z[i][1]<<") ";
//			cout<<endl;
//         
            complex<double> tmp=FLL.get_wf(z)/doubledcfl.get_wf(z);
            overlaps[0]+=tmp;
            overlaps[1]+=norm(tmp);
        	//cout<<"output: "<<doubledcfl.get_wf(z)<<" "<<FLL.get_wf(z)<<endl;
            
        }
        
        for (int l=0; l<2; l++) overlaps[l]/=(1.*nMeas);
        overlaps[0]/=sqrt(overlaps[1]);
        cout<<"nbin="<<nbin<<endl;
        cout<<"overlap="<<abs(overlaps[0])<<endl;
        cout<<"1-|overlap|="<<1-abs(overlaps[0])<<endl;
        cout<<1.*doubledcfl.accepts/doubledcfl.tries<<endl;
        cout<<endl;
    }
}
void ParticleHoleSym2(){
    int Ne, invNu, seed, nMeas, nWarmup, nSteps, nBins; bool testing; string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    //this parameter object will be used to initialize LATTICE
    LATTICE_PARAMS params(Ne/invNu), paramsLL(Ne);
	paramsLL.invNu=1;
	paramsLL.testing=true;
	
    double CFLrescale,LLrescale;
    infile>>CFLrescale;
	params.rescale=CFLrescale;
	infile>>LLrescale;
	paramsLL.rescale=LLrescale;
	paramsLL.type="laughlin";

	
    params.testing=testing;
    params.seed=seed;


    //cfl1 is the wavefunction that we will project into filled landau level.
    //We will see if overlap with cfl2 after projection is close to 1 or not.
    vector<vector<wf_info>> wfs(2);
    wfs[0]=vector<wf_info>(2);
    wfs[0][0]=wf_info(false, false, 0, Ne/invNu, 1);
	wfs[0][0].wf=LATTICE(params);
	params.gs=1;
	wfs[0][1]=wf_info(false, false, Ne/invNu, Ne, -1);
	wfs[0][1].wf=LATTICE(params);
	
	wfs[1]=vector<wf_info>(1);
	wfs[1][0]=wf_info(false, false, 0, Ne, 1);
	wfs[1][0].wf=LATTICE(paramsLL);

	LATTICE_WRAPPER ll(Ne, wfs, seed, testing);

    //monte carlo.
	double denom=0, denom1, denom2;
	complex<double> tmp, num=0;
	complex<double> v1, v2, v3;

    for (int nbin=0; nbin<nBins; nbin++) {
        vector<Eigen::MatrixXcd> overlaps(2, Eigen::MatrixXcd::Zero(invNu, invNu));
        
        num=0; denom1=0; denom2=0;
        ll.reset();
        ll.step(nWarmup);
        for (int nmea=0; nmea<nMeas; nmea++) {
            ll.step(nSteps);

//			tmp=FLL.get_wf(ll.get_zs())/ll.get_wf();
//			num+=tmp;
//			denom+=norm(tmp);			
			v1=ll.get_wf(0,0);
			v2=ll.get_wf(0,1);
			v3=ll.get_wf(1,0);
			num+=v1*v2*conj(v3)/norm(ll.get_wf());
			denom1+=norm(v1*v2/ll.get_wf());
			denom2+=norm(v3/ll.get_wf());
        }

//		num/=(1.*nMeas);
//		denom/=(1.*nMeas);
//		num/=sqrt(denom);
		num/=sqrt(denom1*denom2);

        cout<<abs(num)<<" ";
		cout<<1.-abs(num)*sqrt(comb(Ne,Ne/invNu))<<endl;
    }
    ll.acceptance_rate();
}
//Particle Hole Symmetry (Ne9, maximal symmetric ds).
void Explicit(){
    int Ne, invNu, seed, nMeas, nWarmup, nSteps, nBins; bool testing; string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    //this parameter object will be used to initialize LATTICE
    LATTICE_PARAMS params(Ne/invNu);
    double tempdelta;
    complex<double> wf;
    infile>>tempdelta;
    params.w_delta=complex<double>(tempdelta,0);
    params.testing=testing;
    params.seed=seed;

   	int Ne2=Ne*Ne;

    //cfl1 is the wavefunction that we will project into filled landau level.
    //We will see if overlap with cfl2 after projection is close to 1 or not.
    vector<vector<wf_info>> wfs(2);
    wfs[0]=vector<wf_info>(2);
    wfs[0][0]=wf_info(false, false, 0, Ne/invNu, 1);
	wfs[0][0].wf=LATTICE(params);
	params.gs=1;
	wfs[0][1]=wf_info(false, false, Ne/invNu, Ne, -1);
	wfs[0][1].wf=LATTICE(params);
	
	wfs[1]=vector<wf_info>(1);
	wfs[1][0]=wf_info(false, false, 0, Ne, 1);
	wfs[1][0].wf=LATTICE(Ne, 1, testing, "laughlin", seed, 0);

	LATTICE FLL(Ne, 1, testing, "laughlin", seed, 0);
	LATTICE_WRAPPER ll(Ne, wfs, seed, testing);
	    
    vector< vector<int> > zs(Ne, vector<int>(2));
    vector<int> rawzs(Ne);
	vector<vector<int>> all_tempzs(comb(Ne2,Ne/2), vector<int>(Ne/2));

	vector<int> tempzs(Ne/2);
	for(int i=0;i<Ne/2;i++) tempzs[i]=i;
	int p;
	for(int count=0; count<(signed)all_tempzs.size(); count++){
	
		all_tempzs[count]=tempzs;
//		for(int i=0; i<Ne/2; i++) cout<<tempzs[i]<<" ";
//		cout<<endl;
		if(tempzs[0]==Ne2-Ne/2) break;	

		//find next valid conf
		for(p=0; p<Ne/2-1; p++){
			if(tempzs[p]+1<tempzs[p+1]){
				tempzs[p]++;
				break;
			}
		}
		if (p==Ne/2-1){
			tempzs[p]++;	
		}
		for(int i=0; i<p; i++) tempzs[i]=i;
	}

    int temp;
    complex<double> out=0,v1,v2,v3;
    double norm1=0,norm2=0,norm3=0, total=0;
   	vector< vector<int> >::iterator it;
   	bool duplicate, print;
	//stuff for explicit PH calculation   	
	for(int i1=0;i1<(signed)all_tempzs.size();i1++){
		cout<<i1<<endl;
		copy(all_tempzs[i1].begin(), all_tempzs[i1].end(), rawzs.begin());
		for(int i2=0; i2<(signed)all_tempzs.size();i2++){
			copy(all_tempzs[i2].begin(), all_tempzs[i2].end(), rawzs.begin()+Ne/2);
//			cout<<i1<<" "<<i2<<endl;
//			for(int p=0; p<Ne; p++) cout<<rawzs[p]<<" ";
//			cout<<endl;
		
	//		duplicate=false;
	//		for(int p=0;p<Ne;p++){
	//			temp=(i/pow(Ne2,p));
	//			if( (p>0 and p<Ne/2 and temp<=rawzs[p-1]) or (p>Ne/2 and temp<=rawzs[p-1])){
	//				duplicate=true;
	//				break;
	//			}
	//			rawzs[p]=temp%Ne2;
	//		}
	//		if(duplicate) continue;

			for(int p=0; p<Ne; p++){
				zs[p][0]=rawzs[p]/Ne;
				zs[p][1]=rawzs[p]%Ne;
			}
		
			v1=ll.get_wf(0,0,zs);
			v2=ll.get_wf(0,1,zs);
			v3=ll.get_wf(1,0,zs);
			wf=ll.get_wf(zs);

			if(abs(wf)<1e-12 and (abs(v1*v2)>1e-12 or abs(v3)>1e-12)) print=true;
			else print=false;
			//print=true;
			if(print){
				for(int p=0;p<2*Ne;p++){
					cout<<zs[p/2][p%2]<<" ";
				}
			}
		
			if(abs(wf)>1e-12){
				out+=v3*conj(v1*v2)/norm(wf)*norm(wf);
				norm2+=norm(v1*v2)/norm(wf)*norm(wf);
				norm3+=norm(v3)/norm(wf)*norm(wf);
			}

	//		out+=1./conj(v3)/v1/v2*norm(v1*v2*v3);
			if(print){
				total+=norm(v3);
				 cout<<v1<<" "<<v2<<" "<<v3<<endl;
			}
		}
	}
	norm1=1;
	cout<<"final overlap: "<<sqrt(comb(Ne,Ne/2))*abs(out/sqrt(norm1*norm2*norm3))<<" "<<total/norm3<<endl;
}
//this part of the code specifies all the grid points just outside the circle made up of tempNe electrons
//we will loop through all these positions and add electrons to them
//to go to larger sizes, it will be necessary to add more possible values of tempNe
//it may end up being more convenient to write code to automate this step
void get_dlist(string holes, int tempNe, string kind, vector< vector<int> > &extra_ds, vector< vector<int> > &remove_ds, vector< vector<int> > &old_ds){
    if(holes=="particles"){
        if(tempNe==9){
            extra_ds.push_back(vector<int>{2,0});
            extra_ds.push_back(vector<int>{2,1});
            extra_ds.push_back(vector<int>{1,2});
            extra_ds.push_back(vector<int>{0,2});
            extra_ds.push_back(vector<int>{-1,2});
            extra_ds.push_back(vector<int>{-2,1});
            extra_ds.push_back(vector<int>{-2,0});
            extra_ds.push_back(vector<int>{-2,-1});
            extra_ds.push_back(vector<int>{-1,-2});
            extra_ds.push_back(vector<int>{0,-2});
            extra_ds.push_back(vector<int>{1,-2});
            extra_ds.push_back(vector<int>{2,-1});
        }
        else if (tempNe==4){
            old_ds.clear(); extra_ds.clear();
            old_ds.push_back(vector<int>{0,0});
            old_ds.push_back(vector<int>{0,1});
            old_ds.push_back(vector<int>{1,0});
            old_ds.push_back(vector<int>{1,1});
            
            extra_ds.push_back(vector<int>{2,0});
            extra_ds.push_back(vector<int>{0,-1});
            extra_ds.push_back(vector<int>{-1,1});
            extra_ds.push_back(vector<int>{1,2});
        }
        else if(tempNe==21){
            if (kind=="fermisurface") {
                extra_ds.push_back(vector<int>{3,0});
                extra_ds.push_back(vector<int>{3,1});
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{1,3});
                extra_ds.push_back(vector<int>{0,3});
                extra_ds.push_back(vector<int>{-1,3});
                extra_ds.push_back(vector<int>{-2,2});
                extra_ds.push_back(vector<int>{-3,1});
                extra_ds.push_back(vector<int>{-3,0});
                extra_ds.push_back(vector<int>{-3,-1});
                extra_ds.push_back(vector<int>{-2,-2});
                extra_ds.push_back(vector<int>{-1,-3});
                extra_ds.push_back(vector<int>{0,-3});
                extra_ds.push_back(vector<int>{1,-3});
                extra_ds.push_back(vector<int>{2,-2});
                extra_ds.push_back(vector<int>{3,-1});
            }
            else if (kind=="lowe_21") {
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{-2,2});
                extra_ds.push_back(vector<int>{-2,-2});
                extra_ds.push_back(vector<int>{2,-2});
            }
            else if (kind=="lowe2_21") {
                extra_ds.push_back(vector<int>{3,0});
//                extra_ds.push_back(vector<int>{3,1});
                extra_ds.push_back(vector<int>{2,2});
//                extra_ds.push_back(vector<int>{1,3});
                extra_ds.push_back(vector<int>{0,3});
//                extra_ds.push_back(vector<int>{-1,3});
                extra_ds.push_back(vector<int>{-2,2});
//                extra_ds.push_back(vector<int>{-3,1});
                extra_ds.push_back(vector<int>{-3,0});
//                extra_ds.push_back(vector<int>{-3,-1});
                extra_ds.push_back(vector<int>{-2,-2});
//                extra_ds.push_back(vector<int>{-1,-3});
                extra_ds.push_back(vector<int>{0,-3});
//                extra_ds.push_back(vector<int>{1,-3});
                extra_ds.push_back(vector<int>{2,-2});
//                extra_ds.push_back(vector<int>{3,-1});
            }
            else {cout<<"unrecognized loop."<<endl; exit(0);}
        }else if(tempNe==32){
            extra_ds.push_back(vector<int>{4,0});
            extra_ds.push_back(vector<int>{4,1});
            extra_ds.push_back(vector<int>{4,2});
            extra_ds.push_back(vector<int>{3,3});
            extra_ds.push_back(vector<int>{2,4});
            extra_ds.push_back(vector<int>{1,4});
            extra_ds.push_back(vector<int>{0,4});
            extra_ds.push_back(vector<int>{-1,4});
            extra_ds.push_back(vector<int>{-2,3});
            extra_ds.push_back(vector<int>{-3,2});
            extra_ds.push_back(vector<int>{-3,1});
            extra_ds.push_back(vector<int>{-3,0});
            extra_ds.push_back(vector<int>{-3,-1});
            extra_ds.push_back(vector<int>{-2,-2});
            extra_ds.push_back(vector<int>{-1,-3});
            extra_ds.push_back(vector<int>{0,-3});
            extra_ds.push_back(vector<int>{1,-3});
            extra_ds.push_back(vector<int>{2,-3});
            extra_ds.push_back(vector<int>{3,-2});
            extra_ds.push_back(vector<int>{4,-1});
        }else if(tempNe==37){
            extra_ds.push_back(vector<int>{4,-1});
            extra_ds.push_back(vector<int>{4,0});
            extra_ds.push_back(vector<int>{4,1});
            
            extra_ds.push_back(vector<int>{3,2});
            extra_ds.push_back(vector<int>{2,3});
            
            extra_ds.push_back(vector<int>{1,4});
            extra_ds.push_back(vector<int>{0,4});
            extra_ds.push_back(vector<int>{-1,4});
            
            extra_ds.push_back(vector<int>{-2,3});
            extra_ds.push_back(vector<int>{-3,2});
            
            extra_ds.push_back(vector<int>{-4,1});
            extra_ds.push_back(vector<int>{-4,0});
            extra_ds.push_back(vector<int>{-4,-1});
            
            extra_ds.push_back(vector<int>{-3,-2});
            extra_ds.push_back(vector<int>{-2,-3});
            
            extra_ds.push_back(vector<int>{-1,-4});
            extra_ds.push_back(vector<int>{0,-4});
            extra_ds.push_back(vector<int>{1,-4});
            
            extra_ds.push_back(vector<int>{2,-3});
            extra_ds.push_back(vector<int>{3,-2});
        }else if (tempNe==57){
            extra_ds.push_back(vector<int>{5,-1});
            extra_ds.push_back(vector<int>{5,0});
            extra_ds.push_back(vector<int>{5,1});
            
            extra_ds.push_back(vector<int>{4,2});
            extra_ds.push_back(vector<int>{3,3});
            extra_ds.push_back(vector<int>{2,4});
            
            extra_ds.push_back(vector<int>{1,5});
            extra_ds.push_back(vector<int>{0,5});
            extra_ds.push_back(vector<int>{-1,5});
            
            extra_ds.push_back(vector<int>{-2,4});
            extra_ds.push_back(vector<int>{-3,3});
            extra_ds.push_back(vector<int>{-4,2});
            
            extra_ds.push_back(vector<int>{-5,1});
            extra_ds.push_back(vector<int>{-5,0});
            extra_ds.push_back(vector<int>{-5,-1});
            
            extra_ds.push_back(vector<int>{-4,-2});
            extra_ds.push_back(vector<int>{-3,-3});
            extra_ds.push_back(vector<int>{-2,-4});
            
            extra_ds.push_back(vector<int>{-1,-5});
            extra_ds.push_back(vector<int>{0,-5});
            extra_ds.push_back(vector<int>{1,-5});
            
            extra_ds.push_back(vector<int>{2,-4});
            extra_ds.push_back(vector<int>{3,-3});
            extra_ds.push_back(vector<int>{4,-2});
        }
        else if (tempNe==12){
            if (kind=="Scott1") {
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{-1,2});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{2,-1});
            }
            else if (kind=="Scott2") {
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{1,3});
                extra_ds.push_back(vector<int>{0,3});
                extra_ds.push_back(vector<int>{-1,2});
                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{-2,0});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{0,-2});
                extra_ds.push_back(vector<int>{1,-2});
                extra_ds.push_back(vector<int>{2,-1});
                extra_ds.push_back(vector<int>{3,0});
                extra_ds.push_back(vector<int>{3,1});
            }
            else if (kind=="Scott3") {
                extra_ds.push_back(vector<int>{-1,2});
                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{-2,0});
                extra_ds.push_back(vector<int>{-1,-1});
            }
            else if (kind=="Scott4") {
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{-1,2});
                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{-2,0});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{2,-1});
            }
            else if (kind=="Scott5") {
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{-1,2});
                extra_ds.push_back(vector<int>{-2,0});
                extra_ds.push_back(vector<int>{-2,1});
//                extra_ds.push_back(vector<int>{-2,0});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{2,-1});
            }
            else {cout<<"uncognized kind."<<endl; exit(0);}
        }
        else if (tempNe==69){
            if (kind=="loop1") {
                extra_ds.push_back(vector<int>{5,0});
                extra_ds.push_back(vector<int>{5,1});
                extra_ds.push_back(vector<int>{5,2});
                extra_ds.push_back(vector<int>{4,3});
                extra_ds.push_back(vector<int>{3,4});
                extra_ds.push_back(vector<int>{2,5});
                extra_ds.push_back(vector<int>{1,5});
                extra_ds.push_back(vector<int>{0,5});
                extra_ds.push_back(vector<int>{-1,5});
                extra_ds.push_back(vector<int>{-2,5});
                extra_ds.push_back(vector<int>{-3,4});
                extra_ds.push_back(vector<int>{-4,3});
                extra_ds.push_back(vector<int>{-5,2});
                extra_ds.push_back(vector<int>{-5,1});
                extra_ds.push_back(vector<int>{-5,0});
                extra_ds.push_back(vector<int>{-5,-1});
                extra_ds.push_back(vector<int>{-5,-2});
                extra_ds.push_back(vector<int>{-4,-3});
                extra_ds.push_back(vector<int>{-3,-4});
                extra_ds.push_back(vector<int>{-2,-5});
                extra_ds.push_back(vector<int>{-1,-5});
                extra_ds.push_back(vector<int>{0,-5});
                extra_ds.push_back(vector<int>{1,-5});
                extra_ds.push_back(vector<int>{2,-5});
                extra_ds.push_back(vector<int>{3,-4});
                extra_ds.push_back(vector<int>{4,-3});
                extra_ds.push_back(vector<int>{5,-2});
                extra_ds.push_back(vector<int>{5,-1});
            }
            else if (kind=="loop2") {
                extra_ds.push_back(vector<int>{4,3});
                extra_ds.push_back(vector<int>{3,4});
                extra_ds.push_back(vector<int>{-3,4});
                extra_ds.push_back(vector<int>{-4,3});
                extra_ds.push_back(vector<int>{-4,-3});
                extra_ds.push_back(vector<int>{-3,-4});
                extra_ds.push_back(vector<int>{3,-4});
                extra_ds.push_back(vector<int>{4,-3});
            }
            else if (kind=="loop3") {
                extra_ds.push_back(vector<int>{4,3});
                extra_ds.push_back(vector<int>{3,4});
                extra_ds.push_back(vector<int>{-3,4});
                extra_ds.push_back(vector<int>{-4,-3});
                extra_ds.push_back(vector<int>{-4,3});
//                extra_ds.push_back(vector<int>{-4,-3});
                extra_ds.push_back(vector<int>{-3,-4});
                extra_ds.push_back(vector<int>{3,-4});
                extra_ds.push_back(vector<int>{4,-3});
            }
            else if (kind=="loop4") {
                extra_ds.push_back(vector<int>{4,3});
                extra_ds.push_back(vector<int>{3,4});
                extra_ds.push_back(vector<int>{-3,4});
                extra_ds.push_back(vector<int>{-4,3});
                extra_ds.push_back(vector<int>{-3,-4});
                extra_ds.push_back(vector<int>{-4,-3});
                //                extra_ds.push_back(vector<int>{-3,-4});
                extra_ds.push_back(vector<int>{3,-4});
                extra_ds.push_back(vector<int>{4,-3});
            }
            else if (kind=="loop5") {
                extra_ds.push_back(vector<int>{4,3});
                extra_ds.push_back(vector<int>{3,4});
//                extra_ds.push_back(vector<int>{-3,4});
//                extra_ds.push_back(vector<int>{-4,3});
//                extra_ds.push_back(vector<int>{-4,-3});
//                extra_ds.push_back(vector<int>{-3,-4});
                extra_ds.push_back(vector<int>{3,-4});
                extra_ds.push_back(vector<int>{4,-3});
            }
            else {cout<<"unrecognized kind"<<endl; exit(0);}
            
        }
        else if (tempNe==16){
            if (kind=="triangle1") {
                old_ds.clear(); extra_ds.clear();
                for (int i=-3; i<4; i++) {
                    old_ds.push_back(vector<int>{i,0});
                }
                for (int i=-2; i<3; i++) {
                    old_ds.push_back(vector<int>{i,1});
                }
                for (int i=-1; i<2; i++) {
                    old_ds.push_back(vector<int>{i,2});
                }
                old_ds.push_back(vector<int>{0,3});
                
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{-2,2});
            }
            else if (kind=="triangle2") {
                old_ds.clear(); extra_ds.clear();
                for (int i=-3; i<4; i++) {
                    old_ds.push_back(vector<int>{i,0});
                }
                for (int i=-2; i<3; i++) {
                    old_ds.push_back(vector<int>{i,1});
                }
                for (int i=-1; i<2; i++) {
                    old_ds.push_back(vector<int>{i,2});
                }
                old_ds.push_back(vector<int>{0,3});
                
                extra_ds.push_back(vector<int>{0,4});
                extra_ds.push_back(vector<int>{-1,3});
                extra_ds.push_back(vector<int>{-2,2});
                extra_ds.push_back(vector<int>{-3,1});
                extra_ds.push_back(vector<int>{-4,0});
                extra_ds.push_back(vector<int>{-3,-1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{2,-1});
                extra_ds.push_back(vector<int>{3,-1});
                extra_ds.push_back(vector<int>{4,0});
                extra_ds.push_back(vector<int>{3,1});
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{1,3});
            }
            else if (kind=="triangle3") {
                old_ds.clear(); extra_ds.clear();
                for (int i=-3; i<4; i++) {
                    old_ds.push_back(vector<int>{i,0});
                }
                for (int i=-2; i<3; i++) {
                    old_ds.push_back(vector<int>{i,1});
                }
                for (int i=-1; i<2; i++) {
                    old_ds.push_back(vector<int>{i,2});
                }
                old_ds.push_back(vector<int>{0,3});
                
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{-2,2});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{2,-1});
            }
            else if (kind=="triangle4") {
                old_ds.clear(); extra_ds.clear();
                for (int i=-3; i<4; i++) {
                    old_ds.push_back(vector<int>{i,0});
                }
                for (int i=-2; i<3; i++) {
                    old_ds.push_back(vector<int>{i,1});
                }
                for (int i=-1; i<2; i++) {
                    old_ds.push_back(vector<int>{i,2});
                }
                old_ds.push_back(vector<int>{0,3});
                
                extra_ds.push_back(vector<int>{2,2});
                extra_ds.push_back(vector<int>{-2,2});
                extra_ds.push_back(vector<int>{-2,4});
                extra_ds.push_back(vector<int>{2,4});
            }
            else if (kind=="triangle5") {
                old_ds.clear(); extra_ds.clear();
                for (int i=-3; i<4; i++) {
                    old_ds.push_back(vector<int>{i,0});
                }
                for (int i=-2; i<3; i++) {
                    old_ds.push_back(vector<int>{i,1});
                }
                for (int i=-1; i<2; i++) {
                    old_ds.push_back(vector<int>{i,2});
                }
                old_ds.push_back(vector<int>{0,3});
                
                extra_ds.push_back(vector<int>{3,1});
                extra_ds.push_back(vector<int>{1,3});
                extra_ds.push_back(vector<int>{-1,3});
                extra_ds.push_back(vector<int>{-3,1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{2,-1});
            }
        }
        else{
            cout<<"not set up to deal with "<<tempNe<<" electrons"<<endl;
            exit(0);
        }
        //this code does the same thing as above, but it lists all the positions just inside the fermi surface, where electrons should be removed if we are doing holes
    }
    else if(holes=="holes"){
        if(tempNe==21){
            if (kind=="fermisurface") {
                extra_ds.push_back(vector<int>{2,-1});
//                extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{2,1});

//                extra_ds.push_back(vector<int>{1,2});
                extra_ds.push_back(vector<int>{0,2});
//                extra_ds.push_back(vector<int>{-1,2});

//                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{-2,0});
//                extra_ds.push_back(vector<int>{-2,-1});

                extra_ds.push_back(vector<int>{-1,-2});
//                extra_ds.push_back(vector<int>{0,-2});
                extra_ds.push_back(vector<int>{1,-2});
            }
            else if (kind=="surround0") {
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{0,-1});
            }
            else if (kind=="notsurround0") {
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{-1,-2});
                extra_ds.push_back(vector<int>{0,-2});
                extra_ds.push_back(vector<int>{1,-2});
            }
            else if (kind=="notsurround0_2") {
                extra_ds.push_back(vector<int>{1,-1});
                //extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{-1,-2});
                //extra_ds.push_back(vector<int>{0,-2});
                extra_ds.push_back(vector<int>{1,-2});
            }
            else if (kind=="loop4") {
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{0,-2});
            }
            else if (kind=="loop4test") {
//                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{0,-2});
            }
            else if (kind=="loop5") {
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{1,-2});
            }
            else if (kind=="loop6") {
                extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{1,-1});
            }
            else if (kind=="loop7") {
                extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{2,-1});
            }
            else if (kind=="loop8") {
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{0,0});
                extra_ds.push_back(vector<int>{0,-1});
            }
            else if (kind=="loop9") {
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{0,0});
                extra_ds.push_back(vector<int>{1,-1});
            }
            else if (kind=="flatten1" || kind=="flatten2" || kind=="flatten3" || kind=="flatten4"){
                old_ds.clear(); extra_ds.clear();
                for (int i=-3; i<4; i++) {
                    for (int j=-1; j<2; j++) {
                        old_ds.push_back(vector<int>{i,j});
                    }
                }
                if (kind=="flatten1") {
                    extra_ds.push_back(vector<int>{3,1});
                    extra_ds.push_back(vector<int>{-3,1});
                    extra_ds.push_back(vector<int>{-3,-1});
                    extra_ds.push_back(vector<int>{3,-1});
                }
                else if (kind=="flatten2") {
                    extra_ds.push_back(vector<int>{2,1});
                    extra_ds.push_back(vector<int>{-2,1});
                    extra_ds.push_back(vector<int>{-2,-1});
                    extra_ds.push_back(vector<int>{2,-1});
                }
                else if (kind=="flatten3") {
                    extra_ds.push_back(vector<int>{1,1});
                    extra_ds.push_back(vector<int>{-1,1});
                    extra_ds.push_back(vector<int>{-1,-1});
                    extra_ds.push_back(vector<int>{1,-1});
                }
                else if (kind=="flatten4") {
                    extra_ds.push_back(vector<int>{3,1});
                    extra_ds.push_back(vector<int>{1,1});
                    extra_ds.push_back(vector<int>{-3,1});
                    extra_ds.push_back(vector<int>{-3,-1});
                    extra_ds.push_back(vector<int>{1,-1});
                    extra_ds.push_back(vector<int>{3,-1});
                }
                else {cout<<"unrecoginzed kind."<<endl; exit(0);}
            }
            else if (kind=="cornerloop1") {
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{1,2});
                extra_ds.push_back(vector<int>{-1,2});
                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-1,-2});
                extra_ds.push_back(vector<int>{1,-2});
                extra_ds.push_back(vector<int>{2,-1});
            }
            else if (kind=="cornerloop2") {
                extra_ds.push_back(vector<int>{2,1});
                //                extra_ds.push_back(vector<int>{1,2});
                extra_ds.push_back(vector<int>{-1,2});
                //                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{-2,-1});
                //                extra_ds.push_back(vector<int>{-1,-2});
                extra_ds.push_back(vector<int>{1,-2});
                //                extra_ds.push_back(vector<int>{2,-1});
            }
            else if (kind=="cornerloop3") {
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{1,2});
                extra_ds.push_back(vector<int>{-1,2});
                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-1,-2});
                extra_ds.push_back(vector<int>{1,-2});
//                extra_ds.push_back(vector<int>{2,-1});
            }
            else{
                cout<<"unrecognized loop"<<endl;
                exit(0);
            }
        }
        else if (tempNe==3){
            old_ds.clear(); extra_ds.clear();
            old_ds.push_back(vector<int>{0,0});
            old_ds.push_back(vector<int>{0,1});
            old_ds.push_back(vector<int>{1,0});
            extra_ds.push_back(vector<int>{0,1});
            extra_ds.push_back(vector<int>{1,0});
        }
        else if(tempNe==32){
            if (kind=="fermisurface") {
                extra_ds.push_back(vector<int>{3,0});
                extra_ds.push_back(vector<int>{3,1});
                extra_ds.push_back(vector<int>{3,2});
                extra_ds.push_back(vector<int>{2,3});
                extra_ds.push_back(vector<int>{1,3});
                extra_ds.push_back(vector<int>{0,3});
                extra_ds.push_back(vector<int>{-1,3});
                extra_ds.push_back(vector<int>{-2,2});
                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{-2,0});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-1,-2});
                extra_ds.push_back(vector<int>{-0,-2});
                extra_ds.push_back(vector<int>{1,-2});
                extra_ds.push_back(vector<int>{2,-2});
                extra_ds.push_back(vector<int>{3,-1});
            }
            else if (kind=="loop1") {
                extra_ds.push_back(vector<int>{0,0});
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{1,1});
                extra_ds.push_back(vector<int>{1,0});
            }
            else if (kind=="loop2") {
                extra_ds.push_back(vector<int>{0,2});
                extra_ds.push_back(vector<int>{1,2});
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{-1,1});
            }
            else if (kind=="test_amplitude") {
                extra_ds.push_back(vector<int>{1,1});
                extra_ds.push_back(vector<int>{2,2});
            }
        }
        else if(tempNe==57){
            extra_ds.push_back(vector<int>{4,-1});
            extra_ds.push_back(vector<int>{4,0});
            extra_ds.push_back(vector<int>{4,1});
            
            extra_ds.push_back(vector<int>{3,2});
            extra_ds.push_back(vector<int>{2,3});
            
            extra_ds.push_back(vector<int>{1,4});
            extra_ds.push_back(vector<int>{0,4});
            extra_ds.push_back(vector<int>{-1,4});
            
            extra_ds.push_back(vector<int>{-2,3});
            extra_ds.push_back(vector<int>{-3,2});
            
            extra_ds.push_back(vector<int>{-4,1});
            extra_ds.push_back(vector<int>{-4,0});
            extra_ds.push_back(vector<int>{-4,-1});
            
            extra_ds.push_back(vector<int>{-3,-2});
            extra_ds.push_back(vector<int>{-2,-3});
            
            extra_ds.push_back(vector<int>{-1,-4});
            extra_ds.push_back(vector<int>{0,-4});
            extra_ds.push_back(vector<int>{1,-4});
            
            extra_ds.push_back(vector<int>{2,-3});
            extra_ds.push_back(vector<int>{3,-2});
        }
        else if(tempNe==9){
            if (kind=="loop1") {
                extra_ds.push_back(vector<int>{1,1});
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{-1,1});
            }
            else if (kind=="loop2") {
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{0,1});
            }
            else if (kind=="fullloop") {
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{1,1});
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{-1,1});
                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{1,-1});
            }
            else if (kind=="test") {
                extra_ds.push_back(vector<int>{-1,1});
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{1,1});
            }
            else if (kind=="test2") {
                extra_ds.push_back(vector<int>{-1,1});
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{1,1});
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{-1,-1});
            }
            else if (kind=="forbidline") {
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{0,0});
                extra_ds.push_back(vector<int>{0,-1});
            }
            else {cout<<"unrecognized kind."<<endl; exit(0);}
        }
        else if (tempNe==8){
            if (kind=="Scott1" || kind=="Scott2"){
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{-1,-1});
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{-1,1});
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{2,1});
                old_ds.push_back(vector<int>{2,0});
                old_ds.push_back(vector<int>{2,-1});
                old_ds.push_back(vector<int>{1,-1});
                extra_ds=old_ds;
                
                if (kind=="Scott1") {
                    old_ds.erase(remove(old_ds.begin(),old_ds.end(),vector<int>{1,-1}),old_ds.end());
                    extra_ds.erase(remove(extra_ds.begin(),extra_ds.end(),vector<int>{1,-1}),extra_ds.end());
                }
                else if (kind=="Scott2") {
                    old_ds.erase(remove(old_ds.begin(),old_ds.end(),vector<int>{2,0}),old_ds.end());
                    extra_ds.erase(remove(extra_ds.begin(),extra_ds.end(),vector<int>{2,0}),extra_ds.end());
                }
                else {cout<<"unreco kind."<<endl; exit(0);}
            }
            else if (kind=="missing0loop1") {
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{-1,1});
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{1,0});
                old_ds.push_back(vector<int>{1,-1});
                old_ds.push_back(vector<int>{0,-1});
                old_ds.push_back(vector<int>{-1,-1});
                
                extra_ds=old_ds;
            }
            else if (kind=="missing0loop2") {
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{-1,1});
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{1,0});
                old_ds.push_back(vector<int>{1,-1});
                old_ds.push_back(vector<int>{0,-1});
                old_ds.push_back(vector<int>{-1,-1});
                
                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{0,-1});
            }
            else if (kind=="missing0loop3") {
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{-1,1});
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{1,0});
                old_ds.push_back(vector<int>{1,-1});
                old_ds.push_back(vector<int>{0,-1});
                old_ds.push_back(vector<int>{-1,-1});
                
                extra_ds.push_back(vector<int>{1,1});
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{-1,1});
            }
            else {cout<<"unrecognized kind."<<endl; exit(0);}
        }
        else if(tempNe==5){
            if (kind=="loop1") {
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{0,-1});
                old_ds.push_back(vector<int>{1,0});
                extra_ds=old_ds;
                old_ds.push_back(vector<int>{0,0});
            }
            else if (kind=="loop2") {
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{1,-1});
                old_ds.push_back(vector<int>{-1,-1});
                old_ds.push_back(vector<int>{-1,1});
                extra_ds=old_ds;
                old_ds.push_back(vector<int>{0,0});
            }
        }
        else if (tempNe==12){
            if (kind=="loop1") {
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{-1,1});
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{2,1});
                old_ds.push_back(vector<int>{2,0});
                old_ds.push_back(vector<int>{2,-1});
                old_ds.push_back(vector<int>{1,-1});
                old_ds.push_back(vector<int>{0,-1});
                old_ds.push_back(vector<int>{-1,-1});
                extra_ds=old_ds;
                old_ds.push_back(vector<int>{0,0});
                old_ds.push_back(vector<int>{1,0});
            }
            else if (kind=="loop2") {
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{-1,1});
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{2,1});
                old_ds.push_back(vector<int>{2,0});
                old_ds.push_back(vector<int>{2,-1});
                old_ds.push_back(vector<int>{1,-1});
                old_ds.push_back(vector<int>{0,-1});
                old_ds.push_back(vector<int>{-1,-1});
                old_ds.push_back(vector<int>{0,0});
                old_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{-1,1});
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{2,-1});
                extra_ds.push_back(vector<int>{-1,-1});
            }
            else if(kind=="cross"){
                extra_ds.push_back(vector<int>{0,-1});
                //extra_ds.push_back(vector<int>{1,-1});

                //extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{2,1});

                //extra_ds.push_back(vector<int>{1,2});
                extra_ds.push_back(vector<int>{0,2});

                extra_ds.push_back(vector<int>{-1,1});
                //extra_ds.push_back(vector<int>{-1,0});
            
            }
            else if(kind=="loop1_new"){
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{-1,1});
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{2,1});
                old_ds.push_back(vector<int>{2,0});
                old_ds.push_back(vector<int>{0,2});
                old_ds.push_back(vector<int>{1,-1});
                old_ds.push_back(vector<int>{0,-1});
                old_ds.push_back(vector<int>{1,2});
                old_ds.push_back(vector<int>{0,0});
                old_ds.push_back(vector<int>{1,0});
                
                extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{1,2});
                extra_ds.push_back(vector<int>{0,2});
                extra_ds.push_back(vector<int>{-1,1});
                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{1,-1});
                
            }
            else if(kind=="loop2_new"){
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{-1,1});
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{2,1});
                old_ds.push_back(vector<int>{2,0});
                old_ds.push_back(vector<int>{0,2});
                old_ds.push_back(vector<int>{1,-1});
                old_ds.push_back(vector<int>{0,-1});
                old_ds.push_back(vector<int>{1,2});
                old_ds.push_back(vector<int>{0,0});
                old_ds.push_back(vector<int>{1,0});
                
                extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{2,1});
                //                extra_ds.push_back(vector<int>{1,2});
                extra_ds.push_back(vector<int>{0,2});
                extra_ds.push_back(vector<int>{-1,1});
                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{0,-1});
                //                extra_ds.push_back(vector<int>{1,-1});
                
            }
            else if(kind=="loop3_new"){
                old_ds.clear(); extra_ds.clear();
                old_ds.push_back(vector<int>{-1,0});
                old_ds.push_back(vector<int>{-1,1});
                old_ds.push_back(vector<int>{0,1});
                old_ds.push_back(vector<int>{1,1});
                old_ds.push_back(vector<int>{2,1});
                old_ds.push_back(vector<int>{2,0});
                old_ds.push_back(vector<int>{0,2});
                old_ds.push_back(vector<int>{1,-1});
                old_ds.push_back(vector<int>{0,-1});
                old_ds.push_back(vector<int>{1,2});
                old_ds.push_back(vector<int>{0,0});
                old_ds.push_back(vector<int>{1,0});
                
                extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{1,2});
                extra_ds.push_back(vector<int>{0,2});
                
                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{-1,1});
//                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{1,-1});
                
            }
            else {cout<<"unrecoginzed kind."<<endl; exit(0);}
        }
        else if (tempNe==15){
            old_ds.clear(); extra_ds.clear();
            for (int i=-1; i<2; i++) {
                for (int j=-1; j<2; j++) {
                    old_ds.push_back(vector<int>{i,j});
                }
            }
            old_ds.push_back(vector<int>{2,1});
            old_ds.push_back(vector<int>{2,0});
            old_ds.push_back(vector<int>{2,-1});
            old_ds.push_back(vector<int>{-2,1});
            old_ds.push_back(vector<int>{-2,0});
            old_ds.push_back(vector<int>{-2,-1});
            
            if (kind=="loop1") {
                extra_ds.push_back(vector<int>{1,1});
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{-1,1});
            }
            else if (kind=="loop2") {
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{2,-1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-2,1});
            }
            else if (kind=="loop2_1") {
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{2,-1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-2,1});
            }
            else if (kind=="loop3") {
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{2,-1});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{0,1});
            }
            else if (kind=="loop4") {
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{2,-1});
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-2,1});
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{1,1});
            }
            else if (kind=="forbidline") {
                extra_ds.push_back(vector<int>{2,0});
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{0,0});
                extra_ds.push_back(vector<int>{-1,0});
                extra_ds.push_back(vector<int>{-2,0});
            }
            else {cout<<"unrecognized kind."<<endl; exit(0);}
        }
        else if (tempNe==17){
            old_ds.clear(); extra_ds.clear();
            for (int i=-1; i<2; i++) {
                for (int j=-1; j<2; j++) {
                    old_ds.push_back(vector<int>{i,j});
                }
            }
            old_ds.push_back(vector<int>{2,1});
            old_ds.push_back(vector<int>{2,0});
            old_ds.push_back(vector<int>{2,-1});
            old_ds.push_back(vector<int>{-2,1});
            old_ds.push_back(vector<int>{-2,0});
            old_ds.push_back(vector<int>{-2,-1});
            old_ds.push_back(vector<int>{0,2});
            old_ds.push_back(vector<int>{0,-2});
            
            if (kind=="loop1") {
                extra_ds.push_back(vector<int>{1,1});
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{-1,1});
            }
            else if (kind=="loop2") {
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{2,-1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-2,1});
            }
            else if (kind=="loop3") {
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{-1,0});
            }
            else {cout<<"unrecognized kind."<<endl; exit(0);}
        }
        else if (tempNe==16){
            old_ds.clear(); extra_ds.clear();
            for (int i=-1; i<2; i++) {
                for (int j=-1; j<2; j++) {
                    if (i==0 && j==0) continue;
                    old_ds.push_back(vector<int>{i,j});
                }
            }
            old_ds.push_back(vector<int>{2,1});
            old_ds.push_back(vector<int>{2,0});
            old_ds.push_back(vector<int>{2,-1});
            old_ds.push_back(vector<int>{-2,1});
            old_ds.push_back(vector<int>{-2,0});
            old_ds.push_back(vector<int>{-2,-1});
            old_ds.push_back(vector<int>{0,2});
            old_ds.push_back(vector<int>{0,-2});
            
            if (kind=="loop1") {
                extra_ds.push_back(vector<int>{1,1});
                extra_ds.push_back(vector<int>{1,-1});
                extra_ds.push_back(vector<int>{-1,-1});
                extra_ds.push_back(vector<int>{-1,1});
            }
            else if (kind=="loop2") {
                extra_ds.push_back(vector<int>{2,1});
                extra_ds.push_back(vector<int>{2,-1});
                extra_ds.push_back(vector<int>{-2,-1});
                extra_ds.push_back(vector<int>{-2,1});
            }
            else if (kind=="loop3") {
                extra_ds.push_back(vector<int>{0,1});
                extra_ds.push_back(vector<int>{1,0});
                extra_ds.push_back(vector<int>{0,-1});
                extra_ds.push_back(vector<int>{-1,0});
            }
        }
        else if (tempNe==4){
            old_ds.clear(); extra_ds.clear();
            old_ds.push_back(vector<int>{1,1});
            old_ds.push_back(vector<int>{1,-1});
            old_ds.push_back(vector<int>{-1,-1});
            old_ds.push_back(vector<int>{-1,1});
            extra_ds=old_ds;
        }
        else if (tempNe==18){
            if (kind=="disconnected") {
                old_ds.clear(); extra_ds.clear();
                for (int i=0; i<3; i++) {
                    for (int j=0; j<3; j++) {
                        old_ds.push_back(vector<int>{4+i,j-1});
                        old_ds.push_back(vector<int>{-4-i,j-1});
                    }
                }
                extra_ds.push_back(vector<int>{4,1});
                extra_ds.push_back(vector<int>{4,-1});
                extra_ds.push_back(vector<int>{-4,-1});
                extra_ds.push_back(vector<int>{-4,1});
            }
            if (kind=="disconnected2") {
                old_ds.clear(); extra_ds.clear();
                for (int i=0; i<3; i++) {
                    for (int j=0; j<3; j++) {
                        old_ds.push_back(vector<int>{3+i,j-1});
                        old_ds.push_back(vector<int>{-3-i,j-1});
                    }
                }
                extra_ds.push_back(vector<int>{3,1});
                extra_ds.push_back(vector<int>{3,-1});
                extra_ds.push_back(vector<int>{-3,-1});
                extra_ds.push_back(vector<int>{-3,1});
            }
            else if (kind=="disconnected_test") {
                old_ds.clear(); extra_ds.clear();
                for (int i=0; i<3; i++) {
                    for (int j=0; j<3; j++) {
                        old_ds.push_back(vector<int>{4+i,j-1});
                        old_ds.push_back(vector<int>{-4-i,j-1});
                    }
                }
//                extra_ds.push_back(vector<int>{4,1});
                extra_ds.push_back(vector<int>{4,-1});
                extra_ds.push_back(vector<int>{-4,-1});
//                extra_ds.push_back(vector<int>{-4,1});
            }
        }
        else if (tempNe==19){
            if (kind=="disconnected") {
                old_ds.clear(); extra_ds.clear();
                for (int i=0; i<3; i++) {
                    for (int j=0; j<3; j++) {
                        old_ds.push_back(vector<int>{4+i,j-1});
                        old_ds.push_back(vector<int>{-4-i,j-1});
                    }
                }
                old_ds.push_back(vector<int>{0,0});
                
                extra_ds.push_back(vector<int>{4,1});
                extra_ds.push_back(vector<int>{4,-1});
                extra_ds.push_back(vector<int>{-4,-1});
                extra_ds.push_back(vector<int>{-4,1});
            }
        }
        else{
            cout<<"not set up to deal with "<<tempNe<<" electrons"<<endl;
            exit(0);
        }
    }else if(holes=="holes_and_particles"){
    //holes and particles is used to look at amplitudes at dTheta=0
    	if(tempNe==9){
    		extra_ds.push_back(vector<int>{0,2});
    		remove_ds.push_back(vector<int>{0,1});
 
    		extra_ds.push_back(vector<int>{0,2});
    		remove_ds.push_back(vector<int>{1,1});
 
    		extra_ds.push_back(vector<int>{1,2});
    		remove_ds.push_back(vector<int>{0,1});
 
    		extra_ds.push_back(vector<int>{2,2});
    		remove_ds.push_back(vector<int>{1,1});

    	}else if(tempNe==4){
    		extra_ds.push_back(vector<int>{2,2});
    		remove_ds.push_back(vector<int>{1,1});
    	}else if(tempNe==12){
    		extra_ds.push_back(vector<int>{2,2});
    		remove_ds.push_back(vector<int>{1,1});
    	}else if(tempNe==21){
    		extra_ds.push_back(vector<int>{2,2});
    		remove_ds.push_back(vector<int>{1,1});

    		extra_ds.push_back(vector<int>{2,0});
    		remove_ds.push_back(vector<int>{1,0});
            
            //these next 2 look weird because they are for hex
        }else if(tempNe==8){
            extra_ds.push_back(vector<int>{2,0});
            remove_ds.push_back(vector<int>{1,0});
        }else if(tempNe==13){
            extra_ds.push_back(vector<int>{2,0});
            remove_ds.push_back(vector<int>{1,0});
            
            extra_ds.push_back(vector<int>{0,1});
            remove_ds.push_back(vector<int>{0,2});
        }else if(tempNe==16){
            extra_ds.push_back(vector<int>{3,3});
            remove_ds.push_back(vector<int>{2,2});
        }else if(tempNe==25){
            extra_ds.push_back(vector<int>{3,3});
            remove_ds.push_back(vector<int>{2,2});
            
            extra_ds.push_back(vector<int>{3,0});
            remove_ds.push_back(vector<int>{2,0});
        }else if(tempNe==32){
            extra_ds.push_back(vector<int>{3,3});
        	remove_ds.push_back(vector<int>{2,2});
        }else if(tempNe==36){
        //the default here is not a square, so first 4 make it a square
            extra_ds.push_back(vector<int>{3,-2});
            extra_ds.push_back(vector<int>{-2,3});
            extra_ds.push_back(vector<int>{-2,-2});
            remove_ds.push_back(vector<int>{-3,0});
            remove_ds.push_back(vector<int>{-3,1});
            remove_ds.push_back(vector<int>{0,-3});
            remove_ds.push_back(vector<int>{0,4});
            
            extra_ds.push_back(vector<int>{4,4});
        }else{
            cout<<"not set up to deal with "<<tempNe<<" electrons"<<endl;
            exit(0);
        }
    }
}
//pairamplitude
complex<double> interaction(int m1, int m3, int m4, int No, vector<double> vpseu, string type){
    double L=sqrt(2.0*M_PI*No), gamma=(2*M_PI/L)*(2*M_PI/L);
    complex<double> k=0.;
    
    int round=10;
    for (int q1=-No*round; q1<No*round; q1++)
        for (int q2=-No*round; q2<No*round; q2++)
            if ((q2-(m1-m4))%No==0)
                for (int l=0; l<vpseu.size(); l++) {
                    
                    double temp=0.;
                    double qnorm=1.*(q1*q1+q2*q2)*gamma;
                    
                    if (vpseu[l]==0)
                        continue;
                    
                    if (type=="pa")
                        temp=Laguerrel(l,qnorm)*exp(-0.5*qnorm)/(1.*No);
                    else if (type=="ce") {
                        if (q1%No==0 && q2%No==0)
                            continue;
                        temp=1./sqrt(qnorm)*pow(Laguerrel(l, 0.5*qnorm),2)*exp(-0.5*qnorm)/(1.*No);
                    }
                    else {
                        cout<<"not set up type yet"<<endl;
                        exit(0);
                    }
       
                    k+=vpseu[l]*exp(2.0*M_PI*q1*(m1-m3)/No*complex<double>(0,1))*temp;
                }
    return k;
}
complex<double> latticepp(LATTICE ll, int m1, int m2, int m3, int m4, string type) {
    complex<double> ret=0.;
    if (ll.lat_scalex!=ll.lat_scaleq || ll.lat_scalex!=1) {
        cout<<"lat_scalex != lat_scaleq != 1"<<endl;
        exit(0);
    }
    int N=ll.NPhi*ll.lat_scalex;
    
    for (int i1=0; i1<N; i1++)
        for (int j1=0; j1<N; j1++)
            for (int i2=0; i2<N; i2++)
                for (int j2=0; j2<N; j2++) {
                    
                    if (type=="ce") {
                        ret+=
                        conj( ll.landautable[m1][i1][j1] )*conj( ll.landautable[m2][i2][j2] )
                        *
                        ll.landautable[m3][i2][j2]*ll.landautable[m4][i1][j1]
                        *
                        ll.coulomb_tableHLL[0][abs(i1-i2)][abs(j1-j2)];
                    }
                    else if (type=="pa") {
                        ret+=
                        conj( ll.landautable[m1][i1][j1] )*conj( ll.landautable[m2][i2][j2] )
                        *
                        ll.landautable[m3][i2][j2]*ll.landautable[m4][i1][j1]
                        *
                        ll.PA_table[0][abs(i1-i2)][abs(j1-j2)];
                    }
                    else {cout<<"unrecognized type."<<endl; exit(0);}
                }
    return ret;
}
void testlatticepp(double shift){
    ifstream inf("para");
    int Nphi, m1, m2, m3, m4, m5, m6, m7, m8, Kx, Ky, round, lat_scale;
    double Kq;
    string int_type;
    complex<double> tmp=0.;
    
    inf>>Nphi;
    inf>>m1>>m2>>m3>>m4;
    inf>>Kx>>Ky>>Kq;
    inf>>round;
    inf>>lat_scale;
    inf>>int_type;
    
    LATTICE ll(Nphi, 1, 0, "laughlin", 1, 0, 0.5*M_PI, 1.);
    ll.set_lat_scalex(lat_scale);
    ll.set_lat_scaleq(lat_scale);
    ll.shift_ws(shift);
    ll.setup_landautable();
    
    vector<double> exactret, latsum, latsumcom, latsumnewcom;
    //get V1234 for ED.
    for (int i=0; i<round; i++) {
        vector<double> vp=vector<double>(i,0.); vp.push_back(1.);
        if (i==0) {
            tmp=interaction(m1, m3, m4, Nphi, vp, int_type);
            exactret.push_back(1.);
        }
        else exactret.push_back(real(interaction(m1, m3, m4, Nphi, vp, int_type)/tmp));
    }
    //get V1234 from lattice sum.
    for (int i=0; i<round; i++) {
        ll.setup_tables(vector<int>{i}, int_type);
        if (i==0) {
            tmp=latticepp(ll, m1, m2, m3, m4, int_type);
            latsumnewcom.push_back(1.);
        }
        else latsumnewcom.push_back(real(latticepp(ll, m1, m2, m3, m4, int_type)/tmp));
    }
    
    //output
    cout<<"Nphi="<<Nphi<<" round="<<round<<endl;
    cout<<"m1,m2,m3,m4="<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<endl;
    cout<<"lat_scale="<<lat_scale<<endl;
    cout<<"output V_{m1m2,m3m4} with m1="<<m1<<" m2="<<m2<<" m3="<<m3<<" m4="<<m4<<endl;
    cout<<setw(30)<<" "<<"VED"<<setw(30)<<"V(lat)/VED-1"<<endl;
    for (int i=0; i<round; i++) {
        cout<<"n="<<i<<setw(30)<<setprecision(15)<<exactret[i]<<setw(30)<<latsumnewcom[i]/exactret[i]-1<<endl;
    }
}
vector<double> pairamplitude_ExplicitLatticeSum2(int invNu, double shift1, double shift2, vector<int> PP) {
    int scalex=1, scaleq=1;
    int Ne=2, NPhi=Ne*invNu;
    int Nx=scalex*NPhi, Nq=scaleq*NPhi;
    
    cout<<"Explicit Lattice Sum, Ne=2, invNu="<<invNu<<endl;
    cout<<"lat_scalex="<<scalex<<" lat_scaleq="<<scaleq<<endl;
    LATTICE ll(Ne, invNu, 0, "laughlin", 0, 0);
    
    vector<double> PA(PP.size(),0.);
    
    ll.set_lat_scalex(scalex);
    ll.set_lat_scaleq(scaleq);
    ll.shift_ws(shift1);
    ll.setup_tables(PP,"pa");
    
    vector<vector<complex<double>>> wfTable = vector<vector<complex<double>>> (Nx*Nx, vector<complex<double>>(Nx*Nx, 0.));
    for (int i=0; i<Nx*Nx; i++) {
        for (int j=0; j<Nx*Nx; j++) {
            vector<vector<int>> zs(Ne,vector<int>(2));
            zs[0][0]=i%Nx; zs[0][1]=i/Nx; zs[1][0]=j%Nx; zs[1][1]=j/Nx;
            
            vector<vector<double>> zs_double(Ne,vector<double>(2,0.));
            for (int k=0; k<Ne; k++) {
                zs_double[k][0] = zs[k][0]/(1.*Nx);
                zs_double[k][1] = zs[k][1]/(1.*Nx);
            }
            wfTable[i][j]=ll.get_laughlinwf(zs_double);
        }
    }
    
    double normalization=0.;
    for (int i=0; i<Nx*Nx; i++) {
        for (int j=0; j<Nx*Nx; j++) {
            int x=abs(i%Nx-j%Nx), y=abs(i/Nx-j/Nx);
            
            x=((i%Nx-j%Nx)%Nx+Nx)%Nx;
            y=((i/Nx-j/Nx)%Nx+Nx)%Nx;
            
            for (int k=0; k<PA.size(); k++) {
                PA[k]+=real(conj(wfTable[i][j])*wfTable[i][j]*ll.PA_table[k][x][y]);
            }
            normalization+=real(conj(wfTable[i][j])*wfTable[i][j]);
        }
    }
    
    vector<double> value(PA.size(), 0.);
    for (int i=0; i<PA.size(); i++) {
        value[i] = PA[i]/normalization;
        if (value[i]<1e-10) value[i]=0.;
    }
    return value;
}
vector<double> pairamplitude_ExplicitLatticeSum3(int invNu, double shift1, double shift2, vector<int> PP) {
    int Ne=3, NPhi=Ne*invNu;
    LATTICE ll(Ne, invNu, 0, "laughlin", 0, 0);
    cout<<"Explicit Lattice Sum, Ne=3, invNu="<<invNu<<endl;
    
    vector<double> PA(PP.size(),0.);
    ll.set_lat_scalex(1);
    ll.set_lat_scaleq(1);
    ll.shift_ws(shift1);
    ll.setup_tables(PP,"pa");
    
    vector<vector<vector<complex<double>>>> wfTable=vector<vector<vector<complex<double>>>>(NPhi*NPhi,vector<vector<complex<double>>>(NPhi*NPhi,vector<complex<double>>(NPhi*NPhi,0.)));
    
    for (int i=0; i<NPhi*NPhi; i++) {
        for (int j=0; j<NPhi*NPhi; j++) {
            for (int k=0; k<NPhi*NPhi; k++) {
                vector<vector<int>> zs(Ne,vector<int>(2));
                zs[0][0]=i%NPhi; zs[0][1]=i/NPhi; zs[1][0]=j%NPhi; zs[1][1]=j/NPhi; zs[2][0]=k%NPhi; zs[2][1]=k/NPhi;
                
                vector<vector<double>> zs_double(Ne,vector<double>(2,0.));
                for (int k=0; k<Ne; k++) {
                    zs_double[k][0] = zs[k][0]/(1.*NPhi);
                    zs_double[k][1] = zs[k][1]/(1.*NPhi);
                }
                wfTable[i][j][k]=ll.get_laughlinwf(zs_double);
            }
        }
    }
    
    double normalization=0.;
    for (int i=0; i<NPhi*NPhi; i++) {
        for (int j=0; j<NPhi*NPhi; j++) {
            for (int k=0; k<NPhi*NPhi; k++) {
                
                int x1=((i%NPhi-j%NPhi)%NPhi+NPhi)%NPhi;
                int y1=((i/NPhi-j/NPhi)%NPhi+NPhi)%NPhi;
                int x2=((j%NPhi-k%NPhi)%NPhi+NPhi)%NPhi;
                int y2=((j/NPhi-k/NPhi)%NPhi+NPhi)%NPhi;
                int x3=((k%NPhi-i%NPhi)%NPhi+NPhi)%NPhi;
                int y3=((k/NPhi-i/NPhi)%NPhi+NPhi)%NPhi;
                
                for (int m=0; m<PA.size(); m++) {
                    double val;
                    val=ll.PA_table[m][x1][y1]
                       +ll.PA_table[m][x2][y2]
                       +ll.PA_table[m][x3][y3];
                    
                    PA[m]+=real(conj(wfTable[i][j][k])*wfTable[i][j][k])*val;
                }
                normalization+=real(conj(wfTable[i][j][k])*wfTable[i][j][k]);
            }
        }
    }
    
    vector<double> value(PA.size(), 0.);
    for (int i=0; i<PA.size(); i++) {
        value[i] = PA[i]/normalization;
        if (value[i]<1e-10) value[i]=0.;
    }
    return value;
}
/*
void onebody(int m1, int m2, int NPhi, double shift, int lat_scale){
    LATTICE ll(NPhi, 1, 0, "laughlin", 1, 0, 0.5*M_PI, 1.);
    ll.set_lat_scalex(lat_scale);
    ll.set_lat_scaleq(lat_scale);
    ll.shift_ws(shift);
    ll.setup_landautable();
    
    m1=supermod(m1,NPhi);
    m2=supermod(m2,NPhi);
    
    ll.setup_sb_table(-1.);
    
    complex<double> output=0., normalization=0.;
    for (int x=0; x<lat_scale*NPhi; x++) {
        for (int y=0; y<lat_scale*NPhi; y++) {
            output+=conj(ll.landautable[m1][x][y]*ll.landautable[m2][x][y]*ll.sb_table[x][y]);
            normalization+=real( conj(ll.landautable[m1][x][y])*ll.landautable[m1][x][y] );
        }
    }
    output/=normalization;
    
    cout<<"output="<<output<<endl;
    
//    output=0.;
//    for (int x=0; x<lat_scale*NPhi; x++) {
//        for (int y=0; y<lat_scale*NPhi; y++) {
//            output+=conj(ll.landautable[m1][x][y])*ll.landautable[m2][x][y];
//        }
//    }
//    cout<<"test organ"<<output/normalization<<endl;
    
}
*/
