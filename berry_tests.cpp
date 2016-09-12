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

//test ortho laughlin states via monte carlo.
void check_orthogonality(string type){
    int Ne=8, invNu, nWarmup=5000, nMeas=5000, nSteps=20, seed=0;
    cout<<"--->type = "<<type<<", nMeas="<<nMeas<<", nSteps="<<nSteps<<", Ne="<<Ne;
    
    if (type=="CFL") cout<<". General d, sum of d neq 0."<<endl;
    
    if (type=="CFL") invNu=2;
    if (type=="laughlin") invNu=3;
    
    vector<LATTICE> cfl;
    for (int i=0; i<invNu; i++) {
        if (type=="CFL")
            cfl.push_back(LATTICE(Ne, invNu, false, "CFL", 0, i));
        
        else if (type=="laughlin")
            cfl.push_back(LATTICE(Ne, invNu, false, "laughlin", 0, i));
        
        else {
            cout<<"Type Wrong."<<endl; exit(0);
        }
    }
    
    if (type=="CFL") {
        vector<vector<int> > ds;
        for (int i=0; i<Ne; i++) {
            vector<int> tmp {i,0}; ds.push_back(tmp);
        }
        cfl[0].set_ds(ds); cfl[1].set_ds(ds);
    }
    
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
    for (int nBins=0; nBins<5; nBins++) {
        vector<Eigen::MatrixXcd > overlaps(4, Eigen::MatrixXcd::Zero(invNu, invNu));
        for (int l=0; l<2; l++) {
            cfl[l].reset();
            //            cfl[l].print_ds();
            cfl[l].step(nWarmup);
        }
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<2; i++)
                cfl[i].step(nSteps);
            for (int i=0; i<invNu; i++) {
                for (int j=0; j<invNu; j++) {
                    complex<double> temp;
                    //                temp=cfl[j].get_wf(cfl[i].get_locs())/cfl[i].get_wf(cfl[i].get_locs());
                    //                overlaps[0](i,j)+=temp;
                    //                overlaps[1](i,j)+=norm(temp);
                    temp=cfl[j].get_wf(cfl[i].get_locs())/cfl[i].get_wf(cfl[i].get_locs());
                    overlaps[2](i,j)+=temp;
                    overlaps[3](i,j)+=norm(temp);
                }
            }
        }
        
        for (int l=0; l<4; l++) overlaps[l]/=(1.*nMeas);
        //    overlaps[0]=overlaps[0].array()/overlaps[1].array().sqrt();
        overlaps[2]=overlaps[2].array()/overlaps[3].array().sqrt();
        //        hermitianize(overlaps[2]);
        cout<<"\nnBins="<<nBins<<", overlap matrix = "<<overlaps[2](0,0)<<" "<<overlaps[2](1,1)<<" "<<overlaps[2](0,1)<<" "<<overlaps[2](1,0)<<endl;
    }
}
void single_run(){
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
    
    int gs=0;
//    LATTICE ds_generator(Ne, invNu, testing, type, seed, 0);
//    vector< vector<int> > old_ds=ds_generator.get_ds();
    LATTICE ll(Ne,invNu, testing, type, seed, gs);
//    ll.set_ds(old_ds);
//    ll.print_ds();
    
    ofstream outfile("out"), eout("energy");
    for(int s=0;s<nBins;s++){
//        ll.change_dbar_parameter(s*0.1,s*0.1);
        ll.reset();
        ll.step(nWarmup);
        double E=0,E2=0;
        double P=0,P2=0,three=0;
        double e ,p;
        complex<double> berry_phase(0,0);
        deque<double> e_tracker, p_tracker;
        int Ntrack=50;
        vector<double> autocorr_e(Ntrack,0), autocorr_p(Ntrack,0);
        for(int i=0;i<nMeas;i++){
            ll.step(nSteps);
            e=ll.coulomb_energy();
            E+=e;
            E2+=e*e;
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
//            ll.update_structure_factors();
        }
        outfile<<E/(1.*nMeas*ll.Ne)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))/(1.*ll.Ne)<<" "<<real(berry_phase)/(1.*nMeas)<<" "<<imag(berry_phase)/(1.*nMeas)<<endl;
        cout<<"acceptance rate: "<<(1.*ll.accepts)/(1.*ll.tries)<<endl;
        
        ofstream auto_out("auto");
        for(int j=0;j<Ntrack;j++){
            auto_out<<j+1<<" ";
            auto_out<<autocorr_e[j]/(1.*(nMeas-Ntrack))<<" "<<pow(E/(1.*nMeas),2)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))<<" ";
//            auto_out<<autocorr_p[j]/(1.*(nMeas-Ntrack))<<" "<<pow(P/(1.*nMeas),2)<<" "<<(P2/(1.*nMeas)-pow(P/(1.*nMeas),2))<<" ";
            auto_out<<endl;
        }
    }
    outfile<<endl;
    //    ll.print_structure_factors(nMeas*nBins);
    eout.close();
    outfile.close();
}

//Energetics.
void coul_energy(LATTICE& lattice, int nWarmup, int nMeas, int nSteps, int nBins, string filename){
    ofstream outfile(filename, ios::out| ios::app);//write into the file, write from the end.
    for (int s=0; s<nBins; s++) {
        lattice.reset();
        lattice.step(nWarmup);
        double E=0, E2=0, e;
        for(int i=0;i<nMeas;i++){
            lattice.step(nSteps);
            e=lattice.coulomb_energy();
            E+=e;
            E2+=e*e;
        }
        outfile<<s<<" "<<E/(1.*nMeas)<<" "<<E2/(1.*nMeas)<<" "<<pow(E/(1.*nMeas), 2)-E2/(1.*nMeas)<<endl;
    }
    outfile<<endl;
    outfile.close();
}
void plot_CFL_coule_vsdbar(int grid, int Ne, int nMeas, int nBins){
    ofstream bout("CFLenergy"); bout.close();
    int invNu=2, Nphi=Ne*invNu, seed=0, nWarmup=5000, nSteps=20;
    LATTICE cfl(Ne, invNu, false, "CFL", seed, 0);
    
    for (int i=0; i<grid*Nphi; i++) {
        for (int j=0; j<grid*Nphi; j++) {
            double dbar_parameter[2];
            dbar_parameter[0]=i/(1.*Nphi*grid); dbar_parameter[1]=j/(1.*Nphi*grid);
            cfl.change_dbar_parameter(dbar_parameter[0],dbar_parameter[1]);
            
            ofstream outfile("CFLenergy", ios::out| ios::app);
            outfile<<i/(1.*grid)<<" "<<j/(1.*grid)<<endl;
            outfile.close();
            coul_energy(cfl, nWarmup, nMeas, nSteps, nBins, "CFLenergy");
        }
    }
}
void CFL_ne5_energy_var(int nMeas, int nBins, int num_core){
    int Ne=5, invNu=2, nSteps=20, seed=0, gs=0, nWarmup=5000;
    
    int n=6;
    vector<vector<vector<int>>> ds(n);
    if (Ne==5) {
        ds[0].push_back(vector<int>{0,0});
        ds[0].push_back(vector<int>{0,1});
        ds[0].push_back(vector<int>{1,0});
        ds[0].push_back(vector<int>{-1,0});
        ds[0].push_back(vector<int>{1,-1});
        
        ds[1].push_back(vector<int>{0,0});
        ds[1].push_back(vector<int>{0,1});
        ds[1].push_back(vector<int>{1,0});
        ds[1].push_back(vector<int>{1,1});
        ds[1].push_back(vector<int>{2,0});
        
        ds[2].push_back(vector<int>{0,0});
        ds[2].push_back(vector<int>{0,1});
        ds[2].push_back(vector<int>{1,0});
        ds[2].push_back(vector<int>{0,-1});
        ds[2].push_back(vector<int>{-1,0});
        
        ds[3].push_back(vector<int>{0,0});
        ds[3].push_back(vector<int>{1,1});
        ds[3].push_back(vector<int>{1,-1});
        ds[3].push_back(vector<int>{-1,1});
        ds[3].push_back(vector<int>{-1,-1});
        
        ds[4].push_back(vector<int>{0,0});
        ds[4].push_back(vector<int>{1,0});
        ds[4].push_back(vector<int>{-1,0});
        ds[4].push_back(vector<int>{1,1});
        ds[4].push_back(vector<int>{-1,1});
        
        ds[5].push_back(vector<int>{0,0});
        ds[5].push_back(vector<int>{0,1});
        ds[5].push_back(vector<int>{0,2});
        ds[5].push_back(vector<int>{0,3});
        ds[5].push_back(vector<int>{1,0});
    }
    
    //Initialize cfls. Write things into n different files, to avoid mass up during omp.
    vector<string> filename(n);
    vector<LATTICE> cfls(n);
    for (unsigned i=0; i<n; i++) {
        cfls[i]=LATTICE(Ne, invNu, false, "CFL", seed, gs);
    }
    
    //parallel code begin.
    omp_set_num_threads(num_core);
#pragma omp parallel for
    for (int i=0; i<n; i++) {
        int coren=omp_get_thread_num();
        cfls[i].set_ds(ds[i]);
        cfls[i].reset();
        filename[i]="cflne5evar/cflne5evar_"+to_string((long long)i);
//        cout<<i<<" "<<filename[i]<<endl;
        ofstream outfile(filename[i]);
        outfile<<i<<endl;
        coul_energy(cfls[i], nWarmup, nMeas, nSteps, nBins, filename[i]);
        outfile.close();
    }
    
    //combine outfiles.
    ofstream combine("cflne5evar/cflne5evar"); combine.close();
    for (unsigned i=0; i<n; i++) {
        ofstream combine("cflne5evar/cflne5evar", ios::out| ios::app);
        ifstream infile(filename[i]);
        string str;
        while(getline(infile, str)){
            combine<<str<<endl;
        }
        combine.close(); infile.close();
    }
    
}
void findstate(){
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
    
    ifstream dslist("dslist");
    
    vector<oneconfig> configurations;
    string tmp, ds;
    vector<string> dss;
    bool deorad=0;
    vector<int> d(2);
    int num=0;
    while (getline(dslist, tmp)) {
        oneconfig conf;
        deorad=0;
        istringstream record(tmp);
        while (record>>ds) {
            if (ds=="//") {deorad=1; continue;}
            if (!deorad) {
                d[0]=stoi(ds); record>>ds; d[1]=stoi(ds);
                conf.deletelist.push_back(d);
            }
            else {
                d[0]=stoi(ds); record>>ds; d[1]=stoi(ds);
                conf.addlist.push_back(d);
            }
        }
        conf.num=(num++);
        configurations.push_back(conf);
    }
    //    for (int i=0; i<configurations.size(); i++) {
    //        cout<<"num = "<<configurations[i].num<<endl;
    //        cout<<"deletelist = ";
    //        for (int j=0; j<configurations[i].deletelist.size(); j++) {
    //            cout<<configurations[i].deletelist[j][0]<<" "<<configurations[i].deletelist[j][1]<<", ";
    //        }
    //        cout<<endl<<"addlist = ";
    //        for (int j=0; j<configurations[i].addlist.size(); j++) {
    //            cout<<configurations[i].addlist[j][0]<<" "<<configurations[i].addlist[j][1]<<", ";
    //        }
    //        cout<<endl<<endl;
    //    }
    
    int gs=0;
    LATTICE ds_generator(9,2,testing,type,seed,0);
    vector< vector<int> > old_ds, old_ds2=ds_generator.get_ds();
    vector<int> temp_ds(2);
    
    LATTICE ll(Ne,invNu, testing, type, seed, gs);
    
    ofstream outfile("ds_energys");
    for (int k=0; k<configurations.size(); k++) {
        old_ds=old_ds2;
        vector<vector<int> > deletelist=configurations[k].deletelist;
        vector<vector<int> > addlist=configurations[k].addlist;
        for (int i=0; i<deletelist.size(); i++) {
            temp_ds[0]=deletelist[i][0]; temp_ds[1]=deletelist[i][1];
            old_ds.erase(remove(old_ds.begin(),old_ds.end(),temp_ds),old_ds.end());
            outfile<<deletelist[i][0]<<" "<<deletelist[i][1]<<" ";
        }
        outfile<<"// ";
        for (int i=0; i<addlist.size(); i++) {
            temp_ds[0]=addlist[i][0]; temp_ds[1]=addlist[i][1];
            old_ds.push_back(temp_ds);
            outfile<<addlist[i][0]<<" "<<addlist[i][1]<<" ";
        }
        outfile<<endl;
        
        ll.set_ds(old_ds);
        //    ll.print_ds();
        
        for(int s=0;s<nBins;s++){
            ll.reset();
            ll.step(nWarmup);
            double E=0,E2=0;
            //double P=0,P2=0,three=0;
            double e; //,p;
            complex<double> berry_phase(0,0);
            deque<double> e_tracker, p_tracker;
            int Ntrack=10;
            vector<double> autocorr_e(Ntrack,0), autocorr_p(Ntrack,0);
            for(int i=0;i<nMeas;i++){
                ll.step(nSteps);
                e=ll.coulomb_energy();
                E+=e;
                E2+=e*e;
            }
            outfile<<E/(1.*nMeas*ll.Ne)<<" "<<(E2/(1.*nMeas)-pow(E/(1.*nMeas),2))/(1.*ll.Ne)<<" ";
            //            cout<<"acceptance rate: "<<(1.*ll.accepts)/(1.*ll.tries)<<endl;
        }
        outfile<<endl<<endl;
    }
    outfile.close();
}

//laughlin-hole berry phase.
void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne, int num_core){
    int Ne,Ne_t,invNu,nWarmup,nMeas,nMeas_t,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile("params");
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
    vector<vector<Eigen::MatrixXcd > > overlaps(nds, vector<Eigen::MatrixXcd>(4, Eigen::MatrixXcd::Zero(3, 3)));
    
    omp_set_num_threads(num_core);
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(3)), pp(num_core, vector<LATTICE>(3));//do this to avoid wrong memory access since openmp share memory.
    for (int k=0; k<num_core; k++) for (int i=0; i<3; i++) {ll[k][i]=LATTICE(Ne, invNu, testing, type, seed, i); pp[k][i]=LATTICE(Ne, invNu, testing, type, seed, i);}
    
    //parallel programming begin.
#pragma omp parallel for
    for(int b=0; b<nds; b++) {
        int coren = omp_get_thread_num();
        for (int i=0; i<3; i++) {
            ll[coren][i].set_hole(holes[b]);
            pp[coren][i].set_hole(holes2[b]);
            ll[coren][i].reset(); ll[coren][i].step(nWarmup);
            pp[coren][i].reset();
        }
        
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<3; i++) ll[coren][i].step(nSteps);
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
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
    
    double avephase;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_integral);
    datas[0].ang_trace = arg(berrymatrix_integral.trace());
    datas[0].det = arg(berrymatrix_integral.determinant());
    cout<<"\n\n Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<endl;
    cout<<"phase sum = "; for (int i=0; i<invNu; i++) {cout<<phases[i]<<" "; avephase+=phases[i]/(1.*invNu);} cout<<"\nphase average = "<<avephase<<endl;
    cout<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
    cout<<"amp(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<abs(es.eigenvalues()[i])<<" "; cout<<endl;
    cout<<"arg(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<arg(es.eigenvalues()[i])<<" ";cout<<endl;
    avephase=0.; for (int i=0; i<invNu; i++) avephase+=arg(es.eigenvalues()[i])/(1.*invNu); cout<<"sum arg(berrymatrix_integral.eigenvalue) = "<<avephase<<endl;
    cout<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
    cout<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
    cout<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl;
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
void two_holes(string str, int nmeasurement, data& test){
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
    
    vector<vector<double> > holes;
    double x=0.;
    vector<double> a(2);
    while (x<0.05) {
        a[0]=x; a[1]=0.;
        holes.push_back(a);
        x+=0.005;
    }
    int nds=holes.size();
    
    if (str=="test") {
        nMeas=nmeasurement;
        nds=1;
    }
    
    ofstream bout("twoholelaughlinnew");
    vector<LATTICE> ll(invNu),ll2(invNu);//ll is psi(x), ll2 is psi(x).
    vector<vector<Eigen::MatrixXcd> > overlaps;
    //overlaps[b][0]=<psi(0)|psi(xb)>, overlaps[b][1]=<|<psi(0)|psi(xb)>|^2>.
    for (int b=0; b<nds; b++) {
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
    
    for(int gs1=0;gs1<invNu;gs1++){
        ll[gs1].reset();
        ll[gs1].step(nWarmup);
        for(int i=0;i<nMeas;i++){
            ll[gs1].step(nSteps);
            for(int gs2=0;gs2<invNu;gs2++){
                for(int b=0;b<nds;b++){
                    ll2[gs2].set_hole(holes[b]);
                    ll2[gs2].reset();
                    complex<double> temp=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
                    overlaps[b][0](gs1,gs2)+=temp;
                    overlaps[b][1](gs1,gs2)+=norm(temp);
                }
            }
        }
    }
    
    /*
     for (int b=0; b<nds; b++) {
     for (int gs1=0;gs1<invNu;gs1++) {
     ll[gs1].reset();
     ll[gs1].step(nWarmup);
     for (int gs2=0; gs2<invNu; gs2++) {
     ll2[gs2].set_hole(holes[b]);
     ll2[gs2].reset();
     for (int i=0;i<nMeas;i++) {
     ll[gs1].step(nSteps);
     complex<double> temp=ll2[gs2].get_wf(ll[gs1].get_locs())/ll[gs1].get_wf(ll[gs1].get_locs());
     overlaps[b][0](gs1,gs2)+=temp;
     overlaps[b][1](gs1,gs2)+=norm(temp);
     }
     }
     }
     }
     cout<<"\noverlaps[10][0]=\n"<<overlaps[10][0]<<endl;
     cout<<"\noverlaps[10][1]=\n"<<overlaps[10][1]<<endl;
     */
    
    Eigen::Matrix3cd berryloop = Eigen::Matrix3cd::Identity(3,3);
    for (int b=0; b<nds; b++) {
        for (int l=0; l<2; l++) overlaps[b][l]/=(1.*nMeas);
        for (int gs1=0; gs1<invNu; gs1++) {
            for (int gs2=0; gs2<invNu; gs2++) {
                overlaps[b][0](gs1,gs2)/=sqrt(abs(overlaps[b][1](gs1,gs2)));
            }
        }
        berryloop*=overlaps[b][0];
    }
    
    //    cout<<"\nberryloop=\n"<<berryloop<<endl;
    //    cout<<"determinant=\n"<<arg(berryloop.determinant())<<endl;
    
    for (int b=0; b<nds; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b][2]);
        bout<<holes[b][0]<<" "<<holes[b][1]<<" "<<abs(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[1])<<" "<<abs(es.eigenvalues()[2])<<" "<<arg(es.eigenvalues()[2])<<endl;
        //        for (int i=0; i<3; i++) {
        //            test.amp[i]=abs(es.eigenvalues()[i]);
        //            test.ang[i]=arg(es.eigenvalues()[i]);
        //        }
    }
    
}
void test_error(int ne, double loop, double steplength, int nMea, int ncore, string test, int num_core){
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne, int num_core);
    vector<double> length(2); vector<data> datas;
    length[0]=loop; length[1]=loop;
    
    ofstream bout("test");
    if (test=="steplength") {
        for (double steplength=0.001; steplength<loop; steplength+=0.002) {
            laughlinberryphase(length, steplength, datas, nMea, ne, num_core);
            vector<double> phase(3);
            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
            bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
        }
    }
    if (test=="ne") {
        for (int ne=2; ne<22; ne=ne+2) {
            laughlinberryphase(length, steplength, datas, nMea, ne, num_core);
            vector<double> phase(3);
            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
            bout<<ne<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
        }
    }
    if (test=="nMea") {
        for (int nMeas=10; nMeas<20; nMeas+=2){
            laughlinberryphase(length, steplength, datas, nMea, ne, num_core);
            vector<double> phase(3);
            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
            bout<<nMeas<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
        }
    }
    if (test=="loop") {
        for (double x=0.05; x<0.8; x+=0.05) {
            length[0]=x; length[1]=x;
            laughlinberryphase(length, steplength, datas, nMea, ne, num_core);
            vector<double> phase(3);
            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
            bout<<x<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
        }
    }
    if (test=="normality") {
        for (double steplength=0.01; steplength<0.25; steplength+=0.01) {
            laughlinberryphase(length, steplength, datas, 0, 0, num_core);
            
            bout<<steplength<<endl;
            for (int i=0; i<datas.size(); i++) {
                bout<<datas[i].num<<" "<<datas[i].dfromnorm<<endl;
            }
            bout<<endl;
        }
    }
}

//CFL Berry Phase Calculator, Core Part.
void CFL_berry_phases_parallel(string params_name, string output_name, int num_core, string kind){
    //input file is 'params_name'.
    //For some Ne, kind specifies certain type of berry phase loops.
    vector<data> datas;
    //output will be three files.
    //The first contains phase, amp, energy for each step each bin, and final berry matrix.
    //The second is same as the first, but in a format friend to Mathematica.
    //The third is matrix element in each step each bin.
    output_name+=output_name+"_"+kind;
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
    /*IMPORTANT: every time you try a larger system size than anyone has tried before (the largest I have done is 57)
     *you should rerun with testing=1, and make sure all the columns in the output are the same.
     This is because increasing the size may cause floating point overflow errors
     */
    infile>>testing;
    //  a string which chooses mode to run the code in. right now there are 4 choices:
    //  twod: put two electrons outside the circular fermi surface, move them both around
    //  oned: move one electron outside the fermi surface
    //  mtwod: removes two electrons from the circular fermi surface, i.e. adds two "holes", moves these holes around
    //  oned: removes one electron from the fermi surface
    infile>>type;
    
    //tempNe is the number of electrons in the circular part of the Fermi surface
    //Ne is the total number of electrons including the extra electrons/holes
    if(type=="twod") Ne=tempNe+2;
    else if(type=="oned") Ne=tempNe+1;
    else if(type=="moned") Ne=tempNe-1;
    else if(type=="mtwod") Ne=tempNe-2;
    else{
        cout<<"unrecognized type"<<endl;
        exit(0);
    }
    //holes==true if we are removing electrons, false otherwise
    bool holes=false;
    if(type[0]=='m') holes=true;
    
    //this instance of LATTICE is only to set up the circular fermi surface of tempNe electrons
    LATTICE templl(tempNe, invNu, testing, "CFL", seed, 0);
    //old_ds is the maximal circular fermi surface with tempNe electrons
    vector<vector<int> > old_ds=templl.get_ds(), extra_ds;
    //old_dbar is the center of the circular fermi surface
    vector<double> old_dbar=templl.get_dbar_parameter();
    templl.print_ds();
    
    //this part of the code specifies all the grid points just outside the circle made up of tempNe electrons
    //we will loop through all these positions and add electrons to them
    //to go to larger sizes, it will be necessary to add more possible values of tempNe
    //it may end up being more convenient to write code to automate this step
    if(!holes){
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
        }else if(tempNe==21){
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
        }else{
            cout<<"not set up to deal with "<<tempNe<<" electrons"<<endl;
            exit(0);
        }
        //this code does the same thing as above, but it lists all the positions just inside the fermi surface, where electrons should be removed if we are doing holes
    }
    else if(holes){
        if(tempNe==21){
            if (kind=="fermisurface") {
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
            else{
                cout<<"unrecognized loop"<<endl;
                exit(0);
            }
            
        }
        else if(tempNe==32){
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
            else {cout<<"unrecognized kind."<<endl; exit(0);}
        }
        else if(tempNe==5){
            extra_ds.push_back(vector<int>{1,0});
            extra_ds.push_back(vector<int>{0,1});
            extra_ds.push_back(vector<int>{-1,0});
            extra_ds.push_back(vector<int>{0,-1});
        }else{
            cout<<"not set up to deal with "<<tempNe<<" electrons"<<endl;
            exit(0);
        }
    }
    
    int nds=extra_ds.size();
    int dsteps=nds; //nds
    
    //num_core threads to parallel the code.
    //ll is the object we will do monte carlo on, pp is the object with the electrons (or holes) shifted by one space
    omp_set_num_threads(num_core);
    
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(invNu)), pp(num_core, vector<LATTICE>(invNu));
    //do this to avoid accessing wrong memory since openmp shares memory.
    for (int k=0; k<num_core; k++) for (int i=0; i<invNu; i++) {ll[k][i]=LATTICE(Ne, invNu, testing, "CFL", seed, i); pp[k][i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);}
    
    for (unsigned nbin=0; nbin<nBins; nbin++) {
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
            if(!holes){
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
            if(holes){
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
            
            for (int i=0; i<invNu; i++) {
                ll[coren][i].set_ds(new_ds_ll);
                pp[coren][i].set_ds(new_ds_pp);
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
                
                energy[b]+=ll[coren][0].coulomb_energy();
                
                for (int i=0; i<invNu; i++) {
                    for (int j=0; j<invNu; j++) {
                        complex<double> temp;
                        temp=pp[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                        overlaps[b][0](i,j)+=temp*ll[coren][i].rhoq(dKx,dKy,ll[coren][i].get_locs());// <ll|rhoq|pp>
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
        for (int b=0; b<dsteps; b++) {
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step[b]);
            data tmp;
            berrymatrix_integral *= berrymatrix_step[b];
            for (int i=0; i<invNu; i++) {
                phases[i]+=arg(es.eigenvalues()[i]);
                tmp.num = b; tmp.amp[i] = abs(es.eigenvalues()[i]); tmp.ang[i] = arg(es.eigenvalues()[i]);
            }
            //            // dfromnorm. calculates deviation from normality.
            //            double normeigenvalue=0., normmatrix=0.;
            //            for (int i=0; i<invNu; i++) {normeigenvalue+=sqrt(norm(es.eigenvalues()[i]));}
            //            for (int i=0; i<invNu; i++) {
            //                for (int j=0; j<invNu; j++) {normmatrix+=sqrt(norm(berrymatrix_step[b](i,j)));}
            //            }
            //            tmp.dfromnorm=normmatrix-normeigenvalue;
            datas.push_back(tmp);
        }
        
        double avephase=0.;
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_integral);
        
        //write into outfile.
        outfile<<"----------\nnBin="<<nbin<<endl;
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
        outfile2<<"nBin="<<nbin<<", Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<endl;
        for (int b=0; b<dsteps; b++) {
            for (int i=0; i<invNu; i++) outfile2<<datas[b].ang[i]<<" ";//output phases in each step.
            outfile2<<endl;
            for (int i=0; i<invNu; i++) outfile2<<datas[b].amp[i]<<" ";//output amplitude in each step.
            outfile2<<endl;
            outfile2<<energy[b]/(1.*nMeas*Ne)<<endl;
        }
        
        //write into outfile3. Matrix Element in each step.
        outfile3<<"nBin="<<nbin<<", Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<num_core<<endl;
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
    
    Ne=18; invNu=2;
    int Ne1=9, Ne2=Ne-Ne1;
    
    //cfl1 is the wavefunction that we will project into filled landau level.
    //We will see if overlap with cfl2 after projection is close to 1 or not.
    vector<LATTICE> cfl1(invNu), cfl2(invNu);
    for (unsigned gs=0; gs<invNu; gs++) {
        cfl1[gs]=LATTICE(Ne1, invNu, testing, "CFL", seed, gs);
        cfl2[gs]=LATTICE(Ne2, invNu, testing, "CFL", seed, gs);
    }
    LATTICE FLL(Ne, 1, testing, "FilledLL", seed, 0);//Filled LL Wavefunction.
//    cout<<"testing = "<<testing<<endl;
//    cfl1[1].print_ws();
    
    //monte carlo.
    for (unsigned nbin=0; nbin<nBins; nbin++) {
        vector<Eigen::MatrixXcd> overlaps(2, Eigen::MatrixXcd::Zero(invNu, invNu));
        //overlaps[0][m][n]=<FLL|cfl1[m]*cfl2[n]>, overlaps[1][m][n]=<|FLL|cfl1[m]*cfl2[n]|^2>.
        
        FLL.reset();
        FLL.step(nWarmup);
        for (int nmea=0; nmea<nMeas; nmea++) {
            FLL.step(nSteps);
            vector<vector<int>> z=FLL.get_locs(), z1=z, z2=z;
            z1.resize(Ne1);
            z2.erase(z2.begin(), z2.begin()+Ne1);
            
            for (int m=0; m<invNu; m++) {
                for (int n=0; n<invNu; n++) {
                    complex<double> tmp=cfl1[m].get_wf(z1)*cfl2[n].get_wf(z2)/FLL.get_wf(z);
                    overlaps[0](m,n)+=tmp;
                    overlaps[1](m,n)+=norm(tmp);
                }
            }
            
        }
        
        for (int l=0; l<2; l++) overlaps[l]/=(1.*nMeas);
        overlaps[0]=overlaps[0].array()/overlaps[1].array().sqrt();
        cout<<"nbin="<<nbin<<endl;
        for (int m=0; m<invNu; m++) {
            for (int n=0; n<invNu; n++) {
                cout<<"m="<<m<<" ,n="<<n<<" ,overlap="<<overlaps[0](m,n)<<endl;
                cout<<"m="<<m<<" ,n="<<n<<" ,1-|overlap|="<<1-sqrt(norm(overlaps[0](m,n)))<<endl;
                cout<<endl;
            }
        }
        cout<<endl;
    }
}

void testIQHwf(){
    int Ne, invNu, seed, nMeas, nWarmup, nSteps, nBins; bool testing; string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
    Ne=4, invNu=1;
    LATTICE lat1(Ne, invNu, testing, "laughlin", seed, 0);
    LATTICE lat2(Ne, invNu, testing, "FilledLL", seed, 0);
    lat1.print_ws();
    
    lat1.reset(); lat1.step(nWarmup); complex<double> value=0., value2=0.;;
    for (int nmeas=0; nmeas<nMeas; nmeas++) {
        lat1.step(nSteps);
        complex<double> tmp=lat2.get_wf(lat1.get_locs())/lat1.get_wf(lat1.get_locs());
        value+=tmp;
        value2+=norm(tmp);
    }
    value/=(1.*nMeas); value2/=(1.*nMeas);
    value/=sqrt(value2);
    cout<<"overlap="<<value<<endl;
    
}
