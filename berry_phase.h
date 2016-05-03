#ifndef BERRY_PHASE_H
#define BERRY_PHASE_H


class berry_phase{
    friend class LATTICE;
private:
    int NPhi, invNu, nWarmup, nMeas, nSteps, nBins, seed;
    bool testing;
    string type;
    
    vector <vector <int> > tempds;
    
public:
    berry_phase(int);
    int nds;
    vector<int> dpoint;
    vector<vector<int> > dcenter, d2, d3;
    double energy;
    complex<double> berry2,berry3;
    
    double phasemod(complex<double>);
    
    void two_full_braiding();
    void two_half_braiding();
    void three_full_braiding();
};
berry_phase::berry_phase(int nds_t){
    ifstream infile("params");
    infile>>NPhi>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    
    if(NPhi!=46){
        cout<<"the berry phase calculator only works for 44 flux quanta right now!"<<endl;
        exit(0);
    }
    
    dpoint=vector<int>(2);
    nds=nds_t;
    
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
    int supermod(int k, int n);
    for(int i=0;i<nds;i++){
        d2[supermod(i-1,nds)]=dcenter[i];
        d3[supermod(i+1,nds)]=dcenter[i];
    }
}
double berry_phase::phasemod(complex<double> in){
    double out=arg(in);
//    if(out<0) return out+2*M_PI;
//    else if (out>2*M_PI) return out-2*M_PI;
//    else return out;
    return out;
}
void berry_phase::two_full_braiding(){
    double tot_berry_phase=0.;
    
    ofstream bout("berry");
    LATTICE ll(NPhi,invNu, testing, type, seed);
    LATTICE ll2(NPhi,invNu, testing, type, seed);
//    LATTICE ll3(NPhi,invNu, testing, type, seed);
    
    for(int b=0;b<nds;b++){
        //stuff for initializing with a custom set of ds (for the Berry phase calculation)
        double center_frac[2]={0.,0.}; //might need to change this for different N, but it works for N=22
        ll.make_fermi_surface(center_frac, ll.Ne-2);
        tempds=ll.get_ds();//ds are defined on Ne lattice.
        ll.print_ds();
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
        
//        ll3.make_fermi_surface(center_frac,ll.Ne-2);
//        tempds=ll3.get_ds();
//        dpoint[0]=d3[b][0]; dpoint[1]=d3[b][1];
//        tempds.push_back(dpoint);
//        dpoint[0]=-d3[b][0]; dpoint[1]=-d3[b][1];
//        tempds.push_back(dpoint);
//        ll3.set_ds(tempds);
        
        energy=0;
        berry2=0; berry3=0;
        ll.step(nWarmup);
        for(int i=0;i<nMeas;i++){
            ll.step(nSteps);
            energy+=ll.coulomb_energy();
            berry2+=ll2.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
//            berry3+=ll3.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
        }
        bout<<dcenter[b][0]<<" "<<dcenter[b][1]<<" "<<d2[b][0]<<" "<<d2[b][1]<<"   "<<abs(berry2)/(1.*nMeas)<<"   "<<phasemod(berry2)<<"   "<<energy/(1.*nMeas*ll.Ne)<<endl;
        bout<<dcenter[b][0]<<" "<<dcenter[b][1]<<" "<<d3[b][0]<<" "<<d3[b][1]<<"   "<<abs(berry3)/(1.*nMeas)<<"   "<<phasemod(berry3)<<"   "<<energy/(1.*nMeas*ll.Ne)<<endl;
        tot_berry_phase+=phasemod(berry2);
        cout<<"b="<<b<<endl;
    }
    
    cout<<"total berry phase = "<<tot_berry_phase<<endl;
}



#endif
