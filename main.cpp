using namespace std;
#include "lattice.h"
#include "berry_tests.h"
//#include <iomanip>
bool IsOdd (int i) {
    return ((i%2)==1);
}
int main(){
    //CFL berry phase.
    void CFL_berry_phases_parallel(vector<data> &datas, string params_name, string output_name, int num_core);
    vector<data> datas; int num_core;
    num_core=4;
//    CFL_berry_phases_parallel(datas, "params_ne21", "CFL_berryphase_ne21", num_core);
    void CFL_ne20hole_berryphase(vector<data> &datas, string params_name, string out_name, int num_core, string loop);
//    CFL_ne20hole_berryphase(datas, "params_ne21", "CFL_berryphase_ne21", num_core, "fermisurface");
//    CFL_ne20hole_berryphase(datas, "params_ne21", "CFL_berryphase_ne21", num_core, "surround0");
//    CFL_ne20hole_berryphase(datas, "params_ne21", "CFL_berryphase_ne21", num_core, "notsurround0");
//    CFL_ne20hole_berryphase(datas, "params_ne21", "CFL_berryphase_ne21", num_core, "loop4");
//    CFL_ne20hole_berryphase(datas, "params_ne21_m=4", "CFL_berryphase_ne21_m=4", num_core, "fermisurface");
    
    void onestep(int ne, string output_name);
    onestep(20, "M_20");
    
}
void test_laughlinwf(){
    bool testing=false; string type="laughlin"; int seed=0;
    int Ne=2, invNu=3, NPhi=6;
    
    LATTICE lau(Ne, invNu, testing, type, seed, 0);
    vector<vector<int> > zs(2, vector<int>(2));

    zs[0]=vector<int>{1, 1}; zs[1]=vector<int>{2, 2};
    cout<<"w.f. = "<<lau.get_wf(zs)<<endl;
    
    int nWarmup=5000, nMeas=10000, nSteps=10;
    for (int m=1; m<2; m++) {
        for (int n=m; n<2; n++) {
            LATTICE ll, pp;
            ll=LATTICE(Ne, invNu, testing, type, seed, m);
            pp=LATTICE(Ne, invNu, testing, type, seed, n);
            ll.reset(); ll.step(nWarmup);
            
            vector<complex<double> > overlaps(2);
            for (int k=0; k<nMeas; k++) {
                ll.step(nSteps);
                complex<double> temp=pp.get_wf(ll.get_locs())/ll.get_wf(ll.get_locs());
                overlaps[0]+=temp; overlaps[1]+=norm(temp);
            }
            for (int l=0; l<2; l++) overlaps[l]/=(1.*nMeas);
            overlaps[0]/=sqrt(abs(overlaps[1]));
            cout<<"m="<<m<<" n="<<n<<" overlap="<<overlaps[0]<<endl;
        }
    }
   
}
void CFL_berry_phases(vector<data> &datas){
    int supermod(int k, int n);
	int tempNe,Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
	bool testing;
	string type;

	//get  inputs from the params file
	ifstream infile("params");
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
	//a string which chooses mode to run the code in. right now there are 4 choices:  
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
	vector<vector<int> > old_ds=templl.get_ds(), new_ds_ll, new_ds_pp, extra_ds;
	//old_dbar is the center of the circular fermi surface
	vector<double> old_dbar=templl.get_dbar_parameter();

	//this part of the code specifies all the grid points just outside the circle made up of tempNe electrons
	//we will loop through all these positions and add electrons to them
	//to go to larger sizes, it will be necessary to add more possible values of tempNe
	//it may end up being more convenient to write code to automate this step
	if(!holes){
		if(tempNe==4){
			extra_ds.push_back(vector<int>{2,0});	
			extra_ds.push_back(vector<int>{2,1});	
			extra_ds.push_back(vector<int>{1,2});	
			extra_ds.push_back(vector<int>{0,2});	
			extra_ds.push_back(vector<int>{-1,1});	
			extra_ds.push_back(vector<int>{-1,0});	
			extra_ds.push_back(vector<int>{0,-1});	
			extra_ds.push_back(vector<int>{1,-1});	
		}else if(tempNe==21){
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
	}else{
		if(tempNe==9){
			extra_ds=vector< vector<int> > (4, vector<int>(2));
			extra_ds[0][0]=1; extra_ds[0][1]=1;
			extra_ds[1][0]=1; extra_ds[1][1]=-1;
			extra_ds[2][0]=-1; extra_ds[2][1]=-1;
			extra_ds[3][0]=-1; extra_ds[3][1]=1;
		}else if(tempNe==21){
			extra_ds.push_back(vector<int>{1,0});	
			extra_ds.push_back(vector<int>{0,1});	
			extra_ds.push_back(vector<int>{-1,0});	
			extra_ds.push_back(vector<int>{0,-1});	

//			extra_ds.push_back(vector<int>{2,0});	
//			extra_ds.push_back(vector<int>{2,1});	
//			extra_ds.push_back(vector<int>{1,2});	
//			extra_ds.push_back(vector<int>{0,2});	
//			extra_ds.push_back(vector<int>{-1,2});	
//			extra_ds.push_back(vector<int>{-2,1});	
//			extra_ds.push_back(vector<int>{-2,0});	
//			extra_ds.push_back(vector<int>{-2,-1});	
//			extra_ds.push_back(vector<int>{-1,-2});	
//			extra_ds.push_back(vector<int>{0,-2});	
//			extra_ds.push_back(vector<int>{1,-2});	
//			extra_ds.push_back(vector<int>{2,-1});	
		}else if(tempNe==32){
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
		}else if(tempNe==57){
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
		}else{
			cout<<"not set up to deal with "<<tempNe<<" electrons"<<endl;
			exit(0);
		}
	}
	int nds=extra_ds.size();

	//ll is the object we will do monte carlo on, pp is the object with the electrons (or holes) shifted by one space
	int dsteps=nds; //nds
    vector<LATTICE> ll(invNu), pp(invNu);
    for (int i=0; i<invNu; i++) {
        ll[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
        pp[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
    }

    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
    vector<vector<Eigen::MatrixXcd > > overlaps( dsteps, vector<Eigen::MatrixXcd>(4, Eigen::MatrixXcd::Zero(invNu,invNu) ) );
	complex<double> temp;
    double energy;
    int dKx,dKy;
    for(int b=0; b<dsteps; b++) {
    	new_ds_ll=old_ds;
		new_ds_pp=old_ds;
        
		//depending on the mode, this adds one or two electons just outside the Fermi surface
		if(!holes){
			new_ds_ll.push_back(extra_ds[b]);
			new_ds_pp.push_back(extra_ds[supermod(b+1,nds)]);
			
			dKx=(extra_ds[supermod(b+1,nds)][0]-extra_ds[b][0]);
			dKy=(extra_ds[supermod(b+1,nds)][1]-extra_ds[b][1]);
			if(type=="twod"){
				dKx=0; dKy=0;
				new_ds_ll.push_back(extra_ds[supermod(b+nds/2,nds)]);
				new_ds_pp.push_back(extra_ds[supermod(b+nds/2+1,nds)]);
			}
		//this removes one or two electrons from the list of ds, if we are doing holes
		}else{
			new_ds_ll.erase(remove(new_ds_ll.begin(),new_ds_ll.end(),extra_ds[b]),new_ds_ll.end());
			new_ds_pp.erase(remove(new_ds_pp.begin(),new_ds_pp.end(),extra_ds[supermod(b+1,nds)]),new_ds_pp.end());
			dKx=-2*(extra_ds[supermod(b+1,nds)][0]-extra_ds[b][0]);
			dKy=-2*(extra_ds[supermod(b+1,nds)][1]-extra_ds[b][1]);
			if(type=="mtwod"){
				dKx=0; dKy=0;
				new_ds_ll.erase(remove(new_ds_ll.begin(),new_ds_ll.end(),extra_ds[supermod(b+nds/2,nds)]),new_ds_ll.end());
				new_ds_pp.erase(remove(new_ds_pp.begin(),new_ds_pp.end(),extra_ds[supermod(b+nds/2+1,nds)]),new_ds_pp.end());
			}				
        }
        
        for (int i=0; i<invNu; i++) {
            ll[i].set_ds(new_ds_ll);
            pp[i].set_ds(new_ds_pp);
            ll[i].reset(); pp[i].reset();
            ll[i].step(nWarmup);
        }
		cout<<"warmed up"<<endl;      
        pp[0].print_ds();
        energy=0;
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<invNu; i++) {
                ll[i].step(nSteps);
            }
            energy+=ll[0].coulomb_energy();
            for (int i=0; i<invNu; i++) {
                for (int j=0; j<invNu; j++) {
//                	if(i!=j){
//                		overlaps[b][0](i,j)+=0;
//                		overlaps[b][1](i,j)+=1;
//                		continue;
//                	}
                    temp=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[b][0](i,j)+=temp*ll[i].rhoq(dKx,dKy,ll[i].get_locs());//why times rhoq???
                    overlaps[b][1](i,j)+=norm(temp);
                    temp=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[b][2](i,j)+=temp;
                    overlaps[b][3](i,j)+=norm(temp);
                }
            }
        }
        
        for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
        overlaps[b][0]=overlaps[b][0].array()/overlaps[b][1].array().sqrt();
        overlaps[b][2]=overlaps[b][2].array()/overlaps[b][3].array().sqrt();
        hermitianize(overlaps[b][2]);
        cout<<"energy: "<<energy/(1.*nMeas*Ne)<<endl;
    }
    
    
    vector<Eigen::MatrixXcd> berrymatrix_step(dsteps);
    for (int b=0; b<dsteps; b++) berrymatrix_step[b] = overlaps[b][2].inverse() * overlaps[b][0];
    
    Eigen::MatrixXcd berrymatrix_integral = Eigen::MatrixXcd::Identity(invNu, invNu);
    vector<double> phases(invNu, 0.);
    datas.clear();
    for (int b=0; b<dsteps; b++) {
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
    cout<<"\n Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<1<<endl;
    cout<<"phase sum = "; for (int i=0; i<invNu; i++) {cout<<phases[i]<<" "; avephase+=phases[i]/(1.*invNu);} cout<<"\nphase average = "<<avephase<<endl;
    cout<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
    cout<<"amp(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<abs(es.eigenvalues()[i])<<" "; cout<<endl;
    cout<<"arg(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<arg(es.eigenvalues()[i])<<" ";cout<<endl;
    avephase=0.; for (int i=0; i<invNu; i++) avephase+=arg(es.eigenvalues()[i])/(1.*invNu); cout<<"ave arg(berrymatrix_integral.eigenvalue) = "<<avephase<<endl;
    cout<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
    cout<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
    cout<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl;
}

void CFL_berry_phases_parallel(vector<data> &datas, string params_name, string output_name, int num_core){
    ofstream outfile(output_name.c_str());
    string outfile2name=output_name+"_Mathmatica";//output data in format easy for Mathematica to deal with.
    ofstream outfile2(outfile2name.c_str());
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
    
    if(holes){
        if(tempNe==21){
//            extra_ds.push_back(vector<int>{1,0});
//            extra_ds.push_back(vector<int>{0,1});
//            extra_ds.push_back(vector<int>{-1,0});
//            extra_ds.push_back(vector<int>{0,-1});
            
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
            extra_ds.push_back(vector<int>{1,1});
//            extra_ds.push_back(vector<int>{1,0});
            extra_ds.push_back(vector<int>{1,-1});
//            extra_ds.push_back(vector<int>{0,-1});
            extra_ds.push_back(vector<int>{-1,-1});
//            extra_ds.push_back(vector<int>{-1,0});
            extra_ds.push_back(vector<int>{-1,1});
//            extra_ds.push_back(vector<int>{0,1});
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
    
    //ll is the object we will do monte carlo on, pp is the object with the electrons (or holes) shifted by one space
    int dsteps=nds; //nds
    //num_core threads to parallel the code.
    omp_set_num_threads(num_core);
    
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(invNu)), pp(num_core, vector<LATTICE>(invNu));//do this to avoid accessing wrong memory since openmp shares memory.
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
                        temp=ll[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                        overlaps[b][2](i,j)+=temp;// <ll|ll>
                        overlaps[b][3](i,j)+=norm(temp);
                    }
                }
            }
            
            for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
            overlaps[b][0]=overlaps[b][0].array()/overlaps[b][1].array().sqrt();
            overlaps[b][2]=overlaps[b][2].array()/overlaps[b][3].array().sqrt();
            hermitianize(overlaps[b][2]);
            //        cout<<"energy: "<<energy/(1.*nMeas*Ne)<<endl;
        }
        //parallel programming end.
        
        vector<Eigen::MatrixXcd> berrymatrix_step(dsteps);
        for (int b=0; b<dsteps; b++) berrymatrix_step[b] = overlaps[b][2].inverse() * overlaps[b][0];
        
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
    }
    outfile.close();
    outfile2.close();
}

//calculate ne=20 (one hole) state berry, 3 loop configurations available: fermisurface, surround0, notsurround0. See extra_ds.
void CFL_ne20hole_berryphase(vector<data> &datas, string params_name, string out_name, int num_core, string loop){
    string output_name=out_name+"_"+loop;
    ofstream outfile(output_name.c_str());
    string outfile2name=output_name+"_Mathmatica";//output data in format easy for Mathematica to deal with.
    ofstream outfile2(outfile2name.c_str());
    int supermod(int k, int n);
    int tempNe,Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
    bool testing;
    string type;
    
    //get inputs from the params file
    ifstream infile(params_name.c_str());
    //number of electrons, and inverse of filling fraction
    infile>>tempNe>>invNu;
    infile>>nWarmup>>nMeas>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    
    if (tempNe!=21) {
        cout<<"tempNe!=21"<<endl;
        exit(0);
    }
    Ne=tempNe-1;
    
    //this instance of LATTICE is only to set up the circular fermi surface of tempNe electrons
    LATTICE templl(tempNe, invNu, testing, "CFL", seed, 0);
    //old_ds is the maximal circular fermi surface with tempNe electrons
    vector<vector<int> > old_ds=templl.get_ds(), extra_ds;
    //old_dbar is the center of the circular fermi surface
    vector<double> old_dbar=templl.get_dbar_parameter();
    templl.print_ds();
    
    if (loop=="fermisurface") {
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
    else if (loop=="surround0") {
        extra_ds.push_back(vector<int>{1,0});
        extra_ds.push_back(vector<int>{0,1});
        extra_ds.push_back(vector<int>{-1,0});
        extra_ds.push_back(vector<int>{0,-1});
    }
    else if (loop=="notsurround0") {
        extra_ds.push_back(vector<int>{1,-1});
        extra_ds.push_back(vector<int>{0,-1});
        extra_ds.push_back(vector<int>{-1,-1});
        extra_ds.push_back(vector<int>{-1,-2});
        extra_ds.push_back(vector<int>{0,-2});
        extra_ds.push_back(vector<int>{1,-2});
    }
    else if (loop=="loop4") {
        extra_ds.push_back(vector<int>{0,-1});
        extra_ds.push_back(vector<int>{0,-2});
        extra_ds.push_back(vector<int>{1,-1});
    }
    else if (loop=="loop5"){
        extra_ds.push_back(vector<int>{1,-1});
        extra_ds.push_back(vector<int>{0,-1});
        extra_ds.push_back(vector<int>{1,-2});
    }
    else if (loop=="loop6"){
        extra_ds.push_back(vector<int>{2,0});
        extra_ds.push_back(vector<int>{1,0});
        extra_ds.push_back(vector<int>{2,-1});
    }
    else{
        cout<<"unrecognized loop"<<endl;
        exit(0);
    }
    
    int nds=extra_ds.size();
    
    //ll is the object we will do monte carlo on, pp is the object with the electrons (or holes) shifted by one space
    int dsteps=nds; //nds
    //num_core threads to parallel the code.
    omp_set_num_threads(num_core);
    
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(invNu)), pp(num_core, vector<LATTICE>(invNu));//do this to avoid accessing wrong memory since openmp shares memory.
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
            new_ds_ll.erase(remove(new_ds_ll.begin(),new_ds_ll.end(),extra_ds[b]),new_ds_ll.end());
            new_ds_pp.erase(remove(new_ds_pp.begin(),new_ds_pp.end(),extra_ds[supermod(b+1,nds)]),new_ds_pp.end());
            
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
                        temp=ll[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                        overlaps[b][2](i,j)+=temp;// <ll|ll>
                        overlaps[b][3](i,j)+=norm(temp);
                    }
                }
            }
            
            for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
            overlaps[b][0]=overlaps[b][0].array()/overlaps[b][1].array().sqrt();
            overlaps[b][2]=overlaps[b][2].array()/overlaps[b][3].array().sqrt();
            hermitianize(overlaps[b][2]);
            //        cout<<"energy: "<<energy/(1.*nMeas*Ne)<<endl;
        }
        //parallel programming end.
        
        vector<Eigen::MatrixXcd> berrymatrix_step(dsteps);
        for (int b=0; b<dsteps; b++) berrymatrix_step[b] = overlaps[b][2].inverse() * overlaps[b][0];
        
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
    }
    outfile.close();
    outfile2.close();
}

void laughlin_bp_single_state(int gs, vector<double> length, double steplength, int change_nMeas, vector<data> &datas, int num_core){
    int Ne,invNu,nWarmup,nMeas_t,nSteps,nBins,seed;
    bool testing;
    string type;
    ifstream infile("params");
    infile>>Ne>>invNu;
    infile>>nWarmup>>nMeas_t>>nSteps>>nBins;
    infile>>seed;
    infile>>testing;
    infile>>type;
    //initialize MC object
    
//    ofstream outfile(output_name.c_str());
    
    int nMeas;
    if (change_nMeas==0) nMeas=nMeas_t;
    else nMeas=change_nMeas;
    
    
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
    
    vector<LATTICE> ll(num_core), pp(num_core);
    for (int k=0; k<num_core; k++) {ll[k]=LATTICE(Ne, invNu, testing, type, seed, gs); pp[k]=LATTICE(Ne, invNu, testing, type, seed, gs);}
    
    vector<vector<complex<double> > > overlaps(nds, vector<complex<double>>(2));
    
    omp_set_num_threads(num_core);
#pragma omp parallel for
    for(int b=0;b<nds;b++){
        int coren=omp_get_thread_num();
        ll[coren].set_hole(holes[b]); pp[coren].set_hole(holes2[b]);
        ll[coren].reset(); ll[coren].step(nWarmup);
        pp[coren].reset();
        
        for(int k=0;k<nMeas;k++){
            ll[coren].step(nSteps);
            complex<double> temp=pp[coren].get_wf(ll[coren].get_locs())/ll[coren].get_wf(ll[coren].get_locs());
            overlaps[b][0]+=temp; overlaps[b][1]+=norm(temp);
        }
        for (int l=0; l<2; l++) overlaps[b][l]/=(1.*nMeas);
        
        overlaps[b][0]/=sqrt(abs(overlaps[b][1]));
    }
    
    datas.clear();
    for (int b=0; b<nds; b++) {
        data tmp;
        tmp.position[0]=holes[b][0]; tmp.position[1]=holes[b][1];
        tmp.amp[gs]=abs(overlaps[b][0]); tmp.ang[gs]=arg(overlaps[b][0]);
        datas.push_back(tmp);
//        outfile<<tmp.position[0]<<" "<<tmp.position[1]<<" "<<tmp.amp[gs]<<" "<<tmp.ang[gs]<<endl;
    }
//    outfile.close();
    
}

void check_orthogonality(string type){
    cout<<"--->type = "<<type<<", nMeas=50000, nStep=20, Ne=2. General d, sum of d neq 0."<<endl;
    int Ne=2, invNu, nWarmup=5000, nMeas=50000, nSteps=30, seed=0;
    if (type=="CFL") invNu=2;
    if (type=="laughlin") invNu=3;
    vector<LATTICE> cfl;
    for (int i=0; i<invNu; i++) {
        cfl.push_back(LATTICE(Ne, invNu, false, "CFL", 0, i));
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

void phase_variance(){
    ofstream outfile("M");
    int Ne=8,invNu=2,nWarmup=5000,nMeas=100,nSteps=20,nBins=1000,seed=0;
    outfile<<"Ne=8, invNu=2, nWarmup="<<nWarmup<<", nMeas="<<nMeas<<", nSteps="<<nSteps<<", nBins="<<nBins<<endl;
    bool testing=false;
    int Nphi=Ne*invNu;
    //initialize MC object
    
    //this instance of LATTICE is only to set up the circular fermi surface of tempNe electrons
    LATTICE templl(Ne+1, invNu, testing, "CFL", seed, 0);
    vector<vector<int> > old_ds=templl.get_ds(), extra_ds;
    //old_dbar is the center of the circular fermi surface
    vector<double> old_dbar=templl.get_dbar_parameter();
    
    extra_ds.push_back(vector<int>{1,-1});
    extra_ds.push_back(vector<int>{1,1});
    
    vector<vector<int>> ds0=old_ds, ds1=old_ds;
    ds0.erase(remove(ds0.begin(), ds0.end(), extra_ds[0]), ds0.end());//missing {1,-1}.
    ds1.erase(remove(ds1.begin(), ds1.end(), extra_ds[1]), ds1.end());//missing {1,1}.
    
    vector<LATTICE> ll(invNu), pp(invNu);
    for (int i=0; i<invNu; i++) {
        ll[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
        ll[i].set_ds(ds0);
        pp[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
        pp[i].set_ds(ds1);
    }
    
    //many body K.
    for (int i=0; i<invNu; i++) {
        for (int j=0; j<2; j++) {
            if (ll[i].dsum[j]%invNu!=0 || pp[i].dsum[j]%invNu!=0) {
                cout<<"dsum mod invNu != 0, somewhere wrong!"<<endl;
                exit(0);
            }
        }
    }
    int dKx=ll[0].dsum[0]/invNu-pp[0].dsum[0]/invNu, dKy=ll[0].dsum[1]/invNu-pp[0].dsum[1]/invNu;
    //alphabar = K1 L1/Nphi + K2 L2/Nphi = d/invNu. (K1,K2) are many body momentums.
    //So for d = d1 L1/Ne + d2 L2/Ne => K1=sum d1, K2=sum d2.
    //dKx/y is divided by invNu, because in lattice.cpp dsum is defined on L/Nphi lattice.
    
    //Monte Carlo Part & Output.
    for (unsigned nbin=0; nbin<nBins; nbin++) {
        //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
        vector<Eigen::MatrixXcd> overlaps(4, Eigen::MatrixXcd::Zero(invNu, invNu));
        
        for (int i=0; i<invNu; i++) {
            ll[i].reset();
            pp[i].reset();
            ll[i].step(nWarmup);
        }
        
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<invNu; i++)
                ll[i].step(nSteps);
            
            for (int i=0; i<invNu; i++) {
                for (int j=0; j<invNu; j++) {
                    complex<double> temp;
                    temp=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[0](i,j)+=temp*ll[i].rhoq(dKx,dKy,ll[i].get_locs());
                    overlaps[1](i,j)+=norm(temp);
                    temp=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[2](i,j)+=temp;
                    overlaps[3](i,j)+=norm(temp);
                }
            }
        }
        
        for (int l=0; l<4; l++) overlaps[l]/=(1.*nMeas);
        overlaps[0]=overlaps[0].array()/overlaps[1].array().sqrt();
        overlaps[2]=overlaps[2].array()/overlaps[3].array().sqrt();
        hermitianize(overlaps[2]);
        
        //berry matrix.
        Eigen::MatrixXcd berrymatrix=overlaps[2].inverse() * overlaps[0];
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix);
        outfile<<nbin<<" "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[1])<<endl;
        outfile<<nbin<<" "<<berrymatrix(0,0).real()<<" "<<berrymatrix(0,0).imag()<<" "<<berrymatrix(1,1).real()<<" "<<berrymatrix(1,1).imag()<<" "<<berrymatrix(0,1).real()<<" "<<berrymatrix(0,1).imag()<<" "<<berrymatrix(1,0).real()<<" "<<berrymatrix(1,0).imag()<<endl;
        outfile<<endl;
    }
    outfile.close();
}

void onestep(int ne, string output_name){
    ofstream outfile(output_name.c_str());
    int Ne=ne,invNu=2,nWarmup=5000,nMeas=100,nSteps=20,nBins=100,seed=0,Nphi=Ne*invNu;
    outfile<<"Ne="<<Ne<<", invNu=2, nWarmup="<<nWarmup<<", nMeas="<<nMeas<<", nSteps="<<nSteps<<", nBins="<<nBins<<endl;
    bool testing=false;
    //initialize MC object
    
    //set ds.
    vector<vector<vector<int>>> ds(2);
    if (Ne==2) {
        ds[0].push_back(vector<int>{0, 0}); ds[0].push_back(vector<int>{1, 0});
        ds[1].push_back(vector<int>{0, 0}); ds[1].push_back(vector<int>{0, 1});
    }
    if (Ne==3) {
        ds[0].push_back(vector<int>{ 1, 0});
        ds[0].push_back(vector<int>{-1, 0});
        ds[0].push_back(vector<int>{ 0, 1});
        ds[1].push_back(vector<int>{ 1, 0});
        ds[1].push_back(vector<int>{ 0, 1});
        ds[1].push_back(vector<int>{ 0, 0});
    }
    if (Ne==4) {
        ds[0].push_back(vector<int>{0, 0});
        ds[0].push_back(vector<int>{0, 1});
        ds[0].push_back(vector<int>{0, -1});
        ds[1].push_back(vector<int>{0, 0});
        ds[1].push_back(vector<int>{0, 1});
        ds[1].push_back(vector<int>{0, -1});
        
        ds[0].push_back(vector<int>{-1, 0});
        ds[1].push_back(vector<int>{-1, 1});
    }
    if (Ne==20) {
        LATTICE templl(21, 2, false, "CFL", 0, 0);
        vector<vector<int> > old_ds=templl.get_ds();
        
        ds[0]=old_ds; ds[1]=old_ds;
        ds[0].erase(remove(ds[0].begin(),ds[0].end(),vector<int>{0,-1}),ds[0].end());
        ds[1].erase(remove(ds[1].begin(),ds[1].end(),vector<int>{0,-2}),ds[1].end());
    }
    
    //set cfl wfs.
    vector<LATTICE> ll(invNu), pp(invNu);
    for (int i=0; i<invNu; i++) {
        ll[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
        ll[i].set_ds(ds[1]);
        pp[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
        pp[i].set_ds(ds[0]);
    }
    //many body K.
    for (int i=0; i<invNu; i++) {
        for (int j=0; j<2; j++) {
            if (ll[i].dsum[j]%invNu!=0 || pp[i].dsum[j]%invNu!=0) {
                cout<<"dsum mod invNu != 0, somewhere wrong!"<<endl;
                exit(0);
            }
        }
    }
    int dKx=ll[0].dsum[0]/invNu-pp[0].dsum[0]/invNu, dKy=ll[0].dsum[1]/invNu-pp[0].dsum[1]/invNu;
    //alphabar = K1 L1/Nphi + K2 L2/Nphi = d/invNu. (K1,K2) are many body momentums.
    //So for d = d1 L1/Ne + d2 L2/Ne => K1=sum d1, K2=sum d2.
    //dKx/y is divided by invNu, because in lattice.cpp dsum is defined on L/Nphi lattice.
    
    //Monte Carlo Part & Output.
    for (unsigned nbin=0; nbin<nBins; nbin++) {
        //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
        vector<Eigen::MatrixXcd> overlaps(4, Eigen::MatrixXcd::Zero(invNu, invNu));
        vector<Eigen::MatrixXcd> overlaps_density(4, Eigen::MatrixXcd::Zero(invNu, invNu));
        
        for (int i=0; i<invNu; i++) {
            ll[i].reset();
            pp[i].reset();
            ll[i].step(nWarmup);
        }
        
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<invNu; i++)
                ll[i].step(nSteps);
        
            for (int i=0; i<invNu; i++) {
                for (int j=0; j<invNu; j++) {
                    complex<double> temp;
                    //overlap.
                    temp=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[0](i,j)+=temp;//excluding density operator.
                    overlaps[1](i,j)+=norm(temp);
                    temp=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[2](i,j)+=temp;
                    overlaps[3](i,j)+=norm(temp);
                    //overlap_density.
                    temp=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps_density[0](i,j)+=temp*ll[i].rhoq(dKx,dKy,ll[i].get_locs());//including density operator.
                    overlaps_density[1](i,j)+=norm(temp);
                    temp=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps_density[2](i,j)+=temp;
                    overlaps_density[3](i,j)+=norm(temp);
                }
            }
        }
        
        for (int l=0; l<4; l++) overlaps[l]/=(1.*nMeas);
        overlaps[0]=overlaps[0].array()/overlaps[1].array().sqrt();// <ll|pp>
        overlaps[2]=overlaps[2].array()/overlaps[3].array().sqrt();// <ll|ll>
        hermitianize(overlaps[2]);
        
        for (int l=0; l<4; l++) overlaps_density[l]/=(1.*nMeas);
        overlaps_density[0]=overlaps_density[0].array()/overlaps_density[1].array().sqrt();// <ll|rhoq|pp>
        overlaps_density[2]=overlaps_density[2].array()/overlaps_density[3].array().sqrt();// <ll|ll>
        hermitianize(overlaps_density[2]);
        
//        //berry matrix.
//        Eigen::MatrixXcd berrymatrix=overlaps[2].inverse() * overlaps[0];
//        Eigen::MatrixXcd berrymatrix_density=overlaps_density[2].inverse() * overlaps_density[0];
////        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix);
////        outfile<<nbin<<" "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[1])<<endl;
//        outfile<<nbin<<" "<<berrymatrix(0,0).real()<<" "<<berrymatrix(0,0).imag()<<" "<<berrymatrix(1,1).real()<<" "<<berrymatrix(1,1).imag()<<" "<<berrymatrix(0,1).real()<<" "<<berrymatrix(0,1).imag()<<" "<<berrymatrix(1,0).real()<<" "<<berrymatrix(1,0).imag()<<endl;
//        outfile<<nbin<<" "<<berrymatrix_density(0,0).real()<<" "<<berrymatrix_density(0,0).imag()<<" "<<berrymatrix_density(1,1).real()<<" "<<berrymatrix_density(1,1).imag()<<" "<<berrymatrix_density(0,1).real()<<" "<<berrymatrix_density(0,1).imag()<<" "<<berrymatrix_density(1,0).real()<<" "<<berrymatrix_density(1,0).imag()<<endl;
//        outfile<<endl;
        
        //outputs: <ll|pp> matrix, <ll|rhoq|pp>, <ll|ll>.
        outfile<<nbin<<" ";
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                outfile<<overlaps[0](i,j).real()<<" "<<overlaps[0](i,j).imag()<<" ";
            }
        }
        outfile<<endl;
        outfile<<nbin<<" ";
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                outfile<<overlaps_density[0](i,j).real()<<" "<<overlaps_density[0](i,j).imag()<<" ";
            }
        }
        outfile<<endl;
        outfile<<nbin<<" ";
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                outfile<<overlaps[2](i,j).real()<<" "<<overlaps[2](i,j).imag()<<" ";
            }
        }
        outfile<<endl;
    }
    
    outfile.close();
}
