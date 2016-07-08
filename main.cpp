using namespace std;
//#include <omp.h>
#include "lattice.h"
#include "berry_phase.h"

int main(){
    void CFL_berry_phases(vector<data> &datas);
    void CFL_berry_phases_parallel(vector<data> &datas, int num_core);
    vector<data> datas; int num_core;
    num_core=1;
	CFL_berry_phases(datas);
    cout<<"----------"<<endl;
    CFL_berry_phases_parallel(datas, num_core);

//    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne, int num_core);
//    vector<double> length(2); double steplength; int nMeas, ne, ncore; vector<data> datas;
//    length[0]=0.1; length[1]=0.1; steplength=0.01; nMeas=0; ne=0; ncore=2;
//    laughlinberryphase(length, steplength, datas, nMeas, ne, ncore);
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
		if(tempNe==21){
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
		if(tempNe==21){
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
	int dsteps=3; //nds
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
			
			dKx=2*(extra_ds[supermod(b+1,nds)][0]-extra_ds[b][0]);
			dKy=2*(extra_ds[supermod(b+1,nds)][1]-extra_ds[b][1]);
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
    cout<<"\n\n Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<1<<endl;
    cout<<"phase sum = "; for (int i=0; i<invNu; i++) {cout<<phases[i]<<" "; avephase+=phases[i]/(1.*invNu);} cout<<"\nphase average = "<<avephase<<endl;
    cout<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
    cout<<"amp(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<abs(es.eigenvalues()[i])<<" "; cout<<endl;
    cout<<"arg(berrymatrix_integral.eigenvalue) = "; for (int i=0; i<invNu; i++) cout<<arg(es.eigenvalues()[i])<<" ";cout<<endl;
    avephase=0.; for (int i=0; i<invNu; i++) avephase+=arg(es.eigenvalues()[i])/(1.*invNu); cout<<"sum arg(berrymatrix_integral.eigenvalue) = "<<avephase<<endl;
    cout<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
    cout<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
    cout<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl;
    
}



void CFL_berry_phases_parallel(vector<data> &datas, int num_core){
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
    vector<vector<int> > old_ds=templl.get_ds(), extra_ds;
    //old_dbar is the center of the circular fermi surface
    vector<double> old_dbar=templl.get_dbar_parameter();
    
    //this part of the code specifies all the grid points just outside the circle made up of tempNe electrons
    //we will loop through all these positions and add electrons to them
    //to go to larger sizes, it will be necessary to add more possible values of tempNe
    //it may end up being more convenient to write code to automate this step
    if(!holes){
        if(tempNe==21){
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
        if(tempNe==21){
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
    int dsteps=3; //nds
//    vector<LATTICE> ll(invNu), pp(invNu);
//    for (int i=0; i<invNu; i++) {
//        ll[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
//        pp[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
//    }
    
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
    vector<vector<Eigen::MatrixXcd > > overlaps(nds, vector<Eigen::MatrixXcd>(4, Eigen::MatrixXcd::Zero(invNu,invNu) ) );
    
    omp_set_num_threads(num_core);
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(invNu)), pp(num_core, vector<LATTICE>(invNu));//do this to avoid wrong memory access since openmp share memory.
    for (int k=0; k<num_core; k++) for (int i=0; i<invNu; i++) {ll[k][i]=LATTICE(Ne, invNu, testing, "CFL", seed, i); pp[k][i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);}
    
    //parallel programming begin.
    vector<double> energy(dsteps);
#pragma omp parallel for
    for(int b=0; b<dsteps; b++) {
        
        int coren = omp_get_thread_num();
        
        int dKx,dKy;
        vector<vector<int> > new_ds_ll, new_ds_pp;
        new_ds_ll=old_ds;
        new_ds_pp=old_ds;
        
        //depending on the mode, this adds one or two electons just outside the Fermi surface
        if(!holes){
            new_ds_ll.push_back(extra_ds[b]);
            new_ds_pp.push_back(extra_ds[supermod(b+1,nds)]);
            
            dKx=2*(extra_ds[supermod(b+1,nds)][0]-extra_ds[b][0]);
            dKy=2*(extra_ds[supermod(b+1,nds)][1]-extra_ds[b][1]);
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
            ll[coren][i].set_ds(new_ds_ll);
            pp[coren][i].set_ds(new_ds_pp);
            ll[coren][i].reset(); pp[coren][i].reset();
            ll[coren][i].step(nWarmup);
        }
//        cout<<"warmed up"<<endl;
//        pp[0].print_ds();
//        energy=0;
        energy[b]=0.;
        
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<invNu; i++) ll[coren][i].step(nSteps);
            
            energy[b]+=ll[coren][0].coulomb_energy();
            
            for (int i=0; i<invNu; i++) {
                for (int j=0; j<invNu; j++) {
                    complex<double> temp;
                    temp=pp[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                    overlaps[b][0](i,j)+=temp*ll[coren][i].rhoq(dKx,dKy,ll[coren][i].get_locs());//why times rhoq???
                    overlaps[b][1](i,j)+=norm(temp);
                    temp=ll[coren][j].get_wf(ll[coren][i].get_locs())/ll[coren][i].get_wf(ll[coren][i].get_locs());
                    overlaps[b][2](i,j)+=temp;
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
    
    for (int b=0; b<dsteps; b++) {
        cout<<"b = "<<b<<" energy = "<<energy[b]/(1.*nMeas*Ne)<<endl;
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

