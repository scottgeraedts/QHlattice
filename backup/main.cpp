using namespace std;
#include "lattice.h"
//#include "berry_tests.h"

//#include "berry_phase.h"

struct data{
    int num;
    double position[2];
    double amp[3];
    double ang[3];
    double energy;
    double ang_trace;
    double det;
    double dfromnorm;
};
complex<double> chop(complex<double> input){
    double rea=real(input), ima=imag(input);
    if (rea<pow(10,-5)) {
        rea=0.;
    }
    if (ima<pow(10,-5)) {
        ima=0.;
    }
    return rea+complex<double>(0,1)*ima;
}
Eigen::MatrixXcd chop(Eigen::MatrixXcd mat){
    int m=mat.rows(), n=mat.cols();
    Eigen::MatrixXcd ret = Eigen::MatrixXcd(m,n);
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            ret(i,j)=chop(mat(i,j));
        }
    }
    return ret;
}
void hermitianize(Eigen::MatrixXcd &x){
    for(int i=0;i<x.rows();i++){
        for(int j=i+1;j<x.rows();j++){
            x(i,j)=0.5*(x(i,j)+conj(x(j,i)));
            x(j,i)=conj(x(i,j));
        }
    }
}



int main(){
//	void CFL_berry_phases();
//	CFL_berry_phases();
    
//    testeigen();
    
    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int change_nMeas, int change_Ne, int ncore);
//    void test_error(int ne, double loop, double steplength, int nMea, int ncore, string filename, string test);
    
    vector<double> length(2); vector<data> datas; double loop; double steplength; int ne;
    length[0]=loop; length[1]=loop;
    
    int ncore=10;
    
    /*
    loop=0.5; ne=10; steplength=0.01;
    length[0]=loop; length[1]=loop;
    ofstream bout("test_error_stl_ll_10_new");
    for (int nMeas=10000; nMeas<200000; nMeas+=10000) {
        laughlinberryphase(length, steplength, datas, nMeas, ne, ncore);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
    */
    
    
    loop=0.5; ne=4; steplength=0.01;
    length[0]=loop; length[1]=loop;
    ofstream dout("test_error_stl_ll_4_new");
    for (int nMeas=5000; nMeas<200000; nMeas+=5000) {
        laughlinberryphase(length, steplength, datas, nMeas, ne, ncore);
        vector<double> phase(3);
        for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
        dout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
    }
    
}

//void CFL_berry_phases(){
//	int tempNe,Ne,invNu,nWarmup,nMeas,nSteps,nBins,seed;
//	bool testing;
//	string type;
//
//	//get  inputs from the params file
//	ifstream infile("params");
//	//number of electrons, and inverse of filling fraction
//	infile>>tempNe>>invNu; 
//	//monte carlo parameters
//	infile>>nWarmup>>nMeas>>nSteps>>nBins; 
//	//random seed
//	infile>>seed; 
//	//if 1, will calculate the energy every step and print it.
//	//This is helpful for debugging but slows down the code so should be set to 0 if the code is working
//	/*IMPORTANT: every time you try a larger system size than anyone has tried before (the largest I have done is 57)
//	 *you should rerun with testing=1, and make sure all the columns in the output are the same.
//	 This is because increasing the size may cause floating point overflow errors
//	 */
//	infile>>testing;
//	//a string which chooses mode to run the code in. right now there are 4 choices:  
//	//  twod: put two electrons outside the circular fermi surface, move them both around
//	//  oned: move one electron outside the fermi surface
//	//  mtwod: removes two electrons from the circular fermi surface, i.e. adds two "holes", moves these holes around
//	//  oned: removes one electron from the fermi surface
//	infile>>type; 
//
//	//tempNe is the number of electrons in the circular part of the Fermi surface
//	//Ne is the total number of electrons including the extra electrons/holes
//	if(type=="twod") Ne=tempNe+2;
//	else if(type=="oned") Ne=tempNe+1;
//	else if(type=="moned") Ne=tempNe-1;
//	else if(type=="mtwod") Ne=tempNe-2;
//	else{
//		cout<<"unrecognized type"<<endl;
//		exit(0);
//	}
//	//holes==true if we are removing electrons, false otherwise
//	bool holes=false;
//	if(type[0]=='m') holes=true;
//	
//	//this instance of LATTICE is only to set up the circular fermi surface of tempNe electrons
//	LATTICE templl(tempNe,invNu,testing,"CFL",seed,0); 
//	vector<vector <int> > old_ds=templl.get_ds(), new_ds_ll,new_ds_pp, extra_ds;
//	//old_dbar is the center of the circular fermi surface
//	vector<double> old_dbar=templl.get_dbar_parameter();
//	
//	//this part of the code specifies all the grid points just outside the circle made up of tempNe electrons
//	//we will loop through all these positions and add electrons to them
//	//to go to larger sizes, it will be necessary to add more possible values of tempNe
//	//it may end up being more convenient to write code to automate this step
//	if(!holes){
//		if(tempNe==21){
//			extra_ds.push_back(vector<int>{3,0});	
//			extra_ds.push_back(vector<int>{3,1});	
//			extra_ds.push_back(vector<int>{2,2});	
//			extra_ds.push_back(vector<int>{1,3});	
//			extra_ds.push_back(vector<int>{0,3});	
//			extra_ds.push_back(vector<int>{-1,3});	
//			extra_ds.push_back(vector<int>{-2,2});	
//			extra_ds.push_back(vector<int>{-3,1});	
//			extra_ds.push_back(vector<int>{-3,0});	
//			extra_ds.push_back(vector<int>{-3,-1});	
//			extra_ds.push_back(vector<int>{-2,-2});	
//			extra_ds.push_back(vector<int>{-1,-3});	
//			extra_ds.push_back(vector<int>{0,-3});	
//			extra_ds.push_back(vector<int>{1,-3});	
//			extra_ds.push_back(vector<int>{2,-2});	
//			extra_ds.push_back(vector<int>{3,-1});	
//		}else if(tempNe==32){
//			extra_ds.push_back(vector<int>{4,0});	
//			extra_ds.push_back(vector<int>{4,1});	
//			extra_ds.push_back(vector<int>{4,2});	
//			extra_ds.push_back(vector<int>{3,3});	
//			extra_ds.push_back(vector<int>{2,4});	
//			extra_ds.push_back(vector<int>{1,4});	
//			extra_ds.push_back(vector<int>{0,4});	
//			extra_ds.push_back(vector<int>{-1,4});	
//			extra_ds.push_back(vector<int>{-2,3});	
//			extra_ds.push_back(vector<int>{-3,2});	
//			extra_ds.push_back(vector<int>{-3,1});	
//			extra_ds.push_back(vector<int>{-3,0});	
//			extra_ds.push_back(vector<int>{-3,-1});	
//			extra_ds.push_back(vector<int>{-2,-2});	
//			extra_ds.push_back(vector<int>{-1,-3});	
//			extra_ds.push_back(vector<int>{0,-3});	
//			extra_ds.push_back(vector<int>{1,-3});	
//			extra_ds.push_back(vector<int>{2,-3});	
//			extra_ds.push_back(vector<int>{3,-2});	
//			extra_ds.push_back(vector<int>{4,-1});	
//		}else if(tempNe==37){
//			extra_ds.push_back(vector<int>{4,-1});	
//			extra_ds.push_back(vector<int>{4,0});	
//			extra_ds.push_back(vector<int>{4,1});	
//			
//			extra_ds.push_back(vector<int>{3,2});	
//			extra_ds.push_back(vector<int>{2,3});	
//
//			extra_ds.push_back(vector<int>{1,4});	
//			extra_ds.push_back(vector<int>{0,4});	
//			extra_ds.push_back(vector<int>{-1,4});	
//
//			extra_ds.push_back(vector<int>{-2,3});	
//			extra_ds.push_back(vector<int>{-3,2});	
//
//			extra_ds.push_back(vector<int>{-4,1});	
//			extra_ds.push_back(vector<int>{-4,0});	
//			extra_ds.push_back(vector<int>{-4,-1});	
//
//			extra_ds.push_back(vector<int>{-3,-2});	
//			extra_ds.push_back(vector<int>{-2,-3});	
//
//			extra_ds.push_back(vector<int>{-1,-4});	
//			extra_ds.push_back(vector<int>{0,-4});	
//			extra_ds.push_back(vector<int>{1,-4});	
//
//			extra_ds.push_back(vector<int>{2,-3});	
//			extra_ds.push_back(vector<int>{3,-2});
//		}else if (tempNe==57){
//			extra_ds.push_back(vector<int>{5,-1});
//			extra_ds.push_back(vector<int>{5,0});
//			extra_ds.push_back(vector<int>{5,1});
//
//			extra_ds.push_back(vector<int>{4,2});
//			extra_ds.push_back(vector<int>{3,3});
//			extra_ds.push_back(vector<int>{2,4});
//			
//			extra_ds.push_back(vector<int>{1,5});
//			extra_ds.push_back(vector<int>{0,5});
//			extra_ds.push_back(vector<int>{-1,5});
//
//			extra_ds.push_back(vector<int>{-2,4});
//			extra_ds.push_back(vector<int>{-3,3});
//			extra_ds.push_back(vector<int>{-4,2});
//			
//			extra_ds.push_back(vector<int>{-5,1});
//			extra_ds.push_back(vector<int>{-5,0});
//			extra_ds.push_back(vector<int>{-5,-1});
//			
//			extra_ds.push_back(vector<int>{-4,-2});
//			extra_ds.push_back(vector<int>{-3,-3});
//			extra_ds.push_back(vector<int>{-2,-4});
//			
//			extra_ds.push_back(vector<int>{-1,-5});
//			extra_ds.push_back(vector<int>{0,-5});
//			extra_ds.push_back(vector<int>{1,-5});
//			
//			extra_ds.push_back(vector<int>{2,-4});
//			extra_ds.push_back(vector<int>{3,-3});
//			extra_ds.push_back(vector<int>{4,-2});									
//		}else{
//			cout<<"not set up to deal with "<<tempNe<<" electrons"<<endl;
//			exit(0);
//		}
//	//this code does the same thing as above, but it lists all the positions just inside the fermi surface, where electrons should be removed if we are doing holes
//	}else{
//		if(tempNe==21){
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
//		}else if(tempNe==32){
//			extra_ds.push_back(vector<int>{3,0});	
//			extra_ds.push_back(vector<int>{3,1});	
//			extra_ds.push_back(vector<int>{3,2});	
//			extra_ds.push_back(vector<int>{2,3});	
//			extra_ds.push_back(vector<int>{1,3});	
//			extra_ds.push_back(vector<int>{0,3});	
//			extra_ds.push_back(vector<int>{-1,3});	
//			extra_ds.push_back(vector<int>{-2,2});	
//			extra_ds.push_back(vector<int>{-2,1});	
//			extra_ds.push_back(vector<int>{-2,0});	
//			extra_ds.push_back(vector<int>{-2,-1});	
//			extra_ds.push_back(vector<int>{-1,-2});	
//			extra_ds.push_back(vector<int>{-0,-2});	
//			extra_ds.push_back(vector<int>{1,-2});	
//			extra_ds.push_back(vector<int>{2,-2});	
//			extra_ds.push_back(vector<int>{3,-1});	
//		}else if(tempNe==57){
//			extra_ds.push_back(vector<int>{4,-1});	
//			extra_ds.push_back(vector<int>{4,0});	
//			extra_ds.push_back(vector<int>{4,1});
//				
//			extra_ds.push_back(vector<int>{3,2});	
//			extra_ds.push_back(vector<int>{2,3});	
//
//			extra_ds.push_back(vector<int>{1,4});	
//			extra_ds.push_back(vector<int>{0,4});	
//			extra_ds.push_back(vector<int>{-1,4});	
//
//			extra_ds.push_back(vector<int>{-2,3});	
//			extra_ds.push_back(vector<int>{-3,2});	
//
//			extra_ds.push_back(vector<int>{-4,1});	
//			extra_ds.push_back(vector<int>{-4,0});	
//			extra_ds.push_back(vector<int>{-4,-1});	
//
//			extra_ds.push_back(vector<int>{-3,-2});	
//			extra_ds.push_back(vector<int>{-2,-3});	
//
//			extra_ds.push_back(vector<int>{-1,-4});	
//			extra_ds.push_back(vector<int>{0,-4});	
//			extra_ds.push_back(vector<int>{1,-4});	
//
//			extra_ds.push_back(vector<int>{2,-3});	
//			extra_ds.push_back(vector<int>{3,-2});	
//		}else{
//			cout<<"not set up to deal with "<<tempNe<<" electrons"<<endl;
//			exit(0);
//		}
//	}	
//	int nds=extra_ds.size();
//	//ll is the object we will do monte carlo on, pp is the object with the electrons (or holes) shifted by one space
//	int dsteps=4; //nds
//    vector<LATTICE> ll(invNu), pp(invNu);
//    for (int i=0; i<invNu; i++) {
//        ll[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
//        pp[i]=LATTICE(Ne, invNu, testing, "CFL", seed, i);
//    }
//    
//    vector<vector<Eigen::MatrixXcd > > overlaps( dsteps, vector<Eigen::MatrixXcd>(4, Eigen::MatrixXcd::Zero(invNu,invNu) ) );
//	complex<double> temp;
//    double energy;
//    for(int b=0; b<dsteps; b++) {
//    	new_ds_ll=old_ds;
//		new_ds_pp=old_ds;
//		//depending on the mode, this adds one or two electons just outside the Fermi surface
//		if(!holes){
//			new_ds_ll.push_back(extra_ds[b]);
//            int supermod(int k, int n);
//			new_ds_pp.push_back(extra_ds[supermod(b+1,nds)]);
//			
//			if(type=="twod"){
//				new_ds_ll.push_back(extra_ds[supermod(b+nds/2,nds)]);
//				new_ds_pp.push_back(extra_ds[supermod(b+nds/2+1,nds)]);
//			}
//		//this removes one or two electrons from the list of ds, if we are doing holes
//		}else{
//			new_ds_ll.erase(remove(new_ds_ll.begin(),new_ds_ll.end(),extra_ds[b]),new_ds_ll.end());
//			new_ds_pp.erase(remove(new_ds_pp.begin(),new_ds_pp.end(),extra_ds[supermod(b+1,nds)]),new_ds_pp.end());
//			if(type=="mtwod"){
//				new_ds_ll.erase(remove(new_ds_ll.begin(),new_ds_ll.end(),extra_ds[supermod(b+nds/2,nds)]),new_ds_ll.end());
//				new_ds_pp.erase(remove(new_ds_pp.begin(),new_ds_pp.end(),extra_ds[supermod(b+nds/2+1,nds)]),new_ds_pp.end());
//			}				
//		}	    		
//        for (int i=0; i<invNu; i++) {
//            ll[i].set_ds(new_ds_ll);
//            pp[i].set_ds(new_ds_pp);
//            ll[i].step(nWarmup);
//            
//        }
//		cout<<"warmed up"<<endl;      
//        pp[0].print_ds();
//        energy=0;
//        for (int k=0; k<nMeas; k++) {
//            for (int i=0; i<invNu; i++) {
//                ll[i].step(nSteps);
//            }
//            energy+=ll[0].coulomb_energy();
//            for (int i=0; i<invNu; i++) {
//                for (int j=0; j<invNu; j++) {
//                    temp=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
//                    overlaps[b][0](i,j)+=temp;
//                    overlaps[b][1](i,j)+=norm(temp);
//                    temp=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
//                    overlaps[b][2](i,j)+=temp;
//                    overlaps[b][3](i,j)+=norm(temp);
//                }
//            }
//        }
//        for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
//        overlaps[b][0]=overlaps[b][0].array()/overlaps[b][1].array().sqrt();
//        overlaps[b][2]=overlaps[b][2].array()/overlaps[b][3].array().sqrt();
//        hermitianize(overlaps[b][2]);
//        cout<<"energy: "<<energy/(1.*nMeas*Ne)<<endl;
//        if (b>=1) break;
//    }
//
//    //compensate for vectors not being orthogonal (maybe not necessary)
//    vector<Eigen::MatrixXcd> alphas( dsteps, Eigen::MatrixXcd(invNu,invNu) );
//    for (int b=0; b<nds; b++) {
//        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b][2]);
//        for (int i=0; i<invNu; i++) {
//            Eigen::VectorXcd V = es.eigenvectors().col(i);
//            complex<double> temp = es.eigenvalues()[i]*V.squaredNorm();
//            V/=sqrt(temp);
//            alphas[b].col(i)=V;
//        }
//    }
//    
//    Eigen::MatrixXcd berrymatrix_integral = Eigen::MatrixXcd::Identity(invNu,invNu);
//    Eigen::MatrixXcd berrymatrix_step;
//    for (int b=0; b<dsteps; b++) {
//        berrymatrix_step =  alphas[b]*alphas[b].adjoint();//for orthogonal.
////        Eigen::MatrixXcd berrymatrix_step = overlaps[b][2];//for non-orthogonal.
//        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step);
//        
//        cout<<"-----"<<extra_ds[b][0]<<" "<<extra_ds[b][1]<<endl;
//        cout<<overlaps[b][0]<<endl;
//        cout<<overlaps[b][1]<<endl;
//        cout<<overlaps[b][2]<<endl;
//        cout<<overlaps[b][3]<<endl;
//        cout<<berrymatrix_step<<endl;
//        cout<<alphas[b].adjoint()*alphas[b]<<endl;
//        berrymatrix_integral *= berrymatrix_step * overlaps[b][0];
////        berrymatrix_integral *=  overlaps[b][0];
////        output<<"b = "<<b<<"berrymatrix_step = \n"<<berrymatrix_step<<endl;
//    }
//    cout<<berrymatrix_integral<<endl;
//    cout<<"\ntrace of total berry matrix = "<<berrymatrix_integral.trace()<<" ang(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
////    cout<<"\nberrymatrix_integral = \n"<<berrymatrix_integral<<endl;
//    
//    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es_integral(berrymatrix_integral);
//    for(int i=0;i<invNu;i++) cout<<arg(es_integral.eigenvalues()(i))<<endl;
////    cout<<abs(es_integral.eigenvalues()[0])<<" "<<arg(es_integral.eigenvalues()[0])<<" "<<abs(es_integral.eigenvalues()[1])<<" "<<arg(es_integral.eigenvalues()[1])<<" "<<abs(es_integral.eigenvalues()[2])<<" "<<arg(es_integral.eigenvalues()[2])<<endl;
//    double totalberryphase=0.; for (int i=0; i<invNu; i++) totalberryphase+=arg(es_integral.eigenvalues()[i])/invNu;
//    cout<<"total berry phase = "<<totalberryphase<<endl;
//}



//---------- parallel programming code ------//
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
    
    //    vector<LATTICE> ll(3), pp(3);
    //    for (int i=0; i<3; i++) {ll[i]=LATTICE(Ne, invNu, testing, type, seed, i); pp[i]=LATTICE(Ne, invNu, testing, type, seed, i);}
    
    vector<vector<Eigen::MatrixXcd > > overlaps;
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
    for (int b=0; b<nds; b++) {
        vector<Eigen::MatrixXcd> aa;
        Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(3,3);
        for (int i=0; i<4; i++) aa.push_back(a);
        overlaps.push_back(aa);
    }
    
    omp_set_num_threads(num_core);
    vector<vector<LATTICE> > ll(num_core, vector<LATTICE>(3)), pp(num_core, vector<LATTICE>(3));//do this to avoid wrong memory access since openmp share memory.
    for (int k=0; k<num_core; k++) for (int i=0; i<3; i++) {ll[k][i]=LATTICE(Ne, invNu, testing, type, seed, i); pp[k][i]=LATTICE(Ne, invNu, testing, type, seed, i);}
    
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
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {overlaps[b][0](i,j)/=sqrt(abs(overlaps[b][1](i,j))); overlaps[b][2](i,j)/=sqrt(abs(overlaps[b][3](i,j)));}
        
        hermitianize(overlaps[b][2]);
    }
    //parallel programming end.

    vector<Eigen::MatrixXcd> berrymatrix_step(nds);
    for (int b=0; b<nds; b++) berrymatrix_step[b] = overlaps[b][2].inverse() * overlaps[b][0];
    
    Eigen::Matrix3cd berrymatrix_integral = Eigen::Matrix3cd::Identity(3,3);
    vector<double> phases(3,0.);
    datas.clear();
    for (int b=0; b<nds; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step[b]);
        data tmp;
        berrymatrix_integral *= berrymatrix_step[b];
        for (int i=0; i<3; i++) {
            phases[i]+=arg(es.eigenvalues()[i]);
            tmp.num = b; tmp.amp[i] = abs(es.eigenvalues()[i]); tmp.ang[i] = arg(es.eigenvalues()[i]);
        }
        
        //insert dfromnorm.
        double normeigenvalue=0., normmatrix=0.;
        for (int i=0; i<3; i++) {
            normeigenvalue+=sqrt(norm(es.eigenvalues()[i]));
        }
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                normmatrix+=sqrt(norm(berrymatrix_step[b](i,j)));
            }
        }
        tmp.dfromnorm=normmatrix-normeigenvalue;
        
        datas.push_back(tmp);
    }
    
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_integral);
    datas[0].ang_trace = arg(berrymatrix_integral.trace());
    datas[0].det = arg(berrymatrix_integral.determinant());
    cout<<"\n\n Ne="<<Ne<<" nMea="<<nMeas<<" nStep="<<nSteps<<" ncore="<<omp_get_num_threads()<<endl;
    cout<<"phase sum = "<<phases[0]<<" "<<phases[1]<<" "<<phases[2]<<"\nphase average = "<<(phases[0]+phases[1]+phases[2])/3<<endl;
    cout<<"berrymatrix_integral\n"<<berrymatrix_integral<<endl;
    cout<<"amp(berrymatrix_integral.eigenvalue) = "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[0])<<endl;
    cout<<"arg(berrymatrix_integral.eigenvalue) = "<<arg(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[2])<<endl;
    cout<<"sum arg(berrymatrix_integral.eigenvalue) = "<< ( arg(es.eigenvalues()[0])+arg(es.eigenvalues()[1])+arg(es.eigenvalues()[2]) )/3 <<endl;
    cout<<"arg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
    cout<<"amp(trace) = "<<abs(berrymatrix_integral.trace())<<endl;
    cout<<"arg(det) = "<<arg(berrymatrix_integral.determinant())<<endl;
}

//void test_error(int ne, double loop, double steplength, int nMea, int ncore, string filename, string test){
//    void laughlinberryphase_parallel(vector<double> length, double steplength, vector<data> &datas, int, int, int);
//    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, int, int);
//    vector<double> length(2); vector<data> datas;
//    length[0]=loop; length[1]=loop;
//    ofstream bout("test_error_stl_ll_10");
//    if (test=="steplength") {
//        for (double steplength=0.001; steplength<loop; steplength+=0.002) {
//            if (ncore>1) laughlinberryphase_parallel(length, steplength, datas, nMea, ne, ncore);
//            else laughlinberryphase(length, steplength, datas, nMea, ne);
//            vector<double> phase(3);
//            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
//            bout<<steplength<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
//        }
//    }
//    if (test=="ne") {
//        for (int ne=2; ne<22; ne=ne+2) {
//            if (ncore>1) laughlinberryphase_parallel(length, steplength, datas, nMea, ne, ncore);
//            else laughlinberryphase(length, steplength, datas, nMea, ne);
//            vector<double> phase(3);
//            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
//            bout<<ne<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
//        }
//    }
//    if (test=="nMea") {
//        for (int nMeas=10; nMeas<20; nMeas+=2){
//            if (ncore>1) laughlinberryphase_parallel(length, steplength, datas, nMea, ne, ncore);
//            else laughlinberryphase(length, steplength, datas, nMea, ne);
//            vector<double> phase(3);
//            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
//            bout<<nMeas<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
//        }
//    }
//    if (test=="loop") {
//        for (double x=0.05; x<0.8; x+=0.05) {
//            length[0]=x; length[1]=x;
//            if (ncore>1) laughlinberryphase_parallel(length, steplength, datas, nMea, ne, ncore);
//            else laughlinberryphase(length, steplength, datas, nMea, ne);
//            vector<double> phase(3);
//            for (int i=0; i<datas.size(); i++) for (int j=0; j<3; j++) phase[j]+=datas[i].ang[j];
//            bout<<x<<" "<<phase[0]<<" "<<phase[1]<<" "<<phase[2]<<endl;
//        }
//    }
//    if (test=="normality") {
//        for (double steplength=0.01; steplength<0.25; steplength+=0.01) {
//            if (ncore>1) laughlinberryphase_parallel(length, steplength, datas, 0, 0, ncore);
//            else laughlinberryphase(length, steplength, datas, 0, 0);
//            
//            bout<<steplength<<endl;
//            for (int i=0; i<datas.size(); i++) {
//                bout<<datas[i].num<<" "<<datas[i].dfromnorm<<endl;
//            }
//            bout<<endl;
//        }
//    }
//    
//}



