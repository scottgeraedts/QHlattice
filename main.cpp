#include "lattice.h"
#include "berry_tests.h"

//#include "berry_phase.h"


int main(){
//    plot_CFL_coule_vsdbar(10);
    
//    berry_phase bp(20);
//    bp.two_full_braiding();

//    berry_phase bp(16);
//    bp.two_full_braiding();
    
//	void single_run();
//	single_run();

//    void two_holes_scott();
//    two_holes_scott();
//    string str;
//    void two_holes(string str, int nmeasurement, data&);//str="test" or "".
    
//    ofstream bout("test_may20");
//    for (int i=500; i<10000; i=i+50) {
//        data test;
//        two_holes("test", i, test);
//        bout<<i<<" "<<test.amp[0]<<" "<<test.amp[1]<<" "<<test.amp[2]<<" "<<test.ang[0]<<" "<<test.ang[1]<<" "<<test.ang[2]<<endl;
//    }
    
//    data test;
//    two_holes("", 0, test);
    
    
//    void laughlin_bp_single_state(int gs, vector<double> length, double steplength, vector<data> &datas);//gs = 0, 1, 2, labeling ground state.
//    vector<double> length(2); double steplength = 0.01;
////    vector<vector<data> > datas; for (int i=0; i<3; i++) {vector<data> tmp; datas.push_back(tmp);}
//    vector<data> datas;
////    vector<vector<data> > datass(3);
////    
//    length[0]=0.05; length[1]=0.05;
//    void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas);
//    laughlinberryphase(length, steplength, datas);
    
    
//    ofstream bout("output_modified_orthogonal");
//    for (int i=0; i<datas.size(); i++) {
//        for (int gs=0; gs<3; gs++) {
//            bout<<datas[i].num<<" "<<datas[i].position[0]<<" "<<datas[i].position[1]<<" "<<gs<<" "<<datas[i].amp[gs]<<" "<<datas[i].ang[gs]<<endl;
//        }
//    }
    
//    //use model w.f. rather than berry matrix.
//    int gs=0;
//    length[0]=0.5; length[1]=0.8;
//    for (int gs=0; gs<3; gs++) {
//        laughlin_bp_single_state(gs, length, steplength, datass[gs]);
//        double angle=0.;
//        for (int i=0; i<datass[gs].size(); i++) {
//            angle+=datass[gs][i].ang[gs];
//        }
//        cout<<" gs = "<<gs<<" total phase = "<<angle<<endl;
//    }
    
//    //model w.f.
//    int gs=0;
//    length[0]=0.5; length[1]=0.5;
//    ofstream bout("teststeplength");
//    for (double steplength=0.01; steplength<0.5; steplength=steplength+0.02) {
//        laughlin_bp_single_state(gs, length, steplength, datas);
//        double minamp = datas[0].amp[gs], maxamp = datas[0].amp[gs], ang=0.;
//        for (int i=0; i<datas.size(); i++) {
//            ang+=datas[i].ang[gs];
//            if (datas[i].amp[gs]>maxamp) {maxamp=datas[i].amp[gs];}
//            if (datas[i].amp[gs]<minamp) {minamp=datas[i].amp[gs];}
//        }
//        bout<<steplength<<" "<<maxamp<<" "<<minamp<<" "<<ang<<endl;
//    }
//    ofstream bout("berry_laughlin_single_state");
//    for (int i=0; i<datas.size(); i++) {
//        bout<<datas[i].position[0]<<" "<<datas[i].position[1]<<" "<<datas[i].amp[gs]<<" "<<datas[i].ang[gs]<<endl;
//    }
    

//    testeigen();
    
//    test_largesize();

	void CFL_berry_phases();
	CFL_berry_phases();
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

void laughlinberryphase(vector<double> length, double steplength, vector<data> &datas, string str){
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
    
    vector<LATTICE> ll(3), pp(3);
    for (int i=0; i<3; i++) {
        ll[i]=LATTICE(Ne, invNu, testing, type, seed, i);
        pp[i]=LATTICE(Ne, invNu, testing, type, seed, i);
    }
    
    vector<vector<Eigen::MatrixXcd > > overlaps;
    //overlaps[b][0]=<psi(xb)|psi(xb+1)>, overlaps[b][1]=<|<psi(xb)|psi(xb+1)>|^2>, overlaps[b][2]=<psi(xb)|psi(xb)>, overlaps[b][3]=<|<psi(xb)|psi(xb)>|^2>.
    for (int b=0; b<nds; b++) {
        vector<Eigen::MatrixXcd> aa;
        Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(3,3);
        for (int i=0; i<4; i++) {
            aa.push_back(a);
        }
        overlaps.push_back(aa);
    }
    
    for(int b=0; b<nds; b++) {
        for (int i=0; i<3; i++) {
            ll[i].set_hole(holes[b]);
            pp[i].set_hole(holes2[b]);
            ll[i].reset(); ll[i].step(nWarmup);
            pp[i].reset();
        }
        for (int k=0; k<nMeas; k++) {
            for (int i=0; i<3; i++) {
                ll[i].step(nSteps);
            }
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    vector<complex<double> > temp(2);
                    temp[0]=pp[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[b][0](i,j)+=temp[0];
                    overlaps[b][1](i,j)+=norm(temp[0]);
                    temp[1]=ll[j].get_wf(ll[i].get_locs())/ll[i].get_wf(ll[i].get_locs());
                    overlaps[b][2](i,j)+=temp[1];
                    overlaps[b][3](i,j)+=norm(temp[1]);
                }
            }
        }
        for (int l=0; l<4; l++) overlaps[b][l]/=(1.*nMeas);
        for (int i=0; i<3; i++) for (int j=0; j<3; j++) {overlaps[b][0](i,j)/=sqrt(abs(overlaps[b][1](i,j))); overlaps[b][2](i,j)/=sqrt(abs(overlaps[b][3](i,j)));}

        hermitianize(overlaps[b][2]);
    }

    
    //comment this paragraph if want orthogonal.
    vector<Eigen::Matrix3cd> alphas;
    for (int b=0; b<nds; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b][2]);
        Eigen::Matrix3cd alpha;
        for (int i=0; i<3; i++) {
            Eigen::VectorXcd V = es.eigenvectors().col(i);
            complex<double> temp = es.eigenvalues()[i]*V.squaredNorm();
            V/=sqrt(temp);
            alpha.col(i)=V;
        }
        alphas.push_back(alpha);
    }
    
    if (str == "takeinverse") for (int b=0; b<nds; b++) overlaps[b][0]=overlaps[b][2].inverse() * overlaps[b][0];
    
    
    Eigen::Matrix3cd berrymatrix_integral = Eigen::Matrix3cd::Identity(3,3);
    vector<double> phases(3,0.);
    datas.clear();
    for (int b=0; b<nds; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b][0]);
        data tmp;
        berrymatrix_integral*=overlaps[b][0];
        for (int i=0; i<3; i++) {
            phases[i]+=arg(es.eigenvalues()[i]);
            tmp.num = b;
            tmp.amp[i] = abs(es.eigenvalues()[i]);
            tmp.ang[i] = arg(es.eigenvalues()[i]);
        }
        datas.push_back(tmp);
    }
    
    cout<<"\nstr = "<<str<<endl;
    cout<<"phase sum = "<<phases[0]<<" "<<phases[1]<<" "<<phases[2]<<" phase average = "<<(phases[0]+phases[1]+phases[2])/3<<endl;
    
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_integral);
    cout<<"berrymatrix_integral\n"<<berrymatrix_integral<<" \narg(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
    cout<<"amp of berrymatrix_integral = "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[0])<<" "<<abs(es.eigenvalues()[0])<<endl;
    cout<<"phase of berrymatrix_integral = "<<arg(es.eigenvalues()[0])<<" "<<arg(es.eigenvalues()[1])<<" "<<arg(es.eigenvalues()[2])<<endl;
//    for (int b=0; b<10; b++) {
//        cout<<"b= "<<b<<" inverse\n"<<overlaps[b][2].inverse()<<endl;
//    }
    
    
    
    
    /*
//    ofstream output("berrymatrixstep");
//    output<<"Ne = "<<Ne<<" length[0] = "<<length[0]<<" length[1] = "<<length[1]<<endl;
    Eigen::Matrix3cd berrymatrix_integral = Eigen::Matrix3cd::Identity(3,3);
    datas.clear();
    Eigen::MatrixXcd berrymatrix_step;
    for (int b=0; b<nds; b++) {
        berrymatrix_step =  alphas[b]*alphas[b].adjoint();//for orthogonal.
//        Eigen::MatrixXcd berrymatrix_step = overlaps[b][2];//for non-orthogonal.
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step);
        
        data tmp; tmp.position[0] = holes[b][0]; tmp.position[1] = holes[b][1]; tmp.num=b;
        for (int gs=0; gs<3; gs++) {
            tmp.amp[gs] = abs(es.eigenvalues()[gs]); tmp.ang[gs] = arg(es.eigenvalues()[gs]);
        }
        cout<<"-----"<<holes[b][0]<<" "<<holes[b][1]<<endl;
        cout<<"overlap[0]\n"<<overlaps[b][0]<<endl;
        cout<<"overlap[1]\n"<<overlaps[b][1]<<endl;
        cout<<"overlap[2]\n"<<overlaps[b][2]<<endl;
        cout<<"overlap[3]\n"<<overlaps[b][3]<<endl;
        cout<<"berrymatrix_step:\n"<<berrymatrix_step<<endl;
        cout<<"alpha^+ . alpha\n"<<alphas[b].adjoint()*alphas[b]<<endl;
        datas.push_back(tmp);
        berrymatrix_integral *= berrymatrix_step * overlaps[b][0];
//        berrymatrix_integral *=  overlaps[b][0];
//        output<<"b = "<<b<<"berrymatrix_step = \n"<<berrymatrix_step<<endl;
    }
    cout<<"berrymatrix_integral:\n"<<berrymatrix_integral<<endl;
    cout<<"\ntrace of total berry matrix = "<<berrymatrix_integral.trace()<<" ang(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
//    cout<<"\nberrymatrix_integral = \n"<<berrymatrix_integral<<endl;
    
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es_integral(berrymatrix_integral);
    for(int i=0;i<3;i++) cout<<arg(es_integral.eigenvalues()(i))<<endl;
//    cout<<abs(es_integral.eigenvalues()[0])<<" "<<arg(es_integral.eigenvalues()[0])<<" "<<abs(es_integral.eigenvalues()[1])<<" "<<arg(es_integral.eigenvalues()[1])<<" "<<abs(es_integral.eigenvalues()[2])<<" "<<arg(es_integral.eigenvalues()[2])<<endl;
    double totalberryphase=0.; for (int i=0; i<3; i++) totalberryphase+=arg(es_integral.eigenvalues()[i])/3.;
    cout<<"total berry phase = "<<totalberryphase<<endl;
    
//    Eigen::MatrixXcd test;
//    for (int i=0; i<10; i++) {
//        test = alphas[i].adjoint() * overlaps[i][0] * alphas[i];
//        cout<<" \n i = "<<i<<" test = \n"<<chop(test)<<endl;//expect to print out identity matrix.
//        cout<<"\n i = "<<i<<" alpha.adjoint * alpha = \n"<<alphas[i].adjoint()*alphas[i]<<endl;
//    }
    */
    
}

void CFL_berry_phases(){
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
	LATTICE templl(tempNe,invNu,testing,"CFL",seed,0); 
	vector<vector <int> > old_ds=templl.get_ds(), new_ds_ll,new_ds_pp, extra_ds;
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
                    overlaps[b][0](i,j)+=temp*ll[i].rhoq(dKx,dKy,ll[i].get_locs());
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

    //compensate for vectors not being orthogonal (maybe not necessary)
    vector<Eigen::MatrixXcd> alphas( dsteps, Eigen::MatrixXcd(invNu,invNu) );
    for (int b=0; b<dsteps; b++) {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(overlaps[b][2]);
        for (int i=0; i<invNu; i++) {
            Eigen::VectorXcd V = es.eigenvectors().col(i);
            complex<double> temp = es.eigenvalues()[i]*V.squaredNorm();
            V/=sqrt(temp);
            alphas[b].col(i)=V;
        }
    }
    
    Eigen::MatrixXcd berrymatrix_integral = Eigen::MatrixXcd::Identity(invNu,invNu);
    Eigen::MatrixXcd berrymatrix_step;
    for (int b=0; b<dsteps; b++) {
        berrymatrix_step =  alphas[b]*alphas[b].adjoint();//for orthogonal.
//        Eigen::MatrixXcd berrymatrix_step = overlaps[b][2];//for non-orthogonal.
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(berrymatrix_step);
        
        cout<<"-----"<<extra_ds[b][0]<<" "<<extra_ds[b][1]<<endl;
        cout<<overlaps[b][0]<<endl;
        cout<<overlaps[b][1]<<endl;
        cout<<overlaps[b][2]<<endl;
        cout<<overlaps[b][3]<<endl;
        cout<<berrymatrix_step<<endl;
        cout<<alphas[b].adjoint()*alphas[b]<<endl;
        berrymatrix_integral *= berrymatrix_step * overlaps[b][0];
//        berrymatrix_integral *=  overlaps[b][0];
//        output<<"b = "<<b<<"berrymatrix_step = \n"<<berrymatrix_step<<endl;
    }
    cout<<berrymatrix_integral<<endl;
    cout<<"\ntrace of total berry matrix = "<<berrymatrix_integral.trace()<<" ang(trace) = "<<arg(berrymatrix_integral.trace())<<endl;
//    cout<<"\nberrymatrix_integral = \n"<<berrymatrix_integral<<endl;
    
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es_integral(berrymatrix_integral);
    for(int i=0;i<invNu;i++) cout<<arg(es_integral.eigenvalues()(i))<<endl;
//    cout<<abs(es_integral.eigenvalues()[0])<<" "<<arg(es_integral.eigenvalues()[0])<<" "<<abs(es_integral.eigenvalues()[1])<<" "<<arg(es_integral.eigenvalues()[1])<<" "<<abs(es_integral.eigenvalues()[2])<<" "<<arg(es_integral.eigenvalues()[2])<<endl;
    double totalberryphase=0.; for (int i=0; i<invNu; i++) totalberryphase+=arg(es_integral.eigenvalues()[i])/invNu;
    cout<<"total berry phase = "<<totalberryphase<<endl;
    
	
	
}

