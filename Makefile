CC=g++ -std=c++11
HERE = /Users/jiewang/Google\ Drive/jie_programs/QHlattice
CLIBRARY = /Users/jiewang/Google\ Drive/jie_programs/Scott_Clibrary
EIGEN = /Users/jiewang/Google\ Drive/jie_programs/eigen/
CFLAGS = -I$(HERE) -I$(CLIBRARY) -I$(EIGEN)
LIBS=-L/usr/local/lib/gcc/4.9 -lgfortran 

a.out: main.o lattice.o berry_tests.o lattice_wrapper.o z_function_m.o new_coulomb_m.o wf_tools.o
	$(CC) -O3 $(CFLAGS) -fopenmp -o result z_function_m.o new_coulomb_m.o wf_tools.o lattice.o berry_tests.o lattice_wrapper.o main.o $(LIBS)
	
clean:
	rm -f *~ *.o a.out

%.o:	%.cpp
	$(CC) -w -O3 -fopenmp $(CFLAGS) -c $<

%.o:	%.f90
	gfortran -O3 $(CFLAGS) -c $<
 