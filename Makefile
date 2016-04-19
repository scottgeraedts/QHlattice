# including other makefiles.
CC=icpc
HERE = /Users/jiewang/QHlattice
CLIBRARY=/Users/jiewang/myClibrary
EIGEN=/Users/jiewang/QHlattice/eigen-eigen-07105f7124f9/
CFLAGS=-Wall -I$(HERE) -I$(CLIBRARY) -I$(EIGEN)
LIBS=-L/usr/local/lib/gcc/4.9 -lgfortran 

a.out: main.o lattice.o z_function_m.o coulomb2_m.o wf_tools.o
	$(CC) -O3 $(CFLAGS) -o a.out z_function_m.o coulomb2_m.o wf_tools.o lattice.o main.o $(LIBS)
	
clean:
	rm -f *~ *.o a.out


%.o:	%.cpp
	$(CC) -O3 $(CFLAGS) -c $<

%.o:	%.f90
	gfortran -O3 $(CFLAGS) -c $<
