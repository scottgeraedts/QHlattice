# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/product/simple directory makefile.

# including other makefiles.
CC=g++
HERE = /home/sgeraedt/QHlattice
MYDIR = /home/sgeraedt/myClibrary/
CFLAGS=-Wall -I$(HERE) -I$(MYDIR)
LIBS=  -lgfortran 

a.out: main.o lattice.o z_function_m.o coulomb2_m.o wf_tools.o
	$(CC) -O3 $(CFLAGS) -o a.out z_function_m.o coulomb2_m.o wf_tools.o lattice.o main.o $(LIBS)
	
clean:
	rm -f *~ *.o a.out

%.o:	%.cpp
	$(CC) -O3 $(CFLAGS) -c $<

%.o:	%.f90
	gfortran -O3 $(CFLAGS) -c $<
