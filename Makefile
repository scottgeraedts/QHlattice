# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/product/simple directory makefile.

# including other makefiles.
CC=g++ -g -fopenmp 
HERE = /home/sgeraedt/QHlattice
MYDIR = /home/sgeraedt/myClibrary/
CFLAGS=-std=c++11 -Wall -I$(HERE) -I$(MYDIR)
LIBS=  -lgfortran $(MYDIR)/utils.o
OBJECTS=main.o lattice.o berry_tests.o z_function_m.o wf_tools.o new_coulomb_m.o weir3.o lattice_wrapper.o

a.out: $(OBJECTS)
	$(CC) -O3 $(CFLAGS) -o a.out $(OBJECTS) $(LIBS)
	
clean:
	rm -f *~ *.o *.mod a.out

%.o:	%.cpp
	$(CC) -O3 $(CFLAGS) -c $<

%.o:	%.f90
	gfortran -O3 $(CFLAGS) -c $<
