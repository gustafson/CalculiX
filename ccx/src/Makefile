CFLAGS=-march=native -g -O0 ## -flto
FFLAGS=-march=native -g -O0 ## -flto

SPOOLESLIB=-lspooles -lpthread
SPOOLESINC=-I/usr/include/spooles -I/usr/include/spooles/MT -DSPOOLES
ARPACK=-DARPACK
ARPACKLIB= -larpack
## BLAS=-L/home/pi/gustafson/local/lib -lopenblas
TAUCS= # Currently not linking taucs "-DTAUCS"
TAUCSLIB= # Currently not linking metis "-ltaucs -lmetis".
## only required if using TAUCS LAPACK=-L/opt/ATLAS/lib -llapack -lf77blas -lpthread -lcblas -latlas
## BLAS=-L/opt/ATLAS/lib -lf77blas -latlas -lpthread
## PARDISO=-DPARDISO
## PARDISOLIB=-L/home/pi/gustafson/local/lib -lpardiso -lgfortran -fopenmp -lpthread -lm
USE_MT=-DUSE_MT
USE_MPI= # -DCALCULIX_MPI -I/home/pi/gustafson/local/include/spooles/MPI -fopenmp

CFLAGSLOC = $(CFLAGS) -Wall $(SPOOLESINC) -DARCH="Linux" -DMATRIXSTORAGE $(ARPACK) $(TAUCS) $(PARDISO) $(USE_MT) $(USE_MPI)
FFLAGSLOC = $(FFLAGS) -Wall -fopenmp

CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGSLOC) -c $<
.f.o :
	$(FC) $(FFLAGSLOC) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.6.1.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

DIR=../../../SPOOLES.2.2

LIBS = $(SPOOLESLIB) $(ARPACKLIB) $(TAUCSLIB) $(PARDISOLIB) $(BLAS) $(LAPACK) -lm -lc $(LDFLAGS) -fopenmp

ccx_2.6.1: $(OCCXMAIN) ccx_2.6.1.a
	./date.pl; $(CC) $(CFLAGSLOC) -c ccx_2.6.1.c ; $(FC) -Wall $(FFLAGS) -o $@ $(OCCXMAIN) ccx_2.6.1.a $(LIBS)

ccx_2.6.1.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
