## This file is heavily modified from the Makefile distributed by with
## ccx. 

## pkg-config is used to pull all flags from local enviroment.  You
## must set this up.  Typically this would be defining
## PKG_CONFIG_PATH=$PKG_CONFIG_PATH:~/local/lib64/pkgconfig or
## similar.  You will also have to add .pc file for each library you
## use.  I did this because I have several machines with different
## library enviroments.  You could alternatively use compiler flags
## such as -L/path/to/lib -I/path/to/include -llib.  Examples of this
## on my distribution are in the comments below

CFLAGS=-march=native -g -O0 -pg ## debugging
FFLAGS=-march=native -g -O0 -pg ## debugging
CFLAGS=-march=native -O2 -Wall ## conservative
FFLAGS=-march=native -O2 -Wall ## conservative
CFLAGS=-march=native -O3 -Wall -flto ## aggressive
FFLAGS=-march=native -O3 -Wall -flto ## aggressive

## SPOOLES
## CFLAGS  += -I/usr/include/spooles -I/usr/include/spooles/MT -DSPOOLES
## LDFLAGS += -lspooles -lpthread 
CFLAGS  += `pkg-config --cflags spooles` -DSPOOLES
LDFLAGS += `pkg-config --libs spooles`

## ARPACK
## CFLAGS  += -DARPACK
## LDFLAGS += -larpack
CFLAGS  += `pkg-config --cflags arpack` -DARPACK
LDFLAGS += `pkg-config --libs arpack`

## TAUCS 
## CFLAGS  += -DTAUCS
## LDFLAGS += -ltaucs -lmetis

## LAPACK
## CFLAGS  += -I/usr/include/openblas
## LDFLAGS += -lreflapack -lopenblas
CFLAGS  += `pkg-config --cflags lapack`
LDFLAGS += `pkg-config --libs lapack`

## BLAS
## CFLAGS  += -I/usr/include/openblas
## LDFLAGS += -lopenblas
CFLAGS  += `pkg-config --cflags blas`
LDFLAGS += `pkg-config --libs blas`

## PARDISO
## CFLAGS  += -DPARDISO
## LDFLAGS += -L/home/pete/local/lib64/ -lpardiso -lgfortran -lpthread -lm -fopenmp
CFLAGS  += `pkg-config --cflags pardiso` -DPARDISO
LDFLAGS += `pkg-config --libs pardiso`

## Multi Threaded and MPI
USE_MT=-DUSE_MT
USE_MPI= #-DCALCULIX_MPI -fopenmp

## EXODUSII

## CUDACUSP

## Other CCX Options
CFLAGS  += -DARCH="Linux" -DMATRIXSTORAGE 
LDFLAGS += 

CC=gcc
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

## Link to math and standard c
CFLAGS += -lm -lc

ccx_2.6.1: $(OCCXMAIN) ccx_2.6.1.a
	## ./date.pl; 
	$(CC) $(CFLAGS) -c ccx_2.6.1.c $(LDFLAGS); $(FC) -Wall $(FFLAGS) -o $@ $(OCCXMAIN) ccx_2.6.1.a $(LDFLAGS)

ccx_2.6.1.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

clean:
	rm *.a *.o
