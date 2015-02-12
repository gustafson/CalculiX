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

## CFLAGS=-march=native -g -O0 -pg ## debugging
## FFLAGS=-march=native -g -O0 -pg ## debugging
## CFLAGS=-march=native -O2 -Wall ## conservative
## FFLAGS=-march=native -O2 -Wall ## conservative
CFLAGS=-march=native -O3 -Wall -flto ## aggressive
FFLAGS=-march=native -O3 -Wall -flto ## aggressive

## Integer8.  Note at least Arpack and Pardiso need to be compiled for
## Integer8 as well.  Spooles does not appear to be int8 according to Guido
## LONGLONG = -DLONGLONG
## CFLAGS  += $(LONGLONG)
## FFLAGS  += -fdefault-integer-8
## INTEXT = 64

## SPOOLES
## CFLAGS  += -I/usr/include/spooles -I/usr/include/spooles/MT -DSPOOLES
## LDFLAGS += -lspooles -lpthread 
CFLAGS  += `pkg-config --cflags spooles` -DSPOOLES
LDFLAGS += `pkg-config --libs spooles`

## ARPACK
## CFLAGS  += -DARPACK
## LDFLAGS += -larpack
CFLAGS  += `pkg-config --cflags arpack$(INTEXT)` -DARPACK
LDFLAGS += `pkg-config --libs arpack$(INTEXT)`

## LAPACK
## CFLAGS  += -I/usr/include/openblas
## LDFLAGS += -lreflapack -lopenblas
CFLAGS  += `pkg-config --cflags lapack$(INTEXT)`
LDFLAGS += `pkg-config --libs lapack$(INTEXT)`

## Multi Threaded and MPI
CFLAGS += -DUSE_MT
#CFLAGS += -DCALCULIX_MPI -fopenmp

## Other CCX Options
CFLAGS  += -DARCH="Linux" -DMATRIXSTORAGE 
LDFLAGS += 

## EXODUSII
## CFLAGS  += -DEXODUSII
## LDFLAGS += -lexoIIv2c -lnetcdf
CFLAGS  += `pkg-config --cflags exodusii` -DEXODUSII
LDFLAGS += `pkg-config --libs exodusii`

####

## ## TAUCS 
## ## CFLAGS  += -DTAUCS
## ## LDFLAGS += -ltaucs -lmetis
## 
## ## Flags for the gpu compiler
## ## NVCCCFLAGS = -arch=sm_35
## ## NVCCLDFLAGS = -lstdc++
## 
## ## BLAS
## ## CFLAGS  += -I/usr/include/openblas
## ## LDFLAGS += -lopenblas
## CFLAGS  += `pkg-config --cflags blas$(INTEXT)`
## LDFLAGS += `pkg-config --libs blas$(INTEXT)`

## PARDISO
## CFLAGS  += -DPARDISO
## LDFLAGS += -L/home/pete/local/lib64/ -lpardiso -lgfortran -lpthread -lm -fopenmp
CFLAGS  += `pkg-config --cflags pardiso` -DPARDISO
LDFLAGS += `pkg-config --libs pardiso`

## ## CUDACUSP
## ## This is unique because it a template library rather than binary library
## ## CFLAGS  += -I/home/pete/local/include -L/opt/cuda/lib64 -lcuda -lcudart -DCUDACUSP
## ## LDFLAGS += -L/opt/cuda/lib64 -lcuda -lcudart
## CFLAGS  += `pkg-config --cflags cusp` `pkg-config --libs cusp` -DCUDACUSP
## LDFLAGS += `pkg-config --libs cusp`
## 
## ## CHOLDMOD 
## CFLAGS  += -DSUITESPARSE
## LDFLAGS += `pkg-config --libs cholmod`
## LDFLAGS += `pkg-config --libs metis`
## LDFLAGS += `pkg-config --libs colamd`
## LDFLAGS += `pkg-config --libs ccolamd`
## LDFLAGS += `pkg-config --libs amd`
## LDFLAGS += `pkg-config --libs camd`
## ## ## SUITESPARSE (Also requires CUBLAS)
## ## LDFLAGS += `pkg-config --libs tbb`
## ## LDFLAGS += `pkg-config --libs rbio`
## ## LDFLAGS += `pkg-config --libs ldl`
## ## LDFLAGS += `pkg-config --libs cxsparse`
## ## LDFLAGS += `pkg-config --libs btf`
## ## LDFLAGS += `pkg-config --libs spqr`
## ## ## CUBLAS
## ## CFLAGS  += -I/opt/cuda/include -DGPU_BLAS
## ## LDFLAGS += -L/opt/cuda/lib64 -lcuda -lcudart
## ## CFLAGS  += `pkg-config --cflags cublas` -DGPU_BLAS
## ## LDFLAGS += `pkg-config --libs cublas`

CC=gcc
FC=gfortran
NVCC=nvcc -O3 $(LONGLONG) `pkg-config --cflags cusp` `pkg-config --libs cusp` $(NVCCCFLAGS) # -Xcompiler -fopenmp

## .cu file so not have a default implicit rule.  Define all implicit rules used.
.SUFFIXES: .o .c .cu

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<
.cu.o:
	$(NVCC) -DCUDACUSP -c $<

include Makefile.inc

SCCXMAIN = ccx_2.8p2.c

## Define all the object file rules to identify dependencies
OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

## Link to math and standard c
CFLAGS += -lm -lc

ccx_2.8p2: $(OCCXMAIN) ccx_2.8p2.a
	./date.pl; 
	$(CC) $(CFLAGS) -c ccx_2.8p2.c $(LDFLAGS); $(FC) -Wall $(FFLAGS) -o $@ $(OCCXMAIN) ccx_2.8p2.a $(LDFLAGS)

ccx_2.8p2.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

clean:
	rm *.a *.o
