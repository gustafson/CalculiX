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




################################################################################
################################################################################
##                                                                            ##
## Flags related to CPU, optimizations, etc                                   ##
## GPU based optimization later in the make file                              ##
##                                                                            ##
################################################################################
################################################################################

CC=gcc
FC=gfortran

## CFLAGS=-march=native -g -O0 ## debugging  -pg
## FFLAGS=-march=native -g -O0 ## debugging  -pg
## LDFLAGS += -g
## CFLAGS=-march=native -O2 -Wall ## conservative
## FFLAGS=-march=native -O2 -Wall ## conservative
CFLAGS=-march=native -O3 -Wall ## -flto ## aggressive -fprofile-generate -fprofile-use
FFLAGS=-march=native -O3 -Wall ## -flto ## aggressive -fprofile-generate -fprofile-use

## Integer8.  Note at least Arpack and Pardiso need to be compiled for
## Integer8 as well.  Spooles does not appear to be int8 according to Guido
## LONGLONG = -DLONGLONG
## CFLAGS  += $(LONGLONG)
## FFLAGS  += -fdefault-integer-8
## INTEXT = 64

## Other CCX Options
CFLAGS  += -DARCH="Linux" -DMATRIXSTORAGE
LDFLAGS +=

## Multi Threaded and MPI
CFLAGS += -DUSE_MT
#CFLAGS += -DCALCULIX_MPI -fopenmp

## This is now default for calculix and relates to CFD
CFLAGS += -DNETWORKOUT




################################################################################
################################################################################
##                                                                            ##
## Flags related to CPU based solvers                                         ##
##                                                                            ##
################################################################################
################################################################################

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

## TAUCS
## CFLAGS  += -DTAUCS
## LDFLAGS += -ltaucs -lmetis

## LAPACK
## CFLAGS  += -I/usr/include/openblas
## LDFLAGS += -lreflapack -lopenblas
CFLAGS  += `pkg-config --cflags lapack$(INTEXT)`
LDFLAGS += `pkg-config --libs lapack$(INTEXT)`

## BLAS
## CFLAGS  += -I/usr/include/openblas
## LDFLAGS += -lopenblas
CFLAGS  += `pkg-config --cflags blas$(INTEXT)`
LDFLAGS += `pkg-config --libs blas$(INTEXT)`

## PARDISO
## CFLAGS  += -DPARDISO
## LDFLAGS += -L/home/pete/local/lib64/ -lpardiso -lgfortran -lpthread -lm -fopenmp
CFLAGS  += `pkg-config --cflags pardiso` -DPARDISO
LDFLAGS += `pkg-config --libs pardiso`



################################################################################
################################################################################
##                                                                            ##
## Flags related to GPU based solvers                                         ##
##                                                                            ##
################################################################################
################################################################################

## Flags for the gpu compiler
## NVCCCFLAGS = -arch=sm_35
## NVCCLDFLAGS = -lstdc++
##
## NVCC=nvcc -O3 $(LONGLONG) `pkg-config --cflags cusp` `pkg-config --libs cusp` $(NVCCCFLAGS) # -Xcompiler -fopenmp

## CUDACUSP
## This is unique because it a template library rather than binary library
## CFLAGS  += -I/home/pete/local/include -L/opt/cuda/lib64 -lcuda -lcudart -DCUDACUSP
## LDFLAGS += -L/opt/cuda/lib64 -lcuda -lcudart
## CFLAGS  += `pkg-config --cflags cusp` `pkg-config --libs cusp` -DCUDACUSP
## LDFLAGS += `pkg-config --libs cusp`

## CHOLDMOD
## This is unique because it can be CPU or GPU based, depending on how
## SuiteSparse was compiled.  Here it is assumed that SuiteSparse also
## uses CUDA
## CFLAGS  += -DSUITESPARSE
## LDFLAGS += `pkg-config --libs cublas$(INTEXT)`
## LDFLAGS += `pkg-config --libs cholmod$(INTEXT)`
## LDFLAGS += `pkg-config --libs metis$(INTEXT)`
## LDFLAGS += `pkg-config --libs colamd$(INTEXT)`
## LDFLAGS += `pkg-config --libs ccolamd$(INTEXT)`
## LDFLAGS += `pkg-config --libs amd$(INTEXT)`
## LDFLAGS += `pkg-config --libs camd$(INTEXT)`
## LDFLAGS += `pkg-config --libs rbio$(INTEXT)`
## LDFLAGS += `pkg-config --libs ldl$(INTEXT)`
## LDFLAGS += `pkg-config --libs cxsparse$(INTEXT)`
## LDFLAGS += `pkg-config --libs btf$(INTEXT)`



################################################################################
################################################################################
##                                                                            ##
## Flags related to ExodusII output                                           ##
##                                                                            ##
################################################################################
################################################################################


## EXODUSII
## CFLAGS  += -DEXODUSII
## LDFLAGS += -lexoIIv2c -lnetcdf
## CFLAGS  += `pkg-config --cflags exodusii` -DEXODUSII
## LDFLAGS += `pkg-config --libs exodusii`




################################################################################
################################################################################
##                                                                            ##
## Recipes                                                                    ##
##                                                                            ##
################################################################################
################################################################################

## .cu file so not have a default implicit rule.  Define all implicit rules used.
.SUFFIXES: .o .c .cu

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<
.cu.o:
	$(NVCC) -DCUDACUSP -c $<

include Makefile.inc

SCCXMAIN = ccx_2.17.c

## Define all the object file rules to identify dependencies
OCCXCU = $(SCCXCU:.cu=.o)
OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

## Link to math and standard c
CFLAGS += -lm -lc

ccx_2.17: $(OCCXMAIN) ccx_2.17.a
	./date.pl;
	$(CC) $(CFLAGS) -c ccx_2.17.c $(LDFLAGS); $(FC) -Wall $(FFLAGS) -o $@ $(OCCXMAIN) ccx_2.17.a $(LDFLAGS)

ccx_2.17.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

clean:
	rm *.a *.o
