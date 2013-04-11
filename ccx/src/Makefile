
CFLAGS = -Wall -O  -I ../../../SPOOLES.2.2 -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE
FFLAGS = -Wall -O -fopenmp

CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.4.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

DIR=../../../SPOOLES.2.2

LIBS = \
       $(DIR)/spooles.a \
	../../../ARPACK/libarpack_INTEL.a \
       -lm -lc

ccx_2.4: $(OCCXMAIN) ccx_2.4.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.4.c; $(FC) -Wall -O -o $@ $(OCCXMAIN) ccx_2.4.a  -lpthread  $(LIBS)

ccx_2.4.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
                                                                               
