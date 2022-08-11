#!/bin/bash
## This script written by Peter A. Gustafson (peter.gustafson@wmich.edu)
## Copyright 2022, all rights reserved

[[ `env |grep PREFIX` ]] || PREFIX=~/local
[[ `env |grep MKLROOT` ]] || MKLROOT=/opt/intel/oneapi ## MKLROOT=/appl/intelv2019.6/mkl

echo Calculix will be installed into ${PREFIX} using MKL located at ${MKLROOT}

## sudo zypper addrepo https://yum.repos.intel.com/oneapi oneAPI
## zypper install intel-basekit intel-hpckit 

LIB=${PREFIX}/lib64
BIN=${PREFIX}/bin


[[ -d ${LIB} ]] || mkdir -p ${LIB}
[[ -d ${BIN} ]] || mkdir -p ${BIN}

## Download CalculiX with custom mods
[[ -d CalculiX ]] || git clone gustafson.github.io:gustafson/CalculiX.git
pushd CalculiX
git checkout opensuse


## Build Spooles
[[ -d spooles ]] || mkdir spooles
pushd spooles

if [[ ! -f ${LIB}/libspoolesMT.a ]]; then
    ## Download
    [[ -f spooles.2.2.tgz ]] || wget http://www.netlib.org/linalg/spooles/spooles.2.2.tgz
    tar xavf spooles.2.2.tgz

    ## Large integer correction
    [[ -f ccx_2.19.SPOOLEScorrection.tar.bz2 ]] || wget http://www.dhondt.de/ccx_2.19.SPOOLEScorrection.tar.bz2
    tar xavf ccx_2.19.SPOOLEScorrection.tar.bz2
    mv CalculiX/ccx_2.19/SPOOLES.2.2/I2Ohash/src/util.c I2Ohash/src/util.c
    rm -fr CalculiX

    ## Makefile changes
    sed s/drawTree/draw/ -i ./Tree/src/makeGlobalLib 
    grep MT makefile -B2|sed s/"#cd MT"/"\tcd MT"/ -i makefile 
    sed -i -e s/"# CC = cc"/"  CC = gcc"/ -e s/"# OPTLEVEL = -O3"/"OPTLEVEL = -O3"/ Make.inc

    ## Build
    make -j8 lib
    rsync -avP spooles.a ${LIB}/libspooles.a
    rsync -avP MT/src/spoolesMT.a ${LIB}/libspoolesMT.a
fi
popd


## Build Arpack
if [[ ! -f ${LIB}/libarpack.a ]]; then
   [[ -d arpack ]] || mkdir arpack
   pushd arpack
   [[ -f arpack96.tar.gz ]] || wget https://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz
   [[ -f patch.tar.gz ]] || wget https://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz
   tar xavf arpack96.tar.gz
   tar xavf patch.tar.gz
   cd ARPACK
   sed -i -e s./bin/make./usr/bin/make. \
       -e s/"#DIRS         = \$(UTILdir) \$(SRCdir)"/"DIRS         = \$(UTILdir) \$(SRCdir)"/ \
       -e s/FC.*/FC=gfortran/ \
       -e "s/home =.*/home=./" \
       -e "s/libarpack_.*/libarpack.a/" \
       -e "s/FFLAGS.*/FFLAGS=-O3/" -i ARmake.inc
   sed -i -e s/"      EXTERNAL           ETIME"/"*      EXTERNAL           ETIME"/ UTIL/second.f

   patch << 'EOF'
--- Makefile~	1996-09-24 10:55:31.000000000 -0400
+++ Makefile	2021-12-14 17:27:59.964489948 -0500
@@ -62,7 +62,8 @@
 		$(MAKE) $(PRECISIONS); \
 		$(CD) ..; \
 	done );
-	$(RANLIB) $(ARPACKLIB)
+	ar rv libarpack.a */*.o
+	ranlib libarpack.a
 
 cleantest:
 
EOF
   make lib
   rsync -avP libarpack.a ${LIB}/.
   popd
fi


## Build Exodus (library, not ccx extension)
if [[ ! -f ${LIB}/libexodus.a ]]; then
    git clone https://github.com/gsjaardema/seacas.git
    cd seacas/ && export ACCESS=`pwd`
    ## INSTALL_PATH=${PREFIX}
    NEEDS_ZLIB=YES USE_ZLIB_NG=YES PARMETIS=YES JOBS=10 ./install-tpl.sh
    cd $ACCESS
    mkdir build
    cd build
    ../cmake-exodus
    make -j8
    rsync -avP packages/seacas/libraries/exodus/libexodus* ${LIB}/.
    ## make install
    cd ..
    rsync -avP lib/* ${LIB}/.
    rsync -avP lib64/* ${LIB}/.
fi

## Build OpenBLAS
if [[ ! -f ${LIB}/libopenblas.a ]]; then
    [[ -d OpenBLAS ]] || git clone https://github.com/xianyi/OpenBLAS.git
    pushd OpenBLAS
    git checkout v0.3.19
    [[ -d build ]] || mkdir build
    pushd build
    cmake -D BUILD_RELAPACK=OFF -D BUILD_STATIC_LIBS=OFF -D CMAKE_BUILD_TYPE=Release ..
    make -j8
    cp lib/libopenblas.a ${LIB}/.
    popd
    popd
fi

## Build CalculiX
cd ccx/src
make -f Makefile.opensuse -j8
rsync -avP ccx_2.19 ~/local/bin/.

## Add the library path
[[ `grep LD_LIBRARY_PATH ~/.bashrc | grep local/lib64` ]] || echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:~/local/lib64:${MKLROOT}/mkl/latest/lib/intel64:${MKLROOT}/compiler/latest/linux/compiler/lib/intel64/" >> ~/.bashrc
[[ `grep PATH ~/.bashrc | grep local/bin` ]] || echo "export PATH=\${PATH}:~/local/bin" >> ~/.bashrc


