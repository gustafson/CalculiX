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
    rsync -avP spooles.a ${LIB}/libspooles.a || exit
    rsync -avP MT/src/spoolesMT.a ${LIB}/libspoolesMT.a || exit
fi
popd


## Build Arpack-ng (Next generation maintened arpack)
if [[ ! -f ${LIB}/libarpack.so ]]; then
    [[ -d arpack-ng ]] || git clone git@github.com:opencollab/arpack-ng.git
    [[ -d arpack-ng/build ]] || mkdir arpack-ng/build
    pushd arpack-ng/build
    cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ..
    make -j8 install
    popd
fi


## Build Exodus (library, not ccx extension)
if [[ ! -f ${LIB}/libexodus.so ]]; then
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
    rsync -avP lib/* ${LIB}/. || exit
    rsync -avP lib64/* ${LIB}/. || exit
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
    cp lib/libopenblas.a ${LIB}/. || exit
    popd
    popd
fi

## Determine active gcc.  Get the variable if it exists, or get it
## Use should override this with a specific version if necessary
if [[ ! $CC ]]; then
    CC=`which gcc`
    FC=`echo ${CC} |sed s/gcc/gfortran/`
fi
## export CC=gcc-9.4.0
## export FC=gfortran-9.4.0

## Build CalculiX.
### First set variables
GCCMV=`${CC} -v 2>&1 | grep "gcc version"|cut -f1 -d.|rev |cut -f1 -d" "|rev`
if [[ ${GCCMV} -gt 9 ]]; then
    echo "CalculiX currently builds only with gcc <= 9."
    echo "Identified version is ${GCCMV}."
    echo "Please set the CC and FC variables if another version is available."
    echo "An example command is: export CC=gcc-9.4.0; export FC=gfortran-9.4.0"
    exit
fi

### Then build
cd ccx/src
LD_LIBRARY_PATH=${LIB}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
export LIB
export PREFIX
make -f Makefile.superbuild -j8 || exit
rsync -avP ccx_2.19 ${BIN}/. || exit

## Add the library path
[[ `grep LD_LIBRARY_PATH ~/.bashrc | grep ${LIB}` ]] || echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${LIB}:${MKLROOT}/mkl/latest/lib/intel64:${MKLROOT}/compiler/latest/linux/compiler/lib/intel64/" >> ~/.bashrc
[[ `grep PATH ~/.bashrc | grep ${BIN}` ]] || echo "export PATH=\${PATH}:\${BIN}" >> ~/.bashrc

## Inform
echo Install appears to be successful
echo An attempt has been made to include the local prefix in the path
echo You need to activate the new path by sourcing the ~/.bashrc configuration file: source  ~/.bashrc
echo If libraries are still not found, manually add the library path using a command like: LD_LIBRARY_PATH+=:${PREFIX}/lib64\; export LD_LIBRARY_PATH
