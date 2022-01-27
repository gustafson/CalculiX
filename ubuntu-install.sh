#!/bin/bash

echo "Installing CalculiX and Extras"
echo "The following commands would as assumed to have already been issued:"
echo "    git clone https://github.com/gustafson/CalculiX"
echo "    cd CalculiX"

echo "Cloning CalculiX repo and checkout out ubuntu branch"
git checkout ubuntu

echo "Updating software repos (This could take some time)"
sudo apt-get update
echo "Optional: You should also update your system"
echo "    sudo apt-get upgrade"


function add_intel_apt(){
    ## Add intell for older ubuntu versions
    ## Modified from here but not executed directly due to security risk
    ## https://github.com/eddelbuettel/mkl4deb/blob/master/script.sh


    ## cf https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo

    pushd /tmp
    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    
    ## all products:
    #sudo wget https://apt.repos.intel.com/setup/intelproducts.list -O /etc/apt/sources.list.d/intelproducts.list
    ## just MKL
    sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
    ## other (TBB, DAAL, MPI, ...) listed on page

    sudo apt-get install intel-mkl-64bit-2020.0-088
    ## sudo apt-get install intel-mkl


    ## update alternatives
    sudo update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so     libblas.so-x86_64-linux-gnu      /opt/intel/mkl/lib/intel64/libmkl_rt.so 150
    sudo update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so.3   libblas.so.3-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 150
    sudo update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so   liblapack.so-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 150
    sudo update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so.3 liblapack.so.3-x86_64-linux-gnu  /opt/intel/mkl/lib/intel64/libmkl_rt.so 150

    sudo sh -c 'echo /opt/intel/lib/intel64     >  /etc/ld.so.conf.d/mkl.conf'
    sudo sh -c 'echo /opt/intel/mkl/lib/intel64 >> /etc/ld.so.conf.d/mkl.conf'
    sudo ldconfig

    sudo sh -c 'grep MKL_THREADING_LAYER /etc/environment || echo export MKL_THREADING_LAYER=GNU >> /etc/environment'
    rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    popd

    ## sudo ln -s /opt/intel/compilers_and_libraries/linux/mkl/bin/pkgconfig/mkl-* /usr/lib/pkgconfig/.
    ## sudo sed s.\${prefix}./opt/intel. /usr/lib/pkgconfig/mkl-* -i
}

echo Installing compilers, libraries, and recommended tools
sudo apt-get install libexodusii-dev libspooles-dev libblas3 libblas-dev \
        liblapack3 liblapack-dev libarpack2 libarpack2-dev \
        libparpack2 gfortran gcc unzip make paraview

sudo apt-get install intel-mkl || add_intel_apt

echo Compiling source
cd ccx/src
make -j4 -f Makefile.ubuntu

echo "Attempting install"
[[ -d /usr/local/bin ]] || sudo mkdir -p /usr/local/bin
sudo cp ccx_2.19 /usr/local/bin/
sudo sh -c "rm -f /usr/local/bin/ccx; ln -s /usr/local/bin/ccx_2.19 /usr/local/bin/ccx && echo Install completed successfully"
