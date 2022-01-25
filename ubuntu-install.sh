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

echo Installing compilers, libraries, and recommended tools
sudo apt-get install libexodusii-dev libspooles-dev libblas3 libblas-dev \
        liblapack3 liblapack-dev libarpack2 libarpack2-dev \
        libparpack2 gfortran gcc unzip make \
	paraview intel-mkl

echo Compiling source
cd ccx/src
make -j4 -f Makefile.ubuntu

echo "Install"
sudo cp ccx/src/ccx_2.19 /usr/local/bin/
sudo ln -s /usr/local/bin/ccx_2.19 /usr/local/bin/ccx
