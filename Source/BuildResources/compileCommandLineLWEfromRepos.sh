#!/bin/bash -l

# This will compile the command line version from the repos. Should be self contained,
# but you need to have CUDA and MKL available. By default it loads the modules
# as it is done on the MPCDF clusters. Depending on your GPU, you might need to change
# or add to the CUDA architecture in the cmake line

# NOTE THAT THIS WILL DELETE THE LightwaveExplorer FOLDER WHEN DONE

# load modules (you probably don't need this if you're not on the cluster

module purge
module load gcc/13
module load fftw-serial/3.3.10
module load cuda/12.6
module load cmake/3.30
module load ninja/1.11

echo "Cloning LWE repo... "
git clone https://github.com/NickKarpowicz/LightwaveExplorer >& /dev/null
cd LightwaveExplorer
mkdir build
cd build

echo "Starting to compile, this will take a couple of minutes... "
cmake -DCMAKE_CXX_COMPILER=gcc -DCLI=1 -DUSE_CUDA=1 -DCMAKE_CUDA_ARCHITECTURES=80 .. -GNinja
cmake --build . --config Release
echo "Cleaning up"
cp LightwaveExplorer ../../lwe
cd ..

cp Source/BuildResources/ClusterScripts/lweget.sh ../
cp Source/BuildResources/ClusterScripts/fileget.sh ../
cp Source/BuildResources/ClusterScripts/filesend.sh ../
cd ..
chmod +x lweget.sh
chmod +x fileget.sh
chmod +x filesend.sh
rm -rf LightwaveExplorer
