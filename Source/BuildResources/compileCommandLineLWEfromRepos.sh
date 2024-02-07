#!/bin/bash -l

# This will compile the command line version from the repos. Should be self contained,
# but you need to have CUDA and MKL available. By default it loads the modules
# as it is done on the MPCDF clusters. Depending on your GPU, you might need to change
# or add to the gencodes in the nvcc line
# NOTE THAT THIS WILL DELETE THE dlib AND LightwaveExplorer DIRECTORIES WHEN DONE

# load modules (you probably don't need this if you're not on the cluster
# but you will need to set MKL_HOME to the location of MKL)
module purge
module load gcc/12
module load mkl/2024.0
module load cuda/12.2
module load cmake/3.28

echo "Cloning LWE repo... "
git clone https://github.com/NickKarpowicz/LightwaveExplorer >& /dev/null
cd LightwaveExplorer
git checkout compressedOutputFile
mkdir build
cd build

echo "Starting to compile, this will take a couple of minutes... "
cmake -DCMAKE_CXX_COMPILER=gcc -DCLICUDA=True -DMKL_ROOT=$MKLROOT -DCMAKE_CUDA_ARCHITECTURES=80 .. --fresh
make
echo "Cleaning up"
cp LightwaveExplorer ../../lwe
cd ..

cp CrystalDatabase.txt ../CrystalDatabase.txt
cd ..
rm -rf LightwaveExplorer
