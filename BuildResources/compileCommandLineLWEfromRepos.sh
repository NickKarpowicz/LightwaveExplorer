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
module load mkl/2022.2
module load cuda/11.6
git clone https://github.com/NickKarpowicz/LightwaveExplorer
git clone https://github.com/davisking/dlib

cd LightwaveExplorer
echo "Starting to compile, this will take a couple of minutes... "
gcc -I$MKL_HOME/include -I$MKL_HOME/include/fftw -O3 -std=c++20 -c LightwaveExplorerDevices/LightwaveExplorerUtilities.cpp
nvcc -gencode=arch=compute_70,code=\"sm_70,compute_70\" -allow-unsupported-compiler -gencode=arch=compute_75,code=\"sm_75,compute_75\" -gencode=arch=compute_80,code=\"sm_80,compute_80\" -x cu -I$MKL_HOME/include -I$MKL_HOME/include/fftw -I../dlib --machine 64 -Xcompiler -std=c++17 -use_fast_math -O3 -L$MKL_HOME/lib/intel64 -lcufft -lnvidia-ml -lmkl_sequential -lmkl_core -lmkl_intel_lp64 -o lwe LightwaveExplorerCore.cu LightwaveExplorerUtilities.o LightwaveExplorerDevices/DlibLibraryComponents.cpp

echo "Cleaning up "
cp lwe ../lwe
cp CrystalDatabase.txt ../CrystalDatabase.txt
cd ..
rm -rf dlib
rm -rf LightwaveExplorer