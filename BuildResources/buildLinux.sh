#!/bin/bash -l

#build executable
rm -rf build
mkdir build
cd build
cmake -DONEAPI_ROOT=${ONEAPI_ROOT} -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.0/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=75 ..
cmake -DCMAKE_INSTALL_PREFIX=AppDir/usr -DONEAPI_ROOT=${ONEAPI_ROOT} -DSYCL_INCLUDE_DIR=${ONEAPI_ROOT}/compiler/latest/linux/include/sycl -DSYCL_LIBRARY_DIR=${ONEAPI_ROOT}/compiler/latest/linux/lib -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.0/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=75 -DMKLROOT=${MKLROOT} ..
make
cd ..
