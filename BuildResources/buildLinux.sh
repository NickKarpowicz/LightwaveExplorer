#!/bin/bash -l

#find-and-replace to use fftw, replace std::format .etc with fmt::format since clang doesn't have it yet
sed -i'.bak' 's!<format>!<fmt/format.h>!g ; s/std::format/fmt::format/g ; s/std::vformat/fmt::vformat/g ; s/std::make_format_args/fmt::make_format_args/g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
sed -i'.bak' 's!<format>!<fmt/format.h>!g ; s/std::format/fmt::format/g ; s/std::vformat/fmt::vformat/g ; s/std::make_format_args/fmt::make_format_args/g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp

mv BuildResources/CMakeListsLinux.txt CMakeLists.txt
#build executable
rm -rf build
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.0/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=75 -DMKLROOT=/home/nick/intel/oneapi/mkl/2023.0.0 ..
make
cd ..
rm CMakeLists.txt

#restore the original source and clean up
mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
