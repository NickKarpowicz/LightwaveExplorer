#!/bin/bash -l

#find-and-replace to use fftw, replace std::format .etc with fmt::format since clang doesn't have it yet
sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LightwaveExplorerUtilities.h
sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LWEActiveDeviceCPU.h 
sed -i'.bak' 's!<format>!<fmt/format.h>!g ; s/std::format/fmt::format/g ; s/std::vformat/fmt::vformat/g ; s/std::make_format_args/fmt::make_format_args/g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
sed -i'.bak' 's!<format>!<fmt/format.h>!g ; s/std::format/fmt::format/g ; s/std::vformat/fmt::vformat/g ; s/std::make_format_args/fmt::make_format_args/g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp

mv BuildResources/CMakeListsLinuxCPU.txt CMakeLists.txt
#build executable
rm -rf build
mkdir build
cd build
cmake ..
make
cd ..
rm CMakeLists.txt

#restore the original source and clean up
cp AppImageCPU/COPYING COPYING
tar cf GPLsource.tar COPYING makefile *.cpp *.cu *.h LightwaveExplorerGTK/* DlibLibraryComponents/* MacResources/*
rm COPYING
rm LightwaveExplorerUtilities.h
rm LWEActiveDeviceCPU.h
rm LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
rm LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
mv LightwaveExplorerUtilities.h.bak LightwaveExplorerUtilities.h
mv LWEActiveDeviceCPU.h.bak LWEActiveDeviceCPU.h
mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
