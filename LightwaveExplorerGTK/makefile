CC=g++
DPCPP=/opt/intel/oneapi/compiler/latest/linux/bin/icpx
NVCC=/usr/local/cuda-12.0/bin/nvcc
CUDATARGETS = /usr/local/cuda-12.0/targets/x86_64-linux

CFLAGS=-std=c++20 -use_fast_math -fopenmp --machine 64 -w -O3 -D CPUONLY
INCLUDES=`pkg-config --cflags gtk4` -I${MKLROOT}/include -I${MKLROOT}/include/fftw -I../../dlib
LDFLAGS=`pkg-config --libs gtk4` `pkg-config --libs gtk4` -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
SOURCES=LightwaveExplorerFrontendGTK.cpp ../LightwaveExplorerUtilities.cpp ../LightwaveExplorerCoreCPU.cpp ../DlibLibraryComponents/DlibLibraryComponents.cpp
OBJECTS=-o LightwaveExplorerGTK

CUDAFLAGS= -diag-suppress 1650 -diag-suppress 1217 -x cu -D NOCUDAMAIN
CUDAINCLUDES=-I${MKLROOT}/include -I${MKLROOT}/include/fftw -I../../dlib
CUDAARCH=-gencode=arch=compute_75,code=\"sm_75,compute_75\"
CUDAHOSTFLAGS=--machine 64 -std=c++17 -fconcepts -use_fast_math -fopenmp -O3 -Wno-unused-result  -L/home/nick/repos/LightwaveExplorer/LightwaveExplorerGTK
CUDALDFLAGS=-L. -lcufft -lnvidia-ml ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -lgomp -lpthread -lm -ldl
CUDAOUT = -c
CUDAFILES=../LightwaveExplorerCore.cu

DPCPPFLAGS=-Ofast -w -fsycl -fsycl-unnamed-lambda -fsycl-early-optimizations -DMKL_ILP64 -std=c++20 -fsycl-device-code-split=per_kernel
DPCPPFILES=LightwaveExplorerDPCPPlib.cpp LightwaveExplorerUtilities.o DlibLibraryComponents.o LightwaveExplorerCoreCPU.o LightwaveExplorerFrontendGTK.o LightwaveExplorerCore.o 
DPCPPOUTPUT= -o LightwaveExplorerGTK
DPCPPINCLUDES=-I../ -I${ONEAPIROOT}/compiler/latest/linux/include -I${ONEAPIROOT}/compiler/latest/linux/include/sycl -I${MKLROOT}/include -I${MKLROOT}/include/fftw -I${ONEAPIROOT}/dpl/latest/linux/include -I../../dlib -I../
DPCPPLD=-L${CUDATARGETS}/lib/ -L${CUDATARGETS}/lib/stubs -lcufft -lnvidia-ml -lcudart -lgomp `pkg-config --libs gtk4` `pkg-config --libs gtk4` ${MKLROOT}/lib/intel64/libmkl_sycl.a -Wl,-export-dynamic -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_tbb_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L${TBBROOT}/lib/intel64/gcc4.8 -ltbb -lsycl -lOpenCL -lpthread -lm -ldl

default: clean cuda sycl

cuda: 
	${NVCC} ${CUDAARCH} ${CUDAFLAGS} ${CUDAINCLUDES} -Xcompiler "${CUDAHOSTFLAGS} ${INCLUDES}" ${CUDAOUT} ${CUDAFILES} ${SOURCES} -Xlinker "${CUDALDFLAGS}"

cpuonly:
	${CC} ${CFLAGS} ${INCLUDES} ${OBJECTS} ${SOURCES} ${LDFLAGS}

sycl:
	${DPCPP} ${DPCPPFLAGS} ${DPCPPINCLUDES} ${INCLUDES} ${DPCPPFILES} ${DPCPPOUTPUT} ${DPCPPLD}

clean:
	rm -f LightwaveExplorerGTK
	rm -r -f TestFile*
	rm -rf *.o
