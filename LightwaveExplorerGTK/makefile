CC=g++
MKLROOT=/opt/intel/oneapi/mkl/latest
ONEAPIROOT=/opt/intel/oneapi
CFLAGS=-std=c++20 -use_fast_math -fopenmp --machine 64 -O3 -D LINUXCPUONLY
INCLUDES=`pkg-config --cflags gtk4` -I/opt/intel/oneapi/mkl/latest/include -I/opt/intel/oneapi/mkl/latest/include/fftw -I../../dlib
LDFLAGS=`pkg-config --libs gtk4` `pkg-config --libs gtk4` -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
SOURCES=LightwaveExplorerFrontendGTK.cpp ../LightwaveExplorerUtilities.cpp ../LightwaveExplorerCoreCPU.cpp ../DlibLibraryComponents/DlibLibraryComponents.cpp
OBJECTS=-o LightwaveExplorerGTK
CUDACOMPILER=/usr/local/cuda-12.0/bin/nvcc
CUDAFLAGS=-x cu
CUDAINCLUDES=-I${MKLROOT}/include -I${MKLROOT}/include/fftw -I../../dlib
CUDAARCH=-gencode=arch=compute_75,code=\"sm_75,compute_75\"
CUDAHOSTFLAGS=--machine 64 -std=c++17 -fconcepts -use_fast_math -fopenmp -O3 -L/home/nick/repos/LightwaveExplorer/LightwaveExplorerGTK
CUDALDFLAGS=-L. -lcufft -lnvidia-ml `pkg-config --libs gtk4` ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -lLWESYCL -lgomp -lpthread -lm -ldl
CUDAOUT= -o LightwaveExplorerGTK
CUDAFILES=../LightwaveExplorerCore.cu
TBBROOT=/opt/intel/oneapi/tbb/latest
DPCPP=/opt/intel/oneapi/compiler/latest/linux/bin/icpx
DPCPPFLAGS=-O2 -fsycl -fsycl-unnamed-lambda -fsycl-early-optimizations -DMKL_ILP64 -D "NDEBUG" -D "LINUXCPUONLY" -std=c++20 -fPIC -fsycl-link=image
DPCPPFILES=LightwaveExplorerDPCPPlib.cpp ../LightwaveExplorerUtilities.cpp ../DlibLibraryComponents/DlibLibraryComponents.cpp
DPCPPOUTPUT= -L.
DPCPPINCLUDES=-I../ -I${ONEAPIROOT}/compiler/latest/linux/include -I${ONEAPIROOT}/compiler/latest/linux/include/sycl -I${MKLROOT}/include -I${MKLROOT}/include/fftw -I${ONEAPIROOT}/dpl/latest/linux/include -I../../dlib -I../
DPCPPLD=   -L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 -lmkl_tbb_thread -lmkl_core -lsycl -lOpenCL -lpthread -lm -ldl

default: sycl sycl2 cuda

cuda: 
	${CUDACOMPILER} ${CUDAARCH} ${CUDAFLAGS} ${CUDAINCLUDES} -D NOCUDAMAIN  -Xcompiler "${CUDAHOSTFLAGS} ${INCLUDES}" ${CUDAOUT} ${CUDAFILES} ${SOURCES} -Xlinker "${CUDALDFLAGS}"

cpuonly:
	${CC} ${CFLAGS} ${INCLUDES} ${OBJECTS} ${SOURCES} ${LDFLAGS}

sycl:
	${DPCPP} ${DPCPPFLAGS} ${DPCPPINCLUDES} ${DPCPPFILES} -c

sycl2:
	${DPCPP} -fsycl -fPIC -shared LightwaveExplorerDPCPPlib.o LightwaveExplorerUtilities.o DlibLibraryComponents.o -o libLWESYCL.so -L. ${DPCPPLD}
clean:
	rm -f LightwaveExplorerGTK
	rm -r -f TestFile*
	rm -f libLWESYCL.*
	rm -rf *.o
