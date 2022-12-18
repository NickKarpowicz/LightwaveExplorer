CC=g++
MKLROOT=/opt/intel/oneapi/mkl/latest
CFLAGS=-std=c++17 -fconcepts -use_fast_math -fopenmp -O3 -D LINUXCPUONLY
INCLUDES=`pkg-config --cflags gtk4` -I/opt/intel/oneapi/mkl/latest/include -I/opt/intel/oneapi/mkl/latest/include/fftw -I../../dlib
LDFLAGS=`pkg-config --libs gtk4`  ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -lgomp -lpthread -lm -ldl
SOURCES=LightwaveExplorerFrontendGTK.cpp ../LightwaveExplorerUtilities.cpp ../LightwaveExplorerCoreCPU.cpp ../DlibLibraryComponents/DlibLibraryComponents.cpp
OBJECTS=-o LightwaveExplorerGTK
EXECUTABLES=LightwaveExplorerGTK
CUDACOMPILER=/usr/local/cuda-12.0/bin/nvcc
CUDAFLAGS=-x cu
CUDAINCLUDES=-I${MKLROOT}/include -I${MKLROOT}/include/fftw -I../../dlib
CUDAARCH=-gencode=arch=compute_75,code=\"sm_75,compute_75\"
CUDAHOSTFLAGS=--machine 64 -std=c++17 -fconcepts -use_fast_math -fopenmp -O3
CUDALDFLAGS=-lcufft -lnvidia-ml
CUDAOUT= -o LightwaveExplorerGTK
CUDAFILES=../LightwaveExplorerCore.cu

default: cuda

cuda: 
	${CUDACOMPILER} ${CUDAARCH} ${CUDAFLAGS} ${CUDAINCLUDES} -D NOCUDAMAIN  -Xcompiler "${CUDAHOSTFLAGS} ${INCLUDES}" ${CUDAOUT} ${CUDAFILES} ${SOURCES} ${CUDALDFLAGS} -Xlinker "${LDFLAGS} ${CUDALDFLAGS}"

cpuonly:
	${CC} ${CFLAGS} ${INCLUDES} ${OBJECTS} ${SOURCES} ${LDFLAGS}

clean:
	rm LightwaveExplorerGTK


