CC=g++
CFLAGS=-g -O3 -fopenmp -std=c++20
MKLROOT=/opt/intel/oneapi/mkl/latest
INCLUDES=`pkg-config --cflags gtk4` -I/opt/intel/oneapi/mkl/latest/include -I/opt/intel/oneapi/mkl/latest/include/fftw -I../../dlib
LDFLAGS=`pkg-config --libs gtk4`  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
SOURCES=LightwaveExplorerFrontendGTK.cpp ../LightwaveExplorerUtilities.cpp ../LightwaveExplorerCoreCPU.cpp ../DlibLibraryComponents/DlibLibraryComponents.cpp
EXECUTABLES=LightwaveExplorerGTK

default: lwe

lwe:
	${CC} ${CFLAGS} ${INCLUDES} -o ${EXECUTABLES} ${SOURCES} ${LDFLAGS}

clean:
	rm LightwaveExplorerGTK

