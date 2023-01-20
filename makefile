CC=g++
HBPATH=`brew --prefix`
APPLECC=${HBPATH}/opt/llvm/bin/clang++
DPCPP=icpx
NVCC=/usr/local/cuda-12.0/bin/nvcc
CUDATARGETS=/usr/local/cuda-12.0/targets/x86_64-linux

CFLAGS=-std=c++20 -fopenmp -mtune=native -w -Ofast -D CPUONLY
INCLUDES=`pkg-config --cflags gtk4` -I${MKLROOT}/include -I${MKLROOT}/include/fftw -I../dlib

LDFLAGS=`pkg-config --libs gtk4` `pkg-config --libs gtk4` -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
SOURCES= LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp LightwaveExplorerUtilities.cpp LightwaveExplorerCoreCPU.cpp DlibLibraryComponents/DlibLibraryComponents.cpp
OBJECTS=-o LightwaveExplorer

#mac and CPU-only builds compile with fftw instead of mkl
LDFLAGSF=`pkg-config --libs gtk4` `pkg-config --libs gtk4` /usr/lib/x86_64-linux-gnu/libfftw3.a /usr/local/lib/libfmt.a -lgomp -lpthread -lm -ldl
INCLUDESF=`pkg-config --cflags gtk4` -I../dlib

#CUDA settings
CUDAFLAGS= -diag-suppress 1650 -diag-suppress 1217 -x cu -D NOCUDAMAIN
CUDAINCLUDES=-I${MKLROOT}/include -I${MKLROOT}/include/fftw -I../dlib
CUDAARCH=-gencode=arch=compute_75,code=\"sm_75,compute_75\"
CUDAHOSTFLAGS=--machine 64 -std=c++17 -fconcepts -fopenmp -Ofast -w
CUDALDFLAGS=-L. -lcufft -lnvidia-ml ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -lgomp -lpthread -lm -ldl
CUDAOUT = -c
CUDASOURCE=LightwaveExplorerCore.cu

#DPC++/SYCL settings
DPCPPFLAGS=-Ofast -w -fsycl -fsycl-unnamed-lambda -fsycl-early-optimizations -DMKL_ILP64 -std=c++20
DPCPPFILES=./LightwaveExplorerGTK/LightwaveExplorerDPCPPlib.cpp LightwaveExplorerUtilities.o DlibLibraryComponents.o LightwaveExplorerCoreCPU.o LightwaveExplorerFrontendGTK.o LightwaveExplorerCore.o
DPCPPFILESNOCUDA=./LightwaveExplorerGTK/LightwaveExplorerDPCPPlib.cpp LightwaveExplorerUtilities.cpp DlibLibraryComponents/DlibLibraryComponents.cpp LightwaveExplorerCoreCPU.cpp LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
DPCPPOUTPUT= LightwaveExplorer
DPCPPINCLUDES=-I${ONEAPIROOT}/compiler/latest/linux/include -I${ONEAPIROOT}/compiler/latest/linux/include/sycl -I${MKLROOT}/include -I${MKLROOT}/include/fftw -I${ONEAPIROOT}/dpl/latest/linux/include -I../dlib -I./
DPCPPLD=-L${CUDATARGETS}/lib/ -L${CUDATARGETS}/lib/stubs -lcufft -lnvidia-ml -lcudart /usr/local/lib/libfmt.a -lgomp `pkg-config --libs gtk4` `pkg-config --libs gtk4` ${MKLROOT}/lib/intel64/libmkl_sycl.a -Wl,-export-dynamic -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_tbb_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L${TBBROOT}/lib/intel64/gcc4.8 -ltbb -lsycl -lOpenCL -lpthread -lm -ldl

#Apple compilation

ARMHOMEBREW=/opt/homebrew
ARMHOMEBREWXC=/Users/nick/arm-target
APPLEFLAGS=-std=c++20 -O3 -fopenmp -D CPUONLY -w
APPLEINCLUDES=-I../dlib -I${HBPATH}/include/fmt -I${HBPATH}/opt/libomp/include -I${HBPATH}/include -I${HBPATH}/include/gtk-4.0 -I${HBPATH}/include/pango-1.0 -I${HBPATH}/include/glib-2.0 -I${HBPATH}/include/cairo -I${HBPATH}/lib/glib-2.0/include -I${HBPATH}/include/fontconfig -I${HBPATH}/include/freetype2 -I${HBPATH}/include/gdk-pixbuf-2.0 -I${HBPATH}/include/harfbuzz -I${HBPATH}/include/graphene-1.0 -I${HBPATH}/lib/graphene-1.0/include
APPLELDFLAGS=-L${HBPATH}/lib ${HBPATH}/opt/libomp/lib/libomp.a ${HBPATH}/lib/libfmt.a -lpthread -lm -ldl -lgtk-4 -lgio-2.0 -lpangoft2-1.0 -lgdk_pixbuf-2.0 -lcairo -lpango-1.0 -lfreetype -lfontconfig -lgobject-2.0 -lglib-2.0 -lgthread-2.0 ${HBPATH}/lib/libfftw3.a

LINUXCLANGFLAGS=-std=c++20 -fexperimental-library -stdlib=libc++ -O3 -fopenmp -D CPUONLY -w
LINUXCLANGINCLUDES=-I../dlib -I${HBPATH}/opt/llvm/include -I${HBPATH}/opt/libomp/include -I${HBPATH}/include -I${HBPATH}/include/gtk-4.0 -I${HBPATH}/include/pango-1.0 -I${HBPATH}/include/glib-2.0 -I${HBPATH}/include/cairo -I${HBPATH}/lib/glib-2.0/include -I${HBPATH}/include/fontconfig -I${HBPATH}/include/freetype2 -I${HBPATH}/include/gdk-pixbuf-2.0 -I${HBPATH}/include/harfbuzz -I${HBPATH}/include/graphene-1.0 -I${HBPATH}/lib/graphene-1.0/include
LINUXCLANGLDFLAGS=-L${HBPATH}/lib -L${HBPATH}/opt/llvm/lib/c++ -L${HBPATH}/opt/llvm/lib ${HBPATH}/opt/libomp/lib/libomp.a -lc++ -lc++abi -lpthread -lm -ldl -lgtk-4 -lgio-2.0 -lpangoft2-1.0 -lgdk_pixbuf-2.0 -lcairo -lpango-1.0 -lfreetype -lfontconfig -lgobject-2.0 -lglib-2.0 -lgthread-2.0 ${HBPATH}/lib/libfftw3.a

#cross-compile settings for compiling for Arm64 on Intel
APPLEINCLUDESARMXC=-I../dlib -I${ARMHOMEBREWXC}/opt/llvm/include -I${ARMHOMEBREWXC}/opt/libomp/include -I${ARMHOMEBREWXC}/include -I${ARMHOMEBREWXC}/include/gtk-4.0 -I${ARMHOMEBREWXC}/include/pango-1.0 -I${ARMHOMEBREWXC}/include/glib-2.0 -I${ARMHOMEBREWXC}/include/cairo -I${ARMHOMEBREWXC}/lib/glib-2.0/include -I${ARMHOMEBREWXC}/include/fontconfig -I${ARMHOMEBREWXC}/include/freetype2 -I${ARMHOMEBREWXC}/include/gdk-pixbuf-2.0 -I${ARMHOMEBREWXC}/include/harfbuzz -I${ARMHOMEBREWXC}/include/graphene-1.0 -I${ARMHOMEBREWXC}/lib/graphene-1.0/include
APPLELDFLAGSARMXC=-L${ARMHOMEBREWXC}/opt/llvm/lib ${ARMHOMEBREWXC}/opt/libomp/lib/libomp.a -L${ARMHOMEBREWXC}/lib -lgtk-4 -lgio-2.0 -lpangoft2-1.0 -lgdk_pixbuf-2.0 -lcairo -lpango-1.0 -lfreetype -lfontconfig -lgobject-2.0 -lglib-2.0 -lgthread-2.0 ${ARMHOMEBREWXC}/lib/libfftw3.a

default: clean cuda sycl

cuda: 
	sed -i'.bak' 's#<format>#<fmt/format.h>#g ; s#std::format#fmt::format#g ; s#std::vformat#fmt::vformat#g ; s#std::make_format_args#fmt::make_format_args#g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
	sed -i'.bak' 's#<format>#<fmt/format.h>#g ; s#std::format#fmt::format#g ; s#std::vformat#fmt::vformat#g ; s#std::make_format_args#fmt::make_format_args#g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
	${NVCC} ${CUDAARCH} ${CUDAFLAGS} ${CUDAINCLUDES} -Xcompiler "${CUDAHOSTFLAGS} ${INCLUDES}" ${CUDAOUT} ${CUDASOURCE} ${SOURCES} -Xlinker "${CUDALDFLAGS}"
	rm LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
	rm LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
	mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
	mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h

#Note that the Mac version and cpuonly versions make use of FFTW3 instead of MKL for the Fourier transforms
#as such they are released under GPL3 - the temporary codebase from which they are compiled is saved in the
#GPLsource.tar file along with a copy of the GPLv3.
cpuonly:
	sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LightwaveExplorerUtilities.h
	sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LWEActiveDeviceCPU.h 
	sed -i'.bak' 's#<format>#<fmt/format.h>#g ; s#std::format#fmt::format#g ; s#std::vformat#fmt::vformat#g ; s#std::make_format_args#fmt::make_format_args#g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
	sed -i'.bak' 's#<format>#<fmt/format.h>#g ; s#std::format#fmt::format#g ; s#std::vformat#fmt::vformat#g ; s#std::make_format_args#fmt::make_format_args#g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
	cp AppImageCPU/COPYING COPYING
	${CC} ${CFLAGS} ${INCLUDESF} ${OBJECTS} ${SOURCES} ${LDFLAGSF}
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

macARMonIntel:
	sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LightwaveExplorerUtilities.h
	sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LWEActiveDeviceCPU.h 
	cp AppImageCPU/COPYING COPYING
	${APPLECC} -target arm64-apple-macos12 ${APPLEFLAGS} ${APPLEINCLUDESARMXC} -o LightwaveExplorer_Arm64 ${SOURCES} ${APPLELDFLAGSARMXC}
	${APPLECC} ${APPLEFLAGS} ${APPLEINCLUDES} -o LightwaveExplorer_x86 ${SOURCES} ${APPLELDFLAGS}
	lipo -create -output LightwaveExplorer LightwaveExplorer_x86 LightwaveExplorer_Arm64
	tar cf GPLsource.tar COPYING makefile *.cpp *.cu *.h LightwaveExplorerGTK/* DlibLibraryComponents/* MacResources/*
	rm COPYING
	rm LightwaveExplorerUtilities.h
	rm LWEActiveDeviceCPU.h
	mv LightwaveExplorerUtilities.h.bak LightwaveExplorerUtilities.h
	mv LWEActiveDeviceCPU.h.bak LWEActiveDeviceCPU.h

mac:
	sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LightwaveExplorerUtilities.h
	sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LWEActiveDeviceCPU.h 
	sed -i'.bak' 's/<format>/<format.h>/g ; s/std::format/fmt::format/g ; s/std::vformat/fmt::vformat/g ; s/std::make_format_args/fmt::make_format_args/g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
	sed -i'.bak' 's/<format>/<format.h>/g ; s/std::format/fmt::format/g ; s/std::vformat/fmt::vformat/g ; s/std::make_format_args/fmt::make_format_args/g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
	cp AppImageCPU/COPYING COPYING
	${APPLECC} ${APPLEFLAGS} ${APPLEINCLUDES} ${OBJECTS} ${SOURCES} ${APPLELDFLAGS}
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

cleanmacfail:
	rm -f COPYING
	rm -f LightwaveExplorerUtilities.h
	rm -f LWEActiveDeviceCPU.h
	rm -f LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
	rm -f LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
	mv LightwaveExplorerUtilities.h.bak LightwaveExplorerUtilities.h
	mv LWEActiveDeviceCPU.h.bak LWEActiveDeviceCPU.h
	mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
	mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h

linuxclangfftw:
	sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LightwaveExplorerUtilities.h
	sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LWEActiveDeviceCPU.h 
	cp AppImageCPU/COPYING COPYING
	${APPLECC} ${LINUXCLANGFLAGS} ${LINUXCLANGINCLUDES} ${OBJECTS} ${SOURCES} ${LINUXCLANGLDFLAGS}
	tar cf GPLsource.tar COPYING makefile *.cpp *.cu *.h LightwaveExplorerGTK/* DlibLibraryComponents/* MacResources/*
	rm COPYING
	rm LightwaveExplorerGTK/LightwaveExplorerUtilities.h
	rm LightwaveExplorerGTK/LWEActiveDeviceCPU.h
	mv LightwaveExplorerGTK/LightwaveExplorerUtilities.h.bak LightwaveExplorerGTK/LightwaveExplorerUtilities.h
	mv LightwaveExplorerGTK/LWEActiveDeviceCPU.h.bak LightwaveExplorerGTK/LWEActiveDeviceCPU.h

sycl:
	${DPCPP} ${DPCPPFLAGS} ${DPCPPINCLUDES} ${INCLUDES} ${DPCPPFILES} -o ${DPCPPOUTPUT} ${DPCPPLD}

nocuda:
	${DPCPP} ${DPCPPFLAGS} -D NOCUDA ${DPCPPINCLUDES} ${INCLUDES} ${DPCPPFILESNOCUDA} -o ${DPCPPOUTPUT} ${DPCPPLD}

install:
	cp ${DPCPPOUTPUT} /usr/bin/${DPCPPOUTPUT}
	mkdir -p /usr/share/LightwaveExplorer
	chmod +x ./LWElauncher/LightwaveExplorerLauncher.sh
	cp LWElauncher/LightwaveExplorerLauncher.sh /usr/bin/LightwaveExplorerLauncher.sh
	cp CrystalDatabase.txt /usr/share/LightwaveExplorer/CrystalDatabase.txt
	cp DefaultValues.ini /usr/share/LightwaveExplorer/DefaultValues.ini

clean:
	rm -f LightwaveExplorer
	rm -rf TestFile*
	rm -rf *.o
	rm -rf GPLsource.tar

