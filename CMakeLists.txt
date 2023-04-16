
if(UNIX AND NOT APPLE)
    if(DEFINED ONEAPI_ROOT)
        cmake_minimum_required(VERSION 3.23.0)
        set(SYCL_INCLUDE_DIR "${ONEAPI_ROOT}/compiler/latest/linux/include/sycl")
        set(SYCL_LIBRARY_DIR "${ONEAPI_ROOT}/compiler/latest/linux/lib")
        set(CMAKE_CUDA_ARCHITECTURES 75)
        set(MKLROOT "${ONEAPI_ROOT}/mkl/latest")
        project(LightwaveExplorer LANGUAGES CXX CUDA)
        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        include_directories(build/dlib ../LightwaveExplorer ${MKLROOT}/include/fftw)
        find_package(CUDA REQUIRED)
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -gencode=arch=compute_75,code=\"sm_75,compute_75\" -O3 -ccbin g++ -D NOCUDAMAIN)
        cuda_compile(LightwaveExplorerCudaObjects LightwaveExplorerCore.cu LightwaveExplorerDevices/LightwaveExplorerCoreFP32.cu)

        find_package(PkgConfig REQUIRED)
        find_package(CUDAToolkit)
        find_package(IntelDPCPP)
        find_package(fmt REQUIRED)
        pkg_check_modules(GTK4 gtk4)
        pkg_check_modules(LIOMP5 openmp)
        pkg_check_modules(TBB tbb)
        
        include_directories(${CUDA_INCLUDE_DIRS})
        include_directories(${GTK4_INCLUDE_DIRS})
        link_directories(${CUDA_LIBRARY_DIRS})
        link_directories(${CUDATookit_LIBRARY_DIRS})
        link_directories(${LIOMP5_LIBRARY_DIRS})
        link_directories(${TBB_LIBRARY_DIRS})
        link_directories(${SYCL_LIBRARY_DIR})
        link_directories(${GTK4_LIBRARY_DIRS})

        add_compile_options(-std=c++20 -fp-model=precise -O3 -fopenmp -fsycl -w)
        add_executable(LightwaveExplorer ${LightwaveExplorerCudaObjects} LightwaveExplorerGTK/LightwaveExplorerDPCPPlib.cpp LightwaveExplorerGTK/LightwaveExplorerDPCPPlibFP32.cpp LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp LightwaveExplorerDevices/LightwaveExplorerUtilities.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPU.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPUFP32.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCounter.cpp LightwaveExplorerDevices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES} ${CUDA_LIBRARIES} ${CUDAToolkit_LIBRARIES})
        target_link_libraries(LightwaveExplorer fmt::fmt)
        target_link_libraries(LightwaveExplorer -lcufft -lnvidia-ml -lm)
        target_link_libraries(LightwaveExplorer
            ${MKLROOT}/lib/intel64/libmkl_sycl.a 
            -Wl,--start-group 
            ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKLROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKLROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            -lsycl
            -lOpenCL
        )
        target_link_libraries(LightwaveExplorer ${LIOMP5_LIBRARIES})
        target_link_libraries(LightwaveExplorer ${TBB_LIBRARIES})
        
        install(TARGETS LightwaveExplorer)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
    elseif(DEFINED GPLFLATPAK)
        cmake_minimum_required(VERSION 3.16)
        project(LightwaveExplorer)
        find_package(PkgConfig REQUIRED)
        find_package(OpenMP)
        find_package(fmt)
        pkg_check_modules(GTK4 REQUIRED gtk4)
        pkg_check_modules(FFTW REQUIRED fftw3)
        pkg_check_modules(FFTWF REQUIRED fftw3f)
        include_directories(build/dlib) 
        include_directories(${FFTW_INCLUDE_DIRS})
        include_directories(${GTK4_INCLUDE_DIRS})
        link_directories(${FFTW_LIBRARY_DIRS})
        link_directories(${GTK4_LIBRARY_DIRS})
   

        add_compile_options(-std=c++20 -O3 -fopenmp -D CPUONLY -w)
        execute_process(COMMAND echo "NOTE: Building CPU only version since OneAPI not specified. This version uses FFTW for performing Fourier transforms and is licensed under the GPL v3.0")
        add_executable(LightwaveExplorer LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp LightwaveExplorerDevices/LightwaveExplorerUtilities.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPU.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPUFP32.cpp  LightwaveExplorerDevices/LightwaveExplorerCoreCounter.cpp LightwaveExplorerDevices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
        target_link_libraries(LightwaveExplorer OpenMP::OpenMP_CXX)
        target_link_libraries(LightwaveExplorer ${FFTW_LIBRARIES})
        target_link_libraries(LightwaveExplorer ${FFTWF_LIBRARIES})
        target_link_libraries(LightwaveExplorer fmt::fmt)

        install(TARGETS LightwaveExplorer)
        install(TARGETS LightwaveExplorer)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)

    else()
        cmake_minimum_required(VERSION 3.16)
        project(LightwaveExplorer)
        find_package(PkgConfig REQUIRED)
        find_package(OpenMP)
        pkg_check_modules(GTK4 REQUIRED gtk4)
        pkg_check_modules(FFTW REQUIRED fftw3)
        pkg_check_modules(FFTWF REQUIRED fftw3f)
        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        include_directories(build/dlib) 
        include_directories(${FFTW_INCLUDE_DIRS})
        include_directories(${GTK4_INCLUDE_DIRS})
        link_directories(${FFTW_LIBRARY_DIRS})
        link_directories(${GTK4_LIBRARY_DIRS})
        
        
        include(FetchContent)
        FetchContent_Declare(fmt
            GIT_REPOSITORY https://github.com/fmtlib/fmt.git
            GIT_TAG        9.1.0
        )
        FetchContent_MakeAvailable(fmt)

        add_compile_options(-std=c++20 -O3 -fopenmp -D CPUONLY -w)
        execute_process(COMMAND echo "NOTE: Building CPU only version since OneAPI not specified. This version uses FFTW for performing Fourier transforms and is licensed under the GPL v3.0")
        add_executable(LightwaveExplorer LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp LightwaveExplorerDevices/LightwaveExplorerUtilities.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPU.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPUFP32.cpp  LightwaveExplorerDevices/LightwaveExplorerCoreCounter.cpp LightwaveExplorerDevices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
        target_link_libraries(LightwaveExplorer OpenMP::OpenMP_CXX)
        target_link_libraries(LightwaveExplorer ${FFTW_LIBRARIES})
        target_link_libraries(LightwaveExplorer ${FFTWF_LIBRARIES})
        target_link_libraries(LightwaveExplorer fmt::fmt)

        install(TARGETS LightwaveExplorer)
        install(TARGETS LightwaveExplorer)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
    endif()
endif()

if(APPLE)
    cmake_minimum_required(VERSION 3.16)
    project(LightwaveExplorer)
    execute_process(COMMAND brew --prefix OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE BREW_PREFIX)
    list(APPEND CMAKE_MODULE_PATH "${BREW_PREFIX}/share/cmake/Modules")
    execute_process(COMMAND brew --prefix llvm OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE LLVM_PREFIX)
    set(CMAKE_CXX_COMPILER "${LLVM_PREFIX}/bin/clang++")
    set(CMAKE_C_COMPILER "${LLVM_PREFIX}/bin/clang")

    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GTK4 REQUIRED gtk4)
    pkg_check_modules(FFTW REQUIRED fftw3)
    pkg_check_modules(FFTWF REQUIRED fftw3f)
    find_library(OpenMP_LIBRARY NAMES libomp.a HINTS ${BREW_PREFIX}/opt/libomp/lib/)

    execute_process(COMMAND git clone https://github.com/davisking/dlib)
    include_directories(build/dlib) 
    include_directories(${FFTW_INCLUDE_DIRS})
    include_directories(${GTK4_INCLUDE_DIRS})

    link_directories(${FFTW_LIBRARY_DIRS})
    link_directories(${GTK4_LIBRARY_DIRS})
    
    add_compile_options(-std=c++20 -O3 -fexperimental-library -stdlib=libc++ -fopenmp -D CPUONLY -Wall -Wno-deprecated-declarations)
    add_executable(LightwaveExplorer MACOSX_BUNDLE LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.mm LightwaveExplorerDevices/LightwaveExplorerUtilities.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPUFP32.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPU.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCounter.cpp LightwaveExplorerDevices/DlibLibraryComponents.cpp)
    set_target_properties(LightwaveExplorer  PROPERTIES
        BUNDLE True
        MACOSX_BUNDLE_BUNDLE_NAME LightwaveExplorer
        MACOSX_BUNDLE_VERSION "0.1"
        MACOSX_BUNDLE_VERSION_STRING "0.1"
        MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/BuildResources/macplistbase.plist
        MACOSX_BUNDLE_ICON_FILE ${CMAKE_SOURCE_DIR}/BuildResources/AppIcon.icns
    )
    target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES} ${FFTW_LIBRARIES} ${FFTWF_LIBRARIES})
    target_link_libraries(LightwaveExplorer ${OpenMP_LIBRARY} ${LLVM_PREFIX}/lib/c++/libc++.a ${LLVM_PREFIX}/lib/c++/libc++abi.a ${LLVM_PREFIX}/lib/c++/libc++experimental.a "-framework Cocoa")
endif()