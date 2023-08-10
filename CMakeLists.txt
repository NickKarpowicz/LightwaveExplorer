if(UNIX AND NOT APPLE)
    # cmake --fresh -DCMAKE_CUDA_COMPILER=nvcc -DCMAKE_CUDA_HOST_COMPILER=clang++-15 -DCMAKE_CUDA_ARCHITECTURES=75 ..
    if(DEFINED CMAKE_CUDA_COMPILER)
        cmake_minimum_required(VERSION 3.23.0)
        project(LightwaveExplorer LANGUAGES CXX CUDA)      
        enable_language("CUDA")
        find_package(PkgConfig REQUIRED)
        find_package(OpenMP)
        pkg_check_modules(GTK4 REQUIRED gtk4)
        find_package(CUDAToolkit)
        find_package(MKL)
        include(FetchContent)
        FetchContent_Declare(fmt
            GIT_REPOSITORY https://github.com/fmtlib/fmt.git
            GIT_TAG        9.1.0
        )
        FetchContent_MakeAvailable(fmt)

        configure_file(${CMAKE_SOURCE_DIR}/BuildResources/makeLinuxIcons.sh ${CMAKE_BINARY_DIR}/makeLinuxIcons.sh COPYONLY)
        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        execute_process(COMMAND git clone https://github.com/kthohr/gcem)
        include_directories(build/dlib) 
        include_directories(build/gcem/include)
        include_directories(../LightwaveExplorer)
        include_directories(${CUDA_INCLUDE_DIRS})
        include_directories(${GTK4_INCLUDE_DIRS})
        include_directories(${MKL_ROOT}/include/fftw)
        add_library(LightwaveExplorerCuda LightwaveExplorerCore.cu LightwaveExplorerDevices/LightwaveExplorerCoreFP32.cu)
        add_compile_options(-std=c++20 -O3 -fopenmp -D NOSYCL)
        add_executable(LightwaveExplorer LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp LightwaveExplorerDevices/LightwaveExplorerUtilities.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPU.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPUFP32.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCounter.cpp LightwaveExplorerDevices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
        target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
        target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
        target_link_libraries(LightwaveExplorer fmt::fmt -lm -lgomp)
        target_link_libraries(LightwaveExplorer
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            )
        target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
        add_custom_command(TARGET LightwaveExplorer POST_BUILD COMMAND ${CMAKE_BINARY_DIR}/makeLinuxIcons.sh)
        install(TARGETS LightwaveExplorer)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
        install(DIRECTORY ${CMAKE_BINARY_DIR}/icons DESTINATION share FILES_MATCHING PATTERN "*")
        install(FILES ${CMAKE_BINARY_DIR}/io.github.NickKarpowicz.LightwaveExplorer.desktop DESTINATION share/applications)
        install(CODE "execute_process(COMMAND gtk-update-icon-cache -q -t -f ${CMAKE_INSTALL_PREFIX}/share/icons/hicolor)")
    elseif(DEFINED ONEAPI_ROOT)
    #Linux build; typical usage:
    #cmake --fresh -DONEAPI_ROOT=${ONEAPI_ROOT} -DCMAKE_CXX_COMPILER=icpx ..
        cmake_minimum_required(VERSION 3.23.0)
        project(LightwaveExplorer LANGUAGES CXX)
        #set(SYCL_INCLUDE_DIR "${ONEAPI_ROOT}/compiler/latest/linux/include/sycl")
        #set(SYCL_LIBRARY_DIR "${ONEAPI_ROOT}/compiler/latest/linux/lib")
        set(MKLROOT "${ONEAPI_ROOT}/mkl/latest")

        configure_file(${CMAKE_SOURCE_DIR}/BuildResources/makeLinuxIcons.sh ${CMAKE_BINARY_DIR}/makeLinuxIcons.sh COPYONLY)
        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        execute_process(COMMAND git clone https://github.com/kthohr/gcem)
        include_directories(build/dlib) 
        include_directories(build/gcem/include)
        include_directories(../LightwaveExplorer ${MKLROOT}/include/fftw)

        find_package(PkgConfig REQUIRED)
        find_package(fmt REQUIRED)
        find_package(IntelSYCL REQUIRED)
        pkg_check_modules(GTK4 gtk4)
        pkg_check_modules(LIOMP5 openmp)
        pkg_check_modules(TBB tbb)
        include_directories(../LightwaveExplorer)
        include_directories(${GTK4_INCLUDE_DIRS})
        link_directories(${LIOMP5_LIBRARY_DIRS})
        link_directories(${TBB_LIBRARY_DIRS})
        link_directories(${SYCL_LIBRARY_DIR})
        link_directories(${GTK4_LIBRARY_DIRS})

        add_compile_options(-std=c++20 -D NOCUDA -fp-model=precise -O3 -fopenmp -fsycl -fPIC -w)
        add_executable(LightwaveExplorer LightwaveExplorerGTK/LightwaveExplorerDPCPPlib.cpp LightwaveExplorerGTK/LightwaveExplorerDPCPPlibFP32.cpp LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp LightwaveExplorerDevices/LightwaveExplorerUtilities.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPU.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPUFP32.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCounter.cpp LightwaveExplorerDevices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
        target_link_libraries(LightwaveExplorer fmt::fmt -lm)
        target_link_libraries(LightwaveExplorer
            ${MKLROOT}/lib/intel64/libmkl_sycl.a 
            -Wl,--start-group 
            ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKLROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKLROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            )
        target_link_libraries(LightwaveExplorer -lsycl ${SYCL_LIBRARY_DIR}/libOpenCL.so)
        target_link_libraries(LightwaveExplorer ${LIOMP5_LIBRARIES})
        target_link_libraries(LightwaveExplorer ${TBB_LIBRARIES})
        add_custom_command(TARGET LightwaveExplorer POST_BUILD COMMAND ${CMAKE_BINARY_DIR}/makeLinuxIcons.sh)
        install(TARGETS LightwaveExplorer)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
        install(DIRECTORY ${CMAKE_BINARY_DIR}/icons DESTINATION share FILES_MATCHING PATTERN "*")
        install(FILES ${CMAKE_BINARY_DIR}/io.github.NickKarpowicz.LightwaveExplorer.desktop DESTINATION share/applications)
        install(CODE "execute_process(COMMAND gtk-update-icon-cache -q -t -f ${CMAKE_INSTALL_PREFIX}/share/icons/hicolor)")
    else()
        cmake_minimum_required(VERSION 3.16)
        project(LightwaveExplorer LANGUAGES CXX)
        configure_file(${CMAKE_SOURCE_DIR}/BuildResources/makeLinuxIcons.sh ${CMAKE_BINARY_DIR}/makeLinuxIcons.sh COPYONLY)
        find_package(PkgConfig REQUIRED)
        find_package(OpenMP)
        pkg_check_modules(GTK4 REQUIRED gtk4)
        pkg_check_modules(FFTW REQUIRED fftw3)
        pkg_check_modules(FFTWF REQUIRED fftw3f)
        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        execute_process(COMMAND git clone https://github.com/kthohr/gcem)
        include_directories(build/dlib) 
        include_directories(build/gcem/include)
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

        add_custom_command(TARGET LightwaveExplorer POST_BUILD COMMAND ${CMAKE_BINARY_DIR}/makeLinuxIcons.sh)
        install(TARGETS LightwaveExplorer)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
        install(DIRECTORY ${CMAKE_BINARY_DIR}/icons DESTINATION share FILES_MATCHING PATTERN "*")
        install(FILES ${CMAKE_BINARY_DIR}/io.github.NickKarpowicz.LightwaveExplorer.desktop DESTINATION share/applications)
        install(CODE "execute_process(COMMAND gtk-update-icon-cache -q -t -f ${CMAKE_INSTALL_PREFIX}/share/icons/hicolor)")
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
    pkg_check_modules(FFTW REQUIRED fftw3)
    pkg_check_modules(FFTWF REQUIRED fftw3f)
    pkg_check_modules(GTK4 gtk4)
    find_library(OpenMP_LIBRARY NAMES libomp.a HINTS ${BREW_PREFIX}/opt/libomp/lib/)

    execute_process(COMMAND git clone https://github.com/davisking/dlib)
    execute_process(COMMAND git clone https://github.com/kthohr/gcem)
    include_directories(build/dlib) 
    include_directories(build/gcem/include)
    include_directories(${FFTW_INCLUDE_DIRS})
    include_directories(${GTK4_INCLUDE_DIRS})
    link_directories(${FFTW_LIBRARY_DIRS})
    link_directories(${GTK4_LIBRARY_DIRS})
    add_compile_options(-std=c++20 -O3 -fexperimental-library -stdlib=libc++ -fopenmp -D CPUONLY -Wall -Wno-deprecated-declarations -x objective-c++)
    add_executable(LightwaveExplorer MACOSX_BUNDLE LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp LightwaveExplorerDevices/LightwaveExplorerUtilities.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPUFP32.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCPU.cpp LightwaveExplorerDevices/LightwaveExplorerCoreCounter.cpp LightwaveExplorerDevices/DlibLibraryComponents.cpp)
    set_target_properties(LightwaveExplorer  PROPERTIES
        BUNDLE True
        MACOSX_BUNDLE_BUNDLE_NAME LightwaveExplorer
        MACOSX_BUNDLE_VERSION "0.1"
        MACOSX_BUNDLE_VERSION_STRING "0.1"
        MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/BuildResources/macplistbase.plist
        MACOSX_BUNDLE_ICON_FILE ${CMAKE_SOURCE_DIR}/BuildResources/AppIcon.icns
    )
    target_link_libraries(LightwaveExplorer ${FFTW_LIBRARIES} ${FFTWF_LIBRARIES})
    target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
    target_link_libraries(LightwaveExplorer ${OpenMP_LIBRARY} ${LLVM_PREFIX}/lib/c++/libc++.a ${LLVM_PREFIX}/lib/c++/libc++abi.a ${LLVM_PREFIX}/lib/c++/libc++experimental.a "-framework Cocoa")
endif()