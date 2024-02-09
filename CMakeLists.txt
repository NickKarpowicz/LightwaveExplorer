if(UNIX AND NOT APPLE)
    if(MAKEFULL)
    #full Linux build; typical usage:
    #. /opt/intel/oneapi/setvars.sh
    #cmake --fresh -DMAKEFULL=TRUE -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_HOST_COMPILER=clang++ -DCMAKE_CUDA_COMPILER=nvcc -DCMAKE_CUDA_ARCHITECTURES=75 .. -G Ninja
        cmake_minimum_required(VERSION 3.23.0)
        project(LightwaveExplorer LANGUAGES C CXX CUDA)

        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        execute_process(COMMAND git clone https://github.com/kthohr/gcem)
        execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
        execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        find_package(CUDAToolkit)
        find_package(PkgConfig REQUIRED)
        find_package(fmt REQUIRED)
        find_package(IntelSYCL REQUIRED)
        find_package(TBB REQUIRED)
        find_package(MKL REQUIRED)
        pkg_check_modules(GTK4 gtk4)
        
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${GTK4_INCLUDE_DIRS})
        include_directories(build/dlib) 
        include_directories(build/gcem/include)
        include_directories(build/miniz)
        link_directories(${GTK4_LIBRARY_DIRS})
        add_library(miniz miniz/miniz.c)
        add_library(LightwaveExplorerCuda 
            Source/LightwaveExplorerCore.cu 
            Source/Devices/LightwaveExplorerCoreFP32.cu)
        
        add_compile_options(-std=c++20 -fp-model=precise -O3 -fopenmp -fsycl -w)
        add_link_options(-fsycl)
        add_executable(LightwaveExplorer 
            Source/Devices/LightwaveExplorerSYCLLinux.cpp 
            Source/Devices/LightwaveExplorerSYCLLinuxFP32.cpp 
            Source/Frontend/LightwaveExplorerFrontendGTK.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES} miniz)
        target_link_libraries(LightwaveExplorer fmt::fmt -lm)
        target_link_libraries(LightwaveExplorer
            ${MKL_ROOT}/lib/intel64/libmkl_sycl.a 
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
            -lsycl -lOpenCL)
        target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
        target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
        target_link_libraries(LightwaveExplorer TBB::tbb)

        configure_file(${CMAKE_SOURCE_DIR}/Source/BuildResources/makeLaunchScript.sh ${CMAKE_BINARY_DIR}/makeLaunchScript.sh COPYONLY)
        add_custom_command(TARGET LightwaveExplorer POST_BUILD COMMAND sh ${CMAKE_BINARY_DIR}/makeLaunchScript.sh)
        install(TARGETS LightwaveExplorer)
        install(PROGRAMS ${CMAKE_BINARY_DIR}/LightwaveExplorerLauncher.sh DESTINATION bin)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
        install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source//BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
        install(FILES Source/BuildResources/DesktopFileLauncher 
            DESTINATION share/applications 
            RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
        install(CODE "execute_process(COMMAND gtk-update-icon-cache -q -t -f ${CMAKE_INSTALL_PREFIX}/share/icons/hicolor)")
    elseif(MAKEONEAPI)
    #Oneapi Linux build without CUDA; typical usage:
    #. /opt/intel/oneapi/setvars.sh
    #cmake --fresh -DMAKEONEAPI=TRUE -DCMAKE_CXX_COMPILER=icpx .. -G Ninja
        cmake_minimum_required(VERSION 3.23.0)
        project(LightwaveExplorer LANGUAGES C CXX)

        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        execute_process(COMMAND git clone https://github.com/kthohr/gcem)
        execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
        execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        find_package(PkgConfig REQUIRED)
        find_package(fmt REQUIRED)
        find_package(IntelSYCL REQUIRED)
        find_package(TBB REQUIRED)
        find_package(MKL REQUIRED)
        
        pkg_check_modules(GTK4 gtk4)
        
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${GTK4_INCLUDE_DIRS})
        include_directories(build/dlib) 
        include_directories(build/gcem/include)
        include_directories(build/miniz)
        link_directories(${GTK4_LIBRARY_DIRS})
        add_library(miniz miniz/miniz.c)

        add_compile_options(-std=c++20 -DNOCUDA -fp-model=precise -O3 -fopenmp -fsycl -w -Wall)
        add_link_options(-fsycl)
        add_executable(LightwaveExplorer 
            Source/Devices/LightwaveExplorerSYCLLinux.cpp 
            Source/Devices/LightwaveExplorerSYCLLinuxFP32.cpp 
            Source/Frontend/LightwaveExplorerFrontendGTK.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
        target_link_libraries(LightwaveExplorer miniz)
        target_link_libraries(LightwaveExplorer fmt::fmt -lm)
        target_link_libraries(LightwaveExplorer
            ${MKL_ROOT}/lib/intel64/libmkl_sycl.a 
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
            -lsycl -lOpenCL)
        target_link_libraries(LightwaveExplorer TBB::tbb)

        configure_file(${CMAKE_SOURCE_DIR}/Source/BuildResources/makeLaunchScript.sh ${CMAKE_BINARY_DIR}/makeLaunchScript.sh COPYONLY)
        add_custom_command(TARGET LightwaveExplorer POST_BUILD COMMAND sh ${CMAKE_BINARY_DIR}/makeLaunchScript.sh)
        install(TARGETS LightwaveExplorer)
        install(PROGRAMS ${CMAKE_BINARY_DIR}/LightwaveExplorerLauncher.sh DESTINATION bin)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
        install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source//BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
        install(FILES Source/BuildResources/DesktopFileLauncher 
            DESTINATION share/applications 
            RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
        install(CODE "execute_process(COMMAND gtk-update-icon-cache -q -t -f ${CMAKE_INSTALL_PREFIX}/share/icons/hicolor)")
    elseif(MAKECUDAFLATPAK)
    #Flatpak using CUDA, and CPU calculations with FFTs in MKL. 
    #Current release version of the Flatpak.
        cmake_minimum_required(VERSION 3.23.0)
        project(LightwaveExplorer LANGUAGES CXX)   
        enable_language("CUDA")
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(GTK4 REQUIRED gtk4)
        find_package(CUDAToolkit REQUIRED)
        find_package(TBB PATHS /app/opt/intel/oneapi/tbb/latest/lib/cmake REQUIRED)
        find_package(MKL PATHS /app/opt/intel/oneapi/mkl/latest/lib/cmake REQUIRED)
        find_package(fmt)
        find_package(OpenMP)
        include_directories(${CUDA_INCLUDE_DIRS})
        include_directories(${GTK4_INCLUDE_DIRS})
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${MKL_ROOT}/include)
        add_library(LightwaveExplorerCuda 
            Source/LightwaveExplorerCore.cu 
            Source/Devices/LightwaveExplorerCoreFP32.cu)
        add_executable(LightwaveExplorer 
            Source/Frontend/LightwaveExplorerFrontendGTK.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_compile_options(LightwaveExplorer PRIVATE -std=c++20 -O3 -fopenmp -DNOSYCL)
        target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
        target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
        target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
        target_link_libraries(LightwaveExplorer fmt::fmt)
        target_link_libraries(LightwaveExplorer -lm -lminiz
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            ${MKL_ROOT}/../../compiler/latest/linux/compiler/lib/intel64/libiomp5.a
            )
        target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
        target_link_libraries(LightwaveExplorer TBB::tbb)
        add_executable(LightwaveExplorerNoCuda 
            Source/Frontend/LightwaveExplorerFrontendGTK.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_compile_options(LightwaveExplorerNoCuda PRIVATE -std=c++20 -O3 -fopenmp -DNOSYCL -DNOCUDA)
        target_link_libraries(LightwaveExplorerNoCuda ${GTK4_LIBRARIES})
        target_link_libraries(LightwaveExplorerNoCuda fmt::fmt)
        target_link_libraries(LightwaveExplorerNoCuda -lm
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            ${MKL_ROOT}/../../compiler/latest/linux/compiler/lib/intel64/libiomp5.a
            )
        target_link_libraries(LightwaveExplorerNoCuda TBB::tbb)
        install(TARGETS LightwaveExplorer LightwaveExplorerNoCuda DESTINATION bin)
        install(PROGRAMS ${CMAKE_SOURCE_DIR}/Source/BuildResources/flatpakLauncher.sh DESTINATION bin)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
        install(FILES Source/BuildResources/io.github.NickKarpowicz.LightwaveExplorer.metainfo.xml DESTINATION share/metainfo)
        install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source/BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
        install(FILES Source/BuildResources/DesktopFileFlatpak DESTINATION share/applications RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
        elseif(MAKECUDAINTELFLATPAK)
        #Flatpak using CUDA, SYCL, and CPU calculations with FFTs in MKL. 
        #Pre release version of the Flatpak.
            cmake_minimum_required(VERSION 3.23.0)
            project(LightwaveExplorer LANGUAGES CXX)   
            enable_language("CUDA")
            find_package(PkgConfig REQUIRED)
            pkg_check_modules(GTK4 REQUIRED gtk4)
            find_package(CUDAToolkit REQUIRED)
            find_package(IntelSYCL PATHS /app/opt/intel/oneapi/compiler/latest/linux REQUIRED)
            find_package(TBB PATHS /app/opt/intel/oneapi/tbb/latest/lib/cmake REQUIRED)
            find_package(MKL PATHS /app/opt/intel/oneapi/mkl/latest/lib/cmake /app/opt/intel/oneapi/tbb/latest/lib/cmake REQUIRED)
            find_package(oneDPL PATHS /app/opt/intel/oneapi/dpl/latest/lib/cmake REQUIRED)
            
            find_package(fmt)
            find_package(OpenMP)
            include_directories(../LightwaveExplorer)
            include_directories(${CUDA_INCLUDE_DIRS})
            include_directories(${GTK4_INCLUDE_DIRS})
            include_directories(${MKL_ROOT}/include/fftw)
            include_directories(${MKL_ROOT}/include)
            include_directories(/app/opt/intel/oneapi/dpl/latest/linux/include)
            add_library(LightwaveExplorerCuda 
                Source/LightwaveExplorerCore.cu 
                Source/Devices/LightwaveExplorerCoreFP32.cu)
            add_executable(LightwaveExplorer 
                Source/Frontend/LightwaveExplorerFrontendGTK.cpp 
                Source/LightwaveExplorerUtilities.cpp 
                Source/Devices/LightwaveExplorerCoreCPU.cpp 
                Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
                Source/Devices/LightwaveExplorerCoreCounter.cpp 
                Source/Devices/DlibLibraryComponents.cpp)
            target_compile_options(LightwaveExplorer PRIVATE -std=c++20 -fp-model=precise -O3 -fopenmp -DNOSYCL)
            target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
            target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
            target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
            target_link_libraries(LightwaveExplorer fmt::fmt)
            target_link_libraries(LightwaveExplorer -lm
                -Wl,--start-group 
                ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
                ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
                ${MKL_ROOT}/lib/intel64/libmkl_core.a 
                -Wl,--end-group
                ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
                )
            target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
            target_link_libraries(LightwaveExplorer TBB::tbb)

            add_link_options(-fsycl)
            add_executable(LightwaveExplorerNoCuda 
                Source/Frontend/LightwaveExplorerFrontendGTK.cpp 
                Source/LightwaveExplorerUtilities.cpp 
                Source/Devices/LightwaveExplorerCoreCPU.cpp 
                Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
                Source/Devices/LightwaveExplorerCoreCounter.cpp 
                Source/Devices/DlibLibraryComponents.cpp
                Source/Devices/LightwaveExplorerSYCLLinux.cpp 
                Source/Devices/LightwaveExplorerSYCLLinuxFP32.cpp 
                )
            target_compile_options(LightwaveExplorerNoCuda PRIVATE -std=c++20 -O3 -fp-model=precise -fopenmp -std=c++20 -fp-model=precise -O3 -fopenmp -fsycl -DNOCUDA)
            target_link_libraries(LightwaveExplorerNoCuda ${GTK4_LIBRARIES})
            target_link_libraries(LightwaveExplorerNoCuda fmt::fmt)
            target_link_libraries(LightwaveExplorerNoCuda -lm
                ${MKL_ROOT}/lib/intel64/libmkl_sycl.a 
                -Wl,--start-group 
                ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
                ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
                ${MKL_ROOT}/lib/intel64/libmkl_core.a 
                -Wl,--end-group
                ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
                -lsycl -lOpenCL)
            target_link_libraries(LightwaveExplorerNoCuda TBB::tbb)
            install(TARGETS LightwaveExplorer LightwaveExplorerNoCuda DESTINATION bin)
            install(PROGRAMS ${CMAKE_SOURCE_DIR}/Source/BuildResources/flatpakLauncher.sh DESTINATION bin)
            install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
            install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
            install(FILES Source/BuildResources/io.github.NickKarpowicz.LightwaveExplorer.metainfo.xml DESTINATION share/metainfo)
            install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source/BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
            install(FILES Source/BuildResources/DesktopFileFlatpak DESTINATION share/applications RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
    elseif(CLICPU)
        cmake_minimum_required(VERSION 3.16)
        project(LightwaveExplorer LANGUAGES CXX)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(FFTW REQUIRED fftw3)
        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
        execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        include_directories(build/dlib) 
        include_directories(build/miniz)
        include_directories(${FFTW_INCLUDE_DIRS})
        link_directories(${FFTW_LIBRARY_DIRS})

        add_compile_options(-std=c++20 -O3 -fopenmp -DCPUONLY -DNOSYCL)
        add_executable(LightwaveExplorer 
            Source/LightwaveExplorerCommandLineMain.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer -lgomp miniz)
        target_link_libraries(LightwaveExplorer ${FFTW_LIBRARIES})
    elseif(CLICUDA)
    #full Linux build; typical usage:
    #. /opt/intel/oneapi/setvars.sh
    #cmake --fresh -DMAKEFULL=TRUE -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_HOST_COMPILER=clang++ -DCMAKE_CUDA_COMPILER=nvcc -DCMAKE_CUDA_ARCHITECTURES=75 .. -G Ninja
        cmake_minimum_required(VERSION 3.23.0)
        project(LightwaveExplorer LANGUAGES C CXX CUDA)

        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
        execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        find_package(CUDAToolkit)
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${MKL_ROOT}/include)
        include_directories(build/dlib) 
        include_directories(build/miniz)
        add_library(miniz miniz/miniz.c)
        add_compile_options(-std=c++20 -O3 -DRUNONCUDA -w)
        add_executable(LightwaveExplorer 
        Source/LightwaveExplorerCore.cu
        Source/LightwaveExplorerCommandLineMain.cu
        Source/LightwaveExplorerUtilities.cpp 
        Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer miniz)
        target_link_libraries(LightwaveExplorer -lm)
        target_link_libraries(LightwaveExplorer
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_sequential.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            )
        target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
    else()
    #if nothing specified, build GPL version based on FFTW
        cmake_minimum_required(VERSION 3.16)
        project(LightwaveExplorer LANGUAGES CXX)
        find_package(PkgConfig REQUIRED)
        find_package(OpenMP)
        pkg_check_modules(GTK4 REQUIRED gtk4)
        pkg_check_modules(FFTW REQUIRED fftw3)
        pkg_check_modules(FFTWF REQUIRED fftw3f)
        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        execute_process(COMMAND git clone https://github.com/kthohr/gcem)
        execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
        execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        include_directories(build/dlib) 
        include_directories(build/gcem/include)
        include_directories(build/miniz)
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

        add_compile_options(-std=c++20 -O3 -fopenmp -DCPUONLY -DNOSYCL -DNDEBUG)
        execute_process(COMMAND echo "NOTE: Building CPU only version since OneAPI not specified. This version uses FFTW for performing Fourier transforms and is licensed under the GPL v3.0")
        add_executable(LightwaveExplorer 
            Source/Frontend/LightwaveExplorerFrontendGTK.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
        target_link_libraries(LightwaveExplorer OpenMP::OpenMP_CXX)
        target_link_libraries(LightwaveExplorer ${FFTW_LIBRARIES})
        target_link_libraries(LightwaveExplorer ${FFTWF_LIBRARIES})
        target_link_libraries(LightwaveExplorer fmt::fmt miniz)

        
        install(TARGETS LightwaveExplorer)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
        install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source/BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
        install(FILES Source/BuildResources/DesktopFileNoLauncher DESTINATION share/applications RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
        install(CODE "execute_process(COMMAND gtk-update-icon-cache -q -t -f ${CMAKE_INSTALL_PREFIX}/share/icons/hicolor)")
    endif()
endif()

if(APPLE)
    cmake_minimum_required(VERSION 3.16)
    project(LightwaveExplorer)
    find_package(PkgConfig REQUIRED)
    find_package(fmt CONFIG REQUIRED)
    find_package(gcem CONFIG REQUIRED)
    find_package(miniz CONFIG REQUIRED)
    pkg_check_modules(GTK4 gtk4)
    include_directories(${CMAKE_SOURCE_DIR}/../dlib) 
    include_directories(${CMAKE_SOURCE_DIR}/../fftw/usr/local/include) 
    include_directories(${GTK4_INCLUDE_DIRS})
    link_directories(${GTK4_LIBRARY_DIRS})
    add_compile_options(-std=c++20 -O3 -DCPUONLY -DNOSYCL -Wall -x objective-c++)
    add_executable(LightwaveExplorer MACOSX_BUNDLE 
        Source/Frontend/LightwaveExplorerFrontendGTK.cpp 
        Source/LightwaveExplorerUtilities.cpp 
        Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
        Source/Devices/LightwaveExplorerCoreCPU.cpp 
        Source/Devices/LightwaveExplorerCoreCounter.cpp 
        Source/Devices/DlibLibraryComponents.cpp)
    set_target_properties(LightwaveExplorer  PROPERTIES
        BUNDLE True
        MACOSX_BUNDLE_BUNDLE_NAME LightwaveExplorer
        MACOSX_BUNDLE_VERSION "0.1"
        MACOSX_BUNDLE_VERSION_STRING "0.1"
        MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/Source/BuildResources/macplistbase.plist
        MACOSX_BUNDLE_ICON_FILE ${CMAKE_SOURCE_DIR}/Source/BuildResources/AppIcon.icns
    )
    target_link_libraries(LightwaveExplorer fmt::fmt)
    target_link_libraries(LightwaveExplorer miniz::miniz)
    target_link_libraries(LightwaveExplorer ${CMAKE_SOURCE_DIR}/../fftw/usr/local/lib/libfftw3_threads.a)
    target_link_libraries(LightwaveExplorer ${CMAKE_SOURCE_DIR}/../fftw/usr/local/lib/libfftw3f_threads.a)
        target_link_libraries(LightwaveExplorer ${CMAKE_SOURCE_DIR}/../fftw/usr/local/lib/libfftw3.a)
    target_link_libraries(LightwaveExplorer ${CMAKE_SOURCE_DIR}/../fftw/usr/local/lib/libfftw3f.a)
    target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES})
    target_link_libraries(LightwaveExplorer "-framework Cocoa")
endif()