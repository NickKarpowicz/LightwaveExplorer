cmake_minimum_required(VERSION 3.23.0)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(UNIX AND NOT APPLE)
    if(MAKEFULL)   
        project(LightwaveExplorer LANGUAGES C CXX CUDA)
        find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
        set(CMAKE_AUTOMOC ON)
        find_package(fmt REQUIRED)
        find_package(OpenMP REQUIRED)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(CAIRO REQUIRED cairo)
        find_package(CUDAToolkit REQUIRED)
        find_package(TBB REQUIRED)
        find_package(MKL REQUIRED)

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/dlib)
        message("Using existing dlib clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/davisking/dlib)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/gcem)
            message("Using existing gcem clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/kthohr/gcem)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/miniz)
            message("Using existing miniz download")
        else()
            execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
            execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        endif()
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${CAIRO_INCLUDE_DIRS})
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/dlib) 
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/gcem/include)
        include_directories(${CMAKE_CURRENT_BINARY_DIR})

        add_library(miniz miniz/miniz.c)

        add_library(LightwaveExplorerCuda 
            Source/LightwaveExplorerCore.cu 
            Source/Devices/LightwaveExplorerCoreFP32.cu)

        add_compile_options(-fp-model=precise -O3 ${OpenMP_CXX_FLAGS} -fsycl)
        add_link_options(-fsycl)
        add_executable(LightwaveExplorer 
            Source/Devices/LightwaveExplorerSYCLLinux.cpp 
            Source/Devices/LightwaveExplorerSYCLLinuxFP32.cpp 
            Source/Frontend/LightwaveExplorerFrontendQt.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer Qt6::Widgets Qt6::DBus)
        target_link_libraries(LightwaveExplorer fmt::fmt -lm miniz)
        target_link_libraries(LightwaveExplorer
            ${MKL_ROOT}/lib/intel64/libmkl_sycl.a 
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            -lsycl -lOpenCL)
        target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
        target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
        target_link_libraries(LightwaveExplorer TBB::tbb ${OpenMP_CXX_LIBRARIES})
        target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})
        
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
    elseif(MAKECUDA)   
        project(LightwaveExplorer LANGUAGES C CXX CUDA)
        find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
        set(CMAKE_AUTOMOC ON)
        find_package(fmt REQUIRED)
        find_package(PkgConfig REQUIRED)
        find_package(OpenMP REQUIRED)
        pkg_check_modules(CAIRO REQUIRED cairo)
        find_package(CUDAToolkit REQUIRED)
        find_package(TBB REQUIRED)
        
        if(USEFFTW)
            pkg_check_modules(FFTW REQUIRED fftw3)
            pkg_check_modules(FFTWF REQUIRED fftw3f)
            include_directories(${FFTW_INCLUDE_DIRS})
            link_directories(${FFTW_LIBRARY_DIRS})
        else()
            find_package(MKL REQUIRED)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/dlib)
        message("Using existing dlib clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/davisking/dlib)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/gcem)
            message("Using existing gcem clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/kthohr/gcem)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/miniz)
            message("Using existing miniz download")
        else()
            execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
            execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        endif()
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${CAIRO_INCLUDE_DIRS})
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/dlib) 
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/gcem/include)
        include_directories(${CMAKE_CURRENT_BINARY_DIR})

        add_library(miniz miniz/miniz.c)

        add_library(LightwaveExplorerCuda 
            Source/LightwaveExplorerCore.cu 
            Source/Devices/LightwaveExplorerCoreFP32.cu)
        add_definitions(-DNOSYCL)
        add_compile_options(-O3 ${OpenMP_CXX_FLAGS} -Wpedantic)
        if(USEFFTW)
            add_definitions(-DUSEFFTW)
        endif()
        add_executable(LightwaveExplorer 
            Source/Frontend/LightwaveExplorerFrontendQt.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer Qt6::Widgets Qt6::DBus)
        target_link_libraries(LightwaveExplorer fmt::fmt -lm miniz)
        target_link_libraries(LightwaveExplorer ${OpenMP_CXX_LIBRARIES})

        if(USEFFTW)
            target_link_libraries(LightwaveExplorer ${FFTW_LIBRARIES} ${FFTWF_LIBRARIES})
        else()
            target_link_libraries(LightwaveExplorer
                -Wl,--start-group 
                ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
                ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.a 
                ${MKL_ROOT}/lib/intel64/libmkl_core.a
                -Wl,--end-group)
        endif()

        target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
        target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
        target_link_libraries(LightwaveExplorer TBB::tbb)
        target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})
        
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
    elseif(MAKESYCL)   
        project(LightwaveExplorer LANGUAGES C CXX)
        find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
        set(CMAKE_AUTOMOC ON)
        find_package(fmt REQUIRED)
        find_package(OpenMP REQUIRED)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(CAIRO REQUIRED cairo)
        find_package(TBB REQUIRED)
        find_package(MKL REQUIRED)
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${CAIRO_INCLUDE_DIRS})
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/dlib) 
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/gcem/include)
        include_directories(${CMAKE_CURRENT_BINARY_DIR})

        add_library(miniz miniz/miniz.c)
        add_definitions(-DNOCUDA)
        add_compile_options(-fp-model=precise -O3 ${OpenMP_CXX_FLAGS} -fsycl -w)
        add_link_options(-fsycl)
        add_executable(LightwaveExplorer 
            Source/Devices/LightwaveExplorerSYCLLinux.cpp 
            Source/Devices/LightwaveExplorerSYCLLinuxFP32.cpp 
            Source/Frontend/LightwaveExplorerFrontendQt.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer Qt6::Widgets Qt6::DBus)
        target_link_libraries(LightwaveExplorer fmt::fmt -lm miniz)
        target_link_libraries(LightwaveExplorer
            ${MKL_ROOT}/lib/intel64/libmkl_sycl.a 
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            -lsycl -lOpenCL)
        target_link_libraries(LightwaveExplorer TBB::tbb)
        target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
        
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
    elseif(MAKEGTK)
    #full Linux build; typical usage:
    #. /opt/intel/oneapi/setvars.sh
    #cmake --fresh -DMAKEFULL=TRUE -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_HOST_COMPILER=clang++ -DCMAKE_CUDA_COMPILER=nvcc -DCMAKE_CUDA_ARCHITECTURES=75 .. -G Ninja
        project(LightwaveExplorer LANGUAGES C CXX CUDA)


        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/dlib)
            message("Using existing dlib clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/davisking/dlib)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/gcem)
            message("Using existing gcem clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/kthohr/gcem)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/miniz)
            message("Using existing miniz download")
        else()
            execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
            execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        endif()
        
        find_package(CUDAToolkit)
        find_package(PkgConfig REQUIRED)
        find_package(fmt REQUIRED)
        find_package(IntelSYCL REQUIRED)
        find_package(TBB REQUIRED)
        find_package(MKL REQUIRED)
        find_package(OpenMP REQUIRED)
        pkg_check_modules(GTK4 gtk4)
        
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${GTK4_INCLUDE_DIRS})
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/dlib) 
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/gcem/include)
        include_directories(${CMAKE_CURRENT_BINARY_DIR})
        link_directories(${GTK4_LIBRARY_DIRS})
        add_library(miniz miniz/miniz.c)
        add_library(LightwaveExplorerCuda 
            Source/LightwaveExplorerCore.cu 
            Source/Devices/LightwaveExplorerCoreFP32.cu)
        
        add_compile_options(-fp-model=precise -O3 ${OpenMP_CXX_FLAGS} -fsycl -w)
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
    elseif(MAKECUDAFLATPAK)
    #Flatpak using CUDA, and CPU calculations with FFTs in MKL. 
    #Current release version of the Flatpak.
        project(LightwaveExplorer LANGUAGES CXX)   
        enable_language("CUDA")
        
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(CAIRO REQUIRED cairo)
        find_package(CUDAToolkit REQUIRED)
        find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
        set(CMAKE_AUTOMOC ON)
        find_package(OpenMP REQUIRED)
        find_package(TBB PATHS /app/opt/intel/oneapi/tbb/latest/lib/cmake REQUIRED)
        find_package(MKL PATHS /app/opt/intel/oneapi/mkl/latest/lib/cmake)
        find_package(fmt)
        find_package(miniz)
        include_directories(${CUDA_INCLUDE_DIRS})
        include_directories(${Qt6_INCLUDE_DIRS})
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${MKL_ROOT}/include)
        include_directories(${CAIRO_INCLUDE_DIRS})
        add_definitions(-DLWEFLATPAK)
        add_library(LightwaveExplorerCuda 
            Source/LightwaveExplorerCore.cu 
            Source/Devices/LightwaveExplorerCoreFP32.cu)
        add_executable(LightwaveExplorer 
            Source/Frontend/LightwaveExplorerFrontendQt.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        set_target_properties(LightwaveExplorer PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
        target_compile_options(LightwaveExplorer PRIVATE -O3 ${OpenMP_CXX_FLAGS} -DNOSYCL -DLWEFLATPAK)
        target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
        target_link_libraries(LightwaveExplorer CUDA::cudart_static CUDA::cufft_static CUDA::nvml)
        target_link_libraries(LightwaveExplorer Qt6::Widgets Qt6::DBus)
        target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})
        target_link_libraries(LightwaveExplorer fmt::fmt miniz::miniz)
        target_link_libraries(LightwaveExplorer -lm
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
            )
        target_link_libraries(LightwaveExplorer TBB::tbb)
        add_executable(LightwaveExplorerNoCuda 
            Source/Frontend/LightwaveExplorerFrontendQt.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_compile_options(LightwaveExplorerNoCuda PRIVATE ${OpenMP_CXX_FLAGS} -DNOSYCL -DNOCUDA -DLWEFLATPAK)
        target_link_libraries(LightwaveExplorerNoCuda Qt6::Widgets Qt6::DBus)
        target_link_libraries(LightwaveExplorerNoCuda ${CAIRO_LIBRARIES})
        target_link_libraries(LightwaveExplorerNoCuda fmt::fmt miniz::miniz)
        target_link_libraries(LightwaveExplorerNoCuda -lm
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
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
        #Flatpak using CUDA, and CPU calculations with FFTs in MKL. 
        #Current release version of the Flatpak.
            project(LightwaveExplorer LANGUAGES CXX)   
            enable_language("CUDA")
            
            find_package(PkgConfig REQUIRED)
            pkg_check_modules(CAIRO REQUIRED cairo)
            find_package(CUDAToolkit REQUIRED)
            find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
            set(CMAKE_AUTOMOC ON)
            find_package(OpenMP REQUIRED)
            find_package(TBB PATHS /app/opt/intel/oneapi/tbb/latest/lib/cmake REQUIRED)
            find_package(MKL PATHS /app/opt/intel/oneapi/mkl/latest/lib/cmake REQUIRED)
            find_package(oneDPL PATHS /app/opt/intel/oneapi/dpl/latest/lib/cmake REQUIRED)
            find_package(fmt)
            find_package(miniz)
            include_directories(${CUDA_INCLUDE_DIRS})
            include_directories(${Qt6_INCLUDE_DIRS})
            include_directories(${MKL_ROOT}/include/fftw)
            include_directories(${MKL_ROOT}/include)
            include_directories(${CAIRO_INCLUDE_DIRS})
            add_definitions(-DLWEFLATPAK)
            add_library(LightwaveExplorerCuda 
                Source/LightwaveExplorerCore.cu 
                Source/Devices/LightwaveExplorerCoreFP32.cu)
            add_executable(LightwaveExplorer 
                Source/Frontend/LightwaveExplorerFrontendQt.cpp 
                Source/LightwaveExplorerUtilities.cpp 
                Source/Devices/LightwaveExplorerCoreCPU.cpp 
                Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
                Source/Devices/LightwaveExplorerCoreCounter.cpp 
                Source/Devices/DlibLibraryComponents.cpp)
            set_target_properties(LightwaveExplorer PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
            target_compile_options(LightwaveExplorer PRIVATE -O3 ${OpenMP_CXX_FLAGS} -DNOSYCL -DLWEFLATPAK)
            target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
            target_link_libraries(LightwaveExplorer CUDA::cudart_static CUDA::cufft_static CUDA::nvml)
            target_link_libraries(LightwaveExplorer Qt6::Widgets Qt6::DBus)
            target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})
            target_link_libraries(LightwaveExplorer fmt::fmt miniz::miniz)
            target_link_libraries(LightwaveExplorer -lm
                -Wl,--start-group 
                ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
                ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
                ${MKL_ROOT}/lib/intel64/libmkl_core.a 
                -Wl,--end-group
                ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
                )
            target_link_libraries(LightwaveExplorer TBB::tbb)

            add_link_options(-fsycl)
            add_executable(LightwaveExplorerNoCuda 
                Source/Frontend/LightwaveExplorerFrontendQt.cpp 
                Source/LightwaveExplorerUtilities.cpp 
                Source/Devices/LightwaveExplorerCoreCPU.cpp 
                Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
                Source/Devices/LightwaveExplorerCoreCounter.cpp 
                Source/Devices/DlibLibraryComponents.cpp
                Source/Devices/LightwaveExplorerSYCLLinux.cpp 
                Source/Devices/LightwaveExplorerSYCLLinuxFP32.cpp 
                )
            target_compile_options(LightwaveExplorerNoCuda PRIVATE -O3 ${OpenMP_CXX_FLAGS} -fsycl -DNOCUDA -DLWEFLATPAK)
            target_link_libraries(LightwaveExplorerNoCuda Qt6::Widgets Qt6::DBus)
            target_link_libraries(LightwaveExplorerNoCuda ${CAIRO_LIBRARIES})
            target_link_libraries(LightwaveExplorerNoCuda fmt::fmt miniz::miniz)
            target_link_libraries(LightwaveExplorerNoCuda -lm
                ${MKL_ROOT}/lib/intel64/libmkl_sycl.a 
                -Wl,--start-group 
                ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
                ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
                ${MKL_ROOT}/lib/intel64/libmkl_core.a 
                -Wl,--end-group
                ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
                -lsycl -lOpenCL
                )
            target_link_libraries(LightwaveExplorerNoCuda TBB::tbb)
    
            install(TARGETS LightwaveExplorer LightwaveExplorerNoCuda DESTINATION bin)
            install(PROGRAMS ${CMAKE_SOURCE_DIR}/Source/BuildResources/flatpakLauncher.sh DESTINATION bin)
            install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
            install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
            install(FILES Source/BuildResources/io.github.NickKarpowicz.LightwaveExplorer.metainfo.xml DESTINATION share/metainfo)
            install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source/BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
            install(FILES Source/BuildResources/DesktopFileFlatpak DESTINATION share/applications RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)    
    elseif(CLICPU)
        project(LightwaveExplorer LANGUAGES CXX)
        find_package(PkgConfig REQUIRED)
        find_package(OpenMP REQUIRED)
        pkg_check_modules(FFTW REQUIRED fftw3)
        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/dlib)
            message("Using existing dlib clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/davisking/dlib)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/miniz)
            message("Using existing miniz download")
        else()
            execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
            execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        endif()
        include_directories(build/dlib) 
        include_directories(build/miniz)
        include_directories(${FFTW_INCLUDE_DIRS})
        link_directories(${FFTW_LIBRARY_DIRS})
        add_definitions( -DCPUONLY -DUSEFFTW -DNOSYCL)
        add_compile_options(-O3 ${OpenMP_CXX_FLAGS})
        add_executable(LightwaveExplorer 
            Source/LightwaveExplorerCommandLineMain.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer ${OpenMP_CXX_LIBRARIES} miniz)
        target_link_libraries(LightwaveExplorer ${FFTW_LIBRARIES})
    elseif(CLICUDA)
    #. /opt/intel/oneapi/setvars.sh
    #cmake --fresh -DCLICUDA=1 -DCMAKE_CUDA_HOST_COMPILER=clang++ -DCMAKE_CUDA_COMPILER=nvcc -DCMAKE_CUDA_ARCHITECTURES=75 .. -G Ninja
        project(LightwaveExplorer LANGUAGES C CXX CUDA)

        execute_process(COMMAND git clone https://github.com/davisking/dlib)
        execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
        execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        find_package(CUDAToolkit)
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${MKL_ROOT}/include)
        include_directories(build/dlib) 
        include_directories(build/miniz)
        include_directories(build)
        add_library(miniz miniz/miniz.c)
        add_definitions(-DRUNONCUDA)
        add_compile_options(-O3)
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
    #if nothing specified, build CPU version. specify -DUSEMKL to use MKL.
        project(LightwaveExplorer LANGUAGES C CXX)
        find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
        set(CMAKE_AUTOMOC ON)
        find_package(fmt REQUIRED)
        find_package(TBB REQUIRED)
        find_package(PkgConfig REQUIRED)
        find_package(OpenMP REQUIRED)
        pkg_check_modules(CAIRO REQUIRED cairo)
        include_directories(${CAIRO_INCLUDE_DIRS})
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/dlib) 
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/gcem/include)
        include_directories(${CMAKE_CURRENT_BINARY_DIR})
        if(USEMKL)
            find_package(MKL REQUIRED)
            include_directories(${MKL_ROOT}/include/fftw)    
        else()
            pkg_check_modules(FFTW REQUIRED fftw3)
            pkg_check_modules(FFTWF REQUIRED fftw3f)
            include_directories(${FFTW_INCLUDE_DIRS})
            link_directories(${FFTW_LIBRARY_DIRS})         
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/dlib)
            message("Using existing dlib clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/davisking/dlib)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/gcem)
            message("Using existing gcem clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/kthohr/gcem)
        endif()

        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/miniz)
            message("Using existing miniz download")
        else()
            execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
            execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
        endif()
        
        add_library(miniz miniz/miniz.c)
        add_compile_options(-O3 ${OpenMP_CXX_FLAGS} -Wall)
        add_definitions(-DNOCUDA -DNOSYCL -DNDEBUG)
        if(NOT DEFINED USEMKL)
            add_definitions(-DUSEFFTW)
            execute_process(COMMAND echo "NOTE: Building CPU only version since OneAPI not specified. This version uses FFTW for performing Fourier transforms and is licensed under the GPL v3.0")
        endif()

        add_executable(LightwaveExplorer 
            Source/Frontend/LightwaveExplorerFrontendQt.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(LightwaveExplorer Qt6::Widgets Qt6::DBus)
        target_link_libraries(LightwaveExplorer ${OpenMP_CXX_LIBRARIES})
        if(USEMKL)
            target_link_libraries(LightwaveExplorer
                -Wl,--start-group 
                ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
                ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.a 
                ${MKL_ROOT}/lib/intel64/libmkl_core.a 
                -Wl,--end-group)  
        else()
            target_link_libraries(LightwaveExplorer ${FFTW_LIBRARIES} ${FFTWF_LIBRARIES})          
        endif()
        target_link_libraries(LightwaveExplorer TBB::tbb)
        target_link_libraries(LightwaveExplorer fmt::fmt miniz)
        target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})

        install(TARGETS LightwaveExplorer)
        install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
        install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
        install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source//BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
        if(USEMKL)
            configure_file(${CMAKE_SOURCE_DIR}/Source/BuildResources/makeLaunchScript.sh ${CMAKE_BINARY_DIR}/makeLaunchScript.sh COPYONLY)
            add_custom_command(TARGET LightwaveExplorer POST_BUILD COMMAND sh ${CMAKE_BINARY_DIR}/makeLaunchScript.sh)
            install(PROGRAMS ${CMAKE_BINARY_DIR}/LightwaveExplorerLauncher.sh DESTINATION bin)
            install(FILES Source/BuildResources/DesktopFileLauncher 
                DESTINATION share/applications 
                RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
        else()
            install(FILES Source/BuildResources/DesktopFileNoLauncher 
                DESTINATION share/applications 
                RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
        endif()
    endif()
endif()

if(APPLE)
    project(LightwaveExplorer LANGUAGES C CXX)
    find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
    set(CMAKE_AUTOMOC ON)
    find_package(PkgConfig REQUIRED)
    find_package(fmt CONFIG REQUIRED)
    pkg_check_modules(CAIRO REQUIRED cairo)
    include_directories(${CAIRO_INCLUDE_DIRS})
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/dlib) 
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/gcem/include)
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    include_directories(${CMAKE_SOURCE_DIR}/../fftw/usr/local/include)
    link_directories(${CAIRO_LIBRARY_DIRS})
    link_directories(${CMAKE_SOURCE_DIR}/../fftw/usr/local/lib)
    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/dlib)
        message("Using existing dlib clone")
    else()
        execute_process(COMMAND git clone --depth=1 https://github.com/davisking/dlib)
    endif()

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/gcem)
        message("Using existing gcem clone")
    else()
        execute_process(COMMAND git clone --depth=1 https://github.com/kthohr/gcem)
    endif()

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/miniz)
        message("Using existing miniz download")
    else()
        execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip)
        execute_process(COMMAND unzip -o miniz-3.0.2 -d miniz)
    endif()
    add_library(miniz miniz/miniz.c)

    add_definitions(-DCPUONLY -DUSEFFTW -DNOSYCL -DNDEBUG)
    add_compile_options(-O3 -Wall -x objective-c++)
    add_executable(LightwaveExplorer MACOSX_BUNDLE 
        Source/Frontend/LightwaveExplorerFrontendQt.cpp 
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
    target_link_libraries(LightwaveExplorer miniz)
    target_link_libraries(LightwaveExplorer libfftw3_threads.a)
    target_link_libraries(LightwaveExplorer libfftw3f_threads.a)
    target_link_libraries(LightwaveExplorer libfftw3.a)
    target_link_libraries(LightwaveExplorer libfftw3f.a)
    target_link_libraries(LightwaveExplorer Qt6::Widgets Qt6::DBus)
    target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})
endif()

if(WIN32)
    function(copy_after_build TARGET FILE_PATH)
        add_custom_command(TARGET ${TARGET} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${FILE_PATH}
            $<TARGET_FILE_DIR:${TARGET}>
        )
    endfunction()

    if(MAKESYCL)
        #Build this first to make the LightwaveExplorerSYCL .dll and .lib, then build the main project.   
        #cmake --fresh -DMAKESYCL=1 .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" -T "Intel(R) oneAPI DPC++ Compiler 2024"
        #cmake --build . --config Release
        project(LightwaveExplorer LANGUAGES CXX)
        find_package(TBB REQUIRED)
        find_package(MKL REQUIRED)
        find_package(dlib REQUIRED)
        find_package(miniz REQUIRED)
        find_package(OpenMP REQUIRED)
        include_directories(${CMAKE_CURRENT_BINARY_DIR})
        include_directories(${MKL_ROOT}/include/fftw) 
        include_directories(${MKL_ROOT}/include)  
        include_directories(${CMAKE_SOURCE_DIR}/Source)
        add_definitions(-DLIGHTWAVEEXPLORERSYCL_EXPORTS)
        add_compile_options(-fp:precise -O3 -${OpenMP_CXX_FLAGS} -fsycl)
        add_library(LightwaveExplorerSYCL SHARED
            Source/Devices/LightwaveExplorerSYCL.cpp 
            Source/Devices/LightwaveExplorerSYCLFP32.cpp 
            Source/LightwaveExplorerUtilities.cpp)
        target_link_libraries(LightwaveExplorerSYCL miniz::miniz dlib::dlib)
        target_link_libraries(LightwaveExplorerSYCL
            ${MKL_ROOT}/lib/mkl_sycl.lib 
            ${MKL_ROOT}/lib/mkl_intel_ilp64.lib 
            ${MKL_ROOT}/lib/mkl_tbb_thread.lib 
            ${MKL_ROOT}/lib/mkl_core.lib 
            ${MKL_ROOT}/../../compiler/latest/lib/libiomp5md.lib
            sycl7.lib OpenCL.lib)
        target_link_libraries(LightwaveExplorerSYCL TBB::tbb ${OpenMP_CXX_LIBRARIES})
    else()
        #cmake --fresh .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" -DCMAKE_CUDA_ARCHITECTURES="75;86"
        #cmake --build . --config Release
        project(LightwaveExplorer LANGUAGES CUDA CXX)
        find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
        set(CMAKE_AUTOMOC ON)
        find_package(gcem REQUIRED)
        find_package(miniz REQUIRED)
        find_package(PkgConfig REQUIRED)
        find_package(CUDAToolkit REQUIRED)
        find_package(dlib REQUIRED)
        find_package(MKL REQUIRED)
        find_package(OpenMP REQUIRED)
        pkg_check_modules(CAIRO REQUIRED cairo)
        include_directories(${MKL_ROOT}/include/fftw) 
        include_directories(${MKL_ROOT}/include)                   
        include_directories(${CAIRO_INCLUDE_DIRS})
        link_directories(${CAIRO_LIBRARY_DIRS})
        include_directories(${CMAKE_CURRENT_BINARY_DIR})
        add_compile_options(-O3 ${OpenMP_CXX_FLAGS})
        add_executable(LightwaveExplorer WIN32 
            Source/Frontend/LightwaveExplorerFrontendQt.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp
            Source/LightwaveExplorerCore.cu 
            Source/Devices/LightwaveExplorerCoreFP32.cu
            Source/Frontend/LightwaveExplorerIcon.rc)
        target_link_libraries(LightwaveExplorer Qt6::Widgets Qt6::DBus)
        if(USEFFTW)
            target_link_libraries(LightwaveExplorer ${FFTW_LIBRARIES} ${FFTWF_LIBRARIES})
        else()
            target_link_libraries(LightwaveExplorer
                ${MKL_ROOT}/lib/mkl_intel_ilp64.lib 
                ${MKL_ROOT}/lib/mkl_intel_thread.lib
                ${MKL_ROOT}/lib/mkl_core.lib 
                ${MKL_ROOT}/../../compiler/latest/lib/libiomp5md.lib)
        endif()
        target_link_libraries(LightwaveExplorer miniz::miniz)
        target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})
        target_link_libraries(LightwaveExplorer dlib::dlib)
        target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
        target_link_options(LightwaveExplorer PRIVATE "/DELAYLOAD:nvml.dll")
        target_link_libraries(LightwaveExplorer DelayImp.lib ${OpenMP_CXX_LIBRARIES})
        target_link_libraries(LightwaveExplorer ${CMAKE_CURRENT_BINARY_DIR}/Release/LightwaveExplorerSYCL.lib)
        target_link_options(LightwaveExplorer PRIVATE "/DELAYLOAD:LightwaveExplorerSYCL.dll")
        get_target_property(Qt6Core_LOCATION Qt6::Core LOCATION)
        get_filename_component(Qt6_BIN_DIR ${Qt6Core_LOCATION} DIRECTORY)
        get_filename_component(Qt6_BIN_DIR ${Qt6_BIN_DIR} DIRECTORY)
        add_custom_command(TARGET LightwaveExplorer POST_BUILD
            COMMAND ${Qt6_BIN_DIR}/bin/windeployqt $<TARGET_FILE_DIR:LightwaveExplorer>/LightwaveExplorer.exe --release --no-translations --skip-plugin-types imageformats,iconengines,tls,networkinformation,generic
        )
        copy_after_build(LightwaveExplorer ${CMAKE_SOURCE_DIR}/CrystalDatabase.txt)
        copy_after_build(LightwaveExplorer ${CMAKE_SOURCE_DIR}/Source/BuildResources/DefaultValues.ini)
        copy_after_build(LightwaveExplorer ${CMAKE_SOURCE_DIR}/Source/BuildResources/Licenses.txt)
        copy_after_build(LightwaveExplorer ${MKL_ROOT}/../../compiler/latest/bin/libiomp5md.dll)
        #Include CUDA dlls, whatver their current number is
        file(GLOB CUFFT_DLLS "${CUDAToolkit_BIN_DIR}/cufft64*.dll")
        foreach(CUFFT_DLL ${CUFFT_DLLS})
            add_custom_command(TARGET LightwaveExplorer POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CUFFT_DLL}
                $<TARGET_FILE_DIR:LightwaveExplorer>
            )
        endforeach()
        file(GLOB CUDART_DLLS "${CUDAToolkit_BIN_DIR}/cudart64*.dll")
        foreach(CUDART_DLL ${CUDART_DLLS})
            add_custom_command(TARGET LightwaveExplorer POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CUDART_DLL}
                $<TARGET_FILE_DIR:LightwaveExplorer>
            )
        endforeach()

        #hide the dll files so that one can actually see the thing to click
        add_custom_command(TARGET LightwaveExplorer POST_BUILD
            COMMAND cmd /c "attrib +h $<TARGET_FILE_DIR:LightwaveExplorer>\\*.dll"
        )
        install(TARGETS LightwaveExplorer)
    endif()
endif()
