enable_language(C)
if(USE_CUDA)
    enable_language(CUDA)      
    find_package(CUDAToolkit REQUIRED)
endif()

if(USE_CUDA AND USE_SYCL)   
    find_package(TBB REQUIRED)
    linux_ui_packages()
    conditionally_fetch_dependencies()
    resolve_fft_library()
    add_oneapi_interfaces()

    add_library(LightwaveExplorerCuda STATIC
        Source/LightwaveExplorerCore.cu 
        Source/Devices/LightwaveExplorerCoreFP32.cu)

    set_sycl_compile_flags()

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
    link_fft_library(LightwaveExplorer)
    target_link_libraries(LightwaveExplorer -lsycl -lOpenCL)
    target_link_libraries(LightwaveExplorer TBB::tbb onemkl)
    target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
    target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
    if(OpenMP_FOUND)
        target_link_libraries(LightwaveExplorer ${OpenMP_CXX_LIBRARIES})
    endif()
    target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})

    linux_install_launchscript()
    install(TARGETS onemkl)
elseif(USE_CUDA)   
    find_package(TBB REQUIRED)
    conditionally_fetch_dependencies()
    linux_ui_packages()
    resolve_fft_library()
    add_library(LightwaveExplorerCuda STATIC
        Source/LightwaveExplorerCore.cu 
        Source/Devices/LightwaveExplorerCoreFP32.cu)
    add_definitions(-DNOSYCL)
    add_compile_options(-O3 -Wpedantic)
    if(OpenMP_FOUND)
        add_compile_options(${OpenMP_CXX_FLAGS})
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
    if(OpenMP_FOUND)
        target_link_libraries(LightwaveExplorer ${OpenMP_CXX_LIBRARIES})
    endif()
    link_fft_library(LightwaveExplorer)
    target_link_libraries(LightwaveExplorer LightwaveExplorerCuda)
    target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft CUDA::nvml)
    target_link_libraries(LightwaveExplorer TBB::tbb)
    target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})
    linux_install_nolaunchscript()

elseif(AMDLIBRARY)   
    find_package(TBB REQUIRED)
    add_compile_options(-fPIC)
    resolve_fft_library()
    conditionally_fetch_dependencies()
    add_oneapi_interfaces()
    add_definitions(-DNOCUDA)
    set_sycl_compile_flags()
    set(UTILITY_LIBRARY_NAME "${EXECUTABLE_NAME}_utilities")

    add_library(${UTILITY_LIBRARY_NAME} SHARED 
        Source/LightwaveExplorerUtilities.cpp 
        Source/Devices/DlibLibraryComponents.cpp)
    target_link_libraries(${UTILITY_LIBRARY_NAME} miniz -lm)
    link_fft_library(${UTILITY_LIBRARY_NAME})

    set(LIBRARY_NAME "${EXECUTABLE_NAME}_amd")
    add_library(${LIBRARY_NAME} SHARED 
        Source/Devices/LightwaveExplorerSYCLLinux.cpp 
        Source/Devices/LightwaveExplorerSYCLLinuxFP32.cpp )
    target_link_libraries(${LIBRARY_NAME} ${UTILITY_LIBRARY_NAME})
    target_link_libraries(${LIBRARY_NAME}  onemkl)
    target_link_libraries(${LIBRARY_NAME}  -lsycl -lOpenCL)
    target_link_libraries(${LIBRARY_NAME}  -lm)
    target_link_libraries(${LIBRARY_NAME}  TBB::tbb)

    if(BUILD_EXECUTABLE)
        linux_ui_packages()
        add_executable(${EXECUTABLE_NAME} 
            Source/Frontend/LightwaveExplorerFrontendQt.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
            Source/Devices/LightwaveExplorerCoreCounter.cpp)
        target_link_libraries(${EXECUTABLE_NAME}  ${LIBRARY_NAME} onemkl)
        target_link_libraries(${EXECUTABLE_NAME}  -lsycl -lOpenCL)
        target_link_libraries(${EXECUTABLE_NAME}  Qt6::Widgets Qt6::DBus)
        target_link_libraries(${EXECUTABLE_NAME}  fmt::fmt -lm miniz)
        link_fft_library(${EXECUTABLE_NAME})
        target_link_libraries(${EXECUTABLE_NAME}  TBB::tbb)
        target_link_libraries(${EXECUTABLE_NAME}  ${CAIRO_LIBRARIES})
        if(OpenMP_FOUND)
            target_link_libraries(LightwaveExplorer ${OpenMP_CXX_LIBRARIES})
        endif()
    endif()

elseif(USE_SYCL)   
    find_package(TBB REQUIRED)

    resolve_fft_library()
    conditionally_fetch_dependencies()
    linux_ui_packages()
    add_oneapi_interfaces()
    add_definitions(-DNOCUDA)
    set_sycl_compile_flags()

    add_executable(${EXECUTABLE_NAME} 
        Source/Devices/LightwaveExplorerSYCLLinux.cpp 
        Source/Devices/LightwaveExplorerSYCLLinuxFP32.cpp 
        Source/Frontend/LightwaveExplorerFrontendQt.cpp 
        Source/LightwaveExplorerUtilities.cpp 
        Source/Devices/LightwaveExplorerCoreCPU.cpp 
        Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
        Source/Devices/LightwaveExplorerCoreCounter.cpp 
        Source/Devices/DlibLibraryComponents.cpp)
    target_link_libraries(${EXECUTABLE_NAME}  onemkl)
    target_link_libraries(${EXECUTABLE_NAME}  -lsycl -lOpenCL)
    target_link_libraries(${EXECUTABLE_NAME}  Qt6::Widgets Qt6::DBus)
    target_link_libraries(${EXECUTABLE_NAME}  fmt::fmt -lm miniz)
    link_fft_library(${EXECUTABLE_NAME})
    target_link_libraries(${EXECUTABLE_NAME}  TBB::tbb)
    target_link_libraries(${EXECUTABLE_NAME}  ${CAIRO_LIBRARIES})
    if(OpenMP_FOUND)
        target_link_libraries(LightwaveExplorer ${OpenMP_CXX_LIBRARIES})
    endif()
    
    linux_install_launchscript()
    install(TARGETS onemkl)
   
elseif(CLI)
    if(USE_CUDA)
        #. /opt/intel/oneapi/setvars.sh
        #cmake --fresh -DCLICUDA=1 -DCMAKE_CUDA_HOST_COMPILER=clang++ -DCMAKE_CUDA_COMPILER=nvcc -DCMAKE_CUDA_ARCHITECTURES=75 .. -G Ninja
        conditionally_fetch_dependencies()
        include_directories(${MKL_ROOT}/include/fftw)
        include_directories(${MKL_ROOT}/include)
        add_definitions(-DRUNONCUDA)
        add_compile_options(-O3)
        add_executable(${EXECUTABLE_NAME} 
        Source/LightwaveExplorerCore.cu
        Source/LightwaveExplorerCommandLineMain.cu
        Source/LightwaveExplorerUtilities.cpp 
        Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(${EXECUTABLE_NAME} miniz)
        target_link_libraries(${EXECUTABLE_NAME} -lm)
        target_link_libraries(${EXECUTABLE_NAME}
            -Wl,--start-group 
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
            ${MKL_ROOT}/lib/intel64/libmkl_sequential.a 
            ${MKL_ROOT}/lib/intel64/libmkl_core.a 
            -Wl,--end-group
            )
        target_link_libraries(${EXECUTABLE_NAME} CUDA::cudart CUDA::cufft CUDA::nvml)
    elseif(USE_SYCL)
        find_package(TBB REQUIRED)
        find_package(MKL REQUIRED)
        conditionally_fetch_dependencies()
        add_oneapi_interfaces()
        include_directories(${MKL_ROOT}/include/fftw)
        add_definitions(-DNOCUDA)
        add_definitions(-DRUNONSYCL)
        set_sycl_compile_flags()

        if(FP32)
            add_definitions(-DLWEFLOATINGPOINT=32)
            add_executable(${EXECUTABLE_NAME} 
                Source/Devices/LightwaveExplorerSYCLLinuxFP32.cpp 
                Source/LightwaveExplorerCommandLineMain.cpp 
                Source/LightwaveExplorerUtilities.cpp 
                Source/Devices/DlibLibraryComponents.cpp)
        else()
            add_executable(${EXECUTABLE_NAME} 
                Source/Devices/LightwaveExplorerSYCLLinux.cpp 
                Source/LightwaveExplorerCommandLineMain.cpp 
                Source/LightwaveExplorerUtilities.cpp 
                Source/Devices/DlibLibraryComponents.cpp)
        endif()

        target_link_libraries(${EXECUTABLE_NAME} MKL::MKL_SYCL onemkl)
        target_link_libraries(${EXECUTABLE_NAME} -lsycl -lOpenCL)
        target_link_libraries(${EXECUTABLE_NAME} -lm miniz)
        target_link_libraries(${EXECUTABLE_NAME} TBB::tbb)
        target_link_libraries(${EXECUTABLE_NAME} ${OpenMP_CXX_LIBRARIES})
    else()
        find_package(OpenMP REQUIRED)
        resolve_fft_library()
        conditionally_fetch_dependencies()
        add_definitions( -DCPUONLY -DNOSYCL)
        add_compile_options(-O3 ${OpenMP_CXX_FLAGS})
        add_executable(${EXECUTABLE_NAME} 
            Source/LightwaveExplorerCommandLineMain.cpp 
            Source/LightwaveExplorerUtilities.cpp 
            Source/Devices/LightwaveExplorerCoreCPU.cpp 
            Source/Devices/DlibLibraryComponents.cpp)
        target_link_libraries(${EXECUTABLE_NAME} ${OpenMP_CXX_LIBRARIES} miniz)
        link_fft_library(${EXECUTABLE_NAME})
    endif()
else()
#if nothing specified, build CPU version. specify -DUSE_MKL to use MKL.
    find_package(TBB REQUIRED)
    find_package(OpenMP REQUIRED)
    conditionally_fetch_dependencies()
    linux_ui_packages()
    resolve_fft_library()
    add_compile_options(-O3 ${OpenMP_CXX_FLAGS} -Wall)
    add_definitions(-DNOCUDA -DNOSYCL -DNDEBUG)
    add_executable(${EXECUTABLE_NAME} 
        Source/Frontend/LightwaveExplorerFrontendQt.cpp 
        Source/LightwaveExplorerUtilities.cpp 
        Source/Devices/LightwaveExplorerCoreCPU.cpp 
        Source/Devices/LightwaveExplorerCoreCPUFP32.cpp 
        Source/Devices/LightwaveExplorerCoreCounter.cpp 
        Source/Devices/DlibLibraryComponents.cpp)
    target_link_libraries(${EXECUTABLE_NAME} Qt6::Widgets Qt6::DBus)
    target_link_libraries(${EXECUTABLE_NAME} ${OpenMP_CXX_LIBRARIES})
    link_fft_library(${EXECUTABLE_NAME})
    target_link_libraries(${EXECUTABLE_NAME} TBB::tbb)
    target_link_libraries(${EXECUTABLE_NAME} fmt::fmt miniz)
    target_link_libraries(${EXECUTABLE_NAME} ${CAIRO_LIBRARIES})
    linux_install_nolaunchscript()
endif()