enable_language(C)
if(USE_CUDA)
    enable_language(CUDA)      
    find_package(CUDAToolkit REQUIRED)
endif()

if(USE_CUDA)
    conditionally_fetch_dependencies()
    resolve_fft_library()
    add_definitions(-DRUNONCUDA)
    add_compile_options(-O3)
    add_executable(${EXECUTABLE_NAME} 
    Source/LightwaveExplorerCore.cu
    Source/LightwaveExplorerCommandLineMain.cu
    Source/LightwaveExplorerUtilities.cpp 
    Source/Devices/DlibLibraryComponents.cpp)
    target_link_libraries(${EXECUTABLE_NAME} miniz)
    target_link_libraries(${EXECUTABLE_NAME} -lm)
    link_fft_library(${EXECUTABLE_NAME})
    target_link_libraries(${EXECUTABLE_NAME} CUDA::cudart CUDA::cufft CUDA::nvml)
elseif(USE_SYCL)
    find_package(MKL REQUIRED)
    conditionally_fetch_dependencies()
    add_oneapi_interfaces()
    resolve_fft_library()
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

    target_link_libraries(${EXECUTABLE_NAME} onemkl)
    target_link_libraries(${EXECUTABLE_NAME} -lsycl)
    target_link_libraries(${EXECUTABLE_NAME} -lm miniz)
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
