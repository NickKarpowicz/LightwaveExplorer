function(copy_after_build TARGET FILE_PATH)
    add_custom_command(TARGET ${TARGET} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${FILE_PATH}
        $<TARGET_FILE_DIR:${TARGET}>
    )
endfunction()

macro(conditionally_fetch_dependencies)
    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/dlib)
        message("Using existing dlib clone")
    else()
        execute_process(COMMAND wget https://github.com/davisking/dlib/archive/refs/tags/v19.24.6.zip)
        execute_process(COMMAND unzip -q -o v19.24.6.zip)
        execute_process(COMMAND mv dlib-19.24.6 dlib)
    endif()

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/gcem)
        message("Using existing gcem clone")
    else()
        execute_process(COMMAND wget https://github.com/kthohr/gcem/archive/refs/tags/v1.18.0.zip)
        execute_process(COMMAND unzip -q -o v1.18.0.zip)
        execute_process(COMMAND mv gcem-1.18.0 gcem)
    endif()

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/miniz)
        message("Using existing miniz download")
    else()
        execute_process(COMMAND wget https://github.com/richgel999/miniz/releases/download/3.1.0/miniz-3.1.0.zip)
        execute_process(COMMAND unzip -q -o miniz-3.1.0 -d miniz)
    endif()

    include_directories(${CMAKE_CURRENT_BINARY_DIR}/dlib)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/gcem/include)
    include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/miniz)

endmacro()

if(MAKEDEPENDENCIES)
    project(LightwaveExplorer LANGUAGES C CXX)
    conditionally_fetch_dependencies()
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    add_library(LightwaveExplorerDependencies STATIC
        Source/Devices/DlibLibraryComponents.cpp)
    add_library(miniz STATIC miniz/miniz.c)
elseif(MAKESYCL)
    #Build this first to make the LightwaveExplorerSYCL .dll and .lib, then build the main project.
    #cmake --fresh -DMAKESYCL=1 .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" -T "Intel(R) oneAPI DPC++ Compiler 2024"
    #cmake --build . --config Release
    project(LightwaveExplorer LANGUAGES CXX)
    conditionally_fetch_dependencies()
    find_package(TBB REQUIRED)
    find_package(MKL REQUIRED)

    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    include_directories(${MKL_ROOT}/include/fftw)
    include_directories(${MKL_ROOT}/include)
    include_directories(${CMAKE_SOURCE_DIR}/Source)
    add_definitions(-DLIGHTWAVEEXPLORERSYCL_EXPORTS)
    add_compile_options(-fp:precise -O3 -fsycl)
    add_library(LightwaveExplorerSYCL SHARED
        Source/Devices/LightwaveExplorerSYCL.cpp
        Source/Devices/LightwaveExplorerSYCLFP32.cpp
        Source/LightwaveExplorerUtilities.cpp)
    target_link_libraries(LightwaveExplorerSYCL ${CMAKE_CURRENT_BINARY_DIR}/Release/LightwaveExplorerDependencies.lib)
    target_link_libraries(LightwaveExplorerSYCL ${CMAKE_CURRENT_BINARY_DIR}/Release/miniz.lib)
    target_link_libraries(LightwaveExplorerSYCL
        ${MKL_ROOT}/lib/mkl_sycl.lib
        ${MKL_ROOT}/lib/mkl_intel_ilp64.lib
        ${MKL_ROOT}/lib/mkl_tbb_thread.lib
        ${MKL_ROOT}/lib/mkl_core.lib
        ${MKL_ROOT}/../../compiler/latest/lib/libiomp5md.lib
        sycl8.lib OpenCL.lib)
    target_link_libraries(LightwaveExplorerSYCL TBB::tbb)
else()
    #cmake --fresh .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" -DCMAKE_CUDA_ARCHITECTURES="75;86"
    #cmake --build . --config Release
    project(LightwaveExplorer LANGUAGES CUDA CXX)
    conditionally_fetch_dependencies()
    find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
    set(CMAKE_AUTOMOC ON)
    find_package(PkgConfig REQUIRED)
    find_package(CUDAToolkit REQUIRED)
    find_package(MKL REQUIRED)
    find_package(OpenMP)
    pkg_check_modules(CAIRO REQUIRED cairo)
    include_directories(${MKL_ROOT}/include/fftw)
    include_directories(${MKL_ROOT}/include)
    include_directories(${CAIRO_INCLUDE_DIRS})
    link_directories(${CAIRO_LIBRARY_DIRS})
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    add_executable(LightwaveExplorer WIN32
        Source/Frontend/LightwaveExplorerFrontendQt.cpp
        Source/LightwaveExplorerUtilities.cpp
        Source/Devices/LightwaveExplorerCoreCPU.cpp
        Source/Devices/LightwaveExplorerCoreCPUFP32.cpp
        Source/Devices/LightwaveExplorerCoreCounter.cpp
        Source/LightwaveExplorerCore.cu
        Source/Devices/LightwaveExplorerCoreFP32.cu
        Source/Frontend/LWEVisualizationsCPU.cpp
        Source/Frontend/LightwaveExplorerIcon.rc)
    target_compile_options(LightwaveExplorer
        PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>
    )
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
    target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})
    target_link_libraries(LightwaveExplorer ${CMAKE_CURRENT_BINARY_DIR}/Release/LightwaveExplorerDependencies.lib)
    target_link_libraries(LightwaveExplorer ${CMAKE_CURRENT_BINARY_DIR}/Release/miniz.lib)
    target_link_libraries(LightwaveExplorer CUDA::cudart CUDA::cufft)
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
    file(GLOB CUFFT_DLLS "${CUDAToolkit_BIN_DIR}/x64/cufft64*.dll")
    foreach(CUFFT_DLL ${CUFFT_DLLS})
        add_custom_command(TARGET LightwaveExplorer POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CUFFT_DLL}
            $<TARGET_FILE_DIR:LightwaveExplorer>
        )
    endforeach()
    file(GLOB CUDART_DLLS "${CUDAToolkit_BIN_DIR}/x64/cudart64*.dll")
    foreach(CUDART_DLL ${CUDART_DLLS})
        add_custom_command(TARGET LightwaveExplorer POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CUDART_DLL}
            $<TARGET_FILE_DIR:LightwaveExplorer>
        )
    endforeach()

    #hide the dll files so that one can actually see the thing to click
    add_custom_command(TARGET LightwaveExplorer POST_BUILD
        COMMAND attrib +h $<TARGET_FILE_DIR:LightwaveExplorer>\\*.dll
    )
    install(TARGETS LightwaveExplorer)
endif()
