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
    include_directories(${CMAKE_CURRENT_BINARY_DIR})

    add_library(miniz STATIC miniz/miniz.c)
endmacro()

macro(linux_install_launchscript)
    configure_file(${CMAKE_SOURCE_DIR}/Source/BuildResources/makeLaunchScript.sh ${CMAKE_BINARY_DIR}/makeLaunchScript.sh COPYONLY)
    add_custom_command(TARGET ${EXECUTABLE_NAME} POST_BUILD COMMAND sh ${CMAKE_BINARY_DIR}/makeLaunchScript.sh)
    install(TARGETS ${EXECUTABLE_NAME})
    install(PROGRAMS ${CMAKE_BINARY_DIR}/LightwaveExplorerLauncher.sh DESTINATION bin)
    install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
    install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source//BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
    install(FILES Source/BuildResources/DesktopFileLauncher
        DESTINATION share/applications
        RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
endmacro()

macro(linux_install_nolaunchscript)
    install(TARGETS ${EXECUTABLE_NAME})
    install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
    install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source//BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
    install(FILES Source/BuildResources/DesktopFileNoLauncher
            DESTINATION share/applications
            RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
endmacro()

macro(linux_ui_packages)
    find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
    set(CMAKE_AUTOMOC ON)
    if(NOT DEFINED BACKEND_ROCM)
        find_package(OpenMP)
    endif()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CAIRO REQUIRED cairo)
    include_directories(${CAIRO_INCLUDE_DIRS})
endmacro()

macro(add_oneapi_interfaces)
    set(SYCL_TARGETS "")
    if(BACKEND_ROCM)
        if(NOT BACKEND_ROCM MATCHES "^gfx[0-9]+")
            message(FATAL_ERROR "BACKEND_ROCM is set to an invalid value: ${BACKEND_ROCM}. It should be an AMD GPU architecture, for example, gfx906 or gfx1030. Note that only one architecture can be built in a single binary.")
        endif()
        set(ENABLE_ROCFFT_BACKEND True)
        set(ENABLE_MKLCPU_BACKEND False)
        set(ENABLE_MKLGPU_BACKEND False)
        list(APPEND SYCL_TARGETS "amdgcn-amd-amdhsa")
    endif()
    if(BACKEND_CUDA)
        set(ENABLE_CUFFT_BACKEND True)
        set(ENABLE_MKLCPU_BACKEND False)
        set(ENABLE_MKLGPU_BACKEND False)
        list(APPEND SYCL_TARGETS "nvptx64-nvidia-cuda")
    endif()
    if(NOT SYCL_TARGETS)
        set(BACKEND_INTEL True)
    endif()
    if(BACKEND_INTEL)
        set(ENABLE_MKLCPU_BACKEND True)
        set(ENABLE_MKLGPU_BACKEND True)
        list(APPEND SYCL_TARGETS "spir64")
    endif()


    list(JOIN SYCL_TARGETS "," SYCL_TARGETS_STRING)
    include(FetchContent)
    set(BUILD_FUNCTIONAL_TESTS False)
    set(BUILD_EXAMPLES False)
    FetchContent_Declare(
            onemkl_interface_library
            GIT_REPOSITORY https://github.com/uxlfoundation/oneMath.git
            GIT_TAG v0.9
    )
    FetchContent_MakeAvailable(onemkl_interface_library)

    message(STATUS "SYCL Targets: ${SYCL_TARGETS_STRING}")
endmacro()

macro(set_sycl_compile_flags)
    add_compile_options(
        -ffp-model=precise
        -O3
        -fsycl
        -fsycl-targets=${SYCL_TARGETS_STRING}
        -w)
    if(OpenMP_FOUND)
        add_compile_options(${OpenMP_CXX_FLAGS})
    endif()
    if(BACKEND_ROCM)
        add_compile_options(
            -Xsycl-target-backend=amdgcn-amd-amdhsa --offload-arch=${BACKEND_ROCM}
        )
        if(ROCM_LIB_PATH)
            add_compile_options(--rocm-device-lib-path=${ROCM_LIB_PATH})
        endif()
    endif()

    add_link_options(
        -fsycl
        -fsycl-targets=${SYCL_TARGETS_STRING})

    if(BACKEND_ROCM)
        add_link_options(
            -Xsycl-target-backend=amdgcn-amd-amdhsa --offload-arch=${BACKEND_ROCM}
            --rocm-device-lib-path=${ROCM_LIB_PATH}
        )
        if(ROCM_LIB_PATH)
            add_link_options(--rocm-device-lib-path=${ROCM_LIB_PATH})
        endif()
    endif()
endmacro()

macro(resolve_fft_library)
    find_package(MKL QUIET)
    if(MKL_FOUND)
        find_package(TBB REQUIRED)
        include_directories(${MKL_ROOT}/include/fftw)
        set(USING_MKL True)
        message("Using MKL for FFTs.")
    else()
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(FFTW REQUIRED fftw3)
        pkg_check_modules(FFTWF REQUIRED fftw3f)
        include_directories(${FFTW_INCLUDE_DIRS})
        link_directories(${FFTW_LIBRARY_DIRS})
        set(USING_FFTW True)
        add_definitions(-DUSEFFTW)
        message("Using FFTW for FFTs.")
    endif()
endmacro()

macro(link_fft_library FFT_TARGET)
    if(USING_MKL)
        target_link_libraries(${FFT_TARGET}
            -Wl,--start-group
            ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a
            ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a
            ${MKL_ROOT}/lib/intel64/libmkl_core.a
            -Wl,--end-group)
        target_link_libraries(${FFT_TARGET} TBB::tbb)
    elseif(USING_FFTW)
        target_link_libraries(${FFT_TARGET} ${FFTW_LIBRARIES} ${FFTWF_LIBRARIES})
    else()
        message(FATAL_ERROR "Could not link because neither FFT library was found.")
    endif()
endmacro()
