enable_language(CUDA)
find_package(CUDAToolkit REQUIRED)

if(NOT DEFINED USE_SYCL)
    #Flatpak using CUDA, and CPU calculations with FFTs in MKL.
    #Current release version of the Flatpak.
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CAIRO REQUIRED cairo)
    find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
    set(CMAKE_AUTOMOC ON)
    find_package(OpenMP REQUIRED)
    find_package(TBB PATHS /app/opt/intel/oneapi/tbb/latest/lib/cmake REQUIRED)
    find_package(MKL PATHS /app/opt/intel/oneapi/mkl/latest/lib/cmake)
    find_package(miniz)
    include_directories(${CUDA_INCLUDE_DIRS})
    include_directories(${Qt6_INCLUDE_DIRS})
    include_directories(${MKL_ROOT}/include/fftw)
    include_directories(${MKL_ROOT}/include)
    include_directories(${CAIRO_INCLUDE_DIRS})
    add_definitions(-DLWEFLATPAK)
    add_library(LightwaveExplorerCuda STATIC
        Source/LightwaveExplorerCore.cu
        Source/Devices/LightwaveExplorerCoreFP32.cu)
    add_executable(LightwaveExplorer
        Source/Frontend/LightwaveExplorerFrontendQt.cpp
        Source/LightwaveExplorerUtilities.cpp
        Source/Devices/LightwaveExplorerCoreCPU.cpp
        Source/Devices/LightwaveExplorerCoreCPUFP32.cpp
        Source/Devices/LightwaveExplorerCoreCounter.cpp
        Source/Devices/DlibLibraryComponents.cpp
        Source/Frontend/LWEVisualizationsCPU.cpp)
    set_target_properties(${EXECUTABLE_NAME} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    target_compile_options(${EXECUTABLE_NAME} PRIVATE -O3 ${OpenMP_CXX_FLAGS} -DNOSYCL -DLWEFLATPAK)
    target_link_libraries(${EXECUTABLE_NAME} LightwaveExplorerCuda)
    target_link_libraries(${EXECUTABLE_NAME} CUDA::cudart_static CUDA::cufft_static CUDA::nvml)
    target_link_libraries(${EXECUTABLE_NAME} Qt6::Widgets Qt6::DBus)
    target_link_libraries(${EXECUTABLE_NAME} ${CAIRO_LIBRARIES})
    target_link_libraries(${EXECUTABLE_NAME} miniz::miniz)
    target_link_libraries(${EXECUTABLE_NAME} -lm
        -Wl,--start-group
        ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a
        ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a
        ${MKL_ROOT}/lib/intel64/libmkl_core.a
        -Wl,--end-group
        ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
        )
    target_link_libraries(${EXECUTABLE_NAME} TBB::tbb)
    add_executable(LightwaveExplorerNoCuda
        Source/Frontend/LightwaveExplorerFrontendQt.cpp
        Source/LightwaveExplorerUtilities.cpp
        Source/Devices/LightwaveExplorerCoreCPU.cpp
        Source/Devices/LightwaveExplorerCoreCPUFP32.cpp
        Source/Devices/LightwaveExplorerCoreCounter.cpp
        Source/Devices/DlibLibraryComponents.cpp
        Source/Frontend/LWEVisualizationsCPU.cpp)
    target_compile_options(LightwaveExplorerNoCuda PRIVATE ${OpenMP_CXX_FLAGS} -DNOSYCL -DNOCUDA -DLWEFLATPAK)
    target_link_libraries(LightwaveExplorerNoCuda Qt6::Widgets Qt6::DBus)
    target_link_libraries(LightwaveExplorerNoCuda ${CAIRO_LIBRARIES})
    target_link_libraries(LightwaveExplorerNoCuda miniz::miniz)
    target_link_libraries(LightwaveExplorerNoCuda -lm
        -Wl,--start-group
        ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a
        ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a
        ${MKL_ROOT}/lib/intel64/libmkl_core.a
        -Wl,--end-group
        ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a
        )
    target_link_libraries(LightwaveExplorerNoCuda TBB::tbb)

    install(TARGETS ${EXECUTABLE_NAME} LightwaveExplorerNoCuda DESTINATION bin)
    install(PROGRAMS ${CMAKE_SOURCE_DIR}/Source/BuildResources/flatpakLauncher.sh DESTINATION bin)
    install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
    install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
    install(FILES Source/BuildResources/io.github.NickKarpowicz.LightwaveExplorer.metainfo.xml DESTINATION share/metainfo)
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source/BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
    install(FILES Source/BuildResources/DesktopFileFlatpak DESTINATION share/applications RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)

elseif(USE_SYCL)
    #Flatpak using CUDA, and CPU calculations with FFTs in MKL.
    #Current release version of the Flatpak.
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CAIRO REQUIRED cairo)
    find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
    set(CMAKE_AUTOMOC ON)
    find_package(OpenMP REQUIRED)
    find_package(TBB PATHS /app/opt/intel/oneapi/tbb/latest/lib/cmake REQUIRED)
    find_package(MKL PATHS /app/opt/intel/oneapi/mkl/latest/lib/cmake REQUIRED)
    find_package(oneDPL PATHS /app/opt/intel/oneapi/dpl/latest/lib/cmake REQUIRED)
    find_package(miniz)
    include_directories(${CUDA_INCLUDE_DIRS})
    include_directories(${Qt6_INCLUDE_DIRS})
    include_directories(${MKL_ROOT}/include/fftw)
    include_directories(${MKL_ROOT}/include)
    include_directories(${CAIRO_INCLUDE_DIRS})
    add_definitions(-DLWEFLATPAK)
    add_library(LightwaveExplorerCuda STATIC
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
    target_link_libraries(LightwaveExplorer miniz::miniz)
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
    target_link_options(LightwaveExplorerNoCuda PRIVATE -fsycl)
    target_link_libraries(LightwaveExplorerNoCuda Qt6::Widgets Qt6::DBus)
    target_link_libraries(LightwaveExplorerNoCuda ${CAIRO_LIBRARIES})
    target_link_libraries(LightwaveExplorerNoCuda miniz::miniz)
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
endif()
