enable_language(CUDA)
find_package(CUDAToolkit REQUIRED)


#Flatpak using CUDA
#Current release version of the Flatpak.
find_package(PkgConfig REQUIRED)
pkg_check_modules(CAIRO REQUIRED cairo)
find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
set(CMAKE_AUTOMOC ON)
find_package(miniz)
find_package(TBB REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${Qt6_INCLUDE_DIRS})
include_directories(${MKL_ROOT}/include)
include_directories(${CAIRO_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(/app/include/miniz)
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
target_compile_options(${EXECUTABLE_NAME} PRIVATE -O3 -DNOSYCL -DLWEFLATPAK)
target_link_libraries(${EXECUTABLE_NAME} LightwaveExplorerCuda)
target_link_libraries(${EXECUTABLE_NAME} CUDA::cudart_static CUDA::cufft_static)
target_link_libraries(${EXECUTABLE_NAME} Qt6::Widgets Qt6::DBus)
target_link_libraries(${EXECUTABLE_NAME} ${CAIRO_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} miniz::miniz)
target_link_libraries(${EXECUTABLE_NAME} TBB::tbb)
add_executable(LightwaveExplorerNoCuda
        Source/Frontend/LightwaveExplorerFrontendQt.cpp
        Source/LightwaveExplorerUtilities.cpp
        Source/Devices/LightwaveExplorerCoreCPU.cpp
        Source/Devices/LightwaveExplorerCoreCPUFP32.cpp
        Source/Devices/LightwaveExplorerCoreCounter.cpp
        Source/Devices/DlibLibraryComponents.cpp
        Source/Frontend/LWEVisualizationsCPU.cpp)
target_compile_options(LightwaveExplorerNoCuda PRIVATE -DNOSYCL -DNOCUDA -DLWEFLATPAK)
target_link_libraries(LightwaveExplorerNoCuda Qt6::Widgets Qt6::DBus)
target_link_libraries(LightwaveExplorerNoCuda ${CAIRO_LIBRARIES})
target_link_libraries(LightwaveExplorerNoCuda miniz::miniz)
target_link_libraries(LightwaveExplorerNoCuda TBB::tbb)
install(TARGETS ${EXECUTABLE_NAME} LightwaveExplorerNoCuda DESTINATION bin)
install(PROGRAMS ${CMAKE_SOURCE_DIR}/Source/BuildResources/flatpakLauncher.sh DESTINATION bin)
install(FILES CrystalDatabase.txt DESTINATION share/LightwaveExplorer)
install(FILES Source/BuildResources/DefaultValues.ini DESTINATION share/LightwaveExplorer)
install(FILES Source/BuildResources/io.github.NickKarpowicz.LightwaveExplorer.metainfo.xml DESTINATION share/metainfo)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/Source/BuildResources/icons DESTINATION share FILES_MATCHING PATTERN "*")
install(FILES Source/BuildResources/DesktopFileFlatpak DESTINATION share/applications RENAME io.github.NickKarpowicz.LightwaveExplorer.desktop)
