enable_language(C)
find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
set(CMAKE_AUTOMOC ON)
find_package(PkgConfig REQUIRED)
find_package(fmt CONFIG REQUIRED)
pkg_check_modules(CAIRO REQUIRED cairo)
include_directories(${CAIRO_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/../fftw/usr/local/include)
link_directories(${CAIRO_LIBRARY_DIRS})
link_directories(${CMAKE_SOURCE_DIR}/../fftw/usr/local/lib)
conditionally_fetch_dependencies()

add_definitions(-DCPUONLY -DUSEFFTW -DNOSYCL -DNDEBUG)
add_compile_options(-O3 -Wall)
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