enable_language(C)
find_package(Qt6 COMPONENTS Widgets DBus REQUIRED)
set(CMAKE_AUTOMOC ON)
find_package(PkgConfig REQUIRED)
pkg_check_modules(CAIRO REQUIRED cairo)
include_directories(${CAIRO_INCLUDE_DIRS})

link_directories(${CAIRO_LIBRARY_DIRS})

conditionally_fetch_dependencies()

add_definitions(-DCPUONLY -DNOSYCL -DNDEBUG)
add_compile_options(-O3 -Wall -mmacosx-version-min=13.3 -march=native)
add_executable(LightwaveExplorer MACOSX_BUNDLE
    Source/Frontend/LightwaveExplorerFrontendQt.cpp
    Source/LightwaveExplorerUtilities.cpp
    Source/Devices/LightwaveExplorerCoreCPUFP32.cpp
    Source/Devices/LightwaveExplorerCoreCPU.cpp
    Source/Devices/LightwaveExplorerCoreCounter.cpp
    Source/Devices/DlibLibraryComponents.cpp
    Source/Frontend/LWEVisualizationsCPU.cpp)
set_target_properties(LightwaveExplorer  PROPERTIES
    BUNDLE True
    MACOSX_BUNDLE_BUNDLE_NAME LightwaveExplorer
    MACOSX_BUNDLE_VERSION "0.1"
    MACOSX_BUNDLE_VERSION_STRING "0.1"
    MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/Source/BuildResources/macplistbase.plist
    MACOSX_BUNDLE_ICON_FILE ${CMAKE_SOURCE_DIR}/Source/BuildResources/AppIcon.icns
)

target_link_libraries(LightwaveExplorer miniz)
target_link_libraries(LightwaveExplorer Qt6::Widgets Qt6::DBus)
target_link_libraries(LightwaveExplorer ${CAIRO_LIBRARIES})
