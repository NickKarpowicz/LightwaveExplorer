cmake_minimum_required(VERSION 3.16)
project(LightwaveExplorer)
execute_process(COMMAND brew --prefix OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE BREW_PREFIX)
list(APPEND CMAKE_MODULE_PATH "${BREW_PREFIX}/share/cmake/Modules")
execute_process(COMMAND brew --prefix llvm OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE LLVM_PREFIX)
set(CMAKE_CXX_COMPILER "${LLVM_PREFIX}/bin/clang++")
set(CMAKE_C_COMPILER "${LLVM_PREFIX}/bin/clang")

find_package(PkgConfig REQUIRED)
pkg_check_modules(GTK4 REQUIRED gtk4)
pkg_check_modules(FFTW REQUIRED fftw3)
pkg_check_modules(FMT REQUIRED fmt)

include_directories(${FFTW_INCLUDE_DIRS})
link_directories(${FFTW_LIBRARY_DIRS})

include_directories(${GTK4_INCLUDE_DIRS})
link_directories(${GTK4_LIBRARY_DIRS})
include_directories(${FMT_INCLUDE_DIRS})
link_directories(${FMT_LIBRARY_DIRS})

if(APPLE)
    find_library(OpenMP_LIBRARY NAMES libomp.a HINTS ${BREW_PREFIX}/opt/libomp/lib/)
else()
    find_package(OpenMP)
endif()

include_directories(../dlib)
add_compile_options(-std=c++20 -O3 -fopenmp -D CPUONLY -w)
add_executable(LightwaveExplorer MACOSX_BUNDLE LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.mm LightwaveExplorerUtilities.cpp LightwaveExplorerCoreCPU.cpp DlibLibraryComponents/DlibLibraryComponents.cpp)
set_target_properties(LightwaveExplorer  PROPERTIES
    BUNDLE True
    MACOSX_BUNDLE_BUNDLE_NAME LightwaveExplorer
    MACOSX_BUNDLE_VERSION "0.1"
    MACOSX_BUNDLE_VERSION_STRING "0.1"
    MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/MacResources/macplistbase.plist
    MACOSX_BUNDLE_ICON_FILE ${CMAKE_SOURCE_DIR}/MacResources/AppIcon.icns
)
target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES} ${FFTW_LIBRARIES} ${FMT_LIBRARIES})
if(APPLE)
    target_link_libraries(LightwaveExplorer ${OpenMP_LIBRARY} "-framework Cocoa")
endif()