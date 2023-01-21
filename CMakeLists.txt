cmake_minimum_required(VERSION 3.16)
project(LightwaveExplorer)

set(BREW_PREFIX, "/usr/local")
list(APPEND CMAKE_MODULE_PATH "${BREW_PREFIX}/share/cmake/Modules")
set(LLVM_PREFIX "/usr/local/opt/llvm")
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

include_directories(../dlib)
add_compile_options(-std=c++20 -O3 -D CPUONLY -w)
execute_process(COMMAND sed -i.bak "s/fftw3_mkl.h/fftw3.h/g" ../LightwaveExplorerUtilities.h)
execute_process(COMMAND sed -i.bak "s/fftw3_mkl.h/fftw3.h/g" ../LWEActiveDeviceCPU.h)
execute_process(COMMAND sed -i.bak "s&<format>&<fmt/format.h>&g ; s&std::format&fmt::format&g ; s&std::vformat&fmt::vformat&g ; s&std::make_format_args&fmt::make_format_args&g" ../LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h)
execute_process(COMMAND sed -i.bak "s&<format>&<fmt/format.h>&g ; s&std::format&fmt::format&g ; s&std::vformat&fmt::vformat&g ; s&std::make_format_args&fmt::make_format_args&g" ../LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp)
cp ../LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp ../LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.mm
add_executable(LightwaveExplorer MACOSX_BUNDLE LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.mm LightwaveExplorerUtilities.cpp LightwaveExplorerCoreCPU.cpp DlibLibraryComponents/DlibLibraryComponents.cpp)
set_target_properties(LightwaveExplorer  PROPERTIES
    BUNDLE True
    MACOSX_BUNDLE_BUNDLE_NAME LightwaveExplorer
    MACOSX_BUNDLE_VERSION "0.1"
    MACOSX_BUNDLE_VERSION_STRING "0.1"
    MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/MacResources/macplistbase.plist
)
target_link_libraries(LightwaveExplorer ${GTK4_LIBRARIES} ${FFTW_LIBRARIES} ${FMT_LIBRARIES})