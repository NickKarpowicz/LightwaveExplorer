cmake_minimum_required(VERSION 3.23.0)
project(LightwaveExplorer LANGUAGES C CXX)
if(LWE_VULKAN)
  set(CMAKE_CXX_STANDARD 23)
else()
  set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


if(NOT DEFINED EXECUTABLE_NAME)
    set(EXECUTABLE_NAME "LightwaveExplorer")
endif()
if(LWE_VULKAN)
    find_package(Vulkan REQUIRED)
    find_package(glslang CONFIG REQUIRED)
    find_package(vpp CONFIG REQUIRED)
    include(Source/BuildResources/cmake/linux_macros.cmake)
    conditionally_fetch_dependencies()
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Source/ExternalLibraries)
    include_directories(
      ${CMAKE_CURRENT_SOURCE_DIR}/Source/ExternalLibraries/vkFFT)

    if(TARGET glslang::glslang)
      get_target_property(
        _glslang_incs glslang::glslang
        INTERFACE_INCLUDE_DIRECTORIES)
    endif()

    find_path(GLSLANG_C_INTERFACE_INCLUDE_DIR
      NAMES glslang_c_interface.h
      HINTS ${_glslang_incs}
      PATH_SUFFIXES glslang/Include
      REQUIRED
    )

    add_library(vkfft INTERFACE)
    target_compile_definitions(vkfft INTERFACE
      -DVKFFT_BACKEND=0)
    target_include_directories(vkfft INTERFACE
      ${CMAKE_CURRENT_SOURCE_DIR}/ExternalLibraries/vkFFT
      ${GLSLANG_C_INTERFACE_INCLUDE_DIR}
    )
    add_executable(${EXECUTABLE_NAME} Source/Core.cpp)
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE
      glslang::SPIRV
      glslang::glslang
      Vulkan::Vulkan
      vkfft
      vpp::vpp
      miniz)
elseif(UNIX AND NOT APPLE)
    if(FLATPAK)
        include(Source/BuildResources/cmake/flatpak.cmake)
    elseif(CLI)
        include(Source/BuildResources/cmake/linux_macros.cmake)
        include(Source/BuildResources/cmake/linux_cli.cmake)
    else()
        include(Source/BuildResources/cmake/linux_macros.cmake)
        include(Source/BuildResources/cmake/linux.cmake)
    endif()
elseif(APPLE)
    include(Source/BuildResources/cmake/linux_macros.cmake)
    include(Source/BuildResources/cmake/mac.cmake)
elseif(WIN32)
    include(Source/BuildResources/cmake/windows.cmake)
else()
    message(FATAL_ERROR "Could not resolve target OS.")
endif()
