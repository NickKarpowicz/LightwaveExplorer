// dllmain.cpp : Defines the entry point for the DLL application.
#define RUNONSYCL
#define LWEFLOATINGPOINT 64
#include "../LightwaveExplorerUtilities.h"
#include "../LightwaveExplorerInterfaceClasses.hpp"
#include "LightwaveExplorerSYCL.h"
#include "LightwaveExplorerCore.cu"
#include <string>
#include <vector>
#if defined __linux__ || defined __APPLE__
#include<fmt/format.h>
#define Sformat fmt::format
#define Svformat fmt::vformat
#define Smake_format_args fmt::make_format_args
#else
#include <format>
#define Sformat std::format
#define Svformat std::vformat
#define Smake_format_args std::make_format_args
#endif

void readSYCLDevices(char* deviceArray, char* deviceListCstring) {
    unsigned char cpuCount = 0;
    unsigned char gpuCount = 0;
    std::vector<std::string> namelist;
    std::string deviceList;
    for (const auto& p : sycl::platform::get_platforms()) {
        for (const auto& d : p.get_devices()) {
            //loop through all devices, but only mention the GPUs and CPUs (maybe add accelerators later if there's something
            //useful to target and not emulators)
            if (d.is_cpu()) {
                if (!(std::find(
                    std::begin(namelist),
                    std::end(namelist),
                    d.get_info<sycl::info::device::name>())
                    != std::end(namelist))) {
                    namelist.push_back(d.get_info<sycl::info::device::name>());
                    cpuCount++;
                    deviceList.append(Sformat("SYCL found a CPU:\n   {}\n", d.get_info<sycl::info::device::name>()));
                }
            }
            if (d.is_gpu()) {
                if (!(std::find(
                    std::begin(namelist),
                    std::end(namelist),
                    d.get_info<sycl::info::device::name>())
                    
                    != std::end(namelist))) {
                    
                    namelist.push_back(d.get_info<sycl::info::device::name>());
                    gpuCount++;
                    deviceList.append(Sformat("SYCL found a GPU:\n   {}\n", d.get_info<sycl::info::device::name>()));
                    if(sizeof(deviceFP) == sizeof(double) && d.get_info<sycl::info::device::double_fp_config>().size() == 0) {
                        //gpuCount--;
                        deviceList.append(Sformat("   Warning: doesn't support FP64\n"));
                    }
                    
                }
            }
        }
    }
    deviceArray[0] = cpuCount;
    deviceArray[1] = gpuCount;
    deviceList.copy(deviceListCstring, 1023);
}

#ifdef _WIN32
#include <windows.h>
BOOL DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return true;
}
#endif
