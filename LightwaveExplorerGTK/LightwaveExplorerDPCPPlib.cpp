#define RUNONSYCL
#define isnan(x) std::isnan(x)
#include "LightwaveExplorerUtilities.h"
#include "LightwaveExplorerDPCPPlib.h"
#include "LightwaveExplorerCore.cu"
#include <stdio.h>

extern "C" size_t readSYCLDevices(wchar_t* deviceListString, wchar_t* defaultDeviceString) {
    size_t convertedChars;
    size_t offset = 0;
    unsigned char cpuCount = 0;
    unsigned char gpuCount = 0;

    for (const auto& p : sycl::platform::get_platforms()) {
        for (const auto& d : p.get_devices()) {
            //loop through all devices, but only mention the GPUs and CPUs (maybe add accelerators later if there's something
            //useful to target and not emulators)
            if (d.is_cpu()) {
                cpuCount++;
                #ifndef __linux__
                offset = wcsnlen(deviceListString, MAX_LOADSTRING);
                swprintf(&deviceListString[offset], MAX_LOADSTRING, L"SYCL found a CPU: %s\r\n", d.get_info<sycl::info::device::name>().c_str());
                #endif
            }
            if (d.is_gpu()) {
                gpuCount++;
                #ifndef __linux__
                offset = wcsnlen(deviceListString, MAX_LOADSTRING);
                swprintf(&deviceListString[offset], MAX_LOADSTRING, L"SYCL found a GPU: %s\r\n", d.get_info<sycl::info::device::name>().c_str());
                #endif
            }
        }
    }
    
    size_t deviceCount = 0;
    unsigned char* deviceArray = (unsigned char*)&deviceCount;
    deviceArray[0] = cpuCount;
    deviceArray[1] = gpuCount;
    return deviceCount;
}


void listSyclDevices() {
    printf("Scanning for devices...\n");
    for (const auto& p : sycl::platform::get_platforms()) {
        for (const auto& d : p.get_devices()) {
            printf("Found: %s \n", d.get_info<sycl::info::device::name>().c_str());
        }
    }
    printf("Scan finished\n");
}
