#define RUNONSYCL
#define isnan(x) std::isnan(x)
#include "LightwaveExplorerUtilities.h"
#include "LightwaveExplorerDPCPPlib.h"
#include "LightwaveExplorerCore.cu"
#include <stdio.h>

extern "C" size_t readSYCLDevices(char* deviceListString, char* defaultDeviceString) {
    size_t convertedChars;
    size_t offset = 0;
    unsigned char cpuCount = 0;
    unsigned char gpuCount = 0;
    size_t deviceCount = 0;
    unsigned char* deviceArray = (unsigned char*)&deviceCount;
    //try{
    auto platforms = sycl::platform::get_platforms();
    for (const auto& p : platforms) {
        for (const auto& d : p.get_devices()) {
            //loop through all devices, but only mention the GPUs and CPUs (maybe add accelerators later if there's something
            //useful to target and not emulators)
            if (d.is_cpu()) {
                cpuCount++;

                offset = strnlen(deviceListString, MAX_LOADSTRING);
                sprintf(&deviceListString[offset], "SYCL found a CPU: %s\r\n", d.get_info<sycl::info::device::name>().c_str());

            }
            if (d.is_gpu()) {
                gpuCount++;

                offset = strnlen(deviceListString, MAX_LOADSTRING);
                sprintf(&deviceListString[offset], "SYCL found a GPU: %s\r\n", d.get_info<sycl::info::device::name>().c_str());

            }
        }
    }
    //}
    // catch(int whatever){
    //     printf("Couldn't get a platform...\n");
    // }

    

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
