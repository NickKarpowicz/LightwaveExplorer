#define RUNONSYCL
#define isnan(x) std::isnan(x)
#include "../LightwaveExplorerUtilities.h"
#include "LightwaveExplorerSYCLLinux.h"
#include "../LightwaveExplorerCore.cu"
#include <stdio.h>
#include <string>
#include <algorithm>
#include <vector>

extern "C" void readSYCLDevices(char* deviceArray, char* deviceListCstring) {
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
                    deviceList.append("SYCL found a CPU:\n");
                    deviceList.append(d.get_info<sycl::info::device::name>());
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
                    deviceList.append("SYCL found a GPU:\n");
                    deviceList.append(d.get_info<sycl::info::device::name>());
                }
            }
        }
    }
    deviceArray[0] = cpuCount;
    deviceArray[1] = gpuCount;
    deviceList.copy(deviceListCstring, 1023);
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
