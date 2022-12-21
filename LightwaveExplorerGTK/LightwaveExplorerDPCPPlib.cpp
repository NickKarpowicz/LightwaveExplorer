#define RUNONSYCL
#include "LightwaveExplorerUtilities.h"
bool isnan(double x) {
	return std::isnan(x);
}
#include "LightwaveExplorerDPCPPlib.h"
#include "LightwaveExplorerCore.cu"
#include <stdio.h>

extern "C" size_t readSYCLDevices(wchar_t* deviceListString, wchar_t* defaultDeviceString) {
    printf("Yes, I'm looking\n");
    size_t convertedChars;
    size_t offset = 0;
    unsigned char cpuCount = 0;
    unsigned char gpuCount = 0;
    for (const auto& p : sycl::platform::get_platforms()) {
        for (const auto& d : p.get_devices()) {
            //loop through all devices, but only mention the GPUs and CPUs (maybe add accelerators later if there's something
            //useful to target and not emulators)
            printf("Found something: %s \n", d.get_info<sycl::info::device::name>().c_str());
            if (d.is_cpu()) {
                cpuCount++;
                offset = wcsnlen(deviceListString, MAX_LOADSTRING);
                swprintf(&deviceListString[offset], MAX_LOADSTRING, L"SYCL found a CPU: %s\r\n", d.get_info<sycl::info::device::name>().c_str());
                printf("SYCL found a CPU: %s\r\n", d.get_info<sycl::info::device::name>().c_str());
            }
            if (d.is_gpu()) {
                gpuCount++;
                offset = wcsnlen(deviceListString, MAX_LOADSTRING);
                swprintf(&deviceListString[offset], MAX_LOADSTRING, L"SYCL found a GPU: %s\r\n", d.get_info<sycl::info::device::name>().c_str());
                printf("SYCL found a GPU: %s\r\n", d.get_info<sycl::info::device::name>().c_str());
            }
        }
    }
    size_t deviceCount = 0;
    unsigned char* deviceArray = (unsigned char*)&deviceCount;
    deviceArray[0] = cpuCount;
    deviceArray[1] = gpuCount;
    return deviceCount;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Running default\n");
        char defaultPath[] = "DefaultValues.ini\0";
        mainSYCL(2, defaultPath);
        return 1;
    }
    mainSYCL(argc, argv[1]);
    return 0;
}