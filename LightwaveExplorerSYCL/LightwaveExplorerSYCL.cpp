// dllmain.cpp : Defines the entry point for the DLL application.
#define RUNONSYCL

#include "LightwaveExplorerUtilities.h"
#include "LightwaveExplorerSYCL.h"
#include "LightwaveExplorerCore.cu"

size_t readSYCLDevices(char* deviceListString, char* defaultDeviceString) {
    char deviceCharString[MAX_LOADSTRING] = { 0 };
    size_t offset = 0;
    unsigned char cpuCount = 0;
    unsigned char gpuCount = 0;
    for (const auto& p : cl::sycl::platform::get_platforms()) {
        for (const auto& d : p.get_devices()) {
            //loop through all devices, but only mention the GPUs and CPUs (maybe add accelerators later if there's something
            //useful to target and not emulators)
            if (d.is_cpu()) {
                cpuCount++;
                offset = strnlen_s(deviceListString, MAX_LOADSTRING);
                sprintf_s(&deviceListString[offset], MAX_LOADSTRING, "SYCL found a CPU: %s\r\n", d.get_info<cl::sycl::info::device::name>().c_str());
                memset(deviceCharString, 0, MAX_LOADSTRING * sizeof(char));
            }
            if (d.is_gpu()) {
                gpuCount++;
                offset = strnlen_s(deviceListString, MAX_LOADSTRING);
                sprintf_s(&deviceListString[offset], MAX_LOADSTRING, "SYCL found a GPU: %s\r\n", d.get_info<cl::sycl::info::device::name>().c_str());
                memset(deviceCharString, 0, MAX_LOADSTRING * sizeof(char));
            }
        }
    }
    size_t deviceCount = 0;
    unsigned char* deviceArray = (unsigned char*)&deviceCount;
    deviceArray[0] = cpuCount;
    deviceArray[1] = gpuCount;
    cl::sycl::default_selector ddefault;
    cl::sycl::queue q(ddefault);
    sprintf_s(defaultDeviceString, MAX_LOADSTRING, "SYCL default: %s\r\n", q.get_device().get_info<cl::sycl::info::device::name>().c_str());
    return deviceCount;
}

#ifdef _WIN32
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
    return TRUE;
}
#endif
