// dllmain.cpp : Defines the entry point for the DLL application.
#define RUNONSYCL
#include "pch.h"
#include "LightwaveExplorerUtilities.h"
#include "LightwaveExplorerSYCL.h"
#undef max
#undef min
#include "LightwaveExplorerCore.cu"



int readSYCLDevices(wchar_t* deviceListString, wchar_t* defaultDeviceString) {
    wchar_t deviceCharString[MAX_LOADSTRING] = { 0 };
    size_t convertedChars;
    size_t offset = 0;
    int deviceCount = 0;
    for (const auto& p : cl::sycl::platform::get_platforms()) {
        for (const auto& d : p.get_devices()) {
            //loop through all devices, but only mention the GPUs and CPUs (maybe add accelerators later if there's something
            //useful to target and not emulators)
            if (d.is_cpu()) {
                deviceCount++;
                mbstowcs_s(&convertedChars, deviceCharString, d.get_info<cl::sycl::info::device::name>().c_str(), MAX_LOADSTRING);
                offset = wcsnlen_s(deviceListString, MAX_LOADSTRING);
                swprintf_s(&deviceListString[offset], MAX_LOADSTRING, L"SYCL found a CPU: %ls\r\n", deviceCharString);
                memset(deviceCharString, 0, MAX_LOADSTRING * sizeof(wchar_t));
            }
            if (d.is_gpu()) {
                deviceCount++;
                mbstowcs_s(&convertedChars, deviceCharString, d.get_info<cl::sycl::info::device::name>().c_str(), MAX_LOADSTRING);
                offset = wcsnlen_s(deviceListString, MAX_LOADSTRING);
                swprintf_s(&deviceListString[offset], MAX_LOADSTRING, L"SYCL found a GPU: %ls\r\n", deviceCharString);
                memset(deviceCharString, 0, MAX_LOADSTRING * sizeof(wchar_t));
            }
        }
    }

    cl::sycl::default_selector ddefault;
    cl::sycl::queue q(ddefault);
    mbstowcs_s(&convertedChars, deviceCharString, q.get_device().get_info<cl::sycl::info::device::name>().c_str(), MAX_LOADSTRING);
    swprintf_s(defaultDeviceString, MAX_LOADSTRING, L"SYCL likes this one: %ls\r\n", deviceCharString);
    return deviceCount;
}
BOOL APIENTRY DllMain( HMODULE hModule,
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

