#pragma once
#ifdef LIGHTWAVEEXPLORERSYCL_EXPORTS
#define LIGHTWAVEEXPLORERSYCL_API __declspec(dllexport)
#else
#define LIGHTWAVEEXPLORERSYCL_API __declspec(dllimport)
#endif
extern "C" LIGHTWAVEEXPLORERSYCL_API unsigned long	solveNonlinearWaveEquationSYCL(void* lpParam);
extern "C" LIGHTWAVEEXPLORERSYCL_API unsigned long  solveNonlinearWaveEquationSequenceSYCL(void* lpParam);
extern "C" LIGHTWAVEEXPLORERSYCL_API int            mainSYCL(int argc, char* filepath);
extern "C" LIGHTWAVEEXPLORERSYCL_API unsigned long	runDlibFittingSYCL(simulationParameterSet* sCPU);
extern "C" LIGHTWAVEEXPLORERSYCL_API size_t			readSYCLDevices(wchar_t* deviceListString, wchar_t* defaultDeviceString);