#pragma once
#if defined LIGHTWAVEEXPLORERSYCL_EXPORTS && defined _WIN32
#define LIGHTWAVEEXPLORERSYCL_API extern "C" __declspec(dllexport)
#elif defined _WIN32
#define LIGHTWAVEEXPLORERSYCL_API extern "C" __declspec(dllimport)
#else
#define LIGHTWAVEEXPLORERSYCL_API
#endif
LIGHTWAVEEXPLORERSYCL_API unsigned long	solveNonlinearWaveEquationSYCL(void* lpParam);
LIGHTWAVEEXPLORERSYCL_API unsigned long solveNonlinearWaveEquationSequenceSYCL(void* lpParam);
LIGHTWAVEEXPLORERSYCL_API int           mainSYCL(int argc, char* filepath);
LIGHTWAVEEXPLORERSYCL_API unsigned long	runDlibFittingSYCL(simulationParameterSet* sCPU);
LIGHTWAVEEXPLORERSYCL_API void readSYCLDevices(char* deviceArray, char* deviceListCstring);