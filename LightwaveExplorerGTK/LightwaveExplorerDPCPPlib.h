#pragma once
extern "C" unsigned long	solveNonlinearWaveEquationSYCL(void* lpParam);
extern "C" unsigned long    solveNonlinearWaveEquationSequenceSYCL(void* lpParam);
extern "C" int              mainSYCL(int argc, char* filepath);
extern "C" unsigned long	runDlibFittingSYCL(simulationParameterSet * sCPU);
extern "C" size_t           readSYCLDevices(char* deviceListString, char* defaultDeviceString);
void listSyclDevices();