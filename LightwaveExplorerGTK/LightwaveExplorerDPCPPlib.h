#pragma once
extern "C" unsigned long	solveNonlinearWaveEquationSYCL(void* lpParam);
extern "C" unsigned long    solveNonlinearWaveEquationSequenceSYCL(void* lpParam);
extern "C" unsigned long	runDlibFittingSYCL(simulationParameterSet * sCPU);
extern "C" void             readSYCLDevices(char* deviceListString, char* defaultDeviceString);
void listSyclDevices();