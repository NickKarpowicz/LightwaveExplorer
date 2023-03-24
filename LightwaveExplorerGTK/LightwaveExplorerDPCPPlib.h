#pragma once
extern "C" unsigned long	solveNonlinearWaveEquationSYCL(simulationParameterSet* lpParam);
extern "C" unsigned long    solveNonlinearWaveEquationSequenceSYCL(simulationParameterSet* lpParam);
extern "C" unsigned long	runDlibFittingSYCL(simulationParameterSet * sCPU);
extern "C" void             readSYCLDevices(char* deviceListString, char* defaultDeviceString);
void listSyclDevices();