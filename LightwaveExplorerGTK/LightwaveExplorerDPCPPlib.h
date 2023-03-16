#pragma once
extern "C" unsigned long	solveNonlinearWaveEquationSYCL(void* lpParam);
extern "C" unsigned long    solveNonlinearWaveEquationSequenceSYCL(void* lpParam);
extern "C" int              mainSYCL(int argc, char* filepath);
extern "C" unsigned long	runDlibFittingSYCL(simulationParameterSet * sCPU);
extern "C" unsigned long	solveNonlinearWaveEquationSYCLFP32(void* lpParam);
extern "C" unsigned long    solveNonlinearWaveEquationSequenceSYCLFP32(void* lpParam);
extern "C" int              mainSYCLFP32(int argc, char* filepath);
extern "C" unsigned long	runDlibFittingSYCLFP32(simulationParameterSet * sCPU);
extern "C" void             readSYCLDevices(char* deviceListString, char* defaultDeviceString);
void listSyclDevices();