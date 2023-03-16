#undef RUNONSYCL
#include "LightwaveExplorerUtilities.h"
unsigned long	solveNonlinearWaveEquationCPU(void* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPU(void* lpParam);
int             mainCPU(int argc, char* filepath);
unsigned long	runDlibFittingCPU(simulationParameterSet* sCPU);
unsigned long	solveNonlinearWaveEquationCPUFP32(void* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPUFP32(void* lpParam);
int             mainCPUFP32(int argc, char* filepath);
unsigned long	runDlibFittingCPUFP32(simulationParameterSet* sCPU);