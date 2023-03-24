#undef RUNONSYCL
#include "LightwaveExplorerUtilities.h"
unsigned long	solveNonlinearWaveEquationCPU(simulationParameterSet* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPU(simulationParameterSet* lpParam);
unsigned long	runDlibFittingCPU(simulationParameterSet* sCPU);
unsigned long	solveNonlinearWaveEquationCPUFP32(simulationParameterSet* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPUFP32(simulationParameterSet* lpParam);
unsigned long	runDlibFittingCPUFP32(simulationParameterSet* sCPU);