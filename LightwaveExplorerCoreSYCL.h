#pragma once
#include "LightwaveExplorerUtilities.h"
unsigned long	solveNonlinearWaveEquationSYCL(void* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceSYCL(void* lpParam);
int             mainSYCL(int argc, char* filepath);
unsigned long	runDlibFittingSYCL(simulationParameterSet* sCPU);