#pragma once
#include "LightwaveExplorerUtilities.h"
unsigned long	solveNonlinearWaveEquationCPU(void* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPU(void* lpParam);
unsigned long   runFittingCPU(simulationParameterSet* sCPU);
int             mainCPU(int argc, char* filepath);
