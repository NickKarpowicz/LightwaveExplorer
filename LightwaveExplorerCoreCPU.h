#undef RUNONSYCL
#include "LightwaveExplorerUtilities.h"
unsigned long	solveNonlinearWaveEquationCPU(void* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPU(void* lpParam);
int             mainCPU(int argc, char* filepath);
unsigned long	runDlibFittingCPU(simulationParameterSet* sCPU);