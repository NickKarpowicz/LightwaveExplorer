#undef RUNONSYCL
#include "../LightwaveExplorerUtilities.h"
unsigned long	solveNonlinearWaveEquationCPU(simulationParameterSet* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPU(simulationParameterSet* lpParam);
int             mainCPU(int argc, char* filepath);
unsigned long	runDlibFittingCPU(simulationParameterSet* sCPU);