#undef RUNONSYCL
#include "../LightwaveExplorerUtilities.h"
unsigned long	solveNonlinearWaveEquationCPU(SimulationParameterSet* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPU(SimulationParameterSet* lpParam);
int             mainCPU(int argc, char* filepath);
unsigned long	runDlibFittingCPU(SimulationParameterSet* sCPU);
