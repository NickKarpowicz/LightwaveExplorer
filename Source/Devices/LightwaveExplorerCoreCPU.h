#undef RUNONSYCL
#include "../LightwaveExplorerUtilities.h"
#include "../LightwaveExplorerInterfaceClasses.hpp"
unsigned long	solveNonlinearWaveEquationCPU(SimulationParameterSet* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPU(SimulationParameterSet* lpParam);
unsigned long	runDlibFittingCPU(SimulationParameterSet* sCPU);
unsigned long	solveNonlinearWaveEquationCPUFP32(SimulationParameterSet* lpParam);
unsigned long   solveNonlinearWaveEquationSequenceCPUFP32(SimulationParameterSet* lpParam);
unsigned long	runDlibFittingCPUFP32(SimulationParameterSet* sCPU);
