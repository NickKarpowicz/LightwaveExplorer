#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "LightwaveExplorerUtilities.h"
#include <stdio.h>
#include <complex>
#include "cufft.h"
#include <thrust/complex.h>

#define MAX_LOADSTRING 1024

unsigned long   solveNonlinearWaveEquationSequence(void* lpParam);
unsigned long	solveNonlinearWaveEquation(void* lpParam);
unsigned long   runFitting(simulationParameterSet* sCPU);



