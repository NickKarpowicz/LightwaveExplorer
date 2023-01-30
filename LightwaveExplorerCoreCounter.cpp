#undef __CUDACC__
#undef RUNONSYCL
#define RUNSTEPCOUNTER
#define kernels kernelsCounter
#define deviceFunctions deviceFunctionsCounter
#define hostFunctions hostFunctionsCounter
#include "LightwaveExplorerCoreCounter.h"
#include "LightwaveExplorerCore.cu"
#undef RUNSTEPCOUNTER
#undef kernels
#undef deviceFunctions
#undef hostFunctions
