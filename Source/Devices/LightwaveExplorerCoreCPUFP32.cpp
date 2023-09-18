#undef __CUDACC__
#undef RUNONSYCL
#define LWEFLOATINGPOINT 32
#include "LightwaveExplorerCoreCPU.h"
#include "../LightwaveExplorerCore.cu"
#undef kernels
#undef deviceFunctions
#undef hostFunctions

