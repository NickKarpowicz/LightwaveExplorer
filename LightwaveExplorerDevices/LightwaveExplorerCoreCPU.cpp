#undef __CUDACC__
#undef RUNONSYCL
#define LWEFLOATINGPOINT 64
#include "LightwaveExplorerCoreCPU.h"
#include "../LightwaveExplorerCore.cu"
#undef kernels
#undef deviceFunctions
#undef hostFunctions

