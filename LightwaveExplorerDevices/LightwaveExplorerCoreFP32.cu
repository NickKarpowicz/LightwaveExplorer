#undef RUNONSYCL
#undef CPUONLY
#define LWEFLOATINGPOINT 32
#include "LightwaveExplorerCoreFP32.cuh"
#include "../LightwaveExplorerCore.cu"
#undef LWEFLOATINGPOINT

