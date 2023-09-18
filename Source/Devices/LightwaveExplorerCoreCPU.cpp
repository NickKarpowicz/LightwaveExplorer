#undef __CUDACC__
#undef RUNONSYCL
#define LWEFLOATINGPOINT 64
#include "../LightwaveExplorerCore.cu"
#undef kernels
#undef deviceFunctions
#undef hostFunctions
#undef LWEFLOATINGPOINT

