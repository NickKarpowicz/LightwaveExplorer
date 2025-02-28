#undef __CUDACC__
#undef RUNONSYCL
#define LWEFLOATINGPOINT 32
#include "VisualizationKernels.cpp"
#undef kernels
#undef deviceFunctions
#undef hostFunctions
#undef LWEFLOATINGPOINT