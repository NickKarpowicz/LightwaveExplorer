#undef __CUDACC__
#undef RUNONSYCL
#define RUNSTEPCOUNTER
#define LWEFLOATINGPOINT 64
#include "LightwaveExplorerCoreCounter.h"
#ifdef __linux__
bool isnan(double x){
    return std::isnan(x);
}
#endif
#include "../LightwaveExplorerCore.cu"
#undef RUNSTEPCOUNTER

