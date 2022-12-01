// define a compiler-specific activeDevice class and set of macros
#ifdef __CUDACC__
#include "LWEActiveDeviceCUDA.cuh"
#elif defined RUNONSYCL
#include "LWEActiveDeviceSYCL.h"
#else
#include "LWEActiveDeviceCPU.h"
#endif

#define deviceFFTD2Z 0
#define deviceFFTZ2D 1
#define deviceFFTD2Z1D 2
#define deviceFFTZ2D1D 3
#define deviceFFTD2ZPolarization 4