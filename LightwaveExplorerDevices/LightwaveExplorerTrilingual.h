// define a compiler-specific activeDevice class and set of macros
#ifdef __CUDACC__
#include "LWEActiveDeviceCUDA.cuh"
#define deviceFunction __device__
#define localIndex threadIdx.x + blockIdx.x * blockDim.x
#define deviceFunctions deviceFunctionsCUDA
#define hostFunctions hostFunctionsCUDA

#define kernelLWE(kernelName,...) \
__global__ static void kernelName(__VA_ARGS__)
#if LWEFLOATINGPOINT==64
#ifndef NOCUDAMAIN
#define mainX main
#else
#define mainX mainCUDA
#endif
#define runDlibFittingX runDlibFitting
#define solveNonlinearWaveEquationX solveNonlinearWaveEquation
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequence
#else
#ifdef YESCUDAMAIN32
#define mainX main
#else
#define mainX mainCUDAFP32
#endif
#define runDlibFittingX runDlibFittingFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceFP32
#endif
#elif defined RUNONSYCL
#include "LWEActiveDeviceSYCL.h"
#define kernelLWE(kernelName,...) \
const auto kernelName = [](size_t trilingualLaunchID, __VA_ARGS__)
#define deviceFunction 
#define localIndex trilingualLaunchID

#if LWEFLOATINGPOINT == 64
#define mainX mainSYCL
#define runDlibFittingX runDlibFittingSYCL
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCL
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCL
#else
#define mainX mainSYCLFP32
#define runDlibFittingX runDlibFittingSYCLFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCLFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCLFP32
#endif

#elif defined RUNSTEPCOUNTER
#include "LWEAciveDeviceCounter.h"
#define kernelLWE(kernelName,...) \
static void kernelName(size_t trilingualLaunchID, __VA_ARGS__)
#define deviceFunction 
#define localIndex trilingualLaunchID
#define mainX mainCounter
#define runDlibFittingX runDlibFittingCounter
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCounter
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCounter
#else
#include "LWEActiveDeviceCPU.h"

#define kernelLWE(kernelName,...) \
static void kernelName(size_t trilingualLaunchID, __VA_ARGS__)
#define deviceFunction 
#define localIndex trilingualLaunchID

#if LWEFLOATINGPOINT == 32
#define mainX mainCPUFP32
#define runDlibFittingX runDlibFittingCPUFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPUFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPUFP32
#else
#define mainX mainCPU
#define runDlibFittingX runDlibFittingCPU
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPU
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPU
#endif
#endif

