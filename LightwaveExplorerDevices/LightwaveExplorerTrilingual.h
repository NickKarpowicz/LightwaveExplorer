// define a compiler-specific activeDevice class and set of macros
#ifdef __CUDACC__
#include "LWEActiveDeviceCUDA.cuh"
#define deviceFunction __device__
#define localIndex threadIdx.x + blockIdx.x * blockDim.x
#define activeDevice CUDADevice
#define kernelLWE(kernelName,...) \
__global__ static void kernelName(__VA_ARGS__)
#if LWEFLOATINGPOINT==64

#define runDlibFittingX runDlibFitting
#define solveNonlinearWaveEquationX solveNonlinearWaveEquation
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequence
#else

#define runDlibFittingX runDlibFittingFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceFP32
#endif
#elif defined RUNONSYCL
#include "LWEActiveDeviceSYCL.h"
#define activeDevice SYCLDevice
#define kernelLWE(kernelName,...) \
const auto kernelName = [](size_t trilingualLaunchID, __VA_ARGS__)
#define deviceFunction 
#define localIndex trilingualLaunchID

#if LWEFLOATINGPOINT == 64
#define runDlibFittingX runDlibFittingSYCL
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCL
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCL
#else

#define runDlibFittingX runDlibFittingSYCLFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCLFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCLFP32
#endif

#elif defined RUNSTEPCOUNTER
#include "LWEAciveDeviceCounter.h"
#define activeDevice counterDevice
#define kernelLWE(kernelName,...) \
static void kernelName(size_t trilingualLaunchID, __VA_ARGS__)
#define deviceFunction 
#define localIndex trilingualLaunchID
#define runDlibFittingX runDlibFittingCounter
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCounter
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCounter
#else
#include "LWEActiveDeviceCPU.h"
#define activeDevice CPUDevice
#define kernelLWE(kernelName,...) \
static void kernelName(size_t trilingualLaunchID, __VA_ARGS__)
#define deviceFunction 
#define localIndex trilingualLaunchID

#if LWEFLOATINGPOINT == 32
#define runDlibFittingX runDlibFittingCPUFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPUFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPUFP32
#else
#define runDlibFittingX runDlibFittingCPU
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPU
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPU
#endif
#endif

