// define a compiler-specific activeDevice class and set of macros
#ifdef __CUDACC__
#include "LWEActiveDeviceCUDA.cuh"
#define trilingual __global__ static void
#define withID 
#define asKernel
#define deviceFunction __device__
#define localIndex threadIdx.x + blockIdx.x * blockDim.x
#define withStream 
#define activeDevice deviceCUDA
#define hardwareCheck hardwareCheckCUDA
#define deviceFunctions deviceFunctionsCUDA
#define atomicAddDevice atomicAdd
#define j0Device j0
#define hostFunctions hostFunctionsCUDA

#define mainArgumentX char* argv[]
#define resolveArgv char* filepath = argv[1];
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
#define trilingual const auto 
#define deviceFunction 
#define localIndex trilingualLaunchID
#define asKernel = []
#define withID size_t trilingualLaunchID,
#define withStream , stream
#define activeDevice deviceSYCL
#define atomicAddDevice atomicAddSYCL
#define j0Device j0SYCL
#define hardwareCheck hardwareCheckSYCL

#define deviceFunctions deviceFunctionsSYCL
#define hostFunctions hostFunctionsSYCL
#define mainArgumentX char* filepath
#define resolveArgv
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
#define trilingual static void 
#define withID size_t trilingualLaunchID, 
#define asKernel
#define deviceFunction 
#define localIndex trilingualLaunchID
#define withStream 
#define activeDevice deviceCounter
#define atomicAddDevice atomicAddCounter
#define j0Device j0Counter
#define hardwareCheck hardwareCheckCounter
#define deviceFunctions deviceFunctionsCounter
#define hostFunctions hostFunctionsCounter
#define mainX mainCounter
#define mainArgumentX char* filepath
#define resolveArgv
#define runDlibFittingX runDlibFittingCounter
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCounter
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCounter
#else
#include "LWEActiveDeviceCPU.h"
#define trilingual static void 
#define withID size_t trilingualLaunchID, 
#define asKernel
#define deviceFunction 
#define localIndex trilingualLaunchID
#define withStream 
#define activeDevice deviceCPU
#define atomicAddDevice atomicAddCPU
#define j0Device j0CPU
#define hardwareCheck hardwareCheckCPU
#define deviceFunctions deviceFunctionsCPU
#define hostFunctions hostFunctionsCPU
#define mainArgumentX char* filepath
#define resolveArgv

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

#define deviceFFTD2Z 0
#define deviceFFTZ2D 1
#define deviceFFTD2Z1D 2
#define deviceFFTZ2D1D 3
#define deviceFFTD2ZPolarization 4