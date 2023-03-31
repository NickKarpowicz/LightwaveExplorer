// define a compiler-specific device class and set of macros
#ifdef __CUDACC__
#include "LWEActiveDeviceCUDA.cuh"
#define deviceFunction __device__
#if LWEFLOATINGPOINT==64
#define kernelNamespace CUDA64Kernels
typedef CUDADevice<double, thrust::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef complexLib::complex<double> deviceComplex;
#define runDlibFittingX runDlibFitting
#define solveNonlinearWaveEquationX solveNonlinearWaveEquation
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequence
#else
#define kernelNamespace CUDA32Kernels
typedef CUDADevice<float, thrust::complex<float>> ActiveDevice;
typedef float deviceFP;
typedef complexLib::complex<float> deviceComplex;
#define runDlibFittingX runDlibFittingFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceFP32
#endif
#elif defined RUNONSYCL
#include "LWEActiveDeviceSYCL.h"
#define deviceFunction 
#if LWEFLOATINGPOINT == 64
#define kernelNamespace SYCL64Kernels
typedef SYCLDevice<double, oneapi::dpl::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef complexLib::complex<double> deviceComplex;
#define runDlibFittingX runDlibFittingSYCL
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCL
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCL
#else
#define kernelNamespace SYCL32Kernels
typedef SYCLDevice<float, oneapi::dpl::complex<float>> ActiveDevice;
typedef float deviceFP;
typedef complexLib::complex<float> deviceComplex;
#define runDlibFittingX runDlibFittingSYCLFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCLFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCLFP32
#endif
#elif defined RUNSTEPCOUNTER
#include "LWEAciveDeviceCounter.h"
#define kernelNamespace CounterKernels
#define deviceFunction 
typedef counterDevice<double, std::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef complexLib::complex<double> deviceComplex;
#define runDlibFittingX runDlibFittingCounter
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCounter
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCounter
#else
#include "LWEActiveDeviceCPU.h"
#define deviceFunction 
#if LWEFLOATINGPOINT == 32
#define kernelNamespace CPU32Kernels
typedef CPUDevice<float, std::complex<float>> ActiveDevice;
typedef float deviceFP;
typedef complexLib::complex<float> deviceComplex;
#define runDlibFittingX runDlibFittingCPUFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPUFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPUFP32
#else
#define kernelNamespace CPU64Kernels
typedef CPUDevice<double, std::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef complexLib::complex<double> deviceComplex;
#define runDlibFittingX runDlibFittingCPU
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPU
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPU
#endif
#endif

