// define a compiler-specific device class and set of macros

#ifndef LWEFLOATINGPOINT
#define LWEFLOATINGPOINT 64
#endif

#ifdef __CUDACC__
#include <cufft.h>
#include <thrust/complex.h>
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
#include <oneapi/dpl/complex>
#include <oneapi/dpl/cmath>
#include <oneapi/mkl/dfti.hpp>
#include <sycl/sycl.hpp>
const auto dftPrecision = (LWEFLOATINGPOINT == 32) ? oneapi::mkl::dft::precision::SINGLE : oneapi::mkl::dft::precision::DOUBLE;
#include "LWEActiveDeviceSYCL.h"
#define deviceFunction 
#if LWEFLOATINGPOINT == 64
namespace deviceLib = oneapi::dpl;
namespace complexLib = oneapi::dpl;
#define deviceFPLib
#define kernelNamespace SYCL64Kernels
typedef SYCLDevice<double, oneapi::dpl::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef oneapi::dpl::complex<double> deviceComplex;
#define runDlibFittingX runDlibFittingSYCL
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCL
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCL
#else
namespace deviceLib = deviceLibSYCLFP32;
namespace deviceFPLib = deviceLibSYCLFP32;
namespace complexLib = oneapi::dpl;
#define kernelNamespace SYCL32Kernels
typedef SYCLDevice<float, oneapi::dpl::complex<float>> ActiveDevice;
typedef float deviceFP;
typedef oneapi::dpl::complex<float> deviceComplex;
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
namespace complexLib = std;
#define deviceFunction 
#if LWEFLOATINGPOINT == 32
namespace deviceFPLib = deviceLibCPUFP32;
namespace deviceLib = deviceLibCPUFP32;
#define kernelNamespace CPU32Kernels
typedef CPUDevice<float, std::complex<float>> ActiveDevice;
typedef float deviceFP;
typedef complexLib::complex<float> deviceComplex;
#define runDlibFittingX runDlibFittingCPUFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPUFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPUFP32
#else
namespace deviceFPLib = std;
namespace deviceLib = std;
#define kernelNamespace CPU64Kernels
typedef CPUDevice<double, std::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef complexLib::complex<double> deviceComplex;
#define runDlibFittingX runDlibFittingCPU
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPU
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPU
#endif
#endif

