// define a compiler-specific device class and set of macros
#ifdef __CUDACC__
#include "LWEActiveDeviceCUDA.cuh"
#define deviceFunction __device__
#define localIndex threadIdx.x + blockIdx.x * blockDim.x
#define kernelLWE(kernelName,...) \
__global__ static void kernelName(__VA_ARGS__)
#if LWEFLOATINGPOINT==64
typedef CUDADevice<double, thrust::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef complexLib::complex<double> deviceComplex;
#define runDlibFittingX runDlibFitting
#define solveNonlinearWaveEquationX solveNonlinearWaveEquation
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequence
#else
typedef CUDADevice<float, thrust::complex<float>> ActiveDevice;
typedef float deviceFP;
typedef complexLib::complex<float> deviceComplex;
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
typedef SYCLDevice<double, oneapi::dpl::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef complexLib::complex<double> deviceComplex;
#define runDlibFittingX runDlibFittingSYCL
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCL
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCL
#else
typedef SYCLDevice<float, oneapi::dpl::complex<float>> ActiveDevice;
typedef float deviceFP;
typedef complexLib::complex<float> deviceComplex;
#define runDlibFittingX runDlibFittingSYCLFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCLFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCLFP32
#endif
#elif defined RUNSTEPCOUNTER
#include "LWEAciveDeviceCounter.h"
#define kernelLWE(kernelName,...) \
static void kernelName(size_t trilingualLaunchID, __VA_ARGS__)
#define deviceFunction 
typedef counterDevice<double, std::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef complexLib::complex<double> deviceComplex;
#define localIndex trilingualLaunchID
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
typedef CPUDevice<float, std::complex<float>> ActiveDevice;
typedef float deviceFP;
typedef complexLib::complex<float> deviceComplex;
#define runDlibFittingX runDlibFittingCPUFP32
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPUFP32
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPUFP32
#else
typedef CPUDevice<double, std::complex<double>> ActiveDevice;
typedef double deviceFP;
typedef complexLib::complex<double> deviceComplex;
#define runDlibFittingX runDlibFittingCPU
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPU
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPU
#endif
#endif

