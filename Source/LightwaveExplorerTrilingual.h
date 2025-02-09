#pragma once
#include "MaxwellConfiguration.hpp"
#ifndef LWEFLOATINGPOINT
#define LWEFLOATINGPOINT 64
#endif

#if defined __CUDACC__ || defined RUNONCUDA
	#include <cufft.h>
	#include <thrust/complex.h>
	#include "Devices/LWEActiveDeviceCUDA.cuh"
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
	#include <oneapi/mkl/dft.hpp>
	#include <sycl/sycl.hpp>
	const auto dftPrecision = (LWEFLOATINGPOINT == 32) ? 
		oneapi::mkl::dft::precision::SINGLE 
		: oneapi::mkl::dft::precision::DOUBLE;
	#include "Devices/LWEActiveDeviceSYCL.h"
	#define deviceFunction 
	#if LWEFLOATINGPOINT == 64
		namespace deviceLib = sycl::ext::oneapi::experimental;
		namespace complexLib = sycl::ext::oneapi::experimental;
		#define deviceFPLib
		#define kernelNamespace SYCL64Kernels
		typedef SYCLDevice<double, sycl::ext::oneapi::experimental::complex<double>> ActiveDevice;
		typedef double deviceFP;
		typedef sycl::ext::oneapi::experimental::complex<double> deviceComplex;
		#define runDlibFittingX runDlibFittingSYCL
		#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCL
		#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCL
	#else
		namespace deviceLib = deviceLibSYCLFP32;
		namespace deviceFPLib = deviceLibSYCLFP32;
		namespace complexLib = sycl::ext::oneapi::experimental;
		#define kernelNamespace SYCL32Kernels
		typedef SYCLDevice<float, sycl::ext::oneapi::experimental::complex<float>> ActiveDevice;
		typedef float deviceFP;
		typedef sycl::ext::oneapi::experimental::complex<float> deviceComplex;
		#define runDlibFittingX runDlibFittingSYCLFP32
		#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCLFP32
		#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCLFP32
	#endif
#elif defined RUNSTEPCOUNTER
	#include "Devices/LWEAciveDeviceCounter.h"
	#define kernelNamespace CounterKernels
	#define deviceFunction 
	typedef counterDevice<double, std::complex<double>> ActiveDevice;
	typedef double deviceFP;
	typedef complexLib::complex<double> deviceComplex;
	#define runDlibFittingX runDlibFittingCounter
	#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCounter
	#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCounter
#else
	#include "Devices/LWEActiveDeviceCPU.h"
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

typedef 
maxwellCalculation<deviceFP, maxwellPoint<deviceFP>, maxwellPoint<deviceFP>, oscillator<deviceFP>> 
maxwell3D;