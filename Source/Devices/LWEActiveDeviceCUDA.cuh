#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>
#include <thrust/complex.h>
#include <iostream>
#include "../LightwaveExplorerUtilities.h"

namespace complexLib = thrust;
#if LWEFLOATINGPOINT==64
//In tests this mattered, since Thrust does math between complex and double up casting the double to a complex.
__device__ static thrust::complex<double> operator/(const double& a, const thrust::complex<double>& b) {
	double divByDenominator = a / (b.real() * b.real() + b.imag() * b.imag());
	return thrust::complex<double>(b.real() * divByDenominator, -b.imag() * divByDenominator);
}
__device__ static thrust::complex<double> operator/(
	const thrust::complex<double>& a, 
	const double& b) { 
	return thrust::complex<double>(a.real() / b, a.imag() / b); 
}
__device__ static thrust::complex<double> operator*(
	const double& b, 
	const thrust::complex<double>& a) { 
	return thrust::complex<double>(a.real() * b, a.imag() * b); }
__device__ static thrust::complex<double> operator*(
	thrust::complex<double> a, 
	double b) { 
	return thrust::complex<double>(a.real() * b, a.imag() * b); 
}
__device__ static thrust::complex<double> operator+(
	const double& a, 
	const thrust::complex<double>& b) { 
	return thrust::complex<double>(b.real() + a, b.imag()); 
}
__device__ static thrust::complex<double> operator+(
	const thrust::complex<double>& a, 
	const double& b) { 
	return thrust::complex<double>(a.real() + b, a.imag()); 
}
__device__ static thrust::complex<double> operator-(
	const double& a, 
	const thrust::complex<double>& b) { 
	return thrust::complex<double>(a - b.real(), -b.imag()); 
}
__device__ static thrust::complex<double> operator-(
	const thrust::complex<double>& a, 
	const double& b) { 
	return thrust::complex<double>(a.real() - b, a.imag()); 
}

namespace deviceLib = thrust;
#define deviceFPLib
const auto CUFFT_fwd = CUFFT_D2Z;
const auto CUFFT_bwd = CUFFT_Z2D;
#else
namespace deviceLibCUDAFP32{
	__device__ static float exp(const float x){
		return expf(x);
	}
	__device__ static float abs(const float x){
		return fabs(x);
	}
	__device__ static float sin(const float x){
		return sinf(x);
	}
	__device__ static float cos(const float x){
		return cosf(x);
	}
	__device__ static float atan(const float x){
		return atanf(x);
	}
	__device__ static float sqrt(const float x){
		return sqrtf(x);
	}
	__device__ static float asin(const float x){
		return asinf(x);
	}
	__device__ static float pow(const float x, const float y){
		return powf(x,y);
	}
	__device__ static float atan2(const float x, const float y) {
		return atan2f(x, y);
	}
	__device__ static float acos(const float x) {
		return acosf(x);
	}
	__device__ static thrust::complex<float> pow(const thrust::complex<float> x, const float y){
		return thrust::pow(x,y);
	}
	__device__ static thrust::complex<double> pow(const thrust::complex<double> x, const float y) {
		return thrust::pow(x, (double)y);
	}
	__device__ static thrust::complex<float> exp(const thrust::complex<float> x){
		return thrust::exp(x);
	}
	__device__ static float abs(const thrust::complex<float> x){
		return thrust::abs(x);
	}
	__device__ static thrust::complex<float> sqrt(const thrust::complex<float> x){
		return thrust::sqrt(x);
	}
	__device__ static float hypot(const float x, const float y) {
		return hypotf(x, y);
	}
};

namespace deviceLib = deviceLibCUDAFP32;
namespace deviceFPLib = deviceLibCUDAFP32;
const auto CUFFT_fwd = CUFFT_R2C;
const auto CUFFT_bwd = CUFFT_C2R;
#endif

__device__ static inline float j0Device(float x) {
	return j0(x);
}

__device__ static inline double j0Device(double x) {
	return j0(x);
}

static int hardwareCheck(int* CUDAdeviceCount) {
	int CUDAdevice;
	cudaGetDeviceCount(CUDAdeviceCount);
	cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
	struct cudaDeviceProp activeCUDADeviceProp;
	if (cuErr == cudaSuccess) {
		std::cout << "Found " << *CUDAdeviceCount << "GPU(s) : " << std::endl;
		for (int i = 0; i < *CUDAdeviceCount; ++i) {
			cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
			std::cout << activeCUDADeviceProp.name << std::endl;
			std::cout << " Memory: " << 
				(int)(activeCUDADeviceProp.totalGlobalMem / (1024 * 1024)) <<
				"MB; Multiprocessors: " << 
				activeCUDADeviceProp.multiProcessorCount
				<<  std::endl;
		}
	}
	else {
		std::cout << "No GPU found." << std::endl;
		return 1;
	}
	return 0;
}

template <typename T>
__global__ static void deviceLaunchFunctorWrapper(const T functor) {
	functor(threadIdx.x + blockIdx.x * blockDim.x);
}

template<typename deviceFP, typename deviceComplex>
class CUDADevice : public LWEDevice {
private:
	cudaError_t deviceSetError;
	cudaError_t streamGenerationError;
	cufftHandle fftPlanD2Z;
	cufftHandle fftPlanZ2D;
	cufftHandle fftPlan1DD2Z;
	cufftHandle fftPlan1DZ2D;
	cufftHandle doublePolfftPlan;

	void fftDestroy() {
		cufftDestroy(fftPlanD2Z);
		cufftDestroy(fftPlanZ2D);
		cufftDestroy(fftPlan1DD2Z);
		cufftDestroy(fftPlan1DZ2D);
		if (s->isCylindric) {
			cufftDestroy(doublePolfftPlan);
		}
		configuredFFT = 0;
	}

public:
	cudaStream_t stream;
	deviceParameterSet<deviceFP, deviceComplex>* s;
	deviceParameterSet<deviceFP, deviceComplex>* dParamsDevice;
	std::unique_ptr<UPPEAllocation<deviceFP, deviceComplex>> allocation;
	deviceFP canaryPixel = 0.0;
	using LWEDevice::deviceMemcpy;
	CUDADevice(simulationParameterSet* sCPU) :
		deviceSetError(cudaSetDevice(sCPU->assignedGPU)),
		streamGenerationError(cudaStreamCreate(&stream))
	{
		allocateSet(sCPU);
	}

	~CUDADevice() {
		fftDestroy();
		cudaStreamDestroy(stream);
	}

	bool isTheCanaryPixelNaN(deviceFP* canaryPointer) {
		cudaMemcpyAsync(&canaryPixel, canaryPointer, sizeof(deviceFP), cudaMemcpyDeviceToHost);
		return(isnan(canaryPixel));
	}

	template <typename T>
	void deviceLaunch (
		const unsigned int Nblock, 
		const unsigned int Nthread, 
		const T& functor) const {
		deviceLaunchFunctorWrapper<<<Nblock, Nthread, 0, stream>>>(functor);
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) override {
		int err = cudaMalloc(ptr, N * elementSize);
		cudaMemset(*ptr, 0, N * elementSize);
		return err;
	}
	void deviceMemset(void* ptr, int value, size_t count) override {
		cudaMemset(ptr, value, count);
	}
	void deviceMemcpyImplementation(void* dst, const void* src, size_t count, copyType kind) override {
		cudaMemcpy(dst, src, count, static_cast<cudaMemcpyKind>(static_cast<int>(kind)));
	}

	void inline deviceMemcpy(
		thrust::complex<float>* dst, 
		const std::complex<double>* src, 
		size_t count, 
		copyType kind) {
			deviceMemcpy(reinterpret_cast<std::complex<float>*>(dst), src, count, kind);
		}

	void inline deviceMemcpy(
		std::complex<double>* dst, 
		const thrust::complex<float>* src, 
		size_t count, 
		copyType kind) {
			deviceMemcpy(dst, reinterpret_cast<const std::complex<float>*>(src), count, kind);
		}

	void deviceFree(void* block) override {
		cudaFree(block);
	}

	void fftInitialize() {
		if (configuredFFT) {
			fftDestroy();
		}

		size_t workSize;
		cufftPlan1d(&fftPlan1DD2Z, (int)(*s).Ntime, CUFFT_fwd, 2 * (int)((*s).Nspace * (*s).Nspace2));
		cufftPlan1d(&fftPlan1DZ2D, (int)(*s).Ntime, CUFFT_bwd, 2 * (int)((*s).Nspace * (*s).Nspace2));
		cufftSetStream(fftPlan1DD2Z, stream);
		cufftSetStream(fftPlan1DZ2D, stream);
		if ((*s).is3D) {
			int cufftSizes1[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
			cufftCreate(&fftPlanD2Z);
			cufftGetSizeMany(fftPlanD2Z, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_fwd, 2, &workSize);
			cufftMakePlanMany(fftPlanD2Z, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_fwd, 2, &workSize);

			cufftCreate(&fftPlanZ2D);
			cufftGetSizeMany(fftPlanZ2D, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_bwd, 2, &workSize);
			cufftMakePlanMany(fftPlanZ2D, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_bwd, 2, &workSize);
		}
		else {
			int cufftSizes1[] = { (int)(*s).Nspace, (int)(*s).Ntime };

			cufftCreate(&fftPlanD2Z);
			cufftGetSizeMany(fftPlanD2Z, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_fwd, 2, &workSize);
			cufftMakePlanMany(fftPlanD2Z, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_fwd, 2, &workSize);

			cufftCreate(&fftPlanZ2D);
			cufftGetSizeMany(fftPlanZ2D, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_bwd, 2, &workSize);
			cufftMakePlanMany(fftPlanZ2D, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_bwd, 2, &workSize);

			if ((*s).isCylindric) {
				int cufftSizes2[] = { 2 * (int)(*s).Nspace, (int)(*s).Ntime };
				cufftCreate(&doublePolfftPlan);
				cufftGetSizeMany(
					doublePolfftPlan, 2, cufftSizes2, NULL, 
					0, 0, 0, 0, 0, CUFFT_fwd, 2 + 2 * (*s).hasPlasma, &workSize);
				cufftMakePlanMany(
					doublePolfftPlan, 2, cufftSizes2, NULL, 
					0, 0, 0, 0, 0, CUFFT_fwd, 2 + 2 * (*s).hasPlasma, &workSize);
				cufftSetStream(doublePolfftPlan, stream);
			}
		}
		cufftSetStream(fftPlanD2Z, stream);
		cufftSetStream(fftPlanZ2D, stream);
		configuredFFT = true;
	}

	void fft(const void* input, void* output, deviceFFT type) override {
		if(sizeof(deviceFP) == sizeof(double)){
			switch (type){
			case deviceFFT::D2Z:
				cufftExecD2Z(fftPlanD2Z, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
				break;
			case deviceFFT::Z2D:
				cufftExecZ2D(fftPlanZ2D, (cufftDoubleComplex*)input, (cufftDoubleReal*)output);
				break;
			case deviceFFT::D2Z_1D:
				cufftExecD2Z(fftPlan1DD2Z, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
				break;
			case deviceFFT::Z2D_1D:
				cufftExecZ2D(fftPlan1DZ2D, (cufftDoubleComplex*)input, (cufftDoubleReal*)output);
				break;
			case deviceFFT::D2Z_Polarization:
				cufftExecD2Z(doublePolfftPlan, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
				break;
			}
		}
		else{
			switch(type){
			case deviceFFT::D2Z:
				cufftExecR2C(fftPlanD2Z, (cufftReal*)input, (cufftComplex*)output);
				break;
			case deviceFFT::Z2D:
				cufftExecC2R(fftPlanZ2D, (cufftComplex*)input, (cufftReal*)output);
				break;
			case deviceFFT::D2Z_1D:
				cufftExecR2C(fftPlan1DD2Z, (cufftReal*)input, (cufftComplex*)output);
				break;
			case deviceFFT::Z2D_1D:
				cufftExecC2R(fftPlan1DZ2D, (cufftComplex*)input, (cufftReal*)output);
				break;
			case deviceFFT::D2Z_Polarization:
				cufftExecR2C(doublePolfftPlan, (cufftReal*)input, (cufftComplex*)output);
				break;
			}
		}
	}

	void reset(simulationParameterSet* sCPU) override {
		bool resetFFT = (s->hasPlasma != sCPU->hasPlasma());
		allocation->useNewParameterSet(sCPU);
		if(resetFFT){
			fftInitialize();
		}
	}
	int allocateSet(simulationParameterSet* sCPU) {
		cParams = sCPU;
		if(memoryStatus == 0) allocation = nullptr;
		allocation = std::make_unique<UPPEAllocation<deviceFP, deviceComplex>>(this, sCPU);
		s = &(allocation->parameterSet);
		fftInitialize();
		dParamsDevice = allocation->parameterSet_deviceCopy.device_ptr();
		memoryStatus = 0;
		return 0;
	}

};