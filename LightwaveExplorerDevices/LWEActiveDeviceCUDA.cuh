#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>
#include <nvml.h>
#include <thrust/complex.h>
#include <iostream>
#include "LightwaveExplorerUtilities.h"
#define DeviceToHost cudaMemcpyDeviceToHost
#define HostToDevice cudaMemcpyHostToDevice
#define DeviceToDevice cudaMemcpyDeviceToDevice
#if LWEFLOATINGPOINT==32
	#define deviceLib deviceLibCUDAFP32
	#define deviceFPLib deviceLibCUDAFP32 
	#define complexLib thrust
	#define CUFFT_fwd CUFFT_R2C
	#define CUFFT_bwd CUFFT_C2R
#else
	#define deviceLib thrust
	#define complexLib thrust
	#define deviceFPLib 
	#define CUFFT_fwd CUFFT_D2Z
	#define CUFFT_bwd CUFFT_Z2D
#endif


#ifdef __CUDACC__
#if LWEFLOATINGPOINT==64
//In tests this mattered, since Thrust does math between complex and double up casting the double to a complex.
__device__ static thrust::complex<double> operator/(const double& a, const thrust::complex<double>& b) {
	double divByDenominator = a / (b.real() * b.real() + b.imag() * b.imag());
	return thrust::complex<double>(b.real() * divByDenominator, -b.imag() * divByDenominator);
}
__device__ static thrust::complex<double> operator/(const thrust::complex<double>& a, const double& b) { return thrust::complex<double>(a.real() / b, a.imag() / b); }

//__device__ static thrust::complex<double> operator*(const double& b, const thrust::complex<double>& a) { return thrust::complex<double>(a.real() * b, a.imag() * b); }
__device__ static thrust::complex<double> operator*(thrust::complex<double> a, double b) { return thrust::complex<double>(a.real() * b, a.imag() * b); }

__device__ static thrust::complex<double> operator+(const double& a, const thrust::complex<double>& b) { return thrust::complex<double>(b.real() + a, b.imag()); }
__device__ static thrust::complex<double> operator+(const thrust::complex<double>& a, const double& b) { return thrust::complex<double>(a.real() + b, a.imag()); }

__device__ static thrust::complex<double> operator-(const double& a, const thrust::complex<double>& b) { return thrust::complex<double>(a - b.real(), -b.imag()); }
__device__ static thrust::complex<double> operator-(const thrust::complex<double>& a, const double& b) { return thrust::complex<double>(a.real() - b, a.imag()); }
#endif
#endif


#if LWEFLOATINGPOINT==32
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
};
#endif

static int hardwareCheckCUDA(int* CUDAdeviceCount) {
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

static class deviceCUDA {
	using deviceFP = LWEFLOATINGPOINTTYPE;
	using deviceComplex = thrust::complex<deviceFP>;
private:
#include "LWEActiveDeviceCommon.cpp"
	bool configuredFFT = FALSE;
	bool isCylindric = FALSE;
	deviceFP canaryPixel = 0.0;
	
	cufftHandle fftPlanD2Z;
	cufftHandle fftPlanZ2D;
	cufftHandle fftPlan1DD2Z;
	cufftHandle fftPlan1DZ2D;
	cufftHandle doublePolfftPlan;

	int checkDeviceMemory(simulationParameterSet* sCPU) {
		nvmlDevice_t nvmlDevice = 0;
		nvmlMemory_t nvmlMemoryInfo;
		nvmlInit_v2();
		nvmlDeviceGetHandleByIndex_v2((*sCPU).assignedGPU, &nvmlDevice);
		nvmlDeviceGetMemoryInfo(nvmlDevice, &nvmlMemoryInfo);
		size_t memoryEstimate = sizeof(deviceFP) * ((*s).Ngrid * 2 * 2 + 2 * (*s).NgridC * 6 * 2 + 2 * (*s).isCylindric * 5 * 2 + 2 * (*s).Ntime + 2 * (*s).Nfreq + 81 + 65536);
		nvmlShutdown();
		if (nvmlMemoryInfo.free < memoryEstimate) {
			(*sCPU).memoryError = -1;
			return 1;
		}
		return 0;
	}

	void fftDestroy() {
		cufftDestroy(fftPlanD2Z);
		cufftDestroy(fftPlanZ2D);
		cufftDestroy(fftPlan1DD2Z);
		cufftDestroy(fftPlan1DZ2D);
		if (isCylindric) {
			cufftDestroy(doublePolfftPlan);
		}
		configuredFFT = 0;
	}
public:
	cudaStream_t stream;
	deviceParameterSet<deviceFP, deviceComplex> deviceStruct;
	deviceParameterSet<deviceFP, deviceComplex>* s;
	deviceParameterSet<deviceFP, deviceComplex>* dParamsDevice;
	simulationParameterSet* cParams;
	int memoryStatus;
	bool hasPlasma;

	deviceCUDA(simulationParameterSet* sCPU) {
		s = &deviceStruct;
		memset(s, 0, sizeof(deviceParameterSet<deviceFP, deviceComplex>));
		memoryStatus = -1;
		configuredFFT = 0;
		isCylindric = 0;
		cudaSetDevice((*sCPU).assignedGPU);
		cudaStreamCreate(&stream);
		deviceCalloc((void**)&dParamsDevice, 1, sizeof(deviceParameterSet<deviceFP, deviceComplex>));
		memoryStatus = allocateSet(sCPU);
	}


	~deviceCUDA() {
		fftDestroy();
		deallocateSet();
		deviceFree(dParamsDevice);
		cudaStreamDestroy(stream);
	}

	bool isTheCanaryPixelNaN(deviceFP* canaryPointer) {
		cudaMemcpyAsync(&canaryPixel, canaryPointer, sizeof(deviceFP), DeviceToHost);
		return(isnan(canaryPixel));
	}

	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) const {
		kernel<<<Nblock, Nthread, 0, stream>>>(args...);
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) {
		int err = cudaMalloc(ptr, N * elementSize);
		cudaMemset(*ptr, 0, N * elementSize);
		return err;
	}

	void deviceMemset(void* ptr, int value, size_t count) {
		cudaMemset(ptr, value, count);
	}

	void deviceMemcpy(void* dst, void* src, size_t count, cudaMemcpyKind kind) {
		cudaMemcpy(dst, src, count, kind);
	}
	void deviceMemcpy(double* dst, float* src, size_t count, cudaMemcpyKind kind) {
		float* copyBuffer = new float[count / sizeof(double)];
		cudaMemcpy(copyBuffer, src, count/2, kind);
		for (size_t i = 0; i < count / sizeof(double); i++) {
			dst[i] = copyBuffer[i];
		}
		delete[] copyBuffer;
	}

	void deviceMemcpy(std::complex<double>* dst, thrust::complex<float>* src, size_t count, cudaMemcpyKind kind) {
		thrust::complex<float>* copyBuffer = new thrust::complex<float>[count / sizeof(std::complex<double>)];
		cudaMemcpy(copyBuffer, src, count/2, kind);
		for (size_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			dst[i] = std::complex<double>(copyBuffer[i].real(), copyBuffer[i].imag());
		}
		delete[] copyBuffer;
	}

	void deviceMemcpy(thrust::complex<float>* dst, std::complex<double>* src, size_t count, cudaMemcpyKind kind) {
		thrust::complex<float>* copyBuffer = new thrust::complex<float>[count / sizeof(std::complex<double>)];
		
		for (size_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			copyBuffer[i] = thrust::complex<float>((float)src[i].real(), (float)src[i].imag());
		}
		cudaMemcpy(dst, copyBuffer, count / 2, kind);
		delete[] copyBuffer;
	}

	void deviceMemcpy(float* dst, double* src, size_t count, cudaMemcpyKind kind) {
		float* copyBuffer = new float[count / sizeof(double)];

		for (size_t i = 0; i < count / sizeof(double); i++) {
			copyBuffer[i] = (float)src[i];
		}
		cudaMemcpy(dst, copyBuffer, count / 2, kind);
		delete[] copyBuffer;
	}


	void deviceFree(void* block) {
		cudaFree(block);
	}

	void fftInitialize() {
		if (configuredFFT) {
			fftDestroy();
		}
		isCylindric = 0;
		hasPlasma = (*s).hasPlasma;
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
				isCylindric = 1;
				int cufftSizes2[] = { 2 * (int)(*s).Nspace, (int)(*s).Ntime };
				cufftCreate(&doublePolfftPlan);
				cufftGetSizeMany(doublePolfftPlan, 2, cufftSizes2, NULL, 0, 0, 0, 0, 0, CUFFT_fwd, 2 + 2 * (*s).hasPlasma, &workSize);
				cufftMakePlanMany(doublePolfftPlan, 2, cufftSizes2, NULL, 0, 0, 0, 0, 0, CUFFT_fwd, 2 + 2 * (*s).hasPlasma, &workSize);
				cufftSetStream(doublePolfftPlan, stream);
			}
		}
		cufftSetStream(fftPlanD2Z, stream);
		cufftSetStream(fftPlanZ2D, stream);
		configuredFFT = 1;
	}

	void fft(void* input, void* output, int type) {
		if (sizeof(deviceFP) == sizeof(float)) type += 5;
		switch (type) {
		case 0:
			cufftExecD2Z(fftPlanD2Z, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
			break;
		case 1:
			cufftExecZ2D(fftPlanZ2D, (cufftDoubleComplex*)input, (cufftDoubleReal*)output);
			break;
		case 2:
			cufftExecD2Z(fftPlan1DD2Z, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
			break;
		case 3:
			cufftExecZ2D(fftPlan1DZ2D, (cufftDoubleComplex*)input, (cufftDoubleReal*)output);
			break;
		case 4:
			cufftExecD2Z(doublePolfftPlan, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
			break;
		case 5:
			cufftExecR2C(fftPlanD2Z, (cufftReal*)input, (cufftComplex*)output);
			break;
		case 6:
			cufftExecC2R(fftPlanZ2D, (cufftComplex*)input, (cufftReal*)output);
			break;
		case 7:
			cufftExecR2C(fftPlan1DD2Z, (cufftReal*)input, (cufftComplex*)output);
			break;
		case 8:
			cufftExecC2R(fftPlan1DZ2D, (cufftComplex*)input, (cufftReal*)output);
			break;
		case 9:
			cufftExecR2C(doublePolfftPlan, (cufftReal*)input, (cufftComplex*)output);
			break;
		}
	}
	void deallocateSet() {
		deviceFree((*s).gridETime1);
		deviceFree((*s).workspace1);
		deviceFree((*s).gridEFrequency1);
		deviceFree((*s).gridPropagationFactor1);
		if ((*s).isCylindric) {
			deviceFree((*s).gridPropagationFactor1Rho1);
			deviceFree((*s).gridRadialLaplacian1);
		}
		deviceFree((*s).gridPolarizationFactor1);
		deviceFree((*s).gridEFrequency1Next1);
		deviceFree((*s).k1);
		deviceFree((*s).gridPolarizationTime1);
		deviceFree((*s).expGammaT);
		deviceFree((*s).chiLinear1);
		deviceFree((*s).fieldFactor1);
		deviceFree((*s).inverseChiLinear1);
	}
	void reset(simulationParameterSet* sCPU) {
		if((*s).hasPlasma != ((*sCPU).nonlinearAbsorptionStrength != 0.0)){
			deallocateSet();
			memoryStatus = allocateSet(sCPU);
		}
		else{
			initializeDeviceParameters(sCPU);
		}
		fillRotationMatricies(sCPU);
		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
			if ((*s).hasPlasma) beamExpansionFactor = 4;
		}
		deviceMemset((*s).gridETime1, 0, 2 * (*s).Ngrid * sizeof(deviceFP));
		deviceMemset((*s).gridPolarizationTime1, 0, 2 * (*s).Ngrid * sizeof(deviceFP));
		deviceMemset((*s).workspace1, 0, beamExpansionFactor * 2 * (*s).NgridC * sizeof(deviceComplex));
		deviceMemset((*s).gridEFrequency1, 0, 2 * (*s).NgridC * sizeof(deviceComplex));
		deviceMemset((*s).gridPropagationFactor1, 0,  2 * (*s).NgridC * sizeof(deviceComplex));
		deviceMemset((*s).gridPolarizationFactor1, 0, 2 * (*s).NgridC * sizeof(deviceComplex));
		deviceMemset((*s).gridEFrequency1Next1, 0, 2 * (*s).NgridC * sizeof(deviceComplex));
		deviceMemset((*s).k1, 0, 2 * (*s).NgridC * sizeof(deviceComplex));

		//cylindric sym grids
		if ((*s).isCylindric) {
			deviceMemset((*s).gridPropagationFactor1Rho1, 0, 4 * (*s).NgridC * sizeof(deviceComplex));
			deviceMemset((*s).gridRadialLaplacian1, 0, 2 * beamExpansionFactor * (*s).Ngrid * sizeof(deviceComplex));
		}

		//smaller helper grids
		deviceMemset((*s).expGammaT, 0, 2 * (*s).Ntime * sizeof(deviceFP));
		deviceMemset((*s).chiLinear1, 0, 2 * (*s).Nfreq * sizeof(deviceComplex));
		deviceMemset((*s).fieldFactor1, 0, 2 * (*s).Nfreq * sizeof(deviceFP));
		deviceMemset((*s).inverseChiLinear1, 0, 2 * (*s).Nfreq * sizeof(deviceFP));

		double* expGammaTCPU = new double[2 * (*s).Ntime];
		for (size_t i = 0; i < (*s).Ntime; ++i) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, HostToDevice);
		delete[] expGammaTCPU;

		finishConfiguration(sCPU);
		deviceMemcpy(dParamsDevice, s, sizeof(deviceParameterSet<deviceFP, deviceComplex>), HostToDevice);
	}
	int allocateSet(simulationParameterSet* sCPU) {
		cParams = sCPU;
		cudaSetDevice((*sCPU).assignedGPU);
		initializeDeviceParameters(sCPU);
		fftInitialize();
		if (checkDeviceMemory(sCPU)) {
			fftDestroy();
			return 1;
		}

		int memErrors = 0;
		

		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
			if ((*s).hasPlasma) beamExpansionFactor = 4;
		}

		fillRotationMatricies(sCPU);

		//GPU allocations
		//
		// currently 8 large grids, meaning memory use is approximately
		// 64 bytes per grid point (8 grids x 2 polarizations x 4ouble precision)
		// plus a little bit for additional constants/workspaces/etc
		memErrors += deviceCalloc((void**)&(*s).gridETime1, 2 * (*s).Ngrid, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationTime1, 2 * (*s).Ngrid, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)&(*s).workspace1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationFactor1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1Next1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).k1, 2 * (*s).NgridC, sizeof(deviceComplex));

		//cylindric sym grids
		if ((*s).isCylindric) {
			memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1Rho1, 4 * (*s).NgridC, sizeof(deviceComplex));
			memErrors += deviceCalloc((void**)&(*s).gridRadialLaplacian1, 2 * beamExpansionFactor * (*s).Ngrid, sizeof(deviceComplex));
		}

		//smaller helper grids
		memErrors += deviceCalloc((void**)&(*s).expGammaT, 2 * (*s).Ntime, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).fieldFactor1, 2 * (*s).Nfreq, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)&(*s).inverseChiLinear1, 2 * (*s).Nfreq, sizeof(deviceFP));

		double* expGammaTCPU = new double[2 * (*s).Ntime];
		for (size_t i = 0; i < (*s).Ntime; ++i) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, HostToDevice);
		delete[] expGammaTCPU;

		(*sCPU).memoryError = memErrors;
		if (memErrors > 0) {
			return memErrors;
		}
		finishConfiguration(sCPU);
		deviceMemcpy(dParamsDevice, s, sizeof(deviceParameterSet<deviceFP, deviceComplex>), HostToDevice);
		return 0;
	}

};