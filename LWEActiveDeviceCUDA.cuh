#include "LWEActiveDeviceCommon.cpp"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>
#include <nvml.h>
#include <memory>
#define trilingual __global__ void
#define withID 
#define asKernel
#define deviceFunction __device__
#define RUNTYPE 0
#define localIndex threadIdx.x + blockIdx.x * blockDim.x
#define withStream 
#define activeDevice deviceCUDA
#define DeviceToHost cudaMemcpyDeviceToHost
#define HostToDevice cudaMemcpyHostToDevice
#define DeviceToDevice cudaMemcpyDeviceToDevice

#define deviceFunctions deviceFunctionsCUDA
#define hostFunctions hostFunctionsCUDA
#define mainX main
#define mainArgumentX char* argv[]
#define resolveArgv char* filepath = argv[1];
#define runDlibFittingX runDlibFitting
#define solveNonlinearWaveEquationX solveNonlinearWaveEquation
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequence

class deviceCUDA {
private:
	bool configuredFFT = FALSE;
	bool isCylindric = FALSE;
	double canaryPixel = 0.0;
	cufftHandle fftPlanD2Z;
	cufftHandle fftPlanZ2D;
	cufftHandle fftPlan1DD2Z;
	cufftHandle fftPlan1DZ2D;
	cufftHandle doublePolfftPlan;

	int checkDeviceMemory(simulationParameterSet* sCPU, cudaParameterSet* s) {
		nvmlDevice_t nvmlDevice = 0;
		nvmlMemory_t nvmlMemoryInfo;
		nvmlInit_v2();
		nvmlDeviceGetHandleByIndex_v2((*sCPU).assignedGPU, &nvmlDevice);
		nvmlDeviceGetMemoryInfo(nvmlDevice, &nvmlMemoryInfo);
		size_t memoryEstimate = sizeof(double) * ((*s).Ngrid * 2 * 2 + 2 * (*s).NgridC * 6 * 2 + 2 * (*s).isCylindric * 5 * 2 + 2 * (*s).Ntime + 2 * (*s).Nfreq + 81 + 65536);
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
	cudaParameterSet* dParams;
	cudaParameterSet* dParamsDevice;
	simulationParameterSet* cParams;
	int memoryStatus;
	deviceCUDA() {
		memoryStatus = -1;
		configuredFFT = 0;
		isCylindric = 0;
		cudaStreamCreate(&stream);
		deviceCalloc((void**)&dParamsDevice, 1, sizeof(cudaParameterSet));
	}

	deviceCUDA(simulationParameterSet* sCPU, cudaParameterSet* s) {
		memoryStatus = -1;
		configuredFFT = 0;
		isCylindric = 0;
		cudaStreamCreate(&stream);
		deviceCalloc((void**)&dParamsDevice, 1, sizeof(cudaParameterSet));
		memoryStatus = allocateSet(sCPU, s);
	}


	~deviceCUDA() {
		fftDestroy();
		deviceFree(dParamsDevice);
		cudaStreamDestroy(stream);
	}

	bool isTheCanaryPixelNaN(double* canaryPointer) {
		cudaMemcpyAsync(&canaryPixel, canaryPointer, sizeof(double), DeviceToHost);
		return(isnan(canaryPixel));
	}

	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) {
		kernel << <Nblock, Nthread, 0, stream >> > (args...);
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

	void deviceFree(void* block) {
		cudaFree(block);
	}

	void fftInitialize(cudaParameterSet* s) {
		if (configuredFFT) {
			fftDestroy();
		}
		isCylindric = 0;
		size_t workSize;
		cufftPlan1d(&fftPlan1DD2Z, (int)(*s).Ntime, CUFFT_D2Z, 2 * (int)((*s).Nspace * (*s).Nspace2));
		cufftPlan1d(&fftPlan1DZ2D, (int)(*s).Ntime, CUFFT_Z2D, 2 * (int)((*s).Nspace * (*s).Nspace2));
		cufftSetStream(fftPlan1DD2Z, stream);
		cufftSetStream(fftPlan1DZ2D, stream);
		if ((*s).is3D) {
			int cufftSizes1[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
			cufftCreate(&fftPlanD2Z);
			cufftGetSizeMany(fftPlanD2Z, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
			cufftMakePlanMany(fftPlanD2Z, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);

			cufftCreate(&fftPlanZ2D);
			cufftGetSizeMany(fftPlanZ2D, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);
			cufftMakePlanMany(fftPlanZ2D, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);
		}
		else {
			int cufftSizes1[] = { (int)(*s).Nspace, (int)(*s).Ntime };

			cufftCreate(&fftPlanD2Z);
			cufftGetSizeMany(fftPlanD2Z, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
			cufftMakePlanMany(fftPlanD2Z, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);

			cufftCreate(&fftPlanZ2D);
			cufftGetSizeMany(fftPlanZ2D, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);
			cufftMakePlanMany(fftPlanZ2D, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);

			if ((*s).isCylindric) {
				isCylindric = 1;
				int cufftSizes2[] = { 2 * (int)(*s).Nspace, (int)(*s).Ntime };
				cufftCreate(&doublePolfftPlan);
				cufftGetSizeMany(doublePolfftPlan, 2, cufftSizes2, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
				cufftMakePlanMany(doublePolfftPlan, 2, cufftSizes2, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
				cufftSetStream(doublePolfftPlan, stream);
			}
		}
		cufftSetStream(fftPlanD2Z, stream);
		cufftSetStream(fftPlanZ2D, stream);
		configuredFFT = 1;
	}

	void fft(void* input, void* output, int type) {
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
		}
	}
	void deallocateSet(cudaParameterSet* s) {
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

	int allocateSet(simulationParameterSet* sCPU, cudaParameterSet* s) {
		dParams = s;
		cParams = sCPU;
		cudaSetDevice((*sCPU).assignedGPU);
		initializeDeviceParameters(sCPU, s);
		fftInitialize(s);
		if (checkDeviceMemory(sCPU, s)) {
			fftDestroy();
			return 1;
		}

		int memErrors = 0;
		double* expGammaTCPU = new double[2 * (*s).Ntime];

		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
		}

		fillRotationMatricies(sCPU, s);

		//GPU allocations
		//
		// currently 8 large grids, meaning memory use is approximately
		// 64 bytes per grid point (8 grids x 2 polarizations x 4ouble precision)
		// plus a little bit for additional constants/workspaces/etc
		memErrors += deviceCalloc((void**)&(*s).gridETime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationTime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).workspace1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1Next1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).k1, 2 * (*s).NgridC, sizeof(std::complex<double>));

		//cylindric sym grids
		if ((*s).isCylindric) {
			memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1Rho1, 4 * (*s).NgridC, sizeof(std::complex<double>));
			memErrors += deviceCalloc((void**)&(*s).gridRadialLaplacian1, 4 * (*s).Ngrid, sizeof(std::complex<double>));
		}

		//smaller helper grids
		memErrors += deviceCalloc((void**)&(*s).expGammaT, 2 * (*s).Ntime, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).fieldFactor1, 2 * (*s).Nfreq, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).inverseChiLinear1, 2 * (*s).Nfreq, sizeof(double));
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
		finishConfiguration(sCPU, s);
		deviceMemcpy(dParamsDevice, s, sizeof(cudaParameterSet), HostToDevice);
		return 0;
	}

};