#include "LWEActiveDeviceCommon.cpp"
#include <fftw3_mkl.h>
#include <atomic>
#include <thread>
#define DeviceToHost 2
#define HostToDevice 1
#define DeviceToDevice 3
#define cudaMemcpyKind int
#ifdef _WIN32
const int deviceThreads = maxN(1, std::thread::hardware_concurrency()/2);
#else
const int deviceThreads = std::thread::hardware_concurrency();
#endif
#if defined CPUONLY || NOCUDA
#define isnan(x) std::isnan(x)
#endif

#if defined __APPLE__ || defined __linux__
void atomicAddCPU(double* pulseSum, double pointEnergy) {
	std::atomic<double>* pulseSumAtomic = (std::atomic<double>*)pulseSum;
	double expected = pulseSumAtomic->load();
	while (!std::atomic_compare_exchange_weak(pulseSumAtomic, &expected, expected + pointEnergy));
}
#else
void atomicAddCPU(double* pulseSum, double pointEnergy) {
	std::atomic<double>* pulseSumAtomic = (std::atomic<double>*)pulseSum;
	(*pulseSumAtomic).fetch_add(pointEnergy);
}
#endif
int hardwareCheckCPU(int* CUDAdeviceCount) {
	*CUDAdeviceCount = 1;
	return 0;
}
class deviceCPU {
private:
	bool configuredFFT = FALSE;
	bool isCylindric = FALSE;
	double canaryPixel = 0.0;
	deviceParameterSet dParamslocal;
	fftw_plan fftPlanD2Z;
	fftw_plan fftPlanZ2D;
	fftw_plan fftPlan1DD2Z;
	fftw_plan fftPlan1DZ2D;
	fftw_plan doublePolfftPlan;

	void fftDestroy() {
		fftw_destroy_plan(fftPlanD2Z);
		fftw_destroy_plan(fftPlanZ2D);
		fftw_destroy_plan(fftPlan1DD2Z);
		fftw_destroy_plan(fftPlan1DZ2D);
		if (isCylindric)fftw_destroy_plan(doublePolfftPlan);
		fftw_cleanup();
	}

public:
	int stream;
	int memoryStatus;
	deviceParameterSet* dParams;
	simulationParameterSet* cParams;
	deviceParameterSet* dParamsDevice;
	deviceCPU() {
		memoryStatus = -1;
		stream = 0;
		dParams = NULL;
		cParams = NULL;
		dParamsDevice = NULL;
		configuredFFT = 0;
		isCylindric = 0;
		fftPlanD2Z = 0;
		fftPlanZ2D = 0;
		fftPlan1DD2Z = 0;
		fftPlan1DZ2D = 0;
		doublePolfftPlan = 0;
		dParamsDevice = &dParamslocal;
	}

	deviceCPU(simulationParameterSet* sCPU, deviceParameterSet* s) {
		memoryStatus = -1;
		stream = 0;
		dParams = NULL;
		cParams = NULL;
		dParamsDevice = NULL;
		configuredFFT = 0;
		isCylindric = 0;
		fftPlanD2Z = 0;
		fftPlanZ2D = 0;
		fftPlan1DD2Z = 0;
		fftPlan1DZ2D = 0;
		doublePolfftPlan = 0;
		dParamsDevice = &dParamslocal;
		memoryStatus = allocateSet(sCPU, s);
	}

	~deviceCPU() {
		fftDestroy();
	}
	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) {
#pragma omp parallel for num_threads(deviceThreads)
		for (int i = 0; i < (int)Nthread; ++i) {
			for (unsigned int j = 0u; j < Nblock; ++j) {
				kernel(j + Nblock * (unsigned int)i, args...);
			}
		}
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) {
		(*ptr) = calloc(N, elementSize);
		return (int)((*ptr) == NULL);
	}

	void deviceMemset(void* ptr, int value, size_t count) {
		memset(ptr, value, count);
	}

	void deviceMemcpy(void* dst, void* src, size_t count, cudaMemcpyKind kind) {
		memcpy(dst, src, count);
	}

	void deviceFree(void* block) {
		free(block);
	}

	bool isTheCanaryPixelNaN(double* canaryPointer) {
		deviceMemcpy(&canaryPixel, canaryPointer, sizeof(double), DeviceToHost);
		return(isnan(canaryPixel));
	}
	void fft(void* input, void* output, int type) {
		if (!configuredFFT) return;
		switch (type) {
		case 0:
			fftw_execute_dft_r2c(fftPlanD2Z, (double*)input, (fftw_complex*)output);
			break;
		case 1:
			fftw_execute_dft_c2r(fftPlanZ2D, (fftw_complex*)input, (double*)output);
			break;
		case 2:
			fftw_execute_dft_r2c(fftPlan1DD2Z, (double*)input, (fftw_complex*)output);
			break;
		case 3:
			fftw_execute_dft_c2r(fftPlan1DZ2D, (fftw_complex*)input, (double*)output);
			break;
		case 4:
			fftw_execute_dft_r2c(doublePolfftPlan, (double*)input, (fftw_complex*)output);
			break;
		}
	}

	void fftInitialize(deviceParameterSet* s) {
		if (configuredFFT) fftDestroy();
		isCylindric = (*s).isCylindric;
		double* setupWorkD = new double[(*s).Ngrid * (2 + 2 * (*s).isCylindric)];
		std::complex<double>* setupWorkC = new std::complex<double>[(*s).NgridC * (2 + 2 * (*s).isCylindric)];
		const int fftw1[1] = { (int)(*s).Ntime };
		fftPlan1DD2Z = fftw_plan_many_dft_r2c(1, fftw1, (int)(*s).Nspace * (int)(*s).Nspace2 * 2, setupWorkD, fftw1, 1, (int)(*s).Ntime, (fftw_complex*)setupWorkC, fftw1, 1, (int)(*s).Nfreq, FFTW_ESTIMATE);
		fftPlan1DZ2D = fftw_plan_many_dft_c2r(1, fftw1, (int)(*s).Nspace * (int)(*s).Nspace2 * 2, (fftw_complex*)setupWorkC, fftw1, 1, (int)(*s).Nfreq, setupWorkD, fftw1, 1, (int)(*s).Ntime, FFTW_ESTIMATE);
		if ((*s).is3D) {
			const int fftwSizes[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
			fftPlanD2Z = fftw_plan_many_dft_r2c(3, fftwSizes, 2, setupWorkD, NULL, 1, (int)(*s).Ngrid, (fftw_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, FFTW_MEASURE);
			fftPlanZ2D = fftw_plan_many_dft_c2r(3, fftwSizes, 2, (fftw_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, setupWorkD, NULL, 1, (int)(*s).Ngrid, FFTW_MEASURE);
		}
		else {
			const int fftwSizes[] = { (int)(*s).Nspace, (int)(*s).Ntime };
			fftPlanD2Z = fftw_plan_many_dft_r2c(2, fftwSizes, 2, setupWorkD, NULL, 1, (int)(*s).Ngrid, (fftw_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, FFTW_MEASURE);
			fftPlanZ2D = fftw_plan_many_dft_c2r(2, fftwSizes, 2, (fftw_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, setupWorkD, NULL, 1, (int)(*s).Ngrid, FFTW_MEASURE);

			if ((*s).isCylindric) {
				const int fftwSizesCyl[] = { (int)(2 * (*s).Nspace), (int)(*s).Ntime };
				doublePolfftPlan = fftw_plan_many_dft_r2c(2, fftwSizesCyl, 2, setupWorkD, NULL, 1, 2 * (int)(*s).Ngrid, (fftw_complex*)setupWorkC, NULL, 1, 2 * (int)(*s).NgridC, FFTW_MEASURE);
			}
		}
		delete[] setupWorkC;
		delete[] setupWorkD;
		configuredFFT = 1;
	}
	void deallocateSet(deviceParameterSet* s) {
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
	int allocateSet(simulationParameterSet* sCPU, deviceParameterSet* s) {
		dParams = s;
		cParams = sCPU;
		initializeDeviceParameters(sCPU, s);
		fftInitialize(s);

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
		deviceMemcpy(dParamsDevice, s, sizeof(deviceParameterSet), HostToDevice);
		return 0;
	}
};