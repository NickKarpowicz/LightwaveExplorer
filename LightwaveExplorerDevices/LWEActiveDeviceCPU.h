#include "LWEActiveDeviceCommon.cpp"
#ifdef CPUONLY
#include <fftw3.h>
#else
#include <fftw3_mkl.h>
#endif
#include <atomic>
#include <thread>
#define DeviceToHost 2
#define HostToDevice 1
#define DeviceToDevice 3
#define cudaMemcpyKind int
#if defined _WIN32 || __linux__
const int deviceThreads = maxN(1, std::thread::hardware_concurrency()/2);
#else
const int deviceThreads = std::thread::hardware_concurrency();
#endif

#ifdef __linux__
#define isnan std::isnan
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
double j0CPU(double x) {
	if (x < 8.0) {
		double y = x * x;
		double ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7 +
			y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
		double ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718 +
			y * (59272.64853 + y * (267.8532712 + y))));
		return ans1 / ans2;
	}
	else {
		double z = 8.0 / x;
		double y = z * z;
		double xx = x - 0.785398164;
		double ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 +
			y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
		double ans2 = -0.1562499995e-1 + y * (0.1430488765e-3 +
			y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7)));
		return sqrt(0.636619772 / x) * (cos(xx) * ans1 - z * sin(xx) * ans2);
	}
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
	bool hasPlasma = FALSE;
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
	const void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) const {
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
	const void fft(void* input, void* output, int type) const {
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
		hasPlasma = (*s).hasPlasma;
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
				doublePolfftPlan = fftw_plan_many_dft_r2c(2, fftwSizesCyl, 2 + 2 * (*s).hasPlasma, setupWorkD, NULL, 1, 2 * (int)(*s).Ngrid, (fftw_complex*)setupWorkC, NULL, 1, 2 * (int)(*s).NgridC, FFTW_MEASURE);
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
	void reset(simulationParameterSet* sCPU, deviceParameterSet* s) {
		if((*s).hasPlasma != ((*sCPU).nonlinearAbsorptionStrength != 0.0)){
			deallocateSet(s);
			memoryStatus = allocateSet(sCPU,s);
		}
		else{
			initializeDeviceParameters(sCPU, s);
		}
		
		fillRotationMatricies(sCPU, s);
		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
			if ((*s).hasPlasma) beamExpansionFactor = 4;
		}
		deviceMemset((*s).gridETime1, 0, 2 * (*s).Ngrid * sizeof(double));
		deviceMemset((*s).gridPolarizationTime1, 0, 2 * (*s).Ngrid * sizeof(double));
		deviceMemset((*s).workspace1, 0, beamExpansionFactor * 2 * (*s).NgridC * sizeof(std::complex<double>));
		deviceMemset((*s).gridEFrequency1, 0, 2 * (*s).NgridC * sizeof(std::complex<double>));
		deviceMemset((*s).gridPropagationFactor1, 0, 2 * (*s).NgridC * sizeof(std::complex<double>));
		deviceMemset((*s).gridPolarizationFactor1, 0, 2 * (*s).NgridC * sizeof(std::complex<double>));
		deviceMemset((*s).gridEFrequency1Next1, 0, 2 * (*s).NgridC * sizeof(std::complex<double>));
		deviceMemset((*s).k1, 0, 2 * (*s).NgridC * sizeof(std::complex<double>));

		//cylindric sym grids
		if ((*s).isCylindric) {
			deviceMemset((*s).gridPropagationFactor1Rho1, 0, 4 * (*s).NgridC * sizeof(std::complex<double>));
			deviceMemset((*s).gridRadialLaplacian1, 0, 2 * beamExpansionFactor * (*s).Ngrid * sizeof(std::complex<double>));
		}

		//smaller helper grids
		deviceMemset((*s).expGammaT, 0, 2 * (*s).Ntime * sizeof(double));
		deviceMemset((*s).chiLinear1, 0, 2 * (*s).Nfreq * sizeof(std::complex<double>));
		deviceMemset((*s).fieldFactor1, 0, 2 * (*s).Nfreq * sizeof(double));
		deviceMemset((*s).inverseChiLinear1, 0, 2 * (*s).Nfreq * sizeof(double));

		double* expGammaTCPU = new double[2 * (*s).Ntime];
		for (size_t i = 0; i < (*s).Ntime; ++i) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, HostToDevice);
		delete[] expGammaTCPU;

		finishConfiguration(sCPU, s);
		deviceMemcpy(dParamsDevice, s, sizeof(deviceParameterSet), HostToDevice);
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
			if ((*s).hasPlasma) beamExpansionFactor = 4;
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
			memErrors += deviceCalloc((void**)&(*s).gridRadialLaplacian1, 2 * beamExpansionFactor * (*s).Ngrid, sizeof(std::complex<double>));
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