#include "LightwaveExplorerUtilities.h"
#ifdef CPUONLY
#include <fftw3.h>
#else
#include <fftw3_mkl.h>
#endif
#include <atomic>
#include <thread>
#include <iostream>
#include <cmath>
#include <cstring>
#include <ranges>
#include <execution>
#include <numeric>
#include <iterator>
#if defined __APPLE__ || defined __linux__
template<typename deviceFP>
[[maybe_unused]] void atomicAdd(deviceFP* pulseSum, deviceFP pointEnergy) {
	std::atomic<deviceFP>* pulseSumAtomic = (std::atomic<deviceFP>*)pulseSum;
	deviceFP expected = pulseSumAtomic->load();
	while (!std::atomic_compare_exchange_weak(pulseSumAtomic, &expected, expected + pointEnergy));
}
#ifdef __linux__
static double inline isnan(double x){
	return std::isnan(x);
}
static float inline isnan(float x){
	return std::isnan(x);
}
#endif
#else
template<typename deviceFP>
static void atomicAdd(deviceFP* pulseSum, deviceFP pointEnergy) {
	std::atomic<deviceFP>* pulseSumAtomic = (std::atomic<deviceFP>*)pulseSum;
	(*pulseSumAtomic).fetch_add(pointEnergy);
}
#endif
[[maybe_unused]]static int hardwareCheck(int* CUDAdeviceCount) {
	*CUDAdeviceCount = 1;
	return 0;
}
namespace deviceLibCPUFP32{
	inline float exp(const float x){
		return expf(x);
	}
	inline float abs(const float x){
		return fabs(x);
	}
	inline float sin(const float x){
		return sinf(x);
	}
	inline float cos(const float x){
		return cosf(x);
	}
	inline float atan(const float x){
		return atanf(x);
	}
	inline float sqrt(const float x){
		return sqrtf(x);
	}
	inline float asin(const float x){
		return asinf(x);
	}
	inline float pow(const float x, const float y){
		return powf(x,y);
	}
	inline float atan2(const float x, const float y) {
		return atan2f(x, y);
	}
	inline float acos(const float x) {
		return acosf(x);
	}
	inline std::complex<float> pow(const std::complex<float> x, const float y){
		return std::pow(x,y);
	}
	inline std::complex<float> exp(const std::complex<float> x){
		return std::exp(x);
	}
	inline float abs(const std::complex<float> x){
		return std::abs(x);
	}
	inline std::complex<float> sqrt(const std::complex<float> x){
		return std::sqrt(x);
	}
};

namespace complexLib = std;
#if LWEFLOATINGPOINT==32
namespace deviceFPLib = deviceLibCPUFP32;
namespace deviceLib = deviceLibCPUFP32;
#else
namespace deviceFPLib = std;
namespace deviceLib = std;
#endif

[[maybe_unused]] static std::complex<double> operator+(const float f, const std::complex<double> x) { return std::complex<double>(x.real() + f, x.imag()); }

[[maybe_unused]] static std::complex<double> operator+(const std::complex<double> x, const float f) { return std::complex<double>(x.real() + f, x.imag()); }

[[maybe_unused]] static std::complex<double> operator-(const std::complex<double> x, const float f) { return std::complex<double>(x.real() - f, x.imag()); }

[[maybe_unused]] static std::complex<double> operator*(const float f, const std::complex<double> x) { return std::complex<double>(x.real() * f, x.imag() * f); }

[[maybe_unused]] static std::complex<double> operator*(const std::complex<double> x, const float f) { return std::complex<double>(x.real() * f, x.imag() * f); }

[[maybe_unused]] static std::complex<double> operator/(const std::complex<double> x, const float f) { return std::complex<double>(x.real() / f, x.imag() / f); }

[[maybe_unused]] static double j0Device(double x) {
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

[[maybe_unused]] static float j0Device(float x) {
	if (x < 8.0f) {
		float y = x * x;
		float ans1 = 57568490574.0f + y * (-13362590354.0f + y * (651619640.7f +
			y * (-11214424.18f + y * (77392.33017f + y * (-184.9052456f)))));
		float ans2 = 57568490411.0f + y * (1029532985.0f + y * (9494680.718f +
			y * (59272.64853f + y * (267.8532712f + y))));
		return ans1 / ans2;
	}
	else {
		float z = 8.0f / x;
		float y = z * z;
		float xx = x - 0.785398164f;
		float ans1 = 1.0f + y * (-0.1098628627e-2f + y * (0.2734510407e-4f +
			y * (-0.2073370639e-5f + y * 0.2093887211e-6f)));
		float ans2 = -0.1562499995e-1f + y * (0.1430488765e-3f +
			y * (-0.6911147651e-5f + y * (0.7621095161e-6f - y * 0.934935152e-7f)));
		return sqrt(0.636619772f / x) * (cos(xx) * ans1 - z * sin(xx) * ans2);
	}
}

template <typename deviceFP, typename deviceComplex>
class CPUDevice {
private:
	bool configuredFFT = false;
	bool isCylindric = false;
	deviceFP canaryPixel = 0.0;
	deviceParameterSet<deviceFP, deviceComplex> dParamslocal;
	fftw_plan fftPlanD2Z;
	fftw_plan fftPlanZ2D;
	fftw_plan fftPlan1DD2Z;
	fftw_plan fftPlan1DZ2D;
	fftw_plan doublePolfftPlan;
	fftwf_plan fftPlanD2Z32;
	fftwf_plan fftPlanZ2D32;
	fftwf_plan fftPlan1DD2Z32;
	fftwf_plan fftPlan1DZ2D32;
	fftwf_plan doublePolfftPlan32;

	void fftDestroy() {
		if(sizeof(deviceFP)==sizeof(double)){
			fftw_destroy_plan(fftPlanD2Z);
			fftw_destroy_plan(fftPlanZ2D);
			fftw_destroy_plan(fftPlan1DD2Z);
			fftw_destroy_plan(fftPlan1DZ2D);
			if (isCylindric)fftw_destroy_plan(doublePolfftPlan);
		}
		else{
			fftwf_destroy_plan(fftPlanD2Z32);
			fftwf_destroy_plan(fftPlanZ2D32);
			fftwf_destroy_plan(fftPlan1DD2Z32);
			fftwf_destroy_plan(fftPlan1DZ2D32);
			if (isCylindric)fftwf_destroy_plan(doublePolfftPlan32);
		}
		fftw_cleanup();
	}

public:
	int stream;
	int memoryStatus;
	bool hasPlasma = false;
	deviceParameterSet<deviceFP, deviceComplex> deviceStruct;
	deviceParameterSet<deviceFP, deviceComplex>* s;
	simulationParameterSet* cParams;
	deviceParameterSet<deviceFP, deviceComplex>* dParamsDevice;

	CPUDevice(simulationParameterSet* sCPU) {
		s = &deviceStruct;
		memoryStatus = -1;
		stream = 0;
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
		memoryStatus = allocateSet(sCPU);
	}

	~CPUDevice() {
		fftDestroy();
		deallocateSet();
	}
//	template<typename Function, typename... Args>
//	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, const Function& kernel, Args... args) const {
//#pragma omp parallel for
//		for (int i = 0; i < (int)Nthread; ++i) {
//			for (unsigned int j = 0u; j < Nblock; ++j) {
//				kernel(j + Nblock * (unsigned int)i, args...);
//			}
//		}
//	}

	template<typename T>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, const T& functor) const {
		std::vector<size_t> vec(Nblock*Nthread);
		std::iota(vec.begin(), vec.end(), 0);
		std::for_each(std::execution::par_unseq, vec.begin(), vec.end(), functor);
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize){
		(*ptr) = calloc(N, elementSize);
		return (int)((*ptr) == NULL);
	}

	void deviceMemset(void* ptr, int value, size_t count){
		memset(ptr, value, count);
	}

	template <typename T>
	void deviceMemcpy(T* dst, const T* src, size_t count, copyType kind) {
		std::memcpy(dst, src, count);
	}

	void deviceMemcpy(double* dst, const float* src, size_t count, copyType kind) {
		for (size_t i = 0; i < count / sizeof(double); i++) {
			dst[i] = (double)src[i];
		}
	}

	void deviceMemcpy(std::complex<double>* dst, const std::complex<float>* src, size_t count, copyType kind) {
		for (size_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			dst[i] = std::complex<double>((double)src[i].real(), (double)src[i].imag());
		}
	}

	void deviceMemcpy(std::complex<float>* dst, const std::complex<double>* src, size_t count, copyType kind) {
		for (size_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			dst[i] = std::complex<float>((float)src[i].real(), (float)src[i].imag());
		}
	}

	void deviceMemcpy(float* dst, const double* src, size_t count, copyType kind) {
		for (size_t i = 0; i < count / sizeof(double); i++) {
			dst[i] = (float)src[i];
		}
	}

	void deviceFree(void* block){
		free(block);
	}

	inline bool isTheCanaryPixelNaN(const deviceFP* canaryPointer) {
		return(isnan(*canaryPointer));
	}
	void fft(void* input, void* output, deviceFFT type) const {
		if (!configuredFFT) return;
		switch (static_cast<int>(type) + 5 * (sizeof(deviceFP) == sizeof(float))){
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
		case 5:
			fftwf_execute_dft_r2c(fftPlanD2Z32, (float*)input, (fftwf_complex*)output);
			break;
		case 6:
			fftwf_execute_dft_c2r(fftPlanZ2D32, (fftwf_complex*)input, (float*)output);
			break;
		case 7:
			fftwf_execute_dft_r2c(fftPlan1DD2Z32, (float*)input, (fftwf_complex*)output);
			break;
		case 8:
			fftwf_execute_dft_c2r(fftPlan1DZ2D32, (fftwf_complex*)input, (float*)output);
			break;
		case 9:
			fftwf_execute_dft_r2c(doublePolfftPlan32, (float*)input, (fftwf_complex*)output);
			break;
		}
	}

	void fftInitialize() {
		if (configuredFFT) fftDestroy();
		isCylindric = (*s).isCylindric;
		hasPlasma = (*s).hasPlasma;
		if(sizeof(deviceFP)==sizeof(double)){
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
		}
		else if(sizeof(deviceFP)==sizeof(float)){
			float* setupWorkD = new float[(*s).Ngrid * (2 + 2 * (*s).isCylindric)];
			std::complex<float>* setupWorkC = new std::complex<float>[(*s).NgridC * (2 + 2 * (*s).isCylindric)];
			const int fftw1[1] = { (int)(*s).Ntime };
			fftPlan1DD2Z32 = fftwf_plan_many_dft_r2c(1, fftw1, (int)(*s).Nspace * (int)(*s).Nspace2 * 2, setupWorkD, fftw1, 1, (int)(*s).Ntime, (fftwf_complex*)setupWorkC, fftw1, 1, (int)(*s).Nfreq, FFTW_ESTIMATE);
			fftPlan1DZ2D32 = fftwf_plan_many_dft_c2r(1, fftw1, (int)(*s).Nspace * (int)(*s).Nspace2 * 2, (fftwf_complex*)setupWorkC, fftw1, 1, (int)(*s).Nfreq, setupWorkD, fftw1, 1, (int)(*s).Ntime, FFTW_ESTIMATE);
			if ((*s).is3D) {
				const int fftwSizes[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
				fftPlanD2Z32 = fftwf_plan_many_dft_r2c(3, fftwSizes, 2, setupWorkD, NULL, 1, (int)(*s).Ngrid, (fftwf_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, FFTW_MEASURE);
				fftPlanZ2D32 = fftwf_plan_many_dft_c2r(3, fftwSizes, 2, (fftwf_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, setupWorkD, NULL, 1, (int)(*s).Ngrid, FFTW_MEASURE);
			}
			else {
				const int fftwSizes[] = { (int)(*s).Nspace, (int)(*s).Ntime };
				fftPlanD2Z32 = fftwf_plan_many_dft_r2c(2, fftwSizes, 2, setupWorkD, NULL, 1, (int)(*s).Ngrid, (fftwf_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, FFTW_MEASURE);
				fftPlanZ2D32 = fftwf_plan_many_dft_c2r(2, fftwSizes, 2, (fftwf_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, setupWorkD, NULL, 1, (int)(*s).Ngrid, FFTW_MEASURE);

				if ((*s).isCylindric) {
					const int fftwSizesCyl[] = { (int)(2 * (*s).Nspace), (int)(*s).Ntime };
					doublePolfftPlan32 = fftwf_plan_many_dft_r2c(2, fftwSizesCyl, 2 + 2 * (*s).hasPlasma, setupWorkD, NULL, 1, 2 * (int)(*s).Ngrid, (fftwf_complex*)setupWorkC, NULL, 1, 2 * (int)(*s).NgridC, FFTW_MEASURE);
				}
			}
			delete[] setupWorkC;
			delete[] setupWorkD;
		}
		else{
			configuredFFT = 0;
			return;
		}
		configuredFFT = 1;
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
			sCPU->initializeDeviceParameters(s);
			hasPlasma = s->hasPlasma;
		}
		
		sCPU->fillRotationMatricies(s);
		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
			if ((*s).hasPlasma) beamExpansionFactor = 4;
		}
		deviceMemset((*s).gridETime1, 0, 2 * (*s).Ngrid * sizeof(deviceFP));
		deviceMemset((*s).gridPolarizationTime1, 0, 2 * (*s).Ngrid * sizeof(deviceFP));
		deviceMemset((*s).workspace1, 0, beamExpansionFactor * 2 * (*s).NgridC * sizeof(std::complex<deviceFP>));
		deviceMemset((*s).gridEFrequency1, 0, 2 * (*s).NgridC * sizeof(std::complex<deviceFP>));
		deviceMemset((*s).gridPropagationFactor1, 0, 2 * (*s).NgridC * sizeof(std::complex<deviceFP>));
		deviceMemset((*s).gridPolarizationFactor1, 0, 2 * (*s).NgridC * sizeof(std::complex<deviceFP>));
		deviceMemset((*s).gridEFrequency1Next1, 0, 2 * (*s).NgridC * sizeof(std::complex<deviceFP>));
		deviceMemset((*s).k1, 0, 2 * (*s).NgridC * sizeof(std::complex<deviceFP>));

		//cylindric sym grids
		if ((*s).isCylindric) {
			deviceMemset((*s).gridPropagationFactor1Rho1, 0, 4 * (*s).NgridC * sizeof(std::complex<deviceFP>));
			deviceMemset((*s).gridRadialLaplacian1, 0, 2 * beamExpansionFactor * (*s).Ngrid * sizeof(std::complex<deviceFP>));
		}

		//smaller helper grids
		deviceMemset((*s).expGammaT, 0, 2 * (*s).Ntime * sizeof(deviceFP));
		deviceMemset((*s).chiLinear1, 0, 2 * (*s).Nfreq * sizeof(std::complex<deviceFP>));
		deviceMemset((*s).fieldFactor1, 0, 2 * (*s).Nfreq * sizeof(deviceFP));
		deviceMemset((*s).inverseChiLinear1, 0, 2 * (*s).Nfreq * sizeof(deviceFP));

		double* expGammaTCPU = new double[2 * (*s).Ntime];
		for (size_t i = 0; i < (*s).Ntime; ++i) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, copyType::ToDevice);
		delete[] expGammaTCPU;

		sCPU->finishConfiguration(s);
		deviceMemcpy(dParamsDevice, s, sizeof(deviceParameterSet<deviceFP, deviceComplex>), copyType::ToDevice);
	}
	int allocateSet(simulationParameterSet* sCPU) {
		cParams = sCPU;
		sCPU->initializeDeviceParameters(s);
		hasPlasma = s->hasPlasma;
		fftInitialize();

		int memErrors = 0;
		double* expGammaTCPU = new double[2 * (*s).Ntime];

		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
			if ((*s).hasPlasma) beamExpansionFactor = 4;
		}
		sCPU->fillRotationMatricies(s);
		//GPU allocations
		//
		// currently 8 large grids, meaning memory use is approximately
		// 64 bytes per grid point (8 grids x 2 polarizations x 4ouble precision)
		// plus a little bit for additional constants/workspaces/etc
		memErrors += deviceCalloc((void**)&(*s).gridETime1, 2 * (*s).Ngrid, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationTime1, 2 * (*s).Ngrid, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)&(*s).workspace1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(std::complex<deviceFP>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1, 2 * (*s).NgridC, sizeof(std::complex<deviceFP>));
		memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1, 2 * (*s).NgridC, sizeof(std::complex<deviceFP>));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationFactor1, 2 * (*s).NgridC, sizeof(std::complex<deviceFP>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1Next1, 2 * (*s).NgridC, sizeof(std::complex<deviceFP>));
		memErrors += deviceCalloc((void**)&(*s).k1, 2 * (*s).NgridC, sizeof(std::complex<deviceFP>));

		//cylindric sym grids
		if ((*s).isCylindric) {
			memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1Rho1, 4 * (*s).NgridC, sizeof(std::complex<deviceFP>));
			memErrors += deviceCalloc((void**)&(*s).gridRadialLaplacian1, 2 * beamExpansionFactor * (*s).Ngrid, sizeof(std::complex<deviceFP>));
		}

		//smaller helper grids
		memErrors += deviceCalloc((void**)&(*s).expGammaT, 2 * (*s).Ntime, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(std::complex<deviceFP>));
		memErrors += deviceCalloc((void**)&(*s).fieldFactor1, 2 * (*s).Nfreq, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)&(*s).inverseChiLinear1, 2 * (*s).Nfreq, sizeof(deviceFP));
		for (size_t i = 0; i < (*s).Ntime; ++i) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, copyType::ToDevice);
		delete[] expGammaTCPU;
		(*sCPU).memoryError = memErrors;
		if (memErrors > 0) {
			return memErrors;
		}
		sCPU->finishConfiguration(s);
		deviceMemcpy(dParamsDevice, s, sizeof(deviceParameterSet<deviceFP, deviceComplex>), copyType::ToDevice);
		return 0;
	}
};