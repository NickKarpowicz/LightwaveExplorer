#include "../LightwaveExplorerUtilities.h"
#include "../LightwaveExplorerInterfaceClasses.hpp"
#ifdef USEFFTW
#include <fftw3.h>
#else
#include <fftw3_mkl.h>
#endif

#include <atomic>
#include <thread>
#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <execution>
#include <numeric>
#ifdef __APPLE__
#include <dispatch/dispatch.h>
#endif
static const int LWEThreadCount =
	maxN(static_cast<int>(std::thread::hardware_concurrency() / 2), 2);
#if defined __APPLE__ || defined __linux__
template<typename deviceFP>
[[maybe_unused]] void atomicAdd(deviceFP* pulseSum, deviceFP pointEnergy) {
	std::atomic<deviceFP>* pulseSumAtomic = (std::atomic<deviceFP>*)pulseSum;
	deviceFP expected = pulseSumAtomic->load();
	while (!std::atomic_compare_exchange_weak(pulseSumAtomic, &expected, expected + pointEnergy));
}
#ifdef __linux__
using std::isnan;
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
	inline float hypot(const float x, const float y) {
		return hypotf(x, y);
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

[[maybe_unused]] static std::complex<double> operator+(
	const float f,
	const std::complex<double> x) {
	return std::complex<double>(x.real() + f, x.imag());
}
[[maybe_unused]] static std::complex<double> operator+(
	const std::complex<double> x,
	const float f) {
	return std::complex<double>(x.real() + f, x.imag());
}
[[maybe_unused]] static std::complex<double> operator-(
	const std::complex<double> x,
	const float f) {
	return std::complex<double>(x.real() - f, x.imag()); }
[[maybe_unused]] static std::complex<double> operator*(
	const float f,
	const std::complex<double> x) {
	return std::complex<double>(x.real() * f, x.imag() * f);
}
[[maybe_unused]] static std::complex<double> operator*(
	const std::complex<double> x, const float f) {
	return std::complex<double>(x.real() * f, x.imag() * f);
}
[[maybe_unused]] static std::complex<double> operator/(
	const std::complex<double> x, const float f) {
	return std::complex<double>(x.real() / f, x.imag() / f);
}

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
class CPUDevice : public LWEDevice {
private:
#ifdef __APPLE__
	dispatch_queue_t queue{};
#elif defined _WIN32 || defined __linux__
	std::vector<int64_t> indices;
#endif
	fftw_plan fftPlanD2Z = {};
	fftw_plan fftPlanZ2D = {};
	fftw_plan fftPlan1DD2Z = {};
	fftw_plan fftPlan1DZ2D = {};
	fftw_plan doublePolfftPlan = {};
	fftwf_plan fftPlanD2Z32 = {};
	fftwf_plan fftPlanZ2D32 = {};
	fftwf_plan fftPlan1DD2Z32 = {};
	fftwf_plan fftPlan1DZ2D32 = {};
	fftwf_plan doublePolfftPlan32 = {};

	void fftDestroy() {
		if(sizeof(deviceFP)==sizeof(double)){
			fftw_destroy_plan(fftPlanD2Z);
			fftw_destroy_plan(fftPlanZ2D);
			fftw_destroy_plan(fftPlan1DD2Z);
			fftw_destroy_plan(fftPlan1DZ2D);
			if (s->isCylindric)fftw_destroy_plan(doublePolfftPlan);
		}
		else if(visualizationOnly){
			fftwf_destroy_plan(fftPlan1DD2Z32);
		}
		else{
			fftwf_destroy_plan(fftPlanD2Z32);
			fftwf_destroy_plan(fftPlanZ2D32);
			fftwf_destroy_plan(fftPlan1DD2Z32);
			fftwf_destroy_plan(fftPlan1DZ2D32);
			if (s->isCylindric)fftwf_destroy_plan(doublePolfftPlan32);
		}
		fftw_cleanup();
	}

public:
	deviceParameterSet<deviceFP, deviceComplex>* s;
	deviceParameterSet<deviceFP, deviceComplex>* dParamsDevice;
	std::unique_ptr<UPPEAllocation<deviceFP, deviceComplex>> allocation;
	std::unique_ptr<VisualizationAllocation<deviceComplex>> visualization;
	using LWEDevice::deviceMemcpy;
	CPUDevice(simulationParameterSet* sCPU) {
		memoryStatus = allocateSet(sCPU);
		s = &(allocation->parameterSet);
	#ifdef __APPLE__
		queue = dispatch_queue_create("Kernel", DISPATCH_QUEUE_CONCURRENT);
	#endif
	}

	CPUDevice(int64_t width, int64_t height, simulationParameterSet* sCPU){
		cParams = sCPU;
		visualizationOnly = true;
		visualization = std::make_unique<VisualizationAllocation<deviceComplex>>(this, width, height, sCPU);
		s = &(visualization->parameterSet);
		fftInitialize();
		#if defined _WIN32 || defined __linux__
            int64_t max_launch_size = maxN(
    		    maxN(
    				4 * (*s).NgridC,
    				(*sCPU).Nspace * (*sCPU).Nspace2 * (*sCPU).Npropagation),
       			    height * width
    		);
			indices = std::vector<int64_t>(max_launch_size);
			std::iota(indices.begin(), indices.end(), 0);
		#endif
		#ifdef __APPLE__
			queue = dispatch_queue_create("Kernel", DISPATCH_QUEUE_CONCURRENT);
		#endif
	}

	~CPUDevice() {
		fftDestroy();
	}

	//Use different threading models on different platforms:
	// MacOS: Grand Central Dispatch
	// Window/Linux: C++ std::execution::par or OpenMP
	template<typename T>
#ifdef __APPLE__
	void deviceLaunch(
		const unsigned int Nblock,
		const unsigned int Nthread,
		const T& functor) const {
			for (auto i = 0; i < Nthread; i++) {
				dispatch_async(queue, ^{
					const int64_t offset = i * static_cast<int64_t>(Nblock);
					for(int64_t j = offset; j < static_cast<int64_t>(offset+Nblock); functor(j++)){}
				});
			}
			dispatch_barrier_sync(queue,^{});
		}
#else
	void deviceLaunch(
		const unsigned int Nblock,
		const unsigned int Nthread,
		const T& functor) const {
		if (cParams->useOpenMP) {
#pragma omp parallel for num_threads(LWEThreadCount)
			for (int i = 0; i < static_cast<int>(Nblock); i++) {
				const int64_t offset = i * static_cast<int64_t>(Nthread);
				for (int64_t j = offset; j < static_cast<int64_t>(offset + Nthread); functor(j++)) {}
			}
		}
		else {
#if defined _WIN32 || defined __linux__ && not defined CPUONLY
			std::for_each(
				std::execution::par,
				indices.begin(),
				indices.begin() + static_cast<int64_t>(Nblock) * Nthread,
				functor);
#endif
		}
	}
#endif
	int deviceCalloc(void** ptr, const size_t N, const size_t elementSize) override {
		*ptr = calloc(N, elementSize);
		return static_cast<int>(*ptr == nullptr);
	}
	void deviceMemset(void* ptr, int value, size_t count) override {
		memset(ptr, value, count);
	}
	void deviceMemcpyImplementation(void* dst, const void* src, size_t count, copyType kind) override {
		std::memcpy(dst, src, count);
	}

	void deviceFree(void* block) override {
		free(block);
	}

	inline bool isTheCanaryPixelNaN(const deviceFP* canaryPointer) {
		return isnan(*canaryPointer);
	}

	void fft(const void* input, void* output, deviceFFT type) override {
		if (!configuredFFT) return;
		if(sizeof(deviceFP) == sizeof(double)){
			switch (type){
			case deviceFFT::D2Z:
				fftw_execute_dft_r2c(fftPlanD2Z, (double*)input, (fftw_complex*)output);
				break;
			case deviceFFT::Z2D:
				fftw_execute_dft_c2r(fftPlanZ2D, (fftw_complex*)input, (double*)output);
				break;
			case deviceFFT::D2Z_1D:
				fftw_execute_dft_r2c(fftPlan1DD2Z, (double*)input, (fftw_complex*)output);
				break;
			case deviceFFT::Z2D_1D:
				fftw_execute_dft_c2r(fftPlan1DZ2D, (fftw_complex*)input, (double*)output);
				break;
			case deviceFFT::D2Z_Polarization:
				fftw_execute_dft_r2c(doublePolfftPlan, (double*)input, (fftw_complex*)output);
				break;
			}
		}
		else{
			switch(type){
			case deviceFFT::D2Z:
				fftwf_execute_dft_r2c(fftPlanD2Z32, (float*)input, (fftwf_complex*)output);
				break;
			case deviceFFT::Z2D:
				fftwf_execute_dft_c2r(fftPlanZ2D32, (fftwf_complex*)input, (float*)output);
				break;
			case deviceFFT::D2Z_1D:
				fftwf_execute_dft_r2c(fftPlan1DD2Z32, (float*)input, (fftwf_complex*)output);
				break;
			case deviceFFT::Z2D_1D:
				fftwf_execute_dft_c2r(fftPlan1DZ2D32, (fftwf_complex*)input, (float*)output);
				break;
			case deviceFFT::D2Z_Polarization:
				fftwf_execute_dft_r2c(doublePolfftPlan32, (float*)input, (fftwf_complex*)output);
				break;
			}
		}
	}

	void fftInitialize() {
		if (configuredFFT) fftDestroy();
		if (sizeof(deviceFP) == sizeof(double)) {
			std::vector<double> setupWorkD((*s).Ngrid * (2 + 2 * (*s).isCylindric));
			std::vector<std::complex<double>> setupWorkC((*s).NgridC * (2 + 2 * (*s).isCylindric));
			const int fftw1[1] = { (int)(*s).Ntime };
			fftPlan1DD2Z = fftw_plan_many_dft_r2c(
				1,
				fftw1,
				(int)(*s).Nspace * (int)(*s).Nspace2 * 2,
				setupWorkD.data(),
				fftw1,
				1,
				(int)(*s).Ntime, (fftw_complex*)setupWorkC.data(),
				fftw1,
				1,
				(int)(*s).Nfreq,
				FFTW_ESTIMATE);
			fftPlan1DZ2D = fftw_plan_many_dft_c2r(
				1,
				fftw1,
				(int)(*s).Nspace * (int)(*s).Nspace2 * 2,
				(fftw_complex*)setupWorkC.data(),
				fftw1,
				1,
				(int)(*s).Nfreq,
				setupWorkD.data(),
				fftw1,
				1,
				(int)(*s).Ntime,
				FFTW_ESTIMATE);
			if ((*s).is3D) {
				const int fftwSizes[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
				fftPlanD2Z = fftw_plan_many_dft_r2c(
					3,
					fftwSizes,
					2,
					setupWorkD.data(),
					NULL,
					1,
					(int)(*s).Ngrid,
					(fftw_complex*)setupWorkC.data(),
					NULL, 1,
					(int)(*s).NgridC,
					FFTW_MEASURE);
				fftPlanZ2D = fftw_plan_many_dft_c2r(
					3,
					fftwSizes,
					2,
					(fftw_complex*)setupWorkC.data(),
					NULL,
					1,
					(int)(*s).NgridC,
					setupWorkD.data(),
					NULL,
					1,
					(int)(*s).Ngrid,
					FFTW_MEASURE);
			}
			else {
				const int fftwSizes[] = { (int)(*s).Nspace, (int)(*s).Ntime };
				fftPlanD2Z = fftw_plan_many_dft_r2c(
					2,
					fftwSizes,
					2,
					setupWorkD.data(),
					NULL,
					1,
					(int)(*s).Ngrid,
					(fftw_complex*)setupWorkC.data(),
					NULL,
					1,
					(int)(*s).NgridC,
					FFTW_MEASURE);
				fftPlanZ2D = fftw_plan_many_dft_c2r(
					2,
					fftwSizes,
					2,
					(fftw_complex*)setupWorkC.data(),
					NULL,
					1,
					(int)(*s).NgridC,
					setupWorkD.data(),
					NULL,
					1,
					(int)(*s).Ngrid,
					FFTW_MEASURE);

				if ((*s).isCylindric) {
					const int fftwSizesCyl[] = { (int)(2 * (*s).Nspace), (int)(*s).Ntime };
					doublePolfftPlan = fftw_plan_many_dft_r2c(
						2,
						fftwSizesCyl,
						2 + 2 * (*s).hasPlasma,
						setupWorkD.data(),
						NULL,
						1,
						2 * (int)(*s).Ngrid,
						(fftw_complex*)setupWorkC.data(),
						NULL,
						1,
						2 * (int)(*s).NgridC,
						FFTW_MEASURE);
				}
			}
		}
		else if(sizeof(deviceFP)==sizeof(float)){
			std::vector<float> setupWorkD((*s).Ngrid * (2 + 2 * (*s).isCylindric));
			std::vector<std::complex<float>> setupWorkC((*s).NgridC * (2 + 2 * (*s).isCylindric));
			const int fftw1[1] = { (int)(*s).Ntime };
			fftPlan1DD2Z32 = fftwf_plan_many_dft_r2c(
				1,
				fftw1,
				(int)(*s).Nspace * (int)(*s).Nspace2 * 2,
				setupWorkD.data(),
				fftw1,
				1,
				(int)(*s).Ntime,
				(fftwf_complex*)setupWorkC.data(),
				fftw1,
				1,
				(int)(*s).Nfreq,
				FFTW_ESTIMATE);
			if(!visualizationOnly)fftPlan1DZ2D32 = fftwf_plan_many_dft_c2r(
				1,
				fftw1,
				(int)(*s).Nspace * (int)(*s).Nspace2 * 2,
				(fftwf_complex*)setupWorkC.data(),
				fftw1,
				1,
				(int)(*s).Nfreq,
				setupWorkD.data(),
				fftw1,
				1,
				(int)(*s).Ntime,
				FFTW_ESTIMATE);
			if ((*s).is3D && !visualizationOnly) {
				const int fftwSizes[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
				fftPlanD2Z32 = fftwf_plan_many_dft_r2c(
					3,
					fftwSizes,
					2,
					setupWorkD.data(),
					NULL,
					1,
					(int)(*s).Ngrid,
					(fftwf_complex*)setupWorkC.data(),
					NULL,
					1,
					(int)(*s).NgridC,
					FFTW_MEASURE);
				fftPlanZ2D32 = fftwf_plan_many_dft_c2r(
					3,
					fftwSizes,
					2,
					(fftwf_complex*)setupWorkC.data(),
					NULL,
					1,
					(int)(*s).NgridC,
					setupWorkD.data(),
					NULL,
					1,
					(int)(*s).Ngrid,
					FFTW_MEASURE);
			}
			else {
				const int fftwSizes[] = { (int)(*s).Nspace, (int)(*s).Ntime };
				if(!visualizationOnly)fftPlanD2Z32 = fftwf_plan_many_dft_r2c(
					2,
					fftwSizes,
					2,
					setupWorkD.data(),
					NULL,
					1,
					(int)(*s).Ngrid,
					(fftwf_complex*)setupWorkC.data(),
					NULL,
					1,
					(int)(*s).NgridC,
					FFTW_MEASURE);
				if(!visualizationOnly)fftPlanZ2D32 = fftwf_plan_many_dft_c2r(
					2,
					fftwSizes,
					2,
					(fftwf_complex*)setupWorkC.data(),
					NULL,
					1,
					(int)(*s).NgridC,
					setupWorkD.data(),
					NULL,
					1,
					(int)(*s).Ngrid,
					FFTW_MEASURE);

				if ((*s).isCylindric && !visualizationOnly) {
					const int fftwSizesCyl[] = { (int)(2 * (*s).Nspace), (int)(*s).Ntime };
					doublePolfftPlan32 = fftwf_plan_many_dft_r2c(
						2,
						fftwSizesCyl,
						2 + 2 * (*s).hasPlasma,
						setupWorkD.data(),
						NULL,
						1,
						2 * (int)(*s).Ngrid,
						(fftwf_complex*)setupWorkC.data(),
						NULL,
						1,
						2 * (int)(*s).NgridC,
						FFTW_MEASURE);
				}
			}
		}
		else{
			configuredFFT = false;
			return;
		}
		configuredFFT = true;
	}

	void reset(simulationParameterSet* sCPU) override {
		bool resetFFT = (s->hasPlasma != sCPU->hasPlasma());
		cParams = sCPU;
		if(visualizationOnly){
			resetFFT = visualization->Nx != sCPU->Nspace
			|| visualization->Ny != sCPU->Nspace2
			|| visualization->Nt != sCPU->Ntime;
			visualization->setSimulationDimensions(sCPU);
		}
		else{
			allocation->useNewParameterSet(sCPU);
		}
		if(resetFFT){
			fftInitialize();
		}
	#if defined _WIN32 || defined __linux__
		if(indices.size() < static_cast<std::size_t>(4 * (*s).NgridC) ||
		(visualizationOnly && indices.size() < (static_cast<std::size_t>((*sCPU).Nspace * maxN((*sCPU).Nspace, (*sCPU).Nspace2))))){
			int64_t max_launch_size = maxN(
			    maxN(
					4 * (*s).NgridC,
					(*sCPU).Nspace * (*sCPU).Nspace2 * (*sCPU).Npropagation),
				visualizationOnly * (*sCPU).Nspace * maxN((*sCPU).Nspace, (*sCPU).Nspace2)
			);
		    indices = std::vector<int64_t>(max_launch_size);
			std::iota(indices.begin(), indices.end(), 0);
		}
	#endif
	}

	int allocateSet(simulationParameterSet* sCPU) {
		cParams = sCPU;
		if(memoryStatus == 0) allocation = nullptr;
		allocation = std::make_unique<UPPEAllocation<deviceFP, deviceComplex>>(this, sCPU);
		s = &(allocation->parameterSet);
		fftInitialize();
		dParamsDevice = allocation->parameterSet_deviceCopy.device_ptr();
		#if defined _WIN32 || defined __linux__
            int64_t max_launch_size = maxN(
    		    maxN(
    				4 * (*s).NgridC,
    				(*sCPU).Nspace * (*sCPU).Nspace2 * (*sCPU).Npropagation),
    			visualizationOnly * (*sCPU).Nspace * maxN((*sCPU).Nspace, (*sCPU).Nspace2)
    		);
			indices = std::vector<int64_t>(max_launch_size);
			std::iota(indices.begin(), indices.end(), 0);
		#endif
		return 0;
	}
};
