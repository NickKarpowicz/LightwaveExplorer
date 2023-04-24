#include "LightwaveExplorerUtilities.h"
#include <sycl/sycl.hpp>
#include <sycl/atomic.hpp>
#include <oneapi/mkl/dfti.hpp>
#include <oneapi/dpl/complex>
using std::isnan;

template <typename deviceFP>
static void atomicAdd(deviceFP* pulseSum, deviceFP pointEnergy) {
	sycl::atomic_ref<deviceFP, sycl::memory_order::relaxed, sycl::memory_scope::device> a(*pulseSum);
	a.fetch_add(pointEnergy);
}

[[maybe_unused]] static constexpr oneapi::dpl::complex<double> operator+(const float f, const oneapi::dpl::complex<double> x) { return oneapi::dpl::complex<double>(x.real() + f, x.imag()); }
[[maybe_unused]] static constexpr oneapi::dpl::complex<double> operator+(const oneapi::dpl::complex<double> x, const float f) { return oneapi::dpl::complex<double>(x.real() + f, x.imag()); }
[[maybe_unused]] static constexpr oneapi::dpl::complex<double> operator-(const oneapi::dpl::complex<double> x, const float f) { return oneapi::dpl::complex<double>(x.real() - f, x.imag()); }
[[maybe_unused]] static constexpr oneapi::dpl::complex<double> operator*(const float f, const oneapi::dpl::complex<double> x) { return oneapi::dpl::complex<double>(x.real() * f, x.imag() * f); }
[[maybe_unused]] static constexpr oneapi::dpl::complex<double> operator*(const oneapi::dpl::complex<double> x, const float f) { return oneapi::dpl::complex<double>(x.real() * f, x.imag() * f); }
[[maybe_unused]] static constexpr oneapi::dpl::complex<double> operator/(const oneapi::dpl::complex<double> x, const float f) { return oneapi::dpl::complex<double>(x.real() / f, x.imag() / f); }

namespace deviceLibSYCLFP32{
	static constexpr inline float exp(const float x){
		return expf(x);
	}
	static constexpr inline float abs(const float x){
		return fabs(x);
	}
	static constexpr inline float sin(const float x){
		return sinf(x);
	}
	static constexpr inline float cos(const float x){
		return cosf(x);
	}
	static constexpr inline float atan(const float x){
		return atanf(x);
	}
	static constexpr inline float sqrt(const float x){
		return sqrtf(x);
	}
	static constexpr inline float asin(const float x){
		return asinf(x);
	}
	static constexpr inline float pow(const float x, const float y){
		return powf(x,y);
	}
	static constexpr inline float acos(const float x) {
		return acosf(x);
	}
	static constexpr inline float atan2(const float x, const float y) {
		return atan2f(x, y);
	}
	static inline float hypot(const float x, const float y) {
		return hypotf(x, y);
	}
	static oneapi::dpl::complex<float> pow(const oneapi::dpl::complex<float> x, const float y){
		float r = sqrtf(x.real() * x.real() + x.imag() * x.imag());
		float theta = atan2f(x.imag(), x.real());
		float rn = powf(r, y);
		return oneapi::dpl::complex<float>(rn*cosf(y*theta),rn*sinf(y*theta));
	}
	static inline oneapi::dpl::complex<float> exp(const oneapi::dpl::complex<float> x){
		return oneapi::dpl::exp(x);
	}
	static inline float abs(const oneapi::dpl::complex<float> x){
		return oneapi::dpl::abs(x);
	}
	static inline oneapi::dpl::complex<float> sqrt(const oneapi::dpl::complex<float> x){
		return oneapi::dpl::sqrt(x);
	}
};


static int hardwareCheck(int* CUDAdeviceCount) {
	*CUDAdeviceCount = 1;
	return 0;
}
static constexpr double j0Device(const double x) {
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

static constexpr float j0Device(const float x) {
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
		return sqrtf(0.636619772f / x) * (cosf(xx) * ans1 - z * sinf(xx) * ans2);
	}
}

static constexpr oneapi::dpl::complex<double> operator/(const double a, const oneapi::dpl::complex<double> b) {
	double divByDenominator = a / (b.real() * b.real() + b.imag() * b.imag());
	return oneapi::dpl::complex<double>(b.real() * divByDenominator, -b.imag() * divByDenominator);
}

template<typename deviceFP, typename deviceComplex>
class SYCLDevice {
private:
	bool configuredFFT = false;
	bool isCylindric = false;
	deviceFP canaryPixel = 0.0;
	oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>* fftPlanD2Z;
	oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>* fftPlanZ2D;
	oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>* fftPlan1DD2Z;
	oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>* fftPlan1DZ2D;
	oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>* doublePolfftPlan;

	void fftDestroy() {
		delete fftPlan1DD2Z;
		delete fftPlan1DZ2D;
		delete fftPlanD2Z;
		delete fftPlanZ2D;
		if ((*s).isCylindric) delete doublePolfftPlan;
	}
public:
	sycl::queue stream;
	deviceParameterSet<deviceFP, deviceComplex> deviceStruct;
	deviceParameterSet<deviceFP, deviceComplex>* s;
	deviceParameterSet<deviceFP, deviceComplex>* dParamsDevice;
	simulationParameterSet* cParams;
	int memoryStatus;
	bool hasPlasma = false;

	SYCLDevice(simulationParameterSet* sCPU) {
		memoryStatus = -1;
		configuredFFT = 0;
		isCylindric = 0;
		s = &deviceStruct;
		memoryStatus = allocateSet(sCPU);
	}

	~SYCLDevice() {
		stream.wait();
		fftDestroy();
		deallocateSet();
		deviceFree(dParamsDevice);
	}

	bool isTheCanaryPixelNaN(const deviceFP* canaryPointer) {
		stream.memcpy(&canaryPixel, canaryPointer, sizeof(deviceFP));
		stream.wait();
		return(isnan(canaryPixel));
	}

	template <typename T>
	void deviceLaunch(const unsigned int Nblock, const unsigned int Nthread, const T& functor) {
		size_t i = static_cast<size_t>(Nblock) * static_cast<size_t>(Nthread);
		stream.submit([&](sycl::handler& h) {
			h.parallel_for(i, functor);
			});
}

	int deviceCalloc(void** ptr, const size_t N, const size_t elementSize) {
		(*ptr) = sycl::aligned_alloc_device(2 * sizeof(deviceFP), N * elementSize, stream.get_device(), stream.get_context());
		stream.memset((*ptr), 0, N * elementSize);
		stream.wait();
		return 0;
	}

	void deviceMemset(void* ptr, int value, size_t count) {
		stream.wait();
		stream.memset(ptr, value, count);
	}

	void deviceMemcpy(void* dst, const void* src, const size_t count, const copyType kind) {
		stream.wait();
		stream.memcpy(dst, src, count);
		stream.wait();
	}

	void deviceMemcpy(double* dst, const float* src, const size_t count, const copyType kind) {
		stream.wait();
		float* copyBuffer = new float[count / sizeof(double)]();
		stream.memcpy(copyBuffer, src, count/2);
		stream.wait();
		for (size_t i = 0; i < count / sizeof(double); i++) {
			dst[i] = static_cast<double>(copyBuffer[i]);
		}
		delete[] copyBuffer;
	}

	void deviceMemcpy(std::complex<double>* dst, const oneapi::dpl::complex<float>* src, size_t count, copyType kind) {
		stream.wait();
		std::complex<float>* copyBuffer = new std::complex<float>[count / sizeof(std::complex<double>)]();
		stream.memcpy(copyBuffer, src, count/2);
		stream.wait();
		for (size_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			dst[i] = std::complex<double>(
				static_cast<float>(copyBuffer[i].real()), 
				static_cast<float>(copyBuffer[i].imag()));
		}
		delete[] copyBuffer;
	}

	void deviceMemcpy(oneapi::dpl::complex<float>* dst, const std::complex<double>* src, const size_t count, const copyType kind) {
		stream.wait();
		std::complex<float>* copyBuffer = new std::complex<float>[count / sizeof(std::complex<double>)]();
		for (size_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			copyBuffer[i] = std::complex<float>(
				static_cast<float>(src[i].real()), 
				static_cast<float>(src[i].imag()));
		}
		stream.memcpy(dst, copyBuffer, count / 2);
		stream.wait();
		delete[] copyBuffer;
	}

	void deviceMemcpy(float* dst, const double* src, const size_t count, const copyType kind) {
		stream.wait();
		float* copyBuffer = new float[count / sizeof(double)]();
		for (size_t i = 0; i < count / sizeof(double); i++) {
			copyBuffer[i] = static_cast<float>(src[i]);
		}
		stream.memcpy(dst, copyBuffer, count / 2);
		stream.wait();
		delete[] copyBuffer;
	}

	void deviceFree(void* block) {
		stream.wait();
		sycl::free(block, stream);
	}

	void fftInitialize() {
		if (configuredFFT) {
			fftDestroy();
		}
		isCylindric = 0;
		hasPlasma = (*s).hasPlasma;
		size_t workSize;
		int Ntime = (*s).Ntime;
		int Nspace = (*s).Nspace;
		int Nspace2 = (*s).Nspace2;
		int Nfreq = (*s).Nfreq;

		fftPlan1DD2Z = new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>((int)Ntime);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
		std::int64_t outputStrides[2] = { 0, 1 };
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, outputStrides);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2 * (int)(Nspace * Nspace2));

		fftPlan1DZ2D = new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>(Ntime);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);

		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, outputStrides);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2 * (Nspace * Nspace2));

		if ((*s).is3D) {
			int cufftSizes1[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
			fftPlanD2Z = new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>
				(std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1], cufftSizes1[2]});
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);

			std::int64_t outputStride3D[4] = { 0, Nspace * Nfreq, Nfreq, 1 };
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, outputStride3D);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace * Nspace2);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace * Nspace2);

			fftPlanZ2D = new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>
				(std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1], cufftSizes1[2]});
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);

			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, outputStride3D);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace * Nspace2);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace * Nspace2);

		}
		else {
			int cufftSizes1[] = { (int)(*s).Nspace, (int)(*s).Ntime };

			fftPlanD2Z = new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>(
				std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1]});
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
			std::int64_t outputStride2D[3] = { 0, Nfreq, 1 };
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, outputStride2D);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace);

			fftPlanZ2D = new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>(
				std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1]});
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, outputStride2D);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace);

			if ((*s).isCylindric) {
				isCylindric = 1;
				int cufftSizes2[] = { 2 * (int)(*s).Nspace, (int)(*s).Ntime };
				doublePolfftPlan = new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>(
					std::vector<std::int64_t>{cufftSizes2[0], cufftSizes2[1]});
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2 + 2*(*s).hasPlasma);
				std::int64_t outputStrideCyl[3] = { 0, Nfreq, 1 };
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, outputStrideCyl);
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * 2 * Nspace);
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * 2 * Nspace);
			}
		}
		fftPlan1DD2Z->commit(stream);
		fftPlan1DZ2D->commit(stream);
		fftPlanD2Z->commit(stream);
		fftPlanZ2D->commit(stream);
		if ((*s).isCylindric) doublePolfftPlan->commit(stream);
		configuredFFT = 1;
	}

	void fft(const void* input, void* output, const deviceFFT type) const {
		switch (static_cast<int>(type)) {
		case 0:
			oneapi::mkl::dft::compute_forward(*fftPlanD2Z, (deviceFP*)input, (deviceComplex*)output);
			break;
		case 1:
			oneapi::mkl::dft::compute_backward(*fftPlanZ2D, (deviceComplex*)input, (deviceFP*)output);
			break;
		case 2:
			oneapi::mkl::dft::compute_forward(*fftPlan1DD2Z, (deviceFP*)input, (deviceComplex*)output);
			break;
		case 3:
			oneapi::mkl::dft::compute_backward(*fftPlan1DZ2D, (deviceComplex*)input, (deviceFP*)output);
			break;
		case 4:
			oneapi::mkl::dft::compute_forward(*doublePolfftPlan, (deviceFP*)input, (deviceComplex*)output);
			break;
		}
	}
	void deallocateSet() {
		stream.wait();
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
		deviceMemset((*s).workspace1, 0, beamExpansionFactor * 2 * (*s).NgridC * sizeof(deviceComplex));
		deviceMemset((*s).gridEFrequency1, 0, 2 * (*s).NgridC * sizeof(deviceComplex));
		deviceMemset((*s).gridPropagationFactor1, 0, 2 * (*s).NgridC * sizeof(deviceComplex));
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
		deviceMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, copyType::ToDevice);
		delete[] expGammaTCPU;

		sCPU->finishConfiguration(s);
		deviceMemcpy(dParamsDevice, s, sizeof(deviceParameterSet<deviceFP, deviceComplex>), copyType::ToDevice);
	}
	int allocateSet(simulationParameterSet* sCPU) {
		
		if ((*sCPU).assignedGPU) {
			try {
				sycl::queue gpuStream{ sycl::gpu_selector_v, { sycl::property::queue::in_order() } };
				stream = gpuStream;
			}
			catch (sycl::exception const& e) {
				sycl::queue defaultStream{ sycl::default_selector_v, sycl::property::queue::in_order() };
				stream = defaultStream;
			}
		}
		else if ((*sCPU).runningOnCPU) {
			sycl::queue cpuStream{ sycl::cpu_selector_v, sycl::property::queue::in_order() };
			stream = cpuStream;
		}
		else {
			sycl::queue defaultStream{ sycl::default_selector_v, sycl::property::queue::in_order() };
			if (sizeof(deviceFP) == sizeof(double) && defaultStream.get_device().get_info<cl::sycl::info::device::double_fp_config>().size() == 0) {
				sycl::queue cpuStream{ sycl::cpu_selector_v, sycl::property::queue::in_order() };
				stream = cpuStream;
			}
			else {
				stream = defaultStream;
			}
			
		}
		deviceCalloc((void**)&dParamsDevice, 1, sizeof(deviceParameterSet<deviceFP, deviceComplex>));
		
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
		memErrors += deviceCalloc((void**)&(*s).workspace1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationFactor1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1Next1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)&(*s).k1, 2 * (*s).NgridC, sizeof(deviceComplex));

		//cylindric sym grids
		if ((*s).isCylindric) {
			memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1Rho1, 4 * (*s).NgridC, sizeof(deviceComplex));
			memErrors += deviceCalloc((void**)&(*s).gridRadialLaplacian1, beamExpansionFactor * 2 * (*s).Ngrid, sizeof(deviceComplex));
		}

		//smaller helper grids
		memErrors += deviceCalloc((void**)&(*s).expGammaT, 2 * (*s).Ntime, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(deviceComplex));
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