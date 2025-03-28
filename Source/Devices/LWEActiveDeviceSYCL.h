#include "../LightwaveExplorerUtilities.h"
#define SYCL_EXT_ONEAPI_COMPLEX                                        
#include <sycl/ext/oneapi/experimental/complex/complex.hpp>
#include <sycl/sycl.hpp>   
#include <sycl/atomic.hpp>
#include <oneapi/mkl/dft.hpp>
using std::isnan;

template <typename deviceFP>
static void atomicAdd(deviceFP* pulseSum, deviceFP pointEnergy) {
	sycl::atomic_ref<deviceFP, sycl::memory_order::relaxed, sycl::memory_scope::device> a(*pulseSum);
	a.fetch_add(pointEnergy);
}

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
	[[maybe_unused]] static sycl::ext::oneapi::experimental::complex<float> pow(const sycl::ext::oneapi::experimental::complex<float> x, const float y){
		float theta = atan2f(x.imag(), x.real());
		float rn = powf(x.real() * x.real() + x.imag() * x.imag(), 0.5f*y);
		return sycl::ext::oneapi::experimental::complex<float>(rn*cosf(y*theta),rn*sinf(y*theta));
	}
	static inline sycl::ext::oneapi::experimental::complex<float> exp(const sycl::ext::oneapi::experimental::complex<float> x){
		return sycl::ext::oneapi::experimental::complex<float>{exp(x.real())*cos(x.imag()),exp(x.real())*sin(x.imag())};
	}
	static inline float abs(const sycl::ext::oneapi::experimental::complex<float> x){
		return sqrtf(x.real() * x.real() + x.imag() * x.imag());
	}
	static inline sycl::ext::oneapi::experimental::complex<float> sqrt(const sycl::ext::oneapi::experimental::complex<float> x){
		float h = sqrtf(abs(x) + x.real());
		return 0.70710678118f * sycl::ext::oneapi::experimental::complex<float>{
			h,
			x.imag()/h
		};
	}
};

[[maybe_unused]] static int hardwareCheck(int* CUDAdeviceCount) {
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

template<typename deviceFP, typename deviceComplex>
class SYCLDevice : public LWEDevice {
private:
	deviceFP canaryPixel = 0.0;
	oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>* fftPlanD2Z;
	oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>* fftPlan1DD2Z;
	oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>* doublePolfftPlan;

	void fftDestroy() {
		delete fftPlan1DD2Z;
		delete fftPlanD2Z;
		if (s->isCylindric) delete doublePolfftPlan;
	}
public:
	sycl::queue stream;
	deviceParameterSet<deviceFP, deviceComplex>* s;
	deviceParameterSet<deviceFP, deviceComplex>* dParamsDevice;
	std::unique_ptr<UPPEAllocation<deviceFP, deviceComplex>> allocation;
	using LWEDevice::deviceMemcpy;
	SYCLDevice(simulationParameterSet* sCPU) {
		memoryStatus = allocateSet(sCPU);
		s = &(allocation->parameterSet);
	}

	~SYCLDevice() {
		stream.wait();
		fftDestroy();
	}

	bool isTheCanaryPixelNaN(const deviceFP* canaryPointer) {
		stream.memcpy(&canaryPixel, canaryPointer, sizeof(deviceFP));
		stream.wait();
		return(isnan(canaryPixel));
		return false;
	}

	template <typename T>
	void deviceLaunch(const unsigned int Nblock, const unsigned int Nthread, const T& functor) {
		int64_t i = static_cast<int64_t>(Nblock) * static_cast<int64_t>(Nthread);
		stream.submit([&](sycl::handler& h) {
			h.parallel_for(i, functor);
			});
}

	int deviceCalloc(void** ptr, const size_t N, const size_t elementSize) override {
		(*ptr) = sycl::aligned_alloc_device(
			2 * sizeof(deviceFP), 
			N * elementSize, 
			stream.get_device(), 
			stream.get_context());
		stream.memset((*ptr), 0, N * elementSize);
		stream.wait();
		return 0;
	}

	void deviceMemset(void* ptr, int value, size_t count) override {
		stream.wait();
		stream.memset(ptr, value, count);
	}

	void deviceMemcpyImplementation(
		void* dst, 
		const void* src, 
		const size_t count, 
		const copyType kind) override {
		stream.wait();
		stream.memcpy(dst, src, count);
		stream.wait();
	}

	void deviceMemcpy(
		std::complex<double>* dst, 
		const sycl::ext::oneapi::experimental::complex<float>* src, 
		size_t count, 
		copyType kind) {
		deviceMemcpy(dst, reinterpret_cast<const std::complex<float>*>(src),count,kind);
	}

	void deviceMemcpy(
		sycl::ext::oneapi::experimental::complex<float>* dst, 
		const std::complex<double>* src, 
		const size_t count, 
		const copyType kind) {
		deviceMemcpy(reinterpret_cast<std::complex<float>*>(dst),src,count,kind);
	}

	void deviceFree(void* block) override {
		stream.wait();
		sycl::free(block, stream);
	}

	void fftInitialize() {
		if (configuredFFT) {
			fftDestroy();
		}
		fftPlan1DD2Z = 
			new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>((*s).Ntime);
		fftPlan1DD2Z->set_value(
			oneapi::mkl::dft::config_param::PLACEMENT, oneapi::mkl::dft::config_value::NOT_INPLACE);
		std::int64_t outputStrides[3] = { 0, 1 };
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, outputStrides);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, outputStrides);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, (*s).Ntime);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, (*s).Nfreq);
		fftPlan1DD2Z->set_value(
			oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2 * ((*s).Nspace * (*s).Nspace2));

		if ((*s).is3D) {
			fftPlanD2Z = new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>
				(std::vector<std::int64_t>{(*s).Nspace2, (*s).Nspace, (*s).Ntime});
			fftPlanD2Z->set_value(
				oneapi::mkl::dft::config_param::PLACEMENT, oneapi::mkl::dft::config_value::NOT_INPLACE);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);

			std::int64_t forwardStride3D[5] = { 0, (*s).Nspace * (*s).Ntime, (*s).Nfreq, 1 };
			std::int64_t backwardStride3D[5] = { 0, (*s).Nspace * (*s).Nfreq, (*s).Nfreq, 1 };
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, forwardStride3D);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, backwardStride3D);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, (*s).Ntime * (*s).Nspace * (*s).Nspace2);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, (*s).Nfreq * (*s).Nspace * (*s).Nspace2);
		}
		else {
			fftPlanD2Z = new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>(
				std::vector<std::int64_t>{(*s).Nspace, (*s).Ntime});
			fftPlanD2Z->set_value(
				oneapi::mkl::dft::config_param::PLACEMENT, oneapi::mkl::dft::config_value::NOT_INPLACE);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
			std::int64_t backwardStride2D[4] = { 0, (*s).Nfreq, 1 };
			std::int64_t forwardStride2D[4] = { 0, (*s).Ntime, 1 };
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, forwardStride2D);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, backwardStride2D);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, (*s).Ntime * (*s).Nspace);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, (*s).Nfreq * (*s).Nspace);

			if ((*s).isCylindric) {
				doublePolfftPlan = 
					new oneapi::mkl::dft::descriptor<dftPrecision, oneapi::mkl::dft::domain::REAL>(
					std::vector<std::int64_t>{2 * (*s).Nspace, (*s).Ntime});
				doublePolfftPlan->set_value(
					oneapi::mkl::dft::config_param::PLACEMENT, oneapi::mkl::dft::config_value::NOT_INPLACE);
				doublePolfftPlan->set_value(
					oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2 + 2*(*s).hasPlasma);
				doublePolfftPlan->set_value(
					oneapi::mkl::dft::config_param::BWD_STRIDES, backwardStride2D);
					doublePolfftPlan->set_value(
					oneapi::mkl::dft::config_param::FWD_STRIDES, forwardStride2D);
				doublePolfftPlan->set_value(
					oneapi::mkl::dft::config_param::FWD_DISTANCE, (*s).Ntime * 2 * (*s).Nspace);
				doublePolfftPlan->set_value(
					oneapi::mkl::dft::config_param::BWD_DISTANCE, (*s).Nfreq * 2 * (*s).Nspace);
			}
		}

		fftPlan1DD2Z->commit(stream);
		fftPlanD2Z->commit(stream);
		if (s->isCylindric) doublePolfftPlan->commit(stream);
		configuredFFT = true;
	}

	void fft(const void* input, void* output, const deviceFFT type) override {
		switch (type) {
		case deviceFFT::D2Z:
			oneapi::mkl::dft::compute_forward(*fftPlanD2Z, (deviceFP*)input, (std::complex<deviceFP>*)output);
			break;
		case deviceFFT::Z2D:
			oneapi::mkl::dft::compute_backward(*fftPlanD2Z, (std::complex<deviceFP>*)input, (deviceFP*)output);
			break;
		case deviceFFT::D2Z_1D:
			oneapi::mkl::dft::compute_forward(*fftPlan1DD2Z, (deviceFP*)input, (std::complex<deviceFP>*)output);
			break;
		case deviceFFT::Z2D_1D:
			oneapi::mkl::dft::compute_backward(*fftPlan1DD2Z, (std::complex<deviceFP>*)input, (deviceFP*)output);
			break;
		case deviceFFT::D2Z_Polarization:
			oneapi::mkl::dft::compute_forward(*doublePolfftPlan, (deviceFP*)input, (std::complex<deviceFP>*)output);
			break;
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
			if (sizeof(deviceFP) == sizeof(double) 
				&& defaultStream.get_device().get_info<sycl::info::device::double_fp_config>().size()
				== 0) {
				sycl::queue cpuStream{ sycl::cpu_selector_v, sycl::property::queue::in_order() };
				stream = cpuStream;
			}
			else {
				stream = defaultStream;
			}
			
		}
		cParams = sCPU;
		if(memoryStatus == 0) allocation = nullptr;
		allocation = std::make_unique<UPPEAllocation<deviceFP, deviceComplex>>(this, sCPU);
		s = &(allocation->parameterSet);
		fftInitialize();
		dParamsDevice = allocation->parameterSet_deviceCopy.device_ptr();
		return 0;
	}
};