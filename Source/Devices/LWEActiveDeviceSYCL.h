#include "../LightwaveExplorerUtilities.h"
#define SYCL_EXT_ONEAPI_COMPLEX  
#include <sycl/sycl.hpp> 
#include <dpct/dpct.hpp>
#include <dpct/fft_utils.hpp>                                      
#include <sycl/ext/oneapi/experimental/complex/complex.hpp>  
#include <sycl/atomic.hpp>
#include <oneapi/mkl/dfti.hpp>
using std::isnan;
#if LWEFLOATINGPOINT == 64
const auto MKLFFT_fwd = dpct::fft::fft_type::real_double_to_complex_double;
const auto MKLFFT_bwd = dpct::fft::fft_type::complex_double_to_real_double;
#else
const auto MKLFFT_fwd = dpct::fft::fft_type::real_float_to_complex_float;
const auto MKLFFT_bwd = dpct::fft::fft_type::complex_float_to_real_float;
#endif
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
	static sycl::ext::oneapi::experimental::complex<float> pow(const sycl::ext::oneapi::experimental::complex<float> x, const float y){
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

template<typename deviceFP, typename deviceComplex>
class SYCLDevice {
private:
	bool configuredFFT = false;
	bool isCylindric = false;
	deviceFP canaryPixel = 0.0;
	dpct::fft::fft_engine_ptr fftPlanD2Z;
    dpct::fft::fft_engine_ptr fftPlanZ2D;
    dpct::fft::fft_engine_ptr fftPlan1DD2Z;
    dpct::fft::fft_engine_ptr fftPlan1DZ2D;
    dpct::fft::fft_engine_ptr doublePolfftPlan;

	void fftDestroy() {
		stream.wait();
		dpct::fft::fft_engine::destroy(fftPlanD2Z);
		dpct::fft::fft_engine::destroy(fftPlanZ2D);
		dpct::fft::fft_engine::destroy(fftPlan1DD2Z);
		dpct::fft::fft_engine::destroy(fftPlan1DZ2D);
		if (isCylindric) {
				dpct::fft::fft_engine::destroy(doublePolfftPlan);
		}
		stream.wait();
		configuredFFT = 0;
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
		int64_t i = static_cast<int64_t>(Nblock) * static_cast<int64_t>(Nthread);
		stream.submit([&](sycl::handler& h) {
			h.parallel_for(i, functor);
			});
}

	int deviceCalloc(void** ptr, const int64_t N, const int64_t elementSize) {
		(*ptr) = sycl::aligned_alloc_device(
			2 * sizeof(deviceFP), 
			N * elementSize, 
			stream.get_device(), 
			stream.get_context());
		stream.memset((*ptr), 0, N * elementSize);
		stream.wait();
		return 0;
	}

	void deviceMemset(void* ptr, int value, int64_t count) {
		stream.wait();
		stream.memset(ptr, value, count);
	}

	void deviceMemcpy(
		void* dst, 
		const void* src, 
		const int64_t count, 
		const copyType kind) {
		stream.wait();
		stream.memcpy(dst, src, count);
		stream.wait();
	}

	void deviceMemcpy(
		double* dst, 
		const float* src, 
		const int64_t count, 
		const copyType kind) {
		stream.wait();
		float* copyBuffer = new float[count / sizeof(double)]();
		stream.memcpy(copyBuffer, src, count/2);
		stream.wait();
		for (int64_t i = 0; i < count / sizeof(double); i++) {
			dst[i] = static_cast<double>(copyBuffer[i]);
		}
		delete[] copyBuffer;
	}

	void deviceMemcpy(
		std::complex<double>* dst, 
		const sycl::ext::oneapi::experimental::complex<float>* src, 
		int64_t count, 
		copyType kind) {
		stream.wait();
		std::complex<float>* copyBuffer = new std::complex<float>[count / sizeof(std::complex<double>)]();
		stream.memcpy(copyBuffer, src, count/2);
		stream.wait();
		for (int64_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			dst[i] = std::complex<double>(
				static_cast<float>(copyBuffer[i].real()), 
				static_cast<float>(copyBuffer[i].imag()));
		}
		delete[] copyBuffer;
	}

	void deviceMemcpy(
		sycl::ext::oneapi::experimental::complex<float>* dst, 
		const std::complex<double>* src, 
		const int64_t count, 
		const copyType kind) {
		stream.wait();
		std::complex<float>* copyBuffer = new std::complex<float>[count / sizeof(std::complex<double>)]();
		for (int64_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			copyBuffer[i] = std::complex<float>(
				static_cast<float>(src[i].real()), 
				static_cast<float>(src[i].imag()));
		}
		stream.memcpy(dst, copyBuffer, count / 2);
		stream.wait();
		delete[] copyBuffer;
	}

	void deviceMemcpy(
		float* dst, 
		const double* src, 
		const int64_t count, 
		const copyType kind) {
		stream.wait();
		float* copyBuffer = new float[count / sizeof(double)]();
		for (int64_t i = 0; i < count / sizeof(double); i++) {
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
		fftPlan1DD2Z = dpct::fft::fft_engine::create(
			&dpct::get_in_order_queue(), (int)(*s).Ntime, MKLFFT_fwd,
			2 * (int)((*s).Nspace * (*s).Nspace2));
		fftPlan1DZ2D = dpct::fft::fft_engine::create(
			&dpct::get_in_order_queue(), (int)(*s).Ntime, MKLFFT_bwd,
			2 * (int)((*s).Nspace * (*s).Nspace2));
        if ((*s).is3D) {
			int cufftSizes1[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
            fftPlanD2Z = dpct::fft::fft_engine::create();
			dpct::fft::fft_engine::estimate_size(
				3, cufftSizes1, NULL, 0, 0, 0, 0, 0, MKLFFT_fwd, 2,
				&workSize);
			fftPlanD2Z->commit(&stream, 3,
								cufftSizes1, NULL, 0, 0, 0, 0, 0,
								MKLFFT_fwd, 2, &workSize);

			fftPlanZ2D = dpct::fft::fft_engine::create();
			dpct::fft::fft_engine::estimate_size(
				3, cufftSizes1, NULL, 0, 0, 0, 0, 0, MKLFFT_bwd, 2,
				&workSize);
			fftPlanZ2D->commit(&stream, 3,
								cufftSizes1, NULL, 0, 0, 0, 0, 0,
								MKLFFT_bwd, 2, &workSize);
        }
		else {
			int cufftSizes1[] = { (int)(*s).Nspace, (int)(*s).Ntime };

			fftPlanD2Z = dpct::fft::fft_engine::create();
			dpct::fft::fft_engine::estimate_size(
				2, cufftSizes1, NULL, 0, 0, 0, 0, 0, MKLFFT_fwd, 2,
				&workSize);
			fftPlanD2Z->commit(&stream, 2,
								cufftSizes1, NULL, 0, 0, 0, 0, 0,
								MKLFFT_fwd, 2, &workSize);

			fftPlanZ2D = dpct::fft::fft_engine::create();
			dpct::fft::fft_engine::estimate_size(
				2, cufftSizes1, NULL, 0, 0, 0, 0, 0, MKLFFT_bwd, 2,
				&workSize);
			fftPlanZ2D->commit(&stream, 2,
								cufftSizes1, NULL, 0, 0, 0, 0, 0,
								MKLFFT_bwd, 2, &workSize);

            if ((*s).isCylindric) {
				isCylindric = 1;
				int cufftSizes2[] = { 2 * (int)(*s).Nspace, (int)(*s).Ntime };
            	doublePolfftPlan = dpct::fft::fft_engine::create();
                dpct::fft::fft_engine::estimate_size(
					2, cufftSizes2, NULL, 
					0, 0, 0, 0, 0, MKLFFT_fwd, 2 + 2 * (*s).hasPlasma, &workSize);
				doublePolfftPlan->commit(
					&stream, 2, cufftSizes2, NULL, 
					0, 0, 0, 0, 0, MKLFFT_fwd, 2 + 2 * (*s).hasPlasma, &workSize);
            }
		}
        configuredFFT = 1;
	}

	void fft(const void* input, void* output, deviceFFT type) {
		stream.wait();
		if(sizeof(deviceFP) == sizeof(double)){
			switch (type){
			case deviceFFT::D2Z:
				fftPlanD2Z->compute<double, std::complex<double>>(
						(double *)input,
						(std::complex<double> *)output,
						dpct::fft::fft_direction::forward);
				break;
			case deviceFFT::Z2D:
				fftPlanZ2D->compute<std::complex<double>, double>(
						(std::complex<double> *)input,
						(double *)output,
						dpct::fft::fft_direction::backward);
				break;
			case deviceFFT::D2Z_1D:
				fftPlan1DD2Z->compute<
									double, std::complex<double>>(
					(double *)input, (std::complex<double> *)output,
					dpct::fft::fft_direction::forward);
				break;
			case deviceFFT::Z2D_1D:
				fftPlan1DZ2D->compute<
									std::complex<double>, double>(
					(std::complex<double> *)input, (double *)output,
					dpct::fft::fft_direction::backward);
				break;
			case deviceFFT::D2Z_Polarization:
                doublePolfftPlan->compute<double, std::complex<double>>(
                                    (double *)input, (std::complex<double> *)output,
                                    dpct::fft::fft_direction::forward);
                                break;
			}
		}
		else{
			switch(type){
			case deviceFFT::D2Z:
                fftPlanD2Z->compute<float, std::complex<float>>(
                                        (float *)input, (std::complex<float> *)output,
                                        dpct::fft::fft_direction::forward);
                                break;
			case deviceFFT::Z2D:
                fftPlanZ2D->compute<std::complex<float>, float>(
                                        (std::complex<float> *)input, (float *)output,
                                        dpct::fft::fft_direction::backward);
                                break;
			case deviceFFT::D2Z_1D:
				fftPlan1DD2Z->compute<float, std::complex<float>>(
                                        (float *)input, (std::complex<float> *)output,
                                        dpct::fft::fft_direction::forward);
                                break;
			case deviceFFT::Z2D_1D:
                fftPlan1DZ2D->compute<std::complex<float>, float>(
                                        (std::complex<float> *)input, (float *)output,
                                        dpct::fft::fft_direction::backward);
                                break;
			case deviceFFT::D2Z_Polarization:
                doublePolfftPlan->compute<float,std::complex<float>>(
                                        (float *)input, (std::complex<float> *)output,
                                        dpct::fft::fft_direction::forward);
                                break;
			}
		}
		stream.wait();
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
		int64_t beamExpansionFactor = 1;
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
			deviceMemset(
				(*s).gridPropagationFactor1Rho1, 
				0, 
				4 * (*s).NgridC * sizeof(deviceComplex));
			deviceMemset(
				(*s).gridRadialLaplacian1, 
				0, 
				2 * beamExpansionFactor * (*s).Ngrid * sizeof(deviceComplex));
		}

		//smaller helper grids
		deviceMemset((*s).expGammaT, 0, 2 * (*s).Ntime * sizeof(deviceFP));
		deviceMemset((*s).chiLinear1, 0, 2 * (*s).Nfreq * sizeof(deviceComplex));
		deviceMemset((*s).fieldFactor1, 0, 2 * (*s).Nfreq * sizeof(deviceFP));
		deviceMemset((*s).inverseChiLinear1, 0, 2 * (*s).Nfreq * sizeof(deviceFP));

		std::vector<double> expGammaTCPU(2 * (*s).Ntime);
		for (int64_t i = 0; i < (*s).Ntime; ++i) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU.data(), 2 * sizeof(double) * (*s).Ntime, copyType::ToDevice);

		sCPU->finishConfiguration(s);
		deviceMemcpy(
			dParamsDevice, 
			s, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);
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
		deviceCalloc((void**)&dParamsDevice, 1, sizeof(deviceParameterSet<deviceFP, deviceComplex>));
		
		cParams = sCPU;
		sCPU->initializeDeviceParameters(s);
		hasPlasma = s->hasPlasma;
		fftInitialize();
		int memErrors = 0;

		int64_t beamExpansionFactor = 1;
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
		memErrors += deviceCalloc((void**)
			&(*s).gridETime1, 2 * (*s).Ngrid, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)
			&(*s).gridPolarizationTime1, 2 * (*s).Ngrid, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)
			&(*s).workspace1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)
			&(*s).gridEFrequency1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)
			&(*s).gridPropagationFactor1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)
			&(*s).gridPolarizationFactor1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)
			&(*s).gridEFrequency1Next1, 2 * (*s).NgridC, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)
			&(*s).k1, 2 * (*s).NgridC, sizeof(deviceComplex));

		//cylindric sym grids
		if ((*s).isCylindric) {
			memErrors += deviceCalloc((void**)
				&(*s).gridPropagationFactor1Rho1, 4 * (*s).NgridC, sizeof(deviceComplex));
			memErrors += deviceCalloc((void**)
				&(*s).gridRadialLaplacian1, beamExpansionFactor * 2 * (*s).Ngrid, sizeof(deviceComplex));
		}

		//smaller helper grids
		memErrors += deviceCalloc((void**)
			&(*s).expGammaT, 2 * (*s).Ntime, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)
			&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(deviceComplex));
		memErrors += deviceCalloc((void**)
			&(*s).fieldFactor1, 2 * (*s).Nfreq, sizeof(deviceFP));
		memErrors += deviceCalloc((void**)
			&(*s).inverseChiLinear1, 2 * (*s).Nfreq, sizeof(deviceFP));
		std::vector<double> expGammaTCPU(2 * (*s).Ntime);
		for (int64_t i = 0; i < (*s).Ntime; ++i) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU.data(), 2 * sizeof(double) * (*s).Ntime, copyType::ToDevice);
		(*sCPU).memoryError = memErrors;
		if (memErrors > 0) {
			return memErrors;
		}
		sCPU->finishConfiguration(s);
		deviceMemcpy(
			dParamsDevice, 
			s, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);
		return 0;
	}
};