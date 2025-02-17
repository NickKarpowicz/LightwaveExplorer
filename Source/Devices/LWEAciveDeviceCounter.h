#pragma once
#include "../LightwaveExplorerUtilities.h"
namespace deviceLib = std;
namespace deviceFPLib = std;
namespace complexLib = std;
int hardwareCheck(int* CUDAdeviceCount) {
	*CUDAdeviceCount = 1;
	return 0;
}

[[maybe_unused]] static void atomicAdd(float* pulseSum, float pointEnergy) {
}

[[maybe_unused]] static void atomicAdd(double* pulseSum, double pointEnergy) {
}


[[maybe_unused]]static std::complex<double> operator+(
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
	return std::complex<double>(x.real() - f, x.imag()); 
}

[[maybe_unused]] static std::complex<double> operator*(
	const float f, 
	const std::complex<double> x) { 
	return std::complex<double>(x.real() * f, x.imag() * f); 
}

[[maybe_unused]] static std::complex<double> operator*(
	const std::complex<double> x, 
	const float f) { 
	return std::complex<double>(x.real() * f, x.imag() * f); 
}

[[maybe_unused]] static std::complex<double> operator/(
	const std::complex<double> x, 
	const float f) { 
	return std::complex<double>(x.real() / f, x.imag() / f); 
}


[[maybe_unused]] static float j0Device(float x) {
	return x;
}

[[maybe_unused]] static double j0Device(double x) {
	return x;
}

template <typename deviceFP, typename deviceComplex>
class counterDevice : public LWEDevice {
private:
	deviceParameterSet<deviceFP, deviceComplex> dParamslocal;
public:
	deviceParameterSet<deviceFP, deviceComplex> deviceStruct;
	deviceParameterSet<deviceFP, deviceComplex>* s;
	deviceParameterSet<deviceFP, deviceComplex>* dParamsDevice;
	counterDevice(simulationParameterSet* sCPU) {
		s = &deviceStruct;
		memoryStatus = 0;
		configuredFFT = true;
		cParams = sCPU;
		dParamsDevice = &dParamslocal;
		sCPU->initializeDeviceParameters(s);
	}
	template<typename T>
	void deviceLaunch(const unsigned int Nblock, const unsigned int Nthread, T functor) {}

	int deviceCalloc(void** ptr, const size_t N, const size_t elementSize) override {return 0;}
	void deviceMemset(void* ptr, int value, size_t count) override {}
	void deviceMemcpyImplementation(void* dst, const void* src, size_t count, copyType kind) override {}
	void deviceFree(void* block) override {}
	void fft(const void* input, void* output, deviceFFT type) override {}
	bool isTheCanaryPixelNaN(const deviceFP* canaryPointer) {
		return false;
	}

	void reset(simulationParameterSet* sCPU) override {
		sCPU->initializeDeviceParameters(s);
	}
};