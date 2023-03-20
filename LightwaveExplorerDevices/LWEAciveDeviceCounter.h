#pragma once
#include "LightwaveExplorerUtilities.h"
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


#if LWEFLOATINGPOINT==64

[[maybe_unused]]static std::complex<double> operator+(const float f, const std::complex<double> x) { return std::complex<double>(x.real() + f, x.imag()); }

[[maybe_unused]] static std::complex<double> operator+(const std::complex<double> x, const float f) { return std::complex<double>(x.real() + f, x.imag()); }

[[maybe_unused]] static std::complex<double> operator-(const std::complex<double> x, const float f) { return std::complex<double>(x.real() - f, x.imag()); }

[[maybe_unused]] static std::complex<double> operator*(const float f, const std::complex<double> x) { return std::complex<double>(x.real() * f, x.imag() * f); }

[[maybe_unused]] static std::complex<double> operator*(const std::complex<double> x, const float f) { return std::complex<double>(x.real() * f, x.imag() * f); }

[[maybe_unused]] static std::complex<double> operator/(const std::complex<double> x, const float f) { return std::complex<double>(x.real() / f, x.imag() / f); }

#endif
[[maybe_unused]] static float j0Device(float x) {
	return x;
}

[[maybe_unused]] static double j0Device(double x) {
	return x;
}

template <typename deviceFP, typename deviceComplex>
class counterDevice {
private:
#include "LWEActiveDeviceCommon.cpp"
	bool configuredFFT = 0;
	bool isCylindric = 0;
	deviceParameterSet<deviceFP, deviceComplex> dParamslocal;
	void fftDestroy() {
	}

public:
	int stream;
	int memoryStatus;
	bool hasPlasma;
	deviceParameterSet<deviceFP, deviceComplex> deviceStruct;
	deviceParameterSet<deviceFP, deviceComplex>* s;
	simulationParameterSet* cParams;
	deviceParameterSet<deviceFP, deviceComplex>* dParamsDevice;
	counterDevice(simulationParameterSet* sCPU) {
		s = &deviceStruct;
		memoryStatus = 0;
		stream = 0;
		configuredFFT = 0;
		isCylindric = 0;
		cParams = sCPU;
		dParamsDevice = &dParamslocal;
		initializeDeviceParameters(sCPU);
	}

	~counterDevice() {
	}
	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) {
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) {
		return 0;
	}

	void deviceMemset(void* ptr, int value, size_t count) {
	}

	void deviceMemcpy(void* dst, void* src, size_t count, copyType kind) {
	}

	void deviceFree(void* block) {
	}

	bool isTheCanaryPixelNaN(deviceFP* canaryPointer) {
		return false;
	}

	void fft(void* input, void* output, deviceFFT type) {
	}

	void fftInitialize() {
		hasPlasma = (*s).hasPlasma;
		configuredFFT = 1;
	}
	void deallocateSet() {
	}
	void reset(simulationParameterSet* sCPU) {
		initializeDeviceParameters(sCPU);
	}
	int allocateSet(simulationParameterSet* sCPU) {
		return 0;
	}
};