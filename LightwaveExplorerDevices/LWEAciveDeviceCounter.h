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
		sCPU->initializeDeviceParameters(s);
		hasPlasma = s->hasPlasma;
	}

	~counterDevice() {
	}
	template<typename T>
	void deviceLaunch(const unsigned int Nblock, const unsigned int Nthread, T functor) {
	}

	int deviceCalloc(void** ptr, const size_t N, const size_t elementSize) {
		return 0;
	}

	void deviceMemset(void* ptr, const int value, const size_t count) {
	}

	void deviceMemcpy(void* dst, const void* src, const size_t count, const copyType kind) {
	}

	void deviceFree(const void* block) {
	}

	bool isTheCanaryPixelNaN(const deviceFP* canaryPointer) {
		return false;
	}

	void fft(const void* input, const void* output, const deviceFFT type) {
	}

	void fftInitialize() {
		hasPlasma = (*s).hasPlasma;
		configuredFFT = 1;
	}
	void deallocateSet() {
	}
	void reset(simulationParameterSet* sCPU) {
		sCPU->initializeDeviceParameters(s);
		hasPlasma = s->hasPlasma;
	}
	int allocateSet(simulationParameterSet* sCPU) {
		return 0;
	}
};