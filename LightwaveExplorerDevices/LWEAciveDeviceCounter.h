#pragma once
#include "LightwaveExplorerUtilities.h"
#define DeviceToHost 2
#define HostToDevice 1
#define DeviceToDevice 3
#define cudaMemcpyKind int
#define deviceLib std
#define deviceFPLib std
#define complexLib std
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
class activeDevice {
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
	activeDevice(simulationParameterSet* sCPU) {
		s = &deviceStruct;
		memoryStatus = 0;
		stream = 0;
		configuredFFT = 0;
		isCylindric = 0;
		cParams = sCPU;
		dParamsDevice = &dParamslocal;
		(*s).Ntime = (*sCPU).Ntime;
		(*s).Nspace = (*sCPU).Nspace;
		(*s).Nspace2 = (*sCPU).Nspace2;
		(*s).is3D = (*sCPU).is3D;
		(*s).Nfreq = ((*s).Ntime / 2 + 1);
		(*s).Ngrid = (*s).Ntime * (*s).Nspace * (*s).Nspace2;
		(*s).NgridC = (*s).Nfreq * (*s).Nspace * (*s).Nspace2; //size of the positive frequency side of the grid
		(*s).fftNorm = (deviceFP)1.0 / (*s).Ngrid;
		(*s).dt = (deviceFP)(*sCPU).tStep;
		(*s).dx = (deviceFP)(*sCPU).rStep;
		(*s).dk1 = (deviceFP)(twoPi<double>() / ((*sCPU).Nspace * (*sCPU).rStep));
		(*s).dk2 = (deviceFP)(twoPi<double>() / ((*sCPU).Nspace2 * (*sCPU).rStep));
		(*s).fStep = (deviceFP)(*sCPU).fStep;
		(*s).Nsteps = (size_t)round((*sCPU).crystalThickness / (*sCPU).propagationStep);
		(*s).h = (deviceFP)(*sCPU).crystalThickness / ((*s).Nsteps); //adjust step size so that thickness can be varied continuously by fitting
		(*s).axesNumber = (*sCPU).axesNumber;
		(*s).sellmeierType = (*sCPU).sellmeierType;
		(*s).crystalPhi = (deviceFP)(*sCPU).crystalPhi;
		(*s).crystalTheta = (deviceFP)(*sCPU).crystalTheta;
		(*s).f0 = (deviceFP)(*sCPU).pulse1.frequency;
		(*s).Nthread = THREADS_PER_BLOCK;
		(*s).Nblock = (int)((*s).Ngrid / THREADS_PER_BLOCK);
		(*s).NblockC = (int)((*s).NgridC / THREADS_PER_BLOCK);
		(*s).isCylindric = (*sCPU).isCylindric;
		(*s).forceLinear = (*sCPU).forceLinear;
		(*s).isNonLinear = ((*sCPU).nonlinearSwitches[0] + (*sCPU).nonlinearSwitches[1]) > 0;
		(*s).isUsingMillersRule = ((*sCPU).crystalDatabase[(*sCPU).materialIndex].nonlinearReferenceFrequencies[0]) != 0.0;

		if ((*sCPU).nonlinearAbsorptionStrength > 0.) {
			(*s).hasPlasma = TRUE;
			(*s).isNonLinear = TRUE;
		}
		else {
			(*s).hasPlasma = FALSE;
		}

		if ((*s).forceLinear) {
			(*s).hasPlasma = FALSE;
			(*s).isNonLinear = FALSE;
		}
		hasPlasma = (*s).hasPlasma;
	}

	~activeDevice() {
	}
	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) {
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) {
		return 0;
	}

	void deviceMemset(void* ptr, int value, size_t count) {
	}

	void deviceMemcpy(void* dst, void* src, size_t count, cudaMemcpyKind kind) {
	}

	void deviceFree(void* block) {
	}

	bool isTheCanaryPixelNaN(deviceFP* canaryPointer) {
		return FALSE;
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