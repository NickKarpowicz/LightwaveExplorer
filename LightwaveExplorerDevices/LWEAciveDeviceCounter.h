#pragma once
#include "LightwaveExplorerUtilities.h"
#define DeviceToHost 2
#define HostToDevice 1
#define DeviceToDevice 3
#define cudaMemcpyKind int
#define deviceLib std
#define deviceFPLib std

int hardwareCheckCounter(int* CUDAdeviceCount) {
	*CUDAdeviceCount = 1;
	return 0;
}

void atomicAddCounter(float* pulseSum, float pointEnergy) {
}

void atomicAddCounter(double* pulseSum, double pointEnergy) {
}



double j0Counter(float x) {
	return x;
}

double j0Counter(double x) {
	return x;
}

class deviceCounter {
private:
#include "LWEActiveDeviceCommon.cpp"
	bool configuredFFT = 0;
	bool isCylindric = 0;
	deviceParameterSet dParamslocal;
	void fftDestroy() {
	}

public:
	int stream;
	int memoryStatus;
	bool hasPlasma;
	deviceParameterSet deviceStruct;
	deviceParameterSet* s;
	simulationParameterSet* cParams;
	deviceParameterSet* dParamsDevice;

	deviceCounter(simulationParameterSet* sCPU) {
		s = &deviceStruct;
		memset(s, 0, sizeof(deviceParameterSet));

		memoryStatus = 0;
		stream = 0;
		cParams = NULL;
		dParamsDevice = NULL;
		configuredFFT = 0;
		isCylindric = 0;
		cParams = sCPU;
		dParamsDevice = &dParamslocal;
		initializeDeviceParameters(sCPU);
	}

	~deviceCounter() {
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

	void fft(void* input, void* output, int type) {
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