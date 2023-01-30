#pragma once
#include "LWEActiveDeviceCommon.cpp"

#define DeviceToHost 2
#define HostToDevice 1
#define DeviceToDevice 3
#define cudaMemcpyKind int

int hardwareCheckCounter(int* CUDAdeviceCount) {
	*CUDAdeviceCount = 1;
	return 0;
}

void atomicAddCounter(double* pulseSum, double pointEnergy) {
}

class deviceCounter {
private:
	bool configuredFFT = FALSE;
	bool isCylindric = FALSE;
	double canaryPixel = 0.0;
	deviceParameterSet dParamslocal;
	void fftDestroy() {
	}

public:
	int stream;
	int memoryStatus;
	deviceParameterSet* dParams;
	simulationParameterSet* cParams;
	deviceParameterSet* dParamsDevice;
	deviceCounter() {
		memoryStatus = -1;
		stream = 0;
		dParams = NULL;
		cParams = NULL;
		dParamsDevice = NULL;
		configuredFFT = 0;
		isCylindric = 0;
		dParamsDevice = &dParamslocal;
	}

	deviceCounter(simulationParameterSet* sCPU, deviceParameterSet* s) {
		memoryStatus = 0;
		stream = 0;
		dParams = NULL;
		cParams = NULL;
		dParamsDevice = NULL;
		configuredFFT = 0;
		isCylindric = 0;
		dParams = s;
		cParams = sCPU;
		dParamsDevice = &dParamslocal;
		initializeDeviceParameters(sCPU, s);

	}

	~deviceCounter() {
	}
	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) {
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) {
		return NULL;
	}

	void deviceMemset(void* ptr, int value, size_t count) {
	}

	void deviceMemcpy(void* dst, void* src, size_t count, cudaMemcpyKind kind) {
	}

	void deviceFree(void* block) {
	}

	bool isTheCanaryPixelNaN(double* canaryPointer) {
		return FALSE;
	}

	void fft(void* input, void* output, int type) {
	}

	void fftInitialize(deviceParameterSet* s) {
		configuredFFT = 1;
	}
	void deallocateSet(deviceParameterSet* s) {
	}
	int allocateSet(simulationParameterSet* sCPU, deviceParameterSet* s) {
		return 0;
	}
};