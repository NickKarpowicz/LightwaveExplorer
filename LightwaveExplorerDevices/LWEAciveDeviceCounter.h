#pragma once

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

void atomicAddCounter(deviceFP* pulseSum, deviceFP pointEnergy) {
}

double j0Counter(deviceFP x) {
	return x;
}

namespace {

	void initializeDeviceParameters(simulationParameterSet* sCPU, deviceParameterSet* s) {
		(*s).Ntime = (*sCPU).Ntime;
		(*s).Nspace = (*sCPU).Nspace;
		(*s).Nspace2 = (*sCPU).Nspace2;
		(*s).is3D = (*sCPU).is3D;
		(*s).Nfreq = ((*s).Ntime / 2 + 1);
		(*s).Ngrid = (*s).Ntime * (*s).Nspace * (*s).Nspace2;
		(*s).NgridC = (*s).Nfreq * (*s).Nspace * (*s).Nspace2; //size of the positive frequency side of the grid
		(*s).fftNorm = 1.0 / (*s).Ngrid;
		(*s).dt = (*sCPU).tStep;
		(*s).dx = (*sCPU).rStep;
		(*s).dk1 = TWOPI / ((*sCPU).Nspace * (*sCPU).rStep);
		(*s).dk2 = TWOPI / ((*sCPU).Nspace2 * (*sCPU).rStep);
		(*s).fStep = (*sCPU).fStep;
		(*s).Nsteps = (size_t)round((*sCPU).crystalThickness / (*sCPU).propagationStep);
		(*s).h = (*sCPU).crystalThickness / ((*s).Nsteps); //adjust step size so that thickness can be varied continuously by fitting
	}
}
class deviceCounter {
private:
	bool configuredFFT = FALSE;
	bool isCylindric = FALSE;
	deviceParameterSet dParamslocal;
	void fftDestroy() {
	}

public:
	int stream;
	int memoryStatus;
	bool hasPlasma;
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

	void fftInitialize(deviceParameterSet* s) {
		hasPlasma = (*s).hasPlasma;
		configuredFFT = 1;
	}
	void deallocateSet(deviceParameterSet* s) {
	}
	void reset(simulationParameterSet* sCPU, deviceParameterSet* s) {
		initializeDeviceParameters(sCPU, s);
	}
	int allocateSet(simulationParameterSet* sCPU, deviceParameterSet* s) {
		return 0;
	}
};