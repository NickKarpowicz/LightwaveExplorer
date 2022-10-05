#include <memory>
#include "LightwaveExplorerUtilities.h"

// the anonymous namespace contains a few helper functions that are identical
// for all classes, taken out here to have less duplicated code
namespace {
	int fillRotationMatricies(simulationParameterSet* sCPU, cudaParameterSet* s) {
		double cosT = cos((*sCPU).crystalTheta);
		double sinT = sin((*sCPU).crystalTheta);
		double cosP = cos((*sCPU).crystalPhi);
		double sinP = sin((*sCPU).crystalPhi);
		double forward[9] =
		{ cosP, sinP, 0, -cosT * sinP, cosT * cosP, sinT, sinT * sinP, -sinT * cosP, cosT };

		//reverse direction (same array contents)
		sinT *= -1;
		sinP *= -1;
		double backward[9] =
		{ cosP, sinP, 0, -cosT * sinP, cosT * cosP, sinT, sinT * sinP, -sinT * cosP, cosT };

		memcpy((*s).rotationForward, forward, 9 * sizeof(double));
		memcpy((*s).rotationBackward, backward, 9 * sizeof(double));
		return 0;
	}

	void initializeDeviceParameters(simulationParameterSet* sCPU, cudaParameterSet* s) {
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
		(*s).axesNumber = (*sCPU).axesNumber;
		(*s).sellmeierType = (*sCPU).sellmeierType;
		(*s).f0 = (*sCPU).frequency1;
		(*s).Nthread = THREADS_PER_BLOCK;
		(*s).Nblock = (int)((*s).Ngrid / THREADS_PER_BLOCK);
		(*s).NblockC = (int)((*s).NgridC / THREADS_PER_BLOCK);
		(*s).isCylindric = (*sCPU).isCylindric;
		(*s).forceLinear = (*sCPU).forceLinear;
		(*s).isNonLinear = ((*sCPU).nonlinearSwitches[0] + (*sCPU).nonlinearSwitches[1]) > 0;
		(*s).isUsingMillersRule = ((*sCPU).crystalDatabase[(*sCPU).materialIndex].nonlinearReferenceFrequencies[0]) != 0;
	}

	void finishConfiguration(simulationParameterSet* sCPU, cudaParameterSet* s) {
		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
		}
		//second polarization grids are to pointers within the first polarization
		//to have contiguous memory
		(*s).gridETime2 = (*s).gridETime1 + (*s).Ngrid;
		(*s).workspace2 = (*s).workspace1 + (*s).NgridC;
		(*s).gridPolarizationTime2 = (*s).gridPolarizationTime1 + (*s).Ngrid;
		(*s).workspace2P = (*s).workspace1 + beamExpansionFactor * (*s).NgridC;
		(*s).k2 = (*s).k1 + (*s).NgridC;
		(*s).chiLinear2 = (*s).chiLinear1 + (*s).Nfreq;
		(*s).fieldFactor2 = (*s).fieldFactor1 + (*s).Nfreq;
		(*s).inverseChiLinear2 = (*s).inverseChiLinear1 + (*s).Nfreq;
		(*s).gridRadialLaplacian2 = (*s).gridRadialLaplacian1 + (*s).Ngrid;
		(*s).gridPropagationFactor1Rho2 = (*s).gridPropagationFactor1Rho1 + (*s).NgridC;
		(*s).gridPolarizationFactor2 = (*s).gridPolarizationFactor1 + (*s).NgridC;
		(*s).gridEFrequency1Next2 = (*s).gridEFrequency1Next1 + (*s).NgridC;
		(*s).gridPropagationFactor2 = (*s).gridPropagationFactor1 + (*s).NgridC;
		(*s).gridEFrequency2 = (*s).gridEFrequency1 + (*s).NgridC;

		double firstDerivativeOperation[6] = { -1. / 60.,  3. / 20., -3. / 4.,  3. / 4.,  -3. / 20., 1. / 60. };
		for (size_t i = 0; i < 6; i++) {
			firstDerivativeOperation[i] *= (-2.0 / ((*s).Ngrid * (*s).dx));
		}

		//set nonlinearSwitches[3] to the number of photons needed to overcome bandgap
		(*sCPU).nonlinearSwitches[3] = (int)ceil((*sCPU).bandGapElectronVolts * 241.79893e12 / (*sCPU).frequency1) - 2;
		double plasmaParametersCPU[6] = { 0 };

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
		plasmaParametersCPU[0] = (*sCPU).nonlinearAbsorptionStrength; //nonlinear absorption strength parameter
		plasmaParametersCPU[1] = (*sCPU).drudeGamma; //gamma
		if ((*sCPU).nonlinearAbsorptionStrength > 0.) {
			plasmaParametersCPU[2] = (*sCPU).tStep * (*sCPU).tStep
				* 2.817832e-08 / (1.6022e-19 * (*sCPU).bandGapElectronVolts * (*sCPU).effectiveMass); // (dt^2)*e* e / (m * band gap));
		}
		else {
			plasmaParametersCPU[2] = 0;
		}

		memcpy((*s).chi2Tensor, (*sCPU).chi2Tensor, 18 * sizeof(double));
		for (int j = 0; j < 18; j++) {
			(*s).chi2Tensor[j] *= 2e-12; //go from d in pm/V to chi2 in m/V
			if (j > 8) (*s).chi2Tensor[j] *= 2.0; //multiply cross-terms by 2 for consistency with convention
		}
		memcpy((*s).nonlinearSwitches, (*sCPU).nonlinearSwitches, 4 * sizeof(int));
		memcpy((*s).chi3Tensor, (*sCPU).chi3Tensor, 81 * sizeof(double));
		memcpy((*s).absorptionParameters, (*sCPU).absorptionParameters, 6 * sizeof(double));
		memcpy((*s).plasmaParameters, plasmaParametersCPU, 6 * sizeof(double));
		memcpy((*s).firstDerivativeOperation, firstDerivativeOperation, 6 * sizeof(double));
	}
}

// define a language specific class and set of macros
#ifdef __CUDACC__
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>
#include <nvml.h>
#define trilingual __global__ void
#define withID 
#define asKernel
#define deviceFunction __device__
#define RUNTYPE 0
#define localIndex threadIdx.x + blockIdx.x * blockDim.x
#define withStream 
#define activeDevice deviceCUDA
#define DeviceToHost cudaMemcpyDeviceToHost
#define HostToDevice cudaMemcpyHostToDevice
#define DeviceToDevice cudaMemcpyDeviceToDevice

class deviceCUDA {
private:
	bool configuredFFT;
	bool isCylindric;
	cufftHandle fftPlanD2Z;
	cufftHandle fftPlanZ2D;
	cufftHandle fftPlan1DD2Z;
	cufftHandle fftPlan1DZ2D;
	cufftHandle doublePolfftPlan;

	int checkDeviceMemory(simulationParameterSet *sCPU, cudaParameterSet *s) {
		nvmlDevice_t nvmlDevice = 0;
		nvmlMemory_t nvmlMemoryInfo;
		nvmlInit_v2();
		nvmlDeviceGetHandleByIndex_v2((*sCPU).assignedGPU, &nvmlDevice);
		nvmlDeviceGetMemoryInfo(nvmlDevice, &nvmlMemoryInfo);
		size_t memoryEstimate = sizeof(double) * ((*s).Ngrid * 2 * 2 + 2 * (*s).NgridC * 6 * 2 + 2 * (*s).isCylindric * 5 * 2 + 2 * (*s).Ntime + 2 * (*s).Nfreq + 81 + 65536);
		nvmlShutdown();
		if (nvmlMemoryInfo.free < memoryEstimate) {
			(*sCPU).memoryError = -1;
			return 1;
		}
		return 0;
	}

	void fftDestroy() {
		cufftDestroy(fftPlanD2Z);
		cufftDestroy(fftPlanZ2D);
		cufftDestroy(fftPlan1DD2Z);
		cufftDestroy(fftPlan1DZ2D);
		if (isCylindric) {
			cufftDestroy(doublePolfftPlan);
		}
		configuredFFT = 0;
	}
public:
	cudaStream_t stream;
	cudaParameterSet* dParams;
	cudaParameterSet* dParamsDevice;
	simulationParameterSet* cParams;
	deviceCUDA() {
		configuredFFT = 0;
		isCylindric = 0;
		cudaStreamCreate(&stream);
		deviceCalloc((void**) &dParamsDevice, 1, sizeof(cudaParameterSet));
	}

	~deviceCUDA() {
		fftDestroy();
		deviceFree(dParamsDevice);
		cudaStreamDestroy(stream);
	}

	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) {
		kernel << <Nblock, Nthread, 0, stream >> > (args...);
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) {
		int err = cudaMalloc(ptr, N * elementSize);
		cudaMemset(*ptr, 0, N * elementSize);
		return err;
	}

	void deviceMemset(void* ptr, int value, size_t count) {
		cudaMemset(ptr, value, count);
	}

	void deviceMemcpy(void* dst, void* src, size_t count, cudaMemcpyKind kind) {
		cudaMemcpy(dst, src, count, kind);
	}

	void deviceFree(void* block) {
		cudaFree(block);
	}

	void fftInitialize(cudaParameterSet *s) {
		if (configuredFFT) {
			fftDestroy();
		}
		isCylindric = 0;
		size_t workSize;
		cufftPlan1d(&fftPlan1DD2Z, (int)(*s).Ntime, CUFFT_D2Z, 2 * (int)((*s).Nspace * (*s).Nspace2));
		cufftPlan1d(&fftPlan1DZ2D, (int)(*s).Ntime, CUFFT_Z2D, 2 * (int)((*s).Nspace * (*s).Nspace2));
		cufftSetStream(fftPlan1DD2Z, stream);
		cufftSetStream(fftPlan1DZ2D, stream);
		if ((*s).is3D) {
			int cufftSizes1[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
			cufftCreate(&fftPlanD2Z);
			cufftGetSizeMany(fftPlanD2Z, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
			cufftMakePlanMany(fftPlanD2Z, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);

			cufftCreate(&fftPlanZ2D);
			cufftGetSizeMany(fftPlanZ2D, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);
			cufftMakePlanMany(fftPlanZ2D, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);
		}
		else {
			int cufftSizes1[] = { (int)(*s).Nspace, (int)(*s).Ntime };

			cufftCreate(&fftPlanD2Z);
			cufftGetSizeMany(fftPlanD2Z, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
			cufftMakePlanMany(fftPlanD2Z, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);

			cufftCreate(&fftPlanZ2D);
			cufftGetSizeMany(fftPlanZ2D, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);
			cufftMakePlanMany(fftPlanZ2D, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);

			if ((*s).isCylindric) {
				isCylindric = 1;
				int cufftSizes2[] = { 2 * (int)(*s).Nspace, (int)(*s).Ntime };
				cufftCreate(&doublePolfftPlan);
				cufftGetSizeMany(doublePolfftPlan, 2, cufftSizes2, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
				cufftMakePlanMany(doublePolfftPlan, 2, cufftSizes2, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
				cufftSetStream(doublePolfftPlan, stream);
			}
		}
		cufftSetStream(fftPlanD2Z, stream);
		cufftSetStream(fftPlanZ2D, stream);
		configuredFFT = 1;
	}

	void fft(void* input, void* output, int type) {
		switch (type) {
		case 0:
			cufftExecD2Z(fftPlanD2Z, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
			break;
		case 1:
			cufftExecZ2D(fftPlanZ2D, (cufftDoubleComplex*)input, (cufftDoubleReal*)output);
			break;
		case 2:
			cufftExecD2Z(fftPlan1DD2Z, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
			break;
		case 3:
			cufftExecZ2D(fftPlan1DZ2D, (cufftDoubleComplex*)input, (cufftDoubleReal*)output);
			break;
		case 4:
			cufftExecD2Z(doublePolfftPlan, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
			break;
		}
	}
	void deallocateSet(cudaParameterSet* s) {
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
	int allocateSet(simulationParameterSet* sCPU, cudaParameterSet* s) {
		dParams = s;
		cParams = sCPU;
		cudaSetDevice((*sCPU).assignedGPU);
		initializeDeviceParameters(sCPU, s);
		fftInitialize(s);
		if (checkDeviceMemory(sCPU, s)) {
			fftDestroy();
			return 1;
		}

		int memErrors = 0;
		double* expGammaTCPU = new double[2 * (*s).Ntime];

		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
		}

		fillRotationMatricies(sCPU, s);

		//GPU allocations
		//
		// currently 8 large grids, meaning memory use is approximately
		// 64 bytes per grid point (8 grids x 2 polarizations x 4ouble precision)
		// plus a little bit for additional constants/workspaces/etc
		memErrors += deviceCalloc((void**)&(*s).gridETime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationTime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).workspace1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1Next1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).k1, 2 * (*s).NgridC, sizeof(std::complex<double>));

		//cylindric sym grids
		if ((*s).isCylindric) {
			memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1Rho1, 4 * (*s).NgridC, sizeof(std::complex<double>));
			memErrors += deviceCalloc((void**)&(*s).gridRadialLaplacian1, 4 * (*s).Ngrid, sizeof(std::complex<double>));
		}

		//smaller helper grids
		memErrors += deviceCalloc((void**)&(*s).expGammaT, 2 * (*s).Ntime, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).fieldFactor1, 2 * (*s).Nfreq, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).inverseChiLinear1, 2 * (*s).Nfreq, sizeof(double));
		for (size_t i = 0; i < (*s).Ntime; i++) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, HostToDevice);
		delete[] expGammaTCPU;
		(*sCPU).memoryError = memErrors;
		if (memErrors > 0) {
			return memErrors;
		}
		finishConfiguration(sCPU, s);
		deviceMemcpy(dParamsDevice, s, sizeof(cudaParameterSet), HostToDevice);
		return 0;
	}

};
#elif defined RUNONSYCL
#include <CL/sycl.hpp>
#include <CL/sycl/atomic.hpp>
#include <oneapi/mkl.hpp>
//#include <dpct/dpct.hpp>
#define trilingual const auto 
#define deviceFunction 
#define RUNTYPE 2
#define DeviceToHost 2
#define HostToDevice 1
#define DeviceToDevice 3
#define localIndex trilingualLaunchID
#define cudaMemcpyKind int
#define asKernel = []
#define withID size_t trilingualLaunchID,
#define withStream , stream
#define activeDevice deviceSYCL
class deviceSYCL {
private:
	bool configuredFFT;
	bool isCylindric;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL> *fftPlanD2Z;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL> *fftPlanZ2D;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL> *fftPlan1DD2Z;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL> *fftPlan1DZ2D;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL> *doublePolfftPlan;

	void fftDestroy() {
		delete fftPlan1DD2Z;
		delete fftPlan1DZ2D;
		delete fftPlanD2Z;
		delete fftPlanZ2D;
		if ((*dParams).isCylindric) delete doublePolfftPlan;
	}
public:
	cl::sycl::queue stream;
	cudaParameterSet* dParams;
	cudaParameterSet* dParamsDevice;
	simulationParameterSet* cParams;

	deviceSYCL() {
		
		configuredFFT = 0;
		isCylindric = 0;
		cl::sycl::default_selector syclDefaultSelector;
		cl::sycl::queue initStream(syclDefaultSelector);
		stream = initStream;
		deviceCalloc((void**)&dParamsDevice, 1, sizeof(cudaParameterSet));
	}

	~deviceSYCL() {
		stream.wait();
		fftDestroy();
		deviceFree(dParamsDevice);
	}

	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) {
		auto event = stream.submit([&](cl::sycl::handler& h) {
			h.parallel_for(Nblock * Nthread, [=](auto i) {kernel(i, args...); });
			});
		event.wait();
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) {
		(*ptr) = cl::sycl::malloc_device(N * elementSize, stream.get_device(), stream.get_context());
		auto event = stream.memset((*ptr), 0, N * elementSize);
		event.wait();
		return 0;
	}

	void deviceMemset(void* ptr, int value, size_t count) {
		auto event = stream.memset(ptr, value, count);
		event.wait();
	}

	void deviceMemcpy(void* dst, void* src, size_t count, cudaMemcpyKind kind) {
		auto event = stream.memcpy(dst, src, count);
		event.wait();
	}

	void deviceFree(void* block) {
		cl::sycl::free(block, stream);
	}

	//to do
	void fftInitialize(cudaParameterSet* s) {
		if (configuredFFT) {
			fftDestroy();
		}
		printf("Starting fft init\n");
		isCylindric = 0;
		size_t workSize;
		int Ntime = (*s).Ntime;
		int Nspace = (*s).Nspace;
		int Nspace2 = (*s).Nspace2;
		int Nfreq = (*s).Nfreq;
		//cufftPlan1d(&fftPlan1DD2Z, (int)(*s).Ntime, CUFFT_D2Z, 2 * (int)((*s).Nspace * (*s).Nspace2));
		fftPlan1DD2Z = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>((int)Ntime);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
		std::int64_t output_stride_ct1[2] = { 0, 1 };
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, output_stride_ct1);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2 * (int)(Nspace * Nspace2));

		fftPlan1DZ2D = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>(Ntime);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
		std::int64_t input_stride_ct2[2] = { 0, 1 };
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, input_stride_ct2);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2 * (Nspace * Nspace2));

		if ((*s).is3D) {
			int cufftSizes1[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
			fftPlanD2Z = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>
				(std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1], cufftSizes1[2]});
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
			
			std::int64_t output_stride_ct11[4] = {0, Nspace * Nfreq, Nfreq, 1 };
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, output_stride_ct11);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace * Nspace2);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace * Nspace2);

			fftPlanZ2D = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>
				(std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1],cufftSizes1[2]});
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
			
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, output_stride_ct11);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace * Nspace2);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace * Nspace2);

		}
		else {
			int cufftSizes1[] = { (int)(*s).Nspace, (int)(*s).Ntime };

			fftPlanD2Z = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>(
					std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1]});
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
			std::int64_t output_stride_ct19[3] = { 0, Nfreq, 1 };
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, output_stride_ct19);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace);

			fftPlanZ2D = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>(
					std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1]});
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, output_stride_ct19);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace);

			if ((*s).isCylindric) {
				isCylindric = 1;
				int cufftSizes2[] = { 2 * (int)(*s).Nspace, (int)(*s).Ntime };
				doublePolfftPlan = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>(
						std::vector<std::int64_t>{cufftSizes2[0], cufftSizes2[1]});
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
				std::int64_t output_stride_ct27[3] = { 0, Nfreq, 1 };
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, output_stride_ct27);
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * 2 * Nspace);
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * 2 * Nspace);
			}
		}
		printf("committing plans\n");
		fftPlan1DD2Z->commit(stream);
		fftPlan1DZ2D->commit(stream);
		fftPlanD2Z->commit(stream);
		fftPlanZ2D->commit(stream);
		if((*s).isCylindric) doublePolfftPlan->commit(stream);
		printf("configured\n");
		configuredFFT = 1;
	}

	//to do
	void fft(void* input, void* output, int type) {
		switch (type) {
		case 0:
			oneapi::mkl::dft::compute_forward(*fftPlanD2Z, (double*)input, (double*)(sycl::double2*)output);
			break;
		case 1:
			oneapi::mkl::dft::compute_backward(*fftPlanZ2D, (double*)(sycl::double2*)input, (double*)output);
			break;
		case 2:
			oneapi::mkl::dft::compute_forward(*fftPlan1DD2Z, (double*)input, (double*)(sycl::double2*)output);
			break;
		case 3:
			oneapi::mkl::dft::compute_backward(*fftPlan1DZ2D, (double*)(sycl::double2*)input, (double*)output);
			break;
		case 4:
			oneapi::mkl::dft::compute_forward(*doublePolfftPlan, (double*)input, (double*)(sycl::double2*)output);
			break;
		}
		stream.wait();
	}
	void deallocateSet(cudaParameterSet* s) {
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
	int allocateSet(simulationParameterSet* sCPU, cudaParameterSet* s) {
		dParams = s;
		cParams = sCPU;
		initializeDeviceParameters(sCPU, s);
		fftInitialize(s);
		int memErrors = 0;
		double* expGammaTCPU = new double[2 * (*s).Ntime];

		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
		}

		fillRotationMatricies(sCPU, s);
		//GPU allocations
		//
		// currently 8 large grids, meaning memory use is approximately
		// 64 bytes per grid point (8 grids x 2 polarizations x 4ouble precision)
		// plus a little bit for additional constants/workspaces/etc
		memErrors += deviceCalloc((void**)&(*s).gridETime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationTime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).workspace1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1Next1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).k1, 2 * (*s).NgridC, sizeof(std::complex<double>));

		//cylindric sym grids
		if ((*s).isCylindric) {
			memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1Rho1, 4 * (*s).NgridC, sizeof(std::complex<double>));
			memErrors += deviceCalloc((void**)&(*s).gridRadialLaplacian1, 4 * (*s).Ngrid, sizeof(std::complex<double>));
		}

		//smaller helper grids
		memErrors += deviceCalloc((void**)&(*s).expGammaT, 2 * (*s).Ntime, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).fieldFactor1, 2 * (*s).Nfreq, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).inverseChiLinear1, 2 * (*s).Nfreq, sizeof(double));
		for (size_t i = 0; i < (*s).Ntime; i++) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, HostToDevice);
		delete[] expGammaTCPU;
		(*sCPU).memoryError = memErrors;
		if (memErrors > 0) {
			return memErrors;
		}
		finishConfiguration(sCPU, s);
		deviceMemcpy(dParamsDevice, s, sizeof(cudaParameterSet), HostToDevice);
		return 0;
	}
};
#else
#define trilingual void 
#define withID size_t trilingualLaunchID, 
#define asKernel
#define deviceFunction 
#define RUNTYPE 1
#define DeviceToHost 2
#define HostToDevice 1
#define DeviceToDevice 3
#define localIndex trilingualLaunchID
#define cudaMemcpyKind int
#define withStream 
#define activeDevice deviceCPU

class deviceCPU {
#include <fftw3_mkl.h>
private:
	bool configuredFFT;
	bool isCylindric;
	cudaParameterSet dParamslocal;
	fftw_plan fftPlanD2Z;
	fftw_plan fftPlanZ2D;
	fftw_plan fftPlan1DD2Z;
	fftw_plan fftPlan1DZ2D;
	fftw_plan doublePolfftPlan;

	void fftDestroy() {
		fftw_destroy_plan(fftPlanD2Z);
		fftw_destroy_plan(fftPlanZ2D);
		fftw_destroy_plan(fftPlan1DD2Z);
		fftw_destroy_plan(fftPlan1DZ2D);
		if (isCylindric)fftw_destroy_plan(doublePolfftPlan);
		fftw_cleanup();
	}

public:
	int stream;
	cudaParameterSet* dParams;
	simulationParameterSet* cParams;
	cudaParameterSet* dParamsDevice;
	deviceCPU() {
		stream = 0;
		dParams = NULL;
		cParams = NULL;
		dParamsDevice = NULL;
		configuredFFT = 0;
		isCylindric = 0;
		fftPlanD2Z = 0;
		fftPlanZ2D = 0;
		fftPlan1DD2Z = 0;
		fftPlan1DZ2D = 0;
		doublePolfftPlan = 0;
		dParamsDevice = &dParamslocal;
	}

	~deviceCPU() {
		fftDestroy();
	}
	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) {
#pragma omp parallel for
		for (int i = 0; i < (int)Nthread; i++) {
			for (unsigned int j = 0u; j < Nblock; j++) {
				kernel(j + Nblock * (unsigned int)i, args...);
			}
		}
	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) {
		(*ptr) = calloc(N, elementSize);
		return (int)((*ptr) == NULL);
	}

	void deviceMemset(void* ptr, int value, size_t count) {
		memset(ptr, value, count);
	}

	void deviceMemcpy(void* dst, void* src, size_t count, cudaMemcpyKind kind) {
		memcpy(dst, src, count);
	}

	void deviceFree(void* block) {
		free(block);
	}

	void fft(void* input, void* output, int type) {
		if (!configuredFFT) return;
		switch (type) {
		case 0:
			fftw_execute_dft_r2c(fftPlanD2Z, (double*)input, (fftw_complex*)output);
			break;
		case 1:
			fftw_execute_dft_c2r(fftPlanZ2D, (fftw_complex*)input, (double*)output);
			break;
		case 2:
			fftw_execute_dft_r2c(fftPlan1DD2Z, (double*)input, (fftw_complex*)output);
			break;
		case 3:
			fftw_execute_dft_c2r(fftPlan1DZ2D, (fftw_complex*)input, (double*)output);
			break;
		case 4:
			fftw_execute_dft_r2c(doublePolfftPlan, (double*)input, (fftw_complex*)output);
			break;
		}
	}

	void fftInitialize(cudaParameterSet *s) {
		if (configuredFFT) fftDestroy();
		isCylindric = (*s).isCylindric;
		double* setupWorkD = new double[(*s).Ngrid * (2 + 2 * (*s).isCylindric)];
		std::complex<double>* setupWorkC = new std::complex<double>[(*s).NgridC * (2 + 2 * (*s).isCylindric)];
		const int fftw1[1] = { (int)(*s).Ntime };
		fftPlan1DD2Z = fftw_plan_many_dft_r2c(1, fftw1, (int)(*s).Nspace * (int)(*s).Nspace2 * 2, setupWorkD, fftw1, 1, (int)(*s).Ntime, (fftw_complex*)setupWorkC, fftw1, 1, (int)(*s).Nfreq, FFTW_ESTIMATE);
		fftPlan1DZ2D = fftw_plan_many_dft_c2r(1, fftw1, (int)(*s).Nspace * (int)(*s).Nspace2 * 2, (fftw_complex*)setupWorkC, fftw1, 1, (int)(*s).Nfreq, setupWorkD, fftw1, 1, (int)(*s).Ntime, FFTW_ESTIMATE);
		if ((*s).is3D) {
			const int fftwSizes[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
			fftPlanD2Z = fftw_plan_many_dft_r2c(3, fftwSizes, 2, setupWorkD, NULL, 1, (int)(*s).Ngrid, (fftw_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, FFTW_MEASURE);
			fftPlanZ2D = fftw_plan_many_dft_c2r(3, fftwSizes, 2, (fftw_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, setupWorkD, NULL, 1, (int)(*s).Ngrid, FFTW_MEASURE);
		}
		else {
			const int fftwSizes[] = { (int)(*s).Nspace, (int)(*s).Ntime };
			fftPlanD2Z = fftw_plan_many_dft_r2c(2, fftwSizes, 2, setupWorkD, NULL, 1, (int)(*s).Ngrid, (fftw_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, FFTW_MEASURE);
			fftPlanZ2D = fftw_plan_many_dft_c2r(2, fftwSizes, 2, (fftw_complex*)setupWorkC, NULL, 1, (int)(*s).NgridC, setupWorkD, NULL, 1, (int)(*s).Ngrid, FFTW_MEASURE);

			if ((*s).isCylindric) {
				const int fftwSizesCyl[] = { (int)(2 * (*s).Nspace), (int)(*s).Ntime };
				doublePolfftPlan = fftw_plan_many_dft_r2c(2, fftwSizesCyl, 2, setupWorkD, NULL, 1, 2 * (int)(*s).Ngrid, (fftw_complex*)setupWorkC, NULL, 1, 2 * (int)(*s).NgridC, FFTW_MEASURE);
			}
		}
		delete[] setupWorkC;
		delete[] setupWorkD;
		configuredFFT = 1;
	}
	void deallocateSet(cudaParameterSet *s) {
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
	int allocateSet(simulationParameterSet *sCPU, cudaParameterSet *s){
		dParams = s;
		cParams = sCPU;
		initializeDeviceParameters(sCPU, s);
		fftInitialize(s);

		int memErrors = 0;
		double* expGammaTCPU = new double[2 * (*s).Ntime];

		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
		}

		fillRotationMatricies(sCPU, s);

		//GPU allocations
		//
		// currently 8 large grids, meaning memory use is approximately
		// 64 bytes per grid point (8 grids x 2 polarizations x 4ouble precision)
		// plus a little bit for additional constants/workspaces/etc
		memErrors += deviceCalloc((void**)&(*s).gridETime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationTime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).workspace1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridPolarizationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).gridEFrequency1Next1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).k1, 2 * (*s).NgridC, sizeof(std::complex<double>));

		//cylindric sym grids
		if ((*s).isCylindric) {
			memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1Rho1, 4 * (*s).NgridC, sizeof(std::complex<double>));
			memErrors += deviceCalloc((void**)&(*s).gridRadialLaplacian1, 4 * (*s).Ngrid, sizeof(std::complex<double>));
		}

		//smaller helper grids
		memErrors += deviceCalloc((void**)&(*s).expGammaT, 2 * (*s).Ntime, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).fieldFactor1, 2 * (*s).Nfreq, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).inverseChiLinear1, 2 * (*s).Nfreq, sizeof(double));
		for (size_t i = 0; i < (*s).Ntime; i++) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		deviceMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, HostToDevice);
		delete[] expGammaTCPU;
		(*sCPU).memoryError = memErrors;
		if (memErrors > 0) {
			return memErrors;
		}
		finishConfiguration(sCPU, s);
		deviceMemcpy(dParamsDevice, s, sizeof(cudaParameterSet), HostToDevice);
		return 0;
	}
};
#endif
