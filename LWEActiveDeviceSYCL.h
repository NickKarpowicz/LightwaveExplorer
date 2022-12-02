#include "LWEActiveDeviceCommon.cpp"
#include <CL/sycl.hpp>
#include <CL/sycl/atomic.hpp>
#include <oneapi/mkl/dfti.hpp>
#define DeviceToHost 2
#define HostToDevice 1
#define DeviceToDevice 3
#define cudaMemcpyKind int
void atomicAddSYCL(double* pulseSum, double pointEnergy) {
	cl::sycl::atomic_ref<double, cl::sycl::memory_order::relaxed, cl::sycl::memory_scope::device> a(*pulseSum);
	a.fetch_add(pointEnergy);
}

int hardwareCheckSYCL(int* CUDAdeviceCount) {
	*CUDAdeviceCount = 1;
	return 0;
}

class deviceSYCL {
private:
	bool configuredFFT = FALSE;
	bool isCylindric = FALSE;
	double canaryPixel = 0.0;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>* fftPlanD2Z;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>* fftPlanZ2D;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>* fftPlan1DD2Z;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>* fftPlan1DZ2D;
	oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>* doublePolfftPlan;

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
	int memoryStatus;
	deviceSYCL() {
		memoryStatus = -1;
		configuredFFT = 0;
		isCylindric = 0;
	}

	deviceSYCL(simulationParameterSet* sCPU, cudaParameterSet* s) {
		memoryStatus = -1;
		configuredFFT = 0;
		isCylindric = 0;
		memoryStatus = allocateSet(sCPU, s);
	}

	~deviceSYCL() {
		stream.wait();
		fftDestroy();
		deviceFree(dParamsDevice);
	}

	bool isTheCanaryPixelNaN(double* canaryPointer) {
		deviceMemcpy(&canaryPixel, canaryPointer, sizeof(double), DeviceToHost);
		return(isnan(canaryPixel));
	}

	template<typename Function, typename... Args>
	void deviceLaunch(unsigned int Nblock, unsigned int Nthread, Function kernel, Args... args) {
		stream.submit([&](cl::sycl::handler& h) {
			h.parallel_for(Nblock * Nthread, [=](auto i) {kernel(i, args...); });
			});

	}

	int deviceCalloc(void** ptr, size_t N, size_t elementSize) {
		(*ptr) = cl::sycl::aligned_alloc_device(2 * sizeof(double), N * elementSize, stream.get_device(), stream.get_context());
		stream.memset((*ptr), 0, N * elementSize);
		return 0;
	}

	void deviceMemset(void* ptr, int value, size_t count) {
		stream.wait();
		stream.memset(ptr, value, count);
	}

	void deviceMemcpy(void* dst, void* src, size_t count, cudaMemcpyKind kind) {
		stream.wait();
		stream.memcpy(dst, src, count);
	}

	void deviceFree(void* block) {
		stream.wait();
		cl::sycl::free(block, stream);
	}

	//to do
	void fftInitialize(cudaParameterSet* s) {
		if (configuredFFT) {
			fftDestroy();
		}
		isCylindric = 0;
		size_t workSize;
		int Ntime = (*s).Ntime;
		int Nspace = (*s).Nspace;
		int Nspace2 = (*s).Nspace2;
		int Nfreq = (*s).Nfreq;

		fftPlan1DD2Z = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>((int)Ntime);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
		std::int64_t outputStrides[2] = { 0, 1 };
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, outputStrides);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq);
		fftPlan1DD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2 * (int)(Nspace * Nspace2));

		fftPlan1DZ2D = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>(Ntime);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);

		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, outputStrides);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq);
		fftPlan1DZ2D->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2 * (Nspace * Nspace2));

		if ((*s).is3D) {
			int cufftSizes1[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
			fftPlanD2Z = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>
				(std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1], cufftSizes1[2]});
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);

			std::int64_t outputStride3D[4] = { 0, Nspace * Nfreq, Nfreq, 1 };
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, outputStride3D);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace * Nspace2);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace * Nspace2);

			fftPlanZ2D = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>
				(std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1], cufftSizes1[2]});
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);

			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, outputStride3D);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace * Nspace2);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace * Nspace2);

		}
		else {
			int cufftSizes1[] = { (int)(*s).Nspace, (int)(*s).Ntime };

			fftPlanD2Z = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>(
				std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1]});
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
			std::int64_t outputStride2D[3] = { 0, Nfreq, 1 };
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, outputStride2D);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace);
			fftPlanD2Z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace);

			fftPlanZ2D = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>(
				std::vector<std::int64_t>{cufftSizes1[0], cufftSizes1[1]});
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, outputStride2D);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * Nspace);
			fftPlanZ2D->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * Nspace);

			if ((*s).isCylindric) {
				isCylindric = 1;
				int cufftSizes2[] = { 2 * (int)(*s).Nspace, (int)(*s).Ntime };
				doublePolfftPlan = new oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL>(
					std::vector<std::int64_t>{cufftSizes2[0], cufftSizes2[1]});
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_CONFIG_VALUE::DFTI_NOT_INPLACE);
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 2);
				std::int64_t outputStrideCyl[3] = { 0, Nfreq, 1 };
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, outputStrideCyl);
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, Ntime * 2 * Nspace);
				doublePolfftPlan->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, Nfreq * 2 * Nspace);
			}
		}
		fftPlan1DD2Z->commit(stream);
		fftPlan1DZ2D->commit(stream);
		fftPlanD2Z->commit(stream);
		fftPlanZ2D->commit(stream);
		if ((*s).isCylindric) doublePolfftPlan->commit(stream);
		configuredFFT = 1;
	}

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
	}
	void deallocateSet(cudaParameterSet* s) {
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
	int allocateSet(simulationParameterSet* sCPU, cudaParameterSet* s) {

		cl::sycl::gpu_selector dGPU;
		cl::sycl::cpu_selector dCPU;
		cl::sycl::default_selector d;
		cl::sycl::queue defaultStream(d, { cl::sycl::property::queue::in_order() });
		cl::sycl::queue cpuStream(dCPU, { cl::sycl::property::queue::in_order() });

		if ((*sCPU).assignedGPU == 1) {
			try {
				cl::sycl::queue gpuStream(dGPU, { cl::sycl::property::queue::in_order() });
				stream = gpuStream;
			}
			catch (sycl::exception const& e) {
				stream = cpuStream;
			}
		}
		else if ((*sCPU).runningOnCPU) {
			stream = cpuStream;
		}
		else {
			stream = defaultStream;
		}
		stream.is_in_order();

		deviceCalloc((void**)&dParamsDevice, 1, sizeof(cudaParameterSet));
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
			memErrors += deviceCalloc((void**)&(*s).gridPropagationFactor1Rho1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(std::complex<double>));
			memErrors += deviceCalloc((void**)&(*s).gridRadialLaplacian1, beamExpansionFactor * 2 * (*s).Ngrid, sizeof(std::complex<double>));
		}

		//smaller helper grids
		memErrors += deviceCalloc((void**)&(*s).expGammaT, 2 * (*s).Ntime, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(std::complex<double>));
		memErrors += deviceCalloc((void**)&(*s).fieldFactor1, 2 * (*s).Nfreq, sizeof(double));
		memErrors += deviceCalloc((void**)&(*s).inverseChiLinear1, 2 * (*s).Nfreq, sizeof(double));
		for (size_t i = 0; i < (*s).Ntime; ++i) {
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