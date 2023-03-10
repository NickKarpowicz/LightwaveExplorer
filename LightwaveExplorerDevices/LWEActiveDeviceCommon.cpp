#include "LightwaveExplorerUtilities.h"

namespace {
	int fillRotationMatricies(simulationParameterSet* sCPU, deviceParameterSet* s) {
		double cosT = cos((*sCPU).crystalTheta);
		double sinT = sin((*sCPU).crystalTheta);
		double cosP = cos((*sCPU).crystalPhi);
		double sinP = sin((*sCPU).crystalPhi);
		double forward[9] =
		{ cosT * cosP, sinP, -sinT * cosP, -sinP * cosT, cosP, sinP * sinT, sinT, 0, cosT };

		//reverse direction (different order of operations)
		double backward[9] =
		{ cosT * cosP, -sinP * cosT, sinT, sinP, cosP, 0, -sinT * cosP, sinP * sinT, cosT };

		for (size_t i = 0; i < 9; i++) {
			(*s).rotationForward[i] = forward[i];
			(*s).rotationBackward[i] = backward[i];
		}

		return 0;
	}

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
		(*s).axesNumber = (*sCPU).axesNumber;
		(*s).sellmeierType = (*sCPU).sellmeierType;
		(*s).crystalPhi = (*sCPU).crystalPhi;
		(*s).crystalTheta = (*sCPU).crystalTheta;
		(*s).f0 = (*sCPU).pulse1.frequency;
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

	}

	void finishConfiguration(simulationParameterSet* sCPU, deviceParameterSet* s) {
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
		for (size_t i = 0; i < 6; ++i) {
			firstDerivativeOperation[i] *= (-2.0 / ((*s).dx));
		}

		//set nonlinearSwitches[3] to the number of photons needed to overcome bandgap
		(*sCPU).nonlinearSwitches[3] = (int)ceil((*sCPU).bandGapElectronVolts * 241.79893e12 / (*sCPU).pulse1.frequency) - 2;
		double plasmaParametersCPU[6] = { 0 };


		plasmaParametersCPU[0] = (*sCPU).nonlinearAbsorptionStrength; //nonlinear absorption strength parameter
		plasmaParametersCPU[1] = (*sCPU).drudeGamma; //gamma
		if ((*sCPU).nonlinearAbsorptionStrength > 0.) {
			plasmaParametersCPU[2] = (*sCPU).tStep * (*sCPU).tStep
				* 2.817832e-08 / (1.6022e-19 * (*sCPU).bandGapElectronVolts * (*sCPU).effectiveMass); // (dt^2)*e* e / (m * band gap));
		}
		else {
			plasmaParametersCPU[2] = 0;
		}

		for (int j = 0; j < 18; ++j) {
			(*s).chi2Tensor[j] = 2e-12 * (*sCPU).chi2Tensor[j]; //go from d in pm/V to chi2 in m/V
			if (j > 8) (*s).chi2Tensor[j] *= 2.0; //multiply cross-terms by 2 for consistency with convention
		}
		memcpy((*s).nonlinearSwitches, (*sCPU).nonlinearSwitches, 4 * sizeof(int));

		for (size_t i = 0; i < 81; i++) {
			(*s).chi3Tensor[i] = (*sCPU).chi3Tensor[i];
		}

		for (size_t i = 0; i < 6; i++) {
			(*s).absorptionParameters[i] = (*sCPU).absorptionParameters[i];
		}

		for (size_t i = 0; i < 6; i++) {
			(*s).plasmaParameters[i] = plasmaParametersCPU[i];
		}

		for (size_t i = 0; i < 6; i++) {
			(*s).firstDerivativeOperation[i] = firstDerivativeOperation[i];
		}

	}
}