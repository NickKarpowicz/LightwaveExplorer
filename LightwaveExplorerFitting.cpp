#include "LightwaveExplorerCoreCPU.h"
#include "LightwaveExplorerCore.cuh"
#include "LightwaveExplorerUtilities.h"
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

typedef dlib::matrix<double, 0, 1> column_vector;

simulationParameterSet* fittingSet;


#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#define TWOPI 6.2831853071795862
#define PI 3.1415926535897931
#define DEG2RAD 1.7453292519943295e-02
#define LIGHTC 2.99792458e8
#define EPS0 8.8541878128e-12


double getResidual(const column_vector& x) {

	double multipliers[36] = { 0,
1, 1, 1e12, 1e12,
1e12, 1e12, PI, PI,
1e-15, 1e-15, 1e-30, 1e-30,
1e-45, 1e-45, 1e-6, 1e-6,
1e-6, 1e-6,
1e-6, 1e-6, 1e-6, 1e-6,
DEG2RAD, DEG2RAD, DEG2RAD, DEG2RAD,
1, 1, DEG2RAD, DEG2RAD,
1, 1e12, 1, 1e-6,

1e-9 };
	double result = 0.0;
	double* targets[36] = { 0,
	&(*fittingSet).pulseEnergy1, &(*fittingSet).pulseEnergy2, &(*fittingSet).frequency1, &(*fittingSet).frequency2,
	&(*fittingSet).bandwidth1, &(*fittingSet).bandwidth2, &(*fittingSet).cephase1, &(*fittingSet).cephase2,
	&(*fittingSet).delay1, &(*fittingSet).delay2, &(*fittingSet).gdd1, &(*fittingSet).gdd2,
	&(*fittingSet).tod1, &(*fittingSet).tod2, &(*fittingSet).phaseMaterialThickness1, &(*fittingSet).phaseMaterialThickness2,
	&(*fittingSet).beamwaist1, &(*fittingSet).beamwaist2,
	&(*fittingSet).x01, &(*fittingSet).x02, &(*fittingSet).z01, &(*fittingSet).z02,
	&(*fittingSet).propagationAngle1, &(*fittingSet).propagationAngle2, &(*fittingSet).polarizationAngle1, &(*fittingSet).polarizationAngle2,
	&(*fittingSet).circularity1, &(*fittingSet).circularity2, &(*fittingSet).crystalTheta, &(*fittingSet).crystalPhi,
	&(*fittingSet).nonlinearAbsorptionStrength, &(*fittingSet).drudeGamma, &(*fittingSet).effectiveMass, &(*fittingSet).crystalThickness,
	&(*fittingSet).propagationStep };

	for (int i = 0; i < (*fittingSet).Nfitting; i++) {
		*targets[(int)(*fittingSet).fittingArray[3 * i]] = multipliers[(int)(*fittingSet).fittingArray[3 * i]] * x(i);
	}

	if ((*fittingSet).runningOnCPU) {
		if ((*fittingSet).isInSequence) {
			solveNonlinearWaveEquationSequenceCPU(fittingSet);
		}
		else {
			solveNonlinearWaveEquationCPU(fittingSet);
		}
	}

	else {
		if ((*fittingSet).isInSequence) {
			solveNonlinearWaveEquationSequence(fittingSet);
		}
		else {
			solveNonlinearWaveEquation(fittingSet);
		}
	}


	//maximize total spectrum in ROI
	if ((*fittingSet).fittingMode != 3) {
		for (int i = 0; i < (*fittingSet).fittingROIsize; i++) {
			result += (*fittingSet).totalSpectrum[(*fittingSet).fittingMode * (*fittingSet).Nfreq + (*fittingSet).fittingROIstart + i];
		}
		return result;
	}

	//mode 3: match total spectrum to reference given in ascii file
	double a;
	double maxSim = 0;
	double maxRef = 0;
	double sumSim = 0;
	double sumRef = 0;
	double* simSpec = &(*fittingSet).totalSpectrum[2 * (*fittingSet).Nfreq + (*fittingSet).fittingROIstart];
	double* refSpec = &(*fittingSet).fittingReference[(*fittingSet).fittingROIstart];
	for (int i = 0; i < (*fittingSet).fittingROIsize; i++) {
		maxSim = max(maxSim, simSpec[i]);
		maxRef = max(maxRef, refSpec[i]);
		sumSim += simSpec[i];
		sumRef += refSpec[i];
	}

	if (maxSim == 0) {
		maxSim = 1;
	}
	if (maxRef == 0) {
		maxRef = 1;
	}
	result = 0.0;
	for (int i = 0; i < (*fittingSet).fittingROIsize; i++) {
		a = (refSpec[i] / maxRef) - (simSpec[i] / maxSim);
		result += a * a;
	}
	return sqrt(result);
}


unsigned long runDlibFitting(simulationParameterSet* sCPU) {
	printf("Starting ceres\n");
	fittingSet = (simulationParameterSet*)calloc(1, sizeof(simulationParameterSet));
	if (fittingSet == NULL) return 1;
	memcpy(fittingSet, sCPU, sizeof(simulationParameterSet));

	column_vector parameters;
	parameters.set_size((*sCPU).Nfitting);
	column_vector lowerBounds;
	lowerBounds.set_size((*sCPU).Nfitting);
	column_vector upperBounds;
	upperBounds.set_size((*sCPU).Nfitting);
	double* targets[36] = { 0,
	&(*sCPU).pulseEnergy1, &(*sCPU).pulseEnergy2, &(*sCPU).frequency1, &(*sCPU).frequency2,
	&(*sCPU).bandwidth1, &(*sCPU).bandwidth2, &(*sCPU).cephase1, &(*sCPU).cephase2,
	&(*sCPU).delay1, &(*sCPU).delay2, &(*sCPU).gdd1, &(*sCPU).gdd2,
	&(*sCPU).tod1, &(*sCPU).tod2, &(*sCPU).phaseMaterialThickness1, &(*sCPU).phaseMaterialThickness2,
	&(*sCPU).beamwaist1, &(*sCPU).beamwaist2,
	&(*sCPU).x01, &(*sCPU).x02, &(*sCPU).z01, &(*sCPU).z02,
	&(*sCPU).propagationAngle1, &(*sCPU).propagationAngle2, &(*sCPU).polarizationAngle1, &(*sCPU).polarizationAngle2,
	&(*sCPU).circularity1, &(*sCPU).circularity2, &(*sCPU).crystalTheta, &(*sCPU).crystalPhi,
	&(*sCPU).nonlinearAbsorptionStrength, &(*sCPU).drudeGamma, &(*sCPU).effectiveMass, &(*sCPU).crystalThickness,
	&(*sCPU).propagationStep };

	double multipliers[36] = { 0,
	1, 1, 1e12, 1e12,
	1e12, 1e12, PI, PI,
	1e-15, 1e-15, 1e-30, 1e-30,
	1e-45, 1e-45, 1e-6, 1e-6,
	1e-6, 1e-6,
	1e-6, 1e-6, 1e-6, 1e-6,
	DEG2RAD, DEG2RAD, DEG2RAD, DEG2RAD,
	1, 1, DEG2RAD, DEG2RAD,
	1, 1e12, 1, 1e-6,
	1e-9 };

	for (int i = 0; i < (*sCPU).Nfitting; i++) {
		parameters(i) = *targets[(int)(*sCPU).fittingArray[3 * i]];
		lowerBounds(i) = (*sCPU).fittingArray[3 * i + 1];
		upperBounds(i) = (*sCPU).fittingArray[3 * i + 2];
	}

	dlib::function_evaluation result;

	if ((*sCPU).fittingMode != 3) {
		result = dlib::find_max_global(getResidual, lowerBounds, upperBounds, dlib::max_function_calls((*sCPU).fittingMaxIterations));
	}
	else {
		result = dlib::find_min_global(getResidual, lowerBounds, upperBounds, dlib::max_function_calls((*sCPU).fittingMaxIterations));
	}

	for (int i = 0; i < (*sCPU).Nfitting; i++) {
		*targets[(int)round((*sCPU).fittingArray[3 * i])] = multipliers[(int)round((*sCPU).fittingArray[3 * i])] * result.x(i);
		(*sCPU).fittingResult[i] = result.x(i);
	}

	if ((*sCPU).runningOnCPU) {
		if ((*sCPU).isInSequence) {
			solveNonlinearWaveEquationSequenceCPU(sCPU);
		}
		else {
			solveNonlinearWaveEquationCPU(sCPU);
		}
	}

	else {
		if ((*sCPU).isInSequence) {
			solveNonlinearWaveEquationSequence(sCPU);
		}
		else {
			solveNonlinearWaveEquation(sCPU);
		}
	}
	free(fittingSet);

	return 0;
}

unsigned long runDlibFittingCPU(simulationParameterSet* sCPU) {
	printf("Starting ceres\n");
	fittingSet = (simulationParameterSet*)calloc(1, sizeof(simulationParameterSet));
	if (fittingSet == NULL) return 1;
	memcpy(fittingSet, sCPU, sizeof(simulationParameterSet));

	column_vector parameters;
	parameters.set_size((*sCPU).Nfitting);
	column_vector lowerBounds;
	lowerBounds.set_size((*sCPU).Nfitting);
	column_vector upperBounds;
	upperBounds.set_size((*sCPU).Nfitting);
	double* targets[36] = { 0,
	&(*sCPU).pulseEnergy1, &(*sCPU).pulseEnergy2, &(*sCPU).frequency1, &(*sCPU).frequency2,
	&(*sCPU).bandwidth1, &(*sCPU).bandwidth2, &(*sCPU).cephase1, &(*sCPU).cephase2,
	&(*sCPU).delay1, &(*sCPU).delay2, &(*sCPU).gdd1, &(*sCPU).gdd2,
	&(*sCPU).tod1, &(*sCPU).tod2, &(*sCPU).phaseMaterialThickness1, &(*sCPU).phaseMaterialThickness2,
	&(*sCPU).beamwaist1, &(*sCPU).beamwaist2,
	&(*sCPU).x01, &(*sCPU).x02, &(*sCPU).z01, &(*sCPU).z02,
	&(*sCPU).propagationAngle1, &(*sCPU).propagationAngle2, &(*sCPU).polarizationAngle1, &(*sCPU).polarizationAngle2,
	&(*sCPU).circularity1, &(*sCPU).circularity2, &(*sCPU).crystalTheta, &(*sCPU).crystalPhi,
	&(*sCPU).nonlinearAbsorptionStrength, &(*sCPU).drudeGamma, &(*sCPU).effectiveMass, &(*sCPU).crystalThickness,
	&(*sCPU).propagationStep };

	double multipliers[36] = { 0,
	1, 1, 1e12, 1e12,
	1e12, 1e12, PI, PI,
	1e-15, 1e-15, 1e-30, 1e-30,
	1e-45, 1e-45, 1e-6, 1e-6,
	1e-6, 1e-6,
	1e-6, 1e-6, 1e-6, 1e-6,
	DEG2RAD, DEG2RAD, DEG2RAD, DEG2RAD,
	1, 1, DEG2RAD, DEG2RAD,
	1, 1e12, 1, 1e-6,
	1e-9 };

	for (int i = 0; i < (*sCPU).Nfitting; i++) {
		parameters(i) = *targets[(int)(*sCPU).fittingArray[3 * i]];
		lowerBounds(i) = (*sCPU).fittingArray[3 * i + 1];
		upperBounds(i) = (*sCPU).fittingArray[3 * i + 2];
	}

	dlib::function_evaluation result;
	if ((*sCPU).fittingMode != 3) {
		result = dlib::find_max_global(getResidual, lowerBounds, upperBounds, dlib::max_function_calls((*sCPU).fittingMaxIterations));
	}
	else {
		result = dlib::find_min_global(getResidual, lowerBounds, upperBounds, dlib::max_function_calls((*sCPU).fittingMaxIterations));
	}

	for (int i = 0; i < (*sCPU).Nfitting; i++) {
		*targets[(int)round((*sCPU).fittingArray[3 * i])] = multipliers[(int)round((*sCPU).fittingArray[3 * i])] * result.x(i);
		(*sCPU).fittingResult[i] = result.x(i);
	}
	free(fittingSet);

	return 0;
}