#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "framework.h"
#include <stdio.h>
#include <complex>
#include <cuComplex.h>
#include "cufft.h"



struct cudaLoop {
	cuDoubleComplex* gridETime;
	cuDoubleComplex* gridETime2;
	cuDoubleComplex* gridETemp;
	cuDoubleComplex* gridETemp2;
	cuDoubleComplex* gridEFrequency;
	cuDoubleComplex* gridEFrequency2;
	cuDoubleComplex* gridPropagationFactor;
	cuDoubleComplex* gridPropagationFactorRho1;
	cuDoubleComplex* gridPropagationFactorRho2;
	cuDoubleComplex* gridPolarizationFactor;
	cuDoubleComplex* gridPolarizationFrequency;
	cuDoubleComplex* gridPropagationFactor2;
	cuDoubleComplex* gridPolarizationFactor2;
	cuDoubleComplex* gridPolarizationFrequency2;
	cuDoubleComplex* gridEFrequencyNext1;
	cuDoubleComplex* gridEFrequencyNext2;
	cuDoubleComplex* gridRadialLaplacian1;
	cuDoubleComplex* gridRadialLaplacian2;
	cuDoubleComplex* k1;
	cuDoubleComplex* k2;

	cuDoubleComplex* ne;
	cuDoubleComplex* no;
	bool isCylindric;
	double* gridPolarizationTime; //future optimization: this could be double rather than complex, if I can figure out the data layout of D2Z cufft
	double* gridPolarizationTime2;
	double* firstDerivativeOperation;
	double* plasmaParameters; //[dt^2 * e^2/m * nonlinearAbsorptionStrength, gamma] 
	double* chi2Tensor;
	double* chi3Tensor;
	double* plasmaCurrent1;
	double* plasmaCurrent2;
	double* absorptionParameters;
	int* nonlinearSwitches;
	long long* propagationInts;
	cufftHandle fftPlan;
	cufftHandle polfftPlan;
	BOOL isNonLinear;
	long long Ntime;
	long long Nspace;
	long long Ngrid;
	int axesNumber;
	int sellmeierType;
	double f0;
	double fStep;
	double dt;
	double dx;
	double h;
	long long Nsteps;
	int Nthread;
	int Nblock;
};

DWORD WINAPI	solveNonlinearWaveEquation(LPVOID lpParam);
int				runRK4Step(struct cudaLoop s, int stepNumber);
int				prepareElectricFieldArrays(struct propthread* s, struct cudaLoop* sc);
int				calcEffectiveChi2Tensor(double* defftensor, double* dtensor, double theta, double phi);
int				fftshiftZ(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int				fftshiftAndFilp(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
std::complex<double> sellmeier(std::complex<double>* ne, std::complex<double>* no, double* a, double f, double theta, double phi, int type, int eqn);
int				preparePropagation2DCartesian(struct propthread* s, struct cudaLoop sc);
int				preparePropagation3DCylindric(struct propthread* s, struct cudaLoop sc);
double			findWalkoffAngles(struct propthread* s, double dk, double f, double tol);
int				loadFrogSpeck(char* frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel, int fieldIndex);
int				plotDataXY(double* X, double* Y, double minX, double maxX, double minY, double maxY, int N, int plotSizeX, int plotSizeY, double lineWidth, double markerWidth, double* plotGrid, double* xTicks, int NxTicks, double* yTicks, int NyTicks);