#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "framework.h"
#include <stdio.h>
#include <complex>
#include <cuComplex.h>
#include "cufft.h"


//    cufftHandle fftPlan, polfftPlan;
// cuDoubleComplex* gridETime, * gridEFrequency, * gridPropagationFactor, * gridPolarizationFactor, * gridPolarizationFrequency, * k1, * k2, * k3, * k4;
// double* gridPolarizationTime;
struct cudaLoop {
	cuDoubleComplex* gridETime;
	cuDoubleComplex* gridETime2;
	cuDoubleComplex* gridETemp;
	cuDoubleComplex* gridETemp2;
	cuDoubleComplex* gridEFrequency;
	cuDoubleComplex* gridEFrequency2;
	cuDoubleComplex* gridPropagationFactor;
	cuDoubleComplex* gridPolarizationFactor;
	cuDoubleComplex* gridPolarizationFrequency;
	cuDoubleComplex* gridPropagationFactor2;
	cuDoubleComplex* gridPolarizationFactor2;
	cuDoubleComplex* gridPolarizationFrequency2;
	cuDoubleComplex* k1;
	cuDoubleComplex* k2;
	cuDoubleComplex* k3;
	cuDoubleComplex* k4;
	cuDoubleComplex* k12;
	cuDoubleComplex* k22;
	cuDoubleComplex* k32;
	cuDoubleComplex* k42;
	cuDoubleComplex* gridPolarizationTime; //future optimization: this could be double rather than complex, if I can figure out the data layout of D2Z cufft
	cuDoubleComplex* gridPolarizationTime2;

	double* chi2Tensor;
	double* chi3Tensor;
	double* absorptionParameters;
	int* nonlinearSwitches;
	long long* propagationInts;
	cufftHandle fftPlan;
	cufftHandle polfftPlan;
	BOOL isNonLinear;
	long long Ntime;
	long long Nspace;
	long long Ngrid;
	double dt;
	double dx;
	double h;
	long long Nsteps;
};

DWORD WINAPI	fftTestCode(LPVOID lpParam);
DWORD WINAPI	propagationLoop(LPVOID lpParam);
int				rkstep(struct cudaLoop s, int stepNumber);
int				pulsegenerator(struct propthread* s, struct cudaLoop* sc);
DWORD WINAPI	propagationLoop(LPVOID lpParam);
int				deff(double* defftensor, double* dtensor, double theta, double phi);
int				fftshiftZ(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int				fftshiftZflip(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
std::complex<double> sellmeier(std::complex<double>* ne, std::complex<double>* no, double* a, double f, double theta, double phi, int type, int eqn);
int				preparepropagation2Dcartesian(struct propthread* s, struct cudaLoop* sc);
double			thetasearch(struct propthread* s, double dk, double f, double tol);