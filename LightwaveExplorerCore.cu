#ifdef __CUDACC__
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <thrust/complex.h>
#include <nvml.h>
#include <cufft.h>
#endif
#include "LightwaveExplorerCore.cuh"
#include "LightwaveExplorerCoreCPU.h"
#include "LightwaveExplorerUtilities.h"
#include <complex>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <mkl.h>
#include <thread>

#define _CRT_SECTURE_NO_WARNINGS

#define THREADS_PER_BLOCK 32
#define MIN_GRIDDIM 8
#define ANGLETOLERANCE 1e-12
#define FALSE 0
#define TRUE 1
#define MAX_LOADSTRING 1024
#define TWOPI 6.2831853071795862
#define PI 3.1415926535897931
#define DEG2RAD 1.7453292519943295e-02
#define LIGHTC 2.99792458e8
#define EPS0 8.8541878128e-12
#define SIXTH 0.1666666666666667
#define THIRD 0.3333333333333333
#define KLORENTZIAN 3182.607353999257 //(e * e / (epsilon_o * m_e)
#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#ifndef __CUDACC__
typedef struct uint3 {
	unsigned int x = 0u;
	unsigned int y = 0u;
	unsigned int z = 0u;
} uint3;
#define cudaMemcpyDeviceToHost 0
#define cudaMemcpyHostToDevice 0
#define cudaMemcpyDeviceToDevice 0
#define cudaMemcpyKind int
#endif


//depending on if the compilation pass is being made by nvcc
//or the c++ compiler, change the prefix to global and
//device kernels. This allows them to work properly on CUDA
//as well as being accessible as usual functions for the
//code to run on the CPU, after applying the bilingualLaunch()
//wrapper defined below.
#ifdef __CUDACC__
#define FGLOBAL __global__
#define FDEVICE __device__
#define GKERN
#define RUNTYPE 0
#else
#define FGLOBAL
#define FDEVICE
#define GKERN uint3 blockIdx, uint3 threadIdx, uint3 blockDim,
#define RUNTYPE 1
#endif

#ifdef __CUDACC__
namespace deviceFunctions {
#else
namespace ordinaryFunctions {
#endif

#ifdef __CUDACC__
	//In tests this mattered, since Thrust does math between complex and double up casting the double to a complex.
	FDEVICE deviceComplex operator/(double a, deviceComplex b) {
		double divByDenominator = a / (b.real() * b.real() + b.imag() * b.imag());
		return deviceComplex(b.real() * divByDenominator, -b.imag() * divByDenominator);
	}
	FDEVICE deviceComplex operator/(deviceComplex a, double b) { return deviceComplex(a.real() / b, a.imag() / b); }

	FDEVICE deviceComplex operator*(double b, deviceComplex a) { return deviceComplex(a.real() * b, a.imag() * b); }
	FDEVICE deviceComplex operator*(deviceComplex a, double b) { return deviceComplex(a.real() * b, a.imag() * b); }

	FDEVICE deviceComplex operator+(double a, deviceComplex b) { return deviceComplex(b.real() + a, b.imag()); }
	FDEVICE deviceComplex operator+(deviceComplex a, double b) { return deviceComplex(a.real() + b, a.imag()); }

	FDEVICE deviceComplex operator-(double a, deviceComplex b) { return deviceComplex(a - b.real(), -b.imag()); }
	FDEVICE deviceComplex operator-(deviceComplex a, double b) { return deviceComplex(a.real() - b, a.imag()); }
#endif
	//Inner function for the Sellmeier equation to provide the refractive indicies
	//current equation form:
	//n^2 = a[0] //background (high freq) contribution
	//      + four resonances, purely real contribution
	//      + parametrized low-frequency correction
	//      + 2 complex-valued Lorenzian contribution
	//inputs:
	//a: 22 component array of the coefficients
	//ls: lamda^2 (microns^2)
	//omega: frequency (rad/s)
	//ii: sqrt(-1)
	//kL: 3183.9 i.e. (e * e / (epsilon_o * m_e)
	FDEVICE deviceComplex sellmeierSubfunctionCuda(
		double* a, double ls, double omega) {
		double realPart = a[0]
			+ (a[1] + a[2] * ls) / (ls + a[3])
			+ (a[4] + a[5] * ls) / (ls + a[6])
			+ (a[7] + a[8] * ls) / (ls + a[9])
			+ (a[10] + a[11] * ls) / (ls + a[12])
			+ a[13] * ls
			+ a[14] * ls * ls
			+ a[15] * ls * ls * ls;

		//traditional sellmeier part is not allowed to give complex values because that almost always
		//means it's out of range and causes instability
		if (realPart < 0) realPart = 1;

		return deviceLib::sqrt(realPart
			+ KLORENTZIAN * a[16] / deviceComplex(a[17] - omega * omega, a[18] * omega)
			+ KLORENTZIAN * a[19] / deviceComplex(a[20] - omega * omega, a[21] * omega));
	}

	FDEVICE deviceComplex lorentzianSubfunctionCuda(
		double* a, double ls, double omega) {
		return deviceLib::sqrt(
			+ KLORENTZIAN * a[0] / deviceComplex(a[1] - omega * omega, a[2] * omega)
			+ KLORENTZIAN * a[3] / deviceComplex(a[4] - omega * omega, a[5] * omega)
			+ KLORENTZIAN * a[6] / deviceComplex(a[7] - omega * omega, a[8] * omega)
			+ KLORENTZIAN * a[9] / deviceComplex(a[10] - omega * omega, a[11] * omega)
			+ KLORENTZIAN * a[12] / deviceComplex(a[13] - omega * omega, a[14] * omega)
			+ KLORENTZIAN * a[15] / deviceComplex(a[16] - omega * omega, a[17] * omega)
			+ KLORENTZIAN * a[18] / deviceComplex(a[19] - omega * omega, a[20] * omega));
	}


	//Sellmeier equation for refractive indicies
	FDEVICE deviceComplex sellmeierCuda(
		deviceComplex* ne, deviceComplex* no, double* a, double f, double theta, double phi, int type, int eqn) {
		if (f == 0) return deviceComplex(1.0, 0.0); //exit immediately for f=0

		//pick which equation to use
		deviceComplex(*sellmeierFunc)(double*, double, double) = &sellmeierSubfunctionCuda;
		if (eqn == 1) sellmeierFunc = &lorentzianSubfunctionCuda;

		double ls = 2.99792458e14 / f; //wavelength in microns
		ls *= ls; //only wavelength^2 is ever used
		double omega = TWOPI * abs(f);

		//option 0: isotropic
		if (type == 0) {
			ne[0] = (*sellmeierFunc)(a, ls, omega);
			no[0] = ne[0];
			return ne[0];
		}
		//option 1: uniaxial
		else if (type == 1) {
			deviceComplex na = (*sellmeierFunc)(a, ls, omega);
			deviceComplex nb = (*sellmeierFunc)(&a[22], ls, omega);
			no[0] = na;
			ne[0] = 1.0 / deviceLib::sqrt(cos(theta) * cos(theta) / (na * na) + sin(theta) * sin(theta) / (nb * nb));
			return ne[0];
		}
		else {
			//type == 2: biaxial
			// X. Yin, S. Zhang and Z. Tian, Optics and Laser Technology 39 (2007) 510 - 513.
			// I am sorry if there is a bug and you're trying to find it, i did my best.
			deviceComplex na = (*sellmeierFunc)(a, ls, omega);
			deviceComplex nb = (*sellmeierFunc)(&a[22], ls, omega);
			deviceComplex nc = (*sellmeierFunc)(&a[44], ls, omega);
			double cosTheta = cos(theta);
			double cosTheta2 = cosTheta * cosTheta;
			double sinTheta = sin(theta);
			double sinTheta2 = sinTheta * sinTheta;
			double sinPhi = sin(phi);
			double sinPhi2 = sinPhi * sinPhi;
			double cosPhi = cos(phi);
			double cosPhi2 = cosPhi * cosPhi;
			double realna2 = na.real() * na.real();
			double realnb2 = nb.real() * nb.real();
			deviceComplex na2 = na * na;
			deviceComplex nb2 = nb * nb;
			deviceComplex nc2 = nc * nc;
			double delta = 0.5 * atan(-((1. / realna2 - 1. / realnb2)
				* sin(2 * phi) * cosTheta) / ((cosPhi2 / realna2 + sinPhi2 / realnb2)
					+ ((sinPhi2 / realna2 + cosPhi2 / realnb2)
						* cosTheta2 + sinTheta2 / (nc.real() * nc.real()))));
			double cosDelta = cos(delta);
			double sinDelta = sin(delta);
			ne[0] = 1.0 / deviceLib::sqrt(cosDelta * cosDelta * (cosTheta2 * (cosPhi2 / na2
				+ sinPhi2 / nb2) + sinTheta2 / nc2)
				+ sinDelta * sinDelta * (sinPhi2 / na2 + cosPhi2 / nb2)
				- 0.5 * sin(2 * phi) * cosTheta * sin(2 * delta) * (1. / na2 - 1. / (nb * nb)));

			no[0] = 1.0 / deviceLib::sqrt(sinDelta * sinDelta * (cosTheta2 * (cosPhi2 / na2
				+ sinPhi2 / nb2) + sinTheta2 / nc2)
				+ cosDelta * cosDelta * (sinPhi2 / na2 + cosPhi2 / nb2)
				+ 0.5 * sin(2 * phi) * cosTheta * sin(2 * delta) * (1. / na2 - 1. / nb2));
			return ne[0];
		}
	}

	FDEVICE double cuCModSquared(deviceComplex a) {
		return a.real() * a.real() + a.imag() * a.imag();
	}

	//provide a list of nearest-3 neighbors for taking spatial derivatives
	// exploiting the fact that the radial grid is offset by 1/4 step from 0
	// this means that midpoints are available on the other side of the origin.
	// returns rho at the given index j
	FDEVICE double resolveNeighborsInOffsetRadialSymmetry(
		long long* neighbors, long long N, int j, double dr, long long Ntime, long long h) {
		if (j < N / 2) {
			neighbors[0] = (N - j - 2) * Ntime + h;
			neighbors[1] = (j + 1) * Ntime + h;
			neighbors[2] = (N - j - 1) * Ntime + h;
			neighbors[3] = (N - j) * Ntime + h;
			neighbors[4] = (j - 1) * Ntime + h;
			neighbors[5] = (N - j + 1) * Ntime + h;
			return -(dr * (j - N / 2) + 0.25 * dr);
		}
		else {
			neighbors[0] = (N - j + 1) * Ntime + h;
			neighbors[1] = (j - 1) * Ntime + h;
			neighbors[2] = (N - j) * Ntime + h;
			neighbors[3] = (N - j - 1) * Ntime + h;
			neighbors[4] = (j + 1) * Ntime + h;
			neighbors[5] = (N - j - 2) * Ntime + h;
			return dr * (j - N / 2) + 0.25 * dr;
		}
	}


	//function to find the effective crystal direction (theta,phi) for a given coordinate in k-space
	// trivial in isotropic media (Snell's law) but in birefringent media,
	// the refracted angle depends on refractive index, but refractive
	// index depends on the refracted angle. 
	//  Using conservation of the transverse momentum, we must solve the system of equations:
	//	kx1 = kx2
	//  ky1 = ky2
	//  where 
	//  kx1 = (w/c)sin(alpha_i)
	//  kx2 = (n(theta+alpha,phi+beta)w/c)*sin(alpha)
	//  ky1 = (w/c)sin(beta_i)
	//  ky2 = (n(theta+alpha,phi+beta)w/c)*sin(beta)
	//
	// The k grid is known, meaning kx1 and ky1 are givens; alpha and beta are unknowns
	// minimize n(alpha,beta)*sin(alpha) - kx1*c/w, n(alpha,beta)*sin(beta) - ky1*c/w
	//
	// starting point n = n(alpha=0, beta=0), which is the actual solution for isotropic medium
	// If isotropic, return
	// If uniaxial, solve 1D problem with n(alpha,0)
	// If biaxial, solve 2D problem
	// Use OGM1; D. Kim, J.A. Fessler, Optimized first-order methods for smooth convex minimization, arXiv:1406.5468
	FDEVICE void findBirefringentCrystalIndex(cudaParameterSet* s, double* sellmeierCoefficients, long long i, deviceComplex* n1, deviceComplex* n2) {
		unsigned long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;

		double f = (*s).fStep * h;
		double kx1 = (LIGHTC / (TWOPI * f)) * (j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace));
		double ky1 = (LIGHTC / (TWOPI * f)) * (k * (*s).dk2 - (k >= ((*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2));
		//alpha is deviation from crystal Theta (x2 polarizations)
		//beta is deviation from crystal Phi
		//
		deviceComplex n[4][2];
		deviceComplex nW;
		sellmeierCuda(&n[0][0], &n[0][1], sellmeierCoefficients, f, sellmeierCoefficients[66], sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
		if ((*s).axesNumber == 0) {
			*n1 = n[0][0];
			*n2 = n[0][1];
			return;
		}

		double gradient[2][2];
		double alpha[2] = { asin(kx1 / n[0][0].real()),asin(kx1 / n[0][1].real()) };
		double beta[2] = { asin(ky1 / n[0][0].real()),asin(ky1 / n[0][1].real()) };

		double gradientStep = 1.0e-7;
		double gradientFactor = 0.5 / gradientStep;
		int it;
		int maxiter = 32;
		//emperical testing: 
		// converges to double precision limit in two iterations for BBO
		// converges in 32 iterations in BiBO

		double errArray[4][2];
		if ((*s).axesNumber == 1) {
			maxiter = 4;
			sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] + gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&n[1][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] - gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
			errArray[0][0] = sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
			errArray[1][0] = sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
			gradient[0][0] = gradientFactor * (errArray[0][0] - errArray[1][0]);

			for (it = 0; it < maxiter; it++) {
				if (abs(gradient[0][0]) > 1e-13) alpha[0] -= 0.5 * (errArray[0][0] + errArray[1][0]) / gradient[0][0];

				sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] + gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&n[1][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] - gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
				errArray[0][0] = sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
				errArray[1][0] = sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
				gradient[0][0] = gradientFactor * (errArray[0][0] - errArray[1][0]);
			}
			sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0], sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&nW, &n[1][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1], sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
			*n1 = n[0][0];
			*n2 = n[1][1];
			return;
		}

		if ((*s).axesNumber == 2) {
			sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] + gradientStep, sellmeierCoefficients[67] + beta[0], (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&n[1][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] - gradientStep, sellmeierCoefficients[67] + beta[0], (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&n[2][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0], sellmeierCoefficients[67] + beta[0] + gradientStep, (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&n[3][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0], sellmeierCoefficients[67] + beta[0] - gradientStep, (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&nW, &n[0][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1] + gradientStep, sellmeierCoefficients[67] + beta[1], (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&nW, &n[1][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1] - gradientStep, sellmeierCoefficients[67] + beta[1], (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&nW, &n[2][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1], sellmeierCoefficients[67] + beta[1] + gradientStep, (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&nW, &n[3][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1], sellmeierCoefficients[67] + beta[1] - gradientStep, (*s).axesNumber, (*s).sellmeierType);
			errArray[0][0] = sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
			errArray[1][0] = sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
			errArray[2][0] = sin(beta[0] + gradientStep) * n[2][0].real() - ky1;
			errArray[3][0] = sin(beta[0] - gradientStep) * n[3][0].real() - ky1;
			errArray[0][1] = sin(alpha[1] + gradientStep) * n[0][1].real() - kx1;
			errArray[1][1] = sin(alpha[1] - gradientStep) * n[1][1].real() - kx1;
			errArray[2][1] = sin(beta[1] + gradientStep) * n[2][1].real() - ky1;
			errArray[3][1] = sin(beta[1] - gradientStep) * n[3][1].real() - ky1;
			gradient[0][0] = gradientFactor * (errArray[0][0] - errArray[1][0]);
			gradient[1][0] = gradientFactor * (errArray[2][0] - errArray[3][0]);
			gradient[0][1] = gradientFactor * (errArray[0][1] - errArray[1][1]);
			gradient[1][1] = gradientFactor * (errArray[2][1] - errArray[3][1]);

			for (it = 0; it < maxiter; it++) {
				if (abs(gradient[0][0]) > 1e-13) alpha[0] -= 0.25 * (errArray[0][0] + errArray[1][0]) / gradient[0][0];
				if (abs(gradient[1][0]) > 1e-13) beta[0] -= 0.25 * (errArray[2][0] + errArray[3][0]) / gradient[1][0];
				if (abs(gradient[0][1]) > 1e-13) alpha[1] -= 0.25 * (errArray[0][1] + errArray[1][1]) / gradient[0][1];
				if (abs(gradient[1][1]) > 1e-13) beta[1] -= 0.25 * (errArray[2][1] + errArray[3][1]) / gradient[1][1];
				sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] + gradientStep, sellmeierCoefficients[67] + beta[0], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&n[1][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] - gradientStep, sellmeierCoefficients[67] + beta[0], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&n[2][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0], sellmeierCoefficients[67] + beta[0] + gradientStep, (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&n[3][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0], sellmeierCoefficients[67] + beta[0] - gradientStep, (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&nW, &n[0][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1] + gradientStep, sellmeierCoefficients[67] + beta[1], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&nW, &n[1][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1] - gradientStep, sellmeierCoefficients[67] + beta[1], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&nW, &n[2][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1], sellmeierCoefficients[67] + beta[1] + gradientStep, (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&nW, &n[3][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1], sellmeierCoefficients[67] + beta[1] - gradientStep, (*s).axesNumber, (*s).sellmeierType);
				errArray[0][0] = sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
				errArray[1][0] = sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
				errArray[2][0] = sin(beta[0] + gradientStep) * n[2][0].real() - ky1;
				errArray[3][0] = sin(beta[0] - gradientStep) * n[3][0].real() - ky1;
				errArray[0][1] = sin(alpha[1] + gradientStep) * n[0][1].real() - kx1;
				errArray[1][1] = sin(alpha[1] - gradientStep) * n[1][1].real() - kx1;
				errArray[2][1] = sin(beta[1] + gradientStep) * n[2][1].real() - ky1;
				errArray[3][1] = sin(beta[1] - gradientStep) * n[3][1].real() - ky1;
				gradient[0][0] = gradientFactor * (errArray[0][0] - errArray[1][0]);
				gradient[1][0] = gradientFactor * (errArray[2][0] - errArray[3][0]);
				gradient[0][1] = gradientFactor * (errArray[0][1] - errArray[1][1]);
				gradient[1][1] = gradientFactor * (errArray[2][1] - errArray[3][1]);
			}
			sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0], sellmeierCoefficients[67] + beta[0], (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&nW, &n[1][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1], sellmeierCoefficients[67] + beta[1], (*s).axesNumber, (*s).sellmeierType);
			*n1 = n[0][0];
			*n2 = n[1][1];
			return;
		}

	}

	FDEVICE void findBirefingentCrystalAngle(double* alphaE, double* alphaO, unsigned long long j, double f, double* sellmeierCoefficients, cudaParameterSet* s) {
		//Find walkoff angle, starting from zero
		// in the case of an extraordinary axis, the angle of propagation is related to the transverse
		// momentum in a complicated way:
		// sin(theta) * n(theta) = delta k * c/omega
		// theta depends on the refractive index, and the refractive index depends on theta
		// so we solve numerically
		double dAlpha = 0.1;
		double nePlus, neMinus;
		double err, errPlus, errMinus;
		deviceComplex ne, no;


		deviceComplex ii = deviceComplex(0, 1);
		double crystalTheta = sellmeierCoefficients[66];
		double crystalPhi = sellmeierCoefficients[67];
		double kStep = sellmeierCoefficients[70];
		double tol = sellmeierCoefficients[72];
		double dk = j * kStep - (j >= ((*s).Nspace / 2)) * (kStep * (*s).Nspace); //frequency grid in transverse direction
		double rhs = LIGHTC * dk / (TWOPI * f);

		//if not biaxial, the o-axis can be solved analytically.
		sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
			crystalTheta, crystalPhi, (*s).axesNumber, (*s).sellmeierType);
		*alphaO = asin(rhs / no.real());
		if ((*s).axesNumber == 2) {
			sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
				crystalTheta + *alphaO, crystalPhi, (*s).axesNumber, (*s).sellmeierType);
			nePlus = no.real();
			err = abs(nePlus * sin(*alphaO) - rhs);

			int iters = 0;
			errPlus = 2;
			errMinus = 2;
			while (err > tol && iters < 2048) {
				iters++;

				sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
					crystalTheta + *alphaO + dAlpha, crystalPhi, (*s).axesNumber, (*s).sellmeierType);
				nePlus = no.real();
				errPlus = abs(nePlus * sin(*alphaO + dAlpha) - rhs);

				sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
					crystalTheta + *alphaO - dAlpha, crystalPhi, (*s).axesNumber, (*s).sellmeierType);
				neMinus = no.real();
				errMinus = abs(neMinus * sin(*alphaO - dAlpha) - rhs);

				//Basic hill climbing algorithm
				//calculate the error at theta +/- dTheta
				// if theta + dTheta has lowest error, theta = theta+dTheta, err = errPlus
				// if theta - dTheta has lowest error, theta = theta-dTheta, err = errMinus
				// if theta has lowest error, step size is too large, dTheta /= 2;
				if (errPlus < err && errPlus < errMinus) {
					*alphaO += dAlpha;
					err = errPlus;
				}
				else if (errMinus < err) {
					*alphaO -= dAlpha;
					err = errMinus;
				}
				else {
					dAlpha *= 0.5;
				}

			}
		}

		//find the extraordinary angle if the crystal isn't isotropic
		*alphaE = *alphaO;
		if ((*s).axesNumber > 0) {
			sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
				crystalTheta + *alphaE, crystalPhi, (*s).axesNumber, (*s).sellmeierType);
			nePlus = ne.real();
			err = abs(nePlus * sin(*alphaE) - rhs);

			int iters = 0;
			errPlus = 2;
			errMinus = 2;
			dAlpha = 0.1;
			while (err > tol && iters < 2048) {
				iters++;

				sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
					crystalTheta + *alphaE + dAlpha, crystalPhi, (*s).axesNumber, (*s).sellmeierType);
				nePlus = ne.real();
				errPlus = abs(nePlus * sin(*alphaE + dAlpha) - rhs);

				sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
					crystalTheta + *alphaE - dAlpha, crystalPhi, (*s).axesNumber, (*s).sellmeierType);
				neMinus = ne.real();
				errMinus = abs(neMinus * sin(*alphaE - dAlpha) - rhs);

				//Basic hill climbing algorithm
				//calculate the error at theta +/- dTheta
				// if theta + dTheta has lowest error, theta = theta+dTheta, err = errPlus
				// if theta - dTheta has lowest error, theta = theta-dTheta, err = errMinus
				// if theta has lowest error, step size is too large, dTheta /= 2;
				if (errPlus < err && errPlus < errMinus) {
					*alphaE += dAlpha;
					err = errPlus;
				}
				else if (errMinus < err) {
					*alphaE -= dAlpha;
					err = errMinus;
				}
				else {
					dAlpha *= 0.5;
				}
			}
		}
	}
}

#ifdef __CUDACC__
using namespace deviceFunctions;
#else
using namespace ordinaryFunctions;
#endif
FGLOBAL void millersRuleNormalizationKernel(GKERN cudaParameterSet* s, double* sellmeierCoefficients, double* referenceFrequencies) {
	if (!(*s).isUsingMillersRule) {
		return;
	}
	size_t i;
	double chi11[7];
	double chi12[7];
	deviceComplex ne, no;
	for (i = 0; i < 7; i++) {
		if (referenceFrequencies[i] == 0) {
			chi11[i] = 100000.0;
			chi12[i] = 100000.0;
		}
		else {
			sellmeierCuda(&ne, &no, sellmeierCoefficients, referenceFrequencies[i], sellmeierCoefficients[66], sellmeierCoefficients[67], (int)sellmeierCoefficients[69], 0);
			chi11[i] =ne.real() *ne.real() - 1;
			chi12[i] =no.real() *no.real() - 1;
		}
	}

	//normalize chi2 tensor values
	(*s).chi2Tensor[0] /= chi11[0] * chi11[1] * chi11[2];
	(*s).chi2Tensor[1] /= chi11[0] * chi11[1] * chi12[2];
	(*s).chi2Tensor[2] /= chi11[0] * chi12[1] * chi11[2];
	(*s).chi2Tensor[3] /= chi11[0] * chi12[1] * chi12[2];
	(*s).chi2Tensor[4] /= chi12[0] * chi12[1] * chi11[2];
	(*s).chi2Tensor[5] /= chi12[0] * chi12[1] * chi12[2];

	//normalize chi3 tensor values
	for (char i = 0; i < 81; i++) {
		(*s).chi3Tensor[i] /= chi11[3] * chi11[4] * chi11[5] * chi11[6];
	}
}

FGLOBAL void totalSpectrumKernel(GKERN deviceComplex* fieldGrid1, deviceComplex* fieldGrid2, double gridStep, size_t Ntime, size_t Nspace, double* spectrum) {
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;
	size_t j;
	double beamCenter1 = 0.;
	double beamCenter2 = 0.;
	double beamTotal1 = 0.;
	double beamTotal2 = 0.;
	double a, x;

	//find beam centers
	for (j = 0; j < Nspace; j++) {
		x = gridStep * j;
		a = cuCModSquared(fieldGrid1[i + j * Ntime]);
		beamTotal1 += a;
		beamCenter1 += x * a;
		a = cuCModSquared(fieldGrid2[i + j * Ntime]);
		beamTotal2 += a;
		beamCenter2 += x * a;
	}
	if (beamTotal1 > 0) {
		beamCenter1 /= beamTotal1;
	}
	if (beamTotal2 > 0) {
		beamCenter2 /= beamTotal2;
	}


	//Integrate total beam power, assuming radially-symmetric beam around
	//the center
	beamTotal1 = 0.;
	beamTotal2 = 0.;
	for (j = 0; j < Nspace; j++) {
		x = gridStep * j;
		beamTotal1 += PI * abs(x - beamCenter1) * cuCModSquared(fieldGrid1[i + j * Ntime]);
		beamTotal2 += PI * abs(x - beamCenter2) * cuCModSquared(fieldGrid2[i + j * Ntime]);
	}
	beamTotal1 *= gridStep / Ntime;
	beamTotal2 *= gridStep / Ntime;

	//put the values into the output spectrum
	spectrum[i] = beamTotal1;
	spectrum[i + Ntime] = beamTotal2;
	spectrum[i + 2 * Ntime] = beamTotal1 + beamTotal2;
}

FGLOBAL void totalSpectrum3DKernel(GKERN deviceComplex* fieldGrid1, deviceComplex* fieldGrid2, double gridStep, size_t Ntime, size_t Nspace, double* spectrum) {
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;
	size_t j;

	double beamTotal1 = 0.;
	double beamTotal2 = 0.;
	//Integrate total beam power
	beamTotal1 = 0.;
	beamTotal2 = 0.;
	for (j = 0; j < Nspace; j++) {
		beamTotal1 += cuCModSquared(fieldGrid1[i + j * Ntime]);
		beamTotal2 += cuCModSquared(fieldGrid2[i + j * Ntime]);
	}
	beamTotal1 *= gridStep * gridStep / Ntime;
	beamTotal2 *= gridStep * gridStep / Ntime;

	//put the values into the output spectrum
	spectrum[i] = beamTotal1;
	spectrum[i + Ntime] = beamTotal2;
	spectrum[i + 2 * Ntime] = beamTotal1 + beamTotal2;
}

//rotate the field around the propagation axis (basis change)
FGLOBAL void rotateFieldKernel(GKERN deviceComplex* Ein1, deviceComplex* Ein2, deviceComplex* Eout1,
	deviceComplex* Eout2, double rotationAngle) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	Eout1[i] = cos(rotationAngle) * Ein1[i] - sin(rotationAngle) * Ein2[i];
	Eout2[i] = sin(rotationAngle) * Ein1[i] + cos(rotationAngle) * Ein2[i];
}



FGLOBAL void radialLaplacianKernel(GKERN cudaParameterSet* s) {
	unsigned long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j = i / (*s).Ntime; //spatial coordinate
	long long h = i % (*s).Ntime; //temporal coordinate
	long long neighbors[6];

	//zero at edges of grid
	if (j<3 || j>((long long)(*s).Nspace - 4)) {
		(*s).gridRadialLaplacian1[i] = 0.;
		(*s).gridRadialLaplacian2[i] = 0.;
	}
	else {
		double rho = resolveNeighborsInOffsetRadialSymmetry(neighbors, (*s).Nspace, (int)j, (*s).dx, (*s).Ntime, h);
		rho = -1.0 / rho;
		(*s).gridRadialLaplacian1[i] = rho * ((*s).firstDerivativeOperation[0] * (*s).gridETime1[neighbors[0]]
			+ (*s).firstDerivativeOperation[1] * (*s).gridETime1[neighbors[1]]
			+ (*s).firstDerivativeOperation[2] * (*s).gridETime1[neighbors[2]]
			+ (*s).firstDerivativeOperation[3] * (*s).gridETime1[neighbors[3]]
			+ (*s).firstDerivativeOperation[4] * (*s).gridETime1[neighbors[4]]
			+ (*s).firstDerivativeOperation[5] * (*s).gridETime1[neighbors[5]]);
		(*s).gridRadialLaplacian2[i] = rho * ((*s).firstDerivativeOperation[0] * (*s).gridETime2[neighbors[0]]
			+ (*s).firstDerivativeOperation[1] * (*s).gridETime2[neighbors[1]]
			+ (*s).firstDerivativeOperation[2] * (*s).gridETime2[neighbors[2]]
			+ (*s).firstDerivativeOperation[3] * (*s).gridETime2[neighbors[3]]
			+ (*s).firstDerivativeOperation[4] * (*s).gridETime2[neighbors[4]]
			+ (*s).firstDerivativeOperation[5] * (*s).gridETime2[neighbors[5]]);
	}

}
//Expand the information contained in the radially-symmetric beam in the offset grid
// representation.
// The grid is offset from the origin; rather than ...-2 -1 0 1 2... etc, which would
// contain redundant information (the symmetry means that -1 and -1 are equivalent)
// the grid is at the points -1.75 -0.75 0.25 1.25 2.25, etc.
// the grid spacing is the same, but now the two sides of the origin contain different
// information. This has effectively doubled the resolution of the nonlinear
// polarization. 
// We make use of this by expanding into the full-resolution beam on the grid
// -2.25 -1.75 -1.25 -0.75 -0.25 0.25 0.75 1.25 1.75 2.25...
// after FFT, we can discard the high frequencies. Thus we have downsampled
// in such a way as to avoid aliasing, which inside the simulation is most
// likely the appear (and cause instability) in the nonlinear terms.
FGLOBAL void expandCylindricalBeam(GKERN cudaParameterSet* s, double* polarization1, double* polarization2) {
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;
	size_t j = i / (*s).Ntime; //spatial coordinate
	size_t k = i % (*s).Ntime; //temporal coordinate

	//positions on the expanded grid corresponding the the current index
	size_t pos1 = 2 * ((*s).Nspace - j - 1) * (*s).Ntime + k;
	size_t pos2 = (2 * j + 1) * (*s).Ntime + k;

	//reuse memory allocated for the radial Laplacian, casting complex double
	//to a 2x larger double real grid
	double* expandedBeam1 = (double*)(*s).gridRadialLaplacian1;
	double* expandedBeam2 = expandedBeam1 + 2 * (*s).Ngrid;

	expandedBeam1[pos1] = polarization1[i];
	expandedBeam1[pos2] = polarization1[i];
	expandedBeam2[pos1] = polarization2[i];
	expandedBeam2[pos2] = polarization2[i];
}



//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
//note that the sellmeier coefficients have extra values appended to the end
//to give info about the current simulation
FGLOBAL void applyFresnelLossKernel(GKERN double* sellmeierCoefficients1, double* sellmeierCoefficients2, cudaParameterSet* s) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	double alpha1, alpha2, alphaO1, alphaO2;
	long long j, k;
	long long Ntime = (*s).Ntime;
	int axesNumber = (*s).axesNumber;
	int sellmeierType = (*s).sellmeierType;
	deviceComplex ne1, no1, ne2, no2, n0;
	deviceComplex cuZero = deviceComplex(0, 0);
	j = i / Ntime; //spatial coordinate
	k = i % Ntime; //temporal coordinate
	deviceComplex ii = deviceComplex(0, 1);
	double crystalTheta = sellmeierCoefficients1[66];
	double crystalPhi = sellmeierCoefficients1[67];
	double fStep = sellmeierCoefficients1[71];

	//frequency being resolved by current thread
	double f = k * fStep;


	findBirefingentCrystalAngle(&alpha1, &alphaO1, j, f, sellmeierCoefficients1, s);
	findBirefingentCrystalAngle(&alpha2, &alphaO2, j, f, sellmeierCoefficients2, s);
	//walkoff angle has been found, generate the rest of the grids


	sellmeierCuda(&ne1, &no1, sellmeierCoefficients1, f,
		crystalTheta + 0*alpha1, crystalPhi, axesNumber, sellmeierType);
	sellmeierCuda(&n0, &no1, sellmeierCoefficients1, f,
		crystalTheta + 0*alphaO1, crystalPhi, axesNumber, sellmeierType);
	if (isnan(ne1.real()) || isnan(no1.real())) {
		ne1 = deviceComplex(1, 0);
		no1 = deviceComplex(1, 0);
	}


	sellmeierCuda(&ne2, &no2, sellmeierCoefficients2, f,
		crystalTheta + alpha2, crystalPhi, axesNumber, sellmeierType);
	sellmeierCuda(&n0, &no2, sellmeierCoefficients2, f,
		crystalTheta + alphaO2, crystalPhi, axesNumber, sellmeierType);
	if (isnan(ne2.real()) || isnan(no2.real())) {
		ne2 = deviceComplex(1, 0);
		no2 = deviceComplex(1, 0);
	}

	deviceComplex ts = 2. * ne1 * cos(alpha1) / (ne1 * cos(alpha1) + ne2 * cos(alpha2));
	deviceComplex tp = 2. * ne1 * cos(alpha1) / (ne2 * cos(alpha1) + ne1 * cos(alpha2));
	if (isnan(ts.real()) || isnan(ts.imag())) ts = deviceComplex(0, 0);
	if (isnan(tp.real()) || isnan(tp.imag())) ts = deviceComplex(0, 0);
	(*s).gridEFrequency1[i] = ts * (*s).gridEFrequency1[i];
	(*s).gridEFrequency2[i] = tp * (*s).gridEFrequency2[i];
}


FGLOBAL void apertureKernel(GKERN cudaParameterSet* s, double radius, double activationParameter) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j, k, col;

	col = i / (*s).Ntime;
	j = col % (*s).Nspace;
	k = col / (*s).Nspace;
	double r;
	if ((*s).is3D) {
		double x = ((*s).dx * (j - (*s).Nspace / 2.0));
		double y = ((*s).dx * (k - (*s).Nspace2 / 2.0));
		r = sqrt(x * x + y * y);
	} 
	else {
		r = abs((*s).dx * ((double)j - (*s).Nspace / 2.0) - 0.25 * (*s).dx);
	}

	double a = 1.0 - (1.0 / (1.0 + exp( - activationParameter*(r - radius)/(*s).dx)));

	//if (r>radius) a = 0;
	(*s).gridETime1[i] *= a;
	(*s).gridETime2[i] *= a;
}

FGLOBAL void parabolicMirrorKernel(GKERN cudaParameterSet* s, double focus) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j, k, h, col;
	h = 1 + i % ((*s).Nfreq - 1);
	col = i / ((*s).Nfreq - 1);
	i = h + col * (*s).Nfreq;
	j = col % (*s).Nspace;
	k = col / (*s).Nspace;

	double w = TWOPI * h * (*s).fStep;
	double r;
	if ((*s).is3D) {
		double x = ((*s).dx * (j - (*s).Nspace / 2.0));
		double y = ((*s).dx * (k - (*s).Nspace2 / 2.0));
		r = sqrt(x * x + y * y);
	}
	else {
		r = abs((*s).dx * ((double)j - (*s).Nspace / 2.0) - 0.25 * (*s).dx);
	}


	

	deviceComplex	u = deviceLib::exp(deviceComplex(0.0,
		w * r * r * (0.5 / focus) / LIGHTC));


	(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
	(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
}

FGLOBAL void sphericalMirrorKernel(GKERN cudaParameterSet* s, double ROC) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j, k, h, col;
	h = 1 + i % ((*s).Nfreq - 1);
	col = i / ((*s).Nfreq - 1);
	i = h + col * (*s).Nfreq;
	j = col % (*s).Nspace;
	k = col / (*s).Nspace;

	double w = TWOPI * h * (*s).fStep;
	double r;
	if ((*s).is3D) {
		double x = ((*s).dx * (j - (*s).Nspace / 2.0));
		double y = ((*s).dx * (k - (*s).Nspace2 / 2.0));
		r = sqrt(x * x + y * y);
	}
	else {
		r = abs((*s).dx * ((double)j - (*s).Nspace / 2.0) - 0.25 * (*s).dx);
	}

	bool isNegative = signbit(ROC);
	ROC = abs(ROC);
	deviceComplex u = deviceComplex(0.0, 0.0);
	if (r < ROC) {
		u = deviceLib::exp(deviceComplex(0.0, 
			2.0 * pow(-1,isNegative)*w*ROC*((sqrt(1.0 - r * r / (ROC * ROC))) - 1.0)/LIGHTC));
	}

	(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
	(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
}

FGLOBAL void applyLinearPropagationKernel(GKERN double* sellmeierCoefficients, double thickness, cudaParameterSet *s) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j, h, k, col;
	int axesNumber = (*s).axesNumber;
	int sellmeierType = (*s).sellmeierType;
	deviceComplex ne, no, n0, n0o;
	deviceComplex cuZero = deviceComplex(0, 0);
	h = 1 + i % ((*s).Nfreq - 1);
	col = i / ((*s).Nfreq - 1);
	i = h + col * ((*s).Nfreq);
	j = col % (*s).Nspace;
	k = col / (*s).Nspace;
	deviceComplex ii = deviceComplex(0, 1);
	double crystalTheta = sellmeierCoefficients[66];
	double crystalPhi = sellmeierCoefficients[67];



	//frequency being resolved by current thread
	double f = h * (*s).fStep;
	double omega = TWOPI * f;
	findBirefringentCrystalIndex(s, sellmeierCoefficients, threadIdx.x + blockIdx.x * blockDim.x, &ne, &no);
	double dk1 = j * (*s).dk1 - (j >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace);
	double dk2 = k * (*s).dk2 - (k >= ((long long)(*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2);
	if (!(*s).is3D)dk2 = 0.0;
	//if ((*s).isCylindric) dk2 = dk1;
	sellmeierCuda(&n0, &n0o, sellmeierCoefficients, (*s).f0,
		crystalTheta, crystalPhi, axesNumber, sellmeierType);
	if (isnan(ne.real()) || isnan(no.real())) {
		ne = deviceComplex(1, 0);
		no = deviceComplex(1, 0);
	}

	deviceComplex ke = ne * omega / LIGHTC;
	deviceComplex ko = no * omega / LIGHTC;
	double k0 = (n0 * omega / LIGHTC).real();
	double kze = (deviceLib::sqrt(ke * ke - dk1 * dk1 - dk2 * dk2)).real();
	double kzo = (deviceLib::sqrt(ko * ko - dk1 * dk1 - dk2 * dk2)).real();

	deviceComplex ts = deviceLib::exp(ii * (k0 - kze) * thickness);
	deviceComplex tp = deviceLib::exp(ii * (k0 - kzo) * thickness);
	if (isnan(ts.real()) || isnan(ts.imag())) ts = deviceComplex(0, 0);
	if (isnan(tp.real()) || isnan(tp.imag())) tp = deviceComplex(0, 0);
	(*s).gridEFrequency1[i] = ts * (*s).gridEFrequency1[i];
	(*s).gridEFrequency2[i] = tp * (*s).gridEFrequency2[i];
}


//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
//note that the sellmeier coefficients have extra values appended to the end
//to give info about the current simulation
FGLOBAL void prepareCartesianGridsKernel(GKERN double* sellmeierCoefficients, cudaParameterSet* s) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j, k;
	int axesNumber = (*s).axesNumber;
	int sellmeierType = (*s).sellmeierType;
	deviceComplex ne, no, n0;
	deviceComplex cuZero = deviceComplex(0, 0);
	j = i / ((*s).Nfreq-1); //spatial coordinate
	k = 1 + (i % ((*s).Nfreq-1)); //temporal coordinate
	i = k + j * (*s).Nfreq;
	deviceComplex ii = deviceComplex(0, 1);
	double crystalTheta = sellmeierCoefficients[66];
	double crystalPhi = sellmeierCoefficients[67];
	double kStep = sellmeierCoefficients[70];
	double fStep = sellmeierCoefficients[71];

	//frequency being resolved by current thread
	double f = -k * fStep;

	//transverse wavevector being resolved
	double dk = j * kStep - (j >= ((long long)(*s).Nspace / 2)) * (kStep * (*s).Nspace); //frequency grid in transverse direction
	sellmeierCuda(&n0, &no, sellmeierCoefficients, abs((*s).f0),
		crystalTheta, crystalPhi, axesNumber, sellmeierType);
	findBirefringentCrystalIndex(s, sellmeierCoefficients, threadIdx.x + blockIdx.x * blockDim.x, &ne, &no);

	//walkoff angle has been found, generate the rest of the grids



	if (isnan(ne.real()) || isnan(no.real())) {
		ne = deviceComplex(1, 0);
		no = deviceComplex(1, 0);
	}

	deviceComplex k0 = deviceComplex(TWOPI * n0.real() * f / LIGHTC, 0);
	deviceComplex ke = TWOPI * ne * f / LIGHTC;
	deviceComplex ko = TWOPI * no * f / LIGHTC;


	deviceComplex chi11 = deviceComplex(1.0, 0);
	deviceComplex chi12 = deviceComplex(1.0, 0);
	if ((*s).isUsingMillersRule) {
		chi11 = (*s).chiLinear1[k];
		chi12 = (*s).chiLinear2[k];
	}
	else {
		chi11 = deviceComplex(1, 0);
		chi12 = deviceComplex(1, 0);
	}

	if (abs(dk) < deviceLib::abs(ke) && k < ((long long)(*s).Nfreq-1)) {
		(*s).gridPropagationFactor1[i] = ii * (ke - k0 - dk * dk / (2. * ke.real())) * (*s).h;
		if (isnan(((*s).gridPropagationFactor1[i]).real())) {
			(*s).gridPropagationFactor1[i] = cuZero;
		}

		(*s).gridPropagationFactor2[i] = ii * (ko - k0 - dk * dk / (2. * ko.real())) * (*s).h;
		if (isnan(((*s).gridPropagationFactor2[i]).real())) {
			(*s).gridPropagationFactor2[i] = cuZero;
		}

		(*s).gridPolarizationFactor1[i] = ii * chi11 * (TWOPI * f) / (2. *ne.real() * LIGHTC) * (*s).h;
		(*s).gridPolarizationFactor2[i] = ii * chi12 * (TWOPI * f) / (2. *no.real() * LIGHTC) * (*s).h;
	}

	else {
		(*s).gridPropagationFactor1[i] = cuZero;
		(*s).gridPropagationFactor2[i] = cuZero;
		(*s).gridPolarizationFactor1[i] = cuZero;
		(*s).gridPolarizationFactor2[i] = cuZero;
	}

}

//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
//note that the sellmeier coefficients have extra values appended to the end
//to give info about the current simulation
FGLOBAL void prepare3DGridsKernel(GKERN double* sellmeierCoefficients, cudaParameterSet* s) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long col,j, k, l;
	int axesNumber = (*s).axesNumber;
	int sellmeierType = (*s).sellmeierType;
	deviceComplex ne, no, n0;
	deviceComplex cuZero = deviceComplex(0, 0);
	col = i / ((*s).Nfreq-1); //spatial coordinate
	j = 1+i % ((*s).Nfreq-1); // frequency coordinate
	i = j + col * (*s).Nfreq;
	k = col % (*s).Nspace;
	l = col / (*s).Nspace;

	deviceComplex ii = deviceComplex(0, 1);
	double crystalTheta = sellmeierCoefficients[66];
	double crystalPhi = sellmeierCoefficients[67];

	//frequency being resolved by current thread
	double f = -j * (*s).fStep;

	//transverse wavevector being resolved
	double dk1 = k * (*s).dk1 -(k >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (long long)(*s).Nspace); //frequency grid in x direction
	double dk2 = l * (*s).dk2 - (l >= ((long long)(*s).Nspace2 / 2)) * ((*s).dk2 * (long long)(*s).Nspace2); //frequency grid in y direction
	sellmeierCuda(&n0, &no, sellmeierCoefficients, abs((*s).f0),
		crystalTheta, crystalPhi, axesNumber, sellmeierType);
	findBirefringentCrystalIndex(s, sellmeierCoefficients, threadIdx.x + blockIdx.x * blockDim.x, &ne, &no);



	if (isnan(ne.real()) || isnan(no.real())) {
		ne = deviceComplex(1, 0);
		no = deviceComplex(1, 0);
	}

	deviceComplex k0 = deviceComplex(TWOPI * n0.real() * f / LIGHTC, 0);
	deviceComplex ke = TWOPI * ne * f / LIGHTC;
	deviceComplex ko = TWOPI * no * f / LIGHTC;


	deviceComplex chi11 = deviceComplex(1.0, 0);
	deviceComplex chi12 = deviceComplex(1.0, 0);
	if ((*s).isUsingMillersRule) {
		chi11 = (*s).chiLinear1[j];
		chi12 = (*s).chiLinear2[j];
	}
	else {
		chi11 = deviceComplex(1, 0);
		chi12 = deviceComplex(1, 0);
	}

	if (max(abs(dk1),abs(dk2)) < deviceLib::abs(ke) && j < ((long long)(*s).Nfreq-1)) {
		(*s).gridPropagationFactor1[i] = ii * (ke - k0 - (dk1 * dk1 + dk2 * dk2) / (2. * ke.real())) * (*s).h;
		if (isnan(((*s).gridPropagationFactor1[i].real()))) {
			(*s).gridPropagationFactor1[i] = cuZero;
		}

		(*s).gridPropagationFactor2[i] = ii * (ko - k0 - (dk1 * dk1 + dk2 * dk2) / (2. * ko.real())) * (*s).h;
		if (isnan(((*s).gridPropagationFactor2[i].real()))) {
			(*s).gridPropagationFactor2[i] = cuZero;
		}

		(*s).gridPolarizationFactor1[i] = ii * chi11 * (TWOPI * f) / (2. *ne.real() * LIGHTC) * (*s).h;
		(*s).gridPolarizationFactor2[i] = ii * chi12 * (TWOPI * f) / (2. *no.real() * LIGHTC) * (*s).h;
	}

	else {
		(*s).gridPropagationFactor1[i] = cuZero;
		(*s).gridPropagationFactor2[i] = cuZero;
		(*s).gridPolarizationFactor1[i] = cuZero;
		(*s).gridPolarizationFactor2[i] = cuZero;
	}

}

FGLOBAL void getChiLinearKernel(GKERN cudaParameterSet* s, double* sellmeierCoefficients) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	int axesNumber = (*s).axesNumber;
	int sellmeierType = (*s).sellmeierType;
	deviceComplex cuZero = deviceComplex(0, 0);


	double crystalTheta = sellmeierCoefficients[66];
	double crystalPhi = sellmeierCoefficients[67];
	double fStep = sellmeierCoefficients[71];

	deviceComplex ne, no, n0;

	//frequency being resolved by current thread
	double f = i * fStep;
	sellmeierCuda(&n0, &no, sellmeierCoefficients, abs((*s).f0), crystalTheta, crystalPhi, axesNumber, sellmeierType);
	sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta, crystalPhi, axesNumber, sellmeierType);
	if (isnan(ne.real()) || isnan(no.real())) {
		ne = deviceComplex(1, 0);
		no = deviceComplex(1, 0);
	}


	(*s).chiLinear1[i] = -1. + ne * ne;
	(*s).chiLinear2[i] = -1. + no * no;
	if ((((*s).chiLinear1[i].real()) == 0) || (((*s).chiLinear2[i].real()) == 0) || isnan(((*s).chiLinear1[i].real())) || isnan(((*s).chiLinear2[i].real()))) {
		(*s).chiLinear1[i] = deviceComplex(1, 0);
		(*s).chiLinear2[i] = deviceComplex(1, 0);
	}
	(*s).inverseChiLinear1[i] = 1.0 / (*s).chiLinear1[i].real();
	(*s).inverseChiLinear2[i] = 1.0 / (*s).chiLinear2[i].real();
}
//prepare the propagation constants under the assumption of cylindrical symmetry of the beam
FGLOBAL void prepareCylindricGridsKernel(GKERN double* sellmeierCoefficients, cudaParameterSet* s) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j, k;
	long long Nspace = (*s).Nspace;
	int axesNumber = (*s).axesNumber;
	int sellmeierType = (*s).sellmeierType;
	deviceComplex cuZero = deviceComplex(0, 0);
	j = i / ((*s).Nfreq-1); //spatial coordinate
	k = 1 + i % ((*s).Nfreq-1); //temporal coordinate
	i = k + j * (*s).Nfreq;


	deviceComplex ii = deviceComplex(0, 1);
	double crystalTheta = sellmeierCoefficients[66];
	double crystalPhi = sellmeierCoefficients[67];
	double kStep = sellmeierCoefficients[70];
	double fStep = sellmeierCoefficients[71];

	deviceComplex ne, no, n0;

	//frequency being resolved by current thread
	double f = -k * fStep;

	//transverse wavevector being resolved
	double dk = j * kStep - (j >= (Nspace / 2)) * (kStep * Nspace); //frequency grid in transverse direction
	sellmeierCuda(&n0, &no, sellmeierCoefficients, abs((*s).f0), crystalTheta, crystalPhi, axesNumber, sellmeierType);
	sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta, crystalPhi, axesNumber, sellmeierType);
	if (isnan(ne.real()) || isnan(no.real())) {
		ne = deviceComplex(1, 0);
		no = deviceComplex(1, 0);
	}

	deviceComplex k0 = deviceComplex(TWOPI * n0.real() * f / LIGHTC, 0);
	deviceComplex ke = TWOPI * ne * f / LIGHTC;
	deviceComplex ko = TWOPI * no * f / LIGHTC;

	deviceComplex chi11 = (*s).chiLinear1[k];
	deviceComplex chi12 = (*s).chiLinear2[k];
	if (!(*s).isUsingMillersRule) {
		chi11 = deviceComplex(1, 0);
		chi12 = deviceComplex(1, 0);
	}

	if (abs(dk) <= min(deviceLib::abs(ke), deviceLib::abs(ko)) && k<((long long)(*s).Nfreq-1)) {
		(*s).gridPropagationFactor1[i] = ii * (ke - k0 - dk * dk / (2. * ke.real())) * (*s).h;
		(*s).gridPropagationFactor1Rho1[i] = ii * (1. / (chi11 * 2. * ke.real())) * (*s).h;
		if (isnan(((*s).gridPropagationFactor1[i].real()))) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor1Rho1[i] = cuZero;
		}

		(*s).gridPropagationFactor2[i] = ii * (ko - k0 - dk * dk / (2. * ko.real())) * (*s).h;
		(*s).gridPropagationFactor1Rho2[i] = ii * (1. / (chi12 * 2. * ko.real())) * (*s).h;
		if (isnan(((*s).gridPropagationFactor2[i].real()))) {
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPropagationFactor1Rho2[i] = cuZero;
		}
		//factor of 0.5 comes from doubled grid size in cylindrical symmetry mode after expanding the beam
		(*s).gridPolarizationFactor1[i] = 0.5 * chi11 * ii * (TWOPI * f) / (2. *ne.real() * LIGHTC) * (*s).h;
		(*s).gridPolarizationFactor2[i] = 0.5 * chi12 * ii * (TWOPI * f) / (2. *no.real() * LIGHTC) * (*s).h;


	}

	else {
		(*s).gridPropagationFactor1[i] = cuZero;
		(*s).gridPropagationFactor2[i] = cuZero;
		(*s).gridPolarizationFactor1[i] = cuZero;
		(*s).gridPolarizationFactor2[i] = cuZero;
		(*s).gridPropagationFactor1[i] = cuZero;
		(*s).gridPropagationFactor1Rho2[i] = cuZero;
	}
}

//replaces E with its complex conjugate
FGLOBAL void conjugateKernel(GKERN deviceComplex* E) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	E[i] = deviceLib::conj(E[i]);
}

FGLOBAL void realToComplexKernel(GKERN double* in, deviceComplex* out) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	out[i] = deviceComplex(in[i], 0.0);
}

FGLOBAL void complexToRealKernel(GKERN deviceComplex* in, double* out) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	out[i] = in[i].real();
}

FGLOBAL void materialPhaseKernel(GKERN double df, size_t Ntime, double* a, double f01, double f02, 
	double thickness1, double thickness2, double* phase1, double* phase2) {
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;
	//frequency being resolved by current thread
	double f = i * df;
	if (i >= Ntime / 2) {
		f -= df * Ntime;
	}

	//give phase shift relative to group velocity (approximated 
	// with low-order finite difference) so the pulse doesn't move
	deviceComplex ne, no, no0, n0p, n0m;
	sellmeierCuda(&ne, &no, a, abs(f), 0, 0, 0, 0);
	f *= -TWOPI;
	sellmeierCuda(&ne, &no0, a, f01, 0, 0, 0, 0);
	sellmeierCuda(&ne, &n0p, a, f01 + 1e11, 0, 0, 0, 0);
	sellmeierCuda(&ne, &n0m, a, f01 - 1e11, 0, 0, 0, 0);
	no0 = no0 + f01 * (n0p - n0m) / 2e11;
	phase1[i] = thickness1 * f * (no.real() - no0.real()) / LIGHTC;
	sellmeierCuda(&ne, &no0, a, f02, 0, 0, 0, 0);
	sellmeierCuda(&ne, &n0p, a, f02 + 1e11, 0, 0, 0, 0);
	sellmeierCuda(&ne, &n0m, a, f02 - 1e11, 0, 0, 0, 0);
	no0 = no0 + f02 * (n0p - n0m) / 2e11;
	phase2[i] = thickness2 * f * (no.real() - no0.real()) / LIGHTC;

}
//replaces NaN values with 0
// don't use this in a final version of anything, but can be useful for debugging
// Comment it out when you're done.
//FGLOBAL void fixnanKernel(GKERN deviceComplex* E) {
//	long long i = threadIdx.x + blockIdx.x * blockDim.x;
//	if (isnan(E[i].real()) || isnan(E[i].imag())) {
//		E[i] = deviceComplex(0., 0.);
//	}
//}

//calculate the nonlinear polarization, after FFT to get the field
//in the time domain
FGLOBAL void nonlinearPolarizationKernel(GKERN cudaParameterSet* s) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	double Ex = (*s).fftNorm * (*s).gridETime1[i];
	double Ey = (*s).fftNorm * (*s).gridETime2[i];

	double Ex2 = Ex * Ex;
	double Ey2 = Ey * Ey;
	(*s).gridPolarizationTime1[i] = 0.;
	(*s).gridPolarizationTime2[i] = 0.;

	//The d2eff tensor has the form
	// | d_xxx d_xyx d_yyx |
	// | d_xxy d_xyy d_yyy |
	if ((*s).nonlinearSwitches[0] == 1) {
		(*s).gridPolarizationTime1[i] += (*s).chi2Tensor[0] * Ex2 + (*s).chi2Tensor[2] * Ex * Ey + (*s).chi2Tensor[4] * Ey2;
		(*s).gridPolarizationTime2[i] += (*s).chi2Tensor[1] * Ex2 + (*s).chi2Tensor[3] * Ex * Ey + (*s).chi2Tensor[5] * Ey2;
	}

	//resolve the full chi3 matrix when (*s).nonlinearSwitches[1]==1
	if ((*s).nonlinearSwitches[1] == 1) {

		//rotate field into crystal frame
		double E3[3] = { (*s).rotationForward[0] * Ex + (*s).rotationForward[1] * Ey,
			(*s).rotationForward[3] * Ex + (*s).rotationForward[4] * Ey,
			(*s).rotationForward[6] * Ex + (*s).rotationForward[7] * Ey };

		//loop over tensor element X_abcd
		//i hope the compiler unrolls this, but no way am I writing that out by hand
		unsigned char a, b, c, d;
		double P3[3] = { 0 };
		for (a = 0; a < 3; a++) {
			for (b = 0; b < 3; b++) {
				for (c = 0; c < 3; c++) {
					for (d = 0; d < 3; d++) {
						P3[d] += (*s).chi3Tensor[a + 3 * b + 9 * c + 27 * d] * E3[a] * E3[b] * E3[c];
					}
				}
			}
		}

		//rotate back into simulation frame
		(*s).gridPolarizationTime1[i] += (*s).rotationBackward[0] * P3[0] + (*s).rotationBackward[1] * P3[1] + (*s).rotationBackward[2] * P3[2];
		(*s).gridPolarizationTime2[i] += (*s).rotationBackward[3] * P3[0] + (*s).rotationBackward[4] * P3[1] + (*s).rotationBackward[5] * P3[2];
	}
	//using only one value of chi3, under assumption of centrosymmetry
	if ((*s).nonlinearSwitches[1] == 2) {
		double Esquared = (*s).chi3Tensor[0] * (Ex2 + Ey2);
		(*s).gridPolarizationTime1[i] += Ex * Esquared;
		(*s).gridPolarizationTime2[i] += Ey * Esquared;
	}
}


//Plasma response with time-dependent carrier density
//This polarization needs a different factor in the nonlinear wave equation
//to account for the integration
//plasmaParameters vector:
// 0    e^2/m_eff
// 1    gamma_drude
// 2    ionization rate/E^N
// 3    absorption strength
//equation for the plasma current:
//J_drude(t) = (e/m)*exp(-gamma*t)*\int_-infty^t dt' exp(gamma*t)*N(t)*E(t)
//J_absorption(t) = beta*E^(2*Nphot-2)*E
//plasmaParameters[0] is the nonlinear absorption parameter
	//nonlinearSwitches[3] is Nphotons-2
	//plasmaParameters[2] is the 1/photon energy, translating the loss of power
	//from the field to the number of free carriers
	//extra factor of (dt^2e^2/(m*photon energy*eo) included as it is needed for the amplitude
	//of the plasma current
FGLOBAL void plasmaCurrentKernel_twoStage_A(GKERN cudaParameterSet* s) {
	unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
	double Esquared, Ex, Ey, a;
	unsigned char pMax = (unsigned char)(*s).nonlinearSwitches[3];
	Ex = (*s).gridETime1[i] * (*s).fftNorm;
	Ey = (*s).gridETime2[i] * (*s).fftNorm;

	//save values in workspaces, casting to double
	double* dN = (double*)(*s).workspace1;
	double* dN2 = dN + (*s).Ngrid;
	double* Jx = (*s).gridPolarizationTime1;
	double* Jy = (*s).gridPolarizationTime2;
	Esquared = Ex * Ex + Ey * Ey;
	a = (*s).plasmaParameters[0] * Esquared;
	for (unsigned char p = 0; p < pMax; p++) {
		a *= Esquared;
	}
	Jx[i] = a * Ex;
	Jy[i] = a * Ey;
	dN[i] = (*s).plasmaParameters[2] * (Jx[i] * Ex + Jy[i] * Ey);
	dN2[i] = dN[i];
}

FGLOBAL void plasmaCurrentKernel_twoStage_B(GKERN cudaParameterSet* s) {
	unsigned int j = threadIdx.x + blockIdx.x * blockDim.x;
	j *= (unsigned int)(*s).Ntime;
	double N = 0;
	double integralx = 0;
	double* expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
	double Ex, a;
	double* dN = j + (double*)(*s).workspace1;
	double* E = &(*s).gridETime1[j];
	double* P = &(*s).gridPolarizationTime1[j];
	for (unsigned int k = 0; k < (*s).Ntime; k++) {
		Ex = E[k] * (*s).fftNorm;
		N += dN[k];
		a = N * (*s).expGammaT[k];
		integralx += a * Ex;
		P[k] = P[k] + expMinusGammaT[k] * integralx;
	}
}

FGLOBAL void updateKwithPolarizationKernel(GKERN cudaParameterSet* sP) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
	long long j = i / ((*sP).Nfreq - 1); //spatial coordinate
	i = h + j * ((*sP).Nfreq);
	h += (j + ((*sP).isCylindric * (j > ((long long)(*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;

	(*sP).k1[i] += (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h];
	(*sP).k2[i] += (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h];
}

FGLOBAL void updateKwithPlasmaKernel(GKERN cudaParameterSet* sP) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
	long long j = i / ((*sP).Nfreq - 1); //spatial coordinate
	i = h + j * ((*sP).Nfreq);

	deviceComplex jfac = deviceComplex(0, -1.0 / (h * (*sP).fStep));
	h += (j + ((*sP).isCylindric * (j > ((long long)(*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;

	(*sP).k1[i] += jfac * (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h] * (*sP).inverseChiLinear1[i % ((*sP).Nfreq)];
	(*sP).k2[i] += jfac * (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h] * (*sP).inverseChiLinear2[i % ((*sP).Nfreq)];
}

//Main kernel for RK4 propagation of the field
FGLOBAL void rkKernel(GKERN cudaParameterSet* sP, uint8_t stepNumber) {
	unsigned int iC = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate

	iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
	if (h == 1) {
		(*sP).k1[iC - 1] = deviceComplex(0., 0.);
		(*sP).k2[iC - 1] = deviceComplex(0., 0.);
		(*sP).gridEFrequency1[iC - 1] = deviceComplex(0., 0.);
		(*sP).gridEFrequency2[iC - 1] = deviceComplex(0., 0.);
		(*sP).gridEFrequency1Next1[iC - 1] = deviceComplex(0., 0.);
		(*sP).gridEFrequency1Next2[iC - 1] = deviceComplex(0., 0.);
		(*sP).workspace1[iC - 1] = deviceComplex(0., 0.);
		(*sP).workspace2[iC - 1] = deviceComplex(0., 0.);
	}
	deviceComplex estimate1, estimate2;

	if ((*sP).isCylindric) {
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		(*sP).k2[iC] = (*sP).k2[iC] + (*sP).gridPropagationFactor1Rho2[iC] * (*sP).workspace2[iC];
	}

	//generate the estimates and do the weighted sum to get the grid at the next step
	//with weights determined by the step number
	switch (stepNumber) {
	case 0:
		estimate1 = (*sP).gridEFrequency1[iC] + 0.5 * (*sP).k1[iC];
		estimate2 = (*sP).gridEFrequency2[iC] + 0.5 * (*sP).k2[iC];
		(*sP).gridEFrequency1Next1[iC] = SIXTH * (*sP).k1[iC] + (*sP).gridEFrequency1[iC];
		(*sP).gridEFrequency1Next2[iC] = SIXTH * (*sP).k2[iC] + (*sP).gridEFrequency2[iC];
		if ((*sP).isUsingMillersRule) {
			(*sP).workspace1[iC] = (*sP).chiLinear1[h] * estimate1;
			(*sP).workspace2[iC] = (*sP).chiLinear2[h] * estimate2;
		}
		else {
			(*sP).workspace1[iC] = estimate1;
			(*sP).workspace2[iC] = estimate2;
		}
		(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * estimate1;
		(*sP).k2[iC] = (*sP).gridPropagationFactor2[iC] * estimate2;
		break;
	case 1:
		estimate1 = (*sP).gridEFrequency1[iC] + 0.5 * (*sP).k1[iC];
		estimate2 = (*sP).gridEFrequency2[iC] + 0.5 * (*sP).k2[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + THIRD * (*sP).k1[iC];
		(*sP).gridEFrequency1Next2[iC] = (*sP).gridEFrequency1Next2[iC] + THIRD * (*sP).k2[iC];
		if ((*sP).isUsingMillersRule) {
			(*sP).workspace1[iC] = (*sP).chiLinear1[h] * estimate1;
			(*sP).workspace2[iC] = (*sP).chiLinear2[h] * estimate2;
		}
		else {
			(*sP).workspace1[iC] = estimate1;
			(*sP).workspace2[iC] = estimate2;
		}
		(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * estimate1;
		(*sP).k2[iC] = (*sP).gridPropagationFactor2[iC] * estimate2;
		break;
	case 2:
		estimate1 = (*sP).gridEFrequency1[iC] + (*sP).k1[iC];
		estimate2 = (*sP).gridEFrequency2[iC] + (*sP).k2[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + THIRD * (*sP).k1[iC];
		(*sP).gridEFrequency1Next2[iC] = (*sP).gridEFrequency1Next2[iC] + THIRD * (*sP).k2[iC];
		if ((*sP).isUsingMillersRule) {
			(*sP).workspace1[iC] = (*sP).chiLinear1[h] * estimate1;
			(*sP).workspace2[iC] = (*sP).chiLinear2[h] * estimate2;
		}
		else {
			(*sP).workspace1[iC] = estimate1;
			(*sP).workspace2[iC] = estimate2;
		}
		(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * estimate1;
		(*sP).k2[iC] = (*sP).gridPropagationFactor2[iC] * estimate2;
		break;
	case 3:
		(*sP).gridEFrequency1[iC] = (*sP).gridEFrequency1Next1[iC] + SIXTH * (*sP).k1[iC];
		(*sP).gridEFrequency2[iC] = (*sP).gridEFrequency1Next2[iC] + SIXTH * (*sP).k2[iC];
		if ((*sP).isUsingMillersRule) {
			(*sP).workspace1[iC] = (*sP).chiLinear1[h] * (*sP).gridEFrequency1[iC];
			(*sP).workspace2[iC] = (*sP).chiLinear2[h] * (*sP).gridEFrequency2[iC];
		}
		else {
			(*sP).workspace1[iC] = (*sP).gridEFrequency1[iC];
			(*sP).workspace2[iC] = (*sP).gridEFrequency2[iC];
		}
		(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC];
		(*sP).k2[iC] = (*sP).gridPropagationFactor2[iC] * (*sP).gridEFrequency2[iC];
		break;
	}
}

FGLOBAL void beamNormalizeKernel(GKERN cudaParameterSet* s, double* rawSum, double* pulse, double pulseEnergy) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	double normFactor = sqrt(pulseEnergy / ((*s).Ntime * (*rawSum)));
	pulse[i] *= normFactor;
}

FGLOBAL void addDoubleArraysKernel(GKERN double* A, double* B) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	A[i] += B[i];
}

FGLOBAL void beamGenerationKernel2D(GKERN deviceComplex* pulse, double* pulseSum, cudaParameterSet* s, double frequency, double bandwidth,
	int sgOrder, double cep, double delay, double gdd, double tod,
	bool hasLoadedField, deviceComplex* loadedField, double* materialPhase,
	double w0, double z0, double x0, double beamAngle,
	double polarizationAngle, double circularity,
	double* sellmeierCoefficients, double crystalTheta, double crystalPhi, int sellmeierType) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j, h;
	h = 1 + i % ((*s).Nfreq - 1);
	j = i / ((*s).Nfreq - 1);
	i = h + j * ((*s).Nfreq);
	double f = h * (*s).fStep;
	double w = TWOPI * (f - frequency);

	//supergaussian pulse spectrum, if no input pulse specified
	deviceComplex specfac = deviceComplex(-pow((f - frequency) / bandwidth, sgOrder), 0);

	deviceComplex specphase = deviceComplex(0,
		-(cep
			+ TWOPI * f * (delay - 0.5 * (*s).dt * (*s).Ntime)
			+ 0.5 * gdd * w * w
			+ tod * w * w * w / 6.0
			+ materialPhase[h]));
	specfac = deviceLib::exp(specfac + specphase);

	if (hasLoadedField) {
		specfac = loadedField[h] * deviceLib::exp(specphase);
	}
	deviceComplex ne, no;
	sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta, crystalPhi, sellmeierType, 0);


	double ko = TWOPI * no.real() * f / LIGHTC;
	double zR = PI * w0 * w0 * ne.real() * f / LIGHTC;
	if (f == 0) {
		zR = 1e3;
	}
	double rB = (x0 - (*s).dx * (j - (*s).Nspace / 2.0) - 0.25 * (*s).dx);
	double r = rB * cos(beamAngle) - z0 * sin(beamAngle);
	double z = rB * sin(beamAngle) + z0 * cos(beamAngle);

	double wz = w0 * sqrt(1 + (z * z / (zR * zR)));
	double Rz = z * (1. + (zR * zR / (z * z)));

	if (z == 0) {
		Rz = 1.0e15;
	}
	double phi = atan(z / zR);
	deviceComplex Eb = (w0 / wz) * deviceLib::exp(deviceComplex(0., 1.) * (ko * (z - z0) + ko * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
	Eb = Eb * specfac;
	if (isnan(cuCModSquared(Eb)) || f <= 0) {
		Eb = deviceComplex(0., 0.);
	}

	pulse[i] = deviceComplex(cos(polarizationAngle), -circularity * sin(polarizationAngle)) * Eb;
	pulse[i + (*s).NgridC] = deviceComplex(sin(polarizationAngle), circularity * cos(polarizationAngle)) * Eb;
	double pointEnergy = abs(r) * (cuCModSquared(pulse[i]) + cuCModSquared(pulse[i + (*s).NgridC]));
	pointEnergy *= 2 * PI * LIGHTC * EPS0 * (*s).dx * (*s).dt;
	//two factors of two cancel here - there should be one for the missing frequency plane, but the sum is over x instead of r
	//accordingly we already multiplied by two
#ifdef __CUDACC__
	atomicAdd(pulseSum, pointEnergy);
#else
	*pulseSum += pointEnergy; //NOT THREAD SAFE, RUN CPU CODE ON SINGLE THREAD
#endif
}

FGLOBAL void beamGenerationKernel3D(GKERN deviceComplex* pulse, double* pulseSum, cudaParameterSet* s, double frequency, double bandwidth,
	int sgOrder, double cep, double delay, double gdd, double tod,
	bool hasLoadedField, deviceComplex* loadedField, double* materialPhase,
	double w0, double z0, double y0, double x0, double beamAngle, double beamAnglePhi,
	double polarizationAngle, double circularity,
	double* sellmeierCoefficients, double crystalTheta, double crystalPhi, int sellmeierType) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j, k, h, col;
	h = 1 + i % ((*s).Nfreq - 1);
	col = i / ((*s).Nfreq - 1);
	i = h + col * ((*s).Nfreq);
	j = col % (*s).Nspace;
	k = col / (*s).Nspace;
	double f = h * (*s).fStep;
	double w = TWOPI * (f - frequency);

	//supergaussian pulse spectrum, if no input pulse specified
	deviceComplex specfac = deviceComplex(-pow((f - frequency) / bandwidth, sgOrder), 0);

	deviceComplex specphase = deviceComplex(0,
		-(cep
			+ TWOPI * f * (delay - 0.5 * (*s).dt * (*s).Ntime)
			+ 0.5 * gdd * w * w
			+ tod * w * w * w / 6.0
			+ materialPhase[h]));
	specfac = deviceLib::exp(specfac + specphase);

	if (hasLoadedField) {
		specfac = loadedField[h] * deviceLib::exp(specphase);
	}
	deviceComplex ne, no;
	sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta, crystalPhi, sellmeierType, 0);


	double ko = TWOPI * no.real() * f / LIGHTC;
	double zR = PI * w0 * w0 * ne.real() * f / LIGHTC;
	if (f == 0) {
		zR = 1e3;
	}
	double xo = ((*s).dx * (j - (*s).Nspace / 2.0)) - x0;
	double yo = ((*s).dx * (k - (*s).Nspace2 / 2.0)) - y0;
	double zo = z0;
	double cB = cos(beamAngle);
	double cA = cos(beamAnglePhi);
	double sB = sin(beamAngle);
	double sA = sin(beamAnglePhi);
	double x = cB * xo + sA * sB * yo + sA * sB * zo;
	double y = cA * yo - sA * zo;
	double z = -sB * xo + sA * cB * yo + cA * cB * zo;
	double r = sqrt(x * x + y * y);

	double wz = w0 * sqrt(1 + (z * z / (zR * zR)));
	double Rz = 1.0e15;
	if (z != 0.0) {
		Rz = z * (1. + (zR * zR / (z * z)));
	}

	double phi = atan(z / zR);
	deviceComplex Eb = (w0 / wz) * deviceLib::exp(deviceComplex(0., 1.) * (ko * (z - z0) + ko * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
	Eb = Eb * specfac;
	if (isnan(cuCModSquared(Eb)) || f <= 0) {
		Eb = deviceComplex(0., 0.);
	}

	pulse[i] = deviceComplex(cos(polarizationAngle), -circularity * sin(polarizationAngle)) * Eb;
	pulse[i + (*s).NgridC] = deviceComplex(sin(polarizationAngle), circularity * cos(polarizationAngle)) * Eb;
	double pointEnergy = (cuCModSquared(pulse[i]) + cuCModSquared(pulse[i + (*s).NgridC]));
	pointEnergy *= 2 * LIGHTC * EPS0 * (*s).dx * (*s).dx * (*s).dt;


	//factor 2 accounts for the missing negative frequency plane
#ifdef __CUDACC__
	atomicAdd(pulseSum, pointEnergy);
#else
	*pulseSum += pointEnergy; //NOT THREAD SAFE, RUN CPU CODE ON SINGLE THREAD
#endif
}

//Take absolute value of complex array
FGLOBAL void absKernel(GKERN double* absOut, deviceComplex* complexIn) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	absOut[i] = deviceLib::abs(complexIn[i]);
}

//Apply fft normalization
//Take absolute value of complex array
FGLOBAL void fftNormalizeKernel(GKERN deviceComplex* A, size_t fftSize) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	A[i] = A[i] / (double)fftSize;
}

//Apply fft normalization
 FGLOBAL void multiplyByConstantKernel(GKERN deviceComplex* A, double val) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	A[i] = val * A[i];
}

FGLOBAL void multiplyByConstantKernelD(GKERN double* A, double val) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	A[i] = val * A[i];
}

//element-wise B*A = C;
FGLOBAL void multiplicationKernel(GKERN deviceComplex* A, deviceComplex* B, deviceComplex* C) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	C[i] = B[i] * A[i];
}

FGLOBAL void multiplicationKernelCompactVector(GKERN deviceComplex* A, deviceComplex* B, deviceComplex* C, cudaParameterSet* s) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long h = i % (*s).Nfreq; //temporal coordinate

	C[i] = A[h] * B[i];
}

FGLOBAL void multiplicationKernelCompact(GKERN deviceComplex* A, deviceComplex* B, deviceComplex* C) {
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	C[i] = A[i] * B[i];
}

//anonymous namespace helps to separate functions that have been compiled under nvcc from those
//compiled by the c++ compiler
namespace {
	simulationParameterSet* fittingSet;
	simulationParameterSet* fittingReferenceSet;
	int				runRK4Step(cudaParameterSet* sH, cudaParameterSet* sD, uint8_t stepNumber);
	int				preparePropagation2DCartesian(simulationParameterSet* s, cudaParameterSet sc);
	int				preparePropagation3DCylindric(simulationParameterSet* s, cudaParameterSet sc);
	int             preparePropagation3D(simulationParameterSet* s, cudaParameterSet sc);
	int             getTotalSpectrum(simulationParameterSet* sCPU, cudaParameterSet* sc);
	int				rotateField(simulationParameterSet* s, double rotationAngle);
	void            runFittingIteration(int* m, int* n, double* fittingValues, double* fittingFunction);
	int				prepareElectricFieldArrays(simulationParameterSet* s, cudaParameterSet* sc);
	int             applyLinearPropagation(simulationParameterSet* s, int materialIndex, double thickness);
	int             fillRotationMatricies(simulationParameterSet* sCPU, cudaParameterSet* s);
	int             deallocateCudaParameterSet(cudaParameterSet* s);
	int             initializeCudaParameterSet(simulationParameterSet* sCPU, cudaParameterSet* s);
	//generate the rotation matricies for translating between the beam coordinates and
	//crystal coordinates

//My weird bilingual wrapper template that lets me either call CUDA kernels
//normally on the GPU, or process them on the CPU
//This is why kernel declarations have FGLOBAL in front of them
//instead of the usual __global__ tag
	template<typename Function, typename... Args>
#ifdef __CUDACC__
	void bilingualLaunch(unsigned int Nblock, unsigned int Nthread, cudaStream_t stream, Function func, Args... args) {
#else
	void bilingualLaunch(unsigned int Nblock, unsigned int Nthread, int stream, Function func, Args... args) {
#endif
#ifdef __CUDACC__
		func <<<Nblock, Nthread, 0, stream >>> (args...);
#else
		uint3 bDim;
		bDim.x = Nthread;
		bool isThreaded = Nthread > 1u;
#pragma omp parallel for if(isThreaded)
		for (int j = 0; j < (int)Nblock; j++) {
			uint3 bIdx, tIdx;
			bIdx.x = (unsigned int)j;
			for (tIdx.x = 0u; tIdx.x < Nthread; tIdx.x++) {
				func(bIdx, tIdx, bDim, args...);
			}
		}
#endif
	}

	//memset for either CUDA or c++ depending on which
	//compiler is running
	int bilingualMemset(void* ptr, int value, size_t count) {
#ifdef __CUDACC__
		cudaMemset(ptr, value, count);
#else
		memset(ptr, value, count);
#endif
		return 0;
	}

	//calloc for either CUDA or c++ depending on which
	//compiler is running
	//in CUDA it's just malloc+memset to zero
	int bilingualCalloc(void** ptr, size_t N, size_t elementSize) {
#ifdef __CUDACC__
		int err = cudaMalloc(ptr, N * elementSize);
		bilingualMemset(*ptr, 0, N * elementSize);
		return err;
#else
		(*ptr) = calloc(N, elementSize);
		return (int)((*ptr) == NULL);
#endif
	}

	//free() or cudaFree() depending on what's compiling
	int bilingualFree(void* block) {
#ifdef __CUDACC__
		cudaFree(block);
#else
		free(block);
#endif
		return 0;
	}

	//memcpy() or cudaMemcpy depending on what's running
	int bilingualMemcpy(void* dst, void* src, size_t count, cudaMemcpyKind kind) {
#ifdef __CUDACC__
		cudaMemcpy(dst, src, count, kind);
#else
		memcpy(dst, src, count);
#endif
		return 0;
	}

	//FFT function that can be called from either CUDA or c++ code.
	//the properly initialized cudaParameterSet s will have either
	//cuFFT plans, or MKL dfti handles, for the following transform types:
	// type 0: double to complex (forward)
	// type 1: complex to double (backward)
	// type 2: 1D x spatial grid double to complex (forward, used in beam generation)
	// type 3: 1D x spatial grid complex to double (backward, used in beam generation)
	// type 4: double to complex with doubled spatial grid dimension (forward in cylindrical propagation)
	int combinedFFT(cudaParameterSet* s, void* input, void* output, int type) {
		type += RUNTYPE * 5;
		switch (type) {
#ifdef __CUDACC__
		case 0:
			cufftExecD2Z((*s).fftPlanD2Z, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
			break;
		case 1:
			cufftExecZ2D((*s).fftPlanZ2D, (cufftDoubleComplex*)input, (cufftDoubleReal*)output);
			break;
		case 2:
			cufftExecD2Z((*s).fftPlan1DD2Z, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
			break;
		case 3:
			cufftExecZ2D((*s).fftPlan1DZ2D, (cufftDoubleComplex*)input, (cufftDoubleReal*)output);
			break;
		case 4:
			cufftExecD2Z((*s).doublePolfftPlan, (cufftDoubleReal*)input, (cufftDoubleComplex*)output);
			break;
#endif
		case 5:
			DftiComputeForward((*s).mklPlanD2Z, input, output);
			break;
		case 6:
			DftiComputeBackward((*s).mklPlanZ2D, input, output);
			break;
		case 7:
			DftiComputeForward((*s).mklPlan1DD2Z, input, output);
			break;
		case 8:
			DftiComputeBackward((*s).mklPlan1DZ2D, input, output);
			break;
		case 9:
			DftiComputeForward((*s).mklPlanDoublePolfft, input, output);
			break;
		}

		return 0;
	}

	int prepareElectricFieldArrays(simulationParameterSet* s, cudaParameterSet* sc) {
		
		//run the beam generation single-threaded on CPU to avoid race condition
		unsigned int beamBlocks = (*sc).Nblock / 2;
		unsigned int beamThreads = (*sc).Nthread;
		if (RUNTYPE == 1) {
			beamBlocks = beamBlocks * beamThreads;
			beamThreads = 1;
		}
			
		cudaParameterSet* scDevice;
		bilingualCalloc((void**)&scDevice, 1, sizeof(cudaParameterSet));
		bilingualMemcpy(scDevice, sc, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);
		if ((*s).isFollowerInSequence && !(*s).isReinjecting) {
			bilingualMemcpy((*sc).gridETime1, (*s).ExtOut, 2 * (*s).Ngrid * sizeof(double), cudaMemcpyHostToDevice);
			combinedFFT(sc, (*sc).gridETime1, (*sc).gridEFrequency1, 0);
			//Copy the field into the temporary array
			bilingualMemcpy((*sc).gridEFrequency1Next1, (*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), cudaMemcpyDeviceToDevice);

			if ((*sc).isUsingMillersRule) {
				bilingualLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), (unsigned int)(2u * MIN_GRIDDIM), (*sc).CUDAStream, multiplicationKernelCompactVector, (*sc).chiLinear1, (*sc).gridEFrequency1Next1, (*sc).workspace1, scDevice);
			}
			else {
				bilingualMemcpy((*sc).workspace1, (*sc).gridEFrequency1Next1, 2 * sizeof(deviceComplex) * (*sc).NgridC, cudaMemcpyDeviceToDevice);
			}

			bilingualLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, (*sc).CUDAStream, multiplicationKernelCompact, (*sc).gridPropagationFactor1, (*sc).gridEFrequency1Next1, (*sc).k1);
			bilingualMemcpy((*sc).gridEFrequency1Next1, (*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), cudaMemcpyDeviceToDevice);
			bilingualFree(scDevice);
			return 0;
		}
		double* materialPhase1CUDA, * materialPhase2CUDA;
		deviceComplex* loadedField1, * loadedField2;

		bilingualCalloc((void**)&loadedField1, (*sc).Ntime, sizeof(deviceComplex));
		bilingualCalloc((void**)&loadedField2, (*sc).Ntime, sizeof(deviceComplex));

		//get the material phase
		double* materialCoefficientsCUDA, * sellmeierPropagationMedium;
		//NOTE TO SELF: add second phase material


		if ((*s).field1IsAllocated) {
			bilingualMemcpy(loadedField1, (*s).loadedField1, (*s).Ntime * sizeof(deviceComplex), cudaMemcpyHostToDevice);
		}
		if ((*s).field2IsAllocated) {
			bilingualMemcpy(loadedField2, (*s).loadedField2, (*s).Ntime * sizeof(deviceComplex), cudaMemcpyHostToDevice);
		}
		bilingualCalloc((void**)&materialCoefficientsCUDA, 66, sizeof(double));
		bilingualCalloc((void**)&sellmeierPropagationMedium, 66, sizeof(double));
		bilingualCalloc((void**)&materialPhase1CUDA, (*s).Ntime, sizeof(double));
		bilingualCalloc((void**)&materialPhase2CUDA, (*s).Ntime, sizeof(double));
		bilingualMemcpy(materialCoefficientsCUDA, (*s).crystalDatabase[(*s).phaseMaterialIndex1].sellmeierCoefficients, 66 * sizeof(double), cudaMemcpyHostToDevice);
		bilingualMemcpy(sellmeierPropagationMedium, (*s).crystalDatabase[(*s).materialIndex].sellmeierCoefficients, 66 * sizeof(double), cudaMemcpyHostToDevice);
		bilingualLaunch((unsigned int)(*s).Ntime, 1, (*sc).CUDAStream, materialPhaseKernel, (*s).fStep, (*s).Ntime, materialCoefficientsCUDA, (*s).frequency1, (*s).frequency2, (*s).phaseMaterialThickness1, (*s).phaseMaterialThickness2, materialPhase1CUDA, materialPhase2CUDA);

		double* pulseSum = &materialCoefficientsCUDA[0];
		//calculate pulse 1 and store it in unused memory
		bilingualMemset(pulseSum, 0, sizeof(double));
		bilingualMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		if ((*sc).is3D) {
			bilingualLaunch(beamBlocks, beamThreads, (*sc).CUDAStream, beamGenerationKernel3D,
				(*sc).workspace1, pulseSum, scDevice, (*s).frequency1, (*s).bandwidth1,
				(*s).sgOrder1, (*s).cephase1, (*s).delay1, (*s).gdd1, (*s).tod1,
				(*s).field1IsAllocated, loadedField1, materialPhase1CUDA, (*s).beamwaist1,
				(*s).z01, (*s).y01, (*s).x01, (*s).propagationAngle1, (*s).propagationAnglePhi1, (*s).polarizationAngle1, (*s).circularity1,
				sellmeierPropagationMedium, (*s).crystalTheta, (*s).crystalPhi, (*s).sellmeierType);
		}
		else {
			bilingualLaunch(beamBlocks, beamThreads, (*sc).CUDAStream, beamGenerationKernel2D,
				(*sc).workspace1, pulseSum, scDevice, (*s).frequency1, (*s).bandwidth1,
				(*s).sgOrder1, (*s).cephase1, (*s).delay1, (*s).gdd1, (*s).tod1,
				(*s).field1IsAllocated, loadedField1, materialPhase1CUDA, (*s).beamwaist1,
				(*s).z01, (*s).x01, (*s).propagationAngle1, (*s).polarizationAngle1, (*s).circularity1,
				sellmeierPropagationMedium, (*s).crystalTheta, (*s).crystalPhi, (*s).sellmeierType);
		}
		
		combinedFFT(sc, (*sc).workspace1, (*sc).gridETime1, 3);

		bilingualLaunch(2 * (*sc).Nblock, (*sc).Nthread, (*sc).CUDAStream, beamNormalizeKernel, scDevice, pulseSum, (*sc).gridETime1, (*s).pulseEnergy1);
		bilingualMemcpy((*sc).gridEFrequency1Next1, (*sc).gridETime1, (*sc).Ngrid * 2 * sizeof(double), cudaMemcpyDeviceToDevice);

		//calculate pulse 2
		bilingualMemset(pulseSum, 0, sizeof(double));
		bilingualMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		if ((*sc).is3D) {
			bilingualLaunch(beamBlocks, beamThreads, (*sc).CUDAStream, beamGenerationKernel3D,
				(*sc).workspace1, pulseSum, scDevice, (*s).frequency2, (*s).bandwidth2,
				(*s).sgOrder2, (*s).cephase2, (*s).delay2, (*s).gdd2, (*s).tod2,
				(*s).field2IsAllocated, loadedField2, materialPhase2CUDA, (*s).beamwaist2,
				(*s).z02, (*s).y02, (*s).x02, (*s).propagationAngle2, (*s).propagationAnglePhi2, (*s).polarizationAngle2, (*s).circularity2,
				sellmeierPropagationMedium, (*s).crystalTheta, (*s).crystalPhi, (*s).sellmeierType);
		}
		else {
			bilingualLaunch(beamBlocks, beamThreads, (*sc).CUDAStream, beamGenerationKernel2D,
				(*sc).workspace1, pulseSum, scDevice, (*s).frequency2, (*s).bandwidth2,
				(*s).sgOrder2, (*s).cephase2, (*s).delay2, (*s).gdd2, (*s).tod2,
				(*s).field2IsAllocated, loadedField2, materialPhase2CUDA, (*s).beamwaist2,
				(*s).z02, (*s).x02, (*s).propagationAngle2, (*s).polarizationAngle2, (*s).circularity2,
				sellmeierPropagationMedium, (*s).crystalTheta, (*s).crystalPhi, (*s).sellmeierType);
		}
		

		combinedFFT(sc, (*sc).workspace1, (*sc).gridETime1, 3);

		bilingualLaunch(2 * (*sc).Nblock, (*sc).Nthread, (*sc).CUDAStream, beamNormalizeKernel, scDevice, pulseSum, (*sc).gridETime1, (*s).pulseEnergy2);

		//add the pulses
		bilingualLaunch(2 * (*sc).Nblock, (*sc).Nthread, (*sc).CUDAStream, addDoubleArraysKernel, (*sc).gridETime1, (double*)(*sc).gridEFrequency1Next1);
		if ((*s).isReinjecting) {
			bilingualMemcpy((*sc).workspace1, (*s).ExtOut, 2 * (*s).Ngrid * sizeof(double), cudaMemcpyHostToDevice);
			bilingualLaunch(2 * (*sc).Nblock, (*sc).Nthread, (*sc).CUDAStream, addDoubleArraysKernel, (*sc).gridETime1, (double*)(*sc).workspace1);
		}
		//fft onto frequency grid
		combinedFFT(sc, (*sc).gridETime1, (*sc).gridEFrequency1, 0);

		//Copy the field into the temporary array
		bilingualMemcpy((*sc).gridEFrequency1Next1, (*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), cudaMemcpyDeviceToDevice);

		if ((*sc).isUsingMillersRule && !(*sc).forceLinear) {
			bilingualLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, (*sc).CUDAStream, multiplicationKernelCompactVector, (*sc).chiLinear1, (*sc).gridEFrequency1Next1, (*sc).workspace1, scDevice);
		}
		else {
			bilingualMemcpy((*sc).workspace1, (*sc).gridEFrequency1Next1, 2 * sizeof(deviceComplex) * (*sc).NgridC, cudaMemcpyDeviceToDevice);
		}

		bilingualLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, (*sc).CUDAStream, multiplicationKernelCompact, (*sc).gridPropagationFactor1, (*sc).gridEFrequency1Next1, (*sc).k1);
		bilingualMemcpy((*sc).gridEFrequency1Next1, (*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), cudaMemcpyDeviceToDevice);

		bilingualFree(materialPhase1CUDA);
		bilingualFree(materialPhase2CUDA);
		bilingualFree(materialCoefficientsCUDA);
		bilingualFree(sellmeierPropagationMedium);
		bilingualFree(loadedField1);
		bilingualFree(loadedField2);
		bilingualFree(scDevice);

		return 0;
	}
	int applyFresnelLoss(simulationParameterSet* s, int materialIndex1, int materialIndex2) {
		cudaParameterSet sc;
		initializeCudaParameterSet(s, &sc);
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).crystalDatabase[materialIndex1].sellmeierCoefficients, 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		double* sellmeierCoefficients1;
		double* sellmeierCoefficients2;
		bilingualCalloc((void**)&sellmeierCoefficients1, 74, sizeof(double));
		bilingualCalloc((void**)&sellmeierCoefficients2, 74, sizeof(double));
		bilingualMemcpy(sellmeierCoefficients1, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);

		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).crystalDatabase[materialIndex2].sellmeierCoefficients, 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		bilingualMemcpy(sellmeierCoefficients2, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);


		bilingualMemcpy(sc.gridEFrequency1, (*s).EkwOut, 2 * (*s).NgridC * sizeof(deviceComplex), cudaMemcpyHostToDevice);

		
		//transform final result
		combinedFFT(&sc, sc.gridEFrequency1, sc.gridETime1, 1);
		bilingualLaunch(2 * sc.Nblock, sc.Nthread, sc.CUDAStream, multiplyByConstantKernelD, sc.gridETime1, 1.0 / sc.Ngrid);
		//copy the field arrays from the GPU to CPU memory
		bilingualMemcpy((*s).ExtOut, sc.gridETime1, 2 * (*s).Ngrid * sizeof(double), cudaMemcpyDeviceToHost);
		bilingualMemcpy((*s).EkwOut, sc.gridEFrequency1, 2 * (*s).Ngrid * sizeof(deviceComplex), cudaMemcpyDeviceToHost);

		bilingualFree(sellmeierCoefficients1);
		bilingualFree(sellmeierCoefficients2);
		deallocateCudaParameterSet(&sc);
		return 0;
	}

	int applyAperature(simulationParameterSet* sCPU, double diameter, double activationParameter) {
		cudaParameterSet s;
		initializeCudaParameterSet(sCPU, &s);
		bilingualMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), cudaMemcpyHostToDevice);

		cudaParameterSet* sDevice;
		bilingualCalloc((void**)&sDevice, 1, sizeof(cudaParameterSet));
		bilingualMemcpy(sDevice, &s, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);
		bilingualLaunch(s.Nblock, s.Nthread, s.CUDAStream, apertureKernel, sDevice, 0.5 * diameter, activationParameter);
		combinedFFT(&s, s.gridETime1, s.gridEFrequency1, 0);
		bilingualMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), cudaMemcpyDeviceToHost);
		bilingualMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), cudaMemcpyDeviceToHost);
		getTotalSpectrum(sCPU, &s);
		deallocateCudaParameterSet(&s);
		bilingualFree(sDevice);
		return 0;
	}

	int applySphericalMirror(simulationParameterSet* sCPU, double ROC) {
		cudaParameterSet s;
		initializeCudaParameterSet(sCPU, &s);

		cudaParameterSet* sDevice;
		bilingualCalloc((void**)&sDevice, 1, sizeof(cudaParameterSet));
		bilingualMemcpy(sDevice, &s, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);

		bilingualMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), cudaMemcpyHostToDevice);
		combinedFFT(&s, s.gridETime1, s.gridEFrequency1, 2);
		bilingualLaunch(s.Nblock / 2, s.Nthread, s.CUDAStream, sphericalMirrorKernel, sDevice, ROC);
		combinedFFT(&s, s.gridEFrequency1, s.gridETime1, 3);
		bilingualLaunch(2 * s.Nblock, s.Nthread, s.CUDAStream, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ntime);
		combinedFFT(&s, s.gridETime1, s.gridEFrequency1, 0);
		bilingualMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), cudaMemcpyDeviceToHost);
		bilingualMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), cudaMemcpyDeviceToHost);
		getTotalSpectrum(sCPU, &s);
		deallocateCudaParameterSet(&s);
		bilingualFree(sDevice);

		return 0;
	}

	int applyParabolicMirror(simulationParameterSet* sCPU, double focus) {
		cudaParameterSet s;
		initializeCudaParameterSet(sCPU, &s);

		cudaParameterSet* sDevice;
		bilingualCalloc((void**)&sDevice, 1, sizeof(cudaParameterSet));
		bilingualMemcpy(sDevice, &s, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);

		bilingualMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), cudaMemcpyHostToDevice);
		combinedFFT(&s, s.gridETime1, s.gridEFrequency1, 2);
		bilingualLaunch(s.Nblock / 2, s.Nthread, s.CUDAStream, parabolicMirrorKernel, sDevice, focus);
		combinedFFT(&s, s.gridEFrequency1, s.gridETime1, 3);
		bilingualLaunch(2 * s.Nblock, s.Nthread, s.CUDAStream, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ntime);
		combinedFFT(&s, s.gridETime1, s.gridEFrequency1, 0);
		bilingualMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), cudaMemcpyDeviceToHost);
		bilingualMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), cudaMemcpyDeviceToHost);
		getTotalSpectrum(sCPU, &s);
		deallocateCudaParameterSet(&s);
		bilingualFree(sDevice);
		return 0;
	}

	int applyLinearPropagation(simulationParameterSet* sCPU, int materialIndex, double thickness) {
		cudaParameterSet s;
		initializeCudaParameterSet(sCPU, &s);


		bilingualMemcpy(s.gridEFrequency1, (*sCPU).EkwOut, s.NgridC * 2 * sizeof(deviceComplex), cudaMemcpyHostToDevice);



		double* sellmeierCoefficients = (double*)s.gridEFrequency1Next1;
		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(sellmeierCoefficientsAugmentedCPU, (*sCPU).crystalDatabase[materialIndex].sellmeierCoefficients, 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*sCPU).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*sCPU).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*sCPU).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*sCPU).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*sCPU).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*sCPU).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		bilingualMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);
		s.axesNumber = (*sCPU).crystalDatabase[materialIndex].axisType;
		s.sellmeierType = (*sCPU).crystalDatabase[materialIndex].sellmeierType;
		cudaParameterSet* sDevice;
		bilingualCalloc((void**)&sDevice, 1, sizeof(cudaParameterSet));
		bilingualMemcpy(sDevice, &s, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);



		bilingualLaunch(s.Nblock / 2, s.Nthread, s.CUDAStream, applyLinearPropagationKernel, sellmeierCoefficients, thickness, sDevice);
		bilingualMemcpy((*sCPU).EkwOut, s.gridEFrequency1, s.NgridC * 2 * sizeof(deviceComplex), cudaMemcpyDeviceToHost);
		combinedFFT(&s, s.gridEFrequency1, s.gridETime1, 1);
		bilingualLaunch(2 * s.Nblock, s.Nthread, s.CUDAStream, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);

		bilingualMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), cudaMemcpyDeviceToHost);

		deallocateCudaParameterSet(&s);
		bilingualFree(sDevice);
		return 0;
	}

	int preparePropagation2DCartesian(simulationParameterSet* s, cudaParameterSet sc) {
		//recycle allocated device memory for the grids needed
		double* sellmeierCoefficients = (double*)sc.gridEFrequency1Next1;

		double* referenceFrequencies;
		bilingualCalloc((void**)&referenceFrequencies, 7, sizeof(double));
		bilingualMemcpy(referenceFrequencies, (*s).crystalDatabase[(*s).materialIndex].nonlinearReferenceFrequencies, 7 * sizeof(double), cudaMemcpyHostToDevice);

		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).sellmeierCoefficients, 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		bilingualMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);

		//prepare the propagation grids
		cudaParameterSet* sD;
		bilingualCalloc((void**)&sD, 1, sizeof(cudaParameterSet));
		bilingualMemcpy(sD, &sc, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);
		bilingualLaunch((unsigned int)sc.Nfreq, 1, sc.CUDAStream, getChiLinearKernel, sD, sellmeierCoefficients);
		bilingualLaunch(sc.Nblock / 2, sc.Nthread, sc.CUDAStream, prepareCartesianGridsKernel, sellmeierCoefficients, sD);
		bilingualLaunch(1, 1, sc.CUDAStream, millersRuleNormalizationKernel, sD, sellmeierCoefficients, referenceFrequencies);

		bilingualFree(sD);

		//clean up
		bilingualMemset(sc.gridEFrequency1Next1, 0, 2 * (*s).NgridC * sizeof(deviceComplex));

		bilingualFree(referenceFrequencies);
		return 0;
	}



	int preparePropagation3D(simulationParameterSet* s, cudaParameterSet sc) {
		//recycle allocated device memory for the grids needed
		double* sellmeierCoefficients = (double*)sc.gridEFrequency1Next1;

		double* referenceFrequencies;
		bilingualCalloc((void**)&referenceFrequencies, 7, sizeof(double));
		bilingualMemcpy(referenceFrequencies, (*s).crystalDatabase[(*s).materialIndex].nonlinearReferenceFrequencies, 7 * sizeof(double), cudaMemcpyHostToDevice);

		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).sellmeierCoefficients, 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		bilingualMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);

		//prepare the propagation grids
		cudaParameterSet* sD;
		bilingualCalloc((void**)&sD, 1, sizeof(cudaParameterSet));
		bilingualMemcpy(sD, &sc, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);

		bilingualLaunch((unsigned int)sc.Nfreq, 1u, sc.CUDAStream, getChiLinearKernel, sD, sellmeierCoefficients);
		bilingualLaunch((unsigned int)sc.Nblock / 2u, (unsigned int)sc.Nthread, sc.CUDAStream, prepare3DGridsKernel, sellmeierCoefficients, sD);
		bilingualLaunch(1, 1, sc.CUDAStream, millersRuleNormalizationKernel, sD, sellmeierCoefficients, referenceFrequencies);
		bilingualFree(sD);

		//clean up
		bilingualMemset(sc.gridEFrequency1Next1, 0, 2 * (*s).NgridC * sizeof(deviceComplex));

		bilingualFree(referenceFrequencies);
		return 0;
	}

	int preparePropagation3DCylindric(simulationParameterSet* s, cudaParameterSet sc) {
		//recycle allocated device memory for the grids needed
		double* sellmeierCoefficients = (double*)sc.gridEFrequency1Next1;
		double* referenceFrequencies;
		bilingualCalloc((void**)&referenceFrequencies, 7, sizeof(double));
		bilingualMemcpy(referenceFrequencies, (*s).crystalDatabase[(*s).materialIndex].nonlinearReferenceFrequencies, 7 * sizeof(double), cudaMemcpyHostToDevice);

		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[74];
		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).sellmeierCoefficients, 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		bilingualMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);

		//prepare the propagation grids
		cudaParameterSet* sD;
		bilingualCalloc((void**)&sD, 1, sizeof(cudaParameterSet));
		bilingualMemcpy(sD, &sc, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);
		bilingualLaunch((unsigned int)sc.Nfreq, 1, sc.CUDAStream, getChiLinearKernel, sD, sellmeierCoefficients);
		bilingualLaunch((unsigned int)sc.Nblock / 2u, (unsigned int)sc.Nthread, sc.CUDAStream, prepareCylindricGridsKernel, sellmeierCoefficients, sD);
		bilingualLaunch(1, 1, sc.CUDAStream, millersRuleNormalizationKernel, sD, sellmeierCoefficients, referenceFrequencies);
		bilingualFree(sD);


		//clean up
		bilingualMemset(sc.gridEFrequency1Next1, 0, 2 * (*s).NgridC * sizeof(deviceComplex));
		bilingualFree(referenceFrequencies);
		return 0;
	}





	//Rotate the field on the GPU
	//Allocates memory and copies from CPU, then copies back to CPU and deallocates
	// - inefficient but the general principle is that only the CPU memory is preserved
	// after simulations finish... and this only runs at the end of the simulation
	int rotateField(simulationParameterSet* s, double rotationAngle) {
		cudaParameterSet sc;
		initializeCudaParameterSet(s, &sc);
		deviceComplex* Ein1, * Eout1, * Ein2, * Eout2;
		Ein1 = sc.gridEFrequency1;
		Ein2 = sc.gridEFrequency2;
		Eout1 = sc.gridEFrequency1Next1;
		Eout2 = sc.gridEFrequency1Next2;

		//retrieve/rotate the field from the CPU memory
		bilingualMemcpy(Ein1, (*s).EkwOut, 2 * (*s).NgridC * sizeof(deviceComplex), cudaMemcpyHostToDevice);
		bilingualLaunch((unsigned int)(sc.NgridC / MIN_GRIDDIM), MIN_GRIDDIM, sc.CUDAStream, rotateFieldKernel, Ein1, Ein2, Eout1, Eout2, rotationAngle);
		bilingualMemcpy((*s).EkwOut, Eout1, 2 * (*s).NgridC * sizeof(deviceComplex), cudaMemcpyDeviceToHost);

		//transform back to time
		combinedFFT(&sc, Eout1, sc.gridETime1, 1);
		bilingualLaunch(2 * sc.Nblock, sc.Nthread, sc.CUDAStream, multiplyByConstantKernelD, sc.gridETime1, 1.0 / sc.Ngrid);
		bilingualMemcpy((*s).ExtOut, sc.gridETime1, 2 * (*s).Ngrid * sizeof(double), cudaMemcpyDeviceToHost);

		//update spectrum
		getTotalSpectrum(s, &sc);

		deallocateCudaParameterSet(&sc);
		return 0;
	}
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

	int initializeCudaParameterSet(simulationParameterSet* sCPU, cudaParameterSet* s) {
		//initialize and take values from the struct handed over by the dispatcher
		
		unsigned long long i;
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
		


		//prepare FFT plans
		//explicitly make different plans for GPU or CPU (most other parts of the code can be universal,
		//but not this, since the libraries are different).
#ifdef __CUDACC__
			size_t workSize;
			cufftPlan1d(&(*s).fftPlan1DD2Z, (int)(*s).Ntime, CUFFT_D2Z, 2 * (int)((*s).Nspace * (*s).Nspace2));
			cufftPlan1d(&(*s).fftPlan1DZ2D, (int)(*s).Ntime, CUFFT_Z2D, 2 * (int)((*s).Nspace * (*s).Nspace2));
			cufftSetStream((*s).fftPlan1DD2Z, (*s).CUDAStream);
			cufftSetStream((*s).fftPlan1DZ2D, (*s).CUDAStream);
			if ((*s).is3D) {
				int cufftSizes1[] = { (int)(*s).Nspace2, (int)(*s).Nspace, (int)(*s).Ntime };
				cufftCreate(&(*s).fftPlanD2Z);
				cufftGetSizeMany((*s).fftPlanD2Z, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
				cufftMakePlanMany((*s).fftPlanD2Z, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);

				cufftCreate(&(*s).fftPlanZ2D);
				cufftGetSizeMany((*s).fftPlanZ2D, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);
				cufftMakePlanMany((*s).fftPlanZ2D, 3, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);
			}
			else {
				int cufftSizes1[] = { (int)(*s).Nspace, (int)(*s).Ntime };

				cufftCreate(&(*s).fftPlanD2Z);
				cufftGetSizeMany((*s).fftPlanD2Z, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
				cufftMakePlanMany((*s).fftPlanD2Z, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);

				cufftCreate(&(*s).fftPlanZ2D);
				cufftGetSizeMany((*s).fftPlanZ2D, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);
				cufftMakePlanMany((*s).fftPlanZ2D, 2, cufftSizes1, NULL, 0, 0, 0, 0, 0, CUFFT_Z2D, 2, &workSize);

				if ((*s).isCylindric) {
					int cufftSizes2[] = { 2 * (int)(*s).Nspace, (int)(*s).Ntime };
					cufftCreate(&(*s).doublePolfftPlan);
					cufftGetSizeMany((*s).doublePolfftPlan, 2, cufftSizes2, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
					cufftMakePlanMany((*s).doublePolfftPlan, 2, cufftSizes2, NULL, 0, 0, 0, 0, 0, CUFFT_D2Z, 2, &workSize);
					cufftSetStream((*s).doublePolfftPlan, (*s).CUDAStream);
				}
			}
			cufftSetStream((*s).fftPlanD2Z, (*s).CUDAStream);
			cufftSetStream((*s).fftPlanZ2D, (*s).CUDAStream);
#else
			DftiCreateDescriptor(&(*s).mklPlan1DD2Z, DFTI_DOUBLE, DFTI_REAL, 1, (*s).Ntime);
			DftiSetValue((*s).mklPlan1DD2Z, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
			DftiSetValue((*s).mklPlan1DD2Z, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
			DftiSetValue((*s).mklPlan1DD2Z, DFTI_NUMBER_OF_TRANSFORMS, 2 * (*s).Nspace * (*s).Nspace2);
			DftiSetValue((*s).mklPlan1DD2Z, DFTI_INPUT_DISTANCE, (*s).Ntime);
			DftiSetValue((*s).mklPlan1DD2Z, DFTI_OUTPUT_DISTANCE, (*s).Nfreq);
			DftiCommitDescriptor((*s).mklPlan1DD2Z);

			DftiCreateDescriptor(&(*s).mklPlan1DZ2D, DFTI_DOUBLE, DFTI_REAL, 1, (*s).Ntime);
			DftiSetValue((*s).mklPlan1DZ2D, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
			DftiSetValue((*s).mklPlan1DZ2D, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
			DftiSetValue((*s).mklPlan1DZ2D, DFTI_NUMBER_OF_TRANSFORMS, 2 * (*s).Nspace * (*s).Nspace2);
			DftiSetValue((*s).mklPlan1DZ2D, DFTI_INPUT_DISTANCE, (*s).Nfreq);
			DftiSetValue((*s).mklPlan1DZ2D, DFTI_OUTPUT_DISTANCE, (*s).Ntime);
			DftiCommitDescriptor((*s).mklPlan1DZ2D);

			if ((*s).is3D) {
				MKL_LONG mklSizes[] = { (MKL_LONG)(*s).Nspace2, (MKL_LONG)(*s).Nspace, (MKL_LONG)(*s).Ntime };
				MKL_LONG mklStrides[4] = { 0, (MKL_LONG)((*s).Nspace * (*s).Nfreq), (MKL_LONG)(*s).Nfreq, 1 };
				DftiCreateDescriptor(&(*s).mklPlanD2Z, DFTI_DOUBLE, DFTI_REAL, 3, mklSizes);
				DftiSetValue((*s).mklPlanD2Z, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
				DftiSetValue((*s).mklPlanD2Z, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
				DftiSetValue((*s).mklPlanD2Z, DFTI_OUTPUT_STRIDES, mklStrides);
				DftiSetValue((*s).mklPlanD2Z, DFTI_NUMBER_OF_TRANSFORMS, 2);
				DftiSetValue((*s).mklPlanD2Z, DFTI_INPUT_DISTANCE, (*s).Ngrid);
				DftiSetValue((*s).mklPlanD2Z, DFTI_OUTPUT_DISTANCE, (*s).NgridC);
				DftiCommitDescriptor((*s).mklPlanD2Z);

				DftiCreateDescriptor(&(*s).mklPlanZ2D, DFTI_DOUBLE, DFTI_REAL, 3, mklSizes);
				DftiSetValue((*s).mklPlanZ2D, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
				DftiSetValue((*s).mklPlanZ2D, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
				DftiSetValue((*s).mklPlanZ2D, DFTI_INPUT_STRIDES, mklStrides);
				DftiSetValue((*s).mklPlanZ2D, DFTI_NUMBER_OF_TRANSFORMS, 2);
				DftiSetValue((*s).mklPlanZ2D, DFTI_INPUT_DISTANCE, (*s).NgridC);
				DftiSetValue((*s).mklPlanZ2D, DFTI_OUTPUT_DISTANCE, (*s).Ngrid);
				DftiCommitDescriptor((*s).mklPlanZ2D);
			}
			else {
				MKL_LONG mklSizes[] = { (MKL_LONG)(*s).Nspace, (MKL_LONG)(*s).Ntime };
				MKL_LONG mklStrides[4] = { 0, (MKL_LONG)(*s).Nfreq, 1, 1 };

				DftiCreateDescriptor(&(*s).mklPlanD2Z, DFTI_DOUBLE, DFTI_REAL, 2, mklSizes);
				DftiSetValue((*s).mklPlanD2Z, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
				DftiSetValue((*s).mklPlanD2Z, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
				DftiSetValue((*s).mklPlanD2Z, DFTI_OUTPUT_STRIDES, mklStrides);
				DftiSetValue((*s).mklPlanD2Z, DFTI_NUMBER_OF_TRANSFORMS, 2);
				DftiSetValue((*s).mklPlanD2Z, DFTI_INPUT_DISTANCE, (*s).Ngrid);
				DftiSetValue((*s).mklPlanD2Z, DFTI_OUTPUT_DISTANCE, (*s).NgridC);
				DftiCommitDescriptor((*s).mklPlanD2Z);

				DftiCreateDescriptor(&(*s).mklPlanZ2D, DFTI_DOUBLE, DFTI_REAL, 2, mklSizes);
				DftiSetValue((*s).mklPlanZ2D, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
				DftiSetValue((*s).mklPlanZ2D, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
				DftiSetValue((*s).mklPlanZ2D, DFTI_INPUT_STRIDES, mklStrides);
				DftiSetValue((*s).mklPlanZ2D, DFTI_NUMBER_OF_TRANSFORMS, 2);
				DftiSetValue((*s).mklPlanZ2D, DFTI_INPUT_DISTANCE, (*s).NgridC);
				DftiSetValue((*s).mklPlanZ2D, DFTI_OUTPUT_DISTANCE, (*s).Ngrid);
				DftiCommitDescriptor((*s).mklPlanZ2D);

				if ((*s).isCylindric) {
					MKL_LONG mklSizesD[] = { (MKL_LONG)(2 * (*s).Nspace), (MKL_LONG)(*s).Ntime };
					DftiCreateDescriptor(&(*s).mklPlanDoublePolfft, DFTI_DOUBLE, DFTI_REAL, 2, mklSizesD);
					DftiSetValue((*s).mklPlanDoublePolfft, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
					DftiSetValue((*s).mklPlanDoublePolfft, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
					DftiSetValue((*s).mklPlanDoublePolfft, DFTI_OUTPUT_STRIDES, mklStrides);
					DftiSetValue((*s).mklPlanDoublePolfft, DFTI_NUMBER_OF_TRANSFORMS, 2);
					DftiSetValue((*s).mklPlanDoublePolfft, DFTI_INPUT_DISTANCE, 2 * (*s).Ngrid);
					DftiSetValue((*s).mklPlanDoublePolfft, DFTI_OUTPUT_DISTANCE, 2 * (*s).NgridC);
					DftiCommitDescriptor((*s).mklPlanDoublePolfft);
				}

			}
#endif

		//check if the GPU has enough free memory left for the grids, if not, free the fft plans and quit
#ifdef __CUDACC__
			nvmlDevice_t nvmlDevice = 0;
			nvmlMemory_t nvmlMemoryInfo;
			nvmlInit_v2();
			nvmlDeviceGetHandleByIndex_v2((*sCPU).assignedGPU, &nvmlDevice);
			nvmlDeviceGetMemoryInfo(nvmlDevice, &nvmlMemoryInfo);
			size_t memoryEstimate = sizeof(double) * ((*s).Ngrid * 2 * 2 + 2 * (*s).NgridC * 6 * 2 + 2 * (*s).isCylindric * 5 * 2 + 2 * (*s).Ntime + 2 * (*s).Nfreq + 81 + 65536);
			nvmlShutdown();
			if (nvmlMemoryInfo.free < memoryEstimate) {
				cufftDestroy((*s).fftPlanD2Z);
				cufftDestroy((*s).fftPlanZ2D);
				cufftDestroy((*s).fftPlan1DD2Z);
				cufftDestroy((*s).fftPlan1DZ2D);
				if ((*s).isCylindric) {
					cufftDestroy((*s).doublePolfftPlan);
				}
				(*sCPU).memoryError = -1;
				return 1;
			}	
#endif		


		int memErrors = 0;
		double* expGammaTCPU = (double*)malloc(2 * sizeof(double) * (*s).Ntime);
		if (expGammaTCPU == NULL) return 2;

		size_t beamExpansionFactor = 1;
		if ((*s).isCylindric) {
			beamExpansionFactor = 2;
		}
#ifdef __CUDACC__
		cudaStreamCreate(&(*s).CUDAStream);
#endif
		fillRotationMatricies(sCPU, s);

		//GPU allocations
		//
		// currently 8 large grids, meaning memory use is approximately
		// 64 bytes per grid point (8 grids x 2 polarizations x 4ouble precision)
		// plus a little bit for additional constants/workspaces/etc

		memErrors += bilingualCalloc((void**)&(*s).gridETime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += bilingualCalloc((void**)&(*s).gridPolarizationTime1, 2 * (*s).Ngrid, sizeof(double));
		memErrors += bilingualCalloc((void**)&(*s).workspace1, beamExpansionFactor * 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += bilingualCalloc((void**)&(*s).gridEFrequency1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += bilingualCalloc((void**)&(*s).gridPropagationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += bilingualCalloc((void**)&(*s).gridPolarizationFactor1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += bilingualCalloc((void**)&(*s).gridEFrequency1Next1, 2 * (*s).NgridC, sizeof(std::complex<double>));
		memErrors += bilingualCalloc((void**)&(*s).k1, 2 * (*s).NgridC, sizeof(std::complex<double>));

		//cylindric sym grids
		if ((*s).isCylindric) {
			memErrors += bilingualCalloc((void**)&(*s).gridPropagationFactor1Rho1, 4 * (*s).NgridC, sizeof(std::complex<double>));
			memErrors += bilingualCalloc((void**)&(*s).gridRadialLaplacian1, 4 * (*s).Ngrid, sizeof(std::complex<double>));
		}

		//smaller helper grids
		memErrors += bilingualCalloc((void**)&(*s).expGammaT, 2 * (*s).Ntime, sizeof(double));

		memErrors += bilingualCalloc((void**)&(*s).chiLinear1, 2 * (*s).Nfreq, sizeof(std::complex<double>));
		memErrors += bilingualCalloc((void**)&(*s).inverseChiLinear1, 2 * (*s).Nfreq, sizeof(double));
		for (i = 0; i < (*s).Ntime; i++) {
			expGammaTCPU[i] = exp((*s).dt * i * (*sCPU).drudeGamma);
			expGammaTCPU[i + (*s).Ntime] = exp(-(*s).dt * i * (*sCPU).drudeGamma);
		}
		bilingualMemcpy((*s).expGammaT, expGammaTCPU, 2 * sizeof(double) * (*s).Ntime, cudaMemcpyHostToDevice);
		free(expGammaTCPU);

		memErrors += bilingualCalloc((void**)&(*s).chi3Tensor, 81, sizeof(double));

		(*sCPU).memoryError = memErrors;
		if (memErrors > 0) {
			return memErrors;
		}

		//second polarization grids are to pointers within the first polarization
		//to have contiguous memory
		(*s).gridETime2 = (*s).gridETime1 + (*s).Ngrid;
		(*s).workspace2 = (*s).workspace1 + (*s).NgridC;
		(*s).gridPolarizationTime2 = (*s).gridPolarizationTime1 + (*s).Ngrid;
		(*s).workspace2P = (*s).workspace1 + beamExpansionFactor * (*s).NgridC;
		(*s).k2 = (*s).k1 + (*s).NgridC;
		(*s).chiLinear2 = (*s).chiLinear1 + (*s).Nfreq;
		(*s).inverseChiLinear2 = (*s).inverseChiLinear1 + (*s).Nfreq;
		(*s).gridRadialLaplacian2 = (*s).gridRadialLaplacian1 + (*s).Ngrid;
		(*s).gridPropagationFactor1Rho2 = (*s).gridPropagationFactor1Rho1 + (*s).NgridC;
		(*s).gridPolarizationFactor2 = (*s).gridPolarizationFactor1 + (*s).NgridC;
		(*s).gridEFrequency1Next2 = (*s).gridEFrequency1Next1 + (*s).NgridC;
		(*s).gridPropagationFactor2 = (*s).gridPropagationFactor1 + (*s).NgridC;
		(*s).gridEFrequency2 = (*s).gridEFrequency1 + (*s).NgridC;


		//prepare effective nonlinearity tensors and put them on the GPU

		double firstDerivativeOperation[6] = { -1. / 60.,  3. / 20., -3. / 4.,  3. / 4.,  -3. / 20., 1. / 60. };
		for (i = 0; i < 6; i++) {
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

		calcEffectiveChi2Tensor((*sCPU).deffTensor, (*sCPU).chi2Tensor, (*sCPU).crystalTheta, (*sCPU).crystalPhi);
		memcpy((*s).chi2Tensor, (*sCPU).deffTensor, 9 * sizeof(double));
		memcpy((*s).nonlinearSwitches, (*sCPU).nonlinearSwitches, 4 * sizeof(int));

		bilingualMemcpy((*s).chi3Tensor, (*sCPU).chi3Tensor, 81 * sizeof(double), cudaMemcpyHostToDevice);
		memcpy((*s).absorptionParameters, (*sCPU).absorptionParameters, 6 * sizeof(double));
		memcpy((*s).plasmaParameters, plasmaParametersCPU, 6 * sizeof(double));
		memcpy((*s).firstDerivativeOperation, firstDerivativeOperation, 6 * sizeof(double));


		return 0;
	}

	int deallocateCudaParameterSet(cudaParameterSet* s) {
		bilingualFree((*s).gridETime1);
		bilingualFree((*s).workspace1);
		bilingualFree((*s).gridEFrequency1);
		bilingualFree((*s).gridPropagationFactor1);
		if ((*s).isCylindric) {
			bilingualFree((*s).gridPropagationFactor1Rho1);
			bilingualFree((*s).gridRadialLaplacian1);
		}
		bilingualFree((*s).gridPolarizationFactor1);
		bilingualFree((*s).gridEFrequency1Next1);
		bilingualFree((*s).k1);
		bilingualFree((*s).gridPolarizationTime1);
		bilingualFree((*s).chi3Tensor);
		bilingualFree((*s).expGammaT);
		bilingualFree((*s).chiLinear1);
		bilingualFree((*s).inverseChiLinear1);
#ifdef __CUDACC__
		cufftDestroy((*s).fftPlanD2Z);
		cufftDestroy((*s).fftPlanZ2D);
		cufftDestroy((*s).fftPlan1DD2Z);
		cufftDestroy((*s).fftPlan1DZ2D);
		if ((*s).isCylindric) {
			cufftDestroy((*s).doublePolfftPlan);
		}
		cudaStreamDestroy((*s).CUDAStream);
#else
		DftiFreeDescriptor(&(*s).mklPlan1DD2Z);
		DftiFreeDescriptor(&(*s).mklPlanD2Z);
		DftiFreeDescriptor(&(*s).mklPlanZ2D);
		if ((*s).isCylindric)DftiFreeDescriptor(&(*s).mklPlanDoublePolfft);
#endif
		return 0;
	}

//function to run a RK4 time step
//stepNumber is the sub-step index, from 0 to 3
	int runRK4Step(cudaParameterSet* sH, cudaParameterSet* sD, uint8_t stepNumber) {
		//operations involving FFT
		if ((*sH).isNonLinear || (*sH).isCylindric) {
			//perform inverse FFT to get time-space electric field
			combinedFFT(sH, (*sH).workspace1, (*sH).gridETime1, 1);

			//Nonlinear polarization
			if ((*sH).isNonLinear) {
				bilingualLaunch((*sH).Nblock, (*sH).Nthread, (*sH).CUDAStream, nonlinearPolarizationKernel, sD);
				if ((*sH).isCylindric) {
					bilingualLaunch((*sH).Nblock, (*sH).Nthread, (*sH).CUDAStream, expandCylindricalBeam, sD, (*sH).gridPolarizationTime1, (*sH).gridPolarizationTime2);
					combinedFFT(sH, (*sH).gridRadialLaplacian1, (*sH).workspace1, 4);
				}
				else {
					combinedFFT(sH, (*sH).gridPolarizationTime1, (*sH).workspace1, 0);
				}
				bilingualLaunch((*sH).Nblock / 2, (*sH).Nthread, (*sH).CUDAStream, updateKwithPolarizationKernel, sD);
			}

			//Plasma/multiphoton absorption
			if ((*sH).hasPlasma) {
				//bilingualLaunch((unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / MIN_GRIDDIM), MIN_GRIDDIM, (*sH).CUDAStream, plasmaCurrentKernel, sD);

				bilingualLaunch((*sH).Nblock, (*sH).Nthread, (*sH).CUDAStream, plasmaCurrentKernel_twoStage_A, sD);
				bilingualLaunch((unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / MIN_GRIDDIM), 2*MIN_GRIDDIM, (*sH).CUDAStream, plasmaCurrentKernel_twoStage_B, sD);
				if ((*sH).isCylindric) {
					bilingualLaunch((*sH).Nblock, (*sH).Nthread, (*sH).CUDAStream, expandCylindricalBeam, sD, (*sH).gridPolarizationTime1, (*sH).gridPolarizationTime2);
					combinedFFT(sH, (*sH).gridRadialLaplacian1, (*sH).workspace1, 4);
				}
				else {
					combinedFFT(sH, (*sH).gridPolarizationTime1, (*sH).workspace1, 0);
				}
				bilingualLaunch((*sH).Nblock / 2, (*sH).Nthread, (*sH).CUDAStream, updateKwithPlasmaKernel, sD);
			}

			//Radial Laplacian
			if ((*sH).isCylindric) {
				bilingualLaunch((*sH).Nblock, (*sH).Nthread, (*sH).CUDAStream, radialLaplacianKernel, sD);
				combinedFFT(sH, (*sH).gridRadialLaplacian1, (*sH).workspace1, 0);
			}
		}

		//advance an RK4 step
		bilingualLaunch((*sH).Nblock / 2, (*sH).Nthread, (*sH).CUDAStream, rkKernel, sD, stepNumber);
		return 0;
	}


#ifdef __CUDACC__
int resolveSequence(int currentIndex, simulationParameterSet* s, crystalEntry* db) {
#else
int resolveSequenceCPU(int currentIndex, simulationParameterSet * s, crystalEntry * db) {
#endif

	double* offsetArray = &(*s).sequenceArray[11 * currentIndex];
	int error = 0;
	//sequence format
	//0: step type
	int stepType = (int)offsetArray[0];
	int materialIndex = 0;
	double thickness = 0;
	// 
	// if stepType == 0, normal propagation
	//1: material index
	//2: theta,
	//3: phi, 
	//4: NL absorption
	//5: Band gap
	//6: Drude relaxation
	//7: Effective mass
	//8: Crystal thickness
	//9: Propagation step size
	//10: rotation angle
	//
	// if stepType == 1, linear propagation
	// same parameters as 0, but only 1,2,3,8, and 10 matter
	//
	// if stepType == 2, fresnel loss
	// 1: incidence material index
	// 2: transmission material index
	// other parameters don't matter
	// 
	// if stepType == 3, spherical mirror
	// 1: ROC (m)
	//
	// if stepType == 4, parabolic mirror
	// 1: focus (m)
	// 
	// if stepType == 5, aperture
	// 1: diameter (m)
	// 2: activation parameter p (function is 1 - 1/(1 + exp(-p*(r-radius)))
	//
	// if stepType == 6, loop back to start (handled by solveNonlinearWaveEquationSequence())
	// 1: counter (counts down to zero)
	//
	// if stepType == 7, reinjection, same as 0, but input fields are added to current fields.

	switch (stepType) {
	case 7:
		(*s).isReinjecting = TRUE;
	case 0:
		if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
		if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
		if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
		if ((int)offsetArray[4] != -1) (*s).nonlinearAbsorptionStrength = offsetArray[4];
		if ((int)offsetArray[5] != -1) (*s).bandGapElectronVolts = offsetArray[5];
		if ((int)offsetArray[6] != -1) (*s).drudeGamma = offsetArray[6];
		if ((int)offsetArray[7] != -1) (*s).effectiveMass = offsetArray[7];
		if ((int)offsetArray[8] != -1) (*s).crystalThickness = 1e-6 * offsetArray[8];
		if ((int)offsetArray[9] != -1) (*s).propagationStep = 1e-9 * offsetArray[9];
		if ((int)offsetArray[8] != -1) (*s).Npropagation
			= (size_t)(1e-6 * offsetArray[8] / (*s).propagationStep);
		if (currentIndex > 0) {
			(*s).isFollowerInSequence = TRUE;
		}
		(*s).chi2Tensor = db[(*s).materialIndex].d;
		(*s).chi3Tensor = db[(*s).materialIndex].chi3;
		(*s).nonlinearSwitches = db[(*s).materialIndex].nonlinearSwitches;
		(*s).absorptionParameters = db[(*s).materialIndex].absorptionParameters;
		(*s).sellmeierCoefficients = db[(*s).materialIndex].sellmeierCoefficients;

		(*s).sellmeierType = db[(*s).materialIndex].sellmeierType;
		(*s).axesNumber = db[(*s).materialIndex].axisType;

#ifdef __CUDACC__
		error = solveNonlinearWaveEquation(s);
#else
		error = solveNonlinearWaveEquationCPU(s);
#endif
		if (offsetArray[10] != 0.0) {
			rotateField(s, DEG2RAD * offsetArray[10]);
		}

		if ((*s).memoryError > 0) {
			printf("Warning: device memory error (%i).\n", (*s).memoryError);
		}
		return error;

	case 1:
		if ((*s).isCylindric) {
			if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
			if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
			if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
			if ((int)offsetArray[4] != -1) (*s).nonlinearAbsorptionStrength = offsetArray[4];
			if ((int)offsetArray[5] != -1) (*s).bandGapElectronVolts = offsetArray[5];
			if ((int)offsetArray[6] != -1) (*s).drudeGamma = offsetArray[6];
			if ((int)offsetArray[7] != -1) (*s).effectiveMass = offsetArray[7];
			if ((int)offsetArray[8] != -1) (*s).crystalThickness = 1e-6 * offsetArray[8];
			if ((int)offsetArray[9] != -1) (*s).propagationStep = 1e-9 * offsetArray[9];
			if ((int)offsetArray[8] != -1) (*s).Npropagation
				= (size_t)(1e-6 * offsetArray[8] / (*s).propagationStep);
			if (currentIndex > 0) {
				(*s).isFollowerInSequence = TRUE;
			}
			(*s).chi2Tensor = db[(*s).materialIndex].d;
			(*s).chi3Tensor = db[(*s).materialIndex].chi3;
			(*s).nonlinearSwitches = db[(*s).materialIndex].nonlinearSwitches;
			(*s).absorptionParameters = db[(*s).materialIndex].absorptionParameters;
			(*s).sellmeierCoefficients = db[(*s).materialIndex].sellmeierCoefficients;
			(*s).sellmeierType = db[(*s).materialIndex].sellmeierType;
			(*s).axesNumber = db[(*s).materialIndex].axisType;
			(*s).forceLinear = TRUE;
#ifdef __CUDACC__
			error = solveNonlinearWaveEquation(s);
#else
			error = solveNonlinearWaveEquationCPU(s);
#endif
		}
		else {
			if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
			if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
			if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
			thickness = 1.0e-6 * offsetArray[8];
			if (offsetArray[8] == -1) {
				thickness = (*s).crystalThickness;
			}
			materialIndex = (int)offsetArray[1];
			if (offsetArray[1] == -1) {
				materialIndex = (*s).materialIndex;
			}
			applyLinearPropagation(s, materialIndex, thickness);
		}

		if (offsetArray[10] != 0.0) {
			rotateField(s, DEG2RAD * offsetArray[10]);
		}
		return 0;

	case 2:
		if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
		if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
		if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
		applyFresnelLoss(s, (int)offsetArray[4], (int)offsetArray[5]);
		return 0;
	case 3:
		applySphericalMirror(s, offsetArray[8]);
		if (offsetArray[10] != 0.0) {
			rotateField(s, DEG2RAD * offsetArray[10]);
		}
		return 0;
	case 4:
		applyParabolicMirror(s, offsetArray[8]);
		if (offsetArray[10] != 0.0) {
			rotateField(s, DEG2RAD * offsetArray[10]);
		}
		return 0;
	case 5:
		applyAperature(s, offsetArray[1], offsetArray[2]);
		if (offsetArray[10] != 0.0) {
			rotateField(s, DEG2RAD * offsetArray[10]);
		}
		return 0;
	}

	return 1;
	}


	void runFittingIteration(int* m, int* n, double* fittingValues, double* fittingFunction) {
		int i;
		int fitLocation;
		double referenceValue;
		//pointers to values that can be scanned in batch mode
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

		double* references[36] = { 0,
		&(*fittingReferenceSet).pulseEnergy1, &(*fittingReferenceSet).pulseEnergy2, &(*fittingReferenceSet).frequency1, &(*fittingReferenceSet).frequency2,
		&(*fittingReferenceSet).bandwidth1, &(*fittingReferenceSet).bandwidth2, &(*fittingReferenceSet).cephase1, &(*fittingReferenceSet).cephase2,
		&(*fittingReferenceSet).delay1, &(*fittingReferenceSet).delay2, &(*fittingReferenceSet).gdd1, &(*fittingReferenceSet).gdd2,
		&(*fittingReferenceSet).tod1, &(*fittingReferenceSet).tod2, &(*fittingReferenceSet).phaseMaterialThickness1, &(*fittingReferenceSet).phaseMaterialThickness2,
		&(*fittingReferenceSet).beamwaist1, &(*fittingReferenceSet).beamwaist2,
		&(*fittingReferenceSet).x01, &(*fittingReferenceSet).x02, &(*fittingReferenceSet).z01, &(*fittingReferenceSet).z02,
		&(*fittingReferenceSet).propagationAngle1, &(*fittingReferenceSet).propagationAngle2, &(*fittingReferenceSet).polarizationAngle1, &(*fittingReferenceSet).polarizationAngle2,
		&(*fittingReferenceSet).circularity1, &(*fittingReferenceSet).circularity2, &(*fittingReferenceSet).crystalTheta, &(*fittingReferenceSet).crystalPhi,
		&(*fittingReferenceSet).nonlinearAbsorptionStrength, &(*fittingReferenceSet).drudeGamma, &(*fittingReferenceSet).effectiveMass, &(*fittingReferenceSet).crystalThickness,
		&(*fittingReferenceSet).propagationStep };


		for (i = 0; i < *n; i++) {
			fitLocation = (int)round((*fittingSet).fittingArray[3 * i]);
			if (fitLocation > 35) return;
			referenceValue = *references[fitLocation];
			if (referenceValue == 0.0) {
				referenceValue = 1.;
			}
			*targets[fitLocation] = fittingValues[i] * referenceValue;
		}
#ifdef __CUDACC__
		if ((*fittingSet).isInSequence) {
			solveNonlinearWaveEquationSequence(fittingSet);
			(*fittingSet).isFollowerInSequence = FALSE;
		}
		else {
			solveNonlinearWaveEquation(fittingSet);
		}
#else
		if ((*fittingSet).isInSequence) {
			solveNonlinearWaveEquationSequenceCPU(fittingSet);
			(*fittingSet).isFollowerInSequence = FALSE;
		}
		else {
			solveNonlinearWaveEquationCPU(fittingSet);
		}
#endif

		//mode 0: maximize total spectrum in ROI
		if ((*fittingSet).fittingMode == 0) {
			for (i = 0; i < *m; i++) {
				fittingFunction[i] = (1.0e8 / ((*fittingSet).totalSpectrum[2 * (*fittingSet).Nfreq + (*fittingSet).fittingROIstart + i]));
			}
		}
		//mode 1: maximize s-polarized spectrum in ROI
		if ((*fittingSet).fittingMode == 1) {
			for (i = 0; i < *m; i++) {
				fittingFunction[i] = (1.0e8 / ((*fittingSet).totalSpectrum[(*fittingSet).fittingROIstart + i]));
			}
		}
		//mode 2: maximize p-polarized spectrum in ROI
		if ((*fittingSet).fittingMode == 2) {
			for (i = 0; i < *m; i++) {
				fittingFunction[i] = (1.0e8 / ((*fittingSet).totalSpectrum[(*fittingSet).Nfreq + (*fittingSet).fittingROIstart + i]));
			}
		}
		//mode 3: match total spectrum to reference given in ascii file
		if ((*fittingSet).fittingMode == 3) {
			double maxSim = 0;
			double maxRef = 0;
			double sumSim = 0;
			double sumRef = 0;
			double* simSpec = &(*fittingSet).totalSpectrum[2 * (*fittingSet).Nfreq + (*fittingSet).fittingROIstart];
			double* refSpec = &(*fittingSet).fittingReference[(*fittingSet).fittingROIstart];
			for (i = 0; i < *m; i++) {
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

			double sumFF = 0;
			for (i = 0; i < *m; i++) {
				fittingFunction[i] = log10(1e5 * refSpec[i] / maxRef) - log10(1e5 * simSpec[i] / maxSim);
				sumFF += fittingFunction[i];
				//fittingFunction[i] = 1.0e8 / ((*fittingSet).totalSpectrum[(*fittingSet).Ntime + (*fittingSet).fittingROIstart + i]);
			}
			sumFF /= *m;
			for (i = 0; i < *m; i++) {
				fittingFunction[i] -= sumFF;
			}
		}


		return;
	}

	int getTotalSpectrum(simulationParameterSet* sCPU, cudaParameterSet* sc) {

		bilingualMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		combinedFFT(sc, (*sc).gridETime1, (*sc).workspace1, 2);
		if ((*sc).is3D) {
			bilingualLaunch((unsigned int)(*sCPU).Nfreq, 1u, (*sc).CUDAStream, totalSpectrum3DKernel, (*sc).workspace1, (*sc).workspace2, (*sCPU).rStep, (*sCPU).Ntime / 2 + 1, (*sCPU).Nspace * (*sCPU).Nspace2, (*sc).gridPolarizationTime1);
		}
		else {
			bilingualLaunch((unsigned int)(*sCPU).Nfreq, 1u, (*sc).CUDAStream, totalSpectrumKernel, (*sc).workspace1, (*sc).workspace2, (*sCPU).rStep, (*sCPU).Ntime / 2 + 1, (*sCPU).Nspace, (*sc).gridPolarizationTime1);
		}

		bilingualMemcpy((*sCPU).totalSpectrum, (*sc).gridPolarizationTime1, 3 * (*sCPU).Nfreq * sizeof(double), cudaMemcpyDeviceToHost);


		return 0;
	}


}
//END OF NAMESPACE

#ifdef __CUDACC__
unsigned long solveNonlinearWaveEquation(void* lpParam) {
	simulationParameterSet* sCPU = (simulationParameterSet*)lpParam;
	cudaSetDevice((*sCPU).assignedGPU);
#else
unsigned long solveNonlinearWaveEquationCPU(void* lpParam) {
	simulationParameterSet* sCPU = (simulationParameterSet*)lpParam;
#endif
	size_t i;
	cudaParameterSet* sDevice;
	cudaParameterSet s;
	memset(&s, 0, sizeof(cudaParameterSet));
	if (initializeCudaParameterSet(sCPU, &s)) return 1;

	//prepare the propagation arrays
	if (s.is3D) {
		preparePropagation3D(sCPU, s);
	}
	else if (s.isCylindric) {
		preparePropagation3DCylindric(sCPU, s);
	}
	else {
		preparePropagation2DCartesian(sCPU, s);
	}
	prepareElectricFieldArrays(sCPU, &s);
	double canaryPixel = 0;
	double* canaryPointer = &s.gridETime1[s.Ntime / 2 + s.Ntime * (s.Nspace / 2 + s.Nspace * (s.Nspace2 / 2))];

	bilingualCalloc((void**)&sDevice, 1, sizeof(cudaParameterSet));
	bilingualMemcpy(sDevice, &s, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);

	//Core propagation loop
	for (i = 0; i < s.Nsteps; i++) {

		//RK4
		runRK4Step(&s, sDevice, 0);
		runRK4Step(&s, sDevice, 1);
		runRK4Step(&s, sDevice, 2);
		runRK4Step(&s, sDevice, 3);
#ifdef __CUDACC__
		cudaMemcpyAsync(&canaryPixel, canaryPointer, sizeof(double), cudaMemcpyDeviceToHost);
#else
		canaryPixel = *canaryPointer;
#endif
		if (isnan(canaryPixel)) {
			break;
		}

		if ((*sCPU).imdone[0] == 2) {
			break;
		}

		if ((*sCPU).imdone[0] == 3) {
			//copy the field arrays from the GPU to CPU memory if requested by the UI
			bilingualMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), cudaMemcpyDeviceToHost);
			bilingualMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * (*sCPU).Ngrid * sizeof(deviceComplex), cudaMemcpyDeviceToHost);

			(*sCPU).imdone[0] = 0;
		}
		(*(*sCPU).progressCounter)++;
	}

	////give the result to the CPU
	bilingualMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), cudaMemcpyDeviceToHost);


	combinedFFT(&s, s.gridEFrequency1, s.gridETime1, 1);

	bilingualLaunch((int)(s.Ngrid / MIN_GRIDDIM), 2 * MIN_GRIDDIM, s.CUDAStream, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);
	bilingualMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), cudaMemcpyDeviceToHost);

	getTotalSpectrum(sCPU, &s);

	deallocateCudaParameterSet(&s);
	bilingualFree(sDevice);
	(*sCPU).imdone[0] = 1;
	return isnan(canaryPixel) * 13;
}

#ifdef __CUDACC__
unsigned long solveNonlinearWaveEquationSequence(void* lpParam) {
#else
unsigned long solveNonlinearWaveEquationSequenceCPU(void* lpParam) {
#endif
	simulationParameterSet* sCPU = (simulationParameterSet*)lpParam;
	simulationParameterSet sCPUbackup[1];
	memcpy(sCPUbackup, sCPU, sizeof(simulationParameterSet));
	int k;
	int error = 0;
	for (k = 0; k < (*sCPU).Nsequence; k++) {
		if ((int)round((*sCPU).sequenceArray[k * 11]) == 6
			&& ((int)round((*sCPU).sequenceArray[k * 11 + 1])) > 0) {
			(*sCPUbackup).sequenceArray[k * 11 + 1] -= 1.0;
			(*sCPUbackup).isFollowerInSequence = TRUE;
			k = 0;
		}
#ifdef __CUDACC__
		error = resolveSequence(k, sCPU, (*sCPU).crystalDatabase);
#else
		error = resolveSequenceCPU(k, sCPU, (*sCPU).crystalDatabase);
#endif

		if (error) break;
		memcpy(sCPU, sCPUbackup, sizeof(simulationParameterSet));
	}

	return error;
}


#ifdef __CUDACC__
unsigned long runFitting(simulationParameterSet* sCPU) {
#else
unsigned long runFittingCPU(simulationParameterSet * sCPU) {
#endif
	int n = (int)(*sCPU).Nfitting;
	int m = (int)(*sCPU).fittingROIsize;
	fittingReferenceSet = sCPU;
	fittingSet = (simulationParameterSet*)malloc((*sCPU).Nsims * sizeof(simulationParameterSet));
	memcpy(fittingSet, sCPU, (*sCPU).Nsims * sizeof(simulationParameterSet));

	double commonPrecision = (*sCPU).fittingPrecision;
	const double eps[6] = { commonPrecision,commonPrecision,commonPrecision,commonPrecision,commonPrecision,commonPrecision }; /* set precisions for stop-criteria */
	double jacobianPrecision = commonPrecision;
	double* x = (double*)mkl_malloc(sizeof(double) * n, 64);
	double* fittingValues = (double*)mkl_malloc(sizeof(double) * m, 64);
	double* fjac = (double*)mkl_malloc(sizeof(double) * m * n, 64);
	double* lowerBounds = (double*)mkl_malloc(sizeof(double) * n, 64);
	double* upperBounds = (double*)mkl_malloc(sizeof(double) * n, 64);
	const int maxIterations = max((*sCPU).fittingMaxIterations, 2);
	const int maxTrialIterations = max(maxIterations / 10, 2);
	/* initial step bound */
	double rs = 0.0;
	int RCI_Request;
	int successful;

	int iter;
	int stopCriterion;
	double inputResiduals = 0.0, outputResiduals = 0.0;
	_TRNSPBC_HANDLE_t handle;
	int i;
	int error = 0;

	//initial guess and bounds
	for (i = 0; i < n; i++) {
		x[i] = 1.;
		upperBounds[i] = (*fittingSet).fittingArray[3 * i + 2];
		lowerBounds[i] = (*fittingSet).fittingArray[3 * i + 1];
	}

	//initialize fitting function and jacobian
	for (i = 0; i < m; i++) {
		fittingValues[i] = 0.0;
	}
	for (i = 0; i < m * n; i++) {
		fjac[i] = 0.0;
	}

	error += dtrnlspbc_init(&handle, &n, &m, x, lowerBounds, upperBounds, eps, &maxIterations, &maxTrialIterations, &rs) != TR_SUCCESS;
	size_t currentIteration = 0;
	if (error == 0) {
		RCI_Request = 0;
		successful = 0;

		while (successful == 0 && (*sCPU).imdone[0] != 2 && currentIteration < maxIterations)
		{
			currentIteration++;
			if (dtrnlspbc_solve(&handle, fittingValues, fjac, &RCI_Request) != TR_SUCCESS)
			{
				successful = -1;
			}

			//check convergence
			if (RCI_Request > -7 && RCI_Request < -1) successful = 1;

			//recalculate
			if (RCI_Request == 1)
			{
				runFittingIteration(&m, &n, x, fittingValues);
			}

			//make jacobian
			if (RCI_Request == 2)
			{
				djacobi(runFittingIteration, &n, &m, fjac, x, &jacobianPrecision);
			}
		}
	}


	dtrnlspbc_get(&handle, &iter, &stopCriterion, &inputResiduals, &outputResiduals);
	memcpy(sCPU, fittingSet, (*fittingSet).Nsims * sizeof(simulationParameterSet));
#ifdef __CUDACC__
	solveNonlinearWaveEquation(sCPU);
#endif
	//free memory
	dtrnlspbc_delete(&handle);
	mkl_free(upperBounds);
	mkl_free(lowerBounds);
	mkl_free(fjac);
	mkl_free(fittingValues);
	mkl_free(x);
	MKL_Free_Buffers();
	free(fittingSet);
	return 0;
}



#ifdef __CUDACC__
int main(int argc, char* argv[]) {
	char* filepath = argv[1];
#else
int mainCPU(int argc, char* filepath) {
#endif
	int i, j;

	size_t progressCounter = 0;
	int CUDAdeviceCount = 1;
#ifdef __CUDACC__
	int CUDAdevice;
	cudaGetDeviceCount(&CUDAdeviceCount);
	cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
	struct cudaDeviceProp activeCUDADeviceProp;
	if (cuErr == cudaSuccess) {
		printf("Found %i GPU(s): \n", CUDAdeviceCount);
		for (i = 0; i < CUDAdeviceCount; i++) {
			cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
			printf("%s\r\n", activeCUDADeviceProp.name);
			printf(" Memory: %lli MB; Multiprocessors: %i\n",
				activeCUDADeviceProp.totalGlobalMem / (1024 * 1024), activeCUDADeviceProp.multiProcessorCount);
		}
	}
	else {
		printf("No GPU found.\n");
		return 1;
	}
#endif
	if (argc < 2) {
		printf("no input file specified.\n");
		return 2;
	}

	// allocate databases, main structs
	simulationParameterSet* sCPU = (simulationParameterSet*)calloc(512, sizeof(simulationParameterSet));
	crystalEntry* crystalDatabasePtr = (crystalEntry*)calloc(512, sizeof(crystalEntry));
	(*sCPU).crystalDatabase = crystalDatabasePtr;
	(*sCPU).progressCounter = &progressCounter;
	// read crystal database
	if (readCrystalDatabase(crystalDatabasePtr) == -2) {
		return 11;
	}
	if ((*crystalDatabasePtr).numberOfEntries == 0) {
		printf("Could not read crystal database.\n");
		free(sCPU);
		free(crystalDatabasePtr);
		return 12;
	}
	printf("Read %i crystal database entries:\n", (*crystalDatabasePtr).numberOfEntries);
	for (j = 0; j < (*crystalDatabasePtr).numberOfEntries; j++) {
		printf("Material %i name: %ls", j, crystalDatabasePtr[j].crystalNameW);
	}

	// read from settings file
	if (readInputParametersFile(sCPU, crystalDatabasePtr, filepath) == 1) {
		printf("Could not read input file.\n");
		free(sCPU);
		free(crystalDatabasePtr);
		return 13;
	}

	allocateGrids(sCPU);
	if (loadPulseFiles(sCPU) == 1) {
		printf("Could not read pulse file.\n");
		free((*sCPU).imdone);
		free((*sCPU).deffTensor);
		free((*sCPU).loadedField1);
		free((*sCPU).loadedField2);
		free(sCPU);
		free(crystalDatabasePtr);
		return 14;
	}

	readSequenceString(sCPU);
	printf("Found %i steps in sequence\n", (*sCPU).Nsequence);
	readFittingString(sCPU);
	configureBatchMode(sCPU);

	auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

	// run simulations
	if ((*sCPU).isInFittingMode) {
		if ((*sCPU).fittingMode == 3) {
			if (loadReferenceSpectrum((*sCPU).fittingPath, sCPU)) {
				printf("Could not load reference spectrum!\n");
				free((*sCPU).imdone);
				free((*sCPU).deffTensor);
				free((*sCPU).loadedField1);
				free((*sCPU).loadedField2);
				free((*sCPU).ExtOut);
				free((*sCPU).EkwOut);
				free((*sCPU).totalSpectrum);
				free((*sCPU).fittingReference);
				free(sCPU);
				free(crystalDatabasePtr);
				return 10;
			}
		}
		printf("Running in fitting mode -- I don't know how long this will take!\n");
#ifdef __CUDACC__
		runFitting(sCPU);
#else
		runFittingCPU(sCPU);
#endif
		auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
		printf("Finished after %8.4lf s. \n",
			1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));

		saveDataSet(sCPU, crystalDatabasePtr, (*sCPU).outputBasePath, FALSE);

		free((*sCPU).imdone);
		free((*sCPU).deffTensor);
		free((*sCPU).loadedField1);
		free((*sCPU).loadedField2);
		free((*sCPU).ExtOut);
		free((*sCPU).EkwOut);
		free((*sCPU).totalSpectrum);
		free((*sCPU).fittingReference);
		free(sCPU);
		free(crystalDatabasePtr);

		return 0;
	}
	std::thread* threadBlock = (std::thread*)calloc((*sCPU).Nsims * (*sCPU).Nsims2, sizeof(std::thread));
	size_t maxThreads = min(CUDAdeviceCount, (*sCPU).Nsims * (*sCPU).Nsims2);
	for (j = 0; j < (*sCPU).Nsims * (*sCPU).Nsims2; j++) {

		sCPU[j].assignedGPU = j % CUDAdeviceCount;
		if (j >= maxThreads) {
			if (threadBlock[j - maxThreads].joinable()) {
				threadBlock[j - maxThreads].join();
			}
		}
#ifdef __CUDACC__
		if ((*sCPU).isInSequence) {
			threadBlock[j] = std::thread(solveNonlinearWaveEquationSequence, &sCPU[j]);
		}
		else {
			threadBlock[j] = std::thread(solveNonlinearWaveEquation, &sCPU[j]);
		}
#else
		if ((*sCPU).isInSequence) {
			threadBlock[j] = std::thread(solveNonlinearWaveEquationSequenceCPU, &sCPU[j]);
		}
		else {
			threadBlock[j] = std::thread(solveNonlinearWaveEquationCPU, &sCPU[j]);
		}
#endif
	}

	for (i = 0; i < (*sCPU).Nsims * (*sCPU).Nsims2; i++) {
		if (sCPU[i].memoryError > 0) {
			printf("Warning: device memory error (%i).\n", sCPU[i].memoryError);
		}
		if (threadBlock[i].joinable()) {
			threadBlock[i].join();
		}
	}

	auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
	printf("Finished after %8.4lf s. \n",
		1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));


	saveDataSet(sCPU, crystalDatabasePtr, (*sCPU).outputBasePath, FALSE);
	//free
	free(threadBlock);
	free((*sCPU).imdone);
	free((*sCPU).deffTensor);
	free((*sCPU).loadedField1);
	free((*sCPU).loadedField2);
	free((*sCPU).ExtOut);
	free((*sCPU).EkwOut);
	free((*sCPU).totalSpectrum);
	free(sCPU);
	free(crystalDatabasePtr);
	return 0;
}