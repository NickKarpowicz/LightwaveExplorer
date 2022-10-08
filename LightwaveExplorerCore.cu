#include "LightwaveExplorerCore.cuh"
//#include "LightwaveExplorerCoreCPU.h"
//#include "LightwaveExplorerSYCL.h"
#include "LightwaveExplorerUtilities.h"
#include "LightwaveExplorerTrilingual.h"
#include <stdlib.h>
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

//Name assignments for the result functions compiled under CUDA, SYCL, and c++
#ifdef __CUDACC__
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <nvml.h>
#define deviceFunctions deviceFunctionsCUDA
#define hostFunctions hostFunctionsCUDA
#define mainX main
#define mainArgumentX char* argv[]
#define resolveArgv char* filepath = argv[1];
#define runDlibFittingX runDlibFitting
#define solveNonlinearWaveEquationX solveNonlinearWaveEquation
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequence
#elif defined RUNONSYCL
#define deviceFunctions deviceFunctionsSYCL
#define hostFunctions hostFunctionsSYCL
#define mainX mainSYCL
#define mainArgumentX char* filepath
#define resolveArgv
#define runDlibFittingX runDlibFittingSYCL
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationSYCL
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceSYCL
#else
#define deviceFunctions deviceFunctionsCPU
#define hostFunctions hostFunctionsCPU
#define mainX mainCPU
#define mainArgumentX char* filepath
#define resolveArgv
#define runDlibFittingX runDlibFittingCPU
#define solveNonlinearWaveEquationX solveNonlinearWaveEquationCPU
#define solveNonlinearWaveEquationSequenceX solveNonlinearWaveEquationSequenceCPU
#endif

namespace deviceFunctions {
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
	deviceFunction deviceComplex sellmeierFunc(double ls, double omega, double* a, int eqn) {
		double realPart;
		deviceComplex compPart;
		double omega2 = omega * omega;
		switch (eqn) {
		case 0:
			realPart = a[0]
				+ (a[1] + a[2] * ls) / (ls + a[3])
				+ (a[4] + a[5] * ls) / (ls + a[6])
				+ (a[7] + a[8] * ls) / (ls + a[9])
				+ (a[10] + a[11] * ls) / (ls + a[12])
				+ a[13] * ls
				+ a[14] * ls * ls
				+ a[15] * ls * ls * ls;
			compPart = a[16] / deviceComplex(a[17] - omega2, a[18] * omega)
				+ a[19] / deviceComplex(a[20] - omega2, a[21] * omega);
			return deviceLib::sqrt(maxN(realPart, 0.0) + KLORENTZIAN * compPart);
		case 1:
			compPart = a[0] / deviceComplex(a[1] - omega2, a[2] * omega)
				+ a[3] / deviceComplex(a[4] - omega2, a[5] * omega)
				+ a[6] / deviceComplex(a[7] - omega2, a[8] * omega)
				+ a[9] / deviceComplex(a[10] - omega2, a[11] * omega)
				+ a[12] / deviceComplex(a[13] - omega2, a[14] * omega)
				+ a[15] / deviceComplex(a[16] - omega2, a[17] * omega)
				+ a[18] / deviceComplex(a[19] - omega2, a[20] * omega);
			return deviceLib::sqrt(KLORENTZIAN * compPart);
		}
		return deviceComplex(1.0, 0.0);
	};

	//Sellmeier equation for refractive indicies
	deviceFunction deviceComplex sellmeierCuda(
		deviceComplex* ne, deviceComplex* no, double* a, double f, double theta, double phi, int type, int eqn) {
		if (f == 0) return deviceComplex(1.0, 0.0); //exit immediately for f=0
		double ls = 2.99792458e14 / f; //wavelength in microns
		ls *= ls; //only wavelength^2 is ever used
		double omega = TWOPI * maxN(f,-f);

		//option 0: isotropic
		if (type == 0) {
			*ne = sellmeierFunc(ls, omega, a, eqn);
			*no = *ne;
			return *ne;
		}
		//option 1: uniaxial
		else if (type == 1) {
			deviceComplex na = sellmeierFunc(ls, omega, a, eqn);
			deviceComplex nb = sellmeierFunc(ls, omega, &a[22], eqn);
			*no = na;
			*ne = 1.0 / deviceLib::sqrt(cos(theta) * cos(theta) / (na * na) + sin(theta) * sin(theta) / (nb * nb));
			return *ne;
		}
		else {
			//type == 2: biaxial
			// X. Yin, S. Zhang and Z. Tian, Optics and Laser Technology 39 (2007) 510 - 513.
			// I am sorry if there is a bug and you're trying to find it, i did my best.
			deviceComplex na = sellmeierFunc(ls, omega, a, eqn);
			deviceComplex nb = sellmeierFunc(ls, omega, &a[22], eqn);
			deviceComplex nc = sellmeierFunc(ls, omega, &a[44], eqn);
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
			deviceComplex ina2 = 1./ (na * na);
			deviceComplex inb2 = 1./(nb * nb);
			deviceComplex inc2 = 1./(nc * nc);
			double delta = 0.5 * atan(-((1. / realna2 - 1. / realnb2)
				* sin(2 * phi) * cosTheta) / ((cosPhi2 / realna2 + sinPhi2 / realnb2)
					+ ((sinPhi2 / realna2 + cosPhi2 / realnb2)
						* cosTheta2 + sinTheta2 / (nc.real() * nc.real()))));
			double cosDelta = cos(delta);
			double sinDelta = sin(delta);
			*ne = 1.0 / deviceLib::sqrt(cosDelta * cosDelta * (cosTheta2 * (cosPhi2 * ina2
				+ sinPhi2 * inb2) + sinTheta2  * inc2)
				+ sinDelta * sinDelta * (sinPhi2 * ina2 + cosPhi2 * inb2)
				- 0.5 * sin(2 * phi) * cosTheta * sin(2 * delta) * (ina2 - inb2));

			*no = 1.0 / deviceLib::sqrt(sinDelta * sinDelta * (cosTheta2 * (cosPhi2 * ina2
				+ sinPhi2 * inb2) + sinTheta2 * inc2)
				+ cosDelta * cosDelta * (sinPhi2 * ina2 + cosPhi2 * inb2)
				+ 0.5 * sin(2 * phi) * cosTheta * sin(2 * delta) * (ina2 - inb2));

			return *ne;
		}
	}

	deviceFunction double cuCModSquared(deviceComplex& a) {
		return a.real() * a.real() + a.imag() * a.imag();
	}

	//provide a list of nearest-3 neighbors for taking spatial derivatives
	// exploiting the fact that the radial grid is offset by 1/4 step from 0
	// this means that midpoints are available on the other side of the origin.
	// returns rho at the given index j
	deviceFunction double resolveNeighborsInOffsetRadialSymmetry(
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
	deviceFunction void findBirefringentCrystalIndex(cudaParameterSet* s, double* sellmeierCoefficients, long long i, deviceComplex* n1, deviceComplex* n2) {
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

		double gradientStep = 1.0e-6;
		double gradientFactor = 0.5 / gradientStep;
		int it;
		int maxiter = 64;
		double gradientTol = 1e-1;
		//emperical testing: 
		// converges to double precision limit in two iterations for BBO
		// converges in 32 iterations in BiBO

		double errArray[4][2];
		if ((*s).axesNumber == 1) {
			maxiter = 32;
			sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] + gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&n[1][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] - gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
			if (isnan(n[0][0].real()) || isnan(n[0][0].imag()) || isnan(n[1][0].real()) || isnan(n[1][0].imag())) {
				*n1 = deviceComplex(0.0, 0.0);
				*n2 = deviceComplex(0.0, 0.0);
				return;
			}
			errArray[0][0] = sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
			errArray[1][0] = sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
			gradient[0][0] = gradientFactor * (errArray[0][0] - errArray[1][0]);

			for (it = 0; it < maxiter; it++) {
				if (isnan(n[0][0].real()) || isnan(n[0][0].imag()) || isnan(n[1][0].real()) || isnan(n[1][0].imag())) {
					*n1 = deviceComplex(0.0, 0.0);
					*n2 = deviceComplex(0.0, 0.0);
					return;
				}
				if (abs(gradient[0][0]) > gradientTol) {
					alpha[0] -= 0.5 * (errArray[0][0] + errArray[1][0]) / gradient[0][0];
				} 
				else {
					break;
				}

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
			if (isnan(n[0][0].real()) || isnan(n[0][0].imag()) || isnan(n[1][0].real()) || isnan(n[1][0].imag())) {
				*n1 = n[0][0];
				*n2 = n[0][1];
				return;
			}
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
				if (isnan(n[0][0].real()) || isnan(n[0][0].imag()) || isnan(n[1][0].real()) || isnan(n[1][0].imag())) {
					*n1 = deviceComplex(0.0, 0.0);
					*n2 = deviceComplex(0.0, 0.0);
					return;
				}
				if (abs(gradient[0][0]) > 1e-2) alpha[0] -= 0.25 * (errArray[0][0] + errArray[1][0]) / gradient[0][0];
				if (abs(gradient[1][0]) > 1e-2) beta[0] -= 0.25 * (errArray[2][0] + errArray[3][0]) / gradient[1][0];
				if (abs(gradient[0][1]) > 1e-2) alpha[1] -= 0.25 * (errArray[0][1] + errArray[1][1]) / gradient[0][1];
				if (abs(gradient[1][1]) > 1e-2) beta[1] -= 0.25 * (errArray[2][1] + errArray[3][1]) / gradient[1][1];

				if (maxN(maxN(abs(gradient[0][0]), abs(gradient[1][0])), maxN(abs(gradient[0][1]), abs(gradient[1][1]))) < gradientTol) break;
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

	deviceFunction void findBirefingentCrystalAngle(double* alphaE, double* alphaO, unsigned long long j, double f, double* sellmeierCoefficients, cudaParameterSet* s) {
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
using namespace deviceFunctions;

namespace kernels {

	trilingual millersRuleNormalizationKernel asKernel(withID cudaParameterSet* s, double* sellmeierCoefficients, double* referenceFrequencies) {
		if (!(*s).isUsingMillersRule) {
			return;
		}

		double chi11[7];
		deviceComplex ne, no;
		for (int i = 0; i < 7; i++) {
			if (referenceFrequencies[i] == 0) {
				chi11[i] = 100000.0;
			}
			else {
				sellmeierCuda(&ne, &no, sellmeierCoefficients, referenceFrequencies[i], sellmeierCoefficients[66], sellmeierCoefficients[67], (int)sellmeierCoefficients[68], (int)sellmeierCoefficients[69]);
				chi11[i] = ne.real() * ne.real() - 1.0;
			}
		}

		//normalize chi2 tensor values
		for (int i = 0; i < 18; i++) {
			(*s).chi2Tensor[i] /= chi11[0] * chi11[1] * chi11[2];
		}

		//normalize chi3 tensor values
		for (int i = 0; i < 81; i++) {
			(*s).chi3Tensor[i] /= chi11[3] * chi11[4] * chi11[5] * chi11[6];
		}
	};

	trilingual totalSpectrumKernel asKernel(withID deviceComplex* fieldGrid1, deviceComplex* fieldGrid2, double gridStep, size_t Ntime, size_t Nspace, double* spectrum) {
		size_t i = localIndex;
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
	};

	trilingual totalSpectrum3DKernel asKernel(withID deviceComplex* fieldGrid1, deviceComplex* fieldGrid2, double gridStep, size_t Ntime, size_t Nspace, double* spectrum) {
		size_t i = localIndex;
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
	};

	//rotate the field around the propagation axis (basis change)
	trilingual rotateFieldKernel asKernel(withID deviceComplex* Ein1, deviceComplex* Ein2, deviceComplex* Eout1,
		deviceComplex* Eout2, double rotationAngle) {
		long long i = localIndex;
		Eout1[i] = cos(rotationAngle) * Ein1[i] - sin(rotationAngle) * Ein2[i];
		Eout2[i] = sin(rotationAngle) * Ein1[i] + cos(rotationAngle) * Ein2[i];
	};



	trilingual radialLaplacianKernel asKernel(withID cudaParameterSet* s) {
		unsigned long long i = localIndex;
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

	};
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
	trilingual expandCylindricalBeam asKernel(withID cudaParameterSet* s, double* polarization1, double* polarization2) {
		size_t i = localIndex;
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
	};

	//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
	//note that the sellmeier coefficients have extra values appended to the end
	//to give info about the current simulation
	trilingual applyFresnelLossKernel asKernel(withID double* sellmeierCoefficients1, double* sellmeierCoefficients2, cudaParameterSet* s) {
		long long i = localIndex;
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
			crystalTheta + 0 * alpha1, crystalPhi, axesNumber, sellmeierType);
		sellmeierCuda(&n0, &no1, sellmeierCoefficients1, f,
			crystalTheta + 0 * alphaO1, crystalPhi, axesNumber, sellmeierType);
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
	};

	trilingual apertureKernel asKernel(withID cudaParameterSet* s, double radius, double activationParameter) {
		long long i = localIndex;
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
			r = abs((*s).dx * ((double)j - (*s).Nspace / 2.0) + 0.25 * (*s).dx);
		}

		double a = 1.0 - (1.0 / (1.0 + exp(-activationParameter * (r - radius) / (*s).dx)));

		//if (r>radius) a = 0;
		(*s).gridETime1[i] *= a;
		(*s).gridETime2[i] *= a;
	};

	trilingual parabolicMirrorKernel asKernel(withID cudaParameterSet* s, double focus) {
		long long i = localIndex;
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
			r = abs((*s).dx * ((double)j - (*s).Nspace / 2.0) + 0.25 * (*s).dx);
		}

		deviceComplex	u = deviceLib::exp(deviceComplex(0.0,
			w * r * r * (0.5 / focus) / LIGHTC));

		(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
	};

	trilingual sphericalMirrorKernel asKernel(withID cudaParameterSet* s, double ROC) {
		long long i = localIndex;
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
			r = abs((*s).dx * ((double)j - (*s).Nspace / 2.0) + 0.25 * (*s).dx);
		}

		bool isNegative = signbit(ROC);
		ROC = abs(ROC);
		deviceComplex u = deviceComplex(0.0, 0.0);
		if (r < ROC) {
			u = deviceLib::exp(deviceComplex(0.0,
				2.0 * pow(-1, isNegative) * w * ROC * ((sqrt(1.0 - r * r / (ROC * ROC))) - 1.0) / LIGHTC));
		}

		(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
	};

	trilingual applyLinearPropagationKernel asKernel(withID double* sellmeierCoefficients, double thickness, cudaParameterSet* s) {
		long long i = localIndex;
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
		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);
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
	};


	//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
	//note that the sellmeier coefficients have extra values appended to the end
	//to give info about the current simulation
	trilingual prepareCartesianGridsKernel asKernel(withID double* sellmeierCoefficients, cudaParameterSet* s) {
		long long i = localIndex;
		long long j, k;
		deviceComplex ne, no;
		deviceComplex n0 = (*s).n0;
		deviceComplex cuZero = deviceComplex(0, 0);
		j = i / ((*s).Nfreq - 1); //spatial coordinate
		k = 1 + (i % ((*s).Nfreq - 1)); //temporal coordinate
		i = k + j * (*s).Nfreq;
		deviceComplex ii = deviceComplex(0, 1);
		double kStep = sellmeierCoefficients[70];
		double fStep = sellmeierCoefficients[71];

		//frequency being resolved by current thread
		double f = -k * fStep;

		//transverse wavevector being resolved
		double dk = j * kStep - (j >= ((long long)(*s).Nspace / 2)) * (kStep * (*s).Nspace); //frequency grid in transverse direction
		//sellmeierCuda(&n0, &no, sellmeierCoefficients, abs((*s).f0),
		//	crystalTheta, crystalPhi, axesNumber, sellmeierType);
		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);

		//if the refractive index was returned weird, then the index isn't valid, so set the propagator to zero for that frequency
		if (ne.real() < 0.9) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			return;
		}

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
		
		if (dk*dk < minN(ke.real() * ke.real() + ke.imag()*ke.imag(), ko.real() * ko.real() + ko.imag() * ko.imag()) 
			&& k < ((long long)(*s).Nfreq - 1)) {
			(*s).gridPropagationFactor1[i] = ii * (ke - k0 - dk * dk / (2. * ke.real())) * (*s).h;
			if (isnan(((*s).gridPropagationFactor1[i]).real())) {
				(*s).gridPropagationFactor1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = ii * (ko - k0 - dk * dk / (2. * ko.real())) * (*s).h;
			if (isnan(((*s).gridPropagationFactor2[i]).real())) {
				(*s).gridPropagationFactor2[i] = cuZero;
			}

			(*s).gridPolarizationFactor1[i] = ii * pow((*s).chiLinear1[k] + 1.0, 0.25) * chi11 * (TWOPI * f) / (2. * ne.real() * LIGHTC) * (*s).h;
			(*s).gridPolarizationFactor2[i] = ii * pow((*s).chiLinear2[k] + 1.0, 0.25) * chi12 * (TWOPI * f) / (2. * no.real() * LIGHTC) * (*s).h;
		}
		else {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
		}
	};

	//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
	//note that the sellmeier coefficients have extra values appended to the end
	//to give info about the current simulation
	trilingual prepare3DGridsKernel asKernel(withID double* sellmeierCoefficients, cudaParameterSet* s) {
		long long i = localIndex;
		long long col, j, k, l;
		deviceComplex ne, no;
		deviceComplex n0 = (*s).n0;
		deviceComplex cuZero = deviceComplex(0, 0);
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;
		k = col % (*s).Nspace;
		l = col / (*s).Nspace;

		deviceComplex ii = deviceComplex(0, 1);

		//frequency being resolved by current thread
		double f = -j * (*s).fStep;

		//transverse wavevector being resolved
		double dk1 = k * (*s).dk1 - (k >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (long long)(*s).Nspace); //frequency grid in x direction
		double dk2 = l * (*s).dk2 - (l >= ((long long)(*s).Nspace2 / 2)) * ((*s).dk2 * (long long)(*s).Nspace2); //frequency grid in y direction

		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);

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

		if (maxN(abs(dk1), abs(dk2)) < deviceLib::abs(ke) && j < ((long long)(*s).Nfreq - 1)) {
			(*s).gridPropagationFactor1[i] = ii * (ke - k0 - (dk1 * dk1 + dk2 * dk2) / (2. * ke.real())) * (*s).h;
			if (isnan(((*s).gridPropagationFactor1[i].real()))) {
				(*s).gridPropagationFactor1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = ii * (ko - k0 - (dk1 * dk1 + dk2 * dk2) / (2. * ko.real())) * (*s).h;
			if (isnan(((*s).gridPropagationFactor2[i].real()))) {
				(*s).gridPropagationFactor2[i] = cuZero;
			}

			(*s).gridPolarizationFactor1[i] = ii * pow((*s).chiLinear1[j] + 1.0, 0.25) * chi11 * (TWOPI * f) / (2. * ne.real() * LIGHTC) * (*s).h;
			(*s).gridPolarizationFactor2[i] = ii * pow((*s).chiLinear2[j] + 1.0, 0.25) * chi12 * (TWOPI * f) / (2. * no.real() * LIGHTC) * (*s).h;
		}

		else {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
		}
	};

	trilingual getChiLinearKernel asKernel(withID cudaParameterSet* s, double* sellmeierCoefficients) {
		long long i = localIndex;
		int axesNumber = (*s).axesNumber;
		int sellmeierType = (*s).sellmeierType;
		deviceComplex cuZero = deviceComplex(0, 0);

		double crystalTheta = sellmeierCoefficients[66];
		double crystalPhi = sellmeierCoefficients[67];
		double fStep = sellmeierCoefficients[71];

		deviceComplex ne, no;

		//frequency being resolved by current thread
		double f = i * fStep;
		
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
		(*s).fieldFactor1[i] = 1.0 / pow((*s).chiLinear1[i].real() + 1.0, 0.25); //account for the effective field strength in the medium (1/n)
		(*s).fieldFactor2[i] = 1.0 / pow((*s).chiLinear2[i].real() + 1.0, 0.25);
		if ((*s).isUsingMillersRule) {
			(*s).fieldFactor1[i] *= (*s).chiLinear1[i].real();
			(*s).fieldFactor2[i] *= (*s).chiLinear2[i].real();
		}

		if (i == 81) {
			deviceComplex n0;
			(*s).n0 = sellmeierCuda(&n0, &no, sellmeierCoefficients, abs((*s).f0), crystalTheta, crystalPhi, axesNumber, sellmeierType);
			(*s).chiLinear1[(*s).Ntime / 2] = deviceComplex(1.0, 0.0);
			(*s).chiLinear2[(*s).Ntime / 2] = deviceComplex(1.0, 0.0);
			(*s).fieldFactor1[(*s).Ntime / 2] = 0.0;
			(*s).fieldFactor2[(*s).Ntime / 2] = 0.0;
			(*s).inverseChiLinear2[(*s).Ntime / 2] = 1.0 / (*s).chiLinear2[i].real();
			(*s).inverseChiLinear2[(*s).Ntime / 2] = 1.0 / (*s).chiLinear2[i].real();
		}

		//apply Miller's rule to nonlinear coefficients
			if (!(*s).isUsingMillersRule || i > 80) {
				return;
			}
			double* referenceFrequencies = &sellmeierCoefficients[72];
			double chi11[7];

			for (int im = (i>17)*3; im < 7; im++) {
				if (referenceFrequencies[im] == 0) {
					chi11[im] = 100000.0;
				}
				else {
					sellmeierCuda(&ne, &no, sellmeierCoefficients, referenceFrequencies[im], sellmeierCoefficients[66], sellmeierCoefficients[67], (int)sellmeierCoefficients[68], (int)sellmeierCoefficients[69]);
					chi11[im] = ne.real() * ne.real() - 1.0;
				}
			}

			//normalize chi2 tensor values
			if (i < 18) {
				(*s).chi2Tensor[i] /= chi11[0] * chi11[1] * chi11[2];
			}

			//normalize chi3 tensor values
			(*s).chi3Tensor[i] /= chi11[3] * chi11[4] * chi11[5] * chi11[6];
	};
	//prepare the propagation constants under the assumption of cylindrical symmetry of the beam
	trilingual prepareCylindricGridsKernel asKernel(withID double* sellmeierCoefficients, cudaParameterSet* s) {
		long long i = localIndex;
		long long j, k;
		long long Nspace = (*s).Nspace;
		int axesNumber = (*s).axesNumber;
		int sellmeierType = (*s).sellmeierType;
		deviceComplex cuZero = deviceComplex(0, 0);
		j = i / ((*s).Nfreq - 1); //spatial coordinate
		k = 1 + i % ((*s).Nfreq - 1); //temporal coordinate
		i = k + j * (*s).Nfreq;


		deviceComplex ii = deviceComplex(0, 1);
		double crystalTheta = sellmeierCoefficients[66];
		double crystalPhi = sellmeierCoefficients[67];
		double kStep = sellmeierCoefficients[70];
		double fStep = sellmeierCoefficients[71];

		deviceComplex ne, no;
		deviceComplex n0 = (*s).n0;

		//frequency being resolved by current thread
		double f = -k * fStep;

		//transverse wavevector being resolved
		double dk = j * kStep - (j >= (Nspace / 2)) * (kStep * Nspace); //frequency grid in transverse direction

		sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta, crystalPhi, axesNumber, sellmeierType);

		//if the refractive index was returned weird, then the index isn't valid, so set the propagator to zero for that frequency
		if (ne.real() < 0.9) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			return;
		}
		
		
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

		if (abs(dk) <= minN(deviceLib::abs(ke), deviceLib::abs(ko)) && k < ((long long)(*s).Nfreq - 1)) {
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
			(*s).gridPolarizationFactor1[i] = 0.5 * pow((*s).chiLinear1[k] + 1.0, 0.25) * chi11 * ii * (TWOPI * f) / (2. * ne.real() * LIGHTC) * (*s).h;
			(*s).gridPolarizationFactor2[i] = 0.5 * pow((*s).chiLinear2[k] + 1.0, 0.25) * chi12 * ii * (TWOPI * f) / (2. * no.real() * LIGHTC) * (*s).h;


		}

		else {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor1Rho2[i] = cuZero;
		}
	};


	trilingual realToComplexKernel asKernel(withID double* in, deviceComplex* out) {
		long long i = localIndex;
		out[i] = deviceComplex(in[i], 0.0);
	};

	trilingual complexToRealKernel asKernel(withID deviceComplex* in, double* out) {
		long long i = localIndex;
		out[i] = in[i].real();
	};

	trilingual materialPhaseKernel asKernel(withID double df, size_t Ntime, double* a, double f01, double f02,
		double thickness1, double thickness2, double* phase1, double* phase2) {
		size_t i = localIndex;
		//frequency being resolved by current thread
		double f = i * df;
		if (i >= Ntime / 2) {
			f -= df * Ntime;
		}

		//give phase shift relative to group velocity (approximated 
		// with low-order finite difference) so the pulse doesn't move
		deviceComplex ne, no, no0, n0p, n0m;
		sellmeierCuda(&ne, &no, a, abs(f), 0, 0, 0, 0);
		f *= TWOPI;
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
	};

	//calculate the nonlinear polarization, after FFT to get the field
	//in the time domain
	trilingual nonlinearPolarizationKernel asKernel(withID cudaParameterSet* s) {
		size_t i = localIndex;
		double Ex = (*s).fftNorm * (*s).gridETime1[i];
		double Ey = (*s).fftNorm * (*s).gridETime2[i];

		double Ex2 = Ex * Ex;
		double Ey2 = Ey * Ey;
		(*s).gridPolarizationTime1[i] = 0.;
		(*s).gridPolarizationTime2[i] = 0.;
		//rotate field into crystal frame
		double E3[3] = { (*s).rotationForward[0] * Ex + (*s).rotationForward[1] * Ey,
			(*s).rotationForward[3] * Ex + (*s).rotationForward[4] * Ey,
			(*s).rotationForward[6] * Ex + (*s).rotationForward[7] * Ey };

		if ((*s).nonlinearSwitches[0] == 1) {
			double P2[3] = { 0.0 };
			for (unsigned char a = 0; a < 3; a++) {
				P2[a] += (*s).chi2Tensor[0 + a] * E3[0] * E3[0];
				P2[a] += (*s).chi2Tensor[3 + a] * E3[1] * E3[1];
				P2[a] += (*s).chi2Tensor[6 + a] * E3[2] * E3[2];
				P2[a] += (*s).chi2Tensor[9 + a] * E3[1] * E3[2];
				P2[a] += (*s).chi2Tensor[12 + a] * E3[0] * E3[2];
				P2[a] += (*s).chi2Tensor[15 + a] * E3[0] * E3[1];
			}
			(*s).gridPolarizationTime1[i] += (*s).rotationBackward[0] * P2[0] + (*s).rotationBackward[1] * P2[1] + (*s).rotationBackward[2] * P2[2];
			(*s).gridPolarizationTime2[i] += (*s).rotationBackward[3] * P2[0] + (*s).rotationBackward[4] * P2[1] + (*s).rotationBackward[5] * P2[2];
		}

		//resolve the full chi3 matrix when (*s).nonlinearSwitches[1]==1
		if ((*s).nonlinearSwitches[1] == 1) {
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
	};


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
	trilingual plasmaCurrentKernel_twoStage_A asKernel(withID cudaParameterSet* s) {
		size_t i = localIndex;
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
	};

	trilingual plasmaCurrentKernel_twoStage_B asKernel(withID cudaParameterSet* s) {
		size_t j = localIndex;
		j *= (*s).Ntime;
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
	};

	trilingual updateKwithPolarizationKernel asKernel(withID cudaParameterSet* sP) {
		size_t i = localIndex;
		unsigned int h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		unsigned int j = i / ((*sP).Nfreq - 1); //spatial coordinate
		i = h + j * ((*sP).Nfreq);
		h += (j + ((*sP).isCylindric * (j > ((long long)(*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;

		(*sP).k1[i] += (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h];
		(*sP).k2[i] += (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h];
	};

	trilingual updateKwithPlasmaKernel asKernel(withID cudaParameterSet* sP) {
		size_t i = localIndex;
		unsigned int h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		unsigned int j = i / ((*sP).Nfreq - 1); //spatial coordinate
		i = h + j * ((*sP).Nfreq);

		deviceComplex jfac = deviceComplex(0, -1.0 / (h * (*sP).fStep));
		h += (j + ((*sP).isCylindric * (j > ((long long)(*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;

		(*sP).k1[i] += jfac * (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h] * (*sP).inverseChiLinear1[i % ((*sP).Nfreq)];
		(*sP).k2[i] += jfac * (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h] * (*sP).inverseChiLinear2[i % ((*sP).Nfreq)];
	};

	//Main kernel for RK4 propagation of the field
	trilingual rkKernel asKernel(withID cudaParameterSet* sP, uint8_t stepNumber) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate

		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		if (h == 1) {
			(*sP).k1[iC - 1] = deviceComplex(0., 0.);
			(*sP).gridEFrequency1[iC - 1] = deviceComplex(0., 0.);
			(*sP).gridEFrequency1Next1[iC - 1] = deviceComplex(0., 0.);
			(*sP).workspace1[iC - 1] = deviceComplex(0., 0.);
		}
		deviceComplex estimate1;

		if ((*sP).isCylindric) {
			(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		}

		//generate the estimates and do the weighted sum to get the grid at the next step
		//with weights determined by the step number
		switch (stepNumber) {
		case 0:
			estimate1 = (*sP).gridEFrequency1[iC] + 0.5 * (*sP).k1[iC];
			(*sP).gridEFrequency1Next1[iC] = SIXTH * (*sP).k1[iC] + (*sP).gridEFrequency1[iC];
			(*sP).workspace1[iC] = (*sP).fieldFactor1[h] * estimate1;
			(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * estimate1;
			break;
		case 1:
			estimate1 = (*sP).gridEFrequency1[iC] + 0.5 * (*sP).k1[iC];
			(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + THIRD * (*sP).k1[iC];
			(*sP).workspace1[iC] = (*sP).fieldFactor1[h] * estimate1;
			(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * estimate1;
			break;
		case 2:
			estimate1 = (*sP).gridEFrequency1[iC] + (*sP).k1[iC];
			(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + THIRD * (*sP).k1[iC];
			(*sP).workspace1[iC] = (*sP).fieldFactor1[h] * estimate1;
			(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * estimate1;
			break;
		case 3:
			(*sP).gridEFrequency1[iC] = (*sP).gridEFrequency1Next1[iC] + SIXTH * (*sP).k1[iC];
			(*sP).workspace1[iC] = (*sP).fieldFactor1[h] * (*sP).gridEFrequency1[iC];
			(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC];
			break;
		}
	};

	trilingual beamNormalizeKernel asKernel(withID cudaParameterSet* s, double* rawSum, double* pulse, double pulseEnergy) {
		size_t i = localIndex;
		double normFactor = sqrt(pulseEnergy / ((*s).Ntime * (*rawSum)));
		pulse[i] *= normFactor;
	};

	trilingual addDoubleArraysKernel asKernel(withID double* A, double* B) {
		size_t i = localIndex;
		A[i] += B[i];
	};

	trilingual beamGenerationKernel2D asKernel(withID deviceComplex* pulse, double* pulseSum, cudaParameterSet* s, double frequency, double bandwidth,
		int sgOrder, double cep, double delay, double gdd, double tod,
		bool hasLoadedField, deviceComplex* loadedField, double* materialPhase,
		double w0, double z0, double x0, double beamAngle,
		double polarizationAngle, double circularity,
		double* sellmeierCoefficients, double crystalTheta, double crystalPhi, int sellmeierType) {
		long long i = localIndex;
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
#elif defined(__APPLE__)
		std::atomic<double>* pulseSumAtomic = (std::atomic<double>*)pulseSum;
		double expected = pulseSumAtomic->load();
		while (!std::atomic_compare_exchange_weak(pulseSumAtomic, &expected, expected + pointEnergy));
#elif defined RUNONSYCL
		cl::sycl::atomic_ref<double, cl::sycl::memory_order::relaxed,cl::sycl::memory_scope::device> a(*pulseSum);
		a.fetch_add(pointEnergy);
#elif defined(NOFETCHADD)
		(*pulseSum) += pointEnergy; //YOLO
#else
		std::atomic<double>* pulseSumAtomic = (std::atomic<double>*)pulseSum;
		(*pulseSumAtomic).fetch_add(pointEnergy);
#endif
	};

	//note to self: please make a beamParameters struct
	trilingual beamGenerationKernel3D asKernel(withID deviceComplex* pulse, double* pulseSum, cudaParameterSet* s, double frequency, double bandwidth,
		int sgOrder, double cep, double delay, double gdd, double tod,
		bool hasLoadedField, deviceComplex* loadedField, double* materialPhase,
		double w0, double z0, double y0, double x0, double beamAngle, double beamAnglePhi,
		double polarizationAngle, double circularity,
		double* sellmeierCoefficients, double crystalTheta, double crystalPhi, int sellmeierType) {
		long long i = localIndex;
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
#elif defined(__APPLE__)
		std::atomic<double>* pulseSumAtomic = (std::atomic<double>*)pulseSum;
		double expected = pulseSumAtomic->load();
		while (!std::atomic_compare_exchange_weak(pulseSumAtomic, &expected, expected + pointEnergy));
#elif defined RUNONSYCL
		cl::sycl::atomic_ref<double, cl::sycl::memory_order::relaxed, cl::sycl::memory_scope::device> a(*pulseSum);
		a.fetch_add(pointEnergy);
#elif defined (NOFETCHADD)
		(*pulseSum) += pointEnergy; //YOLO
#else
		std::atomic<double>* pulseSumAtomic = (std::atomic<double>*)pulseSum;
		(*pulseSumAtomic).fetch_add(pointEnergy);
#endif
	};


	trilingual multiplyByConstantKernelD asKernel(withID double* A, double val) {
		long long i = localIndex;
		A[i] = val * A[i];
	};

	trilingual multiplicationKernelCompactVector asKernel(withID deviceComplex* A, deviceComplex* B, deviceComplex* C, cudaParameterSet* s) {
		long long i = localIndex;
		long long h = i % (*s).Nfreq; //temporal coordinate

		C[i] = A[h] * B[i];
	};

	trilingual multiplicationKernelCompact asKernel(withID deviceComplex* A, deviceComplex* B, deviceComplex* C) {
		long long i = localIndex;
		C[i] = A[i] * B[i];
	};
}
using namespace kernels;

namespace hostFunctions{
	class ticker {
		std::chrono::steady_clock::time_point time0;

	public:
		ticker() {
			time0 = std::chrono::high_resolution_clock::now();
		}
		void set() {
			time0 = std::chrono::high_resolution_clock::now();
		}
		void printElapsed() {
			std::chrono::steady_clock::time_point time1 = std::chrono::high_resolution_clock::now();
			printf("Elapsed time: %8.4lf ms\n",
				1.0e-3 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0).count()));
		}

		void printSet() {
			printElapsed();
			set();
		}
	};
	typedef dlib::matrix<double, 0, 1> column_vector;
	simulationParameterSet* fittingSet;
	
	int getTotalSpectrum(activeDevice& d) {
		simulationParameterSet* sCPU = d.cParams;
		cudaParameterSet* sc = d.dParams;
		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.fft((*sc).gridETime1, (*sc).workspace1, 2);
		if ((*sc).is3D) {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum3DKernel, (*sc).workspace1, (*sc).workspace2, (*sCPU).rStep, (*sCPU).Ntime / 2 + 1, (*sCPU).Nspace * (*sCPU).Nspace2, (*sc).gridPolarizationTime1);
		}
		else {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel, (*sc).workspace1, (*sc).workspace2, (*sCPU).rStep, (*sCPU).Ntime / 2 + 1, (*sCPU).Nspace, (*sc).gridPolarizationTime1);
		}
		d.deviceMemcpy((*sCPU).totalSpectrum, (*sc).gridPolarizationTime1, 3 * (*sCPU).Nfreq * sizeof(double), DeviceToHost);
		return 0;
	}
	
	int prepareElectricFieldArrays(activeDevice& d) {

		simulationParameterSet* s = d.cParams;
		cudaParameterSet* sc = d.dParams;
		cudaParameterSet* scDevice = d.dParamsDevice;
		//run the beam generation single-threaded on CPU to avoid race condition
		unsigned int beamBlocks = (*sc).Nblock / 2;
		unsigned int beamThreads = (*sc).Nthread;
		
		d.deviceMemcpy(d.dParamsDevice, sc, sizeof(cudaParameterSet), HostToDevice);
		if ((*s).isFollowerInSequence && !(*s).isReinjecting) {
			d.deviceMemcpy((*sc).gridETime1, (*s).ExtOut, 2 * (*s).Ngrid * sizeof(double), HostToDevice);
			d.fft((*sc).gridETime1, (*sc).gridEFrequency1, 0);
			//Copy the field into the temporary array
			d.deviceMemcpy((*sc).gridEFrequency1Next1, (*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), DeviceToDevice);

			if ((*sc).isUsingMillersRule) {
				d.deviceLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), (unsigned int)(2u * MIN_GRIDDIM), multiplicationKernelCompactVector, (*sc).chiLinear1, (*sc).gridEFrequency1Next1, (*sc).workspace1, scDevice);
			}
			else {
				d.deviceMemcpy((*sc).workspace1, (*sc).gridEFrequency1Next1, 2 * sizeof(deviceComplex) * (*sc).NgridC, DeviceToDevice);
			}

			d.deviceLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplicationKernelCompact, (*sc).gridPropagationFactor1, (*sc).gridEFrequency1Next1, (*sc).k1);
			d.deviceMemcpy((*sc).gridEFrequency1Next1, (*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), DeviceToDevice);

			return 0;
		}
		double* materialPhase1CUDA, * materialPhase2CUDA;
		deviceComplex* loadedField1, * loadedField2;

		d.deviceCalloc((void**)&loadedField1, (*sc).Ntime, sizeof(deviceComplex));
		d.deviceCalloc((void**)&loadedField2, (*sc).Ntime, sizeof(deviceComplex));

		//get the material phase
		double* materialCoefficientsCUDA, * sellmeierPropagationMedium;
		//NOTE TO SELF: add second phase material

		if ((*s).field1IsAllocated) {
			d.deviceMemcpy(loadedField1, (*s).loadedField1, (*s).Ntime * sizeof(deviceComplex), HostToDevice);
		}
		if ((*s).field2IsAllocated) {
			d.deviceMemcpy(loadedField2, (*s).loadedField2, (*s).Ntime * sizeof(deviceComplex), HostToDevice);
		}
		d.deviceCalloc((void**)&materialCoefficientsCUDA, 66, sizeof(double));
		d.deviceCalloc((void**)&sellmeierPropagationMedium, 66, sizeof(double));
		d.deviceCalloc((void**)&materialPhase1CUDA, (*s).Ntime, sizeof(double));
		d.deviceCalloc((void**)&materialPhase2CUDA, (*s).Ntime, sizeof(double));
		d.deviceMemcpy(materialCoefficientsCUDA, (*s).crystalDatabase[(*s).phaseMaterialIndex1].sellmeierCoefficients, 66 * sizeof(double), HostToDevice);
		d.deviceMemcpy(sellmeierPropagationMedium, (*s).crystalDatabase[(*s).materialIndex].sellmeierCoefficients, 66 * sizeof(double), HostToDevice);

		d.deviceLaunch((unsigned int)(*s).Ntime, 1, materialPhaseKernel, (*s).fStep, (*s).Ntime, materialCoefficientsCUDA, (*s).frequency1, (*s).frequency2, (*s).phaseMaterialThickness1, (*s).phaseMaterialThickness2, materialPhase1CUDA, materialPhase2CUDA);

		double* pulseSum = &materialCoefficientsCUDA[0];
		//calculate pulse 1 and store it in unused memory
		d.deviceMemset(pulseSum, 0, sizeof(double));
		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));

		if ((*sc).is3D) {
			d.deviceLaunch(beamBlocks, beamThreads, beamGenerationKernel3D,
				(*sc).workspace1, pulseSum, scDevice, (*s).frequency1, (*s).bandwidth1,
				(*s).sgOrder1, (*s).cephase1, (*s).delay1, (*s).gdd1, (*s).tod1,
				(*s).field1IsAllocated, loadedField1, materialPhase1CUDA, (*s).beamwaist1,
				(*s).z01, (*s).y01, (*s).x01, (*s).propagationAngle1, (*s).propagationAnglePhi1, (*s).polarizationAngle1, (*s).circularity1,
				sellmeierPropagationMedium, (*s).crystalTheta, (*s).crystalPhi, (*s).sellmeierType);
		}
		else {
			d.deviceLaunch(beamBlocks, beamThreads, beamGenerationKernel2D,
				(*sc).workspace1, pulseSum, scDevice, (*s).frequency1, (*s).bandwidth1,
				(*s).sgOrder1, (*s).cephase1, (*s).delay1, (*s).gdd1, (*s).tod1,
				(*s).field1IsAllocated, loadedField1, materialPhase1CUDA, (*s).beamwaist1,
				(*s).z01, (*s).x01, (*s).propagationAngle1, (*s).polarizationAngle1, (*s).circularity1,
				sellmeierPropagationMedium, (*s).crystalTheta, (*s).crystalPhi, (*s).sellmeierType);
		}

		d.fft((*sc).workspace1, (*sc).gridETime1, 3);

		d.deviceLaunch(2 * (*sc).Nblock, (*sc).Nthread, beamNormalizeKernel, scDevice, pulseSum, (*sc).gridETime1, (*s).pulseEnergy1);

		d.deviceMemcpy((*sc).gridEFrequency1Next1, (*sc).gridETime1, (*sc).Ngrid * 2 * sizeof(double), DeviceToDevice);

		//calculate pulse 2
		d.deviceMemset(pulseSum, 0, sizeof(double));
		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		if ((*sc).is3D) {
			d.deviceLaunch(beamBlocks, beamThreads, beamGenerationKernel3D,
				(*sc).workspace1, pulseSum, scDevice, (*s).frequency2, (*s).bandwidth2,
				(*s).sgOrder2, (*s).cephase2, (*s).delay2, (*s).gdd2, (*s).tod2,
				(*s).field2IsAllocated, loadedField2, materialPhase2CUDA, (*s).beamwaist2,
				(*s).z02, (*s).y02, (*s).x02, (*s).propagationAngle2, (*s).propagationAnglePhi2, (*s).polarizationAngle2, (*s).circularity2,
				sellmeierPropagationMedium, (*s).crystalTheta, (*s).crystalPhi, (*s).sellmeierType);
		}
		else {
			d.deviceLaunch(beamBlocks, beamThreads, beamGenerationKernel2D,
				(*sc).workspace1, pulseSum, scDevice, (*s).frequency2, (*s).bandwidth2,
				(*s).sgOrder2, (*s).cephase2, (*s).delay2, (*s).gdd2, (*s).tod2,
				(*s).field2IsAllocated, loadedField2, materialPhase2CUDA, (*s).beamwaist2,
				(*s).z02, (*s).x02, (*s).propagationAngle2, (*s).polarizationAngle2, (*s).circularity2,
				sellmeierPropagationMedium, (*s).crystalTheta, (*s).crystalPhi, (*s).sellmeierType);
		}

		d.fft((*sc).workspace1, (*sc).gridETime1, 3);
		d.deviceLaunch(2 * (*sc).Nblock, (*sc).Nthread, beamNormalizeKernel, scDevice, pulseSum, (*sc).gridETime1, (*s).pulseEnergy2);
		//add the pulses
		d.deviceLaunch(2 * (*sc).Nblock, (*sc).Nthread, addDoubleArraysKernel, (*sc).gridETime1, (double*)(*sc).gridEFrequency1Next1);
		if ((*s).isReinjecting) {
			d.deviceMemcpy((*sc).workspace1, (*s).ExtOut, 2 * (*s).Ngrid * sizeof(double), HostToDevice);
			d.deviceLaunch(2 * (*sc).Nblock, (*sc).Nthread, addDoubleArraysKernel, (*sc).gridETime1, (double*)(*sc).workspace1);
		}
		//fft onto frequency grid
		d.fft((*sc).gridETime1, (*sc).gridEFrequency1, 0);

		//Copy the field into the temporary array
		d.deviceMemcpy((*sc).gridEFrequency1Next1, (*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), DeviceToDevice);

		if ((*sc).isUsingMillersRule && !(*sc).forceLinear) {
			d.deviceLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplicationKernelCompactVector, (*sc).chiLinear1, (*sc).gridEFrequency1Next1, (*sc).workspace1, scDevice);
		}
		else {
			d.deviceMemcpy((*sc).workspace1, (*sc).gridEFrequency1Next1, 2 * sizeof(deviceComplex) * (*sc).NgridC, DeviceToDevice);
		}
		d.deviceLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplicationKernelCompact, (*sc).gridPropagationFactor1, (*sc).gridEFrequency1Next1, (*sc).k1);

		d.deviceMemcpy((*sc).gridEFrequency1Next1, (*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), DeviceToDevice);
		d.deviceFree(materialPhase1CUDA);
		d.deviceFree(materialPhase2CUDA);
		d.deviceFree(materialCoefficientsCUDA);
		d.deviceFree(sellmeierPropagationMedium);
		d.deviceFree(loadedField1);
		d.deviceFree(loadedField2);
		return 0;
	}

	int applyFresnelLoss(simulationParameterSet* s, int materialIndex1, int materialIndex2) {
		activeDevice d;
		cudaParameterSet sc;
		d.allocateSet(s, &sc);
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
		d.deviceCalloc((void**)&sellmeierCoefficients1, 74, sizeof(double));
		d.deviceCalloc((void**)&sellmeierCoefficients2, 74, sizeof(double));
		d.deviceMemcpy(sellmeierCoefficients1, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), HostToDevice);

		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).crystalDatabase[materialIndex2].sellmeierCoefficients, 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		d.deviceMemcpy(sellmeierCoefficients2, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), HostToDevice);
		d.deviceMemcpy(sc.gridEFrequency1, (*s).EkwOut, 2 * (*s).NgridC * sizeof(deviceComplex), HostToDevice);

		//transform final result
		d.fft(sc.gridEFrequency1, sc.gridETime1, 1);
		d.deviceLaunch(2 * sc.Nblock, sc.Nthread, multiplyByConstantKernelD, sc.gridETime1, 1.0 / sc.Ngrid);
		//copy the field arrays from the GPU to CPU memory
		d.deviceMemcpy((*s).ExtOut, sc.gridETime1, 2 * (*s).Ngrid * sizeof(double), DeviceToHost);
		d.deviceMemcpy((*s).EkwOut, sc.gridEFrequency1, 2 * (*s).Ngrid * sizeof(deviceComplex), DeviceToHost);

		d.deviceFree(sellmeierCoefficients1);
		d.deviceFree(sellmeierCoefficients2);
		d.deallocateSet(&sc);
		return 0;
	}

	int applyAperature(simulationParameterSet* sCPU, double diameter, double activationParameter) {
		cudaParameterSet s;
		activeDevice d;
		d.allocateSet(sCPU, &s);

		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);

		cudaParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(cudaParameterSet), HostToDevice);
		d.deviceLaunch(s.Nblock, s.Nthread, apertureKernel, sDevice, 0.5 * diameter, activationParameter);
		d.fft(s.gridETime1, s.gridEFrequency1, 0);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), DeviceToHost);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);
		getTotalSpectrum(d);
		d.deallocateSet(&s);

		return 0;
	}

	int applySphericalMirror(simulationParameterSet* sCPU, double ROC) {
		cudaParameterSet s;
		activeDevice d;
		d.allocateSet(sCPU, &s);

		cudaParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(cudaParameterSet), HostToDevice);

		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		d.fft(s.gridETime1, s.gridEFrequency1, 2);
		d.deviceLaunch(s.Nblock / 2, s.Nthread, sphericalMirrorKernel, sDevice, ROC);
		d.fft(s.gridEFrequency1, s.gridETime1, 3);
		d.deviceLaunch(2 * s.Nblock, s.Nthread, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ntime);
		d.fft(s.gridETime1, s.gridEFrequency1, 0);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), DeviceToHost);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);
		getTotalSpectrum(d);
		d.deallocateSet(&s);
		return 0;
	}

	int applyParabolicMirror(simulationParameterSet* sCPU, double focus) {
		cudaParameterSet s;
		activeDevice d;
		d.allocateSet(sCPU, &s);
		cudaParameterSet* sDevice = d.dParamsDevice;

		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		d.fft(s.gridETime1, s.gridEFrequency1, 2);
		d.deviceLaunch(s.Nblock / 2, s.Nthread, parabolicMirrorKernel, sDevice, focus);
		d.fft(s.gridEFrequency1, s.gridETime1, 3);
		d.deviceLaunch(2 * s.Nblock, s.Nthread, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ntime);
		d.fft(s.gridETime1, s.gridEFrequency1, 0);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), DeviceToHost);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);
		getTotalSpectrum(d);
		d.deallocateSet(&s);
		return 0;
	}

	int applyLinearPropagation(simulationParameterSet* sCPU, int materialIndex, double thickness) {
		cudaParameterSet s;
		activeDevice d;
		d.allocateSet(sCPU, &s);

		d.deviceMemcpy(s.gridEFrequency1, (*sCPU).EkwOut, s.NgridC * 2 * sizeof(deviceComplex), HostToDevice);

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
		d.deviceMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), HostToDevice);
		s.axesNumber = (*sCPU).crystalDatabase[materialIndex].axisType;
		s.sellmeierType = (*sCPU).crystalDatabase[materialIndex].sellmeierType;
		cudaParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(cudaParameterSet), HostToDevice);

		d.deviceLaunch(s.Nblock / 2, s.Nthread, applyLinearPropagationKernel, sellmeierCoefficients, thickness, sDevice);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, s.NgridC * 2 * sizeof(deviceComplex), DeviceToHost);
		d.fft(s.gridEFrequency1, s.gridETime1, 1);
		d.deviceLaunch(2 * s.Nblock, s.Nthread, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);

		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), DeviceToHost);

		d.deallocateSet(&s);

		return 0;
	}

	int preparePropagationGrids(activeDevice& d) {
		cudaParameterSet* sc = d.dParams;
		simulationParameterSet* s = d.cParams;
		double* sellmeierCoefficients = (double*)(*sc).gridEFrequency1Next1;
		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[79];
		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).sellmeierCoefficients, 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		memcpy(sellmeierCoefficientsAugmentedCPU + 72, (*s).crystalDatabase[(*s).materialIndex].nonlinearReferenceFrequencies, 7 * sizeof(double));
		d.deviceMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, 79 * sizeof(double), HostToDevice);

		//prepare the propagation grids
		cudaParameterSet* sD = d.dParamsDevice;
		d.deviceMemcpy(sD, sc, sizeof(cudaParameterSet), HostToDevice);
		d.deviceLaunch((unsigned int)(*sc).Ntime/(2*MIN_GRIDDIM), MIN_GRIDDIM, getChiLinearKernel, sD, sellmeierCoefficients);
		if ((*s).is3D) {
			d.deviceLaunch((unsigned int)(*sc).Nblock / 2u, (unsigned int)(*sc).Nthread, prepare3DGridsKernel, sellmeierCoefficients, sD);
		}
		else if ((*s).isCylindric) {
			d.deviceLaunch((unsigned int)(*sc).Nblock / 2u, (unsigned int)(*sc).Nthread, prepareCylindricGridsKernel, sellmeierCoefficients, sD);
		}
		else {
			d.deviceLaunch((unsigned int)(*sc).Nblock / 2u, (unsigned int)(*sc).Nthread, prepareCartesianGridsKernel, sellmeierCoefficients, sD);
		}
		d.deviceMemcpy(sc, sD, sizeof(cudaParameterSet), DeviceToHost);
		return 0;
	}

	//Rotate the field on the GPU
	//Allocates memory and copies from CPU, then copies back to CPU and deallocates
	// - inefficient but the general principle is that only the CPU memory is preserved
	// after simulations finish... and this only runs at the end of the simulation
	int rotateField(simulationParameterSet* s, double rotationAngle) {
		cudaParameterSet sc;
		activeDevice d;
		d.allocateSet(s, &sc);
		deviceComplex* Ein1 = sc.gridEFrequency1;
		deviceComplex* Ein2 = sc.gridEFrequency2;
		deviceComplex* Eout1 = sc.gridEFrequency1Next1;
		deviceComplex* Eout2 = sc.gridEFrequency1Next2;
		//retrieve/rotate the field from the CPU memory
		d.deviceMemcpy(Ein1, (*s).EkwOut, 2 * (*s).NgridC * sizeof(deviceComplex), HostToDevice);
		d.deviceLaunch((unsigned int)(sc.NgridC / MIN_GRIDDIM), MIN_GRIDDIM, rotateFieldKernel, Ein1, Ein2, Eout1, Eout2, rotationAngle);
		d.deviceMemcpy((*s).EkwOut, Eout1, 2 * (*s).NgridC * sizeof(deviceComplex), DeviceToHost);

		//transform back to time
		d.fft(Eout1, sc.gridETime1, 1);
		d.deviceLaunch(2 * sc.Nblock, sc.Nthread, multiplyByConstantKernelD, sc.gridETime1, 1.0 / sc.Ngrid);
		d.deviceMemcpy((*s).ExtOut, sc.gridETime1, 2 * (*s).Ngrid * sizeof(double), DeviceToHost);

		//update spectrum
		getTotalSpectrum(d);

		d.deallocateSet(&sc);
		return 0;
	}

//function to run a RK4 time step
//stepNumber is the sub-step index, from 0 to 3
	int runRK4Step(activeDevice& d, uint8_t stepNumber) {
		cudaParameterSet* sH = d.dParams; 
		cudaParameterSet* sD = d.dParamsDevice;
		//operations involving FFT
		if ((*sH).isNonLinear || (*sH).isCylindric) {
			//perform inverse FFT to get time-space electric field
			d.fft((*sH).workspace1, (*sH).gridETime1, 1);

			//Nonlinear polarization
			if ((*sH).isNonLinear) {
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, nonlinearPolarizationKernel, sD);
				if ((*sH).isCylindric) {
					d.deviceLaunch((*sH).Nblock, (*sH).Nthread, expandCylindricalBeam, sD, (*sH).gridPolarizationTime1, (*sH).gridPolarizationTime2);
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, 4);
				}
				else {
					d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, 0);
				}
				d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPolarizationKernel, sD);
			}

			//Plasma/multiphoton absorption
			if ((*sH).hasPlasma) {
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, plasmaCurrentKernel_twoStage_A, sD);
				d.deviceLaunch((unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / MIN_GRIDDIM), 2*MIN_GRIDDIM, plasmaCurrentKernel_twoStage_B, sD);
				if ((*sH).isCylindric) {
					d.deviceLaunch((*sH).Nblock, (*sH).Nthread, expandCylindricalBeam, sD, (*sH).gridPolarizationTime1, (*sH).gridPolarizationTime2);
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, 4);
				}
				else {
					d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, 0);
				}
				d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPlasmaKernel, sD);
			}

			//Radial Laplacian
			if ((*sH).isCylindric) {
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, radialLaplacianKernel, sD);
				d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, 0);
			}
		}

		//advance an RK4 step
		d.deviceLaunch((*sH).Nblock, (*sH).Nthread, rkKernel, sD, stepNumber);
		return 0;
	}

	int resolveSequence(int currentIndex, simulationParameterSet * s, crystalEntry * db) {


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

			error = solveNonlinearWaveEquationX(s);

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
				solveNonlinearWaveEquationX(s);
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


		if ((*fittingSet).isInSequence) {
			solveNonlinearWaveEquationSequenceX(fittingSet);
		}
		else {
			solveNonlinearWaveEquationX(fittingSet);
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
			maxSim = maxN(maxSim, simSpec[i]);
			maxRef = maxN(maxRef, refSpec[i]);
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
}
using namespace hostFunctions;

unsigned long solveNonlinearWaveEquationX(void* lpParam) {
	simulationParameterSet* sCPU = (simulationParameterSet*)lpParam;
	size_t i;
	cudaParameterSet s;
	activeDevice d;
	memset(&s, 0, sizeof(cudaParameterSet));
	if (d.allocateSet(sCPU, &s)) return 1;

	//prepare the propagation arrays
	preparePropagationGrids(d);
	prepareElectricFieldArrays(d);

	double canaryPixel = 0;
	double* canaryPointer = &s.gridETime1[s.Ntime / 2 + s.Ntime * (s.Nspace / 2 + s.Nspace * (s.Nspace2 / 2))];
	d.deviceMemcpy(d.dParamsDevice, &s, sizeof(cudaParameterSet), HostToDevice);
	//Core propagation loop
	for (i = 0; i < s.Nsteps; i++) {

		//RK4
		runRK4Step(d, 0);
		runRK4Step(d, 1);
		runRK4Step(d, 2);
		runRK4Step(d, 3);

		//periodically check if the simulation diverged
		if (i % 8 == 0) {
#ifdef __CUDACC__
			cudaMemcpyAsync(&canaryPixel, canaryPointer, sizeof(double), DeviceToHost);
#else
			d.deviceMemcpy(&canaryPixel, canaryPointer, sizeof(double),DeviceToHost);
#endif
			if (isnan(canaryPixel)) {
				break;
			}
		}

		if ((*sCPU).imdone[0] == 2) {
			break;
		}

		if ((*sCPU).imdone[0] == 3) {
			//copy the field arrays from the GPU to CPU memory if requested by the UI
			d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);
			d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * (*sCPU).Ngrid * sizeof(deviceComplex), DeviceToHost);

			(*sCPU).imdone[0] = 0;
		}
		if(!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
	}
	if ((*sCPU).isInFittingMode && !(*sCPU).isInSequence)(*(*sCPU).progressCounter)++;

	//give the result to the CPU
	d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);


	d.fft(s.gridEFrequency1, s.gridETime1, 1);

	d.deviceLaunch((int)(s.Ngrid / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);
	d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);

	getTotalSpectrum(d);

	d.deallocateSet(&s);
	(*sCPU).imdone[0] = 1;
	return isnan(canaryPixel) * 13;
}


unsigned long solveNonlinearWaveEquationSequenceX(void* lpParam) {

	simulationParameterSet* sCPU = (simulationParameterSet*)lpParam;
	simulationParameterSet sCPUbackupValues;
	simulationParameterSet* sCPUbackup = &sCPUbackupValues;
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

		error = resolveSequence(k, sCPU, (*sCPU).crystalDatabase);


		if (error) break;
		memcpy(sCPU, sCPUbackup, sizeof(simulationParameterSet));
	}
	if ((*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
	return error;
}


unsigned long runDlibFittingX(simulationParameterSet* sCPU) {
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

	size_t fitCounter = 0;
	size_t* originalCounter = (*sCPU).progressCounter;
	(*sCPU).progressCounter = &fitCounter;

	if ((*sCPU).isInSequence) {
		solveNonlinearWaveEquationSequenceX(sCPU);
	}
	else {
		solveNonlinearWaveEquationX(sCPU);
	}

	(*sCPU).progressCounter = originalCounter;

	free(fittingSet);
	
	return 0;
}

int mainX(int argc, mainArgumentX){
	resolveArgv;
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
	simulationParameterSet initializationStruct;
	memset(&initializationStruct, 0, sizeof(simulationParameterSet));
	crystalEntry* crystalDatabasePtr = new crystalEntry[512];
	initializationStruct.crystalDatabase = crystalDatabasePtr;
	initializationStruct.progressCounter = &progressCounter;
	// read crystal database
	if (readCrystalDatabase(crystalDatabasePtr) == -2) {
		return 11;
	}
	if ((*crystalDatabasePtr).numberOfEntries == 0) {
		printf("Could not read crystal database.\n");
		delete[] crystalDatabasePtr;
		return 12;
	}
	printf("Read %i crystal database entries:\n", (*crystalDatabasePtr).numberOfEntries);
	for (j = 0; j < (*crystalDatabasePtr).numberOfEntries; j++) {
		printf("Material %i name: %ls", j, crystalDatabasePtr[j].crystalNameW);
	}

	// read from settings file
	if (readInputParametersFile(&initializationStruct, crystalDatabasePtr, filepath) == 1) {
		printf("Could not read input file.\n");
		delete[] crystalDatabasePtr;
		return 13;
	}
	simulationParameterSet* sCPU = new simulationParameterSet[initializationStruct.Nsims];
	memcpy(sCPU, &initializationStruct, sizeof(simulationParameterSet));
	allocateGrids(sCPU);
	if (loadPulseFiles(sCPU) == 1) {
		printf("Could not read pulse file.\n");
		deallocateGrids(sCPU, TRUE);
		delete[] sCPU;
		delete[] crystalDatabasePtr;
		return 14;
	}

	readSequenceString(sCPU);
	printf("Found %i steps in sequence\n", (*sCPU).Nsequence);
	configureBatchMode(sCPU);
	readFittingString(sCPU);
	auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

	if ((*sCPU).Nfitting != 0) {
		printf("Running optimization for %i iterations...\n", (*sCPU).fittingMaxIterations);

		runDlibFittingX(sCPU);

		auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
		printf("Finished after %8.4lf s. \n",
			1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
		saveDataSet(sCPU, crystalDatabasePtr, (*sCPU).outputBasePath, FALSE);

		printf("Optimization result:\n (index, value)\n");
		for (int i = 0; i < (*sCPU).Nfitting; i++) {
			printf("%i,  %lf\r\n", i, (*sCPU).fittingResult[i]);
		}

		deallocateGrids(sCPU, TRUE);
		delete[] sCPU;
		delete[] crystalDatabasePtr;
		return 0;
	}
	// run simulations
	std::thread* threadBlock = new std::thread[(*sCPU).Nsims * (*sCPU).Nsims2];
	size_t maxThreads = minN(CUDAdeviceCount, (*sCPU).Nsims * (*sCPU).Nsims2);
	for (j = 0; j < (*sCPU).Nsims * (*sCPU).Nsims2; j++) {

		sCPU[j].assignedGPU = j % CUDAdeviceCount;
		if (j >= maxThreads) {
			if (threadBlock[j - maxThreads].joinable()) {
				threadBlock[j - maxThreads].join();
			}
		}

		if ((*sCPU).isInSequence) {
			threadBlock[j] = std::thread(solveNonlinearWaveEquationSequenceX, &sCPU[j]);
		}
		else {
			threadBlock[j] = std::thread(solveNonlinearWaveEquationX, &sCPU[j]);
		}

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
	delete[] threadBlock;
	deallocateGrids(sCPU, TRUE);
	delete[] sCPU;
	delete[] crystalDatabasePtr;
	return 0;
}