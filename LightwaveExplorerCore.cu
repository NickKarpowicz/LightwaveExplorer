#include "LightwaveExplorerCore.cuh"
#include "LightwaveExplorerUtilities.h"
#include "LightwaveExplorerTrilingual.h"
#include <stdlib.h>
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

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
		double omega2 = omega * omega;
		double realPart;
		deviceComplex compPart;
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
			return deviceLib::sqrt(maxN(realPart, 0.9) + KLORENTZIAN * compPart);
		case 1:
			compPart = a[0] / deviceComplex(a[1] - omega2, a[2] * omega)
				+ a[3] / deviceComplex(a[4] - omega2, a[5] * omega)
				+ a[6] / deviceComplex(a[7] - omega2, a[8] * omega)
				+ a[9] / deviceComplex(a[10] - omega2, a[11] * omega)
				+ a[12] / deviceComplex(a[13] - omega2, a[14] * omega)
				+ a[15] / deviceComplex(a[16] - omega2, a[17] * omega)
				+ a[18] / deviceComplex(a[19] - omega2, a[20] * omega);
			return deviceLib::sqrt(KLORENTZIAN * compPart);
		case 100:
			realPart = a[0]
				+ (a[1] + a[2] * ls) / (ls + a[3])
				+ (a[4] + a[5] * ls) / (ls + a[6])
				+ (a[7] + a[8] * ls) / (ls + a[9])
				+ (a[10] + a[11] * ls) / (ls + a[12])
				+ a[13] * ls
				+ a[14] * ls * ls
				+ a[15] * ls * ls * ls;
			//"real-valued equation has no business being < 1" - causality
			return deviceComplex(sqrt(maxN(realPart, 0.9)), 0.0);
		}
		
		return deviceComplex(1.0, 0.0);
	};

	//Sellmeier equation for refractive indicies
	deviceFunction deviceComplex sellmeierCuda(
		deviceComplex* ne, deviceComplex* no, double* a, double f, double theta, double phi, int type, int eqn) {
		if (f == 0) {
			*ne = deviceComplex(1.0, 0.0); 
			*no = deviceComplex(1.0, 0.0); 
			return deviceComplex(1.0, 0.0);
		} //exit immediately for f=0
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
			na *= na;
			nb *= nb;
			double cosT = cos(theta);
			cosT *= cosT;
			double sinT = sin(theta);
			sinT *= sinT;

			*ne = deviceLib::sqrt((na * nb) / (nb * cosT + na * sinT));

			return *ne;
		}
		//option 2: biaxial
		else {
			deviceComplex na = sellmeierFunc(ls, omega, a, eqn);
			na *= na;
			deviceComplex nb = sellmeierFunc(ls, omega, &a[22], eqn);
			nb *= nb;
			deviceComplex nc = sellmeierFunc(ls, omega, &a[44], eqn);
			nc *= nc;
			double cp = cos(phi);
			cp *= cp;
			double sp = sin(phi);
			sp *= sp;
			double ct = cos(theta);
			ct *= ct;
			double st = sin(theta);
			st *= st;

			*ne = deviceLib::sqrt(na * nb * nc /
				(na * nb * st + na * nc * sp * ct + nb * nc * cp * ct));
			*no = deviceLib::sqrt(na * nb /
				(na * cp + nb * sp));

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
	deviceFunction void findBirefringentCrystalIndex(deviceParameterSet* s, double* sellmeierCoefficients, long long i, deviceComplex* n1, deviceComplex* n2) {
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
		deviceComplex n[4][2] = { {deviceComplex(0.0, 0.0)} };
		deviceComplex nW = deviceComplex(0.0, 0.0);
		sellmeierCuda(&n[0][0], &n[0][1], sellmeierCoefficients, f, sellmeierCoefficients[66], sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
		if ((*s).axesNumber == 0) {
			*n1 = n[0][0];
			*n2 = n[0][1];
			return;
		}

		double gradient[2][2] = { {0.0} };
		double alpha[2] = { asin(kx1 / n[0][0].real()),asin(kx1 / n[0][1].real()) };
		double beta[2] = { asin(ky1 / n[0][0].real()),asin(ky1 / n[0][1].real()) };

		double gradientStep = 1.0e-6;
		double gradientFactor = 0.5 / gradientStep;
		int it = 0;
		int maxiter = 64;
		double gradientTol = 1e-3;
		//emperical testing: 
		// converges to double precision limit in two iterations for BBO
		// converges in 32 iterations in BiBO

		double errArray[4][2] = { {0.0} };
		if ((*s).axesNumber == 1) {
			maxiter = 64;
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

			for (it = 0; it < maxiter; ++it) {
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

			for (it = 0; it < maxiter; ++it) {
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
}
using namespace deviceFunctions;

namespace kernels {
	//the syntax might look a bit strange due to the "trilingual" mode of operation. In short:
	// CUDA and OpenMP are fine with function pointers being used to launch kernels, but SYCL
	// doesn't allow them. However, SYCL does allow named lambdas, which have a nearly identical
	// structure. The trilingual preprocessor definition handles the return type (void for CUDA,
	// auto for SYCL, and the askernel definition is empty for CUDA, but contains the []( syntax
	// to declare a lambda. The widthID definition gives an additional parameter for openMP and
	// SYCL threads to be passed their ID.
	// The function's closing } has to be followed by a ; to have valid syntax in SYCL.

	//adjust tensor values based on the input and output frequencies given in the crystal database
	trilingual millersRuleNormalizationKernel asKernel(withID deviceParameterSet* s, double* sellmeierCoefficients, double* referenceFrequencies) {
		if (!(*s).isUsingMillersRule) {
			return;
		}

		double chi11[7];
		deviceComplex ne, no;
		for (int i = 0; i < 7; ++i) {
			if (referenceFrequencies[i] == 0) {
				chi11[i] = 100000.0;
			}
			else {
				sellmeierCuda(&ne, &no, sellmeierCoefficients, referenceFrequencies[i], sellmeierCoefficients[66], sellmeierCoefficients[67], s->axesNumber, s->sellmeierType);
				chi11[i] = ne.real() * ne.real() - 1.0;
			}
		}

		//normalize chi2 tensor values
		for (int i = 0; i < 18; ++i) {
			(*s).chi2Tensor[i] /= chi11[0] * chi11[1] * chi11[2];
		}

		//normalize chi3 tensor values
		for (int i = 0; i < 81; ++i) {
			(*s).chi3Tensor[i] /= chi11[3] * chi11[4] * chi11[5] * chi11[6];
		}
	};

	//calculate the total power spectrum of the beam for the 2D modes. Note that for the cartesian one, it will be treated as
	//a round beam instead of an infinite plane wave in the transverse direction. Thus, the 2D Cartesian spectra are approximations.
	trilingual totalSpectrumKernel asKernel(withID deviceComplex* fieldGrid1, deviceComplex* fieldGrid2, double gridStep, size_t Ntime, size_t Nspace, double* spectrum) {
		size_t i = localIndex;
		size_t j;
		double beamCenter1 = 0.;
		double beamCenter2 = 0.;
		double beamTotal1 = 0.;
		double beamTotal2 = 0.;
		double a, x;

		//find beam centers
		for (j = 0; j < Nspace; ++j) {
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
		for (j = 0; j < Nspace; ++j) {
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

	//Calculate the power spectrum after a 3D propagation
	trilingual totalSpectrum3DKernel asKernel(withID deviceComplex* fieldGrid1, deviceComplex* fieldGrid2, double gridStep, size_t Ntime, size_t Nspace, double* spectrum) {
		size_t i = localIndex;
		size_t j;

		double beamTotal1 = 0.;
		double beamTotal2 = 0.;
		//Integrate total beam power
		beamTotal1 = 0.;
		beamTotal2 = 0.;
		for (j = 0; j < Nspace; ++j) {
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

	//calculate the extra term in the Laplacian encountered in cylindrical coordinates (1/r d/drho)
	trilingual radialLaplacianKernel asKernel(withID deviceParameterSet* s) {
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
	trilingual expandCylindricalBeam asKernel(withID deviceParameterSet* s) {
		long long i = localIndex;
		long long j = i / (*s).Ntime; //spatial coordinate
		long long k = i % (*s).Ntime; //temporal coordinate

		//positions on the expanded grid corresponding the the current index
		long long pos1 = 2 * ((*s).Nspace - j - 1) * (*s).Ntime + k;
		long long pos2 = (2 * j + 1) * (*s).Ntime + k;

		//reuse memory allocated for the radial Laplacian, casting complex double
		//to a 2x larger double real grid
		double* expandedBeam1 = (double*)(*s).gridRadialLaplacian1;
		double* expandedBeam2 = expandedBeam1 + 2 * (*s).Ngrid;

		expandedBeam1[pos1] = (*s).gridPolarizationTime1[i];
		expandedBeam1[pos2] = (*s).gridPolarizationTime1[i];
		expandedBeam2[pos1] = (*s).gridPolarizationTime2[i];
		expandedBeam2[pos2] = (*s).gridPolarizationTime2[i];
	};

	//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
	//note that the sellmeier coefficients have extra values appended to the end
	//to give info about the current simulation
	trilingual applyFresnelLossKernel asKernel(withID double* sellmeierCoefficients1, double* sellmeierCoefficients2, deviceParameterSet* s) {
		//the old version was hopelessly broken, make a new one from scratch.
	};

	trilingual apertureFarFieldKernel asKernel(withID deviceParameterSet* s, double radius, double activationParameter, double xOffset, double yOffset) {
		long long i = localIndex;
		long long col, j, k, l;
		deviceComplex cuZero = deviceComplex(0, 0);
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;
		k = col % (*s).Nspace;
		l = col / (*s).Nspace;

		//magnitude of k vector
		double ko = TWOPI*j * (*s).fStep/LIGHTC;

		//transverse wavevector being resolved
		double dk1 = k * (*s).dk1 - (k >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (long long)(*s).Nspace); //frequency grid in x direction
		double dk2 = 0.0;
		if((*s).is3D) dk2 = l * (*s).dk2 - (l >= ((long long)(*s).Nspace2 / 2)) * ((*s).dk2 * (long long)(*s).Nspace2); //frequency grid in y direction

		//light that won't go the the farfield is immediately zero
		if (dk1*dk1 > ko*ko || dk2*dk2 > ko*ko) {
			(*s).gridEFrequency1[i] = cuZero;
			(*s).gridEFrequency2[i] = cuZero;
			return;
		}

		double theta1 = asin(dk1 / ko);
		double theta2 = asin(dk2 / ko);

		theta1 -= (!(*s).isCylindric) *  xOffset;
		theta2 -= (*s).is3D * yOffset;

		double r = sqrt(theta1 * theta1 + theta2 * theta2);

		double a = 1.0 - (1.0 / (1.0 + exp(-activationParameter * (r-radius))));

		(*s).gridEFrequency1[i] *= a;
		(*s).gridEFrequency2[i] *= a;
	};


	//apply a spectral filter to the beam (full details in docs)
	trilingual filterKernel asKernel(withID deviceParameterSet* s, double f0, double bandwidth, double order, double inBandAmplitude, double outOfBandAmplitude) {
		long long i = localIndex;
		long long col, j;
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;

		double f = (*s).fStep * j - f0;
		for (int p = 1; p < (int)order; p++) {
			bandwidth *= bandwidth;
			f *= f;
		}
		double filterFunction = outOfBandAmplitude + inBandAmplitude*exp(-f / (2.0 * bandwidth));
		(*s).gridEFrequency1[i] *= filterFunction;
		(*s).gridEFrequency2[i] *= filterFunction;
	};

	//Apply a (soft, possibly) aperture
	trilingual apertureKernel asKernel(withID deviceParameterSet* s, double radius, double activationParameter) {
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

	//apply a spatial phase corresponding to a parabolic mirror (on-axis)
	trilingual parabolicMirrorKernel asKernel(withID deviceParameterSet* s, double focus) {
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

	//apply a spatial phase corresponding to a spherical mirror (on axis)
	trilingual sphericalMirrorKernel asKernel(withID deviceParameterSet* s, double ROC) {
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

		bool isNegative = ROC < 0.0;
		ROC = abs(ROC);
		deviceComplex u = deviceComplex(0.0, 0.0);
		if (r < ROC) {
			u = deviceLib::exp(deviceComplex(0.0,
				2.0 * pow(-1, isNegative) * w * ROC * ((sqrt(1.0 - r * r / (ROC * ROC))) - 1.0) / LIGHTC));
		}

		(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
	};

	//apply linear propagation through a given medium to the fields
	trilingual applyLinearPropagationKernel asKernel(withID double* sellmeierCoefficients, double thickness, deviceParameterSet* s) {
		long long i = localIndex;
		long long j, h, k, col;
		int axesNumber = (*s).axesNumber;
		int sellmeierType = (*s).sellmeierType;
		deviceComplex ne, no, n0, n0o;
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
	trilingual prepareCartesianGridsKernel asKernel(withID double* sellmeierCoefficients, deviceParameterSet* s) {
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
		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);

		//if the refractive index was returned weird, then the index isn't valid, so set the propagator to zero for that frequency
		if (minN(ne.real(), no.real()) < 0.9 || isnan(ne.real()) || isnan(no.real()) || isnan(ne.imag()) || isnan(no.imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			return;
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
	trilingual prepare3DGridsKernel asKernel(withID double* sellmeierCoefficients, deviceParameterSet* s) {
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
		if (minN(ne.real(), no.real()) < 0.9 || isnan(ne.real()) || isnan(no.real()) || isnan(ne.imag()) || isnan(no.imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			return;
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

	//prepare the chi(1) arrays that will be needed in the simulation
	trilingual getChiLinearKernel asKernel(withID deviceParameterSet* s, double* sellmeierCoefficients) {
		long long i = localIndex;
		int axesNumber = (*s).axesNumber;
		int sellmeierType = (*s).sellmeierType;
		double crystalTheta = sellmeierCoefficients[66];
		double crystalPhi = sellmeierCoefficients[67];
		double fStep = sellmeierCoefficients[71];

		deviceComplex ne, no;

		//frequency being resolved by current thread
		double f = i * fStep;
		
		sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta, crystalPhi, axesNumber, sellmeierType);
		if (isnan(ne.real()) || isnan(no.real())) {
			ne = deviceComplex(1.0, 0);
			no = ne;
		}

		(*s).chiLinear1[i] = ne * ne - 1.0;
		(*s).chiLinear2[i] = no * no - 1.0;
		if ((*s).chiLinear1[i].real() != 0.0 && (*s).chiLinear2[i].real() != 0.0) {
			(*s).inverseChiLinear1[i] = 1.0 / (*s).chiLinear1[i].real();
			(*s).inverseChiLinear2[i] = 1.0 / (*s).chiLinear2[i].real();
		}
		else {
			(*s).inverseChiLinear1[i] = 0.0;
			(*s).inverseChiLinear2[i] = 0.0;
		}

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
			(*s).inverseChiLinear1[(*s).Ntime / 2] = 1.0 / (*s).chiLinear1[i].real();
			(*s).inverseChiLinear2[(*s).Ntime / 2] = 1.0 / (*s).chiLinear2[i].real();
		}

		//apply Miller's rule to nonlinear coefficients
			if (!(*s).isUsingMillersRule || i > 80) {
				return;
			}
			double* referenceFrequencies = &sellmeierCoefficients[72];
			double chi11[7];

			for (int im = (i>17)*3; im < 7; ++im) {
				if (referenceFrequencies[im] == 0) {
					chi11[im] = 100000.0;
				}
				else {
					sellmeierCuda(&ne, &no, sellmeierCoefficients, referenceFrequencies[im], sellmeierCoefficients[66], sellmeierCoefficients[67], axesNumber, sellmeierType);
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
	trilingual prepareCylindricGridsKernel asKernel(withID double* sellmeierCoefficients, deviceParameterSet* s) {
		long long i = localIndex;
		long long j, k;
		long long Nspace = (*s).Nspace;
		deviceComplex cuZero = deviceComplex(0.0, 0.0);
		j = i / ((*s).Nfreq - 1); //spatial coordinate
		k = 1 + i % ((*s).Nfreq - 1); //temporal coordinate
		i = k + j * (*s).Nfreq;
		deviceComplex ii = deviceComplex(0.0, 1.0);
		double kStep = sellmeierCoefficients[70];
		double fStep = sellmeierCoefficients[71];

		deviceComplex ne, no;
		deviceComplex n0 = (*s).n0;

		//frequency being resolved by current thread
		double f = -k * fStep;

		//transverse wavevector being resolved
		double dk = j * kStep - (j >= (Nspace / 2)) * (kStep * Nspace); //frequency grid in transverse direction

		sellmeierCuda(&ne, &no, sellmeierCoefficients,fStep*k, sellmeierCoefficients[66], sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);

		//if the refractive index was returned weird, then the index isn't valid, so set the propagator to zero for that frequency
		if (minN(ne.real(), no.real()) < 1.0 || isnan(ne.real()) || isnan(no.real()) || isnan(ne.imag()) || isnan(no.imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			(*s).gridPropagationFactor1Rho1[i] = cuZero;
			(*s).gridPropagationFactor1Rho2[i] = cuZero;
			return;
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

		if ((dk * dk < minN(ke.real() * ke.real() + ke.imag() * ke.imag(), ko.real() * ko.real() + ko.imag() * ko.imag())) && (*s).fieldFactor1[k] > 0.0 && (*s).fieldFactor2[k] > 0.0) {
			(*s).gridPropagationFactor1[i] = ii * (ke - k0 - dk * dk / (2. * ke.real())) * (*s).h;
			(*s).gridPropagationFactor1Rho1[i] = deviceComplex(0.0, (*s).h / ((*s).fieldFactor1[k] * 2. * ke.real()));
			if (isnan(((*s).gridPropagationFactor1[i].real()))) {
				(*s).gridPropagationFactor1[i] = cuZero;
				(*s).gridPropagationFactor1Rho1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = ii * (ko - k0 - dk * dk / (2. * ko.real())) * (*s).h;
			(*s).gridPropagationFactor1Rho2[i] = deviceComplex(0.0, (*s).h / ((*s).fieldFactor2[k] * 2. * ko.real()));
			if (isnan(((*s).gridPropagationFactor2[i].real()))) {
				(*s).gridPropagationFactor2[i] = cuZero;
				(*s).gridPropagationFactor1Rho2[i] = cuZero;
			}

			//factor of 0.5 comes from doubled grid size in cylindrical symmetry mode after expanding the beam
			(*s).gridPolarizationFactor1[i] = 0.5 * pow((*s).chiLinear1[k] + 1.0, 0.25) * chi11 * ii * (TWOPI * f) / (2. * ne.real() * LIGHTC) * (*s).h;
			(*s).gridPolarizationFactor2[i] = 0.5 * pow((*s).chiLinear2[k] + 1.0, 0.25) * chi12 * ii * (TWOPI * f) / (2. * no.real() * LIGHTC) * (*s).h;
			if (isnan((*s).gridPolarizationFactor1[i].real()) || isnan((*s).gridPolarizationFactor1[i].imag()) || isnan((*s).gridPolarizationFactor2[i].real()) || isnan((*s).gridPolarizationFactor2[i].imag())) {
				(*s).gridPolarizationFactor1[i] = cuZero;
				(*s).gridPolarizationFactor2[i] = cuZero;
			}
		}

		else {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			(*s).gridPropagationFactor1Rho1[i] = cuZero;
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
	trilingual nonlinearPolarizationKernel asKernel(withID deviceParameterSet* s) {
		size_t i = localIndex;
		double Ex = (*s).gridETime1[i];
		double Ey = (*s).gridETime2[i];

		if ((*s).nonlinearSwitches[0] == 1 || (*s).nonlinearSwitches[1] == 1){
			
			double P[3] = {0.0};

			// rotate field into crystal frame
			double E3[3] = {(*s).rotationForward[0] * Ex + (*s).rotationForward[1] * Ey,
							(*s).rotationForward[3] * Ex + (*s).rotationForward[4] * Ey,
							(*s).rotationForward[6] * Ex + (*s).rotationForward[7] * Ey};

			// second order nonlinearity, element-by-element in the reduced tensor
			//note that for historical reasons, chi2 is column-order and chi3 is row-order...
			if ((*s).nonlinearSwitches[0] == 1){
				for (unsigned char a = 0; a < 3; ++a){
					P[a] += (*s).chi2Tensor[0 + a] * E3[0] * E3[0];
					P[a] += (*s).chi2Tensor[3 + a] * E3[1] * E3[1];
					P[a] += (*s).chi2Tensor[6 + a] * E3[2] * E3[2];
					P[a] += (*s).chi2Tensor[9 + a] * E3[1] * E3[2];
					P[a] += (*s).chi2Tensor[12 + a] * E3[0] * E3[2];
					P[a] += (*s).chi2Tensor[15 + a] * E3[0] * E3[1];
				}
			}

			// resolve the full chi3 matrix when (*s).nonlinearSwitches[1]==1
			if ((*s).nonlinearSwitches[1] == 1){
				// loop over tensor element X_abcd
				// i hope the compiler unrolls this, but no way am I writing that out by hand
				for (unsigned char a = 0; a < 3; ++a){
					for (unsigned char b = 0; b < 3; ++b){
						for (unsigned char c = 0; c < 3; ++c){
							for (unsigned char d = 0; d < 3; ++d){
								P[d] += (*s).chi3Tensor[a + 3 * b + 9 * c + 27 * d] * E3[a] * E3[b] * E3[c];
							}
						}
					}
				}
			}

			(*s).gridPolarizationTime1[i] = (*s).rotationBackward[0] * P[0] + (*s).rotationBackward[1] * P[1] + (*s).rotationBackward[2] * P[2];
			(*s).gridPolarizationTime2[i] = (*s).rotationBackward[3] * P[0] + (*s).rotationBackward[4] * P[1] + (*s).rotationBackward[5] * P[2];

			//using only one value of chi3, under assumption of centrosymmetry when nonlinearSwitches[1]==2
			if ((*s).nonlinearSwitches[1] == 2) {
				double Esquared = (*s).chi3Tensor[0] * (Ex*Ex + Ey*Ey);
				(*s).gridPolarizationTime1[i] += Ex * Esquared;
				(*s).gridPolarizationTime2[i] += Ey * Esquared;
			}
		}
		else{
			//case of no chi2, and centrosymmetric chi3
			if ((*s).nonlinearSwitches[1] == 2) {
				double Esquared = (*s).chi3Tensor[0] * (Ex*Ex + Ey*Ey);
				(*s).gridPolarizationTime1[i] = Ex * Esquared;
				(*s).gridPolarizationTime2[i] = Ey * Esquared;
			}
			else{
			//this should never be called: the simulation thinks there's a nonlinearity, but they're all off
			//zero just in case.
				(*s).gridPolarizationTime1[i] = 0.0;
				(*s).gridPolarizationTime2[i] = 0.0;
			}
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
	//applied in 3 parts due to the different patterns of calculating ionization rate (A)
	//and integrating over trajectories (B). Added to RK4 propagation array afterwards.
	trilingual plasmaCurrentKernel_twoStage_A asKernel(withID deviceParameterSet* s) {
		size_t i = localIndex;
		double Esquared, a;
		unsigned char pMax = (unsigned char)(*s).nonlinearSwitches[3];

		//save values in workspaces, casting to double
		double* dN = (double*)(*s).workspace1;
		double* dN2 = dN + (*s).Ngrid;
		double* Jx = (*s).gridPolarizationTime1;
		double* Jy = (*s).gridPolarizationTime2;
		Esquared = (*s).gridETime1[i] * (*s).gridETime1[i] + (*s).gridETime2[i] * (*s).gridETime2[i];
		a = (*s).plasmaParameters[0] * Esquared;
		for (unsigned char p = 0; p < pMax; ++p) {
			a *= Esquared;
		}
		Jx[i] = a * (*s).gridETime1[i];
		Jy[i] = a * (*s).gridETime2[i];
		dN[i] = (*s).plasmaParameters[2] * (Jx[i] * (*s).gridETime1[i] + Jy[i] * (*s).gridETime2[i]);
		dN2[i] = dN[i];
	};

	trilingual plasmaCurrentKernel_twoStage_B asKernel(withID deviceParameterSet* s) {
		size_t j = localIndex;
		j *= (*s).Ntime;
		double N = 0;
		double integralx = 0;
		double* expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
		double* dN = j + (double*)(*s).workspace1;
		double* E = &(*s).gridETime1[j];
		double* P = &(*s).gridPolarizationTime1[j];
		for (unsigned int k = 0; k < (*s).Ntime; ++k) {
			N += dN[k];
			integralx += N * (*s).expGammaT[k] * E[k];
			P[k] += expMinusGammaT[k] * integralx;
		}
	};

	trilingual updateKwithPolarizationKernel asKernel(withID deviceParameterSet* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		h += j * (*sP).Nfreq;
		(*sP).k1[h] += (*sP).gridPolarizationFactor1[h] * (*sP).workspace1[h];
		(*sP).k2[h] += (*sP).gridPolarizationFactor2[h] * (*sP).workspace2P[h];
	};

	trilingual updateKwithPolarizationKernelCylindric asKernel(withID deviceParameterSet* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		i = h + j * ((*sP).Nfreq);
		h += (j + ((j > ((*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;
		(*sP).k1[i] += (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h];
		(*sP).k2[i] += (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h];
	};

	trilingual updateKwithPlasmaKernel asKernel(withID deviceParameterSet* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		deviceComplex jfac = deviceComplex(0.0, -1.0 / (h * (*sP).fStep));
		h += j * (*sP).Nfreq;
		(*sP).k1[h] += jfac * (*sP).gridPolarizationFactor1[h] * (*sP).workspace1[h] * (*sP).inverseChiLinear1[h % ((*sP).Nfreq)];
		(*sP).k2[h] += jfac * (*sP).gridPolarizationFactor2[h] * (*sP).workspace2P[h] * (*sP).inverseChiLinear2[h % ((*sP).Nfreq)];
	};

	trilingual updateKwithPlasmaKernelCylindric asKernel(withID deviceParameterSet* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		i = h + j * ((*sP).Nfreq);
		deviceComplex jfac = deviceComplex(0.0, -1.0 / (h * (*sP).fStep));
		h += (j + ( (j > ((*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;
		(*sP).k1[i] += jfac * (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h] * (*sP).inverseChiLinear1[i % ((*sP).Nfreq)];
		(*sP).k2[i] += jfac * (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h] * (*sP).inverseChiLinear2[i % ((*sP).Nfreq)];
	};

	//Slightly different kernels for the four stages of RK4. They used to be one big kernel with a switch case
	//but this has slightly better utilization.
	trilingual rkKernel0 asKernel(withID deviceParameterSet* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);

		if (h == 1) (*sP).workspace1[iC - 1] = deviceComplex(0., 0.);
		if ((*sP).isCylindric) (*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];

		deviceComplex estimate1 = (*sP).gridEFrequency1[iC] + 0.5 * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = SIXTH * (*sP).k1[iC] + (*sP).gridEFrequency1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * (*sP).fieldFactor1[h] * estimate1;
		(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * estimate1;
	};

	trilingual rkKernel1 asKernel(withID deviceParameterSet* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);

		if (h == 1) (*sP).workspace1[iC - 1] = deviceComplex(0., 0.);
		if ((*sP).isCylindric) (*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];

		deviceComplex estimate1 = (*sP).gridEFrequency1[iC] + 0.5 * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + THIRD * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * (*sP).fieldFactor1[h] * estimate1;
		(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * estimate1;
	};

	trilingual rkKernel2 asKernel(withID deviceParameterSet* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);

		if (h == 1) (*sP).workspace1[iC - 1] = deviceComplex(0., 0.);
		if ((*sP).isCylindric) (*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];

		deviceComplex estimate1 = (*sP).gridEFrequency1[iC] + (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + THIRD * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * (*sP).fieldFactor1[h] * estimate1;
		(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * estimate1;
	};

	trilingual rkKernel3 asKernel(withID deviceParameterSet* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);

		if (h == 1) (*sP).workspace1[iC - 1] = deviceComplex(0., 0.);
		if ((*sP).isCylindric) (*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];

		(*sP).gridEFrequency1[iC] = (*sP).gridEFrequency1Next1[iC] + SIXTH * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * (*sP).fieldFactor1[h] * (*sP).gridEFrequency1[iC];
		(*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC];
	};

	trilingual beamNormalizeKernel asKernel(withID deviceParameterSet* s, double* rawSum, double* pulse, double pulseEnergy) {
		size_t i = localIndex;
		double normFactor = sqrt(pulseEnergy / ((*s).Ntime * (*rawSum)));
		pulse[i] *= normFactor;
	};

	trilingual addDoubleArraysKernel asKernel(withID double* A, double* B) {
		size_t i = localIndex;
		A[i] += B[i];
	};

	//crease a pulse on the grid for the 2D modes. Note that normalization of the 2D mode assumes radial symmetry (i.e. that it's
	//a gaussian beam, not an infinite plane wave, which would have zero amplitude for finite energy).
	trilingual beamGenerationKernel2D asKernel(withID deviceComplex* field, pulse* p, double* pulseSum, deviceParameterSet* s,
		bool hasLoadedField, deviceComplex* loadedField, double* materialPhase, double* sellmeierCoefficients) {
		long long i = localIndex;
		long long j, h;
		h = 1 + i % ((*s).Nfreq - 1);
		j = i / ((*s).Nfreq - 1);
		i = h + j * ((*s).Nfreq);
		double f = h * (*s).fStep;
		double w = TWOPI * (f - (*p).frequency);

		//supergaussian pulse spectrum, if no input pulse specified
		deviceComplex specfac = deviceComplex(-pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder), 0);

		deviceComplex specphase = deviceComplex(0,
			-((*p).cep
				+ TWOPI * f * ((*p).delay - 0.5 * (*s).dt * (*s).Ntime)
				+ 0.5 * (*p).gdd * w * w
				+ (*p).tod * w * w * w / 6.0
				+ materialPhase[h]));
		specfac = deviceLib::exp(specfac + specphase);

		if (hasLoadedField) {
			specfac = loadedField[h] * deviceLib::exp(specphase);
		}
		deviceComplex ne, no;
		sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), (*s).crystalTheta, (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);


		double ko = TWOPI * no.real() * f / LIGHTC;
		double zR = PI * (*p).beamwaist * (*p).beamwaist * ne.real() * f / LIGHTC;
		if (f == 0) {
			zR = 1e3;
		}
		double rB = ((*p).x0 - (*s).dx * (j - (*s).Nspace / 2.0) - 0.25 * (*s).dx);
		double r = rB * cos((*p).beamAngle) - (*p).z0 * sin((*p).beamAngle);
		double z = rB * sin((*p).beamAngle) + (*p).z0 * cos((*p).beamAngle);

		double wz = (*p).beamwaist * sqrt(1 + (z * z / (zR * zR)));
		double Rz = z * (1. + (zR * zR / (z * z)));

		if (z == 0) {
			Rz = 1.0e15;
		}
		double phi = atan(z / zR);
		deviceComplex Eb = ((*p).beamwaist / wz) * deviceLib::exp(deviceComplex(0., 1.) * (ko * (z - (*p).z0) + ko * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
		Eb = Eb * specfac;
		if (isnan(cuCModSquared(Eb)) || f <= 0) {
			Eb = deviceComplex(0., 0.);
		}

		field[i] = deviceComplex(cos((*p).polarizationAngle), -(*p).circularity * sin((*p).polarizationAngle)) * Eb;
		field[i + (*s).NgridC] = deviceComplex(sin((*p).polarizationAngle), (*p).circularity * cos((*p).polarizationAngle)) * Eb;
		double pointEnergy = abs(r) * (cuCModSquared(field[i]) + cuCModSquared(field[i + (*s).NgridC]));
		pointEnergy *= 2 * PI * LIGHTC * EPS0 * (*s).dx * (*s).dt;
		//two factors of two cancel here - there should be one for the missing frequency plane, but the sum is over x instead of r
		//accordingly we already multiplied by two
		atomicAddDevice(pulseSum, pointEnergy);
	};

	//Generate a beam in full 3D mode
	trilingual beamGenerationKernel3D asKernel(withID deviceComplex* field, pulse* p, double* pulseSum, deviceParameterSet* s,
		bool hasLoadedField, deviceComplex* loadedField, double* materialPhase, double* sellmeierCoefficients) {
		long long i = localIndex;
		long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		i = h + col * ((*s).Nfreq);
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;
		double f = h * (*s).fStep;
		double w = TWOPI * (f - (*p).frequency);

		//supergaussian pulse spectrum, if no input pulse specified
		deviceComplex specfac = deviceComplex(-pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder), 0);

		deviceComplex specphase = deviceComplex(0,
			-((*p).cep
				+ TWOPI * f * ((*p).delay - 0.5 * (*s).dt * (*s).Ntime)
				+ 0.5 * (*p).gdd * w * w
				+ (*p).tod * w * w * w / 6.0
				+ materialPhase[h]));
		specfac = deviceLib::exp(specfac + specphase);
		
		if (hasLoadedField) {
			specfac = loadedField[h] * deviceLib::exp(specphase);
		}
		deviceComplex ne, no;
		sellmeierCuda(&ne, &no, sellmeierCoefficients, f, (*s).crystalTheta, (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);

		double ko = TWOPI * no.real() * f / LIGHTC;
		double zR = PI * (*p).beamwaist * (*p).beamwaist * ne.real() * f / LIGHTC;
		if (f == 0) {
			zR = 1e3;
		}
		double xo = ((*s).dx * (j - (*s).Nspace / 2.0)) - (*p).x0;
		double yo = ((*s).dx * (k - (*s).Nspace2 / 2.0)) - (*p).y0;
		if (!(*s).is3D) yo = 0.0;
		double zo = (*p).z0;
		double cB = cos((*p).beamAngle);
		double cA = cos((*p).beamAnglePhi);
		double sB = sin((*p).beamAngle);
		double sA = sin((*p).beamAnglePhi);
		double x = cB * xo + sA * sB * yo + sA * sB * zo;
		double y = cA * yo - sA * zo;
		double z = -sB * xo + sA * cB * yo + cA * cB * zo;
		double r = sqrt(x * x + y * y);

		double wz = (*p).beamwaist * sqrt(1 + (z * z / (zR * zR)));
		double Rz = 1.0e15;
		if (z != 0.0) {
			Rz = z * (1. + (zR * zR / (z * z)));
		}

		double phi = atan(z / zR);
		deviceComplex Eb = ((*p).beamwaist / wz) * deviceLib::exp(deviceComplex(0., 1.) * (ko * (z - (*p).z0) + ko * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
		Eb = Eb * specfac;
		if (isnan(cuCModSquared(Eb)) || f <= 0) {
			Eb = deviceComplex(0., 0.);
		}

		field[i] = deviceComplex(cos((*p).polarizationAngle), -(*p).circularity * sin((*p).polarizationAngle)) * Eb;
		field[i + (*s).NgridC] = deviceComplex(sin((*p).polarizationAngle), (*p).circularity * cos((*p).polarizationAngle)) * Eb;
		double pointEnergy = (cuCModSquared(field[i]) + cuCModSquared(field[i + (*s).NgridC]));
		pointEnergy *= 2 * LIGHTC * EPS0 * (*s).dx * (*s).dx * (*s).dt;

		//factor 2 accounts for the missing negative frequency plane
		atomicAddDevice(pulseSum, pointEnergy);
	};

	trilingual multiplyByConstantKernelD asKernel(withID double* A, double val) {
		long long i = localIndex;
		A[i] = val * A[i];
	};

	trilingual multiplyByConstantKernelDZ asKernel(withID deviceComplex* A, double val) {
		size_t i = localIndex;
		A[i] = val * A[i];

	};

	trilingual multiplicationKernelCompactVector asKernel(withID deviceComplex* A, deviceComplex* B, deviceComplex* C, deviceParameterSet* s) {
		long long i = localIndex;
		long long h = i % (*s).Nfreq; //temporal coordinate

		C[i] = A[h] * B[i];
	};

	trilingual multiplicationKernelCompactDoubleVector asKernel(withID double* A, deviceComplex* B, deviceComplex* C, deviceParameterSet* s) {
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
	typedef dlib::matrix<double, 0, 1> column_vector;
	simulationParameterSet* fittingSet;
	
	int getTotalSpectrum(activeDevice& d) {
		simulationParameterSet* sCPU = d.cParams;
		deviceParameterSet* sc = d.dParams;
		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.fft((*sc).gridETime1, (*sc).workspace1, deviceFFTD2Z1D);
		if ((*sc).is3D) {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum3DKernel, (*sc).workspace1, (*sc).workspace2, (*sCPU).rStep, (*sCPU).Ntime / 2 + 1, (*sCPU).Nspace * (*sCPU).Nspace2, (*sc).gridPolarizationTime1);
		}
		else {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel, (*sc).workspace1, (*sc).workspace2, (*sCPU).rStep, (*sCPU).Ntime / 2 + 1, (*sCPU).Nspace, (*sc).gridPolarizationTime1);
		}
		d.deviceMemcpy((*sCPU).totalSpectrum, (*sc).gridPolarizationTime1, 3 * (*sCPU).Nfreq * sizeof(double), DeviceToHost);
		return 0;
	}

	int addPulseToFieldArrays(activeDevice& d, pulse& pCPU, bool useLoadedField, std::complex<double>* loadedFieldIn) {

		simulationParameterSet* s = d.cParams;
		deviceParameterSet* sc = d.dParams;
		deviceParameterSet* scDevice = d.dParamsDevice;
		pulse* p;
		d.deviceCalloc((void**)&p, 1, sizeof(pulse));
		d.deviceMemcpy(d.dParamsDevice, sc, sizeof(deviceParameterSet), HostToDevice);

		double* materialPhase;
		deviceComplex* loadedField;

		d.deviceCalloc((void**)&loadedField, (*sc).Ntime, sizeof(deviceComplex));

		//get the material phase
		double* materialCoefficients, * sellmeierPropagationMedium;

		d.deviceCalloc((void**)&materialCoefficients, 66, sizeof(double));
		d.deviceCalloc((void**)&sellmeierPropagationMedium, 66, sizeof(double));
		d.deviceCalloc((void**)&materialPhase, (*s).Ntime, sizeof(double));

		d.deviceMemcpy(materialCoefficients, (*s).crystalDatabase[pCPU.phaseMaterial].sellmeierCoefficients, 66 * sizeof(double), HostToDevice);
		d.deviceMemcpy(sellmeierPropagationMedium, (*s).crystalDatabase[(*s).materialIndex].sellmeierCoefficients, 66 * sizeof(double), HostToDevice);
		d.deviceLaunch((unsigned int)(*s).Ntime, 1, materialPhaseKernel, (*s).fStep, (*s).Ntime, materialCoefficients, pCPU.frequency, pCPU.frequency, pCPU.phaseMaterialThickness, pCPU.phaseMaterialThickness, materialPhase, materialPhase);

		double* pulseSum = &materialCoefficients[0];

		if (useLoadedField) {
			d.deviceMemcpy(loadedField, loadedFieldIn, (*s).Ntime * sizeof(deviceComplex), HostToDevice);
		}

		d.deviceMemset(pulseSum, 0, sizeof(double));
		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.deviceMemcpy(p, &pCPU, sizeof(pulse), HostToDevice);
		if ((*sc).is3D) {
			d.deviceLaunch((*sc).Nblock / 2, (*sc).Nthread, beamGenerationKernel3D,
				(*sc).workspace1, p, pulseSum, scDevice, useLoadedField, loadedField, materialPhase,
				sellmeierPropagationMedium);
		}
		else {
			d.deviceLaunch((*sc).Nblock / 2, (*sc).Nthread, beamGenerationKernel2D,
				(*sc).workspace1, p, pulseSum, scDevice, useLoadedField, loadedField, materialPhase,
				sellmeierPropagationMedium);
		}

		d.fft((*sc).workspace1, (*sc).gridPolarizationTime1, deviceFFTZ2D1D);

		d.deviceLaunch(2 * (*sc).Nblock, (*sc).Nthread, beamNormalizeKernel, scDevice, pulseSum, (*sc).gridPolarizationTime1, pCPU.energy);

		//add the pulses
		d.deviceLaunch(2 * (*sc).Nblock, (*sc).Nthread, addDoubleArraysKernel, (*sc).gridETime1, (double*)(*sc).gridPolarizationTime1);

		//fft onto frequency grid
		d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFTD2Z);

		d.deviceFree(materialPhase);
		d.deviceFree(materialCoefficients);
		d.deviceFree(sellmeierPropagationMedium);
		d.deviceFree(loadedField);
		d.deviceFree(p);
		return 0;
	}
	
	int prepareElectricFieldArrays(activeDevice& d) {

		simulationParameterSet* s = d.cParams;
		deviceParameterSet* sc = d.dParams;
		d.deviceMemcpy(d.dParamsDevice, sc, sizeof(deviceParameterSet), HostToDevice);
		deviceParameterSet* scDevice = d.dParamsDevice;
		
		if (!(*s).isFollowerInSequence || (*s).isReinjecting) {
			if (!(*s).isReinjecting) {
				d.deviceMemset((*sc).gridETime1, 0, 2 * (*sc).Ngrid * sizeof(double));
			}
			else {
				d.deviceMemcpy((*sc).gridETime1, (*s).ExtOut, 2 * (*s).Ngrid * sizeof(double), HostToDevice);
				d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFTD2Z);
			}
			addPulseToFieldArrays(d, d.cParams->pulse1, d.cParams->field1IsAllocated, d.cParams->loadedField1);
			addPulseToFieldArrays(d, d.cParams->pulse2, d.cParams->field2IsAllocated, d.cParams->loadedField2);
		}
		else {
			d.deviceMemcpy((*sc).gridETime1, (*s).ExtOut, 2 * (*s).Ngrid * sizeof(double), HostToDevice);
			d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFTD2Z);
		}
		
		//Copy the field into the temporary array
		d.deviceMemcpy((*sc).gridEFrequency1Next1, (*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), DeviceToDevice);

		//set the propagation grids how they should be at the beginning of the next step
		d.deviceLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplicationKernelCompactDoubleVector, (*sc).fieldFactor1, (*sc).gridEFrequency1Next1, (*sc).workspace1, scDevice);
		d.deviceLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplyByConstantKernelDZ, (*sc).workspace1, (*sc).fftNorm);
		d.deviceLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplicationKernelCompact, (*sc).gridPropagationFactor1, (*sc).gridEFrequency1Next1, (*sc).k1);

		return 0;
	}

	int applyFresnelLoss(simulationParameterSet* s, int materialIndex1, int materialIndex2) {
		activeDevice d;
		deviceParameterSet sc;
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
		d.fft(sc.gridEFrequency1, sc.gridETime1, deviceFFTZ2D);
		d.deviceLaunch(2 * sc.Nblock, sc.Nthread, multiplyByConstantKernelD, sc.gridETime1, 1.0 / sc.Ngrid);
		//copy the field arrays from the GPU to CPU memory
		d.deviceMemcpy((*s).ExtOut, sc.gridETime1, 2 * (*s).Ngrid * sizeof(double), DeviceToHost);
		d.deviceMemcpy((*s).EkwOut, sc.gridEFrequency1, 2 * (*s).Ngrid * sizeof(deviceComplex), DeviceToHost);

		d.deviceFree(sellmeierCoefficients1);
		d.deviceFree(sellmeierCoefficients2);
		d.deallocateSet(&sc);
		return 0;
	}

	int applyFilter(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, double f0, double bandwidth, double order, double inBandAmplitude, double outOfBandAmplitude) {

		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z);
		deviceParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet), HostToDevice);
		d.deviceLaunch(s.Nblock / 2, s.Nthread, filterKernel, sDevice, 1.0e12*f0, 1.0e12*bandwidth, order, inBandAmplitude, outOfBandAmplitude);

		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);

		d.fft(s.gridEFrequency1, s.gridETime1, deviceFFTZ2D);

		d.deviceLaunch((int)(s.Ngrid / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);

		getTotalSpectrum(d);

		return 0;
	}

	int applyAperatureFarField(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, double diameter, double activationParameter, double xOffset, double yOffset) {
		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z);
		deviceParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet), HostToDevice);
		d.deviceLaunch(s.Nblock/2, s.Nthread, apertureFarFieldKernel, sDevice, 0.5 * DEG2RAD * diameter, activationParameter, DEG2RAD * xOffset, DEG2RAD * yOffset);

		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);


		d.fft(s.gridEFrequency1, s.gridETime1, deviceFFTZ2D);

		d.deviceLaunch((int)(s.Ngrid / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);

		getTotalSpectrum(d);

		return 0;
	}


	int applyAperature(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s,double diameter, double activationParameter) {
		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);

		deviceParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet), HostToDevice);
		d.deviceLaunch(s.Nblock, s.Nthread, apertureKernel, sDevice, 0.5 * diameter, activationParameter);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), DeviceToHost);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);
		getTotalSpectrum(d);
		return 0;
	}

	int applySphericalMirror(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, double ROC) {

		deviceParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet), HostToDevice);

		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z1D);
		d.deviceLaunch(s.Nblock / 2, s.Nthread, sphericalMirrorKernel, sDevice, ROC);
		d.fft(s.gridEFrequency1, s.gridETime1, deviceFFTZ2D1D);
		d.deviceLaunch(2 * s.Nblock, s.Nthread, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ntime);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), DeviceToHost);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);
		getTotalSpectrum(d);
		return 0;
	}

	int applyParabolicMirror(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, double focus) {
		deviceParameterSet* sDevice = d.dParamsDevice;

		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z1D);
		d.deviceLaunch(s.Nblock / 2, s.Nthread, parabolicMirrorKernel, sDevice, focus);
		d.fft(s.gridEFrequency1, s.gridETime1, deviceFFTZ2D1D);
		d.deviceLaunch(2 * s.Nblock, s.Nthread, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ntime);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), DeviceToHost);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);
		getTotalSpectrum(d);
		return 0;
	}

	int applyLinearPropagation(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, int materialIndex, double thickness) {

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
		deviceParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet), HostToDevice);

		d.deviceLaunch(s.Nblock / 2, s.Nthread, applyLinearPropagationKernel, sellmeierCoefficients, thickness, sDevice);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, s.NgridC * 2 * sizeof(deviceComplex), DeviceToHost);
		d.fft(s.gridEFrequency1, s.gridETime1, deviceFFTZ2D);
		d.deviceLaunch(2 * s.Nblock, s.Nthread, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);

		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), DeviceToHost);

		return 0;
	}

	int preparePropagationGrids(activeDevice& d) {
		deviceParameterSet* sc = d.dParams;
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
		deviceParameterSet* sD = d.dParamsDevice;
		d.deviceMemcpy(sD, sc, sizeof(deviceParameterSet), HostToDevice);
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
		d.deviceMemcpy(sc, sD, sizeof(deviceParameterSet), DeviceToHost);
		return 0;
	}

	//Rotate the field on the GPU
	//Allocates memory and copies from CPU, then copies back to CPU and deallocates
	// - inefficient but the general principle is that only the CPU memory is preserved
	// after simulations finish... and this only runs at the end of the simulation
	int rotateField(simulationParameterSet* s, double rotationAngle) {
		deviceParameterSet sc;
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
		d.fft(Eout1, sc.gridETime1, deviceFFTZ2D);
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
		deviceParameterSet* sH = d.dParams; 
		deviceParameterSet* sD = d.dParamsDevice;
		//operations involving FFT
		
		if ((*sH).isNonLinear || (*sH).isCylindric) {
			//perform inverse FFT to get time-space electric field
			d.fft((*sH).workspace1, (*sH).gridETime1, deviceFFTZ2D);

			//Nonlinear polarization
			if ((*sH).isNonLinear) {
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, nonlinearPolarizationKernel, sD);
				if ((*sH).isCylindric) {
					d.deviceLaunch((*sH).Nblock, (*sH).Nthread, expandCylindricalBeam, sD);
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFTD2ZPolarization);
					d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPolarizationKernelCylindric, sD);
				}
				else {
					d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFTD2Z);
					d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPolarizationKernel, sD);
				}
			}

			//Plasma/multiphoton absorption
			if ((*sH).hasPlasma) {
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, plasmaCurrentKernel_twoStage_A, sD);
				d.deviceLaunch((unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / MIN_GRIDDIM), 2*MIN_GRIDDIM, plasmaCurrentKernel_twoStage_B, sD);
				if ((*sH).isCylindric) {
					d.deviceLaunch((*sH).Nblock, (*sH).Nthread, expandCylindricalBeam, sD);
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFTD2ZPolarization);
					d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPlasmaKernelCylindric, sD);
				}
				else {
					d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFTD2Z);
					d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPlasmaKernel, sD);
				}
			}

			//Radial Laplacian
			if ((*sH).isCylindric) {
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, radialLaplacianKernel, sD);
				d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFTD2Z);
			}
		}

		//advance an RK4 step
		switch (stepNumber) {
		case 0:
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, rkKernel0, sD);
			break;
		case 1:
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, rkKernel1, sD);
			break;
		case 2:
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, rkKernel2, sD);
			break;
		case 3:
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, rkKernel3, sD);
		}
		return 0;
	}



	unsigned long int solveNonlinearWaveEquationWithDevice(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s) {
		//prepare the propagation arrays
		preparePropagationGrids(d);
		prepareElectricFieldArrays(d);
		d.deviceMemcpy(d.dParamsDevice, &s, sizeof(deviceParameterSet), HostToDevice);
		double* canaryPointer = &s.gridETime1[s.Ntime / 2 + s.Ntime * (s.Nspace / 2 + s.Nspace * (s.Nspace2 / 2))];

		//Core propagation loop
		for (size_t i = 0; i < s.Nsteps; ++i) {

			//RK4
			runRK4Step(d, 0);
			runRK4Step(d, 1);
			runRK4Step(d, 2);
			runRK4Step(d, 3);

			//periodically check if the simulation diverged or was cancelled
			if ((*sCPU).statusFlags[0] == 2) break;
			if (i % 10 == 0) if (d.isTheCanaryPixelNaN(canaryPointer)) break;

			//copy the field arrays from the GPU to CPU memory if requested by the UI
			if ((*sCPU).statusFlags[0] == 3) {
				d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);
				d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * (*sCPU).NgridC * sizeof(deviceComplex), DeviceToHost);
				(*sCPU).statusFlags[0] = 0;
			}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
		}
		if ((*sCPU).isInFittingMode && !(*sCPU).isInSequence)(*(*sCPU).progressCounter)++;

		//take final spectra and transfer the results to the CPU
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);
		d.fft(s.gridEFrequency1, s.gridETime1, deviceFFTZ2D);
		d.deviceLaunch((int)(s.Ngrid / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);
		getTotalSpectrum(d);

		int returnval = 13 * d.isTheCanaryPixelNaN(canaryPointer);

		(*sCPU).statusFlags[0] = 1;
		return returnval;
	}

	constexpr unsigned int funHash(const char* s, int off = 0) {
		return (s[off] == 0 || s[off] == '(') ? 7177 : (funHash(s, off + 1) * 31) ^ s[off];
	}

	unsigned int stringHash(std::string& s, int off = 0){
		return (s.length() == off || s.at(off) == '(') ? 7177 : (stringHash(s,off+1) * 31) ^ s.at(off);
	}
	//Dispatcher of the sequence mode. New functions go here, and should have a unique hash (chances of a conflict are small, and 
	// will produce a compile-time error.
	// Functions cannot start with a number or the string "None".
	int interpretCommand(std::string cc, double* iBlock, double* vBlock, activeDevice& d, simulationParameterSet *sCPU, deviceParameterSet& s) {
		crystalEntry* db = (*sCPU).crystalDatabase;
		int error = 0;
		double parameters[32] = {0.0};
		bool defaultMask[32] = {0};

		switch (stringHash(cc)) {
		case funHash("rotate"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			rotateField(sCPU, DEG2RAD * parameters[0]);
			break;
		case funHash("set"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			vBlock[(int)parameters[0]] = parameters[1];
			break;
		case funHash("plasmaReinject"):
			(*sCPU).isReinjecting = TRUE;
		case funHash("plasma"):
			interpretParameters(cc, 9, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = DEG2RAD * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = DEG2RAD * parameters[2];
			if (!defaultMask[3])(*sCPU).nonlinearAbsorptionStrength = parameters[3];
			if (!defaultMask[4])(*sCPU).bandGapElectronVolts = parameters[4];
			if (!defaultMask[5])(*sCPU).drudeGamma = parameters[5];
			if (!defaultMask[6])(*sCPU).effectiveMass = parameters[6];
			if (!defaultMask[7])(*sCPU).crystalThickness = 1e-6 * parameters[7];
			if (!defaultMask[8])(*sCPU).propagationStep = 1e-9 * parameters[8];
			(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d;
			(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3;
			(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches;
			(*sCPU).absorptionParameters = db[(*sCPU).materialIndex].absorptionParameters;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients;

			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;

			error = solveNonlinearWaveEquationWithDevice(d, sCPU, s);
			(*sCPU).isFollowerInSequence = TRUE;
			break;
		case funHash("nonlinear"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = DEG2RAD * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = DEG2RAD * parameters[2];
			if (!defaultMask[3])(*sCPU).crystalThickness = 1e-6 * parameters[3];
			if (!defaultMask[4])(*sCPU).propagationStep = 1e-9 * parameters[4];

			(*sCPU).nonlinearAbsorptionStrength = 0.0;
			(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d;
			(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3;
			(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches;
			(*sCPU).absorptionParameters = db[(*sCPU).materialIndex].absorptionParameters;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients;

			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU, &s);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU, s);
			(*sCPU).isFollowerInSequence = TRUE;
			break;
		case funHash("default"):
			d.reset(sCPU, &s);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU, s);
			break;
		case funHash("init"):
			(*sCPU).materialIndex = 0;
			(*sCPU).crystalTheta = 0.0;
			(*sCPU).crystalPhi = 0.0;
			(*sCPU).crystalThickness = 0;
			(*sCPU).propagationStep = 1e-9;

			(*sCPU).nonlinearAbsorptionStrength = 0.0;
			(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d;
			(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3;
			(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches;
			(*sCPU).absorptionParameters = db[(*sCPU).materialIndex].absorptionParameters;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients;

			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU, &s);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU, s);
			(*sCPU).isFollowerInSequence = TRUE;
			break;

		case funHash("linear"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if ((*sCPU).isCylindric) {
				if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
				if (!defaultMask[1])(*sCPU).crystalTheta = DEG2RAD * parameters[1];
				if (!defaultMask[2])(*sCPU).crystalPhi = DEG2RAD * parameters[2];
				if (!defaultMask[3])(*sCPU).crystalThickness = 1e-6 * parameters[3];
				if (!defaultMask[4])(*sCPU).propagationStep = 1e-9 * parameters[4];
				(*sCPU).nonlinearAbsorptionStrength = 0.0;
				(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d;
				(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3;
				(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches;
				(*sCPU).absorptionParameters = db[(*sCPU).materialIndex].absorptionParameters;
				(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients;
				(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
				(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
				(*sCPU).forceLinear = TRUE;
				d.reset(sCPU, &s);
				error = solveNonlinearWaveEquationWithDevice(d, sCPU, s);
				(*sCPU).isFollowerInSequence = TRUE;
			}
			else {
				if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
				if (!defaultMask[1])(*sCPU).crystalTheta = DEG2RAD * parameters[1];
				if (!defaultMask[2])(*sCPU).crystalPhi = DEG2RAD * parameters[2];
				if (!defaultMask[3])(*sCPU).crystalThickness = 1e-6 * parameters[3];

				applyLinearPropagation(d, sCPU, s, (*sCPU).materialIndex, (*sCPU).crystalThickness);
			}

			break;
		case funHash("fresnelLoss"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = DEG2RAD * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = DEG2RAD * parameters[2];
			d.reset(sCPU, &s);
			applyFresnelLoss(sCPU,
				(int)parameters[4],
				(int)parameters[5]);
			break;
		case funHash("sphericalMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applySphericalMirror(d, sCPU, s, parameters[0]);
			break;
		case funHash("parabolicMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applyParabolicMirror(d, sCPU, s, parameters[0]);
			break;
		case funHash("aperture"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applyAperature(d, sCPU, s,
				parameters[0],
				parameters[1]);
			break;
		case funHash("farFieldAperture"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applyAperatureFarField(d, sCPU, s,
				parameters[0],
				parameters[1],
				parameters[2],
				parameters[3]);
			break;
		case funHash("energy"):
			{
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			int targetVar = (int)parameters[0];
			int spectrumType = (int)parameters[1];
			double energy = 0.0;
			for(int i = 0; i < (*sCPU).Nfreq; i++){
				energy += (*sCPU).totalSpectrum[i + (*sCPU).Nfreq * spectrumType];
			}
			vBlock[targetVar] = energy;
			}
			break;
		case funHash("filter"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applyFilter(d, sCPU, s,
				parameters[0],
				parameters[1],
				parameters[2],
				parameters[3],
				parameters[4]);
			break;
		case funHash("addPulse"):
		{
			interpretParameters(cc, 20, iBlock, vBlock, parameters, defaultMask);
			activeDevice d;
			d.dParams = &s;
			d.cParams = sCPU;
			d.reset(sCPU, &s);

			d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
			d.deviceMemcpy(s.gridEFrequency1, (*sCPU).EkwOut, 2 * s.NgridC * sizeof(deviceComplex), HostToDevice);

			pulse p;
			memcpy(&p, &(sCPU->pulse1), sizeof(pulse));
			p.energy = parameters[0];
			p.frequency = 1e12 * parameters[1];
			p.bandwidth = 1e12 * parameters[2];
			p.sgOrder = (int)parameters[3];
			p.cep = parameters[4] / PI;
			p.delay = 1e-15 * parameters[5];
			p.gdd = 1e-30 * parameters[6];
			p.tod = 1e-45 * parameters[7];
			p.phaseMaterial = (int)parameters[8];
			p.phaseMaterialThickness = 1e-6 * parameters[9];
			p.beamwaist = 1e-6 * parameters[10];
			p.x0 = 1e-6 * parameters[11];
			p.z0 = 1e-6 *parameters[12];
			p.beamAngle = DEG2RAD * parameters[13];
			p.beamAnglePhi = DEG2RAD * parameters[14];
			p.polarizationAngle = DEG2RAD * parameters[15];
			p.circularity = parameters[16];
			(*sCPU).materialIndex = (int)parameters[17];
			(*sCPU).crystalTheta = DEG2RAD * parameters[18];
			(*sCPU).crystalPhi = DEG2RAD * parameters[19];

			addPulseToFieldArrays(d, p, FALSE, NULL);
			d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(deviceComplex), DeviceToHost);
			d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);

			getTotalSpectrum(d);
		}
			break;
		case funHash("for"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			int counter = (int)parameters[0];
			int targetVar = (int)parameters[1];
			std::string currentString = cc.substr(cc.find_first_of('{')+1,std::string::npos);

			std::string forStartString = currentString;
			for (int i = 0; i < counter; i++) {
				while (currentString.at(0) != '}') {
					interpretCommand(currentString, iBlock, vBlock, d, sCPU, s);
					currentString = currentString.substr(currentString.find_first_of(')')+1,std::string::npos);
				}
				++vBlock[targetVar];
				currentString = forStartString;
			}
			break;
		}
		return error;
	}

	//deprecated sequence mode; kept for compatibility with old files, but better to use interpretCommand() instead
	int resolveSequence(int currentIndex, simulationParameterSet * s, crystalEntry * db) {


		//double* offsetArray = &(*s).sequenceArray[11 * currentIndex];
		//int error = 0;
		////sequence format
		////0: step type
		//int stepType = (int)offsetArray[0];
		//int materialIndex = 0;
		//double thickness = 0;
		//// 
		//// if stepType == 0, normal propagation
		////1: material index
		////2: theta,
		////3: phi, 
		////4: NL absorption
		////5: Band gap
		////6: Drude relaxation
		////7: Effective mass
		////8: Crystal thickness
		////9: Propagation step size
		////10: rotation angle
		////
		//// if stepType == 1, linear propagation
		//// same parameters as 0, but only 1,2,3,8, and 10 matter
		////
		//// if stepType == 2, fresnel loss
		//// 1: incidence material index
		//// 2: transmission material index
		//// other parameters don't matter
		//// 
		//// if stepType == 3, spherical mirror
		//// 1: ROC (m)
		////
		//// if stepType == 4, parabolic mirror
		//// 1: focus (m)
		//// 
		//// if stepType == 5, aperture
		//// 1: diameter (m)
		//// 2: activation parameter p (function is 1 - 1/(1 + exp(-p*(r-radius)))
		////
		//// if stepType == 6, loop back to start (handled by solveNonlinearWaveEquationSequence())
		//// 1: counter (counts down to zero)
		////
		//// if stepType == 7, reinjection, same as 0, but input fields are added to current fields.

		//switch (stepType) {
		//case 7:
		//	(*s).isReinjecting = TRUE;
		//case 0:
		//	if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
		//	if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
		//	if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
		//	if ((int)offsetArray[4] != -1) (*s).nonlinearAbsorptionStrength = offsetArray[4];
		//	if ((int)offsetArray[5] != -1) (*s).bandGapElectronVolts = offsetArray[5];
		//	if ((int)offsetArray[6] != -1) (*s).drudeGamma = offsetArray[6];
		//	if ((int)offsetArray[7] != -1) (*s).effectiveMass = offsetArray[7];
		//	if ((int)offsetArray[8] != -1) (*s).crystalThickness = 1e-6 * offsetArray[8];
		//	if ((int)offsetArray[9] != -1) (*s).propagationStep = 1e-9 * offsetArray[9];
		//	if ((int)offsetArray[8] != -1) (*s).Npropagation
		//		= (size_t)(1e-6 * offsetArray[8] / (*s).propagationStep);
		//	if (currentIndex > 0) {
		//		(*s).isFollowerInSequence = TRUE;
		//	}
		//	(*s).chi2Tensor = db[(*s).materialIndex].d;
		//	(*s).chi3Tensor = db[(*s).materialIndex].chi3;
		//	(*s).nonlinearSwitches = db[(*s).materialIndex].nonlinearSwitches;
		//	(*s).absorptionParameters = db[(*s).materialIndex].absorptionParameters;
		//	(*s).sellmeierCoefficients = db[(*s).materialIndex].sellmeierCoefficients;

		//	(*s).sellmeierType = db[(*s).materialIndex].sellmeierType;
		//	(*s).axesNumber = db[(*s).materialIndex].axisType;

		//	error = solveNonlinearWaveEquationX(s);

		//	if (offsetArray[10] != 0.0) {
		//		rotateField(s, DEG2RAD * offsetArray[10]);
		//	}

		//	if ((*s).memoryError > 0) {
		//		printf("Warning: device memory error (%i).\n", (*s).memoryError);
		//	}
		//	return error;

		//case 1:
		//	if ((*s).isCylindric) {
		//		if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
		//		if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
		//		if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
		//		if ((int)offsetArray[4] != -1) (*s).nonlinearAbsorptionStrength = offsetArray[4];
		//		if ((int)offsetArray[5] != -1) (*s).bandGapElectronVolts = offsetArray[5];
		//		if ((int)offsetArray[6] != -1) (*s).drudeGamma = offsetArray[6];
		//		if ((int)offsetArray[7] != -1) (*s).effectiveMass = offsetArray[7];
		//		if ((int)offsetArray[8] != -1) (*s).crystalThickness = 1e-6 * offsetArray[8];
		//		if ((int)offsetArray[9] != -1) (*s).propagationStep = 1e-9 * offsetArray[9];
		//		if ((int)offsetArray[8] != -1) (*s).Npropagation
		//			= (size_t)(1e-6 * offsetArray[8] / (*s).propagationStep);
		//		if (currentIndex > 0) {
		//			(*s).isFollowerInSequence = TRUE;
		//		}
		//		(*s).chi2Tensor = db[(*s).materialIndex].d;
		//		(*s).chi3Tensor = db[(*s).materialIndex].chi3;
		//		(*s).nonlinearSwitches = db[(*s).materialIndex].nonlinearSwitches;
		//		(*s).absorptionParameters = db[(*s).materialIndex].absorptionParameters;
		//		(*s).sellmeierCoefficients = db[(*s).materialIndex].sellmeierCoefficients;
		//		(*s).sellmeierType = db[(*s).materialIndex].sellmeierType;
		//		(*s).axesNumber = db[(*s).materialIndex].axisType;
		//		(*s).forceLinear = TRUE;
		//		solveNonlinearWaveEquationX(s);
		//	}
		//	else {
		//		if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
		//		if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
		//		if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
		//		thickness = 1.0e-6 * offsetArray[8];
		//		if (offsetArray[8] == -1) {
		//			thickness = (*s).crystalThickness;
		//		}
		//		materialIndex = (int)offsetArray[1];
		//		if (offsetArray[1] == -1) {
		//			materialIndex = (*s).materialIndex;
		//		}
		//		applyLinearPropagation(s, materialIndex, thickness);
		//	}

		//	if (offsetArray[10] != 0.0) {
		//		rotateField(s, DEG2RAD * offsetArray[10]);
		//	}
		//	return 0;

		//case 2:
		//	if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
		//	if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
		//	if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
		//	applyFresnelLoss(s, (int)offsetArray[4], (int)offsetArray[5]);
		//	return 0;
		//case 3:
		//	applySphericalMirror(s, offsetArray[8]);
		//	if (offsetArray[10] != 0.0) {
		//		rotateField(s, DEG2RAD * offsetArray[10]);
		//	}
		//	return 0;
		//case 4:
		//	applyParabolicMirror(s, offsetArray[8]);
		//	if (offsetArray[10] != 0.0) {
		//		rotateField(s, DEG2RAD * offsetArray[10]);
		//	}
		//	return 0;
		//case 5:
		//	applyAperature(s, offsetArray[1], offsetArray[2]);
		//	if (offsetArray[10] != 0.0) {
		//		rotateField(s, DEG2RAD * offsetArray[10]);
		//	}
		//	return 0;
		//}

		return 1;
	}

	// helper function for fitting mode, runs the simulation and returns difference from the desired outcome.
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
		&(*fittingSet).pulse1.energy, &(*fittingSet).pulse2.energy, &(*fittingSet).pulse1.frequency, &(*fittingSet).pulse2.frequency,
		&(*fittingSet).pulse1.bandwidth, &(*fittingSet).pulse2.bandwidth, &(*fittingSet).pulse1.cep, &(*fittingSet).pulse2.cep,
		&(*fittingSet).pulse1.delay, &(*fittingSet).pulse2.delay, &(*fittingSet).pulse1.gdd, &(*fittingSet).pulse2.gdd,
		&(*fittingSet).pulse1.tod, &(*fittingSet).pulse2.tod, &(*fittingSet).pulse1.phaseMaterialThickness, &(*fittingSet).pulse2.phaseMaterialThickness,
		&(*fittingSet).pulse1.beamwaist, &(*fittingSet).pulse2.beamwaist,
		&(*fittingSet).pulse1.x0, &(*fittingSet).pulse2.x0, &(*fittingSet).pulse1.z0, &(*fittingSet).pulse2.z0,
		&(*fittingSet).pulse1.beamAngle, &(*fittingSet).pulse2.beamAngle, &(*fittingSet).pulse1.polarizationAngle, &(*fittingSet).pulse2.polarizationAngle,
		&(*fittingSet).pulse1.circularity, &(*fittingSet).pulse2.circularity, &(*fittingSet).crystalTheta, &(*fittingSet).crystalPhi,
		&(*fittingSet).nonlinearAbsorptionStrength, &(*fittingSet).drudeGamma, &(*fittingSet).effectiveMass, &(*fittingSet).crystalThickness,
		&(*fittingSet).propagationStep };

		for (int i = 0; i < (*fittingSet).Nfitting; ++i) {
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
			for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
				result += (*fittingSet).totalSpectrum[(*fittingSet).fittingMode * (*fittingSet).Nfreq + (*fittingSet).fittingROIstart + i];
			}
			return result;
		}

		//mode 3: match total spectrum to reference given in ascii file
		double a;
		double maxSim = 0.0;
		double maxRef = 0.0;
		double* simSpec = &(*fittingSet).totalSpectrum[2 * (*fittingSet).Nfreq + (*fittingSet).fittingROIstart];
		double* refSpec = &(*fittingSet).fittingReference[(*fittingSet).fittingROIstart];
		for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
			maxSim = maxN(maxSim, simSpec[i]);
			maxRef = maxN(maxRef, refSpec[i]);
		}

		if (maxSim == 0) {
			maxSim = 1;
		}
		if (maxRef == 0) {
			maxRef = 1;
		}
		result = 0.0;
		for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
			a = (refSpec[i] / maxRef) - (simSpec[i] / maxSim);
			result += a * a;
		}
		return sqrt(result);
	}

	//the old sequence mode, can be called by the new sequence mode if it encounters
	// something that looks like an old sequence.
	unsigned long solveNonlinearWaveEquationSequenceOldX(void* lpParam) {
		simulationParameterSet sCPUcurrent;
		simulationParameterSet* sCPU = &sCPUcurrent;//(simulationParameterSet*)lpParam;
		memcpy(sCPU, (simulationParameterSet*)lpParam, sizeof(simulationParameterSet));
		simulationParameterSet sCPUbackupValues;
		simulationParameterSet* sCPUbackup = &sCPUbackupValues;
		memcpy(sCPUbackup, sCPU, sizeof(simulationParameterSet));
		int k;
		int error = 0;
		for (k = 0; k < (*sCPU).Nsequence; ++k) {
			if ((int)round((*sCPU).sequenceArray[k * 11]) == 6) {
				if (((int)round((*sCPU).sequenceArray[k * 11 + 1])) > 0) {
					(*sCPUbackup).sequenceArray[k * 11 + 1] -= 1.0;
					(*sCPUbackup).isFollowerInSequence = TRUE;
					k = (int)round((*sCPU).sequenceArray[k * 11 + 2]);
					error = resolveSequence(k, sCPU, (*sCPU).crystalDatabase);
				}
			}
			else {
				error = resolveSequence(k, sCPU, (*sCPU).crystalDatabase);
			}

			if (error) break;
			memcpy(sCPU, sCPUbackup, sizeof(simulationParameterSet));
		}
		if ((*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
		return error;
	}

	
}
using namespace hostFunctions;

//Main (non sequence) solver. New device typeps should have a unique definition of the name
// e.g. solveNonlinearWaveEquationSYCL so that the correct one may be called. That's why it's
// a preprocessor definition here.
unsigned long solveNonlinearWaveEquationX(void* lpParam) {
	simulationParameterSet* sCPU = (simulationParameterSet*)lpParam;
	deviceParameterSet s;
	memset(&s, 0, sizeof(deviceParameterSet));
	activeDevice d(sCPU, &s);
	if (d.memoryStatus) return 1;
	unsigned long returnValue = solveNonlinearWaveEquationWithDevice(d, sCPU, s);
	d.deallocateSet(&s);
	return returnValue;
}


// Main function for running a sequence
unsigned long solveNonlinearWaveEquationSequenceX(void* lpParam) {
	simulationParameterSet sCPUcurrent;
	simulationParameterSet* sCPU = &sCPUcurrent;//(simulationParameterSet*)lpParam;
	memcpy(sCPU, (simulationParameterSet*)lpParam, sizeof(simulationParameterSet));
	deviceParameterSet s;
	memset(&s, 0, sizeof(deviceParameterSet));
	activeDevice d(sCPU, &s);

	//pointers to where the various parameters are in the struct
	double* targets[36] = { 0,
		&(*sCPU).pulse1.energy, &(*sCPU).pulse2.energy, &(*sCPU).pulse1.frequency, &(*sCPU).pulse2.frequency,
		&(*sCPU).pulse1.bandwidth, &(*sCPU).pulse2.bandwidth, &(*sCPU).pulse1.cep, &(*sCPU).pulse2.cep,
		&(*sCPU).pulse1.delay, &(*sCPU).pulse2.delay, &(*sCPU).pulse1.gdd, &(*sCPU).pulse2.gdd,
		&(*sCPU).pulse1.tod, &(*sCPU).pulse2.tod, &(*sCPU).pulse1.phaseMaterialThickness, &(*sCPU).pulse2.phaseMaterialThickness,
		&(*sCPU).pulse1.beamwaist, &(*sCPU).pulse2.beamwaist,
		&(*sCPU).pulse1.x0, &(*sCPU).pulse2.x0, &(*sCPU).pulse1.z0, &(*sCPU).pulse2.z0,
		&(*sCPU).pulse1.beamAngle, &(*sCPU).pulse2.beamAngle, &(*sCPU).pulse1.polarizationAngle, &(*sCPU).pulse2.polarizationAngle,
		&(*sCPU).pulse1.circularity, &(*sCPU).pulse2.circularity, &(*sCPU).crystalTheta, &(*sCPU).crystalPhi,
		&(*sCPU).nonlinearAbsorptionStrength, &(*sCPU).drudeGamma, &(*sCPU).effectiveMass, &(*sCPU).crystalThickness,
		&(*sCPU).propagationStep };

	//unit multipliers from interface units to SI base units.
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

	//if it starts with 0, it's an old sequence; send it there
	if ((*sCPU).sequenceString[0] == '0') {
		//readSequenceString((simulationParameterSet*)lpParam);
		//solveNonlinearWaveEquationSequenceOldX(lpParam);
		return 0;
	}

	//main text interpreter
	simulationParameterSet sCPUbackupValues;
	simulationParameterSet* sCPUbackup = &sCPUbackupValues;
	memcpy(sCPUbackup, sCPU, sizeof(simulationParameterSet));

	double iBlock[100] = { 0.0 };

	for (int k = 1; k < 36; k++) {
		iBlock[k] = *(targets[k])/multipliers[k];
	}

	double vBlock[100] = { 0.0 };
	std::string currentString((*sCPU).sequenceString);

	//shortest command is either for() or init(), if there's only 4 characters left, it can only
	//be whitespace or other trailing symbols
	size_t minLength = 5;
	for (;;) {
		//skip curly braces (for loops should have been handled by interpretCommand() already)
		if (currentString.at(0) == '{') {
			currentString = currentString.substr(currentString.find_first_of('}'),std::string::npos);
			if(currentString.length()<minLength) break; 
			currentString = currentString.substr(1,std::string::npos);
		}
		//skip angle brackets (comments)
		if (currentString.at(0) == '<') {
			currentString = currentString.substr(currentString.find_first_of('>'),std::string::npos);
			if(currentString.length()<minLength) break; 
			currentString = currentString.substr(1,std::string::npos);
		}

		interpretCommand(currentString, iBlock, vBlock, d, sCPU, s);
		currentString = currentString.substr(currentString.find_first_of(')'),std::string::npos);

		if(currentString.length()<minLength) break; 

		currentString = currentString.substr(1,std::string::npos);
		
		(*sCPUbackup).isFollowerInSequence = (*sCPU).isFollowerInSequence;
		memcpy(sCPU, sCPUbackup, sizeof(simulationParameterSet));
	}

	d.deallocateSet(&s);
	return 0;
}

//run in fitting mode
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
	&(*sCPU).pulse1.energy, &(*sCPU).pulse2.energy, &(*sCPU).pulse1.frequency, &(*sCPU).pulse2.frequency,
	&(*sCPU).pulse1.bandwidth, &(*sCPU).pulse2.bandwidth, &(*sCPU).pulse1.cep, &(*sCPU).pulse2.cep,
	&(*sCPU).pulse1.delay, &(*sCPU).pulse2.delay, &(*sCPU).pulse1.gdd, &(*sCPU).pulse2.gdd,
	&(*sCPU).pulse1.tod, &(*sCPU).pulse2.tod, &(*sCPU).pulse1.phaseMaterialThickness, &(*sCPU).pulse2.phaseMaterialThickness,
	&(*sCPU).pulse1.beamwaist, &(*sCPU).pulse2.beamwaist,
	&(*sCPU).pulse1.x0, &(*sCPU).pulse2.x0, &(*sCPU).pulse1.z0, &(*sCPU).pulse2.z0,
	&(*sCPU).pulse1.beamAngle, &(*sCPU).pulse2.beamAngle, &(*sCPU).pulse1.polarizationAngle, &(*sCPU).pulse2.polarizationAngle,
	&(*sCPU).pulse1.circularity, &(*sCPU).pulse2.circularity, &(*sCPU).crystalTheta, &(*sCPU).crystalPhi,
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

	for (int i = 0; i < (*sCPU).Nfitting; ++i) {
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

	for (int i = 0; i < (*sCPU).Nfitting; ++i) {
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

//main function - if included in the GUI, this should have a different name
// than main() - this one applies when running in command line mode (e.g. on
// the clusters.)
int mainX(int argc, mainArgumentX){
	resolveArgv;
	int i, j;

	size_t progressCounter = 0;
	int CUDAdeviceCount = 1;
	if (hardwareCheck(&CUDAdeviceCount)) return 1;
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
	for (j = 0; j < (*crystalDatabasePtr).numberOfEntries; ++j) {
		printf("Material %i name: %s\n", j, crystalDatabasePtr[j].crystalNameW);
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

	if (((*sCPU).sequenceString[0] != 'N') && (*sCPU).sequenceString[0] != 0) (*sCPU).isInSequence = TRUE;
	configureBatchMode(sCPU);
	readFittingString(sCPU);
	auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

	if ((*sCPU).Nfitting != 0) {
		printf("Running optimization for %i iterations...\n", (*sCPU).fittingMaxIterations);

		runDlibFittingX(sCPU);

		auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
		printf("Finished after %8.4lf s. \n",
			1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
		saveDataSet(sCPU);

		printf("Optimization result:\n (index, value)\n");
		for (int i = 0; i < (*sCPU).Nfitting; ++i) {
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
	for (j = 0; j < (*sCPU).Nsims * (*sCPU).Nsims2; ++j) {

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
	for (i = 0; i < (*sCPU).Nsims * (*sCPU).Nsims2; ++i) {
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

	saveDataSet(sCPU);
	delete[] threadBlock;
	deallocateGrids(sCPU, TRUE);
	delete[] sCPU;
	delete[] crystalDatabasePtr;
	return 0;
}