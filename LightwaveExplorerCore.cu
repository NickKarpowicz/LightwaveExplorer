#include "LightwaveExplorerCore.cuh"
#include "LightwaveExplorerDevices/LightwaveExplorerUtilities.h"
#include "LightwaveExplorerDevices/LightwaveExplorerTrilingual.h"
#include <stdlib.h>
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

namespace deviceFunctions {
	//Expand the information contained in the radially-symmetric beam in the offset grid
	// representation.
	// see the expandCylindricalBeam() kernel for more details
	deviceFunction void expandCylindricalBeamDevice(const deviceParameterSet* s, long long i, deviceFP* expandedBeam1, deviceFP* sourceBeam1, deviceFP* sourceBeam2) {
		long long j = i / (*s).Ntime; //spatial coordinate
		long long k = i % (*s).Ntime; //temporal coordinate

		//positions on the expanded grid corresponding the the current index
		long long pos1 = 2 * ((*s).Nspace - j - 1) * (*s).Ntime + k;
		long long pos2 = (2 * j + 1) * (*s).Ntime + k;
		deviceFP* expandedBeam2 = expandedBeam1 + 2 * (*s).Ngrid;
		expandedBeam1[pos1] = sourceBeam1[i];
		expandedBeam1[pos2] = sourceBeam1[i];
		expandedBeam2[pos1] = sourceBeam2[i];
		expandedBeam2[pos2] = sourceBeam2[i];
	}

	//give the Dawson function value at x
	//used for Gaussian-based "Sellmeier" equation, as it is the Hilbert transform of the Gaussian function (susceptibility obeys Kramers-Kronig relations)
	//based on Rybicki G.B., Computers in Physics, 3,85-87 (1989)
	//this is the simplest implementation of the formula he provides, he also suggests speed improvements in case
	//evaluation of this becomes a bottleneck
	deviceFunction deviceFP deviceDawson(const deviceFP& x) {
		//parameters determining accuracy (higher n, smaller h -> more accurate but slower)
		int n = 15;
		deviceFP h = 0.3;

		//series expansion for small x
		if (deviceFPLib::abs(x) < 0.2) {
			deviceFP x2 = x * x;
			deviceFP x4 = x2 * x2;
			return x * (1.0 - 2.0 * x2 / 3.0 + 4.0 * x4 / 15.0 - 8.0 * x2 * x4 / 105.0 + (16.0 / 945) * x4 * x4 - (32.0 / 10395) * x4 * x4 * x2);
		}

		int n0 = 2 * int(round(0.5 * x / h));
		deviceFP x0 = h * n0;
		deviceFP xp = x - x0;
		deviceFP d = 0.0;
		for (int i = -n; i < n; i++) {
			if (i % 2 != 0) {
				d += exp(-(xp - i * h) * (xp - i * h)) / (i + n0);
			}
		}
		return INVSQRTPI * d;
	}

	//Inner function for the Sellmeier equation to provide the refractive indicies
	//current equation form:
	//n^2 = a[0] //background (high freq) contribution
	//      + four resonances, purely real contribution
	//      + parametrized low-frequency correction
	//      + 2 complex-valued Lorenzian contribution
	//inputs:
	//a: 22 component array of the coefficients
	//ls: lambda^2 (microns^2)
	//omega: frequency (rad/s)
	//ii: sqrt(-1)
	//kL: 3183.9 i.e. (e * e / (epsilon_o * m_e)
	deviceFunction deviceComplex sellmeierFunc(deviceFP ls, deviceFP omega,const deviceFP* a, int eqn) {
		deviceFP omega2 = omega * omega;
		deviceFP realPart;
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
			compPart = (deviceFP)KLORENTZIAN*a[16] / deviceComplex(a[17] - omega2, a[18] * omega)
				+ (deviceFP)KLORENTZIAN*a[19] / deviceComplex(a[20] - omega2, a[21] * omega);
			return deviceLib::sqrt(maxN(realPart, 0.9f) + compPart);
		case 1:
			//up to 7 Lorentzian lines
			compPart = a[0] / deviceComplex(a[1] - omega2, a[2] * omega)
				+ a[3] / deviceComplex(a[4] - omega2, a[5] * omega)
				+ a[6] / deviceComplex(a[7] - omega2, a[8] * omega)
				+ a[9] / deviceComplex(a[10] - omega2, a[11] * omega)
				+ a[12] / deviceComplex(a[13] - omega2, a[14] * omega)
				+ a[15] / deviceComplex(a[16] - omega2, a[17] * omega)
				+ a[18] / deviceComplex(a[19] - omega2, a[20] * omega);
			compPart *= KLORENTZIAN;
			return deviceLib::sqrt((deviceFP)1.0f + compPart);
		case 2:
		{
			//Up to 7 complex Gaussian functions
			//the real part is the Hilbert transform of the Gaussian, th
			deviceFP scaledF;
			compPart = deviceComplex(a[0], 0.0);
			for (int i = 0; i < 7; ++i) {
				if (a[3 * i + 1] != 0.0){
					scaledF = (omega - a[1 + 3 * i]) / (deviceFPLib::sqrt(2.0) * a[2 + 3 * i]);
					compPart += a[3 + 3 * i] * deviceComplex(-INVSQRTPI * deviceDawson(scaledF), -deviceFPLib::exp(-scaledF * scaledF));
				}

			}
			//always select branch with imaginary part < 0
			return deviceComplex((deviceLib::sqrt(compPart)).real(), -deviceFPLib::abs((deviceLib::sqrt(compPart)).imag()));
		}
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
			return deviceComplex(deviceFPLib::sqrt(maxN(realPart, 0.9)), 0.0);
		}
		
		return deviceComplex(1.0, 0.0);
	};

	//Sellmeier equation for refractive indicies
	deviceFunction deviceComplex sellmeierCuda(
		deviceComplex* ne, deviceComplex* no, const deviceFP* a, deviceFP f, deviceFP theta, deviceFP phi, int type, int eqn) {
		if (f == 0) {
			*ne = deviceComplex(1.0, 0.0); 
			*no = deviceComplex(1.0, 0.0); 
			return deviceComplex(1.0, 0.0);
		} //exit immediately for f=0
		deviceFP ls = 2.99792458e14 / f; //wavelength in microns
		ls *= ls; //only wavelength^2 is ever used
		deviceFP omega = TWOPI * maxN(f,-f);

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
			deviceFP cosT = deviceFPLib::cos(theta);
			cosT *= cosT;
			deviceFP sinT = deviceFPLib::sin(theta);
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
			deviceFP cp = deviceFPLib::cos(phi);
			cp *= cp;
			deviceFP sp = deviceFPLib::sin(phi);
			sp *= sp;
			deviceFP ct = deviceFPLib::cos(theta);
			ct *= ct;
			deviceFP st = deviceFPLib::sin(theta);
			st *= st;

			*ne = deviceLib::sqrt(na * nb * nc /
				(na * nb * st + na * nc * sp * ct + nb * nc * cp * ct));
			*no = deviceLib::sqrt(na * nb /
				(na * cp + nb * sp));

			return *ne;
		}
	}

	deviceFunction deviceFP cuCModSquared(const deviceComplex& a) {
		return a.real() * a.real() + a.imag() * a.imag();
	}

	//provide a list of nearest-3 neighbors for taking spatial derivatives
	// exploiting the fact that the radial grid is offset by 1/4 step from 0
	// this means that midpoints are available on the other side of the origin.
	// returns rho at the given index j
	deviceFunction deviceFP resolveNeighborsInOffsetRadialSymmetry(
		long long* neighbors, long long N, int j, deviceFP dr, long long Ntime, long long h) {
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
	//provide the position rho in cylindric mode; a simplified
	//version of the resolveNeighbors function above for cases where
	//the neighbors aren't required
	deviceFunction deviceFP rhoInRadialSymmetry(
		long long N, int j, deviceFP dr) {
		if (j < N / 2) {
			return deviceFPLib::abs( - (dr * (j - N / 2) + 0.25 * dr));
		}
		else {
			return deviceFPLib::abs(dr * (j - N / 2) + 0.25 * dr);
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
	deviceFunction void findBirefringentCrystalIndex(const deviceParameterSet* s, const deviceFP* sellmeierCoefficients, long long i, deviceComplex* n1, deviceComplex* n2) {
		unsigned long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;

		deviceFP f = (*s).fStep * h;
		deviceFP kx1 = (LIGHTC / (TWOPI * f)) * (j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace));
		deviceFP ky1 = (LIGHTC / (TWOPI * f)) * (k * (*s).dk2 - (k >= ((*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2));
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

		deviceFP gradient[2][2] = { {0.0} };
		deviceFP alpha[2] = { deviceFPLib::asin(kx1 / n[0][0].real()),deviceFPLib::asin(kx1 / n[0][1].real()) };
		deviceFP beta[2] = { deviceFPLib::asin(ky1 / n[0][0].real()), deviceFPLib::asin(ky1 / n[0][1].real()) };

		deviceFP gradientStep = 1.0e-6;
		deviceFP gradientFactor = 0.5 / gradientStep;
		int it = 0;
		int maxiter = 64;
		deviceFP gradientTol = 1e-3;
		//emperical testing: 
		// converges to double precision limit in two iterations for BBO
		// converges in 32 iterations in BiBO

		deviceFP errArray[4][2] = { {0.0} };
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
				if (deviceFPLib::abs(gradient[0][0]) > gradientTol) {
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
				if (deviceFPLib::abs(gradient[0][0]) > 1e-2) alpha[0] -= 0.25 * (errArray[0][0] + errArray[1][0]) / gradient[0][0];
				if (deviceFPLib::abs(gradient[1][0]) > 1e-2) beta[0] -= 0.25 * (errArray[2][0] + errArray[3][0]) / gradient[1][0];
				if (deviceFPLib::abs(gradient[0][1]) > 1e-2) alpha[1] -= 0.25 * (errArray[0][1] + errArray[1][1]) / gradient[0][1];
				if (deviceFPLib::abs(gradient[1][1]) > 1e-2) beta[1] -= 0.25 * (errArray[2][1] + errArray[3][1]) / gradient[1][1];

				if (maxN(maxN(deviceFPLib::abs(gradient[0][0]), deviceFPLib::abs(gradient[1][0])), maxN(deviceFPLib::abs(gradient[0][1]), deviceFPLib::abs(gradient[1][1]))) < gradientTol) break;
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

	//calculate the total energy spectrum of the beam for the 2D modes. Note that for the 
	//cartesian one, it will be treated as a round beam instead of an infinite plane wave 
	//in the transverse direction. Thus, the 2D Cartesian spectra are approximations.
	trilingual totalSpectrumKernel asKernel(withID const deviceParameterSet* s) {
		size_t i = localIndex;
		size_t j;
		deviceFP beamCenter1 = 0.;
		deviceFP beamCenter2 = 0.;
		deviceFP beamTotal1 = 0.;
		deviceFP beamTotal2 = 0.;
		deviceFP a, x;

		//find beam centers
		if ((*s).isCylindric) {
			beamCenter1 = ((*s).Nspace / 2 * (*s).dx) + 0.25 * (*s).dx;
			beamCenter2 = beamCenter1;
		}
		else{
			for (j = 0; j < (*s).Nspace; ++j) {
				x = (*s).dx * j;
				a = cuCModSquared((*s).workspace1[i + j * (*s).Nfreq]);
				beamTotal1 += a;
				beamCenter1 += x * a;
				a = cuCModSquared((*s).workspace2[i + j * (*s).Nfreq]);
				beamTotal2 += a;
				beamCenter2 += x * a;
			}
			if (beamTotal1 > 0) {
				beamCenter1 /= beamTotal1;
			}
			if (beamTotal2 > 0) {
				beamCenter2 /= beamTotal2;
			}
		}

		//Integrate total beam power, assuming radially-symmetric beam around
		//the center
		beamTotal1 = 0.;
		beamTotal2 = 0.;
		for (j = 0; j < (*s).Nspace; ++j) {
			x = (*s).dx * j;
			beamTotal1 += deviceFPLib::abs(x - beamCenter1) * cuCModSquared((*s).workspace1[i + j * (*s).Nfreq]);
			beamTotal2 += deviceFPLib::abs(x - beamCenter2) * cuCModSquared((*s).workspace2[i + j * (*s).Nfreq]);
		}
		beamTotal1 *= 2 * PI * LIGHTC * EPS0 * (*s).dx * (*s).dt * (*s).dt;
		beamTotal2 *= 2 * PI * LIGHTC * EPS0 * (*s).dx * (*s).dt * (*s).dt;

		//put the values into the output spectrum
		(*s).gridPolarizationTime1[i] = beamTotal1;
		(*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
		(*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
	};

	//Calculate the energy spectrum after a 2D propagation assuming that the beam
	//height in the non-resolved direction is == the grid width (i.e. square grid)
	//More quantitative than the mapping to round beams, but rather specific
	trilingual totalSpectrum2DSquareKernel asKernel(withID const deviceParameterSet* s) {
		size_t i = localIndex;
		size_t j;

		deviceFP beamTotal1 = 0.;
		deviceFP beamTotal2 = 0.;
		//Integrate total beam power
		beamTotal1 = 0.;
		beamTotal2 = 0.;
		for (j = 0; j < (*s).Nspace; ++j) {
			beamTotal1 += cuCModSquared((*s).workspace1[i + j * (*s).Nfreq]);
			beamTotal2 += cuCModSquared((*s).workspace2[i + j * (*s).Nfreq]);
		}
		beamTotal1 *= 2 * LIGHTC * EPS0 * (*s).dx * (*s).dx * (*s).Nspace * (*s).dt * (*s).dt;
		beamTotal2 *= 2 * LIGHTC * EPS0 * (*s).dx * (*s).dx * (*s).Nspace * (*s).dt * (*s).dt;

		//put the values into the output spectrum
		(*s).gridPolarizationTime1[i] = beamTotal1;
		(*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
		(*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
	};

	//Calculate the energy spectrum after a 3D propagation
	trilingual totalSpectrum3DKernel asKernel(withID const deviceParameterSet* s) {
		size_t i = localIndex;
		size_t j;

		deviceFP beamTotal1 = 0.;
		deviceFP beamTotal2 = 0.;
		//Integrate total beam power
		beamTotal1 = 0.;
		beamTotal2 = 0.;
		for (j = 0; j < (*s).Nspace * (*s).Nspace2; ++j) {
			beamTotal1 += cuCModSquared((*s).workspace1[i + j * (*s).Nfreq]);
			beamTotal2 += cuCModSquared((*s).workspace2[i + j * (*s).Nfreq]);
		}
		beamTotal1 *= 2 * LIGHTC * EPS0 * (*s).dx * (*s).dx * (*s).dt * (*s).dt;
		beamTotal2 *= 2 * LIGHTC * EPS0 * (*s).dx * (*s).dx * (*s).dt * (*s).dt;

		//put the values into the output spectrum
		(*s).gridPolarizationTime1[i] = beamTotal1;
		(*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
		(*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
	};

	//perform a Hankel transform by direct quadrature
	//the offset radial grid allows the sum to be done with a midpoint method
	//with no numerical effort and the rho=0 point is excluded from the grid
	//this function is slow and order N^2 as it is not used in the core loop.
	//the numerical accuracy of Hankel transforms that I've seen in relatively
	//low due to Gibbs phenomena and I find the FFT-based propagation implemented
	//below better for nonlinear phenomena. I might later use this for linear propagation
	//in sequences however.
	trilingual hankelKernel asKernel(withID const deviceParameterSet* s, deviceFP* in, deviceFP* out) {
		size_t i = localIndex;
		size_t col = i / (*s).Ntime; //spatial coordinate
		deviceFP dk = 2.0 / (PI * (*s).dx * (*s).Nspace);
		in += i % (*s).Ntime;
		out[i] = 0.0;
		out[i + (*s).Ngrid] = 0.0;
		deviceFP r0;
		deviceFP J0 = 1.0;
		deviceFP k0 = col * dk;
		for (size_t r = 0; r < (*s).Nspace; ++r) {
			r0 = rhoInRadialSymmetry((*s).Nspace, r, (*s).dx);
			J0 = r0 * j0Device(r0 * k0);
			out[i] += J0 * in[r * (*s).Ntime];
			out[i + (*s).Ngrid] += J0 * in[r * (*s).Ntime + (*s).Ngrid];
		}
		out[i] *= (*s).dx;
		out[i + (*s).Ngrid] *= (*s).dx;
	};

	//inverse Hankel transform from the k-space back to the offset spatial grid
	trilingual inverseHankelKernel asKernel(withID const deviceParameterSet* s, deviceFP* in, deviceFP* out) {
		size_t i = localIndex;
		size_t col = i / (*s).Ntime; //spatial coordinate
		deviceFP dk = 2.0 / (PI * (*s).dx * (*s).Nspace);
		in += i % (*s).Ntime;;
		out[i] = 0.0;
		out[i + (*s).Ngrid] = 0.0;
		deviceFP r0 = rhoInRadialSymmetry((*s).Nspace, col, (*s).dx);
		deviceFP J0 = 1.0;
		deviceFP k0 = col * dk;
		for (size_t k = 0; k < (*s).Nspace; ++k) {
			k0 = k * dk;
			J0 = k0*j0Device(r0 * k0);
			out[i] += J0 * in[k * (*s).Ntime];
			out[i + (*s).Ngrid] += J0 * in[k * (*s).Ntime + (*s).Ngrid];
		}
		out[i] *= 0.5 * dk / ((*s).Ntime);
		out[i + (*s).Ngrid] *= 0.5 * dk / ((*s).Ntime);
	};

	//rotate the field around the propagation axis (basis change)
	trilingual rotateFieldKernel asKernel(withID deviceComplex* Ein1, deviceComplex* Ein2, deviceComplex* Eout1,
		deviceComplex* Eout2, deviceFP rotationAngle) {
		long long i = localIndex;
		Eout1[i] = deviceFPLib::cos(rotationAngle) * Ein1[i] - deviceFPLib::sin(rotationAngle) * Ein2[i];
		Eout2[i] = deviceFPLib::sin(rotationAngle) * Ein1[i] + deviceFPLib::cos(rotationAngle) * Ein2[i];
	};

	//calculate the extra term in the Laplacian encountered in cylindrical coordinates (1/r d/drho)
	trilingual radialLaplacianKernel asKernel(withID const deviceParameterSet* s) {
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
			deviceFP rho = resolveNeighborsInOffsetRadialSymmetry(neighbors, (*s).Nspace, (int)j, (*s).dx, (*s).Ntime, h);
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

	//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
	//note that the sellmeier coefficients have extra values appended to the end
	//to give info about the current simulation
	trilingual applyFresnelLossKernel asKernel(withID const deviceFP* sellmeierCoefficients1, const deviceFP* sellmeierCoefficients2, const deviceParameterSet* s) {
		//the old version was hopelessly broken, make a new one from scratch.
	};

	trilingual apertureFarFieldKernel asKernel(withID const deviceParameterSet* s, deviceFP radius, deviceFP activationParameter, deviceFP xOffset, deviceFP yOffset) {
		long long i = localIndex;
		long long col, j, k, l;
		deviceComplex cuZero = deviceComplex(0, 0);
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;
		k = col % (*s).Nspace;
		l = col / (*s).Nspace;

		//magnitude of k vector
		deviceFP ko = TWOPI*j * (*s).fStep/LIGHTC;

		//transverse wavevector being resolved
		deviceFP dk1 = k * (*s).dk1 - (k >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (long long)(*s).Nspace); //frequency grid in x direction
		deviceFP dk2 = 0.0;
		if((*s).is3D) dk2 = l * (*s).dk2 - (l >= ((long long)(*s).Nspace2 / 2)) * ((*s).dk2 * (long long)(*s).Nspace2); //frequency grid in y direction

		//light that won't go the the farfield is immediately zero
		if (dk1*dk1 > ko*ko || dk2*dk2 > ko*ko) {
			(*s).gridEFrequency1[i] = cuZero;
			(*s).gridEFrequency2[i] = cuZero;
			return;
		}

		deviceFP theta1 = deviceFPLib::asin(dk1 / ko);
		deviceFP theta2 = deviceFPLib::asin(dk2 / ko);

		theta1 -= (!(*s).isCylindric) *  xOffset;
		theta2 -= (*s).is3D * yOffset;

		deviceFP r = deviceFPLib::sqrt(theta1 * theta1 + theta2 * theta2);

		deviceFP a = 1.0 - (1.0 / (1.0 + deviceFPLib::exp(-activationParameter * (r-radius))));

		(*s).gridEFrequency1[i] *= a;
		(*s).gridEFrequency2[i] *= a;
	};

	trilingual apertureFarFieldKernelHankel asKernel(withID const deviceParameterSet* s, deviceFP radius, deviceFP activationParameter, deviceFP xOffset, deviceFP yOffset) {
		long long i = localIndex;
		long long col, j, k;
		deviceComplex cuZero = deviceComplex(0, 0);
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;
		k = col % (*s).Nspace;

		//magnitude of k vector
		deviceFP ko = TWOPI * j * (*s).fStep / LIGHTC;

		//transverse wavevector being resolved
		deviceFP dk1 = k * 2.0 / (PI * (*s).dx * (*s).Nspace);; //frequency grid in x direction

		//light that won't go the the farfield is immediately zero
		if (dk1 * dk1 > ko * ko) {
			(*s).gridEFrequency1[i] = cuZero;
			(*s).gridEFrequency2[i] = cuZero;
			return;
		}

		deviceFP theta1 = deviceFPLib::asin(dk1 / ko);

		deviceFP a = 1.0 - (1.0 / (1.0 + deviceFPLib::exp(-activationParameter * (deviceFPLib::abs(theta1) - radius))));

		(*s).gridEFrequency1[i] *= a;
		(*s).gridEFrequency2[i] *= a;
	};


	//apply a spectral filter to the beam (full details in docs)
	trilingual filterKernel asKernel(withID const deviceParameterSet* s, deviceFP f0, deviceFP bandwidth, deviceFP order, deviceFP inBandAmplitude, deviceFP outOfBandAmplitude) {
		long long i = localIndex;
		long long col, j;
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;

		deviceFP f = (*s).fStep * j - f0;
		for (int p = 1; p < (int)order; p++) {
			bandwidth *= bandwidth;
			f *= f;
		}
		deviceFP filterFunction = outOfBandAmplitude + inBandAmplitude*deviceFPLib::exp(-f / (2.0 * bandwidth));
		(*s).gridEFrequency1[i] *= filterFunction;
		(*s).gridEFrequency2[i] *= filterFunction;
	};

	//apply a lorentzian gain or loss in a certain cross-section of the beam.
	// amplitude - strength of the copy of the pulse applied
	// f0 - resonance frequency of the Lorentzian (THz)
	// gamma - linewidth (radHz)
	// radius - radius of the spot (m)
	// order - supergaussian order of the spot shape
	trilingual lorentzianSpotKernel asKernel(withID const deviceParameterSet* s, deviceFP amplitude, deviceFP f0, deviceFP gamma, deviceFP radius, deviceFP order) {
		long long i = localIndex;
		long long j, h, k, col;

		col = i / ((*s).Nfreq - 1);
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;
		deviceFP r, f, x, y;
		if ((*s).is3D) {
			h = 1 + i % ((*s).Nfreq - 1);
			col = i / ((*s).Nfreq - 1);
			i = h + col * ((*s).Nfreq);
			j = col % (*s).Nspace;
			k = col / (*s).Nspace;
			f = h * (*s).fStep;

			x = ((*s).dx * (j - (*s).Nspace / 2.0));
			y = ((*s).dx * (k - (*s).Nspace2 / 2.0));
			r = sqrt(x * x + y * y);
		}
		else {
			h = 1 + i % ((*s).Nfreq - 1);
			j = i / ((*s).Nfreq - 1);
			i = h + j * ((*s).Nfreq);
			f = h * (*s).fStep;
			r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0) + 0.25 * (*s).dx);
		}

		deviceFP w0 = TWOPI * f0;
		deviceFP w = TWOPI * f;
		deviceComplex lorentzian = gamma * w0 * amplitude / (w0 * w0 - w * w + deviceComplex(0.0, gamma * w));
		deviceFP spotFactor = r / radius;
		for (int p = 1; p < (int)order; p++) {
			spotFactor *= spotFactor;
		}
		deviceComplex filterFunction = deviceComplex(0.0, exp(-spotFactor)) * lorentzian;
		(*s).gridEFrequency1[i] += filterFunction * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] += filterFunction * (*s).gridEFrequency2[i];
	};

	//Apply a (soft, possibly) aperture
	trilingual apertureKernel asKernel(withID const deviceParameterSet* s, deviceFP radius, deviceFP activationParameter) {
		long long i = localIndex;
		long long j, k, col;

		col = i / (*s).Ntime;
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;
		deviceFP r;
		if ((*s).is3D) {
			deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0));
			deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0));
			r = sqrt(x * x + y * y);
		}
		else {
			r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0) + 0.25 * (*s).dx);
		}

		deviceFP a = 1.0 - (1.0 / (1.0 + deviceFPLib::exp(-activationParameter * (r - radius) / (*s).dx)));

		//if (r>radius) a = 0;
		(*s).gridETime1[i] *= a;
		(*s).gridETime2[i] *= a;
	};

	//apply a spatial phase corresponding to a parabolic mirror (on-axis)
	trilingual parabolicMirrorKernel asKernel(withID const deviceParameterSet* s, deviceFP focus) {
		long long i = localIndex;
		long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		i = h + col * (*s).Nfreq;
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;

		deviceFP w = TWOPI * h * (*s).fStep;
		deviceFP r;
		if ((*s).is3D) {
			deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0));
			deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0));
			r = sqrt(x * x + y * y);
		}
		else {
			r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0) + 0.25 * (*s).dx);
		}

		deviceComplex	u = deviceLib::exp(deviceComplex(0.0,
			w * r * r * (0.5 / focus) / LIGHTC));

		(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
	};

	//apply a spatial phase corresponding to a spherical mirror (on axis)
	trilingual sphericalMirrorKernel asKernel(withID const deviceParameterSet* s, deviceFP ROC) {
		long long i = localIndex;
		long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		i = h + col * (*s).Nfreq;
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;

		deviceFP w = TWOPI * h * (*s).fStep;
		deviceFP r;
		if ((*s).is3D) {
			deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0));
			deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0));
			r = deviceFPLib::sqrt(x * x + y * y);
		}
		else {
			r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0) + 0.25 * (*s).dx);
		}

		bool isNegative = ROC < 0.0;
		ROC = deviceFPLib::abs(ROC);
		deviceComplex u = deviceComplex(0.0, 0.0);
		if (r < ROC) {
			u = deviceLib::exp(deviceComplex(0.0,
				2.0 * deviceFPLib::pow(-1, isNegative) * w * ROC * ((deviceFPLib::sqrt(1.0 - r * r / (ROC * ROC))) - 1.0) / LIGHTC));
		}

		(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
		if (isnan((*s).gridEFrequency1[i].real()) || isnan((*s).gridEFrequency2[i].real()) || isnan((*s).gridEFrequency1[i].imag()) || isnan((*s).gridEFrequency2[i].imag())) {
			(*s).gridEFrequency1[i] = deviceComplex(0.0, 0.0);
			(*s).gridEFrequency2[i] = deviceComplex(0.0, 0.0);
		}
	};

	//apply linear propagation through a given medium to the fields
	trilingual applyLinearPropagationKernel asKernel(withID const deviceFP* sellmeierCoefficients, deviceFP thickness, const deviceParameterSet* s) {
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
		deviceFP crystalTheta = sellmeierCoefficients[66];
		deviceFP crystalPhi = sellmeierCoefficients[67];

		//frequency being resolved by current thread
		deviceFP f = h * (*s).fStep;
		deviceFP omega = TWOPI * f;
		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);
		deviceFP dk1 = j * (*s).dk1 - (j >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace);
		deviceFP dk2 = k * (*s).dk2 - (k >= ((long long)(*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2);
		if (!(*s).is3D)dk2 = 0.0;
		//if ((*s).isCylindric) dk2 = dk1;
		sellmeierCuda(&n0, &n0o, sellmeierCoefficients, (*s).f0,
			crystalTheta, crystalPhi, axesNumber, sellmeierType);
		if (isnan(ne.real()) || isnan(no.real())) {
			ne = deviceComplex(1, 0);
			no = deviceComplex(1, 0);
		}

		deviceComplex ke = ne * omega / (deviceFP)LIGHTC;
		deviceComplex ko = no * omega / (deviceFP)LIGHTC;
		deviceFP k0 = (n0 * omega / (deviceFP)LIGHTC).real();
		deviceFP kze = ke.real() * ke.real() - dk1 * dk1 - dk2 * dk2;
		deviceFP kzo = ko.real() * ko.real() - dk1 * dk1 - dk2 * dk2;
			//deviceFP kze = (deviceLib::sqrt(ke * ke - dk1 * dk1 - dk2 * dk2)).real();
			//deviceFP kzo = (deviceLib::sqrt(ko * ko - dk1 * dk1 - dk2 * dk2)).real();

		deviceComplex ts = deviceComplex(0.0, 0.0);//deviceLib::exp(ii * (k0 - kze) * thickness);
		deviceComplex tp = deviceComplex(0.0, 0.0);// deviceLib::exp(ii * (k0 - kzo) * thickness);
		if(kze>=0.0) ts = deviceLib::exp(ii * (k0 - deviceFPLib::sqrt(kze)) * thickness);
		if(kzo>=0.0) tp = deviceLib::exp(ii * (k0 - deviceFPLib::sqrt(kzo)) * thickness);
		if (isnan(ts.real()) || isnan(ts.imag())) ts = deviceComplex(0, 0);
		if (isnan(tp.real()) || isnan(tp.imag())) tp = deviceComplex(0, 0);
		(*s).gridEFrequency1[i] = ts * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] = tp * (*s).gridEFrequency2[i];
	};

	//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
	//note that the sellmeier coefficients have extra values appended to the end
	//to give info about the current simulation
	trilingual prepareCartesianGridsKernel asKernel(withID const deviceFP* sellmeierCoefficients, const deviceParameterSet* s) {
		long long i = localIndex;
		long long j, k;
		deviceComplex ne, no;
		deviceComplex n0 = (*s).n0;
		deviceComplex cuZero = deviceComplex(0, 0);
		j = i / ((*s).Nfreq - 1); //spatial coordinate
		k = 1 + (i % ((*s).Nfreq - 1)); //temporal coordinate
		i = k + j * (*s).Nfreq;
		deviceComplex ii = deviceComplex(0, 1);
		deviceFP kStep = sellmeierCoefficients[70];
		deviceFP fStep = sellmeierCoefficients[71];

		//frequency being resolved by current thread
		deviceFP f = k * fStep;

		//transverse wavevector being resolved
		deviceFP dk = j * kStep - (j >= ((long long)(*s).Nspace / 2)) * (kStep * (*s).Nspace); //frequency grid in transverse direction
		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);

		//if the refractive index was returned weird, then the index isn't valid, so set the propagator to zero for that frequency
		if (minN(ne.real(), no.real()) < 0.9 || isnan(ne.real()) || isnan(no.real()) || isnan(ne.imag()) || isnan(no.imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			return;
		}

		deviceComplex k0 = deviceComplex((deviceFP)TWOPI * n0.real() * f / LIGHTC, 0);
		deviceComplex ke = (deviceFP)TWOPI * ne * f / (deviceFP)LIGHTC;
		deviceComplex ko = (deviceFP)TWOPI * no * f / (deviceFP)LIGHTC;

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
		deviceComplex kz1 = deviceLib::sqrt(ke * ke - dk * dk);
		deviceComplex kz2 = deviceLib::sqrt(ko * ko - dk * dk);

		if (kz1.real() > 0.0 && kz2.real() > 0.0){
			(*s).gridPropagationFactor1[i] = deviceLib::exp(-(deviceFP)0.5 * ii * (kz1 - k0) * (*s).h);
			if (isnan(((*s).gridPropagationFactor1[i]).real())) {
				(*s).gridPropagationFactor1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = deviceLib::exp(-(deviceFP)0.5 * ii * (kz2 - k0) * (*s).h);
			if (isnan(((*s).gridPropagationFactor2[i]).real())) {
				(*s).gridPropagationFactor2[i] = cuZero;
			}

			(*s).gridPolarizationFactor1[i] = -ii * deviceLib::pow((deviceComplex)(*s).chiLinear1[k] + (deviceFP)1.0, (deviceFP)0.25) * chi11 * ((deviceFP)TWOPI * (deviceFP)TWOPI * f * f) / ((2 * (deviceFP)LIGHTC * (deviceFP)LIGHTC * kz1)) * (*s).h;
			(*s).gridPolarizationFactor2[i] = -ii * deviceLib::pow((deviceComplex)(*s).chiLinear2[k] + (deviceFP)1.0, (deviceFP)0.25) * chi12 * ((deviceFP)TWOPI * (deviceFP)TWOPI * f * f) / ((2 * (deviceFP)LIGHTC * (deviceFP)LIGHTC * kz2)) * (*s).h;
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
	trilingual prepare3DGridsKernel asKernel(withID const deviceFP* sellmeierCoefficients, const deviceParameterSet* s) {
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
		deviceFP f = j * (*s).fStep;

		//transverse wavevector being resolved
		deviceFP dk1 = k * (*s).dk1 - (k >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (long long)(*s).Nspace); //frequency grid in x direction
		deviceFP dk2 = l * (*s).dk2 - (l >= ((long long)(*s).Nspace2 / 2)) * ((*s).dk2 * (long long)(*s).Nspace2); //frequency grid in y direction

		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);
		if (minN(ne.real(), no.real()) < 0.9 || isnan(ne.real()) || isnan(no.real()) || isnan(ne.imag()) || isnan(no.imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			return;
		}

		deviceComplex k0 = deviceComplex(TWOPI * n0.real() * f / LIGHTC, 0);
		deviceComplex ke = (deviceFP)TWOPI * ne * f / (deviceFP)LIGHTC;
		deviceComplex ko = (deviceFP)TWOPI * no * f / (deviceFP)LIGHTC;

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
		deviceComplex kz1 = deviceLib::sqrt(ke * ke - dk1 * dk1 - dk2 * dk2);
		deviceComplex kz2 = deviceLib::sqrt(ko * ko - dk1 * dk1 - dk2 * dk2);

		if (kz1.real() > 0.0 && kz2.real() > 0.0) {
			(*s).gridPropagationFactor1[i] = deviceLib::exp(-(deviceFP)0.5 * ii * (kz1 - k0) * (*s).h);
			if (isnan(((*s).gridPropagationFactor1[i].real()))) {
				(*s).gridPropagationFactor1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = deviceLib::exp(-(deviceFP)0.5 * ii * (kz2 - k0) * (*s).h);
			if (isnan(((*s).gridPropagationFactor2[i].real()))) {
				(*s).gridPropagationFactor2[i] = cuZero;
			}

			(*s).gridPolarizationFactor1[i] = -ii * deviceLib::pow((deviceComplex)(*s).chiLinear1[j] + (deviceFP)1.0, (deviceFP)0.25) * chi11 * ((deviceFP)TWOPI * (deviceFP)TWOPI * f * f) / (2 * (deviceFP)LIGHTC * (deviceFP)LIGHTC * kz1) * (*s).h;
			(*s).gridPolarizationFactor2[i] = -ii * deviceLib::pow((deviceComplex)(*s).chiLinear2[j] + (deviceFP)1.0, (deviceFP)0.25) * chi12 * ((deviceFP)TWOPI * (deviceFP)TWOPI * f * f) / (2 * (deviceFP)LIGHTC * (deviceFP)LIGHTC * kz2) * (*s).h;
		}

		else {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
		}
	};

	//prepare the chi(1) arrays that will be needed in the simulation
	trilingual getChiLinearKernel asKernel(withID deviceParameterSet* s, const deviceFP* sellmeierCoefficients) {
		long long i = localIndex;
		int axesNumber = (*s).axesNumber;
		int sellmeierType = (*s).sellmeierType;
		deviceFP crystalTheta = sellmeierCoefficients[66];
		deviceFP crystalPhi = sellmeierCoefficients[67];
		deviceFP fStep = sellmeierCoefficients[71];

		deviceComplex ne, no;

		//frequency being resolved by current thread
		deviceFP f = i * fStep;
		
		sellmeierCuda(&ne, &no, sellmeierCoefficients, deviceFPLib::abs(f), crystalTheta, crystalPhi, axesNumber, sellmeierType);
		if (isnan(ne.real()) || isnan(no.real())) {
			ne = deviceComplex(1.0, 0);
			no = ne;
		}

		(*s).chiLinear1[i] = ne * ne - (deviceFP)1.0;
		(*s).chiLinear2[i] = no * no - (deviceFP)1.0;
		if ((*s).chiLinear1[i].real() != 0.0 && (*s).chiLinear2[i].real() != 0.0) {
			(*s).inverseChiLinear1[i] = (deviceFP)1.0 / (*s).chiLinear1[i].real();
			(*s).inverseChiLinear2[i] = (deviceFP)1.0 / (*s).chiLinear2[i].real();
		}
		else {
			(*s).inverseChiLinear1[i] = 0.0;
			(*s).inverseChiLinear2[i] = 0.0;
		}

		(*s).fieldFactor1[i] = 1.0 / deviceFPLib::pow((*s).chiLinear1[i].real() + (deviceFP)1.0, (deviceFP)0.25); //account for the effective field strength in the medium (1/n)
		(*s).fieldFactor2[i] = 1.0 / deviceFPLib::pow((*s).chiLinear2[i].real() + (deviceFP)1.0, (deviceFP)0.25);
		if ((*s).isUsingMillersRule) {
			(*s).fieldFactor1[i] *= (*s).chiLinear1[i].real();
			(*s).fieldFactor2[i] *= (*s).chiLinear2[i].real();
		}
		if ((*s).chiLinear1[i].real() < -(deviceFP)0.1) {
			(*s).fieldFactor1[i] = 0.0;
		}
		if ((*s).chiLinear2[i].real() < -(deviceFP)0.1) {
			(*s).fieldFactor2[i] = 0.0;
		}

		if (i == 81) {
			deviceComplex n0;
			sellmeierCuda(&n0, &no, sellmeierCoefficients, deviceFPLib::abs((*s).f0), crystalTheta, crystalPhi, axesNumber, sellmeierType);
			(*s).n0 = no;
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
			const deviceFP* referenceFrequencies = &sellmeierCoefficients[72];
			deviceFP chi11[7];

			for (int im = (i>17)*3; im < 7; ++im) {
				if (referenceFrequencies[im] == 0) {
					chi11[im] = 100000.0;
				}
				else {
					sellmeierCuda(&ne, &no, sellmeierCoefficients, referenceFrequencies[im], sellmeierCoefficients[66], sellmeierCoefficients[67], axesNumber, sellmeierType);
					chi11[im] = no.real() * no.real() - 1.0;
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
	trilingual prepareCylindricGridsKernel asKernel(withID deviceFP* sellmeierCoefficients, deviceParameterSet* s) {
		long long i = localIndex;
		long long j, k;
		long long Nspace = (*s).Nspace;
		deviceComplex cuZero = deviceComplex(0.0, 0.0);
		j = i / ((*s).Nfreq - 1); //spatial coordinate
		k = 1 + i % ((*s).Nfreq - 1); //temporal coordinate
		i = k + j * (*s).Nfreq;
		deviceComplex ii = deviceComplex(0.0, 1.0);
		deviceFP kStep = sellmeierCoefficients[70];
		deviceFP fStep = sellmeierCoefficients[71];

		deviceComplex ne, no;
		deviceComplex n0 = (*s).n0;

		//frequency being resolved by current thread
		deviceFP f = -k * fStep;

		//transverse wavevector being resolved
		deviceFP dk = j * kStep - (j >= (Nspace / 2)) * (kStep * Nspace); //frequency grid in transverse direction

		sellmeierCuda(&ne, &no, sellmeierCoefficients,fStep*k, sellmeierCoefficients[66], sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);

		//if the refractive index was returned weird, then the index isn't valid, so set the propagator to zero for that frequency
		if (minN(ne.real(), no.real()) < 0.95 || ne.real() > 6.0 || no.real() > 6.0 || isnan(ne.real()) || isnan(no.real()) || isnan(ne.imag()) || isnan(no.imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			(*s).gridPropagationFactor1Rho1[i] = cuZero;
			(*s).gridPropagationFactor1Rho2[i] = cuZero;
			return;
		}

		deviceComplex k0 = deviceComplex(TWOPI * n0.real() * f / LIGHTC, 0);
		deviceComplex ke = (deviceFP)TWOPI * ne * f / (deviceFP)LIGHTC;
		deviceComplex ko = (deviceFP)TWOPI * no * f / (deviceFP)LIGHTC;

		deviceComplex chi11 = deviceComplex(1.0, 0);
		deviceComplex chi12 = deviceComplex(1.0, 0);
		if ((*s).isUsingMillersRule) {
			chi11 = (*s).chiLinear1[k];
			chi12 = (*s).chiLinear2[k];
		}

		if ((dk * dk < minN(ke.real() * ke.real() + ke.imag() * ke.imag(), ko.real() * ko.real() + ko.imag() * ko.imag())) && (*s).fieldFactor1[k] > 0.0 && (*s).fieldFactor2[k] > 0.0) {
			(*s).gridPropagationFactor1[i] = deviceLib::exp((deviceFP)0.5 * ii * (ke - k0 - dk * dk / (2 * ke.real())) * (*s).h);
			(*s).gridPropagationFactor1Rho1[i] = ii * (*s).h / ((*s).fieldFactor1[k] * 2 * ke);
			if (isnan((deviceLib::abs((*s).gridPropagationFactor1Rho1[i]+(*s).gridPropagationFactor1[i])))) {
				(*s).gridPropagationFactor1[i] = cuZero;
				(*s).gridPropagationFactor1Rho1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = deviceLib::exp((deviceFP)0.5 * ii * (ko - k0 - dk * dk / (2 * ko.real())) * (*s).h);
			(*s).gridPropagationFactor1Rho2[i] = ii * (*s).h / ((*s).fieldFactor2[k] * 2 * ko);
			if (isnan((deviceLib::abs((*s).gridPropagationFactor1Rho2[i]+(*s).gridPropagationFactor2[i])))) {
				(*s).gridPropagationFactor2[i] = cuZero;
				(*s).gridPropagationFactor1Rho2[i] = cuZero;
			}

			//factor of 0.5 comes from deviceFPd grid size in cylindrical symmetry mode after expanding the beam
			(*s).gridPolarizationFactor1[i] = (deviceFP)0.5 * deviceLib::pow((deviceComplex)(*s).chiLinear1[k] + (deviceFP)1.0, (deviceFP)0.25) * chi11 * ii * ((deviceFP)TWOPI * f) / (2 * ne.real() * (deviceFP)LIGHTC) * (*s).h;
			(*s).gridPolarizationFactor2[i] = (deviceFP)0.5 * deviceLib::pow((deviceComplex)(*s).chiLinear2[k] + (deviceFP)1.0, (deviceFP)0.25) * chi12 * ii * ((deviceFP)TWOPI * f) / (2 * no.real() * (deviceFP)LIGHTC) * (*s).h;
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


	trilingual realToComplexKernel asKernel(withID deviceFP* in, deviceComplex* out) {
		long long i = localIndex;
		out[i] = deviceComplex(in[i], 0.0);
	};

	trilingual complexToRealKernel asKernel(withID deviceComplex* in, deviceFP* out) {
		long long i = localIndex;
		out[i] = in[i].real();
	};

	trilingual materialPhaseKernel asKernel(withID deviceFP df, size_t Ntime, deviceFP* a, deviceFP f01, deviceFP f02,
		deviceFP thickness1, deviceFP thickness2, deviceFP* phase1, deviceFP* phase2) {
		size_t i = localIndex;
		//frequency being resolved by current thread
		deviceFP f = i * df;
		if (i >= Ntime / 2) {
			f -= df * Ntime;
		}

		//give phase shift relative to group velocity (approximated 
		// with low-order finite difference) so the pulse doesn't move
		deviceComplex ne, no, no0, n0p, n0m;
		sellmeierCuda(&ne, &no, a, deviceFPLib::abs(f), 0, 0, 0, 0);
		f *= TWOPI;
		sellmeierCuda(&ne, &no0, a, f01, 0, 0, 0, 0);
		sellmeierCuda(&ne, &n0p, a, f01 + (deviceFP)1e11, 0, 0, 0, 0);
		sellmeierCuda(&ne, &n0m, a, f01 - (deviceFP)1e11, 0, 0, 0, 0);
		no0 = no0 + f01 * (n0p - n0m) / (deviceFP)2e11;
		phase1[i] = thickness1 * f * (no.real() - no0.real()) / LIGHTC;
		sellmeierCuda(&ne, &no0, a, f02, 0, 0, 0, 0);
		sellmeierCuda(&ne, &n0p, a, f02 + (deviceFP)1e11, 0, 0, 0, 0);
		sellmeierCuda(&ne, &n0m, a, f02 - (deviceFP)1e11, 0, 0, 0, 0);
		no0 = no0 + f02 * (n0p - n0m) / (deviceFP)2e11;
		phase2[i] = thickness2 * f * (no.real() - no0.real()) / LIGHTC;
	};

	//calculate the nonlinear polarization, after FFT to get the field
	//in the time domain
	trilingual nonlinearPolarizationKernel asKernel(withID const deviceParameterSet* s) {
		size_t i = localIndex;
		deviceFP Ex = (*s).gridETime1[i];
		deviceFP Ey = (*s).gridETime2[i];

		if ((*s).nonlinearSwitches[0] == 1 || (*s).nonlinearSwitches[1] == 1){
			
			deviceFP P[3] = {0.0};

			// rotate field into crystal frame
			deviceFP E3[3] = {(*s).rotationForward[0] * Ex + (*s).rotationForward[1] * Ey,
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
				deviceFP Esquared = (*s).chi3Tensor[0] * (Ex*Ex + Ey*Ey);
				(*s).gridPolarizationTime1[i] += Ex * Esquared;
				(*s).gridPolarizationTime2[i] += Ey * Esquared;
			}
		}
		else{
			//case of no chi2, and centrosymmetric chi3
			if ((*s).nonlinearSwitches[1] == 2) {
				deviceFP Esquared = (*s).chi3Tensor[0] * (Ex*Ex + Ey*Ey);
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
		if ((*s).isCylindric)expandCylindricalBeamDevice(s, i, (*s).gridRadialLaplacian1 + (*s).Ngrid * 4 * (*s).hasPlasma, (*s).gridPolarizationTime1, (*s).gridPolarizationTime2);
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
	trilingual plasmaCurrentKernel_twoStage_A asKernel(withID const deviceParameterSet* s) {
		size_t i = localIndex;
		deviceFP Esquared, a;
		unsigned char pMax = (unsigned char)(*s).nonlinearSwitches[3];

		//save values in workspaces, casting to deviceFP
		deviceFP* dN = (deviceFP*)(*s).workspace1;
		deviceFP* dN2 = dN + (*s).Ngrid;
		deviceFP* Jx = (*s).gridPolarizationTime1;
		deviceFP* Jy = (*s).gridPolarizationTime2;
		Esquared = (*s).gridETime1[i] * (*s).gridETime1[i] + (*s).gridETime2[i] * (*s).gridETime2[i];
		Esquared *= (*s).plasmaParameters[0];
		a = Esquared;
		for (unsigned char p = 0; p < pMax; ++p) {
			a *= Esquared;
		}
		Jx[i] = a * (*s).gridETime1[i];
		Jy[i] = a * (*s).gridETime2[i];
		dN[i] = (*s).plasmaParameters[2] * (Jx[i] * (*s).gridETime1[i] + Jy[i] * (*s).gridETime2[i]);
		dN2[i] = dN[i];
	};

	trilingual plasmaCurrentKernel_twoStage_B asKernel(withID const deviceParameterSet* s) {
		size_t j = localIndex;
		j *= (*s).Ntime;
		deviceFP N = 0;
		deviceFP integralx = 0;
		deviceFP* expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
		deviceFP* dN = j + (deviceFP*)(*s).workspace1;
		deviceFP* E = &(*s).gridETime1[j];
		deviceFP* P = &(*s).gridPolarizationTime1[j];
		for (unsigned int k = 0; k < (*s).Ntime; ++k) {
			N += dN[k];
			integralx += N * (*s).expGammaT[k] * E[k];
			P[k] += expMinusGammaT[k] * integralx;
		}
	};

	trilingual plasmaCurrentKernel_twoStage_B_Cylindric asKernel(withID const deviceParameterSet* s) {
		size_t j = localIndex;
		j *= (*s).Ntime;
		deviceFP N = 0.0;
		deviceFP integral1 = 0.0; 
		deviceFP* expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
		deviceFP* dN = j + (deviceFP*)(*s).workspace1;
		deviceFP* E1 = &(*s).gridETime1[j];
		deviceFP* P1 = &(*s).gridPolarizationTime1[j];

		for (unsigned int k = 0; k < (*s).Ntime; ++k) {
			N += dN[k];
			integral1 += N * (*s).expGammaT[k] * E1[k];
			P1[k] += expMinusGammaT[k] * integral1;
			//expandCylindricalBeamDeviceSingle(s, j + k, (*s).gridRadialLaplacian1, (*s).gridPolarizationTime1);
		}
	};

	trilingual updateKwithPolarizationKernel asKernel(withID const deviceParameterSet* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		h += j * (*sP).Nfreq;
		(*sP).k1[h] += (*sP).gridPolarizationFactor1[h] * (*sP).workspace1[h];
		(*sP).k2[h] += (*sP).gridPolarizationFactor2[h] * (*sP).workspace2P[h];
	};

	trilingual updateKwithPolarizationKernelCylindric asKernel(withID const deviceParameterSet* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		i = h + j * ((*sP).Nfreq);
		h += (j + ((j > ((*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;
		(*sP).k1[i] += (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h];
		(*sP).k2[i] += (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h];
	};

	trilingual updateKwithPlasmaKernel asKernel(withID const deviceParameterSet* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		deviceComplex jfac = deviceComplex(0.0, -1.0 / (h * (*sP).fStep));
		h += j * (*sP).Nfreq;
		(*sP).k1[h] += jfac * (*sP).gridPolarizationFactor1[h] * (*sP).workspace1[h] * (*sP).inverseChiLinear1[h % ((*sP).Nfreq)];
		(*sP).k2[h] += jfac * (*sP).gridPolarizationFactor2[h] * (*sP).workspace2P[h] * (*sP).inverseChiLinear2[h % ((*sP).Nfreq)];
	};

	trilingual updateKwithRadialLaplacianKernel asKernel(withID const deviceParameterSet* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
	};

	trilingual updateKwithPlasmaKernelCylindric asKernel(withID const deviceParameterSet* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		i = h + j * ((*sP).Nfreq);
		deviceComplex jfac = deviceComplex(0.0, -(deviceFP)1.0 / (h * (*sP).fStep));
		h += (j + ( (j > ((*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;
		(*sP).k1[i] += jfac * (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h] * (*sP).inverseChiLinear1[i % ((*sP).Nfreq)];
		(*sP).k2[i] += jfac * (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h] * (*sP).inverseChiLinear2[i % ((*sP).Nfreq)];

		h += 4 * (*sP).NgridC;
		(*sP).k1[i] += (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h];
		(*sP).k2[i] += (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h];
	};

	//Slightly different kernels for the four stages of RK4. They used to be one big kernel with a switch case
	//but this has slightly better utilization.
	trilingual rkKernel0 asKernel(withID const deviceParameterSet* sP) {
		size_t iC = localIndex;
		size_t h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if(iC>(*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] += (*sP).gridPolarizationFactor1[iC] * (*sP).workspace1[iC];
		if (h == 1) (*sP).workspace1[iC - 1] = deviceComplex(0.0, 0.0);
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + (deviceFP)0.5 * (*sP).gridPropagationFactor1[iC] * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridPropagationFactor1[iC] * (*sP).gridPropagationFactor1[iC] * ((deviceFP)SIXTH * (*sP).k1[iC] + (*sP).gridEFrequency1[iC]);
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex(0.0, 0.0);
	};

	trilingual rkKernel1 asKernel(withID const deviceParameterSet* sP) {
		size_t iC = localIndex;
		size_t h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] += (*sP).gridPolarizationFactor1[iC] * (*sP).workspace1[iC];
		if (h == 1)(*sP).workspace1[iC - 1] = deviceComplex(0.0, 0.0);
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + (deviceFP)0.5 * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).gridPropagationFactor1[iC] * (deviceFP)THIRD * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex(0.0, 0.0);
	};

	trilingual rkKernel2 asKernel(withID const deviceParameterSet* sP) {
		size_t iC = localIndex;
		size_t h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] += (*sP).gridPolarizationFactor1[iC] * (*sP).workspace1[iC];
		if (h == 1)(*sP).workspace1[iC - 1] = deviceComplex(0.0, 0.0);
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + (*sP).gridPropagationFactor1[iC] * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).gridPropagationFactor1[iC] * (deviceFP)THIRD * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex(0.0, 0.0);
	};

	trilingual rkKernel3 asKernel(withID const deviceParameterSet* sP) {
		size_t iC = localIndex;
		size_t h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] += (*sP).gridPolarizationFactor1[iC] * (*sP).workspace1[iC];
		if (h == 1)(*sP).workspace1[iC - 1] = deviceComplex(0.0, 0.0);
		(*sP).gridEFrequency1[iC] = (*sP).gridEFrequency1Next1[iC] + (deviceFP)SIXTH * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * (*sP).gridEFrequency1[iC];
		(*sP).k1[iC] = deviceComplex(0.0, 0.0);
	};

	//Kernels for symmetry around z axis use a different form, adding the radial Laplacian
	//instead of the nonlinear polarization
	trilingual rkKernel0Cylindric asKernel(withID const deviceParameterSet* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		if (h == 1) (*sP).workspace1[iC-1] = deviceComplex(0.0, 0.0);
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + (deviceFP)0.5 * (*sP).gridPropagationFactor1[iC] * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridPropagationFactor1[iC] * (*sP).gridPropagationFactor1[iC] * ((deviceFP)SIXTH * (*sP).k1[iC] + (*sP).gridEFrequency1[iC]);
		(*sP).workspace1[iC] =  (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex(0.0, 0.0);
	};

	trilingual rkKernel1Cylindric asKernel(withID const deviceParameterSet* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		if (h == 1)(*sP).workspace1[iC-1] = deviceComplex(0.0, 0.0);
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + (deviceFP)0.5 * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).gridPropagationFactor1[iC] * (deviceFP)THIRD * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex(0.0, 0.0);
	};

	trilingual rkKernel2Cylindric asKernel(withID const deviceParameterSet* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		if (h == 1)(*sP).workspace1[iC-1] = deviceComplex(0.0, 0.0);
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + (*sP).gridPropagationFactor1[iC] * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).gridPropagationFactor1[iC] * (deviceFP)THIRD * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex(0.0, 0.0);
	};

	trilingual rkKernel3Cylindric asKernel(withID const deviceParameterSet* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		if (h == 1)(*sP).workspace1[iC-1] = deviceComplex(0.0, 0.0);
		(*sP).gridEFrequency1[iC] = (*sP).gridEFrequency1Next1[iC] + (deviceFP)SIXTH * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * (*sP).gridEFrequency1[iC];
		(*sP).k1[iC] = deviceComplex(0.0, 0.0);
	};

	trilingual beamNormalizeKernel asKernel(withID const deviceParameterSet* s, deviceFP* rawSum, deviceFP* pulse, deviceFP pulseEnergy) {
		size_t i = localIndex;
		deviceFP normFactor = deviceFPLib::sqrt(pulseEnergy / ((*s).Ntime * (*rawSum)));
		pulse[i] *= normFactor;
	};

	trilingual addDoubleArraysKernel asKernel(withID deviceFP* A, deviceFP* B) {
		size_t i = localIndex;
		A[i] += B[i];
	};

	//crease a pulse on the grid for the 2D modes. Note that normalization of the 2D mode assumes radial symmetry (i.e. that it's
	//a gaussian beam, not an infinite plane wave, which would have zero amplitude for finite energy).
	trilingual beamGenerationKernel2D asKernel(withID deviceComplex* field, devicePulse* p, deviceFP* pulseSum, deviceParameterSet* s,
		bool hasLoadedField, deviceComplex* loadedField, deviceFP* materialPhase, deviceFP* sellmeierCoefficients) {
		long long i = localIndex;
		long long j, h;
		h = 1 + i % ((*s).Nfreq - 1);
		j = i / ((*s).Nfreq - 1);
		i = h + j * ((*s).Nfreq);
		deviceFP f = h * (*s).fStep;
		deviceFP w = TWOPI * (f - (*p).frequency);

		//supergaussian pulse spectrum, if no input pulse specified
		deviceComplex specfac = deviceComplex(-deviceFPLib::pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder), 0);

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
		deviceFP ko = TWOPI * no.real() * f / LIGHTC;
		deviceFP zR = PI * (*p).beamwaist * (*p).beamwaist * no.real() * f / LIGHTC;
		if (f == 0) {
			zR = 1e3;
		}
		deviceFP rB = ((*p).x0 - (*s).dx * (j - (*s).Nspace / 2.0) - 0.25 * (*s).dx);
		deviceFP r = rB * deviceFPLib::cos((*p).beamAngle) - (*p).z0 * deviceFPLib::sin((*p).beamAngle);
		deviceFP z = rB * deviceFPLib::sin((*p).beamAngle) + (*p).z0 * deviceFPLib::cos((*p).beamAngle);

		deviceFP wz = (*p).beamwaist * deviceFPLib::sqrt(1 + (z * z / (zR * zR)));
		deviceFP Rz = z * (1. + (zR * zR / (z * z)));

		if (z == 0) {
			Rz = 1.0e15;
		}
		deviceFP phi = deviceFPLib::atan(z / zR);
		deviceComplex Eb = ((*p).beamwaist / wz) * deviceLib::exp(deviceComplex(0., 1.) * (ko * (z - (*p).z0) + ko * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
		Eb = Eb * specfac;
		if (isnan(cuCModSquared(Eb)) || f <= 0) {
			Eb = deviceComplex(0., 0.);
		}

		field[i] = deviceComplex(deviceFPLib::cos((*p).polarizationAngle), -(*p).circularity * deviceFPLib::sin((*p).polarizationAngle)) * Eb;
		field[i + (*s).NgridC] = deviceComplex(deviceFPLib::sin((*p).polarizationAngle), (*p).circularity * deviceFPLib::cos((*p).polarizationAngle)) * Eb;
		deviceFP pointEnergy = deviceFPLib::abs(r) * (cuCModSquared(field[i]) + cuCModSquared(field[i + (*s).NgridC]));
		pointEnergy *= 2 * PI * LIGHTC * EPS0 * (*s).dx * (*s).dt;
		//two factors of two cancel here - there should be one for the missing frequency plane, but the sum is over x instead of r
		//accordingly we already multiplied by two
		atomicAddDevice(pulseSum, pointEnergy);
	};

	//Generate a beam in full 3D mode
	trilingual beamGenerationKernel3D asKernel(withID deviceComplex* field, devicePulse* p, deviceFP* pulseSum, deviceParameterSet* s,
		bool hasLoadedField, deviceComplex* loadedField, deviceFP* materialPhase, deviceFP* sellmeierCoefficients) {
		long long i = localIndex;
		long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		i = h + col * ((*s).Nfreq);
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;
		deviceFP f = h * (*s).fStep;
		deviceFP w = TWOPI * (f - (*p).frequency);

		//supergaussian pulse spectrum, if no input pulse specified
		deviceComplex specfac = deviceComplex(-deviceFPLib::pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder), 0);

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

		deviceFP ko = TWOPI * no.real() * f / LIGHTC;
		deviceFP zR = PI * (*p).beamwaist * (*p).beamwaist * no.real() * f / LIGHTC;
		if (f == 0) {
			zR = 1e3;
		}
		deviceFP xo = ((*s).dx * (j - (*s).Nspace / 2.0)) - (*p).x0;
		deviceFP yo = ((*s).dx * (k - (*s).Nspace2 / 2.0)) - (*p).y0;
		if (!(*s).is3D) yo = 0.0;
		deviceFP zo = (*p).z0;
		deviceFP cB = deviceFPLib::cos((*p).beamAngle);
		deviceFP cA = deviceFPLib::cos((*p).beamAnglePhi);
		deviceFP sB = deviceFPLib::sin((*p).beamAngle);
		deviceFP sA = deviceFPLib::sin((*p).beamAnglePhi);
		deviceFP x = cB * xo + sA * sB * yo + sA * sB * zo;
		deviceFP y = cA * yo - sA * zo;
		deviceFP z = -sB * xo + sA * cB * yo + cA * cB * zo;
		deviceFP r = sqrt(x * x + y * y);

		deviceFP wz = (*p).beamwaist * sqrt(1 + (z * z / (zR * zR)));
		deviceFP Rz = 1.0e15;
		if (z != 0.0) {
			Rz = z * (1. + (zR * zR / (z * z)));
		}

		deviceFP phi = deviceFPLib::atan(z / zR);
		deviceComplex Eb = ((*p).beamwaist / wz) * deviceLib::exp(deviceComplex(0., 1.) * (ko * (z - (*p).z0) + ko * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
		Eb = Eb * specfac;
		if (isnan(cuCModSquared(Eb)) || f <= 0) {
			Eb = deviceComplex(0., 0.);
		}

		field[i] = deviceComplex(deviceFPLib::cos((*p).polarizationAngle), -(*p).circularity * deviceFPLib::sin((*p).polarizationAngle)) * Eb;
		field[i + (*s).NgridC] = deviceComplex(deviceFPLib::sin((*p).polarizationAngle), (*p).circularity * deviceFPLib::cos((*p).polarizationAngle)) * Eb;
		deviceFP pointEnergy = (cuCModSquared(field[i]) + cuCModSquared(field[i + (*s).NgridC]));
		pointEnergy *= 2 * LIGHTC * EPS0 * (*s).dx * (*s).dx * (*s).dt;

		//factor 2 accounts for the missing negative frequency plane
		atomicAddDevice(pulseSum, pointEnergy);
	};

	trilingual multiplyByConstantKernelD asKernel(withID deviceFP* A, deviceFP val) {
		long long i = localIndex;
		A[i] = val * A[i];
	};

	trilingual multiplyByConstantKernelDZ asKernel(withID deviceComplex* A, deviceFP val) {
		size_t i = localIndex;
		A[i] = val * A[i];

	};

	trilingual multiplicationKernelCompactVector asKernel(withID deviceComplex* A, deviceComplex* B, deviceComplex* C, const deviceParameterSet* s) {
		long long i = localIndex;
		long long h = i % (*s).Nfreq; //temporal coordinate

		C[i] = A[h] * B[i];
	};

	trilingual multiplicationKernelCompactDoubleVector asKernel(withID deviceFP* A, deviceComplex* B, deviceComplex* C, const deviceParameterSet* s) {
		long long i = localIndex;
		long long h = i % (*s).Nfreq; //temporal coordinate

		C[i] = A[h] * B[i];
	};

	trilingual multiplicationKernelCompact asKernel(withID deviceComplex* A, deviceComplex* B, deviceComplex* C) {
		long long i = localIndex;
		C[i] = A[i] * B[i];
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
	trilingual expandCylindricalBeam asKernel(withID const deviceParameterSet* s) {
		long long i = localIndex;
		long long j = i / (*s).Ntime; //spatial coordinate
		long long k = i % (*s).Ntime; //temporal coordinate

		//positions on the expanded grid corresponding the the current index
		long long pos1 = 2 * ((*s).Nspace - j - 1) * (*s).Ntime + k;
		long long pos2 = (2 * j + 1) * (*s).Ntime + k;

		//reuse memory allocated for the radial Laplacian, casting complex double
		//to a 2x larger double real grid
		deviceFP* expandedBeam1 = (deviceFP*)(*s).gridRadialLaplacian1;
		deviceFP* expandedBeam2 = expandedBeam1 + 2 * (*s).Ngrid;

		expandedBeam1[pos1] = (*s).gridPolarizationTime1[i];
		expandedBeam1[pos2] = (*s).gridPolarizationTime1[i];
		expandedBeam2[pos1] = (*s).gridPolarizationTime2[i];
		expandedBeam2[pos2] = (*s).gridPolarizationTime2[i];
	};
}
using namespace kernels;

namespace hostFunctions{
	typedef dlib::matrix<deviceFP, 0, 1> column_vector;
	simulationParameterSet* fittingSet;
	activeDevice* dFit;

	int getTotalSpectrum(activeDevice& d) {
		simulationParameterSet* sCPU = d.cParams;
		deviceParameterSet* sc = d.dParams;

		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.fft((*sc).gridETime1, (*sc).workspace1, deviceFFTD2Z1D);
		if ((*sc).is3D) {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum3DKernel, d.dParamsDevice);
		}
		// else if ((*sc).isCylindric) {
		// 	d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel, d.dParamsDevice);
		// }
		else {
			//uncomment and change logic if I want to use the square spectra
			//d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum2DSquareKernel, d.dParamsDevice);
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel, d.dParamsDevice);
		}
		
		d.deviceMemcpy((double*)(*sCPU).totalSpectrum, (deviceFP*)(*sc).gridPolarizationTime1, 3 * (*sCPU).Nfreq * sizeof(double), DeviceToHost);
		return 0;
	}

	int forwardHankel(activeDevice& d, deviceFP* in, deviceComplex* out) {
		deviceParameterSet* sc = d.dParams;
		d.deviceLaunch((*sc).Nblock, (*sc).Nthread, hankelKernel, d.dParamsDevice, in, (deviceFP*)(*sc).workspace1);
		d.fft((*sc).workspace1, out, deviceFFTD2Z1D);
		return 0;
	}
	int backwardHankel(activeDevice& d, deviceComplex* in, deviceFP* out) {
		deviceParameterSet* sc = d.dParams;
		d.fft(in, (*sc).workspace1, deviceFFTZ2D1D);
		d.deviceLaunch((*sc).Nblock, (*sc).Nthread, inverseHankelKernel, d.dParamsDevice, (deviceFP*)(*sc).workspace1, out);
		return 0;
	}

	int addPulseToFieldArrays(activeDevice& d, pulse& pCPU, bool useLoadedField, std::complex<double>* loadedFieldIn) {

		simulationParameterSet* s = d.cParams;
		deviceParameterSet* sc = d.dParams;
		deviceParameterSet* scDevice = d.dParamsDevice;
		devicePulse* p;
		d.deviceCalloc((void**)&p, 1, sizeof(devicePulse));
		devicePulse devpCPU;
		devpCPU.energy = pCPU.energy;
		devpCPU.frequency = pCPU.frequency;
		devpCPU.bandwidth = pCPU.bandwidth;
		devpCPU.sgOrder = pCPU.sgOrder;
		devpCPU.cep = pCPU.cep;
		devpCPU.delay = pCPU.delay;
		devpCPU.gdd = pCPU.gdd;
		devpCPU.tod = pCPU.tod;
		devpCPU.phaseMaterial = pCPU.phaseMaterial;
		devpCPU.phaseMaterialThickness = pCPU.phaseMaterialThickness;
		devpCPU.beamwaist = pCPU.beamwaist;
		devpCPU.x0 = pCPU.x0;
		devpCPU.y0 = pCPU.y0;
		devpCPU.z0 = pCPU.z0;
		devpCPU.beamAngle = pCPU.beamAngle;
		devpCPU.polarizationAngle = pCPU.polarizationAngle;
		devpCPU.beamAnglePhi = pCPU.beamAnglePhi;
		devpCPU.circularity = pCPU.circularity;
		devpCPU.pulseSum = pCPU.pulseSum;
		d.deviceMemcpy(d.dParamsDevice, sc, sizeof(deviceParameterSet), HostToDevice);

		deviceFP* materialPhase;
		deviceComplex* loadedField;

		d.deviceCalloc((void**)&loadedField, (*sc).Ntime, sizeof(deviceComplex));

		//get the material phase
		deviceFP* materialCoefficients, * sellmeierPropagationMedium;

		d.deviceCalloc((void**)&materialCoefficients, 66, sizeof(deviceFP));
		d.deviceCalloc((void**)&sellmeierPropagationMedium, 66, sizeof(deviceFP));
		d.deviceCalloc((void**)&materialPhase, (*s).Ntime, sizeof(deviceFP));
		d.deviceMemcpy(materialCoefficients, (*s).crystalDatabase[pCPU.phaseMaterial].sellmeierCoefficients, 66 * sizeof(double), HostToDevice);
		d.deviceMemcpy(sellmeierPropagationMedium, (*s).crystalDatabase[(*s).materialIndex].sellmeierCoefficients, 66 * sizeof(double), HostToDevice);
		d.deviceLaunch((unsigned int)(*s).Ntime, 1, materialPhaseKernel, (*s).fStep, (*s).Ntime, materialCoefficients, pCPU.frequency, pCPU.frequency, pCPU.phaseMaterialThickness, pCPU.phaseMaterialThickness, materialPhase, materialPhase);

		deviceFP* pulseSum = &materialCoefficients[0];

		if (useLoadedField) {
			d.deviceMemcpy(loadedField, loadedFieldIn, (*s).Ntime * sizeof(std::complex<double>), HostToDevice);
		}
		d.deviceMemset(pulseSum, 0, sizeof(deviceFP));
		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.deviceMemcpy(p, &devpCPU, sizeof(devicePulse), HostToDevice);
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
		d.deviceLaunch(2 * (*sc).Nblock, (*sc).Nthread, addDoubleArraysKernel, (*sc).gridETime1, (deviceFP*)(*sc).gridPolarizationTime1);

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
				d.deviceMemset((*sc).gridETime1, 0, 2 * (*sc).Ngrid * sizeof(deviceFP));
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
		//d.deviceLaunch((unsigned int)((*sc).NgridC / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplicationKernelCompact, (*sc).gridPropagationFactor1, (*sc).gridEFrequency1Next1, (*sc).k1);

		return 0;
	}

	int applyFresnelLoss(activeDevice& d, simulationParameterSet* s, deviceParameterSet& sc, int materialIndex1, int materialIndex2) {
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).crystalDatabase[materialIndex1].sellmeierCoefficients, 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		deviceFP* sellmeierCoefficients1;
		deviceFP* sellmeierCoefficients2;
		d.deviceCalloc((void**)&sellmeierCoefficients1, 74, sizeof(deviceFP));
		d.deviceCalloc((void**)&sellmeierCoefficients2, 74, sizeof(deviceFP));
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
		d.deviceMemcpy(sc.gridEFrequency1, (*s).EkwOut, 2 * (*s).NgridC * sizeof(std::complex<double>), HostToDevice);

		//transform final result
		d.fft(sc.gridEFrequency1, sc.gridETime1, deviceFFTZ2D);
		d.deviceLaunch(2 * sc.Nblock, sc.Nthread, multiplyByConstantKernelD, sc.gridETime1, 1.0 / sc.Ngrid);
		//copy the field arrays from the GPU to CPU memory
		d.deviceMemcpy((*s).ExtOut, sc.gridETime1, 2 * (*s).Ngrid * sizeof(double), DeviceToHost);
		d.deviceMemcpy((*s).EkwOut, sc.gridEFrequency1, 2 * (*s).Ngrid * sizeof(std::complex<double>), DeviceToHost);

		d.deviceFree(sellmeierCoefficients1);
		d.deviceFree(sellmeierCoefficients2);

		return 0;
	}

	int applyFilter(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, double f0, double bandwidth, double order, double inBandAmplitude, double outOfBandAmplitude) {

		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z);
		deviceParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet), HostToDevice);
		d.deviceLaunch(s.Nblock / 2, s.Nthread, filterKernel, sDevice, 1.0e12*f0, 1.0e12*bandwidth, order, inBandAmplitude, outOfBandAmplitude);

		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(std::complex<double>), DeviceToHost);

		d.fft(s.gridEFrequency1, s.gridETime1, deviceFFTZ2D);

		d.deviceLaunch((int)(s.Ngrid / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);

		getTotalSpectrum(d);

		return 0;
	}

	int applyLorenzian(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, double amplitude, double f0, double gamma, double radius, double order) {

		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z1D);
		deviceParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet), HostToDevice);
		d.deviceLaunch(s.Nblock / 2, s.Nthread, lorentzianSpotKernel, sDevice, amplitude, 1.0e12 * f0, 1.0e12 * gamma, radius, order);
		d.fft(s.gridEFrequency1, s.gridETime1, deviceFFTZ2D1D);
		d.deviceLaunch((int)(s.Ngrid / MIN_GRIDDIM), 2 * MIN_GRIDDIM, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ntime);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(std::complex<double>), DeviceToHost);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);

		getTotalSpectrum(d);

		return 0;
	}

	int applyAperatureFarFieldHankel(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, double diameter, double activationParameter, double xOffset, double yOffset) {
		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		forwardHankel(d, s.gridETime1, s.gridEFrequency1);
		deviceParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet), HostToDevice);
		d.deviceLaunch(s.Nblock / 2, s.Nthread, apertureFarFieldKernelHankel, sDevice, 0.5 * DEG2RAD * diameter, activationParameter, DEG2RAD * xOffset, DEG2RAD * yOffset);
		backwardHankel(d, s.gridEFrequency1, s.gridETime1);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z);
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(std::complex<double>), DeviceToHost);
		getTotalSpectrum(d);
		return 0;
	}

	int applyAperatureFarField(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, double diameter, double activationParameter, double xOffset, double yOffset) {
		if ((*sCPU).isCylindric) {
			return applyAperatureFarFieldHankel(d, sCPU, s, diameter, activationParameter, xOffset, yOffset);
		}
		d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
		d.fft(s.gridETime1, s.gridEFrequency1, deviceFFTD2Z);
		deviceParameterSet* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet), HostToDevice);
		d.deviceLaunch(s.Nblock/2, s.Nthread, apertureFarFieldKernel, sDevice, 0.5 * DEG2RAD * diameter, activationParameter, DEG2RAD * xOffset, DEG2RAD * yOffset);

		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(std::complex<double>), DeviceToHost);


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
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(std::complex<double>), DeviceToHost);
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
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(std::complex<double>), DeviceToHost);
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
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(std::complex<double>), DeviceToHost);
		getTotalSpectrum(d);
		return 0;
	}

	int applyLinearPropagation(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, int materialIndex, double thickness) {

		d.deviceMemcpy(s.gridEFrequency1, (*sCPU).EkwOut, s.NgridC * 2 * sizeof(std::complex<double>), HostToDevice);

		deviceFP* sellmeierCoefficients = (deviceFP*)s.gridEFrequency1Next1;
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
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, s.NgridC * 2 * sizeof(std::complex<double>), DeviceToHost);
		d.fft(s.gridEFrequency1, s.gridETime1, deviceFFTZ2D);
		d.deviceLaunch(2 * s.Nblock, s.Nthread, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);

		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * s.Ngrid * sizeof(double), DeviceToHost);

		return 0;
	}

	int preparePropagationGrids(activeDevice& d) {
		deviceParameterSet* sc = d.dParams;
		simulationParameterSet* s = d.cParams;
		deviceFP* sellmeierCoefficients = (deviceFP*)(*sc).gridEFrequency1Next1;
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
	int rotateField(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s, double rotationAngle) {

		deviceComplex* Ein1 = s.gridEFrequency1;
		deviceComplex* Ein2 = s.gridEFrequency2;
		deviceComplex* Eout1 = s.gridEFrequency1Next1;
		deviceComplex* Eout2 = s.gridEFrequency1Next2;

		//retrieve/rotate the field from the CPU memory
		d.deviceMemcpy(Ein1, (*sCPU).EkwOut, 2 * (*sCPU).NgridC * sizeof(std::complex<double>), HostToDevice);
		d.deviceLaunch((unsigned int)(s.NgridC / MIN_GRIDDIM), MIN_GRIDDIM, rotateFieldKernel, Ein1, Ein2, Eout1, Eout2, rotationAngle);
		d.deviceMemcpy((*sCPU).EkwOut, Eout1, 2 * (*sCPU).NgridC * sizeof(std::complex<double>), DeviceToHost);

		//transform back to time
		d.fft(Eout1, s.gridETime1, deviceFFTZ2D);
		d.deviceLaunch(2 * s.Nblock, s.Nthread, multiplyByConstantKernelD, s.gridETime1, 1.0 / s.Ngrid);
		d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);

		//update spectrum
		getTotalSpectrum(d);
		return 0;
	}

//function to run a RK4 time step
//stepNumber is the sub-step index, from 0 to 3
	int runRK4Step(activeDevice& d, uint8_t stepNumber) {
		deviceParameterSet* sH = d.dParams; 
		deviceParameterSet* sD = d.dParamsDevice;

		// Beam with symmetry around z axis:
		// Nonlinear polarization and plasma use expanded grid
		// Radial laplacian uses standard grid
		// two possible FFT shapes (radial Laplacian always performed)
		if((*sH).isCylindric){
			//ifft to time domain 
			d.fft((*sH).workspace1, (*sH).gridETime1, deviceFFTZ2D);

			//Nonlinear polarization and plasma current are fft-ed in a batch
			//from separate (de-interlaced) time-domain grids.
			//assumption: no plasma without other nonlinearities
			if((*sH).isNonLinear){
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, nonlinearPolarizationKernel, sD);
				if((*sH).hasPlasma){
					d.deviceLaunch((*sH).Nblock, (*sH).Nthread, plasmaCurrentKernel_twoStage_A, sD);
					d.deviceLaunch((unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / MIN_GRIDDIM), 2 * MIN_GRIDDIM, plasmaCurrentKernel_twoStage_B, sD);
					d.deviceLaunch((*sH).Nblock, (*sH).Nthread, expandCylindricalBeam, sD);
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFTD2ZPolarization);
					d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPlasmaKernelCylindric, sD);
				}
				else{
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFTD2ZPolarization);
					d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPolarizationKernelCylindric, sD);
				}
			}
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, radialLaplacianKernel, sD);
			d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFTD2Z);
		}
		// 2D and 3D cartesian
		// Only one type of FFT
		// currently nonlinear polarization and plasma ffts are not batched, could give
		// minor speed boost by combining them, but requires additional memory, so
		// nonlinear polarization and plasma are fft-ed separately to accommodate larger
		// systems.
		else if ((*sH).isNonLinear) {
			//perform inverse FFT to get time-space electric field
			d.fft((*sH).workspace1, (*sH).gridETime1, deviceFFTZ2D);
			//Plasma/multiphoton absorption
			if ((*sH).hasPlasma) {
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, plasmaCurrentKernel_twoStage_A, sD);
				d.deviceLaunch((unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / MIN_GRIDDIM), 2 * MIN_GRIDDIM, plasmaCurrentKernel_twoStage_B, sD);
				d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFTD2Z);
				d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPlasmaKernel, sD);
			}
			//Nonlinear polarization
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, nonlinearPolarizationKernel, sD);
			d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFTD2Z);
		}

		//advance an RK4 step
		switch (stepNumber + 4*(*sH).isCylindric) {
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
			break;
		case 4:
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, rkKernel0Cylindric, sD);
			break;
		case 5:
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, rkKernel1Cylindric, sD);
			break;
		case 6:
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, rkKernel2Cylindric, sD);
			break;
		case 7:
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, rkKernel3Cylindric, sD);
			break;
		}
		return 0;
	}



	unsigned long int solveNonlinearWaveEquationWithDevice(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s) {
		if ((d.hasPlasma != s.hasPlasma) && s.isCylindric) {
			d.deallocateSet(&s);
			d.allocateSet(sCPU, &s);
		}
		//prepare the propagation arrays
		preparePropagationGrids(d);
		prepareElectricFieldArrays(d);
		d.deviceMemcpy(d.dParamsDevice, &s, sizeof(deviceParameterSet), HostToDevice);
		deviceFP* canaryPointer = &s.gridETime1[s.Ntime / 2 + s.Ntime * (s.Nspace / 2 + s.Nspace * (s.Nspace2 / 2))];

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
				d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * (*sCPU).NgridC * sizeof(std::complex<double>), DeviceToHost);
				(*sCPU).statusFlags[0] = 0;
			}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
		}
		if ((*sCPU).isInFittingMode && !(*sCPU).isInSequence)(*(*sCPU).progressCounter)++;

		//take final spectra and transfer the results to the CPU
		d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(std::complex<double>), DeviceToHost);
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
			d.reset(sCPU, &s);
			rotateField(d, sCPU, s, DEG2RAD * parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("set"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			vBlock[(int)parameters[0]] = parameters[1];
			break;
		case funHash("plasmaReinject"):
			(*sCPU).isReinjecting = TRUE;
		case funHash("plasma"):
		{
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
			d.reset(sCPU, &s);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU, s);
			(*sCPU).isFollowerInSequence = TRUE;
		}
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
		case funHash("save"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			{
				size_t saveLoc = (size_t)parameters[0];
				if (saveLoc < (*sCPU).Nsims && saveLoc != 0 && (*sCPU).runType != -1) {
					memcpy(&(*sCPU).ExtOut[saveLoc * (*sCPU).Ngrid * 2], (*sCPU).ExtOut, 2 * (*sCPU).Ngrid * sizeof(double));
					memcpy(&(*sCPU).EkwOut[saveLoc * (*sCPU).NgridC * 2], (*sCPU).EkwOut, 2 * (*sCPU).NgridC * sizeof(std::complex<double>));
					memcpy(&(*sCPU).totalSpectrum[saveLoc * 3 * (*sCPU).Nfreq], (*sCPU).totalSpectrum, 3 * (*sCPU).Nfreq * sizeof(double));
				}
			}
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
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
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
				if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			}

			break;
		case funHash("fresnelLoss"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = DEG2RAD * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = DEG2RAD * parameters[2];
			d.reset(sCPU, &s);
			applyFresnelLoss(d, sCPU, s,
				(int)parameters[4],
				(int)parameters[5]);
			break;
		case funHash("sphericalMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applySphericalMirror(d, sCPU, s, parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("parabolicMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applyParabolicMirror(d, sCPU, s, parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("aperture"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applyAperature(d, sCPU, s,
				parameters[0],
				parameters[1]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("farFieldAperture"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applyAperatureFarField(d, sCPU, s,
				parameters[0],
				parameters[1],
				parameters[2],
				parameters[3]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("energy"):
			{
			if ((*sCPU).runType == -1) break;
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			int targetVar = (int)parameters[0];
			int spectrumType = (int)parameters[1];
			double energy = 0.0;
			for(int i = 0; i < (*sCPU).Nfreq; i++){
				energy += (*sCPU).totalSpectrum[i + (*sCPU).Nfreq * spectrumType];
			}
			energy *= (*sCPU).fStep;
			vBlock[targetVar] = energy;
			}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
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
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("lorentzian"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			applyLorenzian(d, sCPU, s, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("addPulse"):
			if ((*sCPU).runType == -1) break;
		{
			interpretParameters(cc, 21, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU, &s);
			d.deviceMemcpy(s.gridETime1, (*sCPU).ExtOut, 2 * s.Ngrid * sizeof(double), HostToDevice);
			d.deviceMemcpy(s.gridEFrequency1, (*sCPU).EkwOut, 2 * s.NgridC * sizeof(std::complex<double>), HostToDevice);

			pulse p;
			memcpy(&p, &(sCPU->pulse1), sizeof(pulse));
			p.energy = parameters[0];
			p.frequency = 1e12 * parameters[1];
			p.bandwidth = 1e12 * parameters[2];
			p.sgOrder = (int)parameters[3];
			p.cep = parameters[4] * PI;
			p.delay = 1e-15 * parameters[5];
			p.gdd = 1e-30 * parameters[6];
			p.tod = 1e-45 * parameters[7];
			p.phaseMaterial = (int)parameters[8];
			p.phaseMaterialThickness = 1e-6 * parameters[9];
			p.beamwaist = 1e-6 * parameters[10];
			p.x0 = 1e-6 * parameters[11];
			p.y0 = 1e-6 * parameters[12];
			p.z0 = 1e-6 *parameters[13];
			p.beamAngle = DEG2RAD * parameters[14];
			p.beamAnglePhi = DEG2RAD * parameters[15];
			p.polarizationAngle = DEG2RAD * parameters[16];
			p.circularity = parameters[17];
			(*sCPU).materialIndex = (int)parameters[18];
			(*sCPU).crystalTheta = DEG2RAD * parameters[19];
			(*sCPU).crystalPhi = DEG2RAD * parameters[20];

			addPulseToFieldArrays(d, p, FALSE, NULL);
			d.deviceMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * s.NgridC * sizeof(std::complex<double>), DeviceToHost);
			d.deviceMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), DeviceToHost);

			getTotalSpectrum(d);
		}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("for"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			int counter = (int)parameters[0];
			int targetVar = (int)parameters[1];
			std::string currentString = cc.substr(cc.find_first_of('{')+1,std::string::npos);
			std::string forStartString = currentString;
			vBlock[targetVar] = 0.0;
			for (int i = 0; i < counter; i++) {
				while (currentString.length() > 0 && currentString.at(0) != '}'){
					if (currentString.at(0) == '<'){
						currentString = currentString.substr(currentString.find_first_of('>'), std::string::npos);
						if(currentString.length()>0) currentString = currentString.substr(1, std::string::npos);
					}
					if (currentString.at(0) == '{') {
						currentString = currentString.substr(1,std::string::npos);
						while(currentString.find_first_of('{') != std::string::npos 
							&& currentString.find_first_of('{') < currentString.find_first_of('}')){
							currentString = currentString.substr(currentString.find_first_of('}'),std::string::npos);
							currentString = currentString.substr(1, std::string::npos);
						}
						currentString = currentString.substr(currentString.find_first_of('}'),std::string::npos);
						if(currentString.length()<5) break; 
						currentString = currentString.substr(1,std::string::npos);
					}
					interpretCommand(currentString, iBlock, vBlock, d, sCPU, s);
					currentString = currentString.substr(currentString.find_first_of(')') + 1, std::string::npos);
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

	int solveSequenceWithDevice(activeDevice& d, simulationParameterSet* sCPU, deviceParameterSet& s) {
		int error = 0;
		//pointers to where the various parameters are in the struct
		double* targets[38] = { 0,
			&(*sCPU).pulse1.energy, &(*sCPU).pulse2.energy, &(*sCPU).pulse1.frequency, &(*sCPU).pulse2.frequency,
			&(*sCPU).pulse1.bandwidth, &(*sCPU).pulse2.bandwidth, &(*sCPU).pulse1.cep, &(*sCPU).pulse2.cep,
			&(*sCPU).pulse1.delay, &(*sCPU).pulse2.delay, &(*sCPU).pulse1.gdd, &(*sCPU).pulse2.gdd,
			&(*sCPU).pulse1.tod, &(*sCPU).pulse2.tod, &(*sCPU).pulse1.phaseMaterialThickness, &(*sCPU).pulse2.phaseMaterialThickness,
			&(*sCPU).pulse1.beamwaist, &(*sCPU).pulse2.beamwaist,
			&(*sCPU).pulse1.x0, &(*sCPU).pulse2.x0, &(*sCPU).pulse1.z0, &(*sCPU).pulse2.z0,
			&(*sCPU).pulse1.beamAngle, &(*sCPU).pulse2.beamAngle, &(*sCPU).pulse1.polarizationAngle, &(*sCPU).pulse2.polarizationAngle,
			&(*sCPU).pulse1.circularity, &(*sCPU).pulse2.circularity, &(*sCPU).crystalTheta, &(*sCPU).crystalPhi,
			&(*sCPU).nonlinearAbsorptionStrength, &(*sCPU).drudeGamma, &(*sCPU).effectiveMass, &(*sCPU).crystalThickness,
			&(*sCPU).propagationStep, &(*sCPU).i37, &(*sCPU).i37 };

		//unit multipliers from interface units to SI base units.
		double multipliers[38] = { 0,
		1, 1, 1e12, 1e12,
		1e12, 1e12, PI, PI,
		1e-15, 1e-15, 1e-30, 1e-30,
		1e-45, 1e-45, 1e-6, 1e-6,
		1e-6, 1e-6,
		1e-6, 1e-6, 1e-6, 1e-6,
		DEG2RAD, DEG2RAD, DEG2RAD, DEG2RAD,
		1, 1, DEG2RAD, DEG2RAD,
		1, 1e12, 1, 1e-6,
		1e-9, 1, 1 };

		//if it starts with 0, it's an old sequence; quit
		if ((*sCPU).sequenceString[0] == '0') {
			return 15;
		}

		//main text interpreter
		simulationParameterSet sCPUbackupValues;
		simulationParameterSet* sCPUbackup = &sCPUbackupValues;
		memcpy(sCPUbackup, sCPU, sizeof(simulationParameterSet));

		double iBlock[100] = { 0.0 };

		for (int k = 1; k < 38; k++) {
			iBlock[k] = *(targets[k]) / multipliers[k];
		}

		double vBlock[100] = { 0.0 };
		std::string currentString((*sCPU).sequenceString);
		//shortest command is either for() or init(), if there's only 4 characters left, it can only
		//be whitespace or other trailing symbols
		size_t minLength = 5;
		while (currentString.length() > minLength) {
			//skip curly braces (for loops should have been handled by interpretCommand() already)
			if (currentString.at(0) == '{') {
				currentString = currentString.substr(1, std::string::npos);
				while (currentString.find_first_of('{') != std::string::npos
					&& currentString.find_first_of('{') < currentString.find_first_of('}')) {
					currentString = currentString.substr(currentString.find_first_of('}'), std::string::npos);
					currentString = currentString.substr(1, std::string::npos);
				}
				currentString = currentString.substr(currentString.find_first_of('}'), std::string::npos);
				if (currentString.length() < minLength) break;
				currentString = currentString.substr(1, std::string::npos);
			}
			//skip angle brackets (comments)
			if (currentString.at(0) == '<') {
				currentString = currentString.substr(currentString.find_first_of('>'), std::string::npos);
				if (currentString.length() < minLength) break;
				currentString = currentString.substr(1, std::string::npos);
			}

			error = interpretCommand(currentString, iBlock, vBlock, d, sCPU, s);
			if (error) break;
			currentString = currentString.substr(currentString.find_first_of(')'), std::string::npos);

			if (currentString.length() < minLength) break;

			currentString = currentString.substr(1, std::string::npos);

			(*sCPUbackup).isFollowerInSequence = (*sCPU).isFollowerInSequence;
			memcpy(sCPU, sCPUbackup, sizeof(simulationParameterSet));
		}
		return error;
	}

	// helper function for fitting mode, runs the simulation and returns difference from the desired outcome.
	double getResidual(const dlib::matrix<double, 0, 1>& x) {

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

		activeDevice& d = *dFit;
		deviceParameterSet& s = *d.dParams;
		d.cParams = fittingSet;
		d.reset(fittingSet, &s);

		if ((*fittingSet).isInSequence) {
			(*fittingSet).isFollowerInSequence = FALSE;
			solveSequenceWithDevice(d, fittingSet, s);
			(*(*fittingSet).progressCounter)++;
		}
		else {
			solveNonlinearWaveEquationWithDevice(d, fittingSet, s);
		}

		//maximize total spectrum in ROI
		if ((*fittingSet).fittingMode < 3) {
			for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
				result += (*fittingSet).totalSpectrum[(*fittingSet).fittingMode * (*fittingSet).Nfreq + (*fittingSet).fittingROIstart + i];
			}
			return result;
		}

		//mode 3 & 4: match total spectrum to reference given in ascii file
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

		if ((*fittingSet).fittingMode == 4) {
			for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
				a = log10(refSpec[i] / maxRef) - log10(simSpec[i] / maxSim);
				if (!isnan(a)) result += a * a;
			}
		}
		else {
			for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
				a = (refSpec[i] / maxRef) - (simSpec[i] / maxSim);
				result += a * a;
			}
		}


		return sqrt(result);
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
	if ((*sCPU).batchIndex == 36 && (*sCPU).batchLoc1 != 0) return 0;
	deviceParameterSet s;
	memset(&s, 0, sizeof(deviceParameterSet));
	activeDevice d(sCPU, &s);

	unsigned long returnValue = solveSequenceWithDevice(d, sCPU, s);

	d.deallocateSet(&s);
	return returnValue;
}

//run in fitting mode
unsigned long runDlibFittingX(simulationParameterSet* sCPU) {
	simulationParameterSet sCPUcurrent;
	fittingSet = &sCPUcurrent;

	memcpy(fittingSet, sCPU, sizeof(simulationParameterSet));
	deviceParameterSet s;
	memset(&s, 0, sizeof(deviceParameterSet));
	activeDevice d(fittingSet, &s);
	dFit = &d;
	dlib::matrix<double, 0, 1> parameters;
	parameters.set_size((*sCPU).Nfitting);
	dlib::matrix<double, 0, 1> lowerBounds;
	lowerBounds.set_size((*sCPU).Nfitting);
	dlib::matrix<double, 0, 1> upperBounds;
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
	d.cParams = sCPU;
	d.reset(sCPU, &s);
	if ((*sCPU).isInSequence) {
		solveSequenceWithDevice(d, sCPU, s);
	}
	else {
		solveNonlinearWaveEquationWithDevice(d, sCPU, s);
	}

	(*sCPU).progressCounter = originalCounter;

	d.deallocateSet(&s);
	dFit = NULL;
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