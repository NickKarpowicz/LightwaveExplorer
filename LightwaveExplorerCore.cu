#include "LightwaveExplorerDevices/LightwaveExplorerTrilingual.h"
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

namespace deviceFunctions {
	//Expand the information contained in the radially-symmetric beam in the offset grid
	// representation.
	// see the expandCylindricalBeam() kernel for more details
	template<typename T, typename U>
	deviceFunction static void expandCylindricalBeamDevice(const deviceParameterSet<T,U>* s, long long i, T* expandedBeam1, const T* sourceBeam1, const T* sourceBeam2) {
		long long j = i / (*s).Ntime; //spatial coordinate
		long long k = i % (*s).Ntime; //temporal coordinate

		//positions on the expanded grid corresponding the the current index
		long long pos1 = 2 * ((*s).Nspace - j - 1) * (*s).Ntime + k;
		long long pos2 = (2 * j + 1) * (*s).Ntime + k;
		T* expandedBeam2 = expandedBeam1 + 2 * (*s).Ngrid;
		expandedBeam1[pos1] = sourceBeam1[i];
		expandedBeam1[pos2] = sourceBeam1[i];
		expandedBeam2[pos1] = sourceBeam2[i];
		expandedBeam2[pos2] = sourceBeam2[i];
	}
	
	//Calculate the fourier optics propagator (e^ik_z*d) for a given set of values of the maknitude of k, transverse k (dk1, dk2)
	//a reference k0 which defines the speed of the moving frame, and distance d over which to propagate
	template<typename real_t, typename complex_t>
	deviceFunction static complex_t fourierPropagator(complex_t k, real_t dk1, real_t dk2, real_t k0, real_t d) {
		if (deviceFPLib::abs(dk2) < 0.1f * k.real() && deviceFPLib::abs(dk1) < 0.1f *  k.real()) {
			return deviceLib::exp(complex_t(0.0,-d)*((k - k0) - (dk1 * dk1) / (2.0f * k.real()) - (dk2 * dk2) / (2.0f * k.real())));
		}
		complex_t kz = (deviceLib::sqrt(-dk2 * dk2 / (k + deviceFPLib::abs(dk1)) + k - deviceFPLib::abs(dk1)) * deviceLib::sqrt(k + deviceFPLib::abs(dk1)) - k0);
		if (kz.imag() > 0.0f) kz = complex_t(kz.real(), -kz.imag());
		return deviceLib::exp(complex_t(0.0f, -d) * kz);
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
	template<typename deviceFP, typename deviceComplex>
	deviceFunction static deviceComplex sellmeierFunc(deviceFP ls, deviceFP omega,const deviceFP* a, int eqn) {
		deviceFP omega2 = omega * omega;
		deviceFP realPart;
		deviceComplex compPart;
		switch (eqn) {
		case 0:
			[[unlikely]]
			if(ls == -a[3] || ls == -a[6] || ls == -a[9] || ls == -a[12]) return deviceComplex(0.0f,0.0f);
			realPart = a[0]
				+ (a[1] + a[2] * ls) / (ls + a[3])
				+ (a[4] + a[5] * ls) / (ls + a[6])
				+ (a[7] + a[8] * ls) / (ls + a[9])
				+ (a[10] + a[11] * ls) / (ls + a[12])
				+ a[13] * ls
				+ a[14] * ls * ls
				+ a[15] * ls * ls * ls;
			compPart = kLorentzian<deviceFP>()*a[16] / deviceComplex(a[17] - omega2, a[18] * omega)
				+ kLorentzian<deviceFP>()*a[19] / deviceComplex(a[20] - omega2, a[21] * omega);
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
			compPart *= kLorentzian<deviceFP>();
			return deviceLib::sqrt(1.0f + compPart);
		case 2:
		{
			//Up to 7 complex Gaussian functions
			//the real part is the Hilbert transform of the Gaussian, th
			deviceFP scaledF;
			compPart = deviceComplex(a[0], 0.0f);
			for (int i = 0; i < 7; ++i) {
				if (a[3 * i + 1] != 0.0f){
					scaledF = (omega - a[1 + 3 * i]) / (deviceFPLib::sqrt(2.0f) * a[2 + 3 * i]);
					compPart += a[3 + 3 * i] * deviceComplex(-invSqrtPi<deviceFP>() *	(scaledF), -deviceFPLib::exp(-scaledF * scaledF));
				}

			}
			//always select branch with imaginary part < 0
			return deviceComplex((deviceLib::sqrt(compPart)).real(), -deviceFPLib::abs((deviceLib::sqrt(compPart)).imag()));
		}
		case 100:
			[[unlikely]]if(ls == -a[3] || ls == -a[6] || ls == -a[9] || ls == -a[12]) return deviceComplex(0.0f,0.0f);
			realPart = a[0]
				+ (a[1] + a[2] * ls) / (ls + a[3])
				+ (a[4] + a[5] * ls) / (ls + a[6])
				+ (a[7] + a[8] * ls) / (ls + a[9])
				+ (a[10] + a[11] * ls) / (ls + a[12])
				+ a[13] * ls
				+ a[14] * ls * ls
				+ a[15] * ls * ls * ls;
			//"real-valued equation has no business being < 1" - causality
			return deviceComplex(deviceFPLib::sqrt(maxN(realPart, 0.9f)), 0.0f);
		}
		return deviceComplex(1.0f, 0.0f);
	};

	//Sellmeier equation for refractive indicies
	template<typename deviceFP, typename deviceComplex>
	deviceFunction static deviceComplex sellmeierCuda(
		deviceComplex* ne, deviceComplex* no, const deviceFP* a, deviceFP f, deviceFP theta, deviceFP phi, int type, int eqn) {
		if (f == 0.0f) {
			*ne = deviceComplex(1.0f, 0.0f); 
			*no = deviceComplex(1.0f, 0.0f); 
			return deviceComplex(1.0f, 0.0f);
		} //exit immediately for f=0
		deviceFP ls = 2.99792458e14f / f; //wavelength in microns
		ls *= ls; //only wavelength^2 is ever used
		deviceFP omega = twoPi<deviceFP>() * maxN(f,-f);

		//option 0: isotropic
		if (type == 0) {
			*ne = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, a, eqn);
			*no = *ne;
			return *ne;
		}
		//option 1: uniaxial
		else if (type == 1) {
			deviceComplex na = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, a, eqn);
			deviceComplex nb = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, &a[22], eqn);

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
			deviceComplex na = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, a, eqn);
			na *= na;
			deviceComplex nb = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, &a[22], eqn);
			nb *= nb;
			deviceComplex nc = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, &a[44], eqn);
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

	deviceFunction static inline deviceFP cuCModSquared(const deviceComplex& a) {
		return a.real() * a.real() + a.imag() * a.imag();
	}

	//provide a list of nearest-3 neighbors for taking spatial derivatives
	// exploiting the fact that the radial grid is offset by 1/4 step from 0
	// this means that midpoints are available on the other side of the origin.
	// returns rho at the given index j
	template<typename deviceFP>
	deviceFunction static deviceFP resolveNeighborsInOffsetRadialSymmetry(
		long long* neighbors, long long N, int j, deviceFP dr, long long Ntime, long long h) {
		if (j < N / 2) {
			neighbors[0] = (N - j - 2) * Ntime + h;
			neighbors[1] = (j + 1) * Ntime + h;
			neighbors[2] = (N - j - 1) * Ntime + h;
			neighbors[3] = (N - j) * Ntime + h;
			neighbors[4] = (j - 1) * Ntime + h;
			neighbors[5] = (N - j + 1) * Ntime + h;
			return -(dr * (j - N / 2) + 0.25f * dr);
		}
		else {
			neighbors[0] = (N - j + 1) * Ntime + h;
			neighbors[1] = (j - 1) * Ntime + h;
			neighbors[2] = (N - j) * Ntime + h;
			neighbors[3] = (N - j - 1) * Ntime + h;
			neighbors[4] = (j + 1) * Ntime + h;
			neighbors[5] = (N - j - 2) * Ntime + h;
			return dr * (j - N / 2) + 0.25f * dr;
		}
	}
	//provide the position rho in cylindric mode; a simplified
	//version of the resolveNeighbors function above for cases where
	//the neighbors aren't required
	template<typename deviceFP>
	deviceFunction static deviceFP rhoInRadialSymmetry(
		long long N, int j, deviceFP dr) {
		if (j < N / 2) {
			return deviceFPLib::abs( - (dr * (j - N / 2) + 0.25f * dr));
		}
		else {
			return deviceFPLib::abs(dr * (j - N / 2) + 0.25f * dr);
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
	template<typename deviceFP, typename deviceComplex>
	deviceFunction static void findBirefringentCrystalIndex(const deviceParameterSet<deviceFP, deviceComplex>* s, const deviceFP* sellmeierCoefficients, long long i, deviceComplex* n1, deviceComplex* n2) {
		unsigned long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;

		deviceFP f = (*s).fStep * h;
		deviceFP kx1 = (lightC<deviceFP>() / (twoPi<deviceFP>() * f)) * (j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace));
		deviceFP ky1 = (lightC<deviceFP>() / (twoPi<deviceFP>() * f)) * (k * (*s).dk2 - (k >= ((*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2));
		//alpha is deviation from crystal Theta (x2 polarizations)
		//beta is deviation from crystal Phi
		//
		deviceComplex n[4][2]{};
		deviceComplex nW = deviceComplex{};
		sellmeierCuda(&n[0][0], &n[0][1], sellmeierCoefficients, f, sellmeierCoefficients[66], sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
		if ((*s).axesNumber == 0) {
			*n1 = n[0][0];
			*n2 = n[0][1];
			return;
		}

		deviceFP gradient[2][2]{};
		deviceFP alpha[2] = { deviceFPLib::asin(kx1 / n[0][0].real()),deviceFPLib::asin(kx1 / n[0][1].real()) };
		deviceFP beta[2] = { deviceFPLib::asin(ky1 / n[0][0].real()), deviceFPLib::asin(ky1 / n[0][1].real()) };

		deviceFP gradientStep = 1.0e-6f;
		deviceFP gradientFactor = 0.5f / gradientStep;
		int it = 0;
		int maxiter = 64;
		deviceFP gradientTol = 1e-3f;
		//emperical testing: 
		// converges to double precision limit in two iterations for BBO
		// converges in 32 iterations in BiBO

		deviceFP errArray[4][2] = { {0.0f} };
		if ((*s).axesNumber == 1) {
			maxiter = 64;
			sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] + gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
			sellmeierCuda(&n[1][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] - gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
			if (isnan(n[0][0].real()) || isnan(n[0][0].imag()) || isnan(n[1][0].real()) || isnan(n[1][0].imag())) {
				*n1 = deviceComplex{};
				*n2 = deviceComplex{};
				return;
			}
			errArray[0][0] = deviceFPLib::sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
			errArray[1][0] = deviceFPLib::sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
			gradient[0][0] = gradientFactor * (errArray[0][0] - errArray[1][0]);

			for (it = 0; it < maxiter; ++it) {
				if (isnan(n[0][0].real()) || isnan(n[0][0].imag()) || isnan(n[1][0].real()) || isnan(n[1][0].imag())) {
					*n1 = deviceComplex{};
					*n2 = deviceComplex{};
					return;
				}
				if (deviceFPLib::abs(gradient[0][0]) > gradientTol) {
					alpha[0] -= 0.5f * (errArray[0][0] + errArray[1][0]) / gradient[0][0];
				} 
				else {
					break;
				}

				sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] + gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&n[1][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] - gradientStep, sellmeierCoefficients[67], (*s).axesNumber, (*s).sellmeierType);
				errArray[0][0] = deviceFPLib::sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
				errArray[1][0] = deviceFPLib::sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
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
			errArray[0][0] = deviceFPLib::sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
			errArray[1][0] = deviceFPLib::sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
			errArray[2][0] = deviceFPLib::sin(beta[0] + gradientStep) * n[2][0].real() - ky1;
			errArray[3][0] = deviceFPLib::sin(beta[0] - gradientStep) * n[3][0].real() - ky1;
			errArray[0][1] = deviceFPLib::sin(alpha[1] + gradientStep) * n[0][1].real() - kx1;
			errArray[1][1] = deviceFPLib::sin(alpha[1] - gradientStep) * n[1][1].real() - kx1;
			errArray[2][1] = deviceFPLib::sin(beta[1] + gradientStep) * n[2][1].real() - ky1;
			errArray[3][1] = deviceFPLib::sin(beta[1] - gradientStep) * n[3][1].real() - ky1;
			gradient[0][0] = gradientFactor * (errArray[0][0] - errArray[1][0]);
			gradient[1][0] = gradientFactor * (errArray[2][0] - errArray[3][0]);
			gradient[0][1] = gradientFactor * (errArray[0][1] - errArray[1][1]);
			gradient[1][1] = gradientFactor * (errArray[2][1] - errArray[3][1]);

			for (it = 0; it < maxiter; ++it) {
				if (isnan(n[0][0].real()) || isnan(n[0][0].imag()) || isnan(n[1][0].real()) || isnan(n[1][0].imag())) {
					*n1 = deviceComplex{};
					*n2 = deviceComplex{};
					return;
				}
				if (deviceFPLib::abs(gradient[0][0]) > 1e-2f) alpha[0] -= 0.25f * (errArray[0][0] + errArray[1][0]) / gradient[0][0];
				if (deviceFPLib::abs(gradient[1][0]) > 1e-2f) beta[0] -= 0.25f * (errArray[2][0] + errArray[3][0]) / gradient[1][0];
				if (deviceFPLib::abs(gradient[0][1]) > 1e-2f) alpha[1] -= 0.25f * (errArray[0][1] + errArray[1][1]) / gradient[0][1];
				if (deviceFPLib::abs(gradient[1][1]) > 1e-2f) beta[1] -= 0.25f * (errArray[2][1] + errArray[3][1]) / gradient[1][1];

				if (maxN(maxN(deviceFPLib::abs(gradient[0][0]), deviceFPLib::abs(gradient[1][0])), maxN(deviceFPLib::abs(gradient[0][1]), deviceFPLib::abs(gradient[1][1]))) < gradientTol) break;
				sellmeierCuda(&n[0][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] + gradientStep, sellmeierCoefficients[67] + beta[0], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&n[1][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0] - gradientStep, sellmeierCoefficients[67] + beta[0], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&n[2][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0], sellmeierCoefficients[67] + beta[0] + gradientStep, (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&n[3][0], &nW, sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[0], sellmeierCoefficients[67] + beta[0] - gradientStep, (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&nW, &n[0][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1] + gradientStep, sellmeierCoefficients[67] + beta[1], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&nW, &n[1][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1] - gradientStep, sellmeierCoefficients[67] + beta[1], (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&nW, &n[2][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1], sellmeierCoefficients[67] + beta[1] + gradientStep, (*s).axesNumber, (*s).sellmeierType);
				sellmeierCuda(&nW, &n[3][1], sellmeierCoefficients, f, sellmeierCoefficients[66] + alpha[1], sellmeierCoefficients[67] + beta[1] - gradientStep, (*s).axesNumber, (*s).sellmeierType);
				errArray[0][0] = deviceFPLib::sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
				errArray[1][0] = deviceFPLib::sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
				errArray[2][0] = deviceFPLib::sin(beta[0] + gradientStep) * n[2][0].real() - ky1;
				errArray[3][0] = deviceFPLib::sin(beta[0] - gradientStep) * n[3][0].real() - ky1;
				errArray[0][1] = deviceFPLib::sin(alpha[1] + gradientStep) * n[0][1].real() - kx1;
				errArray[1][1] = deviceFPLib::sin(alpha[1] - gradientStep) * n[1][1].real() - kx1;
				errArray[2][1] = deviceFPLib::sin(beta[1] + gradientStep) * n[2][1].real() - ky1;
				errArray[3][1] = deviceFPLib::sin(beta[1] - gradientStep) * n[3][1].real() - ky1;
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
	// structure. The kernelLWE() preprocessor definition handles the return type (void for CUDA,
	// auto for SYCL, and the []( syntax to declare a lambda for SYCL. 
	// The function's closing } has to be followed by a ; to have valid syntax in SYCL 
	// (again, it's a lambda, not a function).
	// localIndex will point to the current thread id. I only use 1D data layouts currently; if
	// 2D or higher dimensions are required, this could be modified (probably together with a different
	// preprocessor definition replacing kernelLWE().

	//calculate the total energy spectrum of the beam for the 2D modes. Note that for the 
	//cartesian one, it will be treated as a round beam instead of an infinite plane wave 
	//in the transverse direction. Thus, the 2D Cartesian spectra are approximations.
	kernelLWE(totalSpectrumKernel, const deviceParameterSet<deviceFP, deviceComplex>* s){
	//kernelLWE(totalSpectrumKernel, const deviceParameterSet<deviceFP, deviceComplex>* s) {
		size_t i = localIndex;
		size_t j;
		deviceFP beamCenter1{};
		deviceFP beamCenter2{};
		deviceFP beamTotal1{};
		deviceFP beamTotal2{};
		deviceFP a, x;

		//find beam centers
		if ((*s).isCylindric) {
			beamCenter1 = ((*s).Nspace / 2.0f * (*s).dx) + 0.25f * (*s).dx;
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
			if (beamTotal1 > 0.0f) {
				beamCenter1 /= beamTotal1;
			}
			if (beamTotal2 > 0.0f) {
				beamCenter2 /= beamTotal2;
			}
		}

		//Integrate total beam power, assuming radially-symmetric beam around
		//the center
		beamTotal1 = 0.0f;
		beamTotal2 = 0.0f;
		for (j = 0; j < (*s).Nspace; ++j) {
			x = (*s).dx * j;
			beamTotal1 += deviceFPLib::abs(x - beamCenter1) * cuCModSquared((*s).workspace1[i + j * (*s).Nfreq]);
			beamTotal2 += deviceFPLib::abs(x - beamCenter2) * cuCModSquared((*s).workspace2[i + j * (*s).Nfreq]);
		}
		beamTotal1 *= constProd(vPi<deviceFP>(), 2.0, lightC<deviceFP>(), eps0<deviceFP>()) * (*s).dx * (*s).dt * (*s).dt;
		beamTotal2 *= constProd(vPi<deviceFP>(), 2.0, lightC<deviceFP>(), eps0<deviceFP>()) * (*s).dx * (*s).dt * (*s).dt;

		//put the values into the output spectrum
		(*s).gridPolarizationTime1[i] = beamTotal1;
		(*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
		(*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
	};

	//Calculate the energy spectrum after a 2D propagation assuming that the beam
	//height in the non-resolved direction is == the grid width (i.e. square grid)
	//More quantitative than the mapping to round beams, but rather specific
	// DISABLED IN FAVOR OF ROUND-BEAM APPROXIMATION
	//kernelLWE(totalSpectrum2DSquareKernel, const deviceParameterSet<deviceFP, deviceComplex>* s) {
	//	size_t i = localIndex;
	//	size_t j;

	//	deviceFP beamTotal1 = 0.0f;
	//	deviceFP beamTotal2 = 0.0f;
	//	//Integrate total beam power
	//	beamTotal1 = 0.0f;
	//	beamTotal2 = 0.0f;
	//	for (j = 0; j < (*s).Nspace; ++j) {
	//		beamTotal1 += cuCModSquared((*s).workspace1[i + j * (*s).Nfreq]);
	//		beamTotal2 += cuCModSquared((*s).workspace2[i + j * (*s).Nfreq]);
	//	}
	//	beamTotal1 *= 2.0f * LIGHTC * eps0<deviceFP>() * (*s).dx * (*s).dx * (*s).Nspace * (*s).dt * (*s).dt;
	//	beamTotal2 *= 2.0f * LIGHTC * eps0<deviceFP>() * (*s).dx * (*s).dx * (*s).Nspace * (*s).dt * (*s).dt;

	//	//put the values into the output spectrum
	//	(*s).gridPolarizationTime1[i] = beamTotal1;
	//	(*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
	//	(*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
	//};

	//Calculate the energy spectrum after a 3D propagation
	kernelLWE(totalSpectrum3DKernel, const deviceParameterSet<deviceFP, deviceComplex>* s) {
		size_t i = localIndex;
		size_t j;

		deviceFP beamTotal1{};
		deviceFP beamTotal2{};
		//Integrate total beam power
		for (j = 0; j < (*s).Nspace * (*s).Nspace2; ++j) {
			beamTotal1 += cuCModSquared((*s).workspace1[i + j * (*s).Nfreq]);
			beamTotal2 += cuCModSquared((*s).workspace2[i + j * (*s).Nfreq]);
		}
		beamTotal1 *= constProd(lightC<deviceFP>(), 2.0 * eps0<deviceFP>()) * (*s).dx * (*s).dx * (*s).dt * (*s).dt;
		beamTotal2 *= constProd(lightC<deviceFP>(), 2.0 * eps0<deviceFP>()) * (*s).dx * (*s).dx * (*s).dt * (*s).dt;

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
	kernelLWE(hankelKernel, const deviceParameterSet<deviceFP, deviceComplex>* s, const deviceFP* in, deviceFP* out) {
		size_t i = localIndex;
		size_t col = i / (*s).Ntime; //spatial coordinate
		deviceFP dk = constProd((deviceFP)(1.0/vPi<deviceFP>()), 2.0) / ((*s).dx * (*s).Nspace);
		in += i % (*s).Ntime;
		out[i] = 0.0f;
		out[i + (*s).Ngrid] = 0.0f;
		deviceFP r0;
		deviceFP J0 = 1.0f;
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
	kernelLWE(inverseHankelKernel, const deviceParameterSet<deviceFP, deviceComplex>* s, const deviceFP* in, deviceFP* out) {
		size_t i = localIndex;
		size_t col = i / (*s).Ntime; //spatial coordinate
		deviceFP dk = constProd((deviceFP)(1.0 / vPi<deviceFP>()), 2.0) / ((*s).dx * (*s).Nspace);
		in += i % (*s).Ntime;;
		out[i] = 0.0f;
		out[i + (*s).Ngrid] = 0.0f;
		deviceFP r0 = rhoInRadialSymmetry((*s).Nspace, col, (*s).dx);
		deviceFP J0 = 1.0f;
		deviceFP k0 = col * dk;
		for (size_t k = 0; k < (*s).Nspace; ++k) {
			k0 = k * dk;
			J0 = k0*j0Device(r0 * k0);
			out[i] += J0 * in[k * (*s).Ntime];
			out[i + (*s).Ngrid] += J0 * in[k * (*s).Ntime + (*s).Ngrid];
		}
		out[i] *= 0.5f * dk / ((*s).Ntime);
		out[i + (*s).Ngrid] *= 0.5f * dk / ((*s).Ntime);
	};

	//rotate the field around the propagation axis (basis change)
	kernelLWE(rotateFieldKernel, const deviceComplex* Ein1, const deviceComplex* Ein2, deviceComplex* Eout1,
		deviceComplex* Eout2, deviceFP rotationAngle) {
		long long i = localIndex;
		Eout1[i] = deviceFPLib::cos(rotationAngle) * Ein1[i] - deviceFPLib::sin(rotationAngle) * Ein2[i];
		Eout2[i] = deviceFPLib::sin(rotationAngle) * Ein1[i] + deviceFPLib::cos(rotationAngle) * Ein2[i];
	};

	//calculate the extra term in the Laplacian encountered in cylindrical coordinates (1/r d/drho)
	kernelLWE(radialLaplacianKernel, const deviceParameterSet<deviceFP, deviceComplex>* s) {
		unsigned long long i = localIndex;
		long long j = i / (*s).Ntime; //spatial coordinate
		long long h = i % (*s).Ntime; //temporal coordinate
		long long neighbors[6];

		//zero at edges of grid
		[[unlikely]] if (j<3 || j>((long long)(*s).Nspace - 4)) {
			(*s).gridRadialLaplacian1[i] = 0.0f;
			(*s).gridRadialLaplacian2[i] = 0.0f;
		}
		else {
			deviceFP rho = resolveNeighborsInOffsetRadialSymmetry(neighbors, (*s).Nspace, (int)j, (*s).dx, (*s).Ntime, h);
			rho = -1.0f / rho;
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

	kernelLWE(apertureFarFieldKernel, const deviceParameterSet<deviceFP, deviceComplex>* s, deviceFP radius, deviceFP activationParameter, deviceFP xOffset, deviceFP yOffset) {
		long long i = localIndex;
		long long col, j, k, l;
		deviceComplex cuZero = deviceComplex{};
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;
		k = col % (*s).Nspace;
		l = col / (*s).Nspace;

		//magnitude of k vector
		deviceFP ko = constProd(twoPi<deviceFP>(),1.0/lightC<double>()) * j * (*s).fStep;

		//transverse wavevector being resolved
		deviceFP dk1 = k * (*s).dk1 - (k >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (long long)(*s).Nspace); //frequency grid in x direction
		deviceFP dk2 = 0.0f;
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

		deviceFP a = 1.0f - (1.0f / (1.0f + deviceFPLib::exp(-activationParameter * (r-radius))));

		(*s).gridEFrequency1[i] *= a;
		(*s).gridEFrequency2[i] *= a;
	};

	kernelLWE(apertureFarFieldKernelHankel, const deviceParameterSet<deviceFP, deviceComplex>* s, deviceFP radius, deviceFP activationParameter, deviceFP xOffset, deviceFP yOffset) {
		long long i = localIndex;
		long long col, j, k;
		deviceComplex cuZero = deviceComplex{};
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;
		k = col % (*s).Nspace;

		//magnitude of k vector
		deviceFP ko = constProd(twoPi<deviceFP>(), 1.0 / lightC<double>()) * j * (*s).fStep;

		//transverse wavevector being resolved
		deviceFP dk1 = constProd((deviceFP)2.0,1.0/vPi<double>()) * k / ((*s).dx * (*s).Nspace);; //frequency grid in x direction

		//light that won't go the the farfield is immediately zero
		if (dk1 * dk1 > ko * ko) {
			(*s).gridEFrequency1[i] = cuZero;
			(*s).gridEFrequency2[i] = cuZero;
			return;
		}

		deviceFP theta1 = deviceFPLib::asin(dk1 / ko);

		deviceFP a = 1.0f - (1.0f / (1.0f + deviceFPLib::exp(-activationParameter * (deviceFPLib::abs(theta1) - radius))));

		(*s).gridEFrequency1[i] *= a;
		(*s).gridEFrequency2[i] *= a;
	};


	//apply a spectral filter to the beam (full details in docs)
	kernelLWE(filterKernel, const deviceParameterSet<deviceFP, deviceComplex>* s, deviceFP f0, deviceFP bandwidth, int order, deviceFP inBandAmplitude, deviceFP outOfBandAmplitude) {
		long long i = localIndex;
		long long col, j;
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;

		deviceFP f = ((*s).fStep * j - f0)/bandwidth;
		for (int p = 1; p < order; p++) {
			f *= f;
		}
		deviceFP filterFunction = outOfBandAmplitude + inBandAmplitude*deviceFPLib::exp(-0.5 *f);
		(*s).gridEFrequency1[i] *= filterFunction;
		(*s).gridEFrequency2[i] *= filterFunction;
	};

	//apply a lorentzian gain or loss in a certain cross-section of the beam.
	// amplitude - strength of the copy of the pulse applied
	// f0 - resonance frequency of the Lorentzian (THz)
	// gamma - linewidth (radHz)
	// radius - radius of the spot (m)
	// order - supergaussian order of the spot shape
	kernelLWE(lorentzianSpotKernel, const deviceParameterSet<deviceFP, deviceComplex>* s, deviceFP amplitude, deviceFP f0, deviceFP gamma, deviceFP radius, deviceFP order) {
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

			x = ((*s).dx * (j - (*s).Nspace / 2.0f));
			y = ((*s).dx * (k - (*s).Nspace2 / 2.0f));
			r = deviceFPLib::sqrt(x * x + y * y);
		}
		else {
			h = 1 + i % ((*s).Nfreq - 1);
			j = i / ((*s).Nfreq - 1);
			i = h + j * ((*s).Nfreq);
			f = h * (*s).fStep;
			r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) + 0.25f * (*s).dx);
		}

		deviceFP w0 = twoPi<deviceFP>() * f0;
		deviceFP w = twoPi<deviceFP>() * f;
		deviceComplex lorentzian = gamma * w0 * amplitude / (w0 * w0 - w * w + deviceComplex(0.0f, gamma * w));
		deviceFP spotFactor = r / radius;
		for (int p = 1; p < (int)order; p++) {
			spotFactor *= spotFactor;
		}
		deviceComplex filterFunction = deviceComplex(0.0f, deviceFPLib::exp(-spotFactor)) * lorentzian;
		(*s).gridEFrequency1[i] += filterFunction * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] += filterFunction * (*s).gridEFrequency2[i];
	};

	//Apply a (soft, possibly) aperture
	kernelLWE(apertureKernel, const deviceParameterSet<deviceFP, deviceComplex>* s, deviceFP radius, deviceFP activationParameter) {
		long long i = localIndex;
		long long j, k, col;

		col = i / (*s).Ntime;
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;
		deviceFP r;
		if ((*s).is3D) {
			deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0f));
			deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0f));
			r = deviceFPLib::sqrt(x * x + y * y);
		}
		else {
			r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) + 0.25f * (*s).dx);
		}

		deviceFP a = 1.0f - (1.0f / (1.0f + deviceFPLib::exp(-activationParameter * (r - radius) / (*s).dx)));

		(*s).gridETime1[i] *= a;
		(*s).gridETime2[i] *= a;
	};

	//apply a spatial phase corresponding to a parabolic mirror (on-axis)
	kernelLWE(parabolicMirrorKernel, const deviceParameterSet<deviceFP, deviceComplex>* s, deviceFP focus) {
		long long i = localIndex;
		long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		i = h + col * (*s).Nfreq;
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;

		deviceFP w = twoPi<deviceFP>() * h * (*s).fStep;
		deviceFP r;
		if ((*s).is3D) {
			deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0f));
			deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0f));
			r = deviceFPLib::sqrt(x * x + y * y);
		}
		else {
			r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) + 0.25f * (*s).dx);
		}

		deviceComplex	u = deviceLib::exp(deviceComplex(0.0f,
			w * r * r * (0.5f / focus) / lightC<deviceFP>()));

		(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
	};

	//apply a spatial phase corresponding to a spherical mirror (on axis)
	kernelLWE(sphericalMirrorKernel, const deviceParameterSet<deviceFP, deviceComplex>* s, deviceFP ROC) {
		long long i = localIndex;
		long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		i = h + col * (*s).Nfreq;
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;

		deviceFP w = twoPi<deviceFP>() * h * (*s).fStep;
		deviceFP r;
		if ((*s).is3D) {
			deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0f));
			deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0f));
			r = deviceFPLib::sqrt(x * x + y * y);
		}
		else {
			r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) + 0.25f * (*s).dx);
		}

		bool isNegative = ROC < 0.0f;
		ROC = deviceFPLib::abs(ROC);
		deviceComplex u = deviceComplex{};
		if (r >= ROC) {
			u = deviceComplex{};
		}
		else if (r > 0.5f * ROC) {
			u = deviceLib::exp(deviceComplex(0.0f,
				2.0f * deviceFPLib::pow(-1.0f, isNegative) * w * ROC * 
				((deviceFPLib::sqrt(1.0f - r * r / (ROC * ROC))) - 1.0f) / lightC<deviceFP>()));
		}
		else {
			deviceFP ratio = r / ROC;
			ratio *= ratio;
			u = deviceLib::exp(deviceComplex(0.0f,
				2.0f * deviceFPLib::pow(-1.0f, isNegative) * w * ROC * 
				(-0.5f * ratio - 0.125f * ratio * ratio - 0.0625f * ratio * ratio * ratio) / lightC<deviceFP>()));
		}


		(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
		if (isnan((*s).gridEFrequency1[i].real()) || isnan((*s).gridEFrequency2[i].real()) || isnan((*s).gridEFrequency1[i].imag()) || isnan((*s).gridEFrequency2[i].imag())) {
			(*s).gridEFrequency1[i] = deviceComplex{};
			(*s).gridEFrequency2[i] = deviceComplex{};
		}
	};

	//apply linear propagation through a given medium to the fields
	kernelLWE(applyLinearPropagationKernel, const deviceFP* sellmeierCoefficients, deviceFP thickness, const deviceParameterSet<deviceFP, deviceComplex>* s) {
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
		deviceFP crystalTheta = sellmeierCoefficients[66];
		deviceFP crystalPhi = sellmeierCoefficients[67];

		//frequency being resolved by current thread
		deviceFP f = h * (*s).fStep;
		deviceFP omega = twoPi<deviceFP>() * f;
		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);
		deviceFP dk1 = j * (*s).dk1 - (j >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace);
		deviceFP dk2 = k * (*s).dk2 - (k >= ((long long)(*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2);
		if (!(*s).is3D)dk2 = 0.0f;
		//if ((*s).isCylindric) dk2 = dk1;
		sellmeierCuda(&n0, &n0o, sellmeierCoefficients, (*s).f0,
			crystalTheta, crystalPhi, axesNumber, sellmeierType);
		if (isnan(ne.real()) || isnan(no.real())) {
			ne = deviceComplex(1.0f, 0.0f);
			no = deviceComplex(1.0f, 0.0f);
		}

		deviceComplex ke = ne * omega / lightC<deviceFP>();
		deviceComplex ko = no * omega / lightC<deviceFP>();
		deviceComplex k0 = n0 * omega / lightC<deviceFP>();
		deviceComplex ts = deviceComplex{};
		deviceComplex tp = deviceComplex{};

		ts = fourierPropagator(ke, dk1, dk2, k0.real(), thickness);
		tp = fourierPropagator(ko, dk1, dk2, k0.real(), thickness);
		

		if (isnan(ts.real()) || isnan(ts.imag())) ts = deviceComplex{};
		if (isnan(tp.real()) || isnan(tp.imag())) tp = deviceComplex{};
		(*s).gridEFrequency1[i] = ts * (*s).gridEFrequency1[i];
		(*s).gridEFrequency2[i] = tp * (*s).gridEFrequency2[i];
		if(isnan((*s).gridEFrequency2[i].real()) ||
			isnan((*s).gridEFrequency2[i].imag()) || 
			isnan((*s).gridEFrequency1[i].real()) ||
			isnan((*s).gridEFrequency1[i].imag())) {
			(*s).gridEFrequency1[i] = deviceComplex{};
			(*s).gridEFrequency2[i] = deviceComplex{};

		}
		if (h == 1) {
			(*s).gridEFrequency1[i-1] = deviceComplex{};
			(*s).gridEFrequency2[i - 1] = deviceComplex{};
		}
	};

	//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
	//note that the sellmeier coefficients have extra values appended to the end
	//to give info about the current simulation
	kernelLWE(prepareCartesianGridsKernel, const deviceFP* sellmeierCoefficients, const deviceParameterSet<deviceFP, deviceComplex>* s) {
		long long i = localIndex;
		long long j, k;
		deviceComplex ne, no;
		deviceComplex n0 = (*s).n0;
		deviceComplex cuZero = deviceComplex{};
		j = i / ((*s).Nfreq - 1); //spatial coordinate
		k = 1 + (i % ((*s).Nfreq - 1)); //temporal coordinate
		i = k + j * (*s).Nfreq;
		deviceComplex ii = deviceComplex(0.0f, 1.0f);
		deviceFP kStep = sellmeierCoefficients[70];
		deviceFP fStep = sellmeierCoefficients[71];

		//frequency being resolved by current thread
		deviceFP f = k * fStep;

		//transverse wavevector being resolved
		deviceFP dk = j * kStep - (j >= ((long long)(*s).Nspace / 2)) * (kStep * (*s).Nspace); //frequency grid in transverse direction
		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);

		//if the refractive index was returned weird, then the index isn't valid, so set the propagator to zero for that frequency
		if (minN(ne.real(), no.real()) < 0.9f || isnan(ne.real()) || isnan(no.real()) || isnan(ne.imag()) || isnan(no.imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			return;
		}

		deviceComplex k0 = deviceComplex(twoPi<deviceFP>() * n0.real() * f / lightC<deviceFP>(), 0.0f);
		deviceComplex ke = twoPi<deviceFP>() * ne * f / lightC<deviceFP>();
		deviceComplex ko = twoPi<deviceFP>() * no * f / lightC<deviceFP>();

		deviceComplex chi11 = deviceComplex(1.0f, 0.0f);
		deviceComplex chi12 = deviceComplex(1.0f, 0.0f);
		if ((*s).isUsingMillersRule) {
			chi11 = (*s).chiLinear1[k];
			chi12 = (*s).chiLinear2[k];
		}
		else {
			chi11 = deviceComplex(1.0f, 0.0f);
			chi12 = deviceComplex(1.0f, 0.0f);
		}
		deviceComplex kz1 = deviceLib::sqrt(ke - dk) * deviceLib::sqrt(ke + dk);
		deviceComplex kz2 = deviceLib::sqrt(ko - dk) * deviceLib::sqrt(ko + dk);

		if (kz1.real() > 0.0f && kz2.real() > 0.0f){
			(*s).gridPropagationFactor1[i] = deviceLib::exp(-0.5f * ii * (kz1 - k0) * (*s).h);
			if (isnan(((*s).gridPropagationFactor1[i]).real())) {
				(*s).gridPropagationFactor1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = deviceLib::exp(-0.5f * ii * (kz2 - k0) * (*s).h);
			if (isnan(((*s).gridPropagationFactor2[i]).real())) {
				(*s).gridPropagationFactor2[i] = cuZero;
			}

			(*s).gridPolarizationFactor1[i] = -ii * deviceLib::pow((deviceComplex)(*s).chiLinear1[k] + (deviceFP)1.0f, (deviceFP)0.25f) * chi11 * ((deviceFP)twoPi<deviceFP>() * (deviceFP)twoPi<deviceFP>() * f * f) / ((2.0f * (deviceFP)lightC<deviceFP>() * (deviceFP)lightC<deviceFP>() * kz1)) * (*s).h;
			(*s).gridPolarizationFactor2[i] = -ii * deviceLib::pow((deviceComplex)(*s).chiLinear2[k] + (deviceFP)1.0f, (deviceFP)0.25f) * chi12 * ((deviceFP)twoPi<deviceFP>() * (deviceFP)twoPi<deviceFP>() * f * f) / ((2.0f * (deviceFP)lightC<deviceFP>() * (deviceFP)lightC<deviceFP>() * kz2)) * (*s).h;
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
	kernelLWE(prepare3DGridsKernel, const deviceFP* sellmeierCoefficients, const deviceParameterSet<deviceFP, deviceComplex>* s) {
		long long i = localIndex;
		long long col, j, k, l;
		deviceComplex ne, no;
		deviceComplex n0 = (*s).n0;
		deviceComplex cuZero = deviceComplex{};
		col = i / ((*s).Nfreq - 1); //spatial coordinate
		j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
		i = j + col * (*s).Nfreq;
		k = col % (*s).Nspace;
		l = col / (*s).Nspace;

		deviceComplex ii = deviceComplex(0.0f, 1.0f);

		//frequency being resolved by current thread
		deviceFP f = j * (*s).fStep;

		//transverse wavevector being resolved
		deviceFP dk1 = k * (*s).dk1 - (k >= ((long long)(*s).Nspace / 2)) * ((*s).dk1 * (long long)(*s).Nspace); //frequency grid in x direction
		deviceFP dk2 = l * (*s).dk2 - (l >= ((long long)(*s).Nspace2 / 2)) * ((*s).dk2 * (long long)(*s).Nspace2); //frequency grid in y direction

		findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);
		if (minN(ne.real(), no.real()) < 0.9f || isnan(ne.real()) || isnan(no.real()) || isnan(ne.imag()) || isnan(no.imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			return;
		}

		deviceComplex k0 = twoPi<deviceFP>() * n0 * f / lightC<deviceFP>();
		deviceComplex ke = twoPi<deviceFP>() * ne * f / lightC<deviceFP>();
		deviceComplex ko = (deviceFP)twoPi<deviceFP>() * no * f / lightC<deviceFP>();

		deviceComplex chi11 = deviceComplex(1.0f, 0.0f);
		deviceComplex chi12 = deviceComplex(1.0f, 0.0f);
		if ((*s).isUsingMillersRule) {
			chi11 = (*s).chiLinear1[j];
			chi12 = (*s).chiLinear2[j];
		}
		else {
			chi11 = deviceComplex(1.0f, 0.0f);
			chi12 = deviceComplex(1.0f, 0.0f);
		}

		deviceComplex kz1 = deviceLib::sqrt(ke * ke - dk1 * dk1 - dk2 * dk2);
		deviceComplex kz2 = deviceLib::sqrt(ko * ko - dk1 * dk1 - dk2 * dk2);
		if (kz1.real() > 0.0f && kz2.real() > 0.0f) {
			(*s).gridPropagationFactor1[i] = fourierPropagator(ke, dk1, dk2, k0.real(), 0.5f * (*s).h);
			if (isnan(((*s).gridPropagationFactor1[i].real()))) {
				(*s).gridPropagationFactor1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = fourierPropagator(ko, dk1, dk2, k0.real(), 0.5f * (*s).h);
			if (isnan(((*s).gridPropagationFactor2[i].real()))) {
				(*s).gridPropagationFactor2[i] = cuZero;
			}

			(*s).gridPolarizationFactor1[i] = -ii * deviceLib::pow((deviceComplex)(*s).chiLinear1[j] + 1.0f, 0.25f) * chi11 * ((deviceFP)twoPi<deviceFP>() * (deviceFP)twoPi<deviceFP>() * f * f) / (2.0f * (deviceFP)lightC<deviceFP>() * (deviceFP)lightC<deviceFP>() * kz1) * (*s).h;
			(*s).gridPolarizationFactor2[i] = -ii * deviceLib::pow((deviceComplex)(*s).chiLinear2[j] + 1.0f, 0.25f) * chi12 * ((deviceFP)twoPi<deviceFP>() * (deviceFP)twoPi<deviceFP>() * f * f) / (2.0f * (deviceFP)lightC<deviceFP>() * (deviceFP)lightC<deviceFP>() * kz2) * (*s).h;
		}
		else {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
		}

		if (isnan((*s).gridPropagationFactor1[i].real() + (*s).gridPropagationFactor1[i].imag()) ||
			isnan((*s).gridPropagationFactor2[i].real() + (*s).gridPropagationFactor2[i].imag()) ||
			isnan((*s).gridPolarizationFactor1[i].real() + (*s).gridPolarizationFactor1[i].imag()) ||
			isnan((*s).gridPolarizationFactor2[i].real() + (*s).gridPolarizationFactor2[i].imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
		}
	};

	//prepare the chi(1) arrays that will be needed in the simulation
	kernelLWE(getChiLinearKernel, deviceParameterSet<deviceFP, deviceComplex>* s, const deviceFP* sellmeierCoefficients) {
		long long i = localIndex;
		int axesNumber = (*s).axesNumber;
		int sellmeierType = (*s).sellmeierType;
		deviceFP crystalTheta = sellmeierCoefficients[66];
		deviceFP crystalPhi = sellmeierCoefficients[67];
		deviceFP fStep = sellmeierCoefficients[71];

		deviceComplex ne, no;

		//frequency being resolved by current thread
		deviceFP f = i * fStep;
		
		sellmeierCuda(&ne, &no, sellmeierCoefficients, f, crystalTheta, crystalPhi, axesNumber, sellmeierType);

		(*s).chiLinear1[i] = ne * ne - 1.0f;
		(*s).chiLinear2[i] = no * no - 1.0f;
		if ((*s).chiLinear1[i].real() != 0.0f && (*s).chiLinear2[i].real() != 0.0f) {
			(*s).inverseChiLinear1[i] = 1.0f / (*s).chiLinear1[i].real();
			(*s).inverseChiLinear2[i] = 1.0f / (*s).chiLinear2[i].real();
		}
		else {
			(*s).inverseChiLinear1[i] = 0.0f;
			(*s).inverseChiLinear2[i] = 0.0f;
		}


		(*s).fieldFactor1[i] = 1.0f / deviceFPLib::pow((*s).chiLinear1[i].real() + 1.0f, 0.25f); //account for the effective field strength in the medium (1/n)
		(*s).fieldFactor2[i] = 1.0f / deviceFPLib::pow((*s).chiLinear2[i].real() + 1.0f, 0.25f);
		if ((*s).isUsingMillersRule) {
			(*s).fieldFactor1[i] *= (*s).chiLinear1[i].real();
			(*s).fieldFactor2[i] *= (*s).chiLinear2[i].real();
		}
	
		if (isnan(ne.real()) || isnan(no.real()) || ne.real() < 0.9f || no.real()<0.9f || ne.imag()>0.0f || no.imag() > 0.0f) {
			ne = deviceComplex(1.0f, 0.0f);
			no = ne;
			(*s).fieldFactor1[i] = 0.0f;
			(*s).fieldFactor2[i] = 0.0f;
			(*s).inverseChiLinear1[i] = 0.0f;
			(*s).inverseChiLinear2[i] = 0.0f;
		}

		if (i == 81) {
			deviceComplex n0;
			sellmeierCuda(&n0, &no, sellmeierCoefficients, deviceFPLib::abs((*s).f0), crystalTheta, crystalPhi, axesNumber, sellmeierType);
			(*s).n0 = no;
			(*s).chiLinear1[(*s).Ntime / 2] = deviceComplex(1.0f, 0.0f);
			(*s).chiLinear2[(*s).Ntime / 2] = deviceComplex(1.0f, 0.0f);
			(*s).fieldFactor1[(*s).Ntime / 2] = 0.0f;
			(*s).fieldFactor2[(*s).Ntime / 2] = 0.0f;
			(*s).inverseChiLinear1[(*s).Ntime / 2] = 1.0f / (*s).chiLinear1[i].real();
			(*s).inverseChiLinear2[(*s).Ntime / 2] = 1.0f / (*s).chiLinear2[i].real();
		}

		//apply Miller's rule to nonlinear coefficients
			if (!(*s).isUsingMillersRule || i > 80) {
				return;
			}
			const deviceFP* referenceFrequencies = &sellmeierCoefficients[72];
			deviceFP chi11[7];

			for (int im = (i>17)*3; im < 7; ++im) {
				if (referenceFrequencies[im] == 0.0f) {
					chi11[im] = 100000.0f;
				}
				else {
					sellmeierCuda(&ne, &no, sellmeierCoefficients, referenceFrequencies[im], sellmeierCoefficients[66], sellmeierCoefficients[67], axesNumber, sellmeierType);
					chi11[im] = no.real() * no.real() - 1.0f;
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
	kernelLWE(prepareCylindricGridsKernel, deviceFP* sellmeierCoefficients, deviceParameterSet<deviceFP, deviceComplex>* s) {
		long long i = localIndex;
		long long j, k;
		long long Nspace = (*s).Nspace;
		deviceComplex cuZero = deviceComplex{};
		j = i / ((*s).Nfreq - 1); //spatial coordinate
		k = 1 + i % ((*s).Nfreq - 1); //temporal coordinate
		i = k + j * (*s).Nfreq;
		deviceComplex ii = deviceComplex(0.0f, 1.0f);
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
		if (minN(ne.real(), no.real()) < 0.95f || ne.real() > 6.0f || no.real() > 6.0f || isnan(ne.real()) || isnan(no.real()) || isnan(ne.imag()) || isnan(no.imag())) {
			(*s).gridPropagationFactor1[i] = cuZero;
			(*s).gridPropagationFactor2[i] = cuZero;
			(*s).gridPolarizationFactor1[i] = cuZero;
			(*s).gridPolarizationFactor2[i] = cuZero;
			(*s).gridPropagationFactor1Rho1[i] = cuZero;
			(*s).gridPropagationFactor1Rho2[i] = cuZero;
			return;
		}

		deviceComplex k0 = deviceComplex(twoPi<deviceFP>() * n0.real() * f / lightC<deviceFP>(), 0.0f);
		deviceComplex ke = twoPi<deviceFP>() * ne * f / lightC<deviceFP>();
		deviceComplex ko = twoPi<deviceFP>() * no * f / lightC<deviceFP>();

		deviceComplex chi11 = deviceComplex(1.0f, 0.0f);
		deviceComplex chi12 = deviceComplex(1.0f, 0.0f);
		if ((*s).isUsingMillersRule) {
			chi11 = (*s).chiLinear1[k];
			chi12 = (*s).chiLinear2[k];
		}
		//fine to here
		if ((dk * dk < minN(ke.real() * ke.real() + ke.imag() * ke.imag(), ko.real() * ko.real() + ko.imag() * ko.imag())) 
			&& (*s).fieldFactor1[k] > 0.0f 
			&& (*s).fieldFactor2[k] > 0.0f) {
			(*s).gridPropagationFactor1[i] = deviceLib::exp(0.5f * ii * (ke - k0 - dk * dk / (2.0f * ke.real())) * (*s).h);
			(*s).gridPropagationFactor1Rho1[i] = ii * (*s).h / ((*s).fieldFactor1[k] * 2.0f * ke);
			if (isnan((deviceLib::abs((*s).gridPropagationFactor1Rho1[i]+(*s).gridPropagationFactor1[i])))) {
				(*s).gridPropagationFactor1[i] = cuZero;
				(*s).gridPropagationFactor1Rho1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = deviceLib::exp(0.5f * ii * (ko - k0 - dk * dk / (2.0f * ko.real())) * (*s).h);
			(*s).gridPropagationFactor1Rho2[i] = ii * (*s).h / ((*s).fieldFactor2[k] * 2.0f * ko);
			if (isnan((deviceLib::abs((*s).gridPropagationFactor1Rho2[i]+(*s).gridPropagationFactor2[i])))) {
				(*s).gridPropagationFactor2[i] = cuZero;
				(*s).gridPropagationFactor1Rho2[i] = cuZero;
			}

			//factor of 0.5 comes from deviceFPd grid size in cylindrical symmetry mode after expanding the beam
			(*s).gridPolarizationFactor1[i] = 0.5f * deviceLib::pow((deviceComplex)(*s).chiLinear1[k] + 1.0f, 0.25f) * chi11 * ii * (twoPi<deviceFP>() * f) / (2.0f * ne.real() * lightC<deviceFP>()) * (*s).h;
			(*s).gridPolarizationFactor2[i] = 0.5f * deviceLib::pow((deviceComplex)(*s).chiLinear2[k] + 1.0f, 0.25f) * chi12 * ii * (twoPi<deviceFP>() * f) / (2.0f * no.real() * lightC<deviceFP>()) * (*s).h;
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
		if ((dk * dk < minN(ke.real() * ke.real() + ke.imag() * ke.imag(), ko.real() * ko.real() + ko.imag() * ko.imag())) && (*s).fieldFactor1[k] > 0.0f && (*s).fieldFactor2[k] > 0.0f) {
			(*s).gridPropagationFactor1[i] = deviceLib::exp(0.5f * ii * (ke - k0 - dk * dk / (2.0f * ke.real())) * (*s).h);
			(*s).gridPropagationFactor1Rho1[i] = ii * (*s).h / ((*s).fieldFactor1[k] * 2.0f * ke);
			if (isnan((deviceLib::abs((*s).gridPropagationFactor1Rho1[i]+(*s).gridPropagationFactor1[i])))) {
				(*s).gridPropagationFactor1[i] = cuZero;
				(*s).gridPropagationFactor1Rho1[i] = cuZero;
			}

			(*s).gridPropagationFactor2[i] = deviceLib::exp(0.5f * ii * (ko - k0 - dk * dk / (2.0f * ko.real())) * (*s).h);
			(*s).gridPropagationFactor1Rho2[i] = ii * (*s).h / ((*s).fieldFactor2[k] * 2.0f * ko);
			if (isnan((deviceLib::abs((*s).gridPropagationFactor1Rho2[i]+(*s).gridPropagationFactor2[i])))) {
				(*s).gridPropagationFactor2[i] = cuZero;
				(*s).gridPropagationFactor1Rho2[i] = cuZero;
			}

			//factor of 0.5 comes from deviceFPd grid size in cylindrical symmetry mode after expanding the beam
			(*s).gridPolarizationFactor1[i] = 0.5f * deviceLib::pow((deviceComplex)(*s).chiLinear1[k] + 1.0f, 0.25f) * chi11 * ii * (twoPi<deviceFP>() * f) / (2.0f * ne.real() * lightC<deviceFP>()) * (*s).h;
			(*s).gridPolarizationFactor2[i] = 0.5f * deviceLib::pow((deviceComplex)(*s).chiLinear2[k] + 1.0f, 0.25f) * chi12 * ii * (twoPi<deviceFP>() * f) / (2.0f * no.real() * lightC<deviceFP>()) * (*s).h;
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

	kernelLWE(materialPhaseKernel, deviceFP df, size_t Ntime, deviceFP* a, deviceFP f0,
		deviceFP thickness, deviceFP* phase1) {
		size_t i = localIndex;
		//frequency being resolved by current thread
		deviceFP f = i * df;

		//give phase shift relative to group velocity (approximated 
		// with low-order finite difference) so the pulse doesn't move
		deviceComplex ne, no, no0, n0p, n0m;
		sellmeierCuda<deviceFP,deviceComplex>(&ne, &no, a, f, 0.0f, 0.0f, 0, 0);
		f *= twoPi<deviceFP>();
		sellmeierCuda<deviceFP, deviceComplex>(&ne, &no0, a, f0, 0.0f, 0.0f, 0, 0);
		sellmeierCuda<deviceFP, deviceComplex>(&ne, &n0p, a, f0 + 1.0e11f, 0.0f, 0.0f, 0, 0);
		sellmeierCuda<deviceFP, deviceComplex>(&ne, &n0m, a, f0 - 1.0e11f, 0.0f, 0.0f, 0, 0);
		no0 = no0 + f0 * (n0p - n0m) / 2.0e11f;
		phase1[i] = thickness* f* (no.real() - no0.real()) / lightC<deviceFP>();
	};

	//calculate the nonlinear polarization, after FFT to get the field
	//in the time domain
	kernelLWE(nonlinearPolarizationKernel, const deviceParameterSet<deviceFP, deviceComplex>* s) {
		size_t i = localIndex;
		deviceFP Ex = (*s).gridETime1[i];
		deviceFP Ey = (*s).gridETime2[i];

		if ((*s).nonlinearSwitches[0] == 1 || (*s).nonlinearSwitches[1] == 1){
			
			deviceFP P[3]{};

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
				(*s).gridPolarizationTime1[i] = 0.0f;
				(*s).gridPolarizationTime2[i] = 0.0f;
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
	kernelLWE(plasmaCurrentKernel_twoStage_A, const deviceParameterSet<deviceFP, deviceComplex>* s) {
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

	kernelLWE(plasmaCurrentKernel_twoStage_B, const deviceParameterSet<deviceFP, deviceComplex>* s) {
		size_t j = localIndex;
		j *= (*s).Ntime;
		deviceFP N{};
		deviceFP integralx{};
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

	kernelLWE(updateKwithPolarizationKernelCylindric, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		i = h + j * ((*sP).Nfreq);
		h += (j + ((j > ((*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;
		(*sP).k1[i] += (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h];
		(*sP).k2[i] += (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h];
	};

	kernelLWE(updateKwithPlasmaKernel, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		deviceComplex jfac = deviceComplex(0.0f, -1.0f / (h * (*sP).fStep));
		h += j * (*sP).Nfreq;
		(*sP).k1[h] += jfac * (*sP).gridPolarizationFactor1[h] * (*sP).workspace1[h] * (*sP).inverseChiLinear1[h % ((*sP).Nfreq)];
		(*sP).k2[h] += jfac * (*sP).gridPolarizationFactor2[h] * (*sP).workspace2P[h] * (*sP).inverseChiLinear2[h % ((*sP).Nfreq)];
	};

	kernelLWE(updateKwithPlasmaKernelCylindric, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t i = localIndex;
		size_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
		size_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
		i = h + j * ((*sP).Nfreq);
		deviceComplex jfac = deviceComplex(0.0f, -1.0f / (h * (*sP).fStep));
		h += (j + ( (j > ((*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;
		(*sP).k1[i] += jfac * (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h] * (*sP).inverseChiLinear1[i % ((*sP).Nfreq)];
		(*sP).k2[i] += jfac * (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h] * (*sP).inverseChiLinear2[i % ((*sP).Nfreq)];

		h += 4 * (*sP).NgridC;
		(*sP).k1[i] += (*sP).gridPolarizationFactor1[i] * (*sP).workspace1[h];
		(*sP).k2[i] += (*sP).gridPolarizationFactor2[i] * (*sP).workspace2P[h];
	};

	//Slightly different kernels for the four stages of RK4. They used to be one big kernel with a switch case
	//but this has slightly better utilization.
	kernelLWE(rkKernel0, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t iC = localIndex;
		size_t h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if(iC>(*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] += (*sP).gridPolarizationFactor1[iC] * (*sP).workspace1[iC];
		[[unlikely]] if (h == 1) (*sP).workspace1[iC - 1] = deviceComplex{};
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + 0.5f * (*sP).gridPropagationFactor1[iC] * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridPropagationFactor1[iC] * (*sP).gridPropagationFactor1[iC] * (sixth<deviceFP>() * (*sP).k1[iC] + (*sP).gridEFrequency1[iC]);
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex{};
	};

	kernelLWE(rkKernel1, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t iC = localIndex;
		size_t h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] += (*sP).gridPolarizationFactor1[iC] * (*sP).workspace1[iC];
		[[unlikely]] if (h == 1)(*sP).workspace1[iC - 1] = deviceComplex{};
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + 0.5f * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).gridPropagationFactor1[iC] * (deviceFP)third<deviceFP>() * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex{};
	};

	kernelLWE(rkKernel2, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t iC = localIndex;
		size_t h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] += (*sP).gridPolarizationFactor1[iC] * (*sP).workspace1[iC];
		[[unlikely]] if (h == 1)(*sP).workspace1[iC - 1] = deviceComplex{};
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + (*sP).gridPropagationFactor1[iC] * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).gridPropagationFactor1[iC] * (deviceFP)third<deviceFP>() * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex{};
	};

	kernelLWE(rkKernel3, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t iC = localIndex;
		size_t h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] += (*sP).gridPolarizationFactor1[iC] * (*sP).workspace1[iC];
		[[unlikely]]if (h == 1)(*sP).workspace1[iC - 1] = deviceComplex{};
		(*sP).gridEFrequency1[iC] = (*sP).gridEFrequency1Next1[iC] + sixth<deviceFP>() * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * (*sP).gridEFrequency1[iC];
		(*sP).k1[iC] = deviceComplex{};
	};

	//Kernels for symmetry around z axis use a different form, adding the radial Laplacian
	//instead of the nonlinear polarization
	kernelLWE(rkKernel0Cylindric, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		[[unlikely]] if (h == 1) (*sP).workspace1[iC-1] = deviceComplex{};
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + 0.5f * (*sP).gridPropagationFactor1[iC] * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridPropagationFactor1[iC] * (*sP).gridPropagationFactor1[iC] * (sixth<deviceFP>() * (*sP).k1[iC] + (*sP).gridEFrequency1[iC]);
		(*sP).workspace1[iC] =  (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex{};
	};

	kernelLWE(rkKernel1Cylindric, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		[[unlikely]] if (h == 1)(*sP).workspace1[iC-1] = deviceComplex{};
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + 0.5f * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).gridPropagationFactor1[iC] * third<deviceFP>() * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex{};
	};

	kernelLWE(rkKernel2Cylindric, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		[[unlikely]] if (h == 1)(*sP).workspace1[iC-1] = deviceComplex{};
		deviceComplex estimate1 = (*sP).gridPropagationFactor1[iC] * (*sP).gridPropagationFactor1[iC] * (*sP).gridEFrequency1[iC] + (*sP).gridPropagationFactor1[iC] * (*sP).k1[iC];
		(*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).gridPropagationFactor1[iC] * third<deviceFP>() * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * estimate1;
		(*sP).k1[iC] = deviceComplex{};
	};

	kernelLWE(rkKernel3Cylindric, const deviceParameterSet<deviceFP, deviceComplex>* sP) {
		size_t iC = localIndex;
		unsigned int h = 1 + iC % ((*sP).Nfreq - 1); //frequency coordinate
		iC = h + (iC / ((unsigned int)(*sP).Nfreq - 1)) * ((unsigned int)(*sP).Nfreq);
		deviceFP ff = (*sP).fieldFactor1[h];
		if (iC > (*sP).NgridC)ff = (*sP).fieldFactor2[h];
		(*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1[iC];
		[[unlikely]] if (h == 1)(*sP).workspace1[iC-1] = deviceComplex{};
		(*sP).gridEFrequency1[iC] = (*sP).gridEFrequency1Next1[iC] + sixth<deviceFP>() * (*sP).k1[iC];
		(*sP).workspace1[iC] = (*sP).fftNorm * ff * (*sP).gridEFrequency1[iC];
		(*sP).k1[iC] = deviceComplex{};
	};

	kernelLWE(beamNormalizeKernel, const deviceParameterSet<deviceFP, deviceComplex>* s, const deviceFP* rawSum, deviceFP* field, const deviceFP pulseEnergy) {
		size_t i = localIndex;
		deviceFP normFactor = deviceFPLib::sqrt(pulseEnergy / ((deviceFP)(*s).Ntime * (*rawSum)));
		field[i]  *= normFactor;
	};

	kernelLWE(addDoubleArraysKernel, deviceFP* A, deviceFP* B) {
		size_t i = localIndex;
		A[i] += B[i];
	};

	//crease a pulse on the grid for the 2D modes. Note that normalization of the 2D mode assumes radial symmetry (i.e. that it's
	//a gaussian beam, not an infinite plane wave, which would have zero amplitude for finite energy).
	kernelLWE(beamGenerationKernel2D, deviceComplex* field, pulse<deviceFP>* p, deviceFP* pulseSum, deviceParameterSet<deviceFP, deviceComplex>* s,
		bool hasLoadedField, deviceComplex* loadedField, deviceFP* materialPhase, deviceFP* sellmeierCoefficients) {
		long long i = localIndex;
		long long j, h;
		h = 1 + i % ((*s).Nfreq - 1);
		j = i / ((*s).Nfreq - 1);
		i = h + j * ((*s).Nfreq);
		deviceFP f = h * (*s).fStep;
		deviceFP w = twoPi<deviceFP>() * (f - (*p).frequency);

		//supergaussian pulse spectrum, if no input pulse specified
		deviceComplex specfac = deviceComplex(-deviceFPLib::pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder), 0.0f);

		deviceComplex specphase = deviceComplex(0.0f,
			-((*p).cep
				+ twoPi<deviceFP>() * f * ((*p).delay - 0.5f * (*s).dt * (*s).Ntime)
				+ 0.5f * (*p).gdd * w * w
				+ sixth<deviceFP>()*(*p).tod * w * w * w
				+ materialPhase[h]));
		specfac = deviceLib::exp(specfac + specphase);

		if (hasLoadedField) {
			specfac = loadedField[h] * deviceLib::exp(specphase);
		}
		deviceComplex ne, no;
		sellmeierCuda(&ne, &no, sellmeierCoefficients, f, (*s).crystalTheta, (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);
		deviceFP ko = twoPi<deviceFP>() * no.real() * f / lightC<deviceFP>();
		deviceFP zR = vPi<deviceFP>() * (*p).beamwaist * (*p).beamwaist * no.real() * f / lightC<deviceFP>();
		if (f == 0.0f) {
			zR = 1e3f;
		}
		deviceFP rB = ((*p).x0 - (*s).dx * (j - (*s).Nspace / 2.0f) - 0.25f * (*s).dx);
		deviceFP r = rB * deviceFPLib::cos((*p).beamAngle) - (*p).z0 * deviceFPLib::sin((*p).beamAngle);
		deviceFP z = rB * deviceFPLib::sin((*p).beamAngle) + (*p).z0 * deviceFPLib::cos((*p).beamAngle);

		deviceFP wz = (*p).beamwaist * deviceFPLib::sqrt(1.0f + (z * z / (zR * zR)));
		deviceFP Rz = z * (1.0f + (zR * zR / (z * z)));

		if (z == 0.0f) {
			Rz = 1.0e15f;
		}
		deviceFP phi = deviceFPLib::atan(z / zR);
		deviceComplex Eb = ((*p).beamwaist / wz) * deviceLib::exp(deviceComplex(0.0f, 1.0f) * (ko * (z - (*p).z0) + ko * r * r / (2.0f * Rz) - phi) - r * r / (wz * wz));
		Eb = Eb * specfac;
		if (isnan(cuCModSquared(Eb)) || f <= 0.0f) {
			Eb = deviceComplex{};
		}

		field[i] = deviceComplex(deviceFPLib::cos((*p).polarizationAngle), -(*p).circularity * deviceFPLib::sin((*p).polarizationAngle)) * Eb;
		field[i + (*s).NgridC] = deviceComplex(deviceFPLib::sin((*p).polarizationAngle), (*p).circularity * deviceFPLib::cos((*p).polarizationAngle)) * Eb;
		deviceFP pointEnergy = deviceFPLib::abs(r) * (cuCModSquared(field[i]) + cuCModSquared(field[i + (*s).NgridC]));
		pointEnergy *= 2.0f * vPi<deviceFP>() * lightC<deviceFP>() * eps0<deviceFP>() * (*s).dx * (*s).dt;
		//two factors of two cancel here - there should be one for the missing frequency plane, but the sum is over x instead of r
		//accordingly we already multiplied by two
		atomicAdd(pulseSum, pointEnergy);
	};

	//Generate a beam in full 3D mode
	kernelLWE(beamGenerationKernel3D, deviceComplex* field, pulse<deviceFP>* p, deviceFP* pulseSum, deviceParameterSet<deviceFP, deviceComplex>* s,
		bool hasLoadedField, deviceComplex* loadedField, deviceFP* materialPhase, deviceFP* sellmeierCoefficients) {
		long long i = localIndex;
		long long j, k, h, col;
		h = 1 + i % ((*s).Nfreq - 1);
		col = i / ((*s).Nfreq - 1);
		i = h + col * ((*s).Nfreq);
		j = col % (*s).Nspace;
		k = col / (*s).Nspace;
		deviceFP f = h * (*s).fStep;
		deviceFP w = twoPi<deviceFP>() * (f - (*p).frequency);

		//supergaussian pulse spectrum, if no input pulse specified
		deviceComplex specfac = deviceComplex(-deviceFPLib::pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder), 0.0f);

		deviceComplex specphase = deviceComplex(0.0f,
			-((*p).cep
				+ twoPi<deviceFP>() * f * ((*p).delay - 0.5f * (*s).dt * (*s).Ntime)
				+ 0.5f * (*p).gdd * w * w
				+ (*p).tod * w * w * w / 6.0f
				+ materialPhase[h]));
		specfac = deviceLib::exp(specfac + specphase);
		
		if (hasLoadedField) {
			specfac = loadedField[h] * deviceLib::exp(specphase);
		}
		deviceComplex ne, no;
		sellmeierCuda(&ne, &no, sellmeierCoefficients, f, (*s).crystalTheta, (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);

		deviceFP ko = twoPi<deviceFP>() * no.real() * f / lightC<deviceFP>();
		deviceFP zR = vPi<deviceFP>() * (*p).beamwaist * (*p).beamwaist * no.real() * f / lightC<deviceFP>();
		if (f == 0.0f) {
			zR = 1e3f;
		}
		deviceFP xo = ((*s).dx * (j - (*s).Nspace / 2.0f)) - (*p).x0;
		deviceFP yo = ((*s).dx * (k - (*s).Nspace2 / 2.0f)) - (*p).y0;
		if (!(*s).is3D) yo = 0.0f;
		deviceFP zo = (*p).z0;
		deviceFP cB = deviceFPLib::cos((*p).beamAngle);
		deviceFP cA = deviceFPLib::cos((*p).beamAnglePhi);
		deviceFP sB = deviceFPLib::sin((*p).beamAngle);
		deviceFP sA = deviceFPLib::sin((*p).beamAnglePhi);
		deviceFP x = cB * xo + sA * sB * yo + sA * sB * zo;
		deviceFP y = cA * yo - sA * zo;
		deviceFP z = -sB * xo + sA * cB * yo + cA * cB * zo;
		deviceFP r = deviceFPLib::sqrt(x * x + y * y);

		deviceFP wz = (*p).beamwaist * deviceFPLib::sqrt(1.0f + (z * z / (zR * zR)));
		deviceFP Rz = 1.0e15f;
		if (z != 0.0f) {
			Rz = z * (1.0f + (zR * zR / (z * z)));
		}

		deviceFP phi = deviceFPLib::atan(z / zR);
		deviceComplex Eb = ((*p).beamwaist / wz) * deviceLib::exp(deviceComplex(0.0f, 1.0f) * (ko * (z - (*p).z0) + ko * r * r / (2.0f * Rz) - phi) - r * r / (wz * wz));
		Eb = Eb * specfac;
		if (isnan(cuCModSquared(Eb)) || f <= 0.0f) {
			Eb = deviceComplex{};
		}

		field[i] = deviceComplex(deviceFPLib::cos((*p).polarizationAngle), -(*p).circularity * deviceFPLib::sin((*p).polarizationAngle)) * Eb;
		field[i + (*s).NgridC] = deviceComplex(deviceFPLib::sin((*p).polarizationAngle), (*p).circularity * deviceFPLib::cos((*p).polarizationAngle)) * Eb;
		deviceFP pointEnergy = (cuCModSquared(field[i]) + cuCModSquared(field[i + (*s).NgridC]));
		pointEnergy *= 2.0f * lightC<deviceFP>() * eps0<deviceFP>() * (*s).dx * (*s).dx * (*s).dt;

		//factor 2 accounts for the missing negative frequency plane
		atomicAdd(pulseSum, pointEnergy);
	};

	kernelLWE(multiplyByConstantKernelD, deviceFP* A, deviceFP val) {
		long long i = localIndex;
		A[i] = val * A[i];
	};

	kernelLWE(multiplyByConstantKernelDZ, deviceComplex* A, deviceFP val) {
		size_t i = localIndex;
		A[i] = val * A[i];

	};

	kernelLWE(multiplicationKernelCompactDoubleVector, deviceFP* A, deviceComplex* B, deviceComplex* C, const deviceParameterSet<deviceFP, deviceComplex>* s) {
		long long i = localIndex;
		long long h = i % (*s).Nfreq; //temporal coordinate

		C[i] = A[h] * B[i];
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
	kernelLWE(expandCylindricalBeam, const deviceParameterSet<deviceFP, deviceComplex>* s) {
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
	static simulationParameterSet* fittingSet;
	static ActiveDevice* dFit;

	static int getTotalSpectrum(ActiveDevice& d) {
		simulationParameterSet* sCPU = d.cParams;
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;

		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.fft((*sc).gridETime1, (*sc).workspace1, deviceFFT::D2Z_1D);
		if ((*sc).is3D) {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum3DKernel, d.dParamsDevice);
		}
		else if ((*sc).isCylindric) {
		 	d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel, d.dParamsDevice);
		}
		else {
			//uncomment and change logic if I want to use the square spectra
			//d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum2DSquareKernel, d.dParamsDevice);
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel, d.dParamsDevice);
		}
		
		d.deviceMemcpy((double*)(*sCPU).totalSpectrum, (deviceFP*)(*sc).gridPolarizationTime1, 3 * (*sCPU).Nfreq * sizeof(double), copyType::ToHost);
		return 0;
	}

	static int forwardHankel(ActiveDevice& d, deviceFP* in, deviceComplex* out) {
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		d.deviceLaunch((*sc).Nblock, (*sc).Nthread, hankelKernel, d.dParamsDevice, in, (deviceFP*)(*sc).workspace1);
		d.fft((*sc).workspace1, out, deviceFFT::D2Z_1D);
		return 0;
	}
	static int backwardHankel(ActiveDevice& d, deviceComplex* in, deviceFP* out) {
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		d.fft(in, (*sc).workspace1, deviceFFT::Z2D_1D);
		d.deviceLaunch((*sc).Nblock, (*sc).Nthread, inverseHankelKernel, d.dParamsDevice, (deviceFP*)(*sc).workspace1, out);
		return 0;
	}

	static int addPulseToFieldArrays(ActiveDevice& d, pulse<double>& pCPU, const bool useLoadedField, const std::complex<double>* loadedFieldIn) {

		simulationParameterSet* s = d.cParams;
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		deviceParameterSet<deviceFP, deviceComplex>* scDevice = d.dParamsDevice;
		pulse<deviceFP>* p;
		d.deviceCalloc((void**)&p, 1, sizeof(pulse<deviceFP>));
		pulse<deviceFP> devpCPU = pCPU;

		d.deviceMemcpy(d.dParamsDevice, sc, sizeof(deviceParameterSet<deviceFP, deviceComplex>), copyType::ToDevice);

		deviceFP* materialPhase;
		deviceComplex* loadedField;

		d.deviceCalloc((void**)&loadedField, (*sc).Ntime, sizeof(deviceComplex));

		//get the material phase
		deviceFP* materialCoefficients, * sellmeierPropagationMedium;

		d.deviceCalloc((void**)&materialCoefficients, 66, sizeof(deviceFP));
		d.deviceCalloc((void**)&sellmeierPropagationMedium, 66, sizeof(deviceFP));
		d.deviceCalloc((void**)&materialPhase, (*s).Nfreq, sizeof(deviceFP));
		d.deviceMemcpy(materialCoefficients, (*s).crystalDatabase[pCPU.phaseMaterial].sellmeierCoefficients.data(), 66 * sizeof(double), copyType::ToDevice);
		d.deviceMemcpy(sellmeierPropagationMedium, (*s).crystalDatabase[(*s).materialIndex].sellmeierCoefficients.data(), 66 * sizeof(double), copyType::ToDevice);
		d.deviceLaunch((unsigned int)(*s).Nfreq, 1, materialPhaseKernel, (deviceFP)(*s).fStep, 
			(*s).Ntime, materialCoefficients, (deviceFP)pCPU.frequency, (deviceFP)pCPU.phaseMaterialThickness, 
			materialPhase);

		deviceFP* pulseSum = &materialCoefficients[0];

		if (useLoadedField) {
			d.deviceMemcpy(loadedField, loadedFieldIn, (*s).Ntime * sizeof(std::complex<double>), copyType::ToDevice);
		}
		d.deviceMemset(pulseSum, 0, sizeof(deviceFP));
		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.deviceMemcpy(p, &devpCPU, sizeof(pulse<deviceFP>), copyType::ToDevice);
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

		d.fft((*sc).workspace1, (*sc).gridPolarizationTime1, deviceFFT::Z2D_1D);

		d.deviceLaunch(2 * (*sc).Nblock, (*sc).Nthread, beamNormalizeKernel, scDevice, pulseSum, (*sc).gridPolarizationTime1, (deviceFP)pCPU.energy);

		//add the pulses
		d.deviceLaunch(2 * (*sc).Nblock, (*sc).Nthread, addDoubleArraysKernel, (*sc).gridETime1, (deviceFP*)(*sc).gridPolarizationTime1);

		//fft onto frequency grid
		d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);

		d.deviceFree(materialPhase);
		d.deviceFree(materialCoefficients);
		d.deviceFree(sellmeierPropagationMedium);
		d.deviceFree(loadedField);
		d.deviceFree(p);
		return 0;
	}
	
	static int prepareElectricFieldArrays(ActiveDevice& d) {

		simulationParameterSet* s = d.cParams;
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		
		d.deviceMemcpy(d.dParamsDevice, sc, sizeof(deviceParameterSet<deviceFP, deviceComplex>), copyType::ToDevice);
		deviceParameterSet<deviceFP, deviceComplex>* scDevice = d.dParamsDevice;
		
		if (!(*s).isFollowerInSequence || (*s).isReinjecting) {
			if (!(*s).isReinjecting) {
				d.deviceMemset((*sc).gridETime1, 0, 2 * (*sc).Ngrid * sizeof(deviceFP));
			}
			else {
				d.deviceMemcpy((*sc).gridETime1, (*s).ExtOut, 2 * (*s).Ngrid * sizeof(double), copyType::ToDevice);
				d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);
			}
			
			addPulseToFieldArrays(d, d.cParams->pulse1, d.cParams->field1IsAllocated, d.cParams->loadedField1);
			
			addPulseToFieldArrays(d, d.cParams->pulse2, d.cParams->field2IsAllocated, d.cParams->loadedField2);
		}
		else {
			d.deviceMemcpy((*sc).gridETime1, (*s).ExtOut, 2 * (*s).Ngrid * sizeof(double), copyType::ToDevice);
			d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);
		}
		
		//Copy the field into the temporary array
		d.deviceMemcpy((void*)(*sc).gridEFrequency1Next1, (void*)(*sc).gridEFrequency1, 2 * (*sc).NgridC * sizeof(deviceComplex), copyType::OnDevice);

		//set the propagation grids how they should be at the beginning of the next step
		d.deviceLaunch((unsigned int)((*sc).NgridC / minGridDimension), 2 * minGridDimension, 
			multiplicationKernelCompactDoubleVector, 
			(*sc).fieldFactor1, (*sc).gridEFrequency1Next1, (*sc).workspace1, scDevice);
		d.deviceLaunch((unsigned int)((*sc).NgridC / minGridDimension), 2 * minGridDimension, 
			multiplyByConstantKernelDZ, (*sc).workspace1, (*sc).fftNorm);

		return 0;
	}

	static int applyFresnelLoss(ActiveDevice& d, simulationParameterSet* s, deviceParameterSet<deviceFP, deviceComplex>& sc, int materialIndex1, int materialIndex2) {
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).crystalDatabase[materialIndex1].sellmeierCoefficients.data(), 66 * (sizeof(double)));
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
		d.deviceMemcpy(sellmeierCoefficients1, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), copyType::ToDevice);

		memcpy(sellmeierCoefficientsAugmentedCPU, (*s).crystalDatabase[materialIndex2].sellmeierCoefficients.data(), 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		d.deviceMemcpy(sellmeierCoefficients2, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), copyType::ToDevice);
		d.deviceMemcpy(sc.gridEFrequency1, (*s).EkwOut, 2 * (*s).NgridC * sizeof(std::complex<double>), copyType::ToDevice);

		//transform final result
		d.fft(sc.gridEFrequency1, sc.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(2 * sc.Nblock, sc.Nthread, multiplyByConstantKernelD, 
			sc.gridETime1, (deviceFP)(1.0 / sc.Ngrid));
		//copy the field arrays from the GPU to CPU memory
		d.deviceMemcpy((*s).ExtOut, sc.gridETime1, 2 * (*s).Ngrid * sizeof(double), copyType::ToHost);
		d.deviceMemcpy((*s).EkwOut, sc.gridEFrequency1, 2 * (*s).Ngrid * sizeof(std::complex<double>), copyType::ToHost);

		d.deviceFree(sellmeierCoefficients1);
		d.deviceFree(sellmeierCoefficients2);

		return 0;
	}

	static int applyFilter(ActiveDevice& d, simulationParameterSet* sCPU, double f0, double bandwidth, double order, double inBandAmplitude, double outOfBandAmplitude) {

		d.deviceMemcpy(d.deviceStruct.gridETime1, (*sCPU).ExtOut, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, filterKernel, 
			sDevice, 
			(deviceFP)(1.0e12 * f0),
			(deviceFP)(1.0e12 * bandwidth),
			(int)round(order),
			(deviceFP)inBandAmplitude,
			(deviceFP)outOfBandAmplitude);

		d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToHost);

		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);

		d.deviceLaunch((int)(d.deviceStruct.Ngrid / minGridDimension), 2 * minGridDimension, multiplyByConstantKernelD, 
			d.deviceStruct.gridETime1, (deviceFP)(1.0 / d.deviceStruct.Ngrid));
		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), copyType::ToHost);

		getTotalSpectrum(d);

		return 0;
	}

	static int applyLorenzian(ActiveDevice& d, simulationParameterSet* sCPU, double amplitude, double f0, double gamma, double radius, double order) {

		d.deviceMemcpy(d.deviceStruct.gridETime1, (*sCPU).ExtOut, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, lorentzianSpotKernel, 
			sDevice, 
			(deviceFP)amplitude, 
			(deviceFP)(1.0e12 * f0),
			(deviceFP)(1.0e12 * gamma),
			(deviceFP)radius,
			(deviceFP)order);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.deviceLaunch((int)(d.deviceStruct.Ngrid / minGridDimension), 2 * minGridDimension, multiplyByConstantKernelD, 
			d.deviceStruct.gridETime1, (deviceFP)(1.0 / d.deviceStruct.Ntime));
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToHost);
		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), copyType::ToHost);

		getTotalSpectrum(d);

		return 0;
	}

	static int applyAperatureFarFieldHankel(ActiveDevice& d, simulationParameterSet* sCPU, double diameter, double activationParameter, double xOffset, double yOffset) {
		d.deviceMemcpy(d.deviceStruct.gridETime1, (*sCPU).ExtOut, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToDevice);
		forwardHankel(d, d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, 
			apertureFarFieldKernelHankel, 
			sDevice, 
			(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter), 
			(deviceFP)activationParameter, 
			(deviceFP)(deg2Rad<deviceFP>() * xOffset), 
			(deviceFP)(deg2Rad<deviceFP>() * yOffset));
		backwardHankel(d, d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1);
		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), copyType::ToHost);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyAperatureFarField(ActiveDevice& d, simulationParameterSet* sCPU, double diameter, double activationParameter, double xOffset, double yOffset) {
		if ((*sCPU).isCylindric) {
			return applyAperatureFarFieldHankel(d, sCPU, diameter, activationParameter, xOffset, yOffset);
		}
		d.deviceMemcpy(d.deviceStruct.gridETime1, (*sCPU).ExtOut, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock/2, d.deviceStruct.Nthread, apertureFarFieldKernel, 
			sDevice, 
			(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter), 
			(deviceFP)activationParameter, 
			(deviceFP)(deg2Rad<deviceFP>() * xOffset),
			(deviceFP)(deg2Rad<deviceFP>() * yOffset));

		d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToHost);


		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);

		d.deviceLaunch((int)(d.deviceStruct.Ngrid / minGridDimension), 2 * minGridDimension, multiplyByConstantKernelD, 
			d.deviceStruct.gridETime1, (deviceFP)(1.0 / d.deviceStruct.Ngrid));
		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), copyType::ToHost);

		getTotalSpectrum(d);

		return 0;
	}

	static int applyAperature(ActiveDevice& d, const simulationParameterSet* sCPU, const double diameter, const double activationParameter) {
		d.deviceMemcpy(d.deviceStruct.gridETime1, (*sCPU).ExtOut, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToDevice);

		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock, d.deviceStruct.Nthread, apertureKernel, 
			sDevice, 
			(deviceFP)(0.5 * diameter), 
			(deviceFP)(activationParameter));
		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToHost);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applySphericalMirror(ActiveDevice& d, const simulationParameterSet* sCPU, deviceParameterSet<deviceFP, deviceComplex>& s, const double ROC) {

		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceMemcpy(sDevice, &s, sizeof(deviceParameterSet<deviceFP, deviceComplex>), copyType::ToDevice);

		d.deviceMemcpy(d.deviceStruct.gridETime1, (*sCPU).ExtOut, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, sphericalMirrorKernel, sDevice, (deviceFP)ROC);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.deviceLaunch(2 * d.deviceStruct.Nblock, d.deviceStruct.Nthread, multiplyByConstantKernelD, d.deviceStruct.gridETime1, (deviceFP)(1.0 / d.deviceStruct.Ntime));
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToHost);
		d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyParabolicMirror(ActiveDevice& d, simulationParameterSet* sCPU, deviceParameterSet<deviceFP, deviceComplex>& s, const double focus) {
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;

		d.deviceMemcpy(d.deviceStruct.gridETime1, (*sCPU).ExtOut, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, parabolicMirrorKernel, sDevice, (deviceFP)focus);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.deviceLaunch(2 * d.deviceStruct.Nblock, d.deviceStruct.Nthread, multiplyByConstantKernelD, d.deviceStruct.gridETime1, (deviceFP)(1.0 / d.deviceStruct.Ntime));
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToHost);
		d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyLinearPropagation(ActiveDevice& d, simulationParameterSet* sCPU, const int materialIndex, const double thickness) {

		if (d.hasPlasma) {
			simulationParameterSet sCopy = *sCPU;
			sCopy.nonlinearAbsorptionStrength = 0.0;
			d.reset(&sCopy);
		}

		d.deviceMemcpy(d.deviceStruct.gridEFrequency1, (*sCPU).EkwOut, d.deviceStruct.NgridC * 2 * sizeof(std::complex<double>), copyType::ToDevice);

		deviceFP* sellmeierCoefficients = (deviceFP*)d.deviceStruct.gridEFrequency1Next1;
		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(sellmeierCoefficientsAugmentedCPU, (*sCPU).crystalDatabase[materialIndex].sellmeierCoefficients.data(), 66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*sCPU).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*sCPU).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*sCPU).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*sCPU).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*sCPU).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*sCPU).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		d.deviceMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), copyType::ToDevice);
		d.deviceStruct.axesNumber = (*sCPU).crystalDatabase[materialIndex].axisType;
		d.deviceStruct.sellmeierType = (*sCPU).crystalDatabase[materialIndex].sellmeierType;
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;

		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, applyLinearPropagationKernel, sellmeierCoefficients, (deviceFP)thickness, sDevice);
		d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, d.deviceStruct.NgridC * 2 * sizeof(std::complex<double>), copyType::ToHost);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(2 * d.deviceStruct.Nblock, d.deviceStruct.Nthread, multiplyByConstantKernelD, d.deviceStruct.gridETime1, (deviceFP)(1.0 / d.deviceStruct.Ngrid));

		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToHost);
		getTotalSpectrum(d);

		return 0;
	}

	static int preparePropagationGrids(ActiveDevice& d) {
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
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
		memcpy(sellmeierCoefficientsAugmentedCPU + 72, (*s).crystalDatabase[(*s).materialIndex].nonlinearReferenceFrequencies.data(), 7 * sizeof(double));
		d.deviceMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, 79 * sizeof(double), copyType::ToDevice);

		//prepare the propagation grids
		deviceParameterSet<deviceFP, deviceComplex>* sD = d.dParamsDevice;
		d.deviceMemcpy(sD, sc, sizeof(deviceParameterSet<deviceFP, deviceComplex>), copyType::ToDevice);
		d.deviceLaunch((unsigned int)(*sc).Ntime/(2*minGridDimension), minGridDimension, getChiLinearKernel, sD, sellmeierCoefficients);
		if ((*s).is3D) {
			d.deviceLaunch((unsigned int)(*sc).Nblock / 2u, (unsigned int)(*sc).Nthread, prepare3DGridsKernel, sellmeierCoefficients, sD);
		}
		else if ((*s).isCylindric) {
			d.deviceLaunch((unsigned int)(*sc).Nblock / 2u, (unsigned int)(*sc).Nthread, prepareCylindricGridsKernel, sellmeierCoefficients, sD);
		}
		else {
			d.deviceLaunch((unsigned int)(*sc).Nblock / 2u, (unsigned int)(*sc).Nthread, prepareCartesianGridsKernel, sellmeierCoefficients, sD);
		}
		d.deviceMemcpy(sc, sD, sizeof(deviceParameterSet<deviceFP, deviceComplex>), copyType::ToHost);
		return 0;
	}

	//Rotate the field on the GPU
	//Allocates memory and copies from CPU, then copies back to CPU and deallocates
	// - inefficient but the general principle is that only the CPU memory is preserved
	// after simulations finish... and this only runs at the end of the simulation
	static int rotateField(ActiveDevice& d, simulationParameterSet* sCPU, const double rotationAngle) {

		deviceComplex* Ein1 = d.deviceStruct.gridEFrequency1;
		deviceComplex* Ein2 = d.deviceStruct.gridEFrequency2;
		deviceComplex* Eout1 = d.deviceStruct.gridEFrequency1Next1;
		deviceComplex* Eout2 = d.deviceStruct.gridEFrequency1Next2;

		//retrieve/rotate the field from the CPU memory
		d.deviceMemcpy(Ein1, (*sCPU).EkwOut, 2 * (*sCPU).NgridC * sizeof(std::complex<double>), copyType::ToDevice);
		d.deviceLaunch((unsigned int)(d.deviceStruct.NgridC / minGridDimension), minGridDimension, rotateFieldKernel, Ein1, Ein2, Eout1, Eout2, (deviceFP)rotationAngle);
		d.deviceMemcpy((*sCPU).EkwOut, Eout1, 2 * (*sCPU).NgridC * sizeof(std::complex<double>), copyType::ToHost);

		//transform back to time
		d.fft(Eout1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(2 * d.deviceStruct.Nblock, d.deviceStruct.Nthread, multiplyByConstantKernelD, d.deviceStruct.gridETime1, (deviceFP)(1.0 / d.deviceStruct.Ngrid));
		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), copyType::ToHost);

		//update spectrum
		getTotalSpectrum(d);
		return 0;
	}

//function to run a RK4 time step
//stepNumber is the sub-step index, from 0 to 3
	static int runRK4Step(ActiveDevice& d, uint8_t stepNumber) {
		deviceParameterSet<deviceFP, deviceComplex>* sH = d.s; 
		deviceParameterSet<deviceFP, deviceComplex>* sD = d.dParamsDevice;

		// Beam with symmetry around z axis:
		// Nonlinear polarization and plasma use expanded grid
		// Radial laplacian uses standard grid
		// two possible FFT shapes (radial Laplacian always performed)
		if((*sH).isCylindric){
			//ifft to time domain 
			d.fft((*sH).workspace1, (*sH).gridETime1, deviceFFT::Z2D);

			//Nonlinear polarization and plasma current are fft-ed in a batch
			//from separate (de-interlaced) time-domain grids.
			//assumption: no plasma without other nonlinearities
			if((*sH).isNonLinear){
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, nonlinearPolarizationKernel, sD);
				if((*sH).hasPlasma){
					d.deviceLaunch((*sH).Nblock, (*sH).Nthread, plasmaCurrentKernel_twoStage_A, sD);
					d.deviceLaunch((unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 2 * minGridDimension, plasmaCurrentKernel_twoStage_B, sD);
					d.deviceLaunch((*sH).Nblock, (*sH).Nthread, expandCylindricalBeam, sD);
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFT::D2Z_Polarization);
					d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPlasmaKernelCylindric, sD);
				}
				else{
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFT::D2Z_Polarization);
					d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPolarizationKernelCylindric, sD);
				}
			}
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, radialLaplacianKernel, sD);
			d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFT::D2Z);
		}
		// 2D and 3D cartesian
		// Only one type of FFT
		// currently nonlinear polarization and plasma ffts are not batched, could give
		// minor speed boost by combining them, but requires additional memory, so
		// nonlinear polarization and plasma are fft-ed separately to accommodate larger
		// systems.
		else if ((*sH).isNonLinear) {
			//perform inverse FFT to get time-space electric field
			d.fft((*sH).workspace1, (*sH).gridETime1, deviceFFT::Z2D);
			//Plasma/multiphoton absorption
			if ((*sH).hasPlasma) {
				d.deviceLaunch((*sH).Nblock, (*sH).Nthread, plasmaCurrentKernel_twoStage_A, sD);
				d.deviceLaunch((unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 2 * minGridDimension, plasmaCurrentKernel_twoStage_B, sD);
				d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFT::D2Z);
				d.deviceLaunch((*sH).Nblock / 2, (*sH).Nthread, updateKwithPlasmaKernel, sD);
			}
			//Nonlinear polarization
			d.deviceLaunch((*sH).Nblock, (*sH).Nthread, nonlinearPolarizationKernel, sD);
			d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFT::D2Z);
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



	static unsigned long int solveNonlinearWaveEquationWithDevice(ActiveDevice& d, simulationParameterSet* sCPU) {
		//prepare the propagation arrays
		preparePropagationGrids(d);
		prepareElectricFieldArrays(d);

		deviceFP* canaryPointer = &d.deviceStruct.gridETime1[d.deviceStruct.Ntime / 2 + d.deviceStruct.Ntime * (d.deviceStruct.Nspace / 2 + d.deviceStruct.Nspace * (d.deviceStruct.Nspace2 / 2))];
		//Core propagation loop
		for (size_t i = 0; i < d.deviceStruct.Nsteps; ++i) {

			//RK4
			runRK4Step(d, 0);
			runRK4Step(d, 1);
			runRK4Step(d, 2);
			runRK4Step(d, 3);

			//periodically check if the simulation diverged or was cancelled
			if ((*sCPU).cancellationCalled) break;
			if (i % 10 == 0) if (d.isTheCanaryPixelNaN(canaryPointer)) break;
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
		}
		if ((*sCPU).isInFittingMode && !(*sCPU).isInSequence)(*(*sCPU).progressCounter)++;

		//take final spectra and transfer the results to the CPU
		d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToHost);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch((int)(d.deviceStruct.Ngrid / minGridDimension), 2 * minGridDimension, multiplyByConstantKernelD, d.deviceStruct.gridETime1, (deviceFP)(1.0 / d.deviceStruct.Ngrid));
		d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), copyType::ToHost);
		getTotalSpectrum(d);

		int returnval =  13 * d.isTheCanaryPixelNaN(canaryPointer);
		return returnval;
	}

	static constexpr unsigned int funHash(const char* s, int off = 0) {
		return (s[off] == 0 || s[off] == '(') ? 7177 : (funHash(s, off + 1) * 31) ^ s[off];
	}

	static unsigned int stringHash(std::string& s, int off = 0){
		return (s.length() == off || s.at(off) == '(') ? 7177 : (stringHash(s,off+1) * 31) ^ s.at(off);
	}
	//Dispatcher of the sequence mode. New functions go here, and should have a unique hash (chances of a conflict are small, and 
	// will produce a compile-time error.
	// Functions cannot start with a number or the string "None".
	static int interpretCommand(std::string cc, double* iBlock, double* vBlock, ActiveDevice& d, simulationParameterSet *sCPU) {
		crystalEntry* db = (*sCPU).crystalDatabase;
		int error = 0;
		double parameters[32] = {0.0};
		bool defaultMask[32] = {0};

		switch (stringHash(cc)) {
		case funHash("rotate"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			rotateField(d, sCPU, deg2Rad<deviceFP>() * parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("set"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			vBlock[(int)parameters[0]] = parameters[1];
			break;
		case funHash("plasmaReinject"):
			(*sCPU).isReinjecting = true;
			[[fallthrough]];
		case funHash("plasma"):
		{
			interpretParameters(cc, 9, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
			if (!defaultMask[3])(*sCPU).nonlinearAbsorptionStrength = parameters[3];
			if (!defaultMask[4])(*sCPU).bandGapElectronVolts = parameters[4];
			if (!defaultMask[5])(*sCPU).drudeGamma = parameters[5];
			if (!defaultMask[6])(*sCPU).effectiveMass = parameters[6];
			if (!defaultMask[7])(*sCPU).crystalThickness = 1e-6 * parameters[7];
			if (!defaultMask[8])(*sCPU).propagationStep = 1e-9 * parameters[8];
			(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d.data();
			(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3.data();
			(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches.data();
			(*sCPU).absorptionParameters = db[(*sCPU).materialIndex].absorptionParameters.data();
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();

			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU);
			(*sCPU).isFollowerInSequence = true;
		}
			break;
		case funHash("nonlinear"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
			if (!defaultMask[3])(*sCPU).crystalThickness = 1e-6 * parameters[3];
			if (!defaultMask[4])(*sCPU).propagationStep = 1e-9 * parameters[4];

			(*sCPU).nonlinearAbsorptionStrength = 0.0;
			(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d.data();
			(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3.data();
			(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches.data();
			(*sCPU).absorptionParameters = db[(*sCPU).materialIndex].absorptionParameters.data();
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();

			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU);
			(*sCPU).isFollowerInSequence = true;
			break;
		case funHash("default"):
			d.reset(sCPU);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU);
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
			(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d.data();
			(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3.data();
			(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches.data();
			(*sCPU).absorptionParameters = db[(*sCPU).materialIndex].absorptionParameters.data();
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();

			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			(*sCPU).isFollowerInSequence = true;
			break;

		case funHash("linear"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if ((*sCPU).isCylindric) {
				if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
				if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
				if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
				if (!defaultMask[3])(*sCPU).crystalThickness = 1e-6 * parameters[3];
				if (!defaultMask[4])(*sCPU).propagationStep = 1e-9 * parameters[4];
				(*sCPU).nonlinearAbsorptionStrength = 0.0;
				(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d.data();
				(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3.data();
				(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches.data();
				(*sCPU).absorptionParameters = db[(*sCPU).materialIndex].absorptionParameters.data();
				(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();
				(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
				(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
				(*sCPU).forceLinear = true;
				d.reset(sCPU);
				error = solveNonlinearWaveEquationWithDevice(d, sCPU);
				(*sCPU).isFollowerInSequence = true;
			}
			else {
				if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
				if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
				if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
				if (!defaultMask[3])(*sCPU).crystalThickness = 1e-6 * parameters[3];
				if (d.hasPlasma) {
					(*sCPU).nonlinearAbsorptionStrength = 0.0;
					(*sCPU).forceLinear = true;
					d.reset(sCPU);
				}

				applyLinearPropagation(d, sCPU, (*sCPU).materialIndex, (*sCPU).crystalThickness);
				if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			}

			break;
		case funHash("fresnelLoss"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
			d.reset(sCPU);
			applyFresnelLoss(d, sCPU, d.deviceStruct,
				(int)parameters[4],
				(int)parameters[5]);
			break;
		case funHash("sphericalMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applySphericalMirror(d, sCPU, d.deviceStruct, parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("parabolicMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyParabolicMirror(d, sCPU, d.deviceStruct, parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("aperture"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyAperature(d, sCPU,
				parameters[0],
				parameters[1]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("farFieldAperture"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyAperatureFarField(d, sCPU,
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
			d.reset(sCPU);
			applyFilter(d, sCPU,
				parameters[0],
				parameters[1],
				parameters[2],
				parameters[3],
				parameters[4]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("lorentzian"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyLorenzian(d, sCPU, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case funHash("addPulse"):
			if ((*sCPU).runType == -1) break;
		{
			interpretParameters(cc, 21, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			d.deviceMemcpy(d.deviceStruct.gridETime1, (*sCPU).ExtOut, 2 * d.deviceStruct.Ngrid * sizeof(double), copyType::ToDevice);
			d.deviceMemcpy(d.deviceStruct.gridEFrequency1, (*sCPU).EkwOut, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToDevice);

			pulse<double> p;
			p = sCPU->pulse1;
			p.energy = parameters[0];
			p.frequency = 1e12 * parameters[1];
			p.bandwidth = 1e12 * parameters[2];
			p.sgOrder = (int)parameters[3];
			p.cep = parameters[4] * vPi<deviceFP>();
			p.delay = 1e-15 * parameters[5];
			p.gdd = 1e-30 * parameters[6];
			p.tod = 1e-45 * parameters[7];
			p.phaseMaterial = (int)parameters[8];
			p.phaseMaterialThickness = 1e-6 * parameters[9];
			p.beamwaist = 1e-6 * parameters[10];
			p.x0 = 1e-6 * parameters[11];
			p.y0 = 1e-6 * parameters[12];
			p.z0 = 1e-6 *parameters[13];
			p.beamAngle = deg2Rad<deviceFP>() * parameters[14];
			p.beamAnglePhi = deg2Rad<deviceFP>() * parameters[15];
			p.polarizationAngle = deg2Rad<deviceFP>() * parameters[16];
			p.circularity = parameters[17];
			(*sCPU).materialIndex = (int)parameters[18];
			(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[19];
			(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[20];

			addPulseToFieldArrays(d, p, false, NULL);
			d.deviceMemcpy((*sCPU).EkwOut, d.deviceStruct.gridEFrequency1, 2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), copyType::ToHost);
			d.deviceMemcpy((*sCPU).ExtOut, d.deviceStruct.gridETime1, 2 * (*sCPU).Ngrid * sizeof(double), copyType::ToHost);

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
					error = interpretCommand(currentString, iBlock, vBlock, d, sCPU);
					currentString = currentString.substr(currentString.find_first_of(')') + 1, std::string::npos);
					if (error || (*sCPU).cancellationCalled) break;
				}
				++vBlock[targetVar];
				currentString = forStartString;
				if (error || (*sCPU).cancellationCalled) break;
			}
			break;
		}
		return error;
	}

	static int solveSequenceWithDevice(ActiveDevice& d, simulationParameterSet* sCPU) {
		int error = 0;

		//if it starts with 0, it's an old sequence; quit
		if ((*sCPU).sequenceString[0] == '0') {
			return 15;
		}

		//main text interpreter
		simulationParameterSet sCPUbackupValues;
		simulationParameterSet* sCPUbackup = &sCPUbackupValues;
		sCPUbackupValues = *sCPU;
		double iBlock[100] = { 0.0 };

		for (int k = 1; k < 38; k++) {
			iBlock[k] = (*sCPU).getByNumberWithMultiplier(k);
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

			error = interpretCommand(currentString, iBlock, vBlock, d, sCPU);
			if (error || (*sCPU).cancellationCalled) break;
			currentString = currentString.substr(currentString.find_first_of(')'), std::string::npos);

			if (currentString.length() < minLength) break;

			currentString = currentString.substr(1, std::string::npos);

			(*sCPUbackup).isFollowerInSequence = (*sCPU).isFollowerInSequence;
			//memcpy(sCPU, sCPUbackup, sizeof(simulationParameterSet));
			*sCPU = *sCPUbackup;
		}
		return error;
	}

	// helper function for fitting mode, runs the simulation and returns difference from the desired outcome.
	static double getResidual(const dlib::matrix<double, 0, 1>& x) {
		double result = 0.0;
		for (int i = 0; i < (*fittingSet).Nfitting; ++i) {
			(*fittingSet).setByNumberWithMultiplier((size_t)(*fittingSet).fittingArray[3 * i], x(i));
		}

		ActiveDevice& d = *dFit;
		d.cParams = fittingSet;
		d.reset(fittingSet);

		if ((*fittingSet).isInSequence) {
			(*fittingSet).isFollowerInSequence = false;
			solveSequenceWithDevice(d, fittingSet);
			(*(*fittingSet).progressCounter)++;
		}
		else {
			solveNonlinearWaveEquationWithDevice(d, fittingSet);
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
unsigned long solveNonlinearWaveEquationX(simulationParameterSet* lpParam) {
	simulationParameterSet* sCPU = (simulationParameterSet*)lpParam;
	ActiveDevice d(sCPU);
	if (d.memoryStatus) return 1;
	unsigned long returnValue = solveNonlinearWaveEquationWithDevice(d, sCPU);
	return returnValue;
}

// Main function for running a sequence
unsigned long solveNonlinearWaveEquationSequenceX(simulationParameterSet* lpParam) {
	simulationParameterSet sCPUcurrent = *((simulationParameterSet*)lpParam);
	simulationParameterSet* sCPU = &sCPUcurrent;
	if ((*sCPU).batchIndex == 36 && (*sCPU).batchLoc1 != 0) return 0;
	ActiveDevice d(sCPU);
	unsigned long returnValue = solveSequenceWithDevice(d, sCPU);
	return returnValue;
}

//run in fitting mode
unsigned long runDlibFittingX(simulationParameterSet* sCPU) {
	simulationParameterSet sCPUcurrent = *sCPU;
	fittingSet = &sCPUcurrent;

	ActiveDevice d(fittingSet);
	dFit = &d;
	dlib::matrix<double, 0, 1> parameters;
	parameters.set_size((*sCPU).Nfitting);
	dlib::matrix<double, 0, 1> lowerBounds;
	lowerBounds.set_size((*sCPU).Nfitting);
	dlib::matrix<double, 0, 1> upperBounds;
	upperBounds.set_size((*sCPU).Nfitting);

	for (int i = 0; i < (*sCPU).Nfitting; ++i) {
		parameters(i) = (*sCPU).getByNumber((size_t)round((*sCPU).fittingArray[3 * i]));
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
		(*sCPU).setByNumberWithMultiplier((size_t)round((*sCPU).fittingArray[3 * i]), result.x(i));
		(*sCPU).fittingResult[i] = result.x(i);
	}

	std::atomic_uint32_t fitCounter{ 0 };
	std::atomic_uint32_t* originalCounter = (*sCPU).progressCounter;
	(*sCPU).progressCounter = &fitCounter;
	d.cParams = sCPU;
	d.reset(sCPU);
	if ((*sCPU).isInSequence) {
		solveSequenceWithDevice(d, sCPU);
	}
	else {
		solveNonlinearWaveEquationWithDevice(d, sCPU);
	}

	(*sCPU).progressCounter = originalCounter;

	dFit = nullptr;
	return 0;
}



