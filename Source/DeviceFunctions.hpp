#include "LightwaveExplorerTrilingual.h"

namespace deviceFunctions {
	//strict isnan for complex type
	deviceFunction static inline bool isComplexNaN(const deviceComplex& x){
		return isnan(x.real()) || isnan(x.imag());
	}


	//Expand the information contained in the radially-symmetric beam in the offset grid
	// representation.
	// see the expandCylindricalBeam() kernel for more details
	template<typename T, typename U>
	deviceFunction static void expandCylindricalBeamDevice(
		const deviceParameterSet<T,U>* s, 
		const int64_t i, 
		T* expandedBeam1, 
		const T* sourceBeam1, 
		const T* sourceBeam2) {
		const int64_t j = i / (*s).Ntime; //spatial coordinate
		const int64_t k = i - j * (*s).Ntime; //temporal coordinate

		//positions on the expanded grid corresponding the the current index
		const int64_t pos1 = 2 * ((*s).Nspace - j - 1) * (*s).Ntime + k;
		const int64_t pos2 = (2 * j + 1) * (*s).Ntime + k;
		T* expandedBeam2 = expandedBeam1 + 2 * (*s).Ngrid;
		expandedBeam1[pos1] = sourceBeam1[i];
		expandedBeam1[pos2] = sourceBeam1[i];
		expandedBeam2[pos1] = sourceBeam2[i];
		expandedBeam2[pos2] = sourceBeam2[i];
	}

	template<typename T>
	deviceFunction static T lagrangeInterpolation(T x, T* yData) {
		T result{};
		for (int i = 0; i < 4; i++) {
			T term = yData[i];
			for (int j = 0; j < 4; j++) {
				if (j != i) term *= (x - (j - 1)) / (i - j);
			}
			result += term;
		}
		return result;
	}

	//Calculate the fourier optics propagator (e^ik_z*d) 
	// for a given set of values of the maknitude of k, transverse k (dk1, dk2)
	// a reference k0 which defines the speed of the moving frame, 
	// and distance d over which to propagate
	template<typename real_t, typename complex_t>
	deviceFunction static complex_t fourierPropagator(
		const complex_t k, 
		const real_t dk1, 
		const real_t dk2, 
		const real_t k0, 
		const real_t d) {

		if (isnan(k0) 
		|| isnan(dk1) 
		|| isnan(dk2) 
		|| isComplexNaN(k) 
		|| k.real() == real_t{}) return complex_t{};
		if ( (deviceFPLib::abs(dk2) < 0.1f * k.real()) 
			&& (deviceFPLib::abs(dk1) < 0.1f *  k.real())) {
			deviceFP halfOverKr = (0.5f / k.real()) * (dk1 * dk1 + dk2 * dk2);
			deviceFP kMoving = k.real() - k0;
			complex_t returnVal = complex_t(
					d * k.imag(),
					d * (halfOverKr - kMoving)
				);
			return isComplexNaN(returnVal) ? complex_t{} : deviceLib::exp(returnVal);
		}
		complex_t kz = 
			(deviceLib::sqrt(-dk2 * dk2 / (k + deviceFPLib::abs(dk1)) + k - deviceFPLib::abs(dk1))
				* deviceLib::sqrt(k + deviceFPLib::abs(dk1)) - k0);
		if (kz.imag() > 0.0f) kz = complex_t(kz.real(), -kz.imag());
		kz = complex_t(d*kz.imag(),-d*kz.real());
		if (isComplexNaN(kz)) return complex_t{};
		return deviceLib::exp(kz);
	}

	//give the Dawson function value at x
	//used for Gaussian-based "Sellmeier" equation, 
	// as it is the Hilbert transform of the Gaussian function 
	// (susceptibility obeys Kramers-Kronig relations)
	//based on Rybicki G.B., Computers in Physics, 3,85-87 (1989)
	//this is the simplest implementation of the formula he provides, 
	// he also suggests speed improvements in case
	//evaluation of this becomes a bottleneck
	//macro for the two versions because I don't want to lose precision in the double version
	//but Intel Xe graphics need the constants as floats, and CUDA doesn't work with the
	//[[maybe_unused]] attribute... and I don't like compiler warnings.
	[[maybe_unused]] deviceFunction static float deviceDawson(
		const float x) {

		//parameters determining accuracy (higher n, smaller h -> more accurate but slower)
		const int n = 15;
		const float h = 0.3f;

		//series expansion for small x
		if (deviceFPLib::abs(x) < 0.2f) {
			const float x2 = x * x;
			const float x4 = x2 * x2;
			return x * (1.0f 
				- 2.0f * x2 / 3.0f 
				+ 4.0f * x4 / 15.0f 
				- 8.0f * x2 * x4 / 105.0f 
				+ (16.0f / 945.0f) * x4 * x4 
				- (32.0f / 10395.0f) * x4 * x4 * x2);
		}

		const int n0 = 2 * (int)(round(0.5f * x / h));
		const float x0 = h * n0;
		const float xp = x - x0;
		float d = 0.0f;
		for (int i = -n; i < n; i++) {
			if (i % 2 != 0) {
				d += deviceFPLib::exp(-(xp - i * h) * (xp - i * h)) / (i + n0);
			}
		}
		return invSqrtPi<float>() * d;
	}

	[[maybe_unused]] deviceFunction static double deviceDawson(
		const double x) {

		//parameters determining accuracy (higher n, smaller h -> more accurate but slower)
		const int n = 15;
		const double h = 0.3;

		//series expansion for small x
		if (deviceFPLib::abs(x) < 0.2) {
			const double x2 = x * x;
			const double x4 = x2 * x2;
			return x * (1.0 - 2.0 * x2 / 3.0 
				+ 4.0 * x4 / 15.0 
				- 8.0 * x2 * x4 / 105.0 
				+ (16.0 / 945) * x4 * x4 
				- (32.0 / 10395) * x4 * x4 * x2);
		}

		const int n0 = 2 * (int)(round(0.5 * x / h));
		const double x0 = h * n0;
		const double xp = x - x0;
		double d = 0.0;
		for (int i = -n; i < n; i++) {
			if (i % 2 != 0) {
				d += deviceFPLib::exp(-(xp - i * h) * (xp - i * h)) / (i + n0);
			}
		}
		return invSqrtPi<double>() * d;
	}
	//Inner function for the Sellmeier equation to provide the refractive indices
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
	deviceFunction static deviceComplex sellmeierFunc(
		deviceFP ls, 
		const deviceFP omega, 
		const deviceFP* a, 
		const int eqn,
		const bool takeSqrt = true) {

		const deviceFP omega2 = omega * omega;
		deviceFP realPart;
		deviceComplex compPart;
		switch (eqn) {
		case 0:
			[[unlikely]]
			if (ls == -a[3] || ls == -a[6] || ls == -a[9] || ls == -a[12]) return deviceComplex{};
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
			return takeSqrt ? deviceLib::sqrt(maxN(realPart, 0.9f) + compPart)
				: maxN(realPart, 0.9f) + compPart;
		case 1:
			//up to 7 Lorentzian lines
			compPart = a[1] / deviceComplex(a[2] - omega2, a[3] * omega)
				+ a[4] / deviceComplex(a[5] - omega2, a[6] * omega)
				+ a[7] / deviceComplex(a[8] - omega2, a[9] * omega)
				+ a[10] / deviceComplex(a[11] - omega2, a[12] * omega)
				+ a[13] / deviceComplex(a[14] - omega2, a[15] * omega)
				+ a[16] / deviceComplex(a[17] - omega2, a[18] * omega)
				+ a[19] / deviceComplex(a[20] - omega2, a[21] * omega);
			compPart *= kLorentzian<deviceFP>();
			compPart += a[0];
			return takeSqrt ?
				deviceComplex((
					deviceLib::sqrt(compPart)).real(),
					-deviceFPLib::abs((deviceLib::sqrt(compPart)).imag()))
				: compPart;
		case 2:
		{
			//Up to 7 complex Gaussian functions
			//the real part is the Hilbert transform of the Gaussian, th
			deviceFP scaledF;
			compPart = deviceComplex(a[0], 0.0f);
			for (int i = 0; i < 7; ++i) {
				if (a[3 * i + 1] != 0.0f){
					scaledF = (omega - a[1 + 3 * i]) / (sqrtTwo<deviceFP>() * a[2 + 3 * i]);
					compPart += a[3 + 3 * i] * deviceComplex(
						-invSqrtPi<deviceFP>() * deviceDawson(scaledF), 
						-deviceFPLib::exp(-scaledF * scaledF));
				}

			}
			//always select branch with imaginary part < 0
			return takeSqrt ?
				deviceComplex((
					deviceLib::sqrt(compPart)).real(),
					-deviceFPLib::abs((deviceLib::sqrt(compPart)).imag()))
				: compPart;

		}
		case 100:
		{
			[[unlikely]] if (ls == -a[3] 
				|| ls == -a[6] 
				|| ls == -a[9] 
				|| ls == -a[12]) return deviceComplex{};
		}
			realPart = a[0]
				+ (a[1] + a[2] * ls) / (ls + a[3])
				+ (a[4] + a[5] * ls) / (ls + a[6])
				+ (a[7] + a[8] * ls) / (ls + a[9])
				+ (a[10] + a[11] * ls) / (ls + a[12])
				+ a[13] * ls
				+ a[14] * ls * ls
				+ a[15] * ls * ls * ls;
			//"real-valued equation has no business being < 1" - causality
			return takeSqrt ?
				deviceComplex(
					deviceFPLib::sqrt(maxN(realPart, 0.9f)),
					0.0f)
				: maxN(realPart, 0.9f);
		default: return cOne<deviceComplex>();
		}
		return cOne<deviceComplex>();
	}
	
	deviceFunction static inline deviceFP square(const deviceFP x){
		return x*x;
	}

	//Sellmeier equation for refractive indices
	template<typename deviceFP, typename deviceComplex>
	deviceFunction static deviceFP sellmeierDevice(
		deviceComplex* ne, 
		deviceComplex* no, 
		const deviceFP* a, 
		const deviceFP f, 
		const deviceFP theta, 
		const deviceFP phi, 
		const int type, 
		const int eqn,
		const bool takeSqrt = true) {
		
		if (f == 0.0f) {
			*ne = cOne<deviceComplex>(); 
			*no = cOne<deviceComplex>(); 
			return deviceFP{};
		} //exit immediately for f=0
		deviceFP ls = 2.99792458e14f / f; //wavelength in microns
		ls *= ls; //only wavelength^2 is ever used
		const deviceFP omega = twoPi<deviceFP>() * maxN(f,-f);

		//option 0: isotropic
		if (type == 0) {
			*ne = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, a, eqn, takeSqrt);
			*no = *ne;
			return deviceFP{};
		}
		
		//option 1: uniaxial
		else if (type == 1) {
			deviceComplex na = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, a, eqn, false);
			deviceComplex nb = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, &a[22], eqn, false);

			
			deviceFP cosT = deviceFPLib::cos(theta);
			cosT *= cosT;
			deviceFP sinT = deviceFPLib::sin(theta);
			sinT *= sinT;
			if (takeSqrt) {
				*no = deviceLib::sqrt(na);
				*ne = deviceLib::sqrt((na * nb) / (nb * cosT + na * sinT));
			}
			else {
				*no = na;
				*ne = (na * nb) / (nb * cosT + na * sinT);
			}
			return deviceFP{};
		}
		//option 2: biaxial
		else {
			deviceComplex na = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, a, eqn, false);
			deviceComplex nb = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, &a[22], eqn, false);
			deviceComplex nc = sellmeierFunc<deviceFP, deviceComplex>(ls, omega, &a[44], eqn, false);
			deviceFP cp = deviceFPLib::cos(phi);
			deviceFP sp = deviceFPLib::sin(phi);
			deviceFP ct = deviceFPLib::cos(theta);
			deviceFP st = deviceFPLib::sin(theta);

			deviceFP d = (na.real() == nb.real()) ? 
			deviceFP{} : 
			0.5f * deviceFPLib::atan(deviceFPLib::sin(2*phi) * ct / 
			(((1.0f/nb.real() - 1.0f/nc.real())/(1.0f/na.real() - 1.0f/nb.real())) * st*st - cp*cp * ct*ct + sp*sp));
			deviceFP cd = deviceFPLib::cos(d);
			deviceFP sd = deviceFPLib::sin(d);

			deviceComplex nTop = na * nb * nc;

			*ne = nTop / 
			(na*nb*st*st*cd*cd 
			+ na*nc * square(sd*cp + sp*cd*ct) 
			+ nb*nc * square(sd*sp - cd*cp*ct));

			*no = nTop / 
			(na*nb*st*st*sd*sd
			+ na*nc * square(sd*sp*ct - cd*cp)
			+ nb*nc * square(sd*cp*ct + sp*cd));



			if (takeSqrt) {
				*ne = deviceLib::sqrt(*ne);
				*no = deviceLib::sqrt(*no);
			}
			return d;
		}
	}

	deviceFunction static inline deviceFP modulusSquared(
		const deviceComplex& a) {
		return a.real() * a.real() + a.imag() * a.imag();
	}

	//provide the position rho in cylindric mode
	template<typename deviceFP>
	deviceFunction static deviceFP rhoInRadialSymmetry(
		const int64_t N, 
		const int64_t j, 
		const deviceFP dr) {

		if (j < N / 2) {
			return deviceFPLib::abs((dr * (N / 2 - j) - 0.25f * dr));
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
	// Use OGM1; D. Kim, J.A. Fessler, 
	// Optimized first-order methods for smooth convex minimization, arXiv:1406.5468
	template<typename deviceFP, typename deviceComplex>
	deviceFunction static deviceFP findBirefringentCrystalIndex(
		const deviceParameterSet<deviceFP, deviceComplex>* s, 
		const deviceFP* sellmeierCoefficients, 
		const int64_t i, 
		deviceComplex* n1, 
		deviceComplex* n2) {

		int64_t h = 1 + i % ((*s).Nfreq - 1);
		int64_t col = i / ((*s).Nfreq - 1);
		int64_t k = col / (*s).Nspace;
		int64_t j = col - k * (*s).Nspace;
		

		deviceFP f = (*s).fStep * h;
		deviceFP kx1 = (lightC<deviceFP>() 
			/ (twoPi<deviceFP>() * f)) 
			* (j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace));
		deviceFP ky1 = (*s).is3D ? 
			(lightC<deviceFP>() 
			/ (twoPi<deviceFP>() * f)) 
			* (k * (*s).dk2 - (k >= ((*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2))
			:
			deviceFP{};
		//alpha is deviation from crystal Theta (x2 polarizations)
		//beta is deviation from crystal Phi
		//
		deviceComplex n[4][2]{};
		deviceComplex nW = deviceComplex{};
		sellmeierDevice(
			&n[0][0], 
			&n[0][1], 
			sellmeierCoefficients, 
			f, 
			(*s).crystalTheta, 
			(*s).crystalPhi, 
			(*s).axesNumber, 
			(*s).sellmeierType);
		if ((*s).axesNumber == 0) {
			*n1 = n[0][0];
			*n2 = n[0][1];
			return deviceFP{};
		}

		deviceFP gradient[2][2]{};
		deviceFP alpha[2] = 
		{ deviceFPLib::asin(kx1 / n[0][0].real()),deviceFPLib::asin(kx1 / n[0][1].real()) };
		deviceFP beta[2] = 
		{ deviceFPLib::asin(ky1 / n[0][0].real()), deviceFPLib::asin(ky1 / n[0][1].real()) };

		deviceFP gradientStep = 1.0e-6f;
		deviceFP gradientFactor = 0.5f / gradientStep;
		int it = 0;
		int maxiter = 64;
		deviceFP gradientTol = 1e-3f;
		//emperical testing: 
		// converges to double precision limit in two iterations for BBO
		// converges in 32 iterations in BiBO

		deviceFP errArray[4][2]{};
		if ((*s).axesNumber == 1) {
			maxiter = 64;
			sellmeierDevice(
				&n[0][0], 
				&nW, 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[0] + gradientStep, 
				(*s).crystalPhi, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&n[1][0], 
				&nW, 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[0] - gradientStep, 
				(*s).crystalPhi, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			if (isComplexNaN(n[0][0]) 
				|| isComplexNaN(n[1][0])) {
				*n1 = deviceComplex{};
				*n2 = deviceComplex{};
				return deviceFP{};
			}
			errArray[0][0] = deviceFPLib::sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
			errArray[1][0] = deviceFPLib::sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
			gradient[0][0] = gradientFactor * (errArray[0][0] - errArray[1][0]);

			for (it = 0; it < maxiter; ++it) {
				if (isComplexNaN(n[0][0]) 
					|| isComplexNaN(n[1][0])) {
					*n1 = deviceComplex{};
					*n2 = deviceComplex{};
					return deviceFP{};
				}
				if (deviceFPLib::abs(gradient[0][0]) > gradientTol) {
					alpha[0] -= 0.5f * (errArray[0][0] + errArray[1][0]) / gradient[0][0];
				} 
				else {
					break;
				}

				sellmeierDevice(
					&n[0][0], 
					&nW, 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[0] + gradientStep, 
					(*s).crystalPhi, 
					(*s).axesNumber, 
					(*s).sellmeierType);
				sellmeierDevice(
					&n[1][0], 
					&nW, 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[0] - gradientStep, 
					(*s).crystalPhi, (*s).axesNumber, 
					(*s).sellmeierType);
				errArray[0][0] = deviceFPLib::sin(alpha[0] + gradientStep) * n[0][0].real() - kx1;
				errArray[1][0] = deviceFPLib::sin(alpha[0] - gradientStep) * n[1][0].real() - kx1;
				gradient[0][0] = gradientFactor * (errArray[0][0] - errArray[1][0]);
			}
			sellmeierDevice(
				&n[0][0], 
				&nW, 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[0], 
				(*s).crystalPhi, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&nW, 
				&n[1][1], 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[1], 
				(*s).crystalPhi, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			*n1 = n[0][0];
			*n2 = n[1][1];
			return deviceFP{};
		}

		if ((*s).axesNumber == 2) {
			sellmeierDevice(
				&n[0][0], 
				&nW, 
				sellmeierCoefficients, 
				f, (*s).crystalTheta + alpha[0] + gradientStep, 
				(*s).crystalPhi + beta[0], 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&n[1][0], 
				&nW, 
				sellmeierCoefficients, 
				f, (*s).crystalTheta + alpha[0] - gradientStep, 
				(*s).crystalPhi + beta[0], 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&n[2][0], 
				&nW, 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[0], 
				(*s).crystalPhi + beta[0] + gradientStep, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&n[3][0], 
				&nW, 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[0], 
				(*s).crystalPhi + beta[0] - gradientStep, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&nW, 
				&n[0][1], 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[1] + gradientStep, 
				(*s).crystalPhi + beta[1], 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&nW, 
				&n[1][1], 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[1] - gradientStep, 
				(*s).crystalPhi + beta[1], 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&nW, 
				&n[2][1], 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[1], 
				(*s).crystalPhi + beta[1] + gradientStep, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&nW, 
				&n[3][1], 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[1], 
				(*s).crystalPhi + beta[1] - gradientStep, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			if (isComplexNaN(n[0][0]) 
				|| isComplexNaN(n[1][0])) {
				*n1 = n[0][0];
				*n2 = n[0][1];
				return deviceFP{};
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
				if (isComplexNaN(n[0][0])
					|| isComplexNaN(n[1][0])) {
					*n1 = deviceComplex{};
					*n2 = deviceComplex{};
					return deviceFP{};
				}
				if (deviceFPLib::abs(gradient[0][0]) > 1e-2f) alpha[0] 
					-= 0.25f * (errArray[0][0] + errArray[1][0]) / gradient[0][0];
				if (deviceFPLib::abs(gradient[1][0]) > 1e-2f) beta[0] 
					-= 0.25f * (errArray[2][0] + errArray[3][0]) / gradient[1][0];
				if (deviceFPLib::abs(gradient[0][1]) > 1e-2f) alpha[1] 
					-= 0.25f * (errArray[0][1] + errArray[1][1]) / gradient[0][1];
				if (deviceFPLib::abs(gradient[1][1]) > 1e-2f) beta[1] 
					-= 0.25f * (errArray[2][1] + errArray[3][1]) / gradient[1][1];

				if (
					maxN(
						maxN(
							deviceFPLib::abs(gradient[0][0]), 
							deviceFPLib::abs(gradient[1][0])), 
						maxN(
							deviceFPLib::abs(gradient[0][1]), 
							deviceFPLib::abs(gradient[1][1]))) < gradientTol) break;
				sellmeierDevice(
					&n[0][0], 
					&nW, 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[0] + gradientStep, 
					(*s).crystalPhi + beta[0], 
					(*s).axesNumber, 
					(*s).sellmeierType);
				sellmeierDevice(
					&n[1][0], 
					&nW, 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[0] - gradientStep, 
					(*s).crystalPhi + beta[0], 
					(*s).axesNumber, 
					(*s).sellmeierType);
				sellmeierDevice(
					&n[2][0], 
					&nW, 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[0], 
					(*s).crystalPhi + beta[0] + gradientStep, 
					(*s).axesNumber, 
					(*s).sellmeierType);
				sellmeierDevice(
					&n[3][0], 
					&nW, 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[0], 
					(*s).crystalPhi + beta[0] - gradientStep, 
					(*s).axesNumber, 
					(*s).sellmeierType);
				sellmeierDevice(
					&nW, 
					&n[0][1], 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[1] + gradientStep, 
					(*s).crystalPhi + beta[1], 
					(*s).axesNumber, 
					(*s).sellmeierType);
				sellmeierDevice(
					&nW, 
					&n[1][1], 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[1] - gradientStep, 
					(*s).crystalPhi + beta[1], 
					(*s).axesNumber, 
					(*s).sellmeierType);
				sellmeierDevice(
					&nW, 
					&n[2][1], 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[1], 
					(*s).crystalPhi + beta[1] + gradientStep, 
					(*s).axesNumber, 
					(*s).sellmeierType);
				sellmeierDevice(
					&nW, 
					&n[3][1], 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta + alpha[1], 
					(*s).crystalPhi + beta[1] - gradientStep, 
					(*s).axesNumber, 
					(*s).sellmeierType);
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
			deviceFP d = sellmeierDevice(
				&n[0][0], 
				&nW, 
				sellmeierCoefficients, 
				f, (*s).crystalTheta + alpha[0], 
				(*s).crystalPhi + beta[0], 
				(*s).axesNumber, 
				(*s).sellmeierType);
			sellmeierDevice(
				&nW, 
				&n[1][1], 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta + alpha[1], 
				(*s).crystalPhi + beta[1], 
				(*s).axesNumber, 
				(*s).sellmeierType);
			*n1 = n[0][0];
			*n2 = n[1][1];
			return d;
		}
		return deviceFP{};
	}
}
