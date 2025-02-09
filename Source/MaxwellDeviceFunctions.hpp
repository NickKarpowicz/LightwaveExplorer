#include "LightwaveExplorerTrilingual.h"
namespace deviceFunctions{

	//define math operators for maxwell grid points
	deviceFunction static inline deviceFP dotProduct(
		const maxwellPoint<deviceFP>& a, 
		const maxwellPoint<deviceFP>& b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	deviceFunction static maxwellPoint<deviceFP> rotateMaxwellPoint(
		const maxwell3D* s, 
		const maxwellPoint<deviceFP>& in, 
		const bool backwards) {
		if (backwards) {
			return maxwellPoint<deviceFP>{
				(*s).rotateBackward[0][0] * in.x 
					+ (*s).rotateBackward[1][0] * in.y 
					+ (*s).rotateBackward[2][0] * in.z,
				(*s).rotateBackward[3][0] * in.x 
					+ (*s).rotateBackward[4][0] * in.y 
					+ (*s).rotateBackward[5][0] * in.z,
				(*s).rotateBackward[6][0] * in.x 
					+ (*s).rotateBackward[7][0] * in.y 
					+ (*s).rotateBackward[8][0] * in.z};
		}
		return maxwellPoint<deviceFP>{
			(*s).rotateForward[0][0] * in.x 
				+ (*s).rotateForward[1][0] * in.y 
				+ (*s).rotateForward[2][0] * in.z,
			(*s).rotateForward[3][0] * in.x 
				+ (*s).rotateForward[4][0] * in.y 
				+ (*s).rotateForward[5][0] * in.z,
			(*s).rotateForward[6][0] * in.x 
				+ (*s).rotateForward[7][0] * in.y 
				+ (*s).rotateForward[8][0] * in.z};
	}
	//determine if the current grid location is close to the boundary
	//5 means its a normal point
	// 0..3 first points in the grid
	// -4..-1 last points in grid
	deviceFunction static inline int64_t boundaryFinder(
		const int64_t i, 
		const int64_t N) {
		//normal points
		if (i < (N - 5) && i > 3) return 5;
		//LHS
		else if (i < 4) return i;
		//RHS
		else return i - N;
	}

	[[maybe_unused]] deviceFunction static inline float firstDerivativeSixthOrder(
		const float M3, 
		const float M2, 
		const float M1, 
		const float P1, 
		const float P2, 
		const float P3) {
		return -0.0046875f * M3 
			+ 0.065104166666666667f * M2 
			- 1.171875f * M1 
			+ 1.171875f * P1 
			- 0.065104166666666667f * P2 
			+ 0.0046875f * P3;
	}
	[[maybe_unused]] deviceFunction static inline double firstDerivativeSixthOrder(
		const double M3, 
		const double M2, 
		const double M1, 
		const double P1, 
		const double P2, 
		const double P3) {
		return -0.0046875 * M3 
			+ 0.065104166666666667 * M2 
			- 1.171875 * M1 
			+ 1.171875 * P1 
			- 0.065104166666666667 * P2 
			+ 0.0046875 * P3;
	}


	//apply special boundary conditions to the end of a FDTD grid
	//they are derived from the finite difference generator of Fornberg
	//and depend only on interior points. They thus allow a "one-way" boundary
	//that has low reflection.
	// cases 0..2: solve for dE/dz at z[0..2] for mapping onto the E grid locations
	//   for the other end of the grid, reverse the order of the inputs and negate.
	// cases 11..13: solve for dH/dz but substitute the associated E into the inputs.
	//   the resulting quantity should additionally be divided by Zo.
	[[maybe_unused]] deviceFunction static inline double resolveMaxwellEndpoints(
		const int64_t type, 
		const double in00, 
		const double in01, 
		const double in02, 
		const double in03, 
		const double in04, 
		const double in05, 
		const double in06, 
		const double in07, 
		const double in08, 
		const double in09, 
		const double in10) {
		switch (type) {
		case 0: //reverse and negate for +z end
			return (-2.928968253968253 * in00
				+ 10. * in01
				- 22.5 * in02
				+ 40. * in03
				- 52.5 * in04
				+ 50.4 * in05
				- 35 * in06
				+ 17.142857142857142 * in07
				- 5.625 * in08
				+ 1.1111111111111111 * in09
				- 0.1 * in10);
		case 1:
			return (-0.7517466711619544 * in00
				- 0.4695846315414186 * in01
				+ 4.22831798735119 * in02
				- 7.892969912574405 * in03
				+ 10.470316569010416 * in04
				- 10.08553466796875 * in05
				+ 7.012410481770833 * in06
				- 3.436110723586309 * in07
				+ 1.127578880673363 * in08
				- 0.22271970718625989 * in09
				+ 0.02004239521329365 * in10);
		case 2:
			return (0.02004239521329365 * in00
				- 0.9722130185081846 * in01
				+ 0.6327471051897322 * in02
				+ 0.9213227771577381 * in03
				- 1.2789794921875 * in04
				+ 1.21072998046875 * in05
				- 0.8259480794270833 * in06
				+ 0.3984200613839286 * in07
				- 0.12911551339285712 * in08
				+ 0.025247143942212297 * in09
				- 0.0022533598400297614 * in10);

			//field derivative onto field grid
		case 11:
			return (-0.1 * in00
				- 1.8289682539682544 * in01
				+ 4.5 * in02
				- 6. * in03
				+ 7. * in04
				- 6.3 * in05
				+ 4.2 * in06
				- 2. * in07
				+ 0.6428571428571428 * in08
				- 0.125 * in09
				+ 0.011111111111111111 * in10);
		case 12:
			return (0.01111111111111111 * in00
				- 0.2222222222222222 * in01
				- 1.2178571428571427 * in02
				+ 2.666666666666667 * in03
				- 2.3333333333333333 * in04
				+ 1.8666666666666667 * in05
				- 1.1666666666666667 * in06
				+ 0.5333333333333333 * in07
				- 0.166666666666666667 * in08
				+ 0.031746031746031744 * in09
				- 0.002777777777777778 * in10);
		case 13:
			return (-2.7777777777777778e-03 * in00
				+ 4.16666666666666667e-02 * in01
				- 3.7500000000000000e-01 * in02
				- 7.5952380952380949e-01 * in03
				+ 1.7500000000000000e+00 * in04
				- 1.0500000000000000e+00 * in05
				+ 5.8333333333333333e-01 * in06
				- 2.5000000000000000e-01 * in07
				+ 7.5e-02 * in08
				- 1.38888888888888889e-02 * in09
				+ 1.1904761904761906e-03 * in10);
		}
		return 0.0;
	}

	[[maybe_unused]] deviceFunction static inline float resolveMaxwellEndpoints(
		const int64_t type, 
		const float in00, 
		const float in01, 
		const float in02, 
		const float in03, 
		const float in04, 
		const float in05, 
		const float in06, 
		const float in07, 
		const float in08, 
		const float in09, 
		const float in10) {

		switch (type) {
		case 0: //reverse and negate for +z end
			return (-2.928968253968253f * in00
				+ 10.f * in01
				- 22.5f * in02
				+ 40.f * in03
				- 52.5f * in04
				+ 50.4f * in05
				- 35.0f * in06
				+ 17.142857142857142f * in07
				- 5.625f * in08
				+ 1.1111111111111111f * in09
				- 0.1f * in10);
		case 1:
			return (-0.7517466711619544f * in00
				- 0.4695846315414186f * in01
				+ 4.22831798735119f * in02
				- 7.892969912574405f * in03
				+ 10.470316569010416f * in04
				- 10.08553466796875f * in05
				+ 7.012410481770833f * in06
				- 3.436110723586309f * in07
				+ 1.127578880673363f * in08
				- 0.22271970718625989f * in09
				+ 0.02004239521329365f * in10);
		case 2:
			return (0.02004239521329365f * in00
				- 0.9722130185081846f * in01
				+ 0.6327471051897322f * in02
				+ 0.9213227771577381f * in03
				- 1.2789794921875f * in04
				+ 1.21072998046875f * in05
				- 0.8259480794270833f * in06
				+ 0.3984200613839286f * in07
				- 0.12911551339285712f * in08
				+ 0.025247143942212297f * in09
				- 0.0022533598400297614f * in10);

			//field derivative onto field grid
		case 11:
			return (-0.1f * in00
				- 1.8289682539682544f * in01
				+ 4.5f * in02
				- 6.f * in03
				+ 7.f * in04
				- 6.3f * in05
				+ 4.2f * in06
				- 2.f * in07
				+ 0.6428571428571428f * in08
				- 0.125f * in09
				+ 0.011111111111111111f * in10);
		case 12:
			return (0.01111111111111111f * in00
				- 0.2222222222222222f * in01
				- 1.2178571428571427f * in02
				+ 2.666666666666667f * in03
				- 2.3333333333333333f * in04
				+ 1.8666666666666667f * in05
				- 1.1666666666666667f * in06
				+ 0.5333333333333333f * in07
				- 0.166666666666666667f * in08
				+ 0.031746031746031744f * in09
				- 0.002777777777777778f * in10);
		case 13:
			return (-2.7777777777777778e-03f * in00
				+ 4.16666666666666667e-02f * in01
				- 3.7500000000000000e-01f * in02
				- 7.5952380952380949e-01f * in03
				+ 1.7500000000000000e+00f * in04
				- 1.0500000000000000e+00f * in05
				+ 5.8333333333333326e-01f * in06
				- 2.5000000000000000e-01f * in07
				+ 7.5e-02f * in08
				- 1.38888888888888889e-02f * in09
				+ 1.1904761904761906e-03f * in10);
		}
		return 0.0f;
	}

	template<typename real_t, typename complex_t>
	deviceFunction static void fourierInterpolation(
		const real_t timeShift, 
		const complex_t* fourierData1, 
		const complex_t* fourierData2, 
		const int64_t dataSize, 
		const real_t dOmega, 
		const int64_t frequencyLimit, 
		maxwellKPoint<real_t>& k) {

		real_t result1{};
		real_t result2{};
		const real_t dOmegaT = dOmega * timeShift;
		for (int i = 1; i < frequencyLimit; i++) {
			const real_t w = static_cast<real_t>(i) * dOmegaT;
			const real_t sinw = deviceFPLib::sin(w);
			const real_t cosw = deviceFPLib::cos(w);
			result1 += static_cast<real_t>(i) * (cosw * fourierData1[i].imag() + sinw * fourierData1[i].real());
			result2 += static_cast<real_t>(i) * (cosw * fourierData2[i].imag() + sinw * fourierData2[i].real());
		}
		
		const real_t normalizationFactor = (460.76f * dOmega) / (static_cast<real_t>(dataSize + 1));
		result1 *= normalizationFactor;
		result2 *= normalizationFactor;
		
		k.kE.x += result1;
		k.kE.y += result2;
		
		k.kH.x += inverseZo<real_t>() * result2;
		k.kH.y -= inverseZo<real_t>() * result1;
	}

	deviceFunction static maxwellKPoint<deviceFP> maxwellDerivativeTerms(
		const maxwell3D* s,
		const int64_t i,
		const maxwellPoint<deviceFP>* EgridIn,
		const maxwellPoint<deviceFP>* HgridIn) {

		const int64_t zIndex = i % s->Nz;
		const int64_t xIndex = (i / s->Nz) % s->Nx;
		const int64_t yIndex = (i / (s->Nz * s->Nx));
		const int64_t zEnd = (s->Nz) * (xIndex + 1) + (s->Nz * s->Nx) * yIndex;
		const int64_t zStart = (s->Nz) * xIndex + (s->Nz * s->Nx) * yIndex;
		const int64_t zyOffset = zIndex + yIndex * (s->Nz * s->Nx);
		deviceFP dExDy{};
		deviceFP dExDz{};
		deviceFP dEyDx{};
		deviceFP dEyDz{};
		deviceFP dEzDx{};
		deviceFP dEzDy{};
		deviceFP dHxDy{};
		deviceFP dHxDz{};
		deviceFP dHyDx{};
		deviceFP dHyDz{};
		deviceFP dHzDx{};
		deviceFP dHzDy{};

		//y-derivatives (only done in 3D mode)
		if (s->Ny > 1) {
			const int64_t zxOffset = zIndex + xIndex * s->Nz;
			const int64_t pageSize = s->Nz * s->Nx;
			switch (boundaryFinder(yIndex, s->Ny)) {
			[[likely]] case 5:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (yIndex - 2) * pageSize].x,
					EgridIn[zxOffset + (yIndex - 1) * pageSize].x,
					EgridIn[zxOffset + (yIndex)*pageSize].x,
					EgridIn[zxOffset + (yIndex + 1) * pageSize].x,
					EgridIn[zxOffset + (yIndex + 2) * pageSize].x,
					EgridIn[zxOffset + (yIndex + 3) * pageSize].x);
				dEzDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (yIndex - 2) * pageSize].z,
					EgridIn[zxOffset + (yIndex - 1) * pageSize].z,
					EgridIn[zxOffset + (yIndex)*pageSize].z,
					EgridIn[zxOffset + (yIndex + 1) * pageSize].z,
					EgridIn[zxOffset + (yIndex + 2) * pageSize].z,
					EgridIn[zxOffset + (yIndex + 3) * pageSize].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (yIndex - 3) * pageSize].x,
					HgridIn[zxOffset + (yIndex - 2) * pageSize].x,
					HgridIn[zxOffset + (yIndex - 1) * pageSize].x,
					HgridIn[zxOffset + (yIndex)*pageSize].x,
					HgridIn[zxOffset + (yIndex + 1) * pageSize].x,
					HgridIn[zxOffset + (yIndex + 2) * pageSize].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (yIndex - 3) * pageSize].z,
					HgridIn[zxOffset + (yIndex - 2) * pageSize].z,
					HgridIn[zxOffset + (yIndex - 1) * pageSize].z,
					HgridIn[zxOffset + (yIndex)*pageSize].z,
					HgridIn[zxOffset + (yIndex + 1) * pageSize].z,
					HgridIn[zxOffset + (yIndex + 2) * pageSize].z);
				break;
			case 0:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 2) * pageSize].x,
					EgridIn[zxOffset + (s->Ny - 1) * pageSize].x,
					EgridIn[zxOffset].x,
					EgridIn[zxOffset + pageSize].x,
					EgridIn[zxOffset + 2 * pageSize].x,
					EgridIn[zxOffset + 3 * pageSize].x);
				dEzDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 2) * pageSize].z,
					EgridIn[zxOffset + (s->Ny - 1) * pageSize].z,
					EgridIn[zxOffset].z,
					EgridIn[zxOffset + pageSize].z,
					EgridIn[zxOffset + 2 * pageSize].z,
					EgridIn[zxOffset + 3 * pageSize].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 3) * pageSize].x,
					HgridIn[zxOffset + (s->Ny - 2) * pageSize].x,
					HgridIn[zxOffset + (s->Ny - 1) * pageSize].x,
					HgridIn[zxOffset].x,
					HgridIn[zxOffset + pageSize].x,
					HgridIn[zxOffset + 2 * pageSize].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 3) * pageSize].z,
					HgridIn[zxOffset + (s->Ny - 2) * pageSize].z,
					HgridIn[zxOffset + (s->Ny - 1) * pageSize].z,
					HgridIn[zxOffset].z,
					HgridIn[zxOffset + pageSize].z,
					HgridIn[zxOffset + 2 * pageSize].z);
				break;
			case 1:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 1) * pageSize].x,
					EgridIn[zxOffset].x,
					EgridIn[zxOffset + 1 * (pageSize)].x,
					EgridIn[zxOffset + 2 * (pageSize)].x,
					EgridIn[zxOffset + 3 * (pageSize)].x,
					EgridIn[zxOffset + 4 * (pageSize)].x);
				dEzDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 1) * pageSize].z,
					EgridIn[zxOffset].z,
					EgridIn[zxOffset + 1 * (pageSize)].z,
					EgridIn[zxOffset + 2 * (pageSize)].z,
					EgridIn[zxOffset + 3 * (pageSize)].z,
					EgridIn[zxOffset + 4 * (pageSize)].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 2) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].x,
					HgridIn[zxOffset].x,
					HgridIn[zxOffset + 1 * (pageSize)].x,
					HgridIn[zxOffset + 2 * (pageSize)].x,
					HgridIn[zxOffset + 3 * (pageSize)].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 2) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].z,
					HgridIn[zxOffset].z,
					HgridIn[zxOffset + 1 * (pageSize)].z,
					HgridIn[zxOffset + 2 * (pageSize)].z,
					HgridIn[zxOffset + 3 * (pageSize)].z);
				break;
			case 2:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset].x,
					EgridIn[zxOffset + 1 * pageSize].x,
					EgridIn[zxOffset + 2 * (pageSize)].x,
					EgridIn[zxOffset + 3 * (pageSize)].x,
					EgridIn[zxOffset + 4 * (pageSize)].x,
					EgridIn[zxOffset + 5 * (pageSize)].x);
				dEzDx = firstDerivativeSixthOrder(
					EgridIn[zxOffset].z,
					EgridIn[zxOffset + pageSize].z,
					EgridIn[zxOffset + 2 * (pageSize)].z,
					EgridIn[zxOffset + 3 * (pageSize)].z,
					EgridIn[zxOffset + 4 * (pageSize)].z,
					EgridIn[zxOffset + 5 * (pageSize)].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].x,
					HgridIn[zxOffset].x,
					HgridIn[zxOffset + pageSize].x,
					HgridIn[zxOffset + 2 * (pageSize)].x,
					HgridIn[zxOffset + 3 * (pageSize)].x,
					HgridIn[zxOffset + 4 * (pageSize)].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].z,
					HgridIn[zxOffset].z,
					HgridIn[zxOffset + pageSize].z,
					HgridIn[zxOffset + 2 * (pageSize)].z,
					HgridIn[zxOffset + 3 * (pageSize)].z,
					HgridIn[zxOffset + 4 * (pageSize)].z);
				break;
			case 3:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (yIndex - 2) * pageSize].x,
					EgridIn[zxOffset + (yIndex - 1) * pageSize].x,
					EgridIn[zxOffset + (yIndex)*pageSize].x,
					EgridIn[zxOffset + (yIndex + 1) * pageSize].x,
					EgridIn[zxOffset + (yIndex + 2) * pageSize].x,
					EgridIn[zxOffset + (yIndex + 3) * pageSize].x);
				dEzDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (yIndex - 2) * pageSize].z,
					EgridIn[zxOffset + (yIndex - 1) * pageSize].z,
					EgridIn[zxOffset + (yIndex)*pageSize].z,
					EgridIn[zxOffset + (yIndex + 1) * pageSize].z,
					EgridIn[zxOffset + (yIndex + 2) * pageSize].z,
					EgridIn[zxOffset + (yIndex + 3) * pageSize].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (yIndex - 3) * pageSize].x,
					HgridIn[zxOffset + (yIndex - 2) * pageSize].x,
					HgridIn[zxOffset + (yIndex - 1) * pageSize].x,
					HgridIn[zxOffset + (yIndex)*pageSize].x,
					HgridIn[zxOffset + (yIndex + 1) * pageSize].x,
					HgridIn[zxOffset + (yIndex + 2) * pageSize].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (yIndex - 3) * pageSize].z,
					HgridIn[zxOffset + (yIndex - 2) * pageSize].z,
					HgridIn[zxOffset + (yIndex - 1) * pageSize].z,
					HgridIn[zxOffset + (yIndex)*pageSize].z,
					HgridIn[zxOffset + (yIndex + 1) * pageSize].z,
					HgridIn[zxOffset + (yIndex + 2) * pageSize].z);
				break;
			case -1:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 3) * (pageSize)].x,
					EgridIn[zxOffset + (s->Ny - 2) * (pageSize)].x,
					EgridIn[zxOffset + (s->Ny - 1) * (pageSize)].x,
					EgridIn[zxOffset].x,
					EgridIn[zxOffset + (pageSize)].x,
					EgridIn[zxOffset + 2 * (pageSize)].x);
				dEzDx = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 3) * (pageSize)].z,
					EgridIn[zxOffset + (s->Ny - 2) * (pageSize)].z,
					EgridIn[zxOffset + (s->Ny - 1) * (pageSize)].z,
					EgridIn[zxOffset].z,
					EgridIn[zxOffset + (pageSize)].z,
					EgridIn[zxOffset + 2 * (pageSize)].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 4) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 3) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 2) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].x,
					HgridIn[zxOffset].x,
					HgridIn[zxOffset + (pageSize)].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 4) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 3) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 2) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].z,
					HgridIn[zxOffset].z,
					HgridIn[zxOffset + (pageSize)].z);
				break;
			case -2:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 4) * (pageSize)].x,
					EgridIn[zxOffset + (s->Ny - 3) * (pageSize)].x,
					EgridIn[zxOffset + (s->Ny - 2) * (pageSize)].x,
					EgridIn[zxOffset + (s->Ny - 1) * (pageSize)].x,
					EgridIn[zxOffset].x,
					EgridIn[zxOffset + (pageSize)].x);
				dEzDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 4) * (pageSize)].z,
					EgridIn[zxOffset + (s->Ny - 3) * (pageSize)].z,
					EgridIn[zxOffset + (s->Ny - 2) * (pageSize)].z,
					EgridIn[zxOffset + (s->Ny - 1) * (pageSize)].z,
					EgridIn[zxOffset].z,
					EgridIn[zxOffset + (pageSize)].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 5) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 4) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 3) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 2) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].x,
					HgridIn[zxOffset].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 5) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 4) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 3) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 2) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].z,
					HgridIn[zxOffset].z);
				break;
			case -3:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 5) * (pageSize)].x,
					EgridIn[zxOffset + (s->Ny - 4) * (pageSize)].x,
					EgridIn[zxOffset + (s->Ny - 3) * (pageSize)].x,
					EgridIn[zxOffset + (s->Ny - 2) * (pageSize)].x,
					EgridIn[zxOffset + (s->Ny - 1) * (pageSize)].x,
					EgridIn[zxOffset].x);
				dEzDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (s->Ny - 5) * (pageSize)].z,
					EgridIn[zxOffset + (s->Ny - 4) * (pageSize)].z,
					EgridIn[zxOffset + (s->Ny - 3) * (pageSize)].z,
					EgridIn[zxOffset + (s->Ny - 2) * (pageSize)].z,
					EgridIn[zxOffset + (s->Ny - 1) * (pageSize)].z,
					EgridIn[zxOffset].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 6) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 5) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 4) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 3) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 2) * (pageSize)].x,
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (s->Ny - 6) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 5) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 4) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 3) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 2) * (pageSize)].z,
					HgridIn[zxOffset + (s->Ny - 1) * (pageSize)].z);
				break;
			case -4:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (yIndex - 2) * pageSize].x,
					EgridIn[zxOffset + (yIndex - 1) * pageSize].x,
					EgridIn[zxOffset + (yIndex)*pageSize].x,
					EgridIn[zxOffset + (yIndex + 1) * pageSize].x,
					EgridIn[zxOffset + (yIndex + 2) * pageSize].x,
					EgridIn[zxOffset + (yIndex + 3) * pageSize].x);
				dEzDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (yIndex - 2) * pageSize].z,
					EgridIn[zxOffset + (yIndex - 1) * pageSize].z,
					EgridIn[zxOffset + (yIndex)*pageSize].z,
					EgridIn[zxOffset + (yIndex + 1) * pageSize].z,
					EgridIn[zxOffset + (yIndex + 2) * pageSize].z,
					EgridIn[zxOffset + (yIndex + 3) * pageSize].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (yIndex - 3) * pageSize].x,
					HgridIn[zxOffset + (yIndex - 2) * pageSize].x,
					HgridIn[zxOffset + (yIndex - 1) * pageSize].x,
					HgridIn[zxOffset + (yIndex)*pageSize].x,
					HgridIn[zxOffset + (yIndex + 1) * pageSize].x,
					HgridIn[zxOffset + (yIndex + 2) * pageSize].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (yIndex - 3) * pageSize].z,
					HgridIn[zxOffset + (yIndex - 2) * pageSize].z,
					HgridIn[zxOffset + (yIndex - 1) * pageSize].z,
					HgridIn[zxOffset + (yIndex)*pageSize].z,
					HgridIn[zxOffset + (yIndex + 1) * pageSize].z,
					HgridIn[zxOffset + (yIndex + 2) * pageSize].z);
				break;
			case -5:
				dExDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (yIndex - 2) * pageSize].x,
					EgridIn[zxOffset + (yIndex - 1) * pageSize].x,
					EgridIn[zxOffset + (yIndex)*pageSize].x,
					EgridIn[zxOffset + (yIndex + 1) * pageSize].x,
					EgridIn[zxOffset + (yIndex + 2) * pageSize].x,
					EgridIn[zxOffset + (yIndex + 3) * pageSize].x);
				dEzDy = firstDerivativeSixthOrder(
					EgridIn[zxOffset + (yIndex - 2) * pageSize].z,
					EgridIn[zxOffset + (yIndex - 1) * pageSize].z,
					EgridIn[zxOffset + (yIndex)*pageSize].z,
					EgridIn[zxOffset + (yIndex + 1) * pageSize].z,
					EgridIn[zxOffset + (yIndex + 2) * pageSize].z,
					EgridIn[zxOffset + (yIndex + 3) * pageSize].z);
				dHxDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (yIndex - 3) * pageSize].x,
					HgridIn[zxOffset + (yIndex - 2) * pageSize].x,
					HgridIn[zxOffset + (yIndex - 1) * pageSize].x,
					HgridIn[zxOffset + (yIndex)*pageSize].x,
					HgridIn[zxOffset + (yIndex + 1) * pageSize].x,
					HgridIn[zxOffset + (yIndex + 2) * pageSize].x);
				dHzDy = firstDerivativeSixthOrder(
					HgridIn[zxOffset + (yIndex - 3) * pageSize].z,
					HgridIn[zxOffset + (yIndex - 2) * pageSize].z,
					HgridIn[zxOffset + (yIndex - 1) * pageSize].z,
					HgridIn[zxOffset + (yIndex)*pageSize].z,
					HgridIn[zxOffset + (yIndex + 1) * pageSize].z,
					HgridIn[zxOffset + (yIndex + 2) * pageSize].z);
				break;
			default:
				break;
			}
		}

		//x-derivatives
		switch (boundaryFinder(xIndex, s->Nx)) {
		[[likely]] case 5:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (xIndex - 2) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex - 1) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 1) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 2) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 3) * (s->Nz)].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (xIndex - 2) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex - 1) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 1) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 2) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 3) * (s->Nz)].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (xIndex - 3) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex - 2) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex - 1) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex + 1) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex + 2) * (s->Nz)].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (xIndex - 3) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex - 2) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex - 1) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex + 1) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex + 2) * (s->Nz)].z);
			break;
		case 0:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				EgridIn[zyOffset].y,
				EgridIn[zyOffset + (s->Nz)].y,
				EgridIn[zyOffset + 2 * (s->Nz)].y,
				EgridIn[zyOffset + 3 * (s->Nz)].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				EgridIn[zyOffset].z,
				EgridIn[zyOffset + (s->Nz)].z,
				EgridIn[zyOffset + 2 * (s->Nz)].z,
				EgridIn[zyOffset + 3 * (s->Nz)].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				HgridIn[zyOffset + 0].y,
				HgridIn[zyOffset + (s->Nz)].y,
				HgridIn[zyOffset + 2 * (s->Nz)].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				HgridIn[zyOffset + 0].z,
				HgridIn[zyOffset + (s->Nz)].z,
				HgridIn[zyOffset + 2 * (s->Nz)].z);
			break;
		case 1:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				EgridIn[zyOffset].y,
				EgridIn[zyOffset + (s->Nz)].y,
				EgridIn[zyOffset + 2 * (s->Nz)].y,
				EgridIn[zyOffset + 3 * (s->Nz)].y,
				EgridIn[zyOffset + 4 * (s->Nz)].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				EgridIn[zyOffset].z,
				EgridIn[zyOffset + (s->Nz)].z,
				EgridIn[zyOffset + 2 * (s->Nz)].z,
				EgridIn[zyOffset + 3 * (s->Nz)].z,
				EgridIn[zyOffset + 4 * (s->Nz)].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				HgridIn[zyOffset].y,
				HgridIn[zyOffset + (s->Nz)].y,
				HgridIn[zyOffset + 2 * (s->Nz)].y,
				HgridIn[zyOffset + 3 * (s->Nz)].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				HgridIn[zyOffset].z,
				HgridIn[zyOffset + (s->Nz)].z,
				HgridIn[zyOffset + 2 * (s->Nz)].z,
				HgridIn[zyOffset + 3 * (s->Nz)].z);
			break;
		case 2:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset].y,
				EgridIn[zyOffset + (s->Nz)].y,
				EgridIn[zyOffset + 2 * (s->Nz)].y,
				EgridIn[zyOffset + 3 * (s->Nz)].y,
				EgridIn[zyOffset + 4 * (s->Nz)].y,
				EgridIn[zyOffset + 5 * (s->Nz)].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset].z,
				EgridIn[zyOffset + (s->Nz)].z,
				EgridIn[zyOffset + 2 * (s->Nz)].z,
				EgridIn[zyOffset + 3 * (s->Nz)].z,
				EgridIn[zyOffset + 4 * (s->Nz)].z,
				EgridIn[zyOffset + 5 * (s->Nz)].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				HgridIn[zyOffset].y,
				HgridIn[zyOffset + (s->Nz)].y,
				HgridIn[zyOffset + 2 * (s->Nz)].y,
				HgridIn[zyOffset + 3 * (s->Nz)].y,
				HgridIn[zyOffset + 4 * (s->Nz)].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				HgridIn[zyOffset].z,
				HgridIn[zyOffset + (s->Nz)].z,
				HgridIn[zyOffset + 2 * (s->Nz)].z,
				HgridIn[zyOffset + 3 * (s->Nz)].z,
				HgridIn[zyOffset + 4 * (s->Nz)].z);
			break;
		case 3:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (xIndex - 2) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex - 1) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 1) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 2) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 3) * (s->Nz)].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (xIndex - 2) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex - 1) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 1) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 2) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 3) * (s->Nz)].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (xIndex - 3) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex - 2) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex - 1) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex + 1) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex + 2) * (s->Nz)].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (xIndex - 3) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex - 2) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex - 1) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex + 1) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex + 2) * (s->Nz)].z);
			break;
		case -1:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				EgridIn[zyOffset].y,
				EgridIn[zyOffset + (s->Nz)].y,
				EgridIn[zyOffset + 2 * (s->Nz)].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				EgridIn[zyOffset].z,
				EgridIn[zyOffset + (s->Nz)].z,
				EgridIn[zyOffset + 2 * (s->Nz)].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				HgridIn[zyOffset].y,
				HgridIn[zyOffset + (s->Nz)].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				HgridIn[zyOffset].z,
				HgridIn[zyOffset + (s->Nz)].z);
			break;
		case -2:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				EgridIn[zyOffset].y,
				EgridIn[zyOffset + (s->Nz)].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				EgridIn[zyOffset].z,
				EgridIn[zyOffset + (s->Nz)].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 5) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				HgridIn[zyOffset].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 5) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				HgridIn[zyOffset].z);
			break;
		case -3:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 5) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].y,
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y,
				EgridIn[zyOffset].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (s->Nx - 5) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].z,
				EgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z,
				EgridIn[zyOffset].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 6) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 5) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].y,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (s->Nx - 6) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 5) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 4) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 3) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 2) * (s->Nz)].z,
				HgridIn[zyOffset + (s->Nx - 1) * (s->Nz)].z);
			break;
		case -4:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (xIndex - 2) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex - 1) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 1) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 2) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 3) * (s->Nz)].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (xIndex - 2) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex - 1) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 1) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 2) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 3) * (s->Nz)].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (xIndex - 3) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex - 2) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex - 1) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex + 1) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex + 2) * (s->Nz)].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (xIndex - 3) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex - 2) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex - 1) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex + 1) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex + 2) * (s->Nz)].z);
			break;
		case -5:
			dEyDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (xIndex - 2) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex - 1) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 1) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 2) * (s->Nz)].y,
				EgridIn[zyOffset + (xIndex + 3) * (s->Nz)].y);
			dEzDx = firstDerivativeSixthOrder(
				EgridIn[zyOffset + (xIndex - 2) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex - 1) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 1) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 2) * (s->Nz)].z,
				EgridIn[zyOffset + (xIndex + 3) * (s->Nz)].z);
			dHyDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (xIndex - 3) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex - 2) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex - 1) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex + 1) * (s->Nz)].y,
				HgridIn[zyOffset + (xIndex + 2) * (s->Nz)].y);
			dHzDx = firstDerivativeSixthOrder(
				HgridIn[zyOffset + (xIndex - 3) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex - 2) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex - 1) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex + 1) * (s->Nz)].z,
				HgridIn[zyOffset + (xIndex + 2) * (s->Nz)].z);
			break;
		default:
			break;
		}

		//z-derivatives
		switch (boundaryFinder(zIndex, s->Nz)) {
		[[likely]] case 5:
			dExDz = firstDerivativeSixthOrder(
				EgridIn[i - 3].x,
				EgridIn[i - 2].x,
				EgridIn[i - 1].x,
				EgridIn[i].x,
				EgridIn[i + 1].x,
				EgridIn[i + 2].x);
			dEyDz = firstDerivativeSixthOrder(
				EgridIn[i - 3].y,
				EgridIn[i - 2].y,
				EgridIn[i - 1].y,
				EgridIn[i].y,
				EgridIn[i + 1].y,
				EgridIn[i + 2].y);
			dHxDz = firstDerivativeSixthOrder(
				HgridIn[i - 2].x,
				HgridIn[i - 1].x,
				HgridIn[i].x,
				HgridIn[i + 1].x,
				HgridIn[i + 2].x,
				HgridIn[i + 3].x);
			dHyDz = firstDerivativeSixthOrder(
				HgridIn[i - 2].y,
				HgridIn[i - 1].y,
				HgridIn[i].y,
				HgridIn[i + 1].y,
				HgridIn[i + 2].y,
				HgridIn[i + 3].y);
			break;
		case 0:
			dExDz = resolveMaxwellEndpoints(0,
				EgridIn[zStart].x,
				EgridIn[zStart + 1].x,
				EgridIn[zStart + 2].x,
				EgridIn[zStart + 3].x,
				EgridIn[zStart + 4].x,
				EgridIn[zStart + 5].x,
				EgridIn[zStart + 6].x,
				EgridIn[zStart + 7].x,
				EgridIn[zStart + 8].x,
				EgridIn[zStart + 9].x,
				EgridIn[zStart + 10].x);
			dEyDz = resolveMaxwellEndpoints(0,
				EgridIn[zStart].y,
				EgridIn[zStart + 1].y,
				EgridIn[zStart + 2].y,
				EgridIn[zStart + 3].y,
				EgridIn[zStart + 4].y,
				EgridIn[zStart + 5].y,
				EgridIn[zStart + 6].y,
				EgridIn[zStart + 7].y,
				EgridIn[zStart + 8].y,
				EgridIn[zStart + 9].y,
				EgridIn[zStart + 10].y);
			dHxDz = -inverseZo<deviceFP>() * dEyDz;
			dHyDz = inverseZo<deviceFP>() * dExDz;
			break;
		case 1:
			dExDz = resolveMaxwellEndpoints(1,
				EgridIn[zStart].x,
				EgridIn[zStart + 1].x,
				EgridIn[zStart + 2].x,
				EgridIn[zStart + 3].x,
				EgridIn[zStart + 4].x,
				EgridIn[zStart + 5].x,
				EgridIn[zStart + 6].x,
				EgridIn[zStart + 7].x,
				EgridIn[zStart + 8].x,
				EgridIn[zStart + 9].x,
				EgridIn[zStart + 10].x);
			dEyDz = resolveMaxwellEndpoints(1,
				EgridIn[zStart].y,
				EgridIn[zStart + 1].y,
				EgridIn[zStart + 2].y,
				EgridIn[zStart + 3].y,
				EgridIn[zStart + 4].y,
				EgridIn[zStart + 5].y,
				EgridIn[zStart + 6].y,
				EgridIn[zStart + 7].y,
				EgridIn[zStart + 8].y,
				EgridIn[zStart + 9].y,
				EgridIn[zStart + 10].y);
			dHxDz = -inverseZo<deviceFP>() * resolveMaxwellEndpoints(11,
				EgridIn[zStart].y,
				EgridIn[zStart + 1].y,
				EgridIn[zStart + 2].y,
				EgridIn[zStart + 3].y,
				EgridIn[zStart + 4].y,
				EgridIn[zStart + 5].y,
				EgridIn[zStart + 6].y,
				EgridIn[zStart + 7].y,
				EgridIn[zStart + 8].y,
				EgridIn[zStart + 9].y,
				EgridIn[zStart + 10].y);
			dHyDz = inverseZo<deviceFP>() * resolveMaxwellEndpoints(11,
				EgridIn[zStart].x,
				EgridIn[zStart + 1].x,
				EgridIn[zStart + 2].x,
				EgridIn[zStart + 3].x,
				EgridIn[zStart + 4].x,
				EgridIn[zStart + 5].x,
				EgridIn[zStart + 6].x,
				EgridIn[zStart + 7].x,
				EgridIn[zStart + 8].x,
				EgridIn[zStart + 9].x,
				EgridIn[zStart + 10].x);
			break;
		case 2:
			dExDz = resolveMaxwellEndpoints(2,
				EgridIn[zStart].x,
				EgridIn[zStart + 1].x,
				EgridIn[zStart + 2].x,
				EgridIn[zStart + 3].x,
				EgridIn[zStart + 4].x,
				EgridIn[zStart + 5].x,
				EgridIn[zStart + 6].x,
				EgridIn[zStart + 7].x,
				EgridIn[zStart + 8].x,
				EgridIn[zStart + 9].x,
				EgridIn[zStart + 10].x);
			dEyDz = resolveMaxwellEndpoints(2,
				EgridIn[zStart].y,
				EgridIn[zStart + 1].y,
				EgridIn[zStart + 2].y,
				EgridIn[zStart + 3].y,
				EgridIn[zStart + 4].y,
				EgridIn[zStart + 5].y,
				EgridIn[zStart + 6].y,
				EgridIn[zStart + 7].y,
				EgridIn[zStart + 8].y,
				EgridIn[zStart + 9].y,
				EgridIn[zStart + 10].y);
			dHxDz = -inverseZo<deviceFP>() * resolveMaxwellEndpoints(12,
				EgridIn[zStart].y,
				EgridIn[zStart + 1].y,
				EgridIn[zStart + 2].y,
				EgridIn[zStart + 3].y,
				EgridIn[zStart + 4].y,
				EgridIn[zStart + 5].y,
				EgridIn[zStart + 6].y,
				EgridIn[zStart + 7].y,
				EgridIn[zStart + 8].y,
				EgridIn[zStart + 9].y,
				EgridIn[zStart + 10].y);
			dHyDz = inverseZo<deviceFP>() * resolveMaxwellEndpoints(12,
				EgridIn[zStart].x,
				EgridIn[zStart + 1].x,
				EgridIn[zStart + 2].x,
				EgridIn[zStart + 3].x,
				EgridIn[zStart + 4].x,
				EgridIn[zStart + 5].x,
				EgridIn[zStart + 6].x,
				EgridIn[zStart + 7].x,
				EgridIn[zStart + 8].x,
				EgridIn[zStart + 9].x,
				EgridIn[zStart + 10].x);
			break;
		case 3:
			dExDz = firstDerivativeSixthOrder(
				EgridIn[zStart].x,
				EgridIn[zStart + 1].x,
				EgridIn[zStart + 2].x,
				EgridIn[zStart + 3].x,
				EgridIn[zStart + 4].x,
				EgridIn[zStart + 5].x);
			dEyDz = firstDerivativeSixthOrder(
				EgridIn[zStart].y,
				EgridIn[zStart + 1].y,
				EgridIn[zStart + 2].y,
				EgridIn[zStart + 3].y,
				EgridIn[zStart + 4].y,
				EgridIn[zStart + 5].y);
			dHxDz = -inverseZo<deviceFP>() * resolveMaxwellEndpoints(13,
				EgridIn[zStart].y,
				EgridIn[zStart + 1].y,
				EgridIn[zStart + 2].y,
				EgridIn[zStart + 3].y,
				EgridIn[zStart + 4].y,
				EgridIn[zStart + 5].y,
				EgridIn[zStart + 6].y,
				EgridIn[zStart + 7].y,
				EgridIn[zStart + 8].y,
				EgridIn[zStart + 9].y,
				EgridIn[zStart + 10].y);
			dHyDz = inverseZo<deviceFP>() * resolveMaxwellEndpoints(13,
				EgridIn[zStart].x,
				EgridIn[zStart + 1].x,
				EgridIn[zStart + 2].x,
				EgridIn[zStart + 3].x,
				EgridIn[zStart + 4].x,
				EgridIn[zStart + 5].x,
				EgridIn[zStart + 6].x,
				EgridIn[zStart + 7].x,
				EgridIn[zStart + 8].x,
				EgridIn[zStart + 9].x,
				EgridIn[zStart + 10].x);
			break;
		case -1:
			//note: last point in E grid excluded from calculation, don't touch it!
			dExDz = -resolveMaxwellEndpoints(0,
				EgridIn[zEnd - 2].x,
				EgridIn[zEnd - 3].x,
				EgridIn[zEnd - 4].x,
				EgridIn[zEnd - 5].x,
				EgridIn[zEnd - 6].x,
				EgridIn[zEnd - 7].x,
				EgridIn[zEnd - 8].x,
				EgridIn[zEnd - 9].x,
				EgridIn[zEnd - 10].x,
				EgridIn[zEnd - 11].x,
				EgridIn[zEnd - 12].x);
			dEyDz = -resolveMaxwellEndpoints(0,
				EgridIn[zEnd - 2].y,
				EgridIn[zEnd - 3].y,
				EgridIn[zEnd - 4].y,
				EgridIn[zEnd - 5].y,
				EgridIn[zEnd - 6].y,
				EgridIn[zEnd - 7].y,
				EgridIn[zEnd - 8].y,
				EgridIn[zEnd - 9].y,
				EgridIn[zEnd - 10].y,
				EgridIn[zEnd - 11].y,
				EgridIn[zEnd - 12].y);
			break;
		case -2:
			dExDz = -resolveMaxwellEndpoints(1,
				EgridIn[zEnd - 2].x,
				EgridIn[zEnd - 3].x,
				EgridIn[zEnd - 4].x,
				EgridIn[zEnd - 5].x,
				EgridIn[zEnd - 6].x,
				EgridIn[zEnd - 7].x,
				EgridIn[zEnd - 8].x,
				EgridIn[zEnd - 9].x,
				EgridIn[zEnd - 10].x,
				EgridIn[zEnd - 11].x,
				EgridIn[zEnd - 12].x);
			dEyDz = -resolveMaxwellEndpoints(1,
				EgridIn[zEnd - 2].y,
				EgridIn[zEnd - 3].y,
				EgridIn[zEnd - 4].y,
				EgridIn[zEnd - 5].y,
				EgridIn[zEnd - 6].y,
				EgridIn[zEnd - 7].y,
				EgridIn[zEnd - 8].y,
				EgridIn[zEnd - 9].y,
				EgridIn[zEnd - 10].y,
				EgridIn[zEnd - 11].y,
				EgridIn[zEnd - 12].y);
			dHxDz = -inverseZo<deviceFP>() * resolveMaxwellEndpoints(0,
				EgridIn[zEnd - 2].y,
				EgridIn[zEnd - 3].y,
				EgridIn[zEnd - 4].y,
				EgridIn[zEnd - 5].y,
				EgridIn[zEnd - 6].y,
				EgridIn[zEnd - 7].y,
				EgridIn[zEnd - 8].y,
				EgridIn[zEnd - 9].y,
				EgridIn[zEnd - 10].y,
				EgridIn[zEnd - 11].y,
				EgridIn[zEnd - 12].y);
			dHyDz = inverseZo<deviceFP>() * resolveMaxwellEndpoints(0,
				EgridIn[zEnd - 2].x,
				EgridIn[zEnd - 3].x,
				EgridIn[zEnd - 4].x,
				EgridIn[zEnd - 5].x,
				EgridIn[zEnd - 6].x,
				EgridIn[zEnd - 7].x,
				EgridIn[zEnd - 8].x,
				EgridIn[zEnd - 9].x,
				EgridIn[zEnd - 10].x,
				EgridIn[zEnd - 11].x,
				EgridIn[zEnd - 12].x);
			break;
		case -3:
			dExDz = -resolveMaxwellEndpoints(2,
				EgridIn[zEnd - 2].x,
				EgridIn[zEnd - 3].x,
				EgridIn[zEnd - 4].x,
				EgridIn[zEnd - 5].x,
				EgridIn[zEnd - 6].x,
				EgridIn[zEnd - 7].x,
				EgridIn[zEnd - 8].x,
				EgridIn[zEnd - 9].x,
				EgridIn[zEnd - 10].x,
				EgridIn[zEnd - 11].x,
				EgridIn[zEnd - 12].x);
			dEyDz = -resolveMaxwellEndpoints(2,
				EgridIn[zEnd - 2].y,
				EgridIn[zEnd - 3].y,
				EgridIn[zEnd - 4].y,
				EgridIn[zEnd - 5].y,
				EgridIn[zEnd - 6].y,
				EgridIn[zEnd - 7].y,
				EgridIn[zEnd - 8].y,
				EgridIn[zEnd - 9].y,
				EgridIn[zEnd - 10].y,
				EgridIn[zEnd - 11].y,
				EgridIn[zEnd - 12].y);
			dHxDz = -inverseZo<deviceFP>() * resolveMaxwellEndpoints(11,
				EgridIn[zEnd - 2].y,
				EgridIn[zEnd - 3].y,
				EgridIn[zEnd - 4].y,
				EgridIn[zEnd - 5].y,
				EgridIn[zEnd - 6].y,
				EgridIn[zEnd - 7].y,
				EgridIn[zEnd - 8].y,
				EgridIn[zEnd - 9].y,
				EgridIn[zEnd - 10].y,
				EgridIn[zEnd - 11].y,
				EgridIn[zEnd - 12].y);
			dHyDz = inverseZo<deviceFP>() * resolveMaxwellEndpoints(11,
				EgridIn[zEnd - 2].x,
				EgridIn[zEnd - 3].x,
				EgridIn[zEnd - 4].x,
				EgridIn[zEnd - 5].x,
				EgridIn[zEnd - 6].x,
				EgridIn[zEnd - 7].x,
				EgridIn[zEnd - 8].x,
				EgridIn[zEnd - 9].x,
				EgridIn[zEnd - 10].x,
				EgridIn[zEnd - 11].x,
				EgridIn[zEnd - 12].x);
			break;
		case -4:
			dExDz = firstDerivativeSixthOrder(
				EgridIn[i - 3].x,
				EgridIn[i - 2].x,
				EgridIn[i - 1].x,
				EgridIn[i].x,
				EgridIn[i + 1].x,
				EgridIn[i + 2].x);
			dEyDz = firstDerivativeSixthOrder(
				EgridIn[i - 3].y,
				EgridIn[i - 2].y,
				EgridIn[i - 1].y,
				EgridIn[i].y,
				EgridIn[i + 1].y,
				EgridIn[i + 2].y);
			dHxDz = -inverseZo<deviceFP>() * resolveMaxwellEndpoints(12,
				EgridIn[zEnd - 2].y,
				EgridIn[zEnd - 3].y,
				EgridIn[zEnd - 4].y,
				EgridIn[zEnd - 5].y,
				EgridIn[zEnd - 6].y,
				EgridIn[zEnd - 7].y,
				EgridIn[zEnd - 8].y,
				EgridIn[zEnd - 9].y,
				EgridIn[zEnd - 10].y,
				EgridIn[zEnd - 11].y,
				EgridIn[zEnd - 12].y);
			dHyDz = inverseZo<deviceFP>() * resolveMaxwellEndpoints(12,
				EgridIn[zEnd - 2].x,
				EgridIn[zEnd - 3].x,
				EgridIn[zEnd - 4].x,
				EgridIn[zEnd - 5].x,
				EgridIn[zEnd - 6].x,
				EgridIn[zEnd - 7].x,
				EgridIn[zEnd - 8].x,
				EgridIn[zEnd - 9].x,
				EgridIn[zEnd - 10].x,
				EgridIn[zEnd - 11].x,
				EgridIn[zEnd - 12].x);
			break;
		case -5:
			dExDz = firstDerivativeSixthOrder(
				EgridIn[i - 3].x,
				EgridIn[i - 2].x,
				EgridIn[i - 1].x,
				EgridIn[i].x,
				EgridIn[i + 1].x,
				EgridIn[i + 2].x);
			dEyDz = firstDerivativeSixthOrder(
				EgridIn[i - 3].y,
				EgridIn[i - 2].y,
				EgridIn[i - 1].y,
				EgridIn[i].y,
				EgridIn[i + 1].y,
				EgridIn[i + 2].y);
			dHxDz = -inverseZo<deviceFP>() * resolveMaxwellEndpoints(13,
				EgridIn[zEnd - 2].y,
				EgridIn[zEnd - 3].y,
				EgridIn[zEnd - 4].y,
				EgridIn[zEnd - 5].y,
				EgridIn[zEnd - 6].y,
				EgridIn[zEnd - 7].y,
				EgridIn[zEnd - 8].y,
				EgridIn[zEnd - 9].y,
				EgridIn[zEnd - 10].y,
				EgridIn[zEnd - 11].y,
				EgridIn[zEnd - 12].y);
			dHyDz = inverseZo<deviceFP>() * resolveMaxwellEndpoints(13,
				EgridIn[zEnd - 2].x,
				EgridIn[zEnd - 3].x,
				EgridIn[zEnd - 4].x,
				EgridIn[zEnd - 5].x,
				EgridIn[zEnd - 6].x,
				EgridIn[zEnd - 7].x,
				EgridIn[zEnd - 8].x,
				EgridIn[zEnd - 9].x,
				EgridIn[zEnd - 10].x,
				EgridIn[zEnd - 11].x,
				EgridIn[zEnd - 12].x);
			break;
		default:
			break;
		}
		//return the compound result (structure of two field points, first one for the advancement of E,
		//second for the advancement of H)
		return maxwellKPoint<deviceFP>{
			maxwellPoint<deviceFP>{
				inverseEps0<deviceFP>()* (s->inverseZStep* dHyDz - s->inverseXyStep * dHzDy),
					inverseEps0<deviceFP>()* (s->inverseXyStep* dHzDx - s->inverseZStep * dHxDz),
					inverseEps0<deviceFP>()* s->inverseXyStep* (dHxDy - dHyDx)},
				maxwellPoint<deviceFP>{
				inverseMu0<deviceFP>()* (s->inverseXyStep* dEzDy - s->inverseZStep * dEyDz),
					inverseMu0<deviceFP>()* (-s->inverseXyStep * dEzDx + s->inverseZStep * dExDz),
					inverseMu0<deviceFP>()* s->inverseXyStep* (dEyDx - dExDy)}
		};
	}

	deviceFunction static void maxwellCurrentTerms(
		const maxwell3D* s,
		const int64_t i,
		const int64_t t,
		const bool isAtMidpoint,
		const maxwellPoint<deviceFP>* gridIn,
		const oscillator<deviceFP>* currentGridIn,
		maxwellKPoint<deviceFP>& k,
		const int rkIndex) {
		const int64_t xIndex = i / s->Nz;
		const int64_t zIndex = i - s->Nz * xIndex;
		bool solveMaterialEquations = s->hasMaterialMap ?
			s->materialMap[i] > 0
			: zIndex >= s->materialStart && zIndex < s->materialStop;

		if (solveMaterialEquations) {
			const int64_t oscillatorIndex = s->hasMaterialMap ?
				s->oscillatorIndexMap[i] * s->Noscillators
				: (zIndex - s->materialStart) * s->Noscillators 
				+ xIndex * (s->materialStop - s->materialStart) * s->Noscillators;
			const int oscillatorType = s->hasMaterialMap ?
				s->materialMap[i] - 1
				: 0;
			
			//rotate the field and the currently-active derivative term into the crystal coordinates
			maxwellPoint<deviceFP> crystalField = rotateMaxwellPoint(s, gridIn[i], false);
			maxwellPoint<deviceFP> kE = rotateMaxwellPoint(s, k.kE, false);
			maxwellPoint<deviceFP> chiInstant = s->sellmeierEquations[0][oscillatorType] - 1.0f;
			
			//get the total dipole current and polarization of the oscillators
			maxwellPoint<deviceFP> J{};
			maxwellPoint<deviceFP> P = chiInstant * crystalField;
			for (int j = 0; j < (s->Noscillators - s->hasPlasma[oscillatorType]); j++) {
				J += currentGridIn[oscillatorIndex + j].J;
				P += currentGridIn[oscillatorIndex + j].P;
			}
			//update dEdt (kE) with the dipole current and divide by the instantaneous
			//part of the dielectric constant
			kE += J * inverseEps0<deviceFP>();
			kE /= s->sellmeierEquations[0][oscillatorType];

			//Use dEdt to calculate dPdt
			maxwellPoint<deviceFP> dPdt = (kE * chiInstant) * eps0<deviceFP>();
			dPdt += J;
			dPdt *= 2.0f;
			P *= 2.0f;

			//Calculate the nonlinearity in two terms, an instantaneous
			//one that gets added directly to the propagation and
			//a nonlinear driver which gets added to the driving force
			//of the oscillators. This allows the dispersion of the
			//nonlinearities to follow Miller's rule
			maxwellPoint<deviceFP> instNonlin{};
			maxwellPoint<deviceFP> nonlinearDriver{};

			//calculate the chi2 nonlinearity
			if (s->hasChi2[oscillatorType]) {
				nonlinearDriver += s->chi2[0][oscillatorType] * (P.x * P.x);
				instNonlin += s->chi2[0][oscillatorType] * (2.0f * dPdt.x * P.x);

				nonlinearDriver += s->chi2[1][oscillatorType] * (P.y * P.y);
				instNonlin += s->chi2[1][0] * (2.0f * dPdt.y * P.y);

				nonlinearDriver += s->chi2[2][oscillatorType] * (P.z * P.z);
				instNonlin += s->chi2[2][oscillatorType] * (2.0f * P.z * dPdt.z);

				nonlinearDriver += s->chi2[3][oscillatorType] * (P.y * P.z);
				instNonlin += s->chi2[3][oscillatorType] * (P.y * dPdt.z + dPdt.y * P.z);

				nonlinearDriver += s->chi2[4][oscillatorType] * (P.x * P.z);
				instNonlin += s->chi2[4][oscillatorType] * (P.x * dPdt.z + dPdt.x * P.z);

				nonlinearDriver += s->chi2[5][oscillatorType] * (P.x * P.y);
				instNonlin += s->chi2[5][oscillatorType] * (P.x * dPdt.y + dPdt.x * P.y);
			}

			//calculate kerr nonlinearity for scalar chi3 (assuming centrosymmetry)
			if (s->hasSingleChi3[oscillatorType]) {
				deviceFP fieldSquaredSum = dotProduct(P, P);
				deviceFP dByDtfieldSquaredSum = 2.0f * dotProduct(dPdt, P);
				nonlinearDriver += P * (s->chi3[0][oscillatorType].x * fieldSquaredSum);
				instNonlin += (dPdt * fieldSquaredSum + P * dByDtfieldSquaredSum) * (s->chi3[0][oscillatorType].x);
			}

			//calculate Chi3 nonlinearity with full tensor
			if (s->hasFullChi3[oscillatorType]) {
				for (auto a = 0; a < 3; ++a) {
					for (auto b = 0; b < 3; ++b) {
						for (auto c = 0; c < 3; ++c) {
							nonlinearDriver += s->chi3[a + 3 * b + 9 * c][oscillatorType]
								* P(a) * P(b) * P(c);
							instNonlin += s->chi3[a + 3 * b + 9 * c][oscillatorType] * (
								dPdt(a) * P(b) * P(c)
								+ P(a) * dPdt(b) * P(c)
								+ P(a) * P(b) * dPdt(c));
						}
					}
				}
			}
			nonlinearDriver *= 2.0f;
			instNonlin /= s->sellmeierEquations[0][oscillatorType];
			kE += (instNonlin * chiInstant) * (2.0f * inverseEps0<deviceFP>());

			//resolve the plasma nonlinearity
			deviceFP absorptionCurrent = (s->hasPlasma[oscillatorType]) ?
				deviceFPLib::pow(
					dotProduct(P, P) * s->kNonlinearAbsorption[oscillatorType],
					s->nonlinearAbsorptionOrder[oscillatorType])
				: 0.0f;
			if (s->hasPlasma[oscillatorType]) {
				maxwellPoint<deviceFP> absorption = -0.25f * inverseEps0<deviceFP>() * absorptionCurrent * crystalField;
				absorption += currentGridIn[oscillatorIndex + s->Noscillators - 1].J * inverseEps0<deviceFP>();
				absorption /= s->sellmeierEquations[0][oscillatorType];
				kE += absorption;
			}

			//rotate the updated dEdt back to the beam coordinates
			k.kE = rotateMaxwellPoint(s, kE, true);

			//Update and advance the material oscillators
			for (int j = 0; j < s->Noscillators; j++) {
				oscillator<deviceFP> kOsc = (j < (s->Noscillators - s->hasPlasma[oscillatorType])) ?
					oscillator<deviceFP>{
					(-eps0<deviceFP>() * kLorentzian<deviceFP>())*
						s->sellmeierEquations[1 + j * 3][oscillatorType] * (crystalField+nonlinearDriver)
						- s->sellmeierEquations[2 + j * 3][oscillatorType] * currentGridIn[oscillatorIndex+j].P
						- s->sellmeierEquations[3 + j * 3][oscillatorType] * currentGridIn[oscillatorIndex+j].J,
						currentGridIn[oscillatorIndex + j].J} :
				oscillator<deviceFP>{
					currentGridIn[oscillatorIndex + j].P.x * s->kDrude[oscillatorType] * crystalField
					- s->gammaDrude[oscillatorType] * currentGridIn[oscillatorIndex + j].J,
					maxwellPoint<deviceFP>{2.0f * absorptionCurrent* dotProduct(P,P)* s->kCarrierGeneration[oscillatorType],
					deviceFP{},
					deviceFP{} } };
				//note that k.P.x is used to store the carrier density

				switch (rkIndex) {
				case 0:
					s->materialGridEstimate[oscillatorIndex + j] = 
						s->materialGrid[oscillatorIndex + j] + kOsc * (s->tStep * 0.5f);
					s->materialGridNext[oscillatorIndex + j] = 
						s->materialGrid[oscillatorIndex + j] + kOsc * (sixth<deviceFP>() * s->tStep);
					break;
				case 1:
					s->materialGridEstimate2[oscillatorIndex + j] = 
						s->materialGrid[oscillatorIndex + j] + kOsc * (s->tStep * 0.5f);
					s->materialGridNext[oscillatorIndex + j] += kOsc * (third<deviceFP>() * s->tStep);
					break;
				case 2:
					s->materialGridEstimate[oscillatorIndex + j] = 
						s->materialGrid[oscillatorIndex + j] + kOsc * (s->tStep);
					s->materialGridNext[oscillatorIndex + j] += kOsc * (third<deviceFP>() * s->tStep);
					break;
				case 3:
					s->materialGrid[oscillatorIndex + j] = 
						s->materialGridNext[oscillatorIndex + j] + kOsc * (sixth<deviceFP>() * s->tStep);
					break;
				}
			}
			return;
		}
		//at the front of the grid, run the injection routine
		else if (zIndex == 0 && (t < s->Ninjection)) {
			deviceFP tCurrent = s->tStep * t + 0.5f * isAtMidpoint * s->tStep;
			deviceComplex* fftDataY = (deviceComplex*)s->inputEyFFT;
			deviceComplex* fftDataX = (deviceComplex*)s->inputExFFT;
			fourierInterpolation(tCurrent, 
				&fftDataX[xIndex * s->fftSize], 
				&fftDataY[xIndex * s->fftSize], 
				s->fftSize, 
				s->omegaStep, 
				s->frequencyLimit, 
				k);
			return;
		}
		return;
	}
}