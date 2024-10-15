#include "LightwaveExplorerTrilingual.h"
#include <dlib/global_optimization.h>

typedef 
maxwellCalculation<deviceFP, maxwellPoint<deviceFP>, maxwellPoint<deviceFP>, oscillator<deviceFP>> 
maxwell3D;

namespace deviceFunctions {
	//strict isnan for complex type
	deviceFunction static inline bool isComplexNaN(const deviceComplex& x){
		return isnan(x.real()) || isnan(x.imag());
	}

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
				2.0f * deviceFPLib::pow(
					dotProduct(P, P) * s->kNonlinearAbsorption[oscillatorType],
					s->nonlinearAbsorptionOrder[oscillatorType])
				: 0.0f;
			if (s->hasPlasma[oscillatorType]) {
				maxwellPoint<deviceFP> absorption = -twoPi<deviceFP>() * absorptionCurrent * crystalField;
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
					maxwellPoint<deviceFP>{absorptionCurrent* dotProduct(P,P)* s->kCarrierGeneration[oscillatorType],
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
		else if (zIndex == 0) {
			deviceFP tCurrent = s->tStep * t + 0.5f * isAtMidpoint * s->tStep;
			deviceComplex* fftDataY = (deviceComplex*)s->inputEyFFT;
			deviceComplex* fftDataX = (deviceComplex*)s->inputExFFT;
			int64_t fftSize = s->NtIO / 2 + 1;
			fourierInterpolation(tCurrent, 
				&fftDataX[xIndex * fftSize], 
				&fftDataY[xIndex * fftSize], 
				fftSize, 
				s->omegaStep, 
				s->frequencyLimit, 
				k);
			return;
		}
		return;
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
	};

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
			na *= na;
			nb *= nb;
			nc *= nc;
			deviceFP cp = deviceFPLib::cos(phi);
			deviceFP sp = deviceFPLib::sin(phi);
			deviceFP ct = deviceFPLib::cos(theta);
			deviceFP st = deviceFPLib::sin(theta);

			deviceFP d = (na == nb) ? 
			deviceFP{} : 
			0.5f * deviceFPLib::atan(deviceFPLib::sin(2*phi) * ct / 
			(((1.0f/nb.real() - 1.0f/nc.real())/(1.0f/na.real() - 1.0f/nb.real())) * st*st - cp*cp * ct*ct + sp*sp));
			deviceFP cd = deviceFPLib::cos(d);
			deviceFP sd = deviceFPLib::sin(d);

			deviceComplex nTop = na * nb * nc;

			*ne = nTop / 
			(na*nb*st*st*cd*cd 
			+ na*nc * deviceFPLib::pow(sd*cp + sp*cd*ct,2) 
			+ nb*nc * deviceFPLib::pow(sd*sp - cd*cp*ct,2));

			*no = nTop / 
			(na*nb*st*st*sd*sd
			+ na*nc * deviceFPLib::pow(sd*sp*ct - cd*cp,2)
			+ nb*nc * deviceFPLib::pow(sd*cp*ct + sp*cd,2));


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
using namespace deviceFunctions;

namespace kernelNamespace{
	//Kernels are written in the form of functors, which can be passed to the three main
	// platforms (CUDA, SYCL, c++) without using macros or much in the way of wrappers.
	// All members are public so that I can do aggregate initliazation without writing
	// constructors for all of them. (Could therefore be slightly simplified by using structs
	// instead of classes, but I think a class makes more sense in this context).
	// A macro is used to name the namespace between different modes to avoid ODR violations.

	//calculate the total energy spectrum of the beam for the 2D modes. Note that for the 
	//cartesian one, it will be treated as a round beam instead of an infinite plane wave 
	//in the transverse direction. Thus, the 2D Cartesian spectra are approximations.
	class totalSpectrumKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t i) const {
			deviceFP beamCenter1{};
			deviceFP beamCenter2{};
			deviceFP beamTotal1{};
			deviceFP beamTotal2{};

			//find beam centers
			if ((*s).isCylindric) {
				beamCenter1 = ((*s).Nspace / 2.0f * (*s).dx) + 0.25f * (*s).dx;
				beamCenter2 = beamCenter1;
			}
			else {
				for (int64_t j{}; j < (*s).Nspace; ++j) {
					deviceFP x = (*s).dx * j;
					deviceFP a = modulusSquared((*s).workspace1[i + j * (*s).Nfreq]);
					beamTotal1 += a;
					beamCenter1 += x * a;
					a = modulusSquared((*s).workspace2[i + j * (*s).Nfreq]);
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
			for (int64_t j{}; j < (*s).Nspace; ++j) {
				deviceFP x = (*s).dx * j;
				beamTotal1 += deviceFPLib::abs(x - beamCenter1) 
					* modulusSquared((*s).workspace1[i + j * (*s).Nfreq]);
				beamTotal2 += deviceFPLib::abs(x - beamCenter2) 
					* modulusSquared((*s).workspace2[i + j * (*s).Nfreq]);
			}

			//put the values into the output spectrum
			(*s).gridPolarizationTime1[i] = beamTotal1;
			(*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
			(*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
		}
	};

	//Calculate the energy spectrum after a 2D propagation assuming that the beam
	//height in the non-resolved direction is == the grid width (i.e. square grid)
	//More quantitative than the mapping to round beams, but rather specific
	// DISABLED IN FAVOR OF ROUND-BEAM APPROXIMATION
	//	int64_t j;
	//	deviceFP beamTotal1 = 0.0f;
	//	deviceFP beamTotal2 = 0.0f;
	//	//Integrate total beam power
	//	beamTotal1 = 0.0f;
	//	beamTotal2 = 0.0f;
	//	for (j = 0; j < (*s).Nspace; ++j) {
	//		beamTotal1 += modulusSquared((*s).workspace1[i + j * (*s).Nfreq]);
	//		beamTotal2 += modulusSquared((*s).workspace2[i + j * (*s).Nfreq]);
	//	}
	//	beamTotal1 *= 2.0f * LIGHTC * eps0<deviceFP>() * (*s).dx * (*s).dx * (*s).Nspace * (*s).dt * (*s).dt;
	//	beamTotal2 *= 2.0f * LIGHTC * eps0<deviceFP>() * (*s).dx * (*s).dx * (*s).Nspace * (*s).dt * (*s).dt;
	//	//put the values into the output spectrum
	//	(*s).gridPolarizationTime1[i] = beamTotal1;
	//	(*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
	//	(*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
	//};

	//Calculate the energy spectrum after a 3D propagation
	class totalSpectrum3DKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(int64_t i) const {
			deviceFP beamTotal1{};
			deviceFP beamTotal2{};
			//Integrate total beam power
			for (int64_t j{}; j < (*s).Nspace * (*s).Nspace2; ++j) {
				beamTotal1 += modulusSquared((*s).workspace1[i + j * (*s).Nfreq]);
				beamTotal2 += modulusSquared((*s).workspace2[i + j * (*s).Nfreq]);
			}

			//put the values into the output spectrum
			(*s).gridPolarizationTime1[i] = beamTotal1;
			(*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
			(*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
		}
	};

	//perform a Hankel transform by direct quadrature
	//the offset radial grid allows the sum to be done with a midpoint method
	//with no numerical effort and the rho=0 point is excluded from the grid
	//this function is slow and order N^2 as it is not used in the core loop.
	//the numerical accuracy of Hankel transforms that I've seen is relatively
	//low due to Gibbs phenomena and I find the FFT-based propagation implemented
	//below better for nonlinear phenomena. I might later use this for linear propagation
	//in sequences however.
	class hankelKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP* in;
		deviceFP* out;
		deviceFunction void operator()(const int64_t i) const {
			int64_t col = i / (*s).Ntime; //spatial coordinate
			deviceFP dk = constProd((deviceFP)(1.0 / vPi<deviceFP>()), 2) / ((*s).dx * (*s).Nspace);
			out[i] = {};
			out[i + (*s).Ngrid] = {};
			deviceFP r0;
			deviceFP J0 = 1.0f;
			deviceFP k0 = col * dk;
			for (int64_t r{}; r < (*s).Nspace; ++r) {
				r0 = rhoInRadialSymmetry((*s).Nspace, r, (*s).dx);
				J0 = r0 * j0Device(r0 * k0);
				out[i] += J0 * in[r * (*s).Ntime + i % (*s).Ntime];
				out[i + (*s).Ngrid] += J0 * in[r * (*s).Ntime + (*s).Ngrid + i % (*s).Ntime];
			}
			out[i] *= (*s).dx;
			out[i + (*s).Ngrid] *= (*s).dx;
		}
	};

	//inverse Hankel transform from the k-space back to the offset spatial grid
	class inverseHankelKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP* in;
		deviceFP* out;
		deviceFunction void operator()(const int64_t i) const {
			int64_t col = i / (*s).Ntime; //spatial coordinate
			deviceFP dk = constProd((deviceFP)(1.0 / vPi<deviceFP>()), 2.0) / ((*s).dx * (*s).Nspace);
			out[i] = {};
			out[i + (*s).Ngrid] = {};
			deviceFP r0 = rhoInRadialSymmetry((*s).Nspace, col, (*s).dx);
			deviceFP J0 = 1.0f;
			deviceFP k0 = col * dk;
			for (int64_t k{}; k < (*s).Nspace; ++k) {
				k0 = k * dk;
				J0 = k0 * j0Device(r0 * k0);
				out[i] += J0 * in[k * (*s).Ntime + i % (*s).Ntime];
				out[i + (*s).Ngrid] += J0 * in[k * (*s).Ntime + (*s).Ngrid + i % (*s).Ntime];
			}
			out[i] *= 0.5f * dk / ((*s).Ntime);
			out[i + (*s).Ngrid] *= 0.5f * dk / ((*s).Ntime);
		}
	};

	//rotate the field around the propagation axis (basis change)
	class rotateFieldKernel {
	public:
		const deviceComplex* Ein1;
		const deviceComplex* Ein2;
		deviceComplex* Eout1;
		deviceComplex* Eout2;
		const deviceFP rotationAngle;
		deviceFunction void operator()(const int64_t i) const {
			Eout1[i] = deviceFPLib::cos(rotationAngle) * Ein1[i] - deviceFPLib::sin(rotationAngle) * Ein2[i];
			Eout2[i] = deviceFPLib::sin(rotationAngle) * Ein1[i] + deviceFPLib::cos(rotationAngle) * Ein2[i];
		}
	};

	//calculate the extra term in the Laplacian encountered in cylindrical coordinates (1/rho d/drho)
	class radialLaplacianKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t j = i / (*s).Ntime;
			
			//zero at edges of grid
			[[unlikely]] if (j < 3 || j > ((*s).Nspace - 4)) {
				(*s).gridRadialLaplacian1[i] = {};
				(*s).gridRadialLaplacian2[i] = {};
			}
			else {
				const int64_t h = i - j * (*s).Ntime;
				const deviceFP* E1 = (*s).gridETime1 + h;
				const deviceFP* E2 = (*s).gridETime2 + h;
				if (j < (*s).Nspace / 2) {
					const deviceFP rhoFac = 2.0f / ((*s).dx * (*s).dx * (static_cast<deviceFP>((*s).Nspace / 2 - j) - 0.25f));
					const int64_t neighbors[6]{
						((*s).Nspace - j - 2) * (*s).Ntime,
						(j + 1) * (*s).Ntime,
						((*s).Nspace - j - 1) * (*s).Ntime,
						((*s).Nspace - j) * (*s).Ntime,
						(j - 1) * (*s).Ntime,
						((*s).Nspace - j + 1) * (*s).Ntime
					};
					(*s).gridRadialLaplacian1[i] = rhoFac
						* (firstDerivativeStencil<deviceFP>(-3) * E1[neighbors[0]]
							+ firstDerivativeStencil<deviceFP>(-2) * E1[neighbors[1]]
							+ firstDerivativeStencil<deviceFP>(-1) * E1[neighbors[2]]
							+ firstDerivativeStencil<deviceFP>(1) * E1[neighbors[3]]
							+ firstDerivativeStencil<deviceFP>(2) * E1[neighbors[4]]
							+ firstDerivativeStencil<deviceFP>(3) * E1[neighbors[5]]);
					(*s).gridRadialLaplacian2[i] = rhoFac
						* (firstDerivativeStencil<deviceFP>(-3) * E2[neighbors[0]]
							+ firstDerivativeStencil<deviceFP>(-2) * E2[neighbors[1]]
							+ firstDerivativeStencil<deviceFP>(-1) * E2[neighbors[2]]
							+ firstDerivativeStencil<deviceFP>(1) * E2[neighbors[3]]
							+ firstDerivativeStencil<deviceFP>(2) * E2[neighbors[4]]
							+ firstDerivativeStencil<deviceFP>(3) * E2[neighbors[5]]);
				}
				else {
					const deviceFP rhoFac = 2.0f / ((*s).dx * (*s).dx * (static_cast<deviceFP>(j - (*s).Nspace / 2) + 0.25f));
					const int64_t neighbors[6]{
						((*s).Nspace - j + 1) * (*s).Ntime,
						(j - 1) * (*s).Ntime,
						((*s).Nspace - j) * (*s).Ntime,
						((*s).Nspace - j - 1) * (*s).Ntime,
						(j + 1) * (*s).Ntime,
						((*s).Nspace - j - 2) * (*s).Ntime
					};
					(*s).gridRadialLaplacian1[i] = rhoFac
						* (firstDerivativeStencil<deviceFP>(-3) * E1[neighbors[0]]
							+ firstDerivativeStencil<deviceFP>(-2) * E1[neighbors[1]]
							+ firstDerivativeStencil<deviceFP>(-1) * E1[neighbors[2]]
							+ firstDerivativeStencil<deviceFP>(1) * E1[neighbors[3]]
							+ firstDerivativeStencil<deviceFP>(2) * E1[neighbors[4]]
							+ firstDerivativeStencil<deviceFP>(3) * E1[neighbors[5]]);
					(*s).gridRadialLaplacian2[i] = rhoFac
						* (firstDerivativeStencil<deviceFP>(-3) * E2[neighbors[0]]
							+ firstDerivativeStencil<deviceFP>(-2) * E2[neighbors[1]]
							+ firstDerivativeStencil<deviceFP>(-1) * E2[neighbors[2]]
							+ firstDerivativeStencil<deviceFP>(1) * E2[neighbors[3]]
							+ firstDerivativeStencil<deviceFP>(2) * E2[neighbors[4]]
							+ firstDerivativeStencil<deviceFP>(3) * E2[neighbors[5]]);
				}
			}
		}
	};

	class apertureFarFieldKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP radius;
		const deviceFP activationParameter;
		const deviceFP xOffset;
		const deviceFP yOffset;
		deviceFunction void operator()(int64_t i) const {
			const int64_t col = i / ((*s).Nfreq - 1); //spatial coordinate
			const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
			i = j + col * (*s).Nfreq;
			const int64_t l = col / (*s).Nspace;
			const int64_t k = col - l * (*s).Nspace;
			

			//magnitude of k vector
			const deviceFP ko = constProd(twoPi<deviceFP>(), 1.0 / lightC<double>()) * j * (*s).fStep;

			//transverse wavevector being resolved
			const deviceFP dk1 = k * (*s).dk1 - (k >= ((int64_t)(*s).Nspace / 2)) 
				* ((*s).dk1 * (int64_t)(*s).Nspace); //frequency grid in x direction
			deviceFP dk2 = {};
			if ((*s).is3D) dk2 = l * (*s).dk2 - (l >= ((int64_t)(*s).Nspace2 / 2)) 
				* ((*s).dk2 * (int64_t)(*s).Nspace2); //frequency grid in y direction

			//light that won't go the the farfield is immediately zero
			if (dk1 * dk1 > ko * ko || dk2 * dk2 > ko * ko) {
				(*s).gridEFrequency1[i] = {};
				(*s).gridEFrequency2[i] = {};
				return;
			}

			deviceFP theta1 = deviceFPLib::asin(dk1 / ko);
			deviceFP theta2 = deviceFPLib::asin(dk2 / ko);

			theta1 -= (!(*s).isCylindric) * xOffset;
			theta2 -= (*s).is3D * yOffset;

			const deviceFP r = deviceFPLib::hypot(theta1, theta2);
			const deviceFP a = 1.0f 
				- (1.0f / (1.0f + deviceFPLib::exp(-activationParameter * (r - radius))));
			(*s).gridEFrequency1[i] *= a;
			(*s).gridEFrequency2[i] *= a;
			if (j == 1) {
				(*s).gridEFrequency1[i - 1] = deviceComplex{};
				(*s).gridEFrequency2[i - 1] = deviceComplex{};
			}
		}
	};

	class inverseApertureFarFieldKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP radius;
		const deviceFP activationParameter;
		const deviceFP xOffset;
		const deviceFP yOffset;
		deviceFunction void operator()(int64_t i) const {
			const int64_t col = i / ((*s).Nfreq - 1); //spatial coordinate
			const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
			i = j + col * (*s).Nfreq;
			const int64_t l = col / (*s).Nspace;
			const int64_t k = col - l * (*s).Nspace;
			

			//magnitude of k vector
			const deviceFP ko = constProd(twoPi<deviceFP>(), 1.0 / lightC<double>()) * j * (*s).fStep;

			//transverse wavevector being resolved
			const deviceFP dk1 = k * (*s).dk1 - (k >= ((int64_t)(*s).Nspace / 2)) 
				* ((*s).dk1 * (int64_t)(*s).Nspace); //frequency grid in x direction
			deviceFP dk2 = {};
			if ((*s).is3D) dk2 = l * (*s).dk2 - (l >= ((int64_t)(*s).Nspace2 / 2)) 
				* ((*s).dk2 * (int64_t)(*s).Nspace2); //frequency grid in y direction

			//light that won't go the the farfield is immediately zero
			if (dk1 * dk1 > ko * ko || dk2 * dk2 > ko * ko) {
				(*s).gridEFrequency1[i] = {};
				(*s).gridEFrequency2[i] = {};
				return;
			}

			deviceFP theta1 = deviceFPLib::asin(dk1 / ko);
			deviceFP theta2 = deviceFPLib::asin(dk2 / ko);

			theta1 -= (!(*s).isCylindric) * xOffset;
			theta2 -= (*s).is3D * yOffset;

			const deviceFP r = deviceFPLib::hypot(theta1, theta2);
			const deviceFP a = (1.0f / (1.0f + deviceFPLib::exp(-activationParameter * (r - radius))));
			(*s).gridEFrequency1[i] *= a;
			(*s).gridEFrequency2[i] *= a;
			if (j == 1) {
				(*s).gridEFrequency1[i - 1] = deviceComplex{};
				(*s).gridEFrequency2[i - 1] = deviceComplex{};
			}
		}
	};

	class apertureFarFieldKernelHankel { 
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP radius;
		const deviceFP activationParameter; 
		const deviceFP xOffset;
		const deviceFP yOffset;
		deviceFunction void operator()(int64_t i) const {
			const int64_t col = i / ((*s).Nfreq - 1); //spatial coordinate
			const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
			i = j + col * (*s).Nfreq;
			const int64_t k = col % (*s).Nspace;

			//magnitude of k vector
			deviceFP ko = constProd(twoPi<deviceFP>(), 1.0 / lightC<double>()) 
				* j * (*s).fStep;

			//transverse wavevector being resolved
			deviceFP dk1 = constProd((deviceFP)2.0, 1.0 / vPi<double>()) 
				* k / ((*s).dx * (*s).Nspace);; //frequency grid in x direction

			//light that won't go the the farfield is immediately zero
			if (dk1 * dk1 > ko * ko) {
				(*s).gridEFrequency1[i] = {};
				(*s).gridEFrequency2[i] = {};
				return;
			}

			const deviceFP theta1 = deviceFPLib::asin(dk1 / ko);
			const deviceFP a = 1.0f - (1.0f / (1.0f + 
				deviceFPLib::exp(-activationParameter * (deviceFPLib::abs(theta1) - radius))));
			(*s).gridEFrequency1[i] *= a;
			(*s).gridEFrequency2[i] *= a;
			if (j == 1) {
				(*s).gridEFrequency1[i - 1] = deviceComplex{};
				(*s).gridEFrequency2[i - 1] = deviceComplex{};
			}
		}
	};

	class inverseApertureFarFieldKernelHankel { 
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP radius;
		const deviceFP activationParameter; 
		const deviceFP xOffset;
		const deviceFP yOffset;
		deviceFunction void operator()(int64_t i) const {
			const int64_t col = i / ((*s).Nfreq - 1); //spatial coordinate
			const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
			i = j + col * (*s).Nfreq;
			const int64_t k = col % (*s).Nspace;

			//magnitude of k vector
			deviceFP ko = constProd(twoPi<deviceFP>(), 1.0 / lightC<double>()) 
				* j * (*s).fStep;

			//transverse wavevector being resolved
			deviceFP dk1 = constProd((deviceFP)2.0, 1.0 / vPi<double>()) 
				* k / ((*s).dx * (*s).Nspace);; //frequency grid in x direction

			//light that won't go the the farfield is immediately zero
			if (dk1 * dk1 > ko * ko) {
				(*s).gridEFrequency1[i] = {};
				(*s).gridEFrequency2[i] = {};
				return;
			}

			const deviceFP theta1 = deviceFPLib::asin(dk1 / ko);
			const deviceFP a = (1.0f / (1.0f + 
				deviceFPLib::exp(-activationParameter * (deviceFPLib::abs(theta1) - radius))));
			(*s).gridEFrequency1[i] *= a;
			(*s).gridEFrequency2[i] *= a;
			if (j == 1) {
				(*s).gridEFrequency1[i - 1] = deviceComplex{};
				(*s).gridEFrequency2[i - 1] = deviceComplex{};
			}
		}
	};


	//apply a spectral filter to the beam (full details in docs)
	class filterKernel { 
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFP f0;
		const deviceFP bandwidth; 
		const int order; const deviceFP inBandAmplitude;
		const deviceFP outOfBandAmplitude;
		deviceFunction void operator()(int64_t i) const {
			const int64_t col = i / ((*s).Nfreq - 1); //spatial coordinate
			const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
			i = j + col * (*s).Nfreq;

			deviceFP f = ((*s).fStep * j - f0) / bandwidth;
			for (int p = 1; p < order; p++) {
				f *= ((*s).fStep * j - f0) / bandwidth;
			}
			const deviceFP filterFunction = 
				outOfBandAmplitude 
				+ inBandAmplitude * deviceFPLib::exp(-0.5f * f);
			(*s).gridEFrequency1[i] *= filterFunction;
			(*s).gridEFrequency2[i] *= filterFunction;
		}
	};

	class applyOpticKernel { 
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceComplex* complexReflectivity;
		const bool applyX;
		const bool applyY;
		deviceFunction void operator()(int64_t i) const {
			const int64_t col = i / ((*s).Nfreq - 1); //spatial coordinate
			const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
			i = j + col * (*s).Nfreq;
			if(applyX) (*s).gridEFrequency1[i] *= complexReflectivity[j]/(*s).Ntime;
			else (*s).gridEFrequency1[i] /= (*s).Ntime;
			if(applyY) (*s).gridEFrequency2[i] *= complexReflectivity[j]/(*s).Ntime;
			else (*s).gridEFrequency2[i] /= (*s).Ntime;
			if(j==1){
				(*s).gridEFrequency1[i-1] = {};
				(*s).gridEFrequency2[i-1] = {};
			}
		}
	};

	//apply a lorentzian gain or loss in a certain cross-section of the beam.
	// amplitude - strength of the copy of the pulse applied
	// f0 - resonance frequency of the Lorentzian (THz)
	// gamma - linewidth (radHz)
	// radius - radius of the spot (m)
	// order - supergaussian order of the spot shape
	class lorentzianSpotKernel { 
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP amplitude;
		const deviceFP f0;
		const deviceFP gamma;
		const deviceFP radius;
		const deviceFP order;
		deviceFunction void operator()(int64_t i) const {
			int64_t col = i / ((*s).Nfreq - 1);
			int64_t j = col % (*s).Nspace;
			int64_t k = col / (*s).Nspace;
			int64_t h = 1 + i % ((*s).Nfreq - 1);
			deviceFP r, f, x, y;
			if ((*s).is3D) {
				col = i / ((*s).Nfreq - 1);
				i = h + col * ((*s).Nfreq);
				j = col % (*s).Nspace;
				k = col / (*s).Nspace;
				f = h * (*s).fStep;
				x = ((*s).dx * (j - (*s).Nspace / 2.0f));
				y = ((*s).dx * (k - (*s).Nspace2 / 2.0f));
				r = deviceFPLib::hypot(x, y);
			}
			else {
				j = i / ((*s).Nfreq - 1);
				i = h + j * ((*s).Nfreq);
				f = h * (*s).fStep;
				r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) + 0.25f * (*s).dx);
			}

			deviceFP w0 = twoPi<deviceFP>() * f0;
			deviceFP w = twoPi<deviceFP>() * f;
			deviceComplex lorentzian = 
				gamma * w0 * amplitude / (w0 * w0 - w * w + deviceComplex(0.0f, gamma * w));
			deviceFP spotFactor = r / radius;
			for (int p = 1; p < (int)order; p++) {
				spotFactor *= spotFactor;
			}
			deviceComplex filterFunction = 
				deviceComplex(0.0f, deviceFPLib::exp(-spotFactor)) * lorentzian;
			(*s).gridEFrequency1[i] += filterFunction * (*s).gridEFrequency1[i];
			(*s).gridEFrequency2[i] += filterFunction * (*s).gridEFrequency2[i];
		}
	};

	//Apply a (soft, possibly) aperture
	class apertureKernel { 
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP radius;
		const deviceFP activationParameter;
		deviceFunction void operator()(int64_t i) const {
			const int64_t col = i / (*s).Ntime;
			const int64_t j = col % (*s).Nspace;
			const int64_t k = col / (*s).Nspace;
			deviceFP r;
			if ((*s).is3D) {
				const deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0f));
				const deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0f));
				r = deviceFPLib::hypot(x, y);
			}
			else {
				r = deviceFPLib::abs((*s).dx 
					* ((deviceFP)j - (*s).Nspace / 2.0f) + 0.25f * (*s).dx);
			}

			const deviceFP a = 1.0f - (1.0f / (1.0f 
				+ deviceFPLib::exp(-activationParameter * (r - radius) / (*s).dx)));
			(*s).gridETime1[i] *= a;
			(*s).gridETime2[i] *= a;
		}
	};

	//apply a spatial phase corresponding to a parabolic mirror (on-axis)
	class parabolicMirrorKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP focus;
		deviceFunction void operator()(int64_t i) const {
			const int64_t h = 1 + i % ((*s).Nfreq - 1);
			const int64_t col = i / ((*s).Nfreq - 1);
			i = h + col * (*s).Nfreq;
			const int64_t j = col % (*s).Nspace;
			const int64_t k = col / (*s).Nspace;

			const deviceFP w = twoPi<deviceFP>() * h * (*s).fStep;
			deviceFP r;
			if ((*s).is3D) {
				const deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0f));
				const deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0f));
				r = deviceFPLib::hypot(x, y);
			}
			else {
				r = deviceFPLib::abs((*s).dx 
					* ((deviceFP)j - (*s).Nspace / 2.0f) + 0.25f * (*s).dx);
			}

			const deviceComplex u = deviceLib::exp(deviceComplex(0.0f,
				w * r * r * (0.5f / focus) / lightC<deviceFP>()));

			(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
			(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
		}
	};

	//apply a spatial phase corresponding to a spherical mirror (on axis)
	class sphericalMirrorKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFP ROC_in;
		deviceFunction void operator()(int64_t i) const {
			const int64_t h = 1 + i % ((*s).Nfreq - 1);
			const int64_t col = i / ((*s).Nfreq - 1);
			i = h + col * (*s).Nfreq;
			const int64_t j = col % (*s).Nspace;
			const int64_t k = col / (*s).Nspace;

			const deviceFP w = twoPi<deviceFP>() * h * (*s).fStep;
			deviceFP r;
			if ((*s).is3D) {
				const deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0f));
				const deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0f));
				r = deviceFPLib::hypot(x, y);
			}
			else {
				r = deviceFPLib::abs((*s).dx 
					* ((deviceFP)j - (*s).Nspace / 2.0f) + 0.25f * (*s).dx);
			}
			const bool isNegative = ROC_in < 0.0f;
			const deviceFP ROC = deviceFPLib::abs(ROC_in);
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
					(-0.5f * ratio - 0.125f * ratio * ratio 
						- 0.0625f * ratio * ratio * ratio) 
					/ lightC<deviceFP>()));
			}


			(*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
			(*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
			if (isComplexNaN((*s).gridEFrequency1[i]) 
				|| isComplexNaN((*s).gridEFrequency2[i])) {
				(*s).gridEFrequency1[i] = deviceComplex{};
				(*s).gridEFrequency2[i] = deviceComplex{};
			}
		}
	};

	//apply linear propagation through a given medium to the fields
	class correctFDTDAmplitudesKernel { 
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t localIndex) const {
			int64_t i = localIndex;
			const int64_t h = 1 + i % ((*s).Nfreq - 1);
			const int64_t col = i / ((*s).Nfreq - 1);
			i = h + col * ((*s).Nfreq);
			const int64_t k = col / (*s).Nspace;
			const int64_t j = col - k * (*s).Nspace;
			
			const deviceFP kMagnitude = (h * (*s).fStep * twoPi<deviceFP>()) / lightC<deviceFP>();
			const deviceFP dk1 = j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace);
			const deviceFP dk2 = k * (*s).dk2 - (k >= ((*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2);
			if(dk2*dk2 + dk1*dk1 < kMagnitude*kMagnitude){
				s->gridEFrequency1[i] *= (*s).fftNorm * kMagnitude/deviceFPLib::sqrt((kMagnitude - dk1)*(kMagnitude + dk1));
				s->gridEFrequency2[i] *= (*s).fftNorm * kMagnitude/deviceFPLib::sqrt((kMagnitude - dk2)*(kMagnitude + dk2));
			}
			else{
				s->gridEFrequency1[i] = {};
				s->gridEFrequency2[i] = {};
			}
			if(h==1){
				s->gridEFrequency1[i-1] = {};
				s->gridEFrequency2[i-1] = {};
			}
		}
	};

		class correctFDTDAmplitudesKernel2D { 
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t localIndex) const {
			int64_t i = localIndex;
			const int64_t h = 1 + i % ((*s).Nfreq - 1);
			const int64_t col = i / ((*s).Nfreq - 1);
			i = h + col * ((*s).Nfreq);
			const int64_t j = col % (*s).Nspace;
			
			const deviceFP kMagnitude = (h * (*s).fStep * twoPi<deviceFP>()) / lightC<deviceFP>();
			const deviceFP dk1 = j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace);
			if(dk1*dk1 < kMagnitude*kMagnitude && h > 2){
				s->gridEFrequency1[i] *= (*s).fftNorm * kMagnitude/deviceFPLib::sqrt((kMagnitude - dk1)*(kMagnitude + dk1));
				s->gridEFrequency2[i] *= (*s).fftNorm;
			}
			else{
				s->gridEFrequency1[i] = {};
				s->gridEFrequency2[i] = {};
			}
			if(h==1){
				s->gridEFrequency1[i-1] = {};
				s->gridEFrequency2[i-1] = {};
			}
		}
	};
	
	//apply linear propagation through a given medium to the fields
	class applyLinearPropagationKernel { 
	public:
		const deviceFP* sellmeierCoefficients;
		const deviceFP thickness;
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t localIndex) const {
			int64_t i = localIndex;
			const int axesNumber = (*s).axesNumber;
			const int sellmeierType = (*s).sellmeierType;
			deviceComplex ne, no, n0, n0o;
			const int64_t h = 1 + i % ((*s).Nfreq - 1);
			const int64_t col = i / ((*s).Nfreq - 1);
			i = h + col * ((*s).Nfreq);
			const int64_t j = col % (*s).Nspace;
			const int64_t k = col / (*s).Nspace;
			const deviceFP crystalTheta = (*s).crystalTheta;
			const deviceFP crystalPhi = (*s).crystalPhi;

			//frequency being resolved by current thread
			const deviceFP f = h * (*s).fStep;
			const deviceFP omega = twoPi<deviceFP>() * f;
			deviceFP delta = findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);
			if(s->axesNumber==2 && i<s->NgridC) s->gridBiaxialDelta[i] = delta;

			const deviceFP dk1 = j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace);
			deviceFP dk2 = k * (*s).dk2 - (k >= ((*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2);
			if (!(*s).is3D)dk2 = 0.0f;
			sellmeierDevice(&n0, &n0o, sellmeierCoefficients, (*s).f0,
				crystalTheta, crystalPhi, axesNumber, sellmeierType);
			if (isComplexNaN(ne) || isComplexNaN(no)) {
				ne = cOne<deviceComplex>();
				no = cOne<deviceComplex>();
			}

			const deviceComplex ke = ne * omega / lightC<deviceFP>();
			const deviceComplex ko = no * omega / lightC<deviceFP>();
			const deviceComplex k0 = n0 * omega / lightC<deviceFP>();

			deviceComplex ts = fourierPropagator(ke, dk1, dk2, k0.real(), thickness);
			deviceComplex tp = fourierPropagator(ko, dk1, dk2, k0.real(), thickness);
			if (((dk1 * dk1 + dk2 * dk2 > modulusSquared(ke))
				|| (dk1 * dk1 + dk2 * dk2 > modulusSquared(ko)))
				&& thickness < 0.0f) {
				ts = deviceComplex{};
				tp = deviceComplex{};
			}
			if (isComplexNaN(ts)) ts = deviceComplex{};
			if (isComplexNaN(tp)) tp = deviceComplex{};
			(*s).gridEFrequency1[i] = ts * (*s).gridEFrequency1[i];
			(*s).gridEFrequency2[i] = tp * (*s).gridEFrequency2[i];
			if (isComplexNaN((*s).gridEFrequency1[i])
			|| isComplexNaN((*s).gridEFrequency2[i])){
				(*s).gridEFrequency1[i] = deviceComplex{};
				(*s).gridEFrequency2[i] = deviceComplex{};
			}
			if (h == 1) {
				(*s).gridEFrequency1[i - 1] = deviceComplex{};
				(*s).gridEFrequency2[i - 1] = deviceComplex{};
			}
		}
	};

	//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
	//note that the sellmeier coefficients have extra values appended to the end
	//to give info about the current simulation
	class prepareCartesianGridsKernel {
	public:
		const deviceFP* sellmeierCoefficients;
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t localIndex) const {
			const int64_t j = localIndex / ((*s).Nfreq - 1); //spatial coordinate
			const int64_t k = 1 + (localIndex % ((*s).Nfreq - 1)); //temporal coordinate
			const int64_t i = k + j * (*s).Nfreq;

			//frequency being resolved by current thread
			const deviceFP f = k * (*s).fStep;

			//transverse wavevector being resolved
			const deviceFP dk = j * (*s).dk1 - (j >= ((*s).Nspace / 2)) 
				* ((*s).dk1 * (*s).Nspace); //frequency grid in transverse direction
			deviceComplex ne, no;
			deviceFP d = findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);
			if(s->axesNumber==2 && i < s->NgridC) s->gridBiaxialDelta[i] = d;
			//if the refractive index was returned weird, 
			//then the index isn't valid, so set the propagator to zero for that frequency
			if (minN(ne.real(), no.real()) < 0.9f 
				|| isComplexNaN(ne) 
				|| isComplexNaN(no)) {
				(*s).gridPropagationFactor1[i] = {};
				(*s).gridPropagationFactor2[i] = {};
				(*s).gridPolarizationFactor1[i] = {};
				(*s).gridPolarizationFactor2[i] = {};
				return;
			}

			const deviceComplex k0 = deviceComplex(
				twoPi<deviceFP>() * (*s).n0.real() * f / lightC<deviceFP>(), 
				0.0f);
			const deviceComplex ke = twoPi<deviceFP>() * ne * f / lightC<deviceFP>();
			const deviceComplex ko = twoPi<deviceFP>() * no * f / lightC<deviceFP>();

			//chi11 factor, also multiplied by -sqrt(-1)!
			const deviceComplex chi11 = ((*s).isUsingMillersRule) ?
				deviceComplex((*s).chiLinear1[k].imag(), -(*s).chiLinear1[k].real())
				: cOne<deviceComplex>();
			const deviceComplex chi12 = ((*s).isUsingMillersRule) ?
				deviceComplex((*s).chiLinear2[k].imag(), -(*s).chiLinear2[k].real())
				: cOne<deviceComplex>();

			const deviceComplex kz1 = deviceLib::sqrt(ke - dk) * deviceLib::sqrt(ke + dk);
			const deviceComplex kz2 = deviceLib::sqrt(ko - dk) * deviceLib::sqrt(ko + dk);

			if (kz1.real() > 0.0f 
			&& kz2.real() > 0.0f
			&& !isComplexNaN(k0)
			&& !isComplexNaN(kz1)
			&& !isComplexNaN(kz2)) {
				(*s).gridPropagationFactor1[i] = fourierPropagator(ke, dk, deviceFP{}, k0.real(), 0.5f * (*s).h);
				(*s).gridPropagationFactor2[i] = fourierPropagator(ko, dk, deviceFP{}, k0.real(), 0.5f * (*s).h);
				(*s).gridPolarizationFactor1[i] = 
					deviceLib::pow((*s).chiLinear1[k] + 1.0f, (deviceFP)0.25f) 
					* chi11 * (twoPi<deviceFP>() * twoPi<deviceFP>() * f * f) 
					/ ((2.0f * lightC<deviceFP>() * lightC<deviceFP>() * kz1)) * (*s).h;
				(*s).gridPolarizationFactor2[i] = 
					deviceLib::pow((*s).chiLinear2[k] + 1.0f, (deviceFP)0.25f) 
					* chi12 * (twoPi<deviceFP>() * twoPi<deviceFP>() * f * f) 
					/ ((2.0f * lightC<deviceFP>() * lightC<deviceFP>() * kz2)) * (*s).h;
			}
			else {
				(*s).gridPropagationFactor1[i] = {};
				(*s).gridPropagationFactor2[i] = {};
				(*s).gridPolarizationFactor1[i] = {};
				(*s).gridPolarizationFactor2[i] = {};
			}

			if(isComplexNaN((*s).gridPropagationFactor1[i])
			|| isComplexNaN((*s).gridPropagationFactor2[i])
			|| isComplexNaN((*s).gridPolarizationFactor1[i])
			|| isComplexNaN((*s).gridPolarizationFactor2[i])){
				(*s).gridPropagationFactor1[i] = {};
				(*s).gridPropagationFactor2[i] = {};
				(*s).gridPolarizationFactor1[i] = {};
				(*s).gridPolarizationFactor2[i] = {};
			}
		}
	};

	//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
	//note that the sellmeier coefficients have extra values appended to the end
	//to give info about the current simulation
	class prepare3DGridsKernel {
	public:
		const deviceFP* sellmeierCoefficients;
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t localIndex) const {
			
			const int64_t col = localIndex / ((*s).Nfreq - 1); //spatial coordinate
			const int64_t j = 1 + localIndex % ((*s).Nfreq - 1); // frequency coordinate
			const int64_t i = j + col * (*s).Nfreq;
			const int64_t k = col % (*s).Nspace;
			const int64_t l = col / (*s).Nspace;

			//frequency being resolved by current thread
			const deviceFP f = j * (*s).fStep;

			//transverse wavevector being resolved
			const deviceFP dk1 = k * (*s).dk1 - (k >= ((*s).Nspace / 2)) 
				* ((*s).dk1 * (*s).Nspace); //frequency grid in x direction
			const deviceFP dk2 = l * (*s).dk2 - (l >= ((*s).Nspace2 / 2)) 
				* ((*s).dk2 * (*s).Nspace2); //frequency grid in y direction

			deviceComplex ne, no;
			deviceFP d = findBirefringentCrystalIndex(s, sellmeierCoefficients, localIndex, &ne, &no);
			if(s->axesNumber==2 && i<s->NgridC) s->gridBiaxialDelta[i] = d;
			if (minN(ne.real(), no.real()) < 0.9f 
				|| isComplexNaN(ne) 
				|| isComplexNaN(no)) {
				(*s).gridPropagationFactor1[i] = {};
				(*s).gridPropagationFactor2[i] = {};
				(*s).gridPolarizationFactor1[i] = {};
				(*s).gridPolarizationFactor2[i] = {};
				return;
			}

			deviceComplex k0 = twoPi<deviceFP>() * (*s).n0 * f / lightC<deviceFP>();
			deviceComplex ke = twoPi<deviceFP>() * ne * f / lightC<deviceFP>();
			deviceComplex ko = (deviceFP)twoPi<deviceFP>() * no * f / lightC<deviceFP>();

			//chi11 factor, also multiplied by -sqrt(-1)!
			const deviceComplex chi11 = ((*s).isUsingMillersRule) ?
				deviceComplex((*s).chiLinear1[j].imag(), -(*s).chiLinear1[j].real())
				: cOne<deviceComplex>();
			const deviceComplex chi12 = ((*s).isUsingMillersRule) ?
				deviceComplex((*s).chiLinear2[j].imag(), -(*s).chiLinear2[j].real())
				: cOne<deviceComplex>();

			deviceComplex kz1 = deviceLib::sqrt(ke * ke - dk1 * dk1 - dk2 * dk2);
			deviceComplex kz2 = deviceLib::sqrt(ko * ko - dk1 * dk1 - dk2 * dk2);
			if (kz1.real() > 0.0f && kz2.real() > 0.0f) {
				(*s).gridPropagationFactor1[i] = fourierPropagator(ke, dk1, dk2, k0.real(), 0.5f * (*s).h);
				if (isnan(((*s).gridPropagationFactor1[i].real()))) {
					(*s).gridPropagationFactor1[i] = {};
				}

				(*s).gridPropagationFactor2[i] = fourierPropagator(ko, dk1, dk2, k0.real(), 0.5f * (*s).h);
				if (isnan(((*s).gridPropagationFactor2[i].real()))) {
					(*s).gridPropagationFactor2[i] = {};
				}

				(*s).gridPolarizationFactor1[i] = 
					deviceLib::pow((*s).chiLinear1[j] + 1.0f, 0.25f) 
					* chi11 * (twoPi<deviceFP>() * twoPi<deviceFP>() * f * f) 
					/ (2.0f * lightC<deviceFP>() * lightC<deviceFP>() * kz1) * (*s).h;
				(*s).gridPolarizationFactor2[i] = 
					deviceLib::pow((deviceComplex)(*s).chiLinear2[j] + 1.0f, 0.25f) 
					* chi12 * (twoPi<deviceFP>() * twoPi<deviceFP>() * f * f) 
					/ (2.0f * lightC<deviceFP>() * lightC<deviceFP>() * kz2) * (*s).h;
			}
			else {
				(*s).gridPropagationFactor1[i] = {};
				(*s).gridPropagationFactor2[i] = {};
				(*s).gridPolarizationFactor1[i] = {};
				(*s).gridPolarizationFactor2[i] = {};
			}

			if (isComplexNaN((*s).gridPropagationFactor1[i])
				|| isComplexNaN((*s).gridPropagationFactor2[i])
				|| isComplexNaN((*s).gridPolarizationFactor1[i])
				|| isComplexNaN((*s).gridPolarizationFactor2[i])) {
				(*s).gridPropagationFactor1[i] = {};
				(*s).gridPropagationFactor2[i] = {};
				(*s).gridPolarizationFactor1[i] = {};
				(*s).gridPolarizationFactor2[i] = {};
			}
		}
	};

	//prepare the chi(1) arrays that will be needed in the simulation
	class getChiLinearKernel {
	public:
		deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP* sellmeierCoefficients;
		deviceFunction void operator()(const int64_t i) const {
			//frequency being resolved by current thread
			const deviceFP f = static_cast<deviceFP>(i) * (*s).fStep;
			deviceComplex ne, no;
			if (i < (*s).Ntime / 2) {
				sellmeierDevice(
					&ne, 
					&no, 
					sellmeierCoefficients, 
					f, 
					(*s).crystalTheta, 
					(*s).crystalPhi, 
					(*s).axesNumber, 
					(*s).sellmeierType,
					false); //note: false to applySqrt, so it returns n^2

				(*s).chiLinear1[i] = ne - 1.0f;
				(*s).chiLinear2[i] = no - 1.0f;
				if ((*s).chiLinear1[i].real() != 0.0f && (*s).chiLinear2[i].real() != 0.0f) {
					(*s).inverseChiLinear1[i] = 1.0f / (*s).chiLinear1[i].real();
					(*s).inverseChiLinear2[i] = 1.0f / (*s).chiLinear2[i].real();
				}
				else {
					(*s).inverseChiLinear1[i] = {};
					(*s).inverseChiLinear2[i] = {};
				}
				
				(*s).fieldFactor1[i] = deviceFPLib::pow((*s).chiLinear1[i].real() + 1.0f, -0.25f); 
				//account for the effective field strength in the medium (1/n)
				(*s).fieldFactor2[i] = deviceFPLib::pow((*s).chiLinear2[i].real() + 1.0f, -0.25f);
				if ((*s).isUsingMillersRule) {
					(*s).fieldFactor1[i] *= (*s).chiLinear1[i].real();
					(*s).fieldFactor2[i] *= (*s).chiLinear2[i].real();
				}
				(*s).fieldFactor1[i] *= (*s).fftNorm;
				(*s).fieldFactor2[i] *= (*s).fftNorm;
				if (isComplexNaN(ne) 
					|| isComplexNaN(no)
					|| isComplexNaN((*s).chiLinear1[i])
					|| isComplexNaN((*s).chiLinear2[i])
					|| isnan((*s).fieldFactor1[i])
					|| isnan((*s).fieldFactor2[i])
					|| isnan((*s).inverseChiLinear1[i])
					|| isnan((*s).inverseChiLinear2[i])
					|| ne.real() < 0.9f 
					|| no.real() < 0.9f 
					|| ne.imag() > 0.0f 
					|| no.imag() > 0.0f) {
					ne = cOne<deviceComplex>();
					no = ne;
					(*s).fieldFactor1[i] = {};
					(*s).fieldFactor2[i] = {};
					(*s).inverseChiLinear1[i] = {};
					(*s).inverseChiLinear2[i] = {};
					(*s).chiLinear1[i] = {};
					(*s).chiLinear2[i] = {};
				}
			}
			

			if (i == 81) {
				deviceComplex n0;
				sellmeierDevice(
					&n0, 
					&no, 
					sellmeierCoefficients, 
					deviceFPLib::abs((*s).f0), 
					(*s).crystalTheta, 
					(*s).crystalPhi, 
					(*s).axesNumber, 
					(*s).sellmeierType);
				(*s).n0 = no;
				(*s).chiLinear1[(*s).Ntime / 2] = cOne<deviceComplex>();
				(*s).chiLinear2[(*s).Ntime / 2] = cOne<deviceComplex>();
				(*s).fieldFactor1[(*s).Ntime / 2] = {};
				(*s).fieldFactor2[(*s).Ntime / 2] = {};
				(*s).inverseChiLinear1[(*s).Ntime / 2] = 1.0f / (*s).chiLinear1[(*s).Ntime / 2].real();
				(*s).inverseChiLinear2[(*s).Ntime / 2] = 1.0f / (*s).chiLinear2[(*s).Ntime / 2].real();
			}

			//apply Miller's rule to nonlinear coefficients
			if (!(*s).isUsingMillersRule || i > 80) {
				return;
			}
			const deviceFP* referenceFrequencies = &sellmeierCoefficients[72];
			deviceFP chi11[7];

			for (int im = (i > 17) * 3; im < 7; ++im) {
				if (referenceFrequencies[im] == 0.0f) {
					chi11[im] = 100000.0f;
				}
				else {
					sellmeierDevice(
						&ne, 
						&no, 
						sellmeierCoefficients, 
						referenceFrequencies[im], 
						(*s).crystalTheta, 
						(*s).crystalPhi, 
						(*s).axesNumber, 
						(*s).sellmeierType,
						false); //note: false to applySqrt, so it returns n^2
					chi11[im] = no.real() - 1.0f;
				}
			}

			//normalize chi2 tensor values
			if (i < 18) {
				(*s).chi2Tensor[i] /= chi11[0] * chi11[1] * chi11[2];
				if(isnan((*s).chi2Tensor[i])) (*s).chi2Tensor[i] = {};
			}

			//normalize chi3 tensor values
			(*s).chi3Tensor[i] /= chi11[3] * chi11[4] * chi11[5] * chi11[6];
			if(isnan((*s).chi3Tensor[i])) (*s).chi3Tensor[i] = {};
		}
	};

	//prepare the propagation constants under the assumption of cylindrical symmetry of the beam
	class prepareCylindricGridsKernel { 
	public:
		deviceFP* sellmeierCoefficients;
		deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(int64_t i) const {
			const int64_t j = i / ((*s).Nfreq - 1); //spatial coordinate
			const int64_t k = 1 + i % ((*s).Nfreq - 1); //temporal coordinate
			i = static_cast<deviceFP>(k) + static_cast<deviceFP>(j) * (*s).Nfreq;
			const deviceComplex ii = deviceComplex(0.0f, 1.0f);
			//frequency being resolved by current thread
			const deviceFP f = -(static_cast<deviceFP>(k) * (*s).fStep);
			//transverse wavevector being resolved
			const deviceFP dk = static_cast<deviceFP>(j) * (*s).dk1 - static_cast<deviceFP>(j >= ((*s).Nspace / 2)) 
				* ((*s).dk1 * (*s).Nspace); //frequency grid in transverse direction
			deviceComplex ne, no;
			deviceFP delta = sellmeierDevice(
				&ne, 
				&no, 
				sellmeierCoefficients, 
				(*s).fStep * k, 
				(*s).crystalTheta, 
				(*s).crystalPhi, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			if((*s).axesNumber==2) (*s).gridBiaxialDelta[i] = delta;
			//if the refractive index was returned weird, then the index isn't valid, 
			//so set the propagator to zero for that frequency
			if (minN(ne.real(), no.real()) < 0.95f 
				|| ne.real() > 6.0f 
				|| no.real() > 6.0f 
				|| isComplexNaN(ne) 
				|| isComplexNaN(no)) {
				(*s).gridPropagationFactor1[i] = {};
				(*s).gridPropagationFactor2[i] = {};
				(*s).gridPolarizationFactor1[i] = {};
				(*s).gridPolarizationFactor2[i] = {};
				(*s).gridPropagationFactor1Rho1[i] = {};
				(*s).gridPropagationFactor1Rho2[i] = {};
				return;
			}

			const deviceComplex k0 = deviceComplex(
				twoPi<deviceFP>() * (*s).n0.real() * f / lightC<deviceFP>(), 
				0.0f);
			const deviceComplex ke = twoPi<deviceFP>() * ne * f / lightC<deviceFP>();
			const deviceComplex ko = twoPi<deviceFP>() * no * f / lightC<deviceFP>();

			const deviceComplex chi11 = ((*s).isUsingMillersRule) ? 
				(*s).chiLinear1[k] : cOne<deviceComplex>();
			const deviceComplex chi12 = ((*s).isUsingMillersRule) ? 
				(*s).chiLinear2[k] : cOne<deviceComplex>();
			
			if ((dk * dk 
				< minN(ke.real() * ke.real() + ke.imag() * ke.imag(), 
					ko.real() * ko.real() + ko.imag() * ko.imag())) 
				&& (*s).fieldFactor1[k] > 0.0f 
				&& (*s).fieldFactor2[k] > 0.0f) {
				(*s).gridPropagationFactor1[i] = 0.5f * ii * (ke - k0 - dk * dk / (2.0f * ke.real())) * (*s).h;
				(*s).gridPropagationFactor1[i] = isComplexNaN((*s).gridPropagationFactor1[i]) ? 
					deviceComplex{} : deviceLib::exp((*s).gridPropagationFactor1[i]);
				(*s).gridPropagationFactor1Rho1[i] = 
					(*s).fftNorm * ii * (*s).h / ((*s).fieldFactor1[k] * 2.0f * ke);
				if (isnan(
					(deviceLib::abs((*s).gridPropagationFactor1Rho1[i] 
						+ (*s).gridPropagationFactor1[i])))) {
					(*s).gridPropagationFactor1[i] = {};
					(*s).gridPropagationFactor1Rho1[i] = {};
				}

				(*s).gridPropagationFactor2[i] = 0.5f * ii * (ko - k0 - dk * dk 
						/ (2.0f * ko.real())) * (*s).h;
				(*s).gridPropagationFactor2[i] = isComplexNaN((*s).gridPropagationFactor2[i]) ?	
					deviceComplex{} : deviceLib::exp((*s).gridPropagationFactor2[i]);
				(*s).gridPropagationFactor1Rho2[i] = 
					(*s).fftNorm * ii * (*s).h 
					/ ((*s).fieldFactor2[k] * 2.0f * ko);
				if (isnan(
					(deviceLib::abs((*s).gridPropagationFactor1Rho2[i] 
						+ (*s).gridPropagationFactor2[i])))) {
					(*s).gridPropagationFactor2[i] = {};
					(*s).gridPropagationFactor1Rho2[i] = {};
				}

				//factor of 0.5 comes from deviceFPd grid size in cylindrical symmetry 
				//mode after expanding the beam
				(*s).gridPolarizationFactor1[i] = 
					0.5f * deviceLib::pow((deviceComplex)(*s).chiLinear1[k] + 1.0f, 0.25f) 
					* chi11 * ii * (twoPi<deviceFP>() * f) 
					/ (2.0f * ne.real() * lightC<deviceFP>()) * (*s).h;
				(*s).gridPolarizationFactor2[i] = 
					0.5f * deviceLib::pow((deviceComplex)(*s).chiLinear2[k] + 1.0f, 0.25f) 
					* chi12 * ii * (twoPi<deviceFP>() * f) 
					/ (2.0f * no.real() * lightC<deviceFP>()) * (*s).h;
				if (isComplexNaN((*s).gridPolarizationFactor1[i]) 
					|| isComplexNaN((*s).gridPolarizationFactor2[i])) {
					(*s).gridPolarizationFactor1[i] = {};
					(*s).gridPolarizationFactor2[i] = {};
				}
			}
			else {
				(*s).gridPropagationFactor1[i] = {};
				(*s).gridPropagationFactor2[i] = {};
				(*s).gridPolarizationFactor1[i] = {};
				(*s).gridPolarizationFactor2[i] = {};
				(*s).gridPropagationFactor1Rho1[i] = {};
				(*s).gridPropagationFactor1Rho2[i] = {};
			}
		}
	};

	class materialPhaseKernel { 
	public:
		const deviceFP df;
		const int64_t Ntime;
		const deviceFP* a;
		const deviceFP f0;
		const deviceFP thickness;
		const int sellmeierEquation;
		deviceFP* phase1;
		deviceFunction void operator()(const int64_t i) const {
			//frequency being resolved by current thread
			const deviceFP f = i * df;
			//give phase shift relative to group velocity (approximated 
			// with low-order finite difference) so the pulse doesn't move
			deviceComplex ne, no, no0, n0p, n0m;
			sellmeierDevice<deviceFP, deviceComplex>(
				&ne, &no, a, f, {}, {}, 0, sellmeierEquation);
			sellmeierDevice<deviceFP, deviceComplex>(
				&ne, &no0, a, f0, {}, {}, 0, sellmeierEquation);
			sellmeierDevice<deviceFP, deviceComplex>(
				&ne, &n0p, a, f0 + 1.0e11f, {}, {}, 0, sellmeierEquation);
			sellmeierDevice<deviceFP, deviceComplex>(
				&ne, &n0m, a, f0 - 1.0e11f, {}, {}, 0, sellmeierEquation);
			no0 = no0 + f0 * (n0p - n0m) / 2.0e11f;
			phase1[i] = thickness * twoPi<deviceFP>() * f 
				* (no.real() - no0.real()) / lightC<deviceFP>();
		}
	};

	//calculate the nonlinear polarization, after FFT to get the field
	//in the time domain
	class nonlinearPolarizationKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t i) const {
			const deviceFP Ex = (*s).gridETime1[i];
			const deviceFP Ey = (*s).gridETime2[i];

			if ((*s).nonlinearSwitches[0] == 1 || (*s).nonlinearSwitches[1] == 1) {

				deviceFP P[3]{};

				// rotate field into crystal frame
				const deviceFP E3[3] = { (*s).rotationForward[0] * Ex + (*s).rotationForward[1] * Ey,
								(*s).rotationForward[3] * Ex + (*s).rotationForward[4] * Ey,
								(*s).rotationForward[6] * Ex + (*s).rotationForward[7] * Ey };

				// second order nonlinearity, element-by-element in the reduced tensor
				//note that for historical reasons, chi2 is column-order and chi3 is row-order...
				if ((*s).nonlinearSwitches[0] == 1) {
					for (auto a = 0; a < 3; ++a) {
						P[a] += (*s).chi2Tensor[0 + a] * E3[0] * E3[0];
						P[a] += (*s).chi2Tensor[3 + a] * E3[1] * E3[1];
						P[a] += (*s).chi2Tensor[6 + a] * E3[2] * E3[2];
						P[a] += (*s).chi2Tensor[9 + a] * E3[1] * E3[2];
						P[a] += (*s).chi2Tensor[12 + a] * E3[0] * E3[2];
						P[a] += (*s).chi2Tensor[15 + a] * E3[0] * E3[1];
					}
				}

				// resolve the full chi3 matrix when (*s).nonlinearSwitches[1]==1
				if ((*s).nonlinearSwitches[1] == 1) {
					// loop over tensor element X_abcd
					// i hope the compiler unrolls this, but no way am I writing that out by hand
					for (auto a = 0; a < 3; ++a) {
						for (auto b = 0; b < 3; ++b) {
							for (auto c = 0; c < 3; ++c) {
								for (auto d = 0; d < 3; ++d) {
									P[d] += (*s).chi3Tensor[a + 3 * b + 9 * c + 27 * d] 
										* E3[a] * E3[b] * E3[c];
								}
							}
						}
					}
				}

				(*s).gridPolarizationTime1[i] = 
					(*s).rotationBackward[0] * P[0] 
					+ (*s).rotationBackward[1] * P[1] 
					+ (*s).rotationBackward[2] * P[2];
				(*s).gridPolarizationTime2[i] = 
					(*s).rotationBackward[3] * P[0] 
					+ (*s).rotationBackward[4] * P[1] 
					+ (*s).rotationBackward[5] * P[2];

				//using only one value of chi3, under assumption of centrosymmetry when nonlinearSwitches[1]==2
				if ((*s).nonlinearSwitches[1] == 2) {
					const deviceFP Esquared = (*s).chi3Tensor[0] * (Ex * Ex + Ey * Ey);
					(*s).gridPolarizationTime1[i] += Ex * Esquared;
					(*s).gridPolarizationTime2[i] += Ey * Esquared;
				}
			}
			else {
				const deviceFP Esquared = (*s).chi3Tensor[0] * (Ex * Ex + Ey * Ey);
				(*s).gridPolarizationTime1[i] = Ex * Esquared;
				(*s).gridPolarizationTime2[i] = Ey * Esquared;
			}
			if ((*s).isCylindric)expandCylindricalBeamDevice(
				s, 
				i, 
				(*s).gridRadialLaplacian1 + (*s).Ngrid * 4 * (*s).hasPlasma, 
				(*s).gridPolarizationTime1, 
				(*s).gridPolarizationTime2);
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
	class plasmaCurrentKernel_twoStage_A {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t i) const {
			const int pMax = static_cast<int>((*s).nonlinearSwitches[3]);

			//save values in workspaces, casting to deviceFP
			deviceFP* dN = (deviceFP*)(*s).workspace1;
			const deviceFP Esquared = 
				(*s).plasmaParameters[0] 
				* ((*s).gridETime1[i] * (*s).gridETime1[i] 
					+ (*s).gridETime2[i] * (*s).gridETime2[i]);
			deviceFP EtoThe2N = Esquared;
			for (int p = 0; p < pMax; ++p) {
				EtoThe2N *= Esquared;
			}
			(*s).gridPolarizationTime1[i] = EtoThe2N * (*s).gridETime1[i];
			(*s).gridPolarizationTime2[i] = EtoThe2N * (*s).gridETime2[i];
			dN[i] = (*s).plasmaParameters[2] * (
				(*s).gridPolarizationTime1[i] * (*s).gridETime1[i] 
				+ (*s).gridPolarizationTime2[i] * (*s).gridETime2[i]);
		}
	};

	class plasmaCurrentKernel_twoStage_B {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t j = i * (*s).Ntime;
			deviceFP N{};
			deviceFP integralx{};
			const deviceFP* expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
			const deviceFP* dN = (j % (*s).Ngrid) + (deviceFP*)(*s).workspace1;
			const deviceFP* E = &(*s).gridETime1[j];
			deviceFP* P = &(*s).gridPolarizationTime1[j];
			for (int64_t k{}; k < (*s).Ntime; ++k) {
				N += dN[k];
				integralx += N * (*s).expGammaT[k] * E[k];
				P[k] += expMinusGammaT[k] * integralx;
			}
		}
	};

	class plasmaCurrentKernel_twoStage_B_simultaneous {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t j = (i) * (*s).Ntime;
			deviceFP N{};
			deviceFP integralx{};
			deviceFP integraly{};
			const deviceFP* expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
			const deviceFP* dN = j + (deviceFP*)(*s).workspace1;
			const deviceFP* E = &(*s).gridETime1[j];
			const deviceFP* Ey = E + (*s).Ngrid;
			deviceFP* P = &(*s).gridPolarizationTime1[j];
			deviceFP* P2 = P + (*s).Ngrid;
			for (int64_t k{}; k < (*s).Ntime; ++k) {
				N += dN[k];
				const deviceFP plasmaFactor = N * (*s).expGammaT[k];
				integralx += plasmaFactor * E[k];
				integraly += plasmaFactor * Ey[k];
				P[k] += expMinusGammaT[k] * integralx;
				P2[k] += expMinusGammaT[k] * integraly;
			}
		}
	};

	class plasmaCurrentKernel_SaveOutput {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t j = (i) * (*s).Ntime;
			deviceFP N{};
			deviceFP integralx{};
			deviceFP integraly{};
			const deviceFP* expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
			const deviceFP* dN = j + reinterpret_cast<deviceFP*>((*s).workspace1);
			const deviceFP* E = &(*s).gridETime1[j];
			const deviceFP* Ey = E + (*s).Ngrid;
			deviceFP* P = &(*s).gridPolarizationTime1[j];
			deviceFP* P2 = P + (*s).Ngrid;
			deviceFP* Ndest = reinterpret_cast<deviceFP*>((*s).gridEFrequency1);
			Ndest += j;
			for (auto k = 0; k < (*s).Ntime; ++k) {
				N += dN[k];
				integralx += N * (*s).expGammaT[k] * E[k];
				integraly += N * (*s).expGammaT[k] * Ey[k];
				P[k] += expMinusGammaT[k] * integralx;
				P2[k] += expMinusGammaT[k] * integraly;
				Ndest[k] = N/(*s).dt;
			}
		}
	};

	class updateKwithPolarizationKernelCylindric {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
			const int64_t radIndex = i / ((*sP).Nfreq - 1);
			const int64_t gridIndex = freqIndex + radIndex * ((*sP).Nfreq);
			const int64_t polIndex = freqIndex + 
				(radIndex + ((radIndex > ((*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;
			(*sP).k1[gridIndex] += 
				(*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[polIndex];
			(*sP).k2[gridIndex] += 
				(*sP).gridPolarizationFactor2[gridIndex] * (*sP).workspace2P[polIndex];
		}
	};

	class updateKwithPlasmaKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			int64_t h = 1 + i % ((*sP).Nfreq - 1); //temporal coordinate
			const int64_t j = i / ((*sP).Nfreq - 1); //spatial coordinate
			const deviceFP jfac = -1.0f / (h * (*sP).fStep);
			h += j * (*sP).Nfreq;
			(*sP).k1[h] += deviceComplex{-jfac * (*sP).gridPolarizationFactor1[h].imag(), jfac * (*sP).gridPolarizationFactor1[h].real()} 
				* (*sP).workspace1[h] * (*sP).inverseChiLinear1[h % ((*sP).Nfreq)];
			(*sP).k2[h] += deviceComplex{-jfac * (*sP).gridPolarizationFactor2[h].imag(), jfac * (*sP).gridPolarizationFactor2[h].real()} 
				* (*sP).workspace2P[h] * (*sP).inverseChiLinear2[h % ((*sP).Nfreq)];
		}
	};

	class updateKwithPlasmaKernelCylindric {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t spaceIndex = i / ((*sP).Nfreq - 1);
			int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
			const deviceFP jfac = -1.0f/(freqIndex * (*sP).fStep);
			const int64_t gridIndex = freqIndex + spaceIndex * ((*sP).Nfreq);
			int64_t fftIndex = freqIndex +
				(spaceIndex + ((spaceIndex > ((*sP).Nspace / 2))) * (*sP).Nspace) * (*sP).Nfreq;
			freqIndex = gridIndex % (*sP).Nfreq;
			(*sP).k1[gridIndex] += 
				deviceComplex{
					-jfac * (*sP).gridPolarizationFactor1[gridIndex].imag(), 
					jfac * (*sP).gridPolarizationFactor1[gridIndex].real()}
				* (*sP).workspace1[fftIndex] * (*sP).inverseChiLinear1[freqIndex];
			(*sP).k2[gridIndex] += 
				deviceComplex{
					-jfac * (*sP).gridPolarizationFactor2[gridIndex].imag(), 
					jfac * (*sP).gridPolarizationFactor2[gridIndex].real()}
				* (*sP).workspace2P[fftIndex] * (*sP).inverseChiLinear2[freqIndex];

			fftIndex += 4 * (*sP).NgridC;
			(*sP).k1[gridIndex] += 
				(*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[fftIndex];
			(*sP).k2[gridIndex] += 
				(*sP).gridPolarizationFactor2[gridIndex] * (*sP).workspace2P[fftIndex];
		}
	};

	class biaxialRotationKernel {
		public:
			const deviceParameterSet<deviceFP, deviceComplex>* sP;
			deviceComplex* x;
			bool backwards;
			deviceFunction void operator()(const int64_t gridIndex) const {
				deviceComplex* y = x + sP->NgridC;
				const int64_t freqIndex = 1 + gridIndex % ((*sP).Nfreq - 1); //frequency coordinate
				const int64_t i = freqIndex + (gridIndex / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);

				deviceFP cosDelta = deviceFPLib::cos(sP->gridBiaxialDelta[i]);
				deviceFP sinDelta = deviceFPLib::sin(sP->gridBiaxialDelta[i]);
				if(backwards) sinDelta = -sinDelta;

				deviceComplex newX = cosDelta * x[i] - sinDelta * y[i];
				y[i] = sinDelta * x[i] + cosDelta * y[i];
				x[i] = newX;

			}
	};

	//Slightly different kernels for the four stages of RK4. 
	//They used to be one big kernel with a switch case
	//but this has slightly better utilization.
	class rkKernel0 {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(int64_t gridIndex) const {
			const int64_t freqIndex = 1 + gridIndex % ((*sP).Nfreq - 1); //frequency coordinate
			gridIndex = freqIndex + (gridIndex / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
			(*sP).k1[gridIndex] += (*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[gridIndex];
			(*sP).gridEFrequency1Next1[gridIndex] = (*sP).gridPropagationFactor1[gridIndex] 
				* (*sP).gridPropagationFactor1[gridIndex] 
				* (sixth<deviceFP>() * (*sP).k1[gridIndex] + (*sP).gridEFrequency1[gridIndex]);
			const deviceFP ff = (gridIndex > (*sP).NgridC) ? 
				(*sP).fieldFactor2[freqIndex] 
				: (*sP).fieldFactor1[freqIndex];
			(*sP).workspace1[gridIndex] = ff * (*sP).gridPropagationFactor1[gridIndex] 
				* ((*sP).gridEFrequency1[gridIndex] + 0.5f * (*sP).k1[gridIndex]);
			(*sP).k1[gridIndex] = {};
			[[unlikely]] if (freqIndex == 1) (*sP).workspace1[gridIndex - 1] = {};
		}
	};

	class rkKernel1 {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1); //frequency coordinate
			const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
			(*sP).k1[gridIndex] += 
				(*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[gridIndex];
			(*sP).gridEFrequency1Next1[gridIndex] = (*sP).gridEFrequency1Next1[gridIndex] 
				+ (*sP).gridPropagationFactor1[gridIndex] 
				* (deviceFP)third<deviceFP>() * (*sP).k1[gridIndex];
			const deviceFP ff = (gridIndex > (*sP).NgridC) ?
				(*sP).fieldFactor2[freqIndex] 
				: (*sP).fieldFactor1[freqIndex];
			(*sP).workspace1[gridIndex] = ff * ((*sP).gridPropagationFactor1[gridIndex] 
				* (*sP).gridEFrequency1[gridIndex] + 0.5f * (*sP).k1[gridIndex]);
			(*sP).k1[gridIndex] = {};
			[[unlikely]] if (freqIndex == 1) (*sP).workspace1[gridIndex - 1] = {};
		}
	};

	class rkKernel2 {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1); //frequency coordinate
			const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
			(*sP).k1[gridIndex] += 
				(*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[gridIndex];
			(*sP).gridEFrequency1Next1[gridIndex] = (*sP).gridEFrequency1Next1[gridIndex] 
				+ (*sP).gridPropagationFactor1[gridIndex] 
				* (deviceFP)third<deviceFP>() * (*sP).k1[gridIndex];
			const deviceFP ff = (gridIndex > (*sP).NgridC) ? 
				(*sP).fieldFactor2[freqIndex] 
				: (*sP).fieldFactor1[freqIndex];
			(*sP).workspace1[gridIndex] = ff * ((*sP).gridPropagationFactor1[gridIndex] 
				* ((*sP).gridPropagationFactor1[gridIndex] * (*sP).gridEFrequency1[gridIndex] + (*sP).k1[gridIndex]));
			(*sP).k1[gridIndex] = {};
			[[unlikely]] if (freqIndex == 1) (*sP).workspace1[gridIndex - 1] = {};
		}
	};

	class rkKernel3 {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1); //frequency coordinate
			const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
			(*sP).k1[gridIndex] += 
				(*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[gridIndex];
			(*sP).gridEFrequency1[gridIndex] = (*sP).gridEFrequency1Next1[gridIndex] 
				+ sixth<deviceFP>() * (*sP).k1[gridIndex];
			const deviceFP ff = (gridIndex > (*sP).NgridC) ? 
				(*sP).fieldFactor2[freqIndex] 
				: (*sP).fieldFactor1[freqIndex];
			(*sP).workspace1[gridIndex] = ff * (*sP).gridEFrequency1[gridIndex];
			(*sP).k1[gridIndex] = {};
			[[unlikely]] if (freqIndex == 1) (*sP).workspace1[gridIndex - 1] = {};
		}
	};

	//Kernels for symmetry around z axis use a different form, adding the radial Laplacian
	//instead of the nonlinear polarization
	class rkKernel0Cylindric {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
			const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * (*sP).Nfreq;
			(*sP).k1[gridIndex] += 
				(*sP).gridPropagationFactor1Rho1[gridIndex] * (*sP).workspace1[gridIndex];
			(*sP).gridEFrequency1Next1[gridIndex] = (*sP).gridPropagationFactor1[gridIndex] 
				* (*sP).gridPropagationFactor1[gridIndex] 
				* (sixth<deviceFP>() * (*sP).k1[gridIndex] + (*sP).gridEFrequency1[gridIndex]);
			const deviceFP ff = (gridIndex > (*sP).NgridC) ? 
				(*sP).fieldFactor2[freqIndex] : (*sP).fieldFactor1[freqIndex];
			(*sP).workspace1[gridIndex] = ff * ((*sP).gridPropagationFactor1[gridIndex] 
				* ((*sP).gridEFrequency1[gridIndex] + 0.5f * (*sP).k1[gridIndex]));
			(*sP).k1[gridIndex] = {};
			[[unlikely]] if (freqIndex == 1) (*sP).workspace1[gridIndex - 1] = {};
		}
	};

	class rkKernel1Cylindric {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
			const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * (*sP).Nfreq;
			(*sP).k1[gridIndex] += 
				(*sP).gridPropagationFactor1Rho1[gridIndex] * (*sP).workspace1[gridIndex];
			(*sP).gridEFrequency1Next1[gridIndex] = (*sP).gridEFrequency1Next1[gridIndex] 
				+ (*sP).gridPropagationFactor1[gridIndex] 
				* third<deviceFP>() * (*sP).k1[gridIndex];
			const deviceFP ff = (gridIndex > (*sP).NgridC) ? 
				(*sP).fieldFactor2[freqIndex] : (*sP).fieldFactor1[freqIndex];
			(*sP).workspace1[gridIndex] = ff * ((*sP).gridPropagationFactor1[gridIndex] 
				* (*sP).gridEFrequency1[gridIndex] + 0.5f * (*sP).k1[gridIndex]);
			(*sP).k1[gridIndex] = {};
			[[unlikely]] if (freqIndex == 1) (*sP).workspace1[gridIndex - 1] = {};
		}
	};

	class rkKernel2Cylindric {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
			const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * (*sP).Nfreq;
			(*sP).k1[gridIndex] += 
				(*sP).gridPropagationFactor1Rho1[gridIndex] * (*sP).workspace1[gridIndex];
			(*sP).gridEFrequency1Next1[gridIndex] = (*sP).gridEFrequency1Next1[gridIndex] 
				+ (*sP).gridPropagationFactor1[gridIndex] 
				* third<deviceFP>() * (*sP).k1[gridIndex];
			const deviceFP ff = (gridIndex > (*sP).NgridC) ? 
				(*sP).fieldFactor2[freqIndex] : (*sP).fieldFactor1[freqIndex];
			(*sP).workspace1[gridIndex] = ff * ((*sP).gridPropagationFactor1[gridIndex] 
				* ((*sP).gridPropagationFactor1[gridIndex] * (*sP).gridEFrequency1[gridIndex] 
					+ (*sP).k1[gridIndex]));
			(*sP).k1[gridIndex] = {};
			[[unlikely]] if (freqIndex == 1) (*sP).workspace1[gridIndex - 1] = {};
		}
	};

	class rkKernel3Cylindric {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* sP;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
			const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * (*sP).Nfreq;
			(*sP).k1[gridIndex] += 
				(*sP).gridPropagationFactor1Rho1[gridIndex] * (*sP).workspace1[gridIndex];
			(*sP).gridEFrequency1[gridIndex] = (*sP).gridEFrequency1Next1[gridIndex] 
				+ sixth<deviceFP>() * (*sP).k1[gridIndex];
			const deviceFP ff = (gridIndex > (*sP).NgridC) ? 
				(*sP).fieldFactor2[freqIndex] : (*sP).fieldFactor1[freqIndex];
			(*sP).workspace1[gridIndex] = ff * (*sP).gridEFrequency1[gridIndex];
			(*sP).k1[gridIndex] = {};
			[[unlikely]] if (freqIndex == 1) (*sP).workspace1[gridIndex - 1] = {};
		}
	};

	class maxwellRKkernel0 {
	public:
		const maxwell3D* s;
		const int64_t t;
		deviceFunction void operator()(const int64_t i) const {
			maxwellKPoint<deviceFP> k = maxwellDerivativeTerms(s, i, s->Egrid, s->Hgrid);
			maxwellCurrentTerms(s, i, t, false, s->Egrid, s->materialGrid, k, 0);
			s->EgridEstimate[i] = s->Egrid[i] + k.kE * (s->tStep * 0.5f);
			s->EgridNext[i] = s->Egrid[i] + k.kE * (sixth<deviceFP>() * s->tStep);
			s->HgridEstimate[i] = s->Hgrid[i] + k.kH * (s->tStep * 0.5f);
			s->HgridNext[i] = s->Hgrid[i] + k.kH * (sixth<deviceFP>() * s->tStep);
		}
	};
	class maxwellRKkernel1 {
	public:
		const maxwell3D* s;
		const int64_t t;
		deviceFunction void operator()(const int64_t i) const {
			maxwellKPoint<deviceFP> k = 
				maxwellDerivativeTerms(s, i, s->EgridEstimate, s->HgridEstimate);
			maxwellCurrentTerms(s, i, t, true, s->EgridEstimate, s->materialGridEstimate, k, 1);
			s->EgridEstimate2[i] = s->Egrid[i] + k.kE * (s->tStep * 0.5f);
			s->EgridNext[i] += k.kE * (third<deviceFP>() * s->tStep);
			s->HgridEstimate2[i] = s->Hgrid[i] + k.kH * (s->tStep * 0.5f);
			s->HgridNext[i] += k.kH * (third<deviceFP>() * s->tStep);
		}
	};
	class maxwellRKkernel2 {
	public:
		const maxwell3D* s;
		const int64_t t;
		deviceFunction void operator()(const int64_t i) const {
			maxwellKPoint<deviceFP> k = 
				maxwellDerivativeTerms(s, i, s->EgridEstimate2, s->HgridEstimate2);
			maxwellCurrentTerms(s, i, t, true, s->EgridEstimate2, s->materialGridEstimate2, k, 2);
			s->EgridEstimate[i] = s->Egrid[i] + k.kE * s->tStep;
			s->EgridNext[i] += k.kE * (third<deviceFP>() * s->tStep);
			s->HgridEstimate[i] = s->Hgrid[i] + k.kH * s->tStep;
			s->HgridNext[i] += k.kH * (third<deviceFP>() * s->tStep);
		}
	};
	class maxwellRKkernel3 {
	public:
		const maxwell3D* s;
		const int64_t t;
		deviceFunction void operator()(const int64_t i) const {
			maxwellKPoint<deviceFP> k = 
				maxwellDerivativeTerms(s, i, s->EgridEstimate, s->HgridEstimate);
			maxwellCurrentTerms(s, i, t + 1, false, s->EgridEstimate, s->materialGridEstimate, k, 3);
			s->Egrid[i] = s->EgridNext[i] + k.kE * (sixth<deviceFP>() * s->tStep);
			s->Hgrid[i] = s->HgridNext[i] + k.kH * (sixth<deviceFP>() * s->tStep);
		}
	};

	//store the field in the observation plane in the in/out Ex and Ey arrays
	class maxwellSampleGrid {
	public:
		maxwell3D* s;
		int64_t time;
		deviceFunction void operator()(int64_t i) const {
			int64_t gridIndex = i * s->Nz + s->observationPoint;
			s->inOutEx[i * s->NtIO + time] = s->Egrid[gridIndex].x;
			s->inOutEy[i * s->NtIO + time] = s->Egrid[gridIndex].y;
		}
	};

	class beamNormalizeKernel {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		const deviceFP* rawSum;
		deviceFP* field;
		const deviceFP pulseEnergy;
		deviceFunction void operator()(const int64_t i) const {
			field[i] *= deviceFPLib::sqrt(pulseEnergy / ((deviceFP)(*s).Ntime * (*rawSum)));
		}
	};

	class addDoubleArraysKernel {
	public:
		deviceFP* A;
		deviceFP* B;
		deviceFunction void operator()(const int64_t i) const {
			A[i] += B[i];
		}
	};

	//crease a pulse on the grid for the 2D modes. 
	// Note that normalization of the 2D mode assumes radial symmetry (i.e. that it's
	//a gaussian beam, not an infinite plane wave, which would have zero amplitude for finite energy).
	class beamGenerationKernel2D {
	public:
		deviceComplex* field;
		const pulse<deviceFP>* p;
		deviceFP* pulseSum;
		deviceParameterSet<deviceFP, deviceComplex>* s;
		const bool hasLoadedField;
		const deviceComplex* loadedField;
		const deviceFP* materialPhase;
		const deviceFP* sellmeierCoefficients;
		deviceFunction void operator()(int64_t i) const {
			const int64_t h = 1 + i % ((*s).Nfreq - 1);
			const int64_t j = i / ((*s).Nfreq - 1);
			i = h + j * ((*s).Nfreq);
			const deviceFP f = h * (*s).fStep;
			const deviceFP w = twoPi<deviceFP>() * (f - (*p).frequency);

			//supergaussian pulse spectrum, if no input pulse specified
			deviceComplex specfac = deviceComplex(
				-deviceFPLib::pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder), 
				0.0f);

			if(isComplexNaN(materialPhase[h])){
				field[i] = deviceComplex{};
				field[i + (*s).NgridC] = deviceComplex{};
				return;
			}
			deviceComplex specphase = deviceComplex(0.0f,
				-((*p).cep
					+ twoPi<deviceFP>() * f * ((*p).delay - 0.5f * (*s).dt * (*s).Ntime)
					+ 0.5f * (*p).gdd * w * w
					+ sixth<deviceFP>() * (*p).tod * w * w * w
					+ materialPhase[h]));
			specfac = deviceLib::exp(specfac + specphase);

			if (hasLoadedField) {
				specfac = loadedField[h] * deviceLib::exp(specphase);
			}
			deviceComplex ne, no;
			sellmeierDevice(
				&ne, 
				&no, 
				sellmeierCoefficients, 
				f, (*s).crystalTheta, 
				(*s).crystalPhi, 
				(*s).axesNumber, 
				(*s).sellmeierType);

			if(isComplexNaN(ne) || isComplexNaN(no)){
				field[i] = deviceComplex{};
				field[i + (*s).NgridC] = deviceComplex{};
				return;
			}

			const deviceFP ko = twoPi<deviceFP>() * no.real() * f / lightC<deviceFP>();
			const deviceFP zR = vPi<deviceFP>() 
				* (*p).beamwaist * (*p).beamwaist 
				* no.real() * f / lightC<deviceFP>();

			const deviceFP rB = ((*p).x0 
				- (*s).dx * (-0.5f * (*s).Nspace + (deviceFP)j) - 0.25f * (*s).dx);
			const deviceFP r = rB * deviceFPLib::cos((*p).beamAngle) 
				- (*p).z0 * deviceFPLib::sin((*p).beamAngle);
			const deviceFP z = rB * deviceFPLib::sin((*p).beamAngle) 
				+ (*p).z0 * deviceFPLib::cos((*p).beamAngle);

			const deviceFP wz = (*p).beamwaist * deviceFPLib::sqrt(1.0f + (z * z / (zR * zR)));
			const deviceFP Rz = (z != 0.0f) ? z * (1.0f + (zR * zR / (z * z))) : 1.0e15f;
			const deviceFP phi = deviceFPLib::atan(z / zR);
			deviceComplex Eb = deviceComplex(0.0f, 1.0f) 
					* (ko * (z - (*p).z0) + ko * r * r / (2.0f * Rz) - phi) - r * r / (wz * wz);
			Eb = isComplexNaN(Eb) ? deviceComplex{} : ((*p).beamwaist / wz) * deviceLib::exp(Eb);
			Eb = Eb * specfac;
			
			field[i] = deviceComplex(
				deviceFPLib::cos((*p).polarizationAngle), 
				-(*p).circularity * deviceFPLib::sin((*p).polarizationAngle)) * Eb;
			field[i + (*s).NgridC] = deviceComplex(
				deviceFPLib::sin((*p).polarizationAngle), 
				(*p).circularity * deviceFPLib::cos((*p).polarizationAngle)) * Eb;
			if(isComplexNaN(field[i])) field[i] = deviceComplex{};
			if(isComplexNaN(field[i + (*s).NgridC])) field[i + (*s).NgridC] = deviceComplex{};
			deviceFP pointEnergy = deviceFPLib::abs(r) 
				* (modulusSquared(field[i]) + modulusSquared(field[i + (*s).NgridC]));
			pointEnergy *= 2.0f * vPi<deviceFP>() * lightC<deviceFP>() * eps0<deviceFP>() 
				* (*s).dx * (*s).dt;
			//two factors of two cancel here - there should be one for the missing 
			//frequency plane, but the sum is over x instead of r
			//accordingly we already multiplied by two
			atomicAdd(pulseSum, pointEnergy);
		}
	};

	//Generate a beam in full 3D mode
	class beamGenerationKernel3D {
	public:
		deviceComplex* field;
		const pulse<deviceFP>* p;
		deviceFP* pulseSum;
		deviceParameterSet<deviceFP, deviceComplex>* s;
		const bool hasLoadedField;
		const deviceComplex* loadedField;
		const deviceFP* materialPhase;
		const deviceFP* sellmeierCoefficients;
		deviceFunction void operator()(int64_t i) const {
			const int64_t h = 1 + i % ((*s).Nfreq - 1);
			const int64_t col = i / ((*s).Nfreq - 1);
			i = h + col * ((*s).Nfreq);
			const int64_t j = col % (*s).Nspace;
			const int64_t k = col / (*s).Nspace;
			const deviceFP f = h * (*s).fStep;
			const deviceFP w = twoPi<deviceFP>() * (f - (*p).frequency);
			if(isComplexNaN(materialPhase[h])){
				field[i] = deviceComplex{};
				field[i + (*s).NgridC] = deviceComplex{};
				return;
			}
			//supergaussian pulse spectrum, if no input pulse specified
			deviceComplex specfac = deviceComplex(
				-deviceFPLib::pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder), 
				0.0f);

			const deviceComplex specphase = deviceComplex(0.0f,
				-((*p).cep
					+ twoPi<deviceFP>() * f * ((*p).delay - 0.5f * (*s).dt * (*s).Ntime)
					+ 0.5f * (*p).gdd * w * w
					+ (*p).tod * w * w * w / 6.0f
					+ materialPhase[h]));
			specfac = isComplexNaN(specphase) ? deviceComplex{} : deviceLib::exp(specfac + specphase);

			if (hasLoadedField) {
				specfac = loadedField[h] * deviceLib::exp(specphase);
			}
			deviceComplex ne, no;
			sellmeierDevice(
				&ne, 
				&no, 
				sellmeierCoefficients, 
				f, 
				(*s).crystalTheta, 
				(*s).crystalPhi, 
				(*s).axesNumber, 
				(*s).sellmeierType);
			if(isComplexNaN(ne) || isComplexNaN(no)){
				field[i] = deviceComplex{};
				field[i + (*s).NgridC] = deviceComplex{};
				return;
			}
			const deviceFP ko = twoPi<deviceFP>() * no.real() * f 
				/ lightC<deviceFP>();
			const deviceFP zR = vPi<deviceFP>() * (*p).beamwaist * (*p).beamwaist 
				* no.real() * f / lightC<deviceFP>();

			const deviceFP xo = ((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f)) - (*p).x0;
			const deviceFP yo = ((*s).dx * ((deviceFP)k - (*s).Nspace2 / 2.0f)) - (*p).y0;
			const deviceFP zo = (*p).z0;
			const deviceFP cB = deviceFPLib::cos((*p).beamAngle);
			const deviceFP cA = deviceFPLib::cos((*p).beamAnglePhi);
			const deviceFP sB = deviceFPLib::sin((*p).beamAngle);
			const deviceFP sA = deviceFPLib::sin((*p).beamAnglePhi);
			const deviceFP x = cB * xo + sA * sB * yo + sA * sB * zo;
			const deviceFP y = cA * yo - sA * zo;
			const deviceFP z = -sB * xo + sA * cB * yo + cA * cB * zo;
			const deviceFP r = deviceFPLib::hypot(x, y);

			const deviceFP wz = (*p).beamwaist * deviceFPLib::sqrt(1.0f + (z * z / (zR * zR)));
			const deviceFP Rz = (z != 0.0f) ? z * (1.0f + (zR * zR / (z * z))) : 1.0e15f;
			const deviceFP phi = deviceFPLib::atan(z / zR);

			deviceComplex Eb = ((*p).beamwaist / wz) 
				* deviceLib::exp(deviceComplex(0.0f, 1.0f) 
					* (ko * (z - (*p).z0) 
						+ ko * r * r / (2.0f * Rz) - phi) 
					- r * r / (wz * wz));
			Eb = Eb * specfac;
			if (isComplexNaN(Eb) || f <= 0.0f) {
				Eb = deviceComplex{};
			}

			field[i] = deviceComplex(
				deviceFPLib::cos((*p).polarizationAngle), 
				-(*p).circularity * deviceFPLib::sin((*p).polarizationAngle)) * Eb;
			field[i + (*s).NgridC] = deviceComplex(
				deviceFPLib::sin((*p).polarizationAngle), 
				(*p).circularity * deviceFPLib::cos((*p).polarizationAngle)) * Eb;
			deviceFP pointEnergy = 
				(modulusSquared(field[i]) + modulusSquared(field[i + (*s).NgridC]));
			pointEnergy *= 
				constProd(lightC<deviceFP>() * eps0<deviceFP>(), 2) * (*s).dx * (*s).dx * (*s).dt;

			//factor 2 accounts for the missing negative frequency plane
			atomicAdd(pulseSum, pointEnergy);
		}
	};

	class multiplyByConstantKernelD {
	public:
		deviceFP* A;
		const deviceFP val;
		deviceFunction void operator()(const int64_t i) const {
			A[i] *= val;
		}
	};

	class multiplicationKernelCompactDoubleVector {
	public:
		const deviceFP* A;
		const deviceComplex* B;
		deviceComplex* C;
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		deviceFunction void operator()(const int64_t i) const {
			const int64_t h = i % (*s).Nfreq; //temporal coordinate
			C[i] = A[h] * B[i];
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
	class expandCylindricalBeam {
	public:
		const deviceParameterSet<deviceFP, deviceComplex>* s;
		
		deviceFunction void operator()(const int64_t i) const {
			const int64_t j = i / (*s).Ntime; //spatial coordinate
			const int64_t k = i % (*s).Ntime; //temporal coordinate

			//positions on the expanded grid corresponding the the current index
			const int64_t pos1 = 2 * ((*s).Nspace - j - 1) * (*s).Ntime + k;
			const int64_t pos2 = (2 * j + 1) * (*s).Ntime + k;

			//reuse memory allocated for the radial Laplacian, casting complex double
			//to a 2x larger double real grid
			deviceFP* expandedBeam1 = (*s).gridRadialLaplacian1;
			deviceFP* expandedBeam2 = expandedBeam1 + 2 * (*s).Ngrid;

			expandedBeam1[pos1] = (*s).gridPolarizationTime1[i];
			expandedBeam1[pos2] = (*s).gridPolarizationTime1[i];
			expandedBeam2[pos1] = (*s).gridPolarizationTime2[i];
			expandedBeam2[pos2] = (*s).gridPolarizationTime2[i];
		}
	};
}
using namespace kernelNamespace;

namespace hostFunctions{
	static simulationParameterSet* fittingSet;
	static ActiveDevice* dFit;

	//forward declarations when useful
	static unsigned long int solveFDTD(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		int64_t tFactor, 
		deviceFP dz, 
		deviceFP frontBuffer, 
		deviceFP backBuffer,
		deviceFP observationPoint = 0.0,
		const std::string& materialMapPath = "",
		bool preserveNearField = false);

	static std::complex<double> hostSellmeierFunc(
		double ls, 
		const double omega, 
		const double* a, 
		const int eqn) {
		const double omega2 = omega * omega;
		double realPart;
		std::complex<double> compPart;
		switch (eqn) {
		case 0:
			if (ls == -a[3] || ls == -a[6] || ls == -a[9] || ls == -a[12]) return std::complex<double>{};
			realPart = a[0]
				+ (a[1] + a[2] * ls) / (ls + a[3])
				+ (a[4] + a[5] * ls) / (ls + a[6])
				+ (a[7] + a[8] * ls) / (ls + a[9])
				+ (a[10] + a[11] * ls) / (ls + a[12])
				+ a[13] * ls
				+ a[14] * ls * ls
				+ a[15] * ls * ls * ls;
			compPart = kLorentzian<double>() * a[16] 
				/ std::complex<double>(a[17] - omega2, a[18] * omega)
				+ kLorentzian<double>() * a[19] 
				/ std::complex<double>(a[20] - omega2, a[21] * omega);
			return std::sqrt(maxN(realPart, 0.9f) + compPart);
		case 1:
			//up to 7 Lorentzian lines
			compPart = a[1] / std::complex<double>(a[2] - omega2, a[3] * omega)
				+ a[4] / std::complex<double>(a[5] - omega2, a[6] * omega)
				+ a[7] / std::complex<double>(a[8] - omega2, a[9] * omega)
				+ a[10] / std::complex<double>(a[11] - omega2, a[12] * omega)
				+ a[13] / std::complex<double>(a[14] - omega2, a[15] * omega)
				+ a[16] / std::complex<double>(a[17] - omega2, a[18] * omega)
				+ a[19] / std::complex<double>(a[20] - omega2, a[21] * omega);
			compPart *= kLorentzian<double>();
			compPart += a[0];
			return std::complex<double>((
				std::sqrt(compPart)).real(), 
				-std::abs((std::sqrt(compPart)).imag()));

		case 100:
		{
			if (ls == -a[3] 
				|| ls == -a[6] 
				|| ls == -a[9] 
				|| ls == -a[12]) return std::complex<double>{};
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
		return std::complex<double>(std::sqrt(maxN(realPart, 0.9f)), 0.0f);
		}
		return std::complex<double>(1.0,0.0);
	};
	static int getTotalSpectrum(ActiveDevice& d) {
		simulationParameterSet* sCPU = d.cParams;
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;

		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.fft((*sc).gridETime1, (*sc).workspace1, deviceFFT::D2Z_1D);
		if ((*sc).is3D) {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum3DKernel{ d.dParamsDevice });
		}
		else if ((*sc).isCylindric) {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel{ d.dParamsDevice });
		}
		else {
			//uncomment and change logic if I want to use the square spectra
			//d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum2DSquareKernel, d.dParamsDevice);
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel{ d.dParamsDevice });
		}
		
		d.deviceMemcpy((double*)(*sCPU).totalSpectrum, 
			(deviceFP*)(*sc).gridPolarizationTime1, 
			3 * (*sCPU).Nfreq * sizeof(double), copyType::ToHost);

		//apply normalization to result of 3D calculation for numerical precision (value may not be
		//represtentable as a float)
		if ((*sCPU).runType != runTypes::counter) {
			double volumeElement = constProd(lightC<double>(), 2 * eps0<double>())
					* (*sCPU).rStep * (*sCPU).tStep * (*sCPU).tStep;
			if((*sCPU).is3D){
				volumeElement *= (*sCPU).rStep;
			} 
			else{
				volumeElement *= vPi<double>();
			}
			for (int64_t i = 0; i < 3 * (*sCPU).Nfreq; i++) {
				(*sCPU).totalSpectrum[i] *= volumeElement;
			}
		}

		return 0;
	}

	static int forwardHankel(ActiveDevice& d, deviceFP* in, deviceComplex* out) {
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		d.deviceLaunch(
			(*sc).Nblock, 
			(*sc).Nthread, 
			hankelKernel{ 
				d.dParamsDevice, 
				in, 
				(deviceFP*)(*sc).workspace1 });
		d.fft((*sc).workspace1, out, deviceFFT::D2Z_1D);
		return 0;
	}
	static int backwardHankel(ActiveDevice& d, deviceComplex* in, deviceFP* out) {
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		d.fft(in, (*sc).workspace1, deviceFFT::Z2D_1D);
		d.deviceLaunch(
			(*sc).Nblock, 
			(*sc).Nthread, 
			inverseHankelKernel{ 
				d.dParamsDevice, 
				(deviceFP*)(*sc).workspace1, 
				out });
		return 0;
	}

	static int addPulseToFieldArrays(
		ActiveDevice& d, 
		pulse<double>& pCPU, 
		const bool useLoadedField, 
		const std::complex<double>* loadedFieldIn) {
		if(pCPU.energy == 0.0) return 0;
		simulationParameterSet* s = d.cParams;
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		deviceParameterSet<deviceFP, deviceComplex>* scDevice = d.dParamsDevice;
		pulse<deviceFP>* p;
		d.deviceCalloc((void**)&p, 1, sizeof(pulse<deviceFP>));
		pulse<deviceFP> devpCPU = pCPU;

		d.deviceMemcpy(
			d.dParamsDevice, 
			sc, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);

		deviceFP* materialPhase;
		deviceComplex* loadedField;

		d.deviceCalloc((void**)&loadedField, (*sc).Ntime, sizeof(deviceComplex));

		//get the material phase
		deviceFP* materialCoefficients, * sellmeierPropagationMedium;

		d.deviceCalloc((void**)&materialCoefficients, 66, sizeof(deviceFP));
		d.deviceCalloc((void**)&sellmeierPropagationMedium, 66, sizeof(deviceFP));
		d.deviceCalloc((void**)&materialPhase, (*s).Nfreq, sizeof(deviceFP));
		d.deviceMemcpy(
			materialCoefficients, 
			(*s).crystalDatabase[pCPU.phaseMaterial].sellmeierCoefficients.data(), 
			66 * sizeof(double), 
			copyType::ToDevice);
		d.deviceMemcpy(
			sellmeierPropagationMedium, 
			(*s).crystalDatabase[(*s).materialIndex].sellmeierCoefficients.data(), 
			66 * sizeof(double), 
			copyType::ToDevice);
		d.deviceLaunch(
			(unsigned int)(*s).Nfreq, 
			1, 
			materialPhaseKernel { 
				(deviceFP)(*s).fStep,
				(*s).Ntime, 
				materialCoefficients, 
				(deviceFP)pCPU.frequency, 
				(deviceFP)pCPU.phaseMaterialThickness,
				(*s).crystalDatabase[pCPU.phaseMaterial].sellmeierType,
				materialPhase });

		deviceFP* pulseSum = &materialCoefficients[0];

		if (useLoadedField) {
			d.deviceMemcpy(
				loadedField, 
				loadedFieldIn, 
				(*s).Nfreq * sizeof(std::complex<double>), 
				copyType::ToDevice);
		}
		d.deviceMemset(pulseSum, 0, sizeof(deviceFP));
		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.deviceMemcpy(
			p, 
			&devpCPU, 
			sizeof(pulse<deviceFP>), 
			copyType::ToDevice);
		if ((*sc).is3D) {
			d.deviceLaunch(
				(*sc).Nblock / 2, 
				(*sc).Nthread, 
				beamGenerationKernel3D{
					(*sc).workspace1, 
					p, 
					pulseSum, 
					scDevice, 
					useLoadedField, 
					loadedField, 
					materialPhase,
					sellmeierPropagationMedium });
		}
		else {
			d.deviceLaunch(
				(*sc).Nblock / 2, 
				(*sc).Nthread, 
				beamGenerationKernel2D{
					(*sc).workspace1, 
					p, 
					pulseSum, 
					scDevice, 
					useLoadedField, 
					loadedField, 
					materialPhase,
					sellmeierPropagationMedium });
		}

		d.fft((*sc).workspace1, (*sc).gridPolarizationTime1, deviceFFT::Z2D_1D);

		d.deviceLaunch(
			2 * (*sc).Nblock, 
			(*sc).Nthread, 
			beamNormalizeKernel{ 
				scDevice, 
				pulseSum, 
				(*sc).gridPolarizationTime1, 
				(deviceFP)pCPU.energy });

		//add the pulses
		d.deviceLaunch(
			2 * (*sc).Nblock, 
			(*sc).Nthread, 
			addDoubleArraysKernel{ 
				(*sc).gridETime1, 
				(deviceFP*)(*sc).gridPolarizationTime1 });

		//fft onto frequency grid
		d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);

		d.deviceFree(materialPhase);
		d.deviceFree(materialCoefficients);
		d.deviceFree(sellmeierPropagationMedium);
		d.deviceFree(loadedField);
		d.deviceFree(p);
		return 0;
	}

	static int addPreviousGridToFieldArrays(
		ActiveDevice& d, 
		double* loadedField) {
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		
		d.deviceMemcpy(
			(deviceFP*)(*sc).gridPolarizationTime1, 
			loadedField, 
			(*sc).Ngrid * 2 * sizeof(double), 
			copyType::ToDevice);

		//add the pulses
		d.deviceLaunch(
			2 * (*sc).Nblock, 
			(*sc).Nthread, 
			addDoubleArraysKernel{ 
				(*sc).gridETime1, 
				(deviceFP*)(*sc).gridPolarizationTime1 });

		//fft onto frequency grid
		d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);

		return 0;
	}
	
	static int prepareElectricFieldArrays(ActiveDevice& d) {

		simulationParameterSet* s = d.cParams;
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		
		d.deviceMemcpy(
			d.dParamsDevice, 
			sc, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);
		deviceParameterSet<deviceFP, deviceComplex>* scDevice = d.dParamsDevice;
		
		if (!(*s).isFollowerInSequence || (*s).isReinjecting) {
			if (!(*s).isReinjecting) {
				d.deviceMemset((*sc).gridETime1, 0, 2 * (*sc).Ngrid * sizeof(deviceFP));
			}
			else {
				d.deviceMemcpy(
					(*sc).gridETime1, 
					(*s).ExtOut, 
					2 * (*s).Ngrid * sizeof(double), 
					copyType::ToDevice);
				d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);
			}
			
			if (d.cParams->pulse1FileType == 3) {
				addPreviousGridToFieldArrays(d, d.cParams->loadedFullGrid1);
			}
			else {
				addPulseToFieldArrays(
					d, 
					d.cParams->pulse1, 
					d.cParams->field1IsAllocated, 
					d.cParams->loadedField1);
			}
			
			if (d.cParams->pulse2FileType == 3) {
				addPreviousGridToFieldArrays(d, d.cParams->loadedFullGrid2);
			}
			else {
				addPulseToFieldArrays(
					d, 
					d.cParams->pulse2, 
					d.cParams->field2IsAllocated, 
					d.cParams->loadedField2);
			}
			
		}
		else {
			d.deviceMemcpy(
				(*sc).gridETime1, 
				(*s).ExtOut, 
				2 * (*s).Ngrid * sizeof(double), 
				copyType::ToDevice);
			d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);
		}
		
		if((*sc).axesNumber==2){
			d.deviceLaunch(
				(*sc).Nblock / 2,
				(*sc).Nthread,
				biaxialRotationKernel {d.dParamsDevice,(*sc).gridEFrequency1,true}
			);
		}

		//Copy the field into the temporary array
		d.deviceMemcpy(
			(void*)(*sc).gridEFrequency1Next1, 
			(void*)(*sc).gridEFrequency1, 
			2 * (*sc).NgridC * sizeof(deviceComplex), 
			copyType::OnDevice);


		

		//set the propagation grids how they should be at the beginning of the next step
		d.deviceLaunch(
			(unsigned int)((*sc).NgridC / minGridDimension), 
			2 * minGridDimension, 
			multiplicationKernelCompactDoubleVector{
				(*sc).fieldFactor1, 
				(*sc).gridEFrequency1Next1, 
				(*sc).workspace1, 
				scDevice });

		return 0;
	}

	static int applyFresnelLoss(
		ActiveDevice& d, 
		simulationParameterSet* s, 
		deviceParameterSet<deviceFP, 
		deviceComplex>& sc, 
		int materialIndex1, 
		int materialIndex2) {
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(
			sellmeierCoefficientsAugmentedCPU, 
			(*s).crystalDatabase[materialIndex1].sellmeierCoefficients.data(), 
			66 * (sizeof(double)));
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
		d.deviceMemcpy(
			sellmeierCoefficients1, 
			sellmeierCoefficientsAugmentedCPU, 
			(66 + 8) * sizeof(double), 
			copyType::ToDevice);

		memcpy(
			sellmeierCoefficientsAugmentedCPU, 
			(*s).crystalDatabase[materialIndex2].sellmeierCoefficients.data(), 
			66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		d.deviceMemcpy(
			sellmeierCoefficients2, 
			sellmeierCoefficientsAugmentedCPU, 
			(66 + 8) * sizeof(double), 
			copyType::ToDevice);
		d.deviceMemcpy(
			sc.gridEFrequency1, 
			(*s).EkwOut, 
			2 * (*s).NgridC * sizeof(std::complex<double>), 
			copyType::ToDevice);

		//transform final result
		d.fft(sc.gridEFrequency1, sc.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(2 * sc.Nblock, sc.Nthread, multiplyByConstantKernelD{
			sc.gridETime1, static_cast<deviceFP>(1.0 / sc.Ngrid) });
		//copy the field arrays from the GPU to CPU memory
		d.deviceMemcpy(
			(*s).ExtOut, 
			sc.gridETime1, 
			2 * (*s).Ngrid * sizeof(double), 
			copyType::ToHost);
		d.deviceMemcpy(
			(*s).EkwOut, 
			sc.gridEFrequency1, 
			2 * (*s).Ngrid * sizeof(std::complex<double>), 
			copyType::ToHost);

		d.deviceFree(sellmeierCoefficients1);
		d.deviceFree(sellmeierCoefficients2);

		return 0;
	}

	static int applyFilter(ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		const double f0, 
		const double bandwidth, 
		const double order, 
		const double inBandAmplitude, 
		const double outOfBandAmplitude) {

		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			filterKernel{
				sDevice,
				static_cast<deviceFP>(1.0e12 * f0),
				static_cast<deviceFP>(1.0e12 * bandwidth),
				static_cast<int>(round(order)),
				static_cast<deviceFP>(inBandAmplitude),
				static_cast<deviceFP>(outOfBandAmplitude) });

		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);

		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);

		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		getTotalSpectrum(d);

		return 0;
	}

	static int applyOptic(
		ActiveDevice& d,
		simulationParameterSet& sCPU,
		const int index,
		const bool applyX,
		const bool applyY){
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			sCPU.ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		std::vector<std::complex<deviceFP>> complexReflectivity = 
			sCPU.optics[index].toComplexSpectrum<deviceFP>(sCPU.Nfreq,sCPU.fStep);
		d.deviceMemcpy(
			d.deviceStruct.gridPropagationFactor1, 
			complexReflectivity.data(), 
			sCPU.Nfreq*sizeof(deviceComplex), 
			copyType::ToDevice
		);
		d.deviceLaunch(
			d.deviceStruct.Nblock/2,
			d.deviceStruct.Nthread,
			applyOpticKernel{d.dParamsDevice,d.deviceStruct.gridPropagationFactor1, applyX, applyY}
		);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			sCPU.EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		d.deviceMemcpy(
			sCPU.ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * sCPU.Ngrid * sizeof(double), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyLorenzian(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		const double amplitude, 
		const double f0, 
		const double gamma, 
		const double radius, 
		const double order) {

		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, lorentzianSpotKernel{
			sDevice,
			(deviceFP)amplitude,
			(deviceFP)(1.0e12 * f0),
			(deviceFP)(1.0e12 * gamma),
			(deviceFP)radius,
			(deviceFP)order });
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ntime) });
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		getTotalSpectrum(d);

		return 0;
	}

	static int applyAperatureFarFieldHankel(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		double diameter, 
		double activationParameter, 
		double xOffset, 
		double yOffset) {
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		forwardHankel(d, d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			apertureFarFieldKernelHankel{
				sDevice,
				(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter),
				(deviceFP)activationParameter,
				(deviceFP)(deg2Rad<deviceFP>() * xOffset),
				(deviceFP)(deg2Rad<deviceFP>() * yOffset) });
		backwardHankel(d, d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyInverseAperatureFarFieldHankel(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		double diameter, 
		double activationParameter, 
		double xOffset, 
		double yOffset) {
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		forwardHankel(d, d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			inverseApertureFarFieldKernelHankel{
				sDevice,
				(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter),
				(deviceFP)activationParameter,
				(deviceFP)(deg2Rad<deviceFP>() * xOffset),
				(deviceFP)(deg2Rad<deviceFP>() * yOffset) });
		backwardHankel(d, d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyAperatureFarField(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		double diameter, 
		double activationParameter, 
		double xOffset, 
		double yOffset) {
		if ((*sCPU).isCylindric) {
			return applyAperatureFarFieldHankel(d, sCPU, diameter, activationParameter, xOffset, yOffset);
		}
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, apertureFarFieldKernel{
			sDevice,
			(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter),
			(deviceFP)activationParameter,
			(deviceFP)(deg2Rad<deviceFP>() * xOffset),
			(deviceFP)(deg2Rad<deviceFP>() * yOffset) });

		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);

		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);

		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		getTotalSpectrum(d);
		return 0;
	}

	static int applyInverseAperatureFarField(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		double diameter, 
		double activationParameter, 
		double xOffset, 
		double yOffset) {
		if ((*sCPU).isCylindric) {
			return applyInverseAperatureFarFieldHankel(d, sCPU, diameter, activationParameter, xOffset, yOffset);
		}
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, inverseApertureFarFieldKernel{
			sDevice,
			(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter),
			(deviceFP)activationParameter,
			(deviceFP)(deg2Rad<deviceFP>() * xOffset),
			(deviceFP)(deg2Rad<deviceFP>() * yOffset) });

		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);

		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);

		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		getTotalSpectrum(d);
		return 0;
	}

	static int applyAperature(
		ActiveDevice& d, 
		const simulationParameterSet* sCPU, 
		const double diameter, 
		const double activationParameter) {
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);

		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(
			d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			apertureKernel{
				sDevice,
				(deviceFP)(0.5 * diameter),
				(deviceFP)(activationParameter) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applySphericalMirror(
		ActiveDevice& d, 
		const simulationParameterSet* sCPU, 
		deviceParameterSet<deviceFP, deviceComplex>& s, 
		const double ROC) {

		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceMemcpy(
			sDevice, 
			&s, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);

		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			sphericalMirrorKernel{ sDevice, (deviceFP)ROC });
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.deviceLaunch(
			2 * d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			multiplyByConstantKernelD{ 
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ntime )});
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyParabolicMirror(ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		deviceParameterSet<deviceFP, deviceComplex>& s, 
		const double focus) {

		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			parabolicMirrorKernel{ sDevice, (deviceFP)focus });
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.deviceLaunch(
			2 * d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			multiplyByConstantKernelD { 
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ntime) });
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyLinearPropagation(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		const int materialIndex, 
		const double thickness) {

		if (d.hasPlasma) {
			simulationParameterSet sCopy = *sCPU;
			sCopy.nonlinearAbsorptionStrength = 0.0;
			d.reset(&sCopy);
		}

		d.deviceMemcpy(
			d.deviceStruct.gridEFrequency1, 
			(*sCPU).EkwOut, 
			d.deviceStruct.NgridC * 2 * sizeof(std::complex<double>), 
			copyType::ToDevice);

		deviceFP* sellmeierCoefficients = (deviceFP*)d.deviceStruct.gridEFrequency1Next1;
		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(
			sellmeierCoefficientsAugmentedCPU, 
			(*sCPU).crystalDatabase[materialIndex].sellmeierCoefficients.data(), 
			66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*sCPU).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*sCPU).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*sCPU).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*sCPU).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*sCPU).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*sCPU).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		d.deviceMemcpy(
			sellmeierCoefficients, 
			sellmeierCoefficientsAugmentedCPU, 
			(66 + 8) * sizeof(double), 
			copyType::ToDevice);
		d.deviceStruct.axesNumber = (*sCPU).crystalDatabase[materialIndex].axisType;
		d.deviceStruct.sellmeierType = (*sCPU).crystalDatabase[materialIndex].sellmeierType;
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;

		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			applyLinearPropagationKernel{ 
				sellmeierCoefficients, 
				(deviceFP)thickness, 
				sDevice });
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			d.deviceStruct.NgridC * 2 * sizeof(std::complex<double>), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(
			2 * d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			multiplyByConstantKernelD{ 
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ngrid) });

		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);
		getTotalSpectrum(d);

		return 0;
	}

	static int preparePropagationGrids(
		ActiveDevice& d) {

		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		simulationParameterSet* s = d.cParams;
		deviceFP* sellmeierCoefficients = (deviceFP*)(*sc).gridEFrequency1Next1;
		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[79];
		memcpy(
			sellmeierCoefficientsAugmentedCPU, 
			(*s).sellmeierCoefficients, 
			66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		memcpy(
			sellmeierCoefficientsAugmentedCPU + 72, 
			(*s).crystalDatabase[(*s).materialIndex].nonlinearReferenceFrequencies.data(), 
			7 * sizeof(double));
		d.deviceMemcpy(
			sellmeierCoefficients, 
			sellmeierCoefficientsAugmentedCPU, 
			79 * sizeof(double), 
			copyType::ToDevice);

		//prepare the propagation grids
		deviceParameterSet<deviceFP, deviceComplex>* sD = d.dParamsDevice;
		d.deviceMemcpy(
			sD, 
			sc, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);
		d.deviceLaunch(
			(unsigned int)maxN((*sc).Ntime / 2, 82), 
			1, 
			getChiLinearKernel{ sD, sellmeierCoefficients });
		if ((*s).is3D) {
			d.deviceLaunch(
				static_cast<unsigned int>((*sc).Nblock / 2u), 
				static_cast<unsigned int>((*sc).Nthread), 
				prepare3DGridsKernel{ sellmeierCoefficients, sD });
		}
		else if ((*s).isCylindric) {
			d.deviceLaunch(
				static_cast<unsigned int>((*sc).Nblock / 2u), 
				static_cast<unsigned int>((*sc).Nthread), 
				prepareCylindricGridsKernel{ sellmeierCoefficients, sD });
		}
		else {
			d.deviceLaunch(
				static_cast<unsigned int>((*sc).Nblock / 2u), 
				static_cast<unsigned int>((*sc).Nthread), 
				prepareCartesianGridsKernel{ sellmeierCoefficients, sD });
		}
		d.deviceMemcpy(
			sc, 
			sD, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToHost);
		return 0;
	}

	//Rotate the field on the GPU
	//Allocates memory and copies from CPU, then copies back to CPU and deallocates
	// - inefficient but the general principle is that only the CPU memory is preserved
	// after simulations finish... and this only runs at the end of the simulation
	static int rotateField(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		const double rotationAngle) {

		deviceComplex* Ein1 = d.deviceStruct.gridEFrequency1;
		deviceComplex* Ein2 = d.deviceStruct.gridEFrequency2;
		deviceComplex* Eout1 = d.deviceStruct.gridEFrequency1Next1;
		deviceComplex* Eout2 = d.deviceStruct.gridEFrequency1Next2;

		//retrieve/rotate the field from the CPU memory
		d.deviceMemcpy(
			Ein1, 
			(*sCPU).EkwOut, 
			2 * (*sCPU).NgridC * sizeof(std::complex<double>), 
			copyType::ToDevice);
		d.deviceLaunch(
			static_cast<unsigned int>(d.deviceStruct.NgridC / minGridDimension), 
			minGridDimension, 
			rotateFieldKernel{ Ein1, Ein2, Eout1, Eout2, (deviceFP)rotationAngle });
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			Eout1, 
			2 * (*sCPU).NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);

		//transform back to time
		d.fft(Eout1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(
			2 * d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			multiplyByConstantKernelD{ 
				d.deviceStruct.gridETime1, 
				static_cast<deviceFP>(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		//update spectrum
		getTotalSpectrum(d);
		return 0;
	}

	static void savePlasma(
		ActiveDevice& d, 
		int64_t saveloc, 
		int densityLoc){

		const deviceParameterSet<deviceFP, deviceComplex>* sH = d.s; 
		const deviceParameterSet<deviceFP, deviceComplex>* sD = d.dParamsDevice;
		//prepare the grids
		preparePropagationGrids(d);
		prepareElectricFieldArrays(d);

		//Clear out the memory to which the data will be saved
		d.deviceMemset(d.deviceStruct.gridEFrequency1,0,4*sizeof(deviceFP)*d.cParams->NgridC);

		//run the plasma kernels
		d.deviceLaunch(
			(*sH).Nblock, 
			(*sH).Nthread, 
			plasmaCurrentKernel_twoStage_A{ sD });
		d.deviceLaunch(
			(unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 
			minGridDimension, 
			plasmaCurrentKernel_SaveOutput{ sD });

		//save to memory
		d.deviceMemcpy(
			d.cParams->ExtOut + 2*d.cParams->Ngrid*saveloc, 
			d.deviceStruct.gridPolarizationTime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);

		if (densityLoc == 1) {
			d.deviceMemcpy(
				d.cParams->ExtOut + 2 * d.cParams->Ngrid * saveloc, 
				reinterpret_cast<deviceFP*>(d.deviceStruct.gridEFrequency1), 
				d.deviceStruct.Ngrid * sizeof(double), 
				copyType::ToHost);
		}
		else if (densityLoc == 2) {
			d.deviceMemcpy(
				d.cParams->ExtOut + 2 * d.cParams->Ngrid * saveloc + d.cParams->Ngrid, 
				reinterpret_cast<deviceFP*>(d.deviceStruct.gridEFrequency1), 
				d.deviceStruct.Ngrid * sizeof(double), 
				copyType::ToHost);
		}
		
	}
//function to run a RK4 time step
//stepNumber is the sub-step index, from 0 to 3
	static int runRK4Step(
		ActiveDevice& d, 
		const uint8_t stepNumber){

		const deviceParameterSet<deviceFP, deviceComplex>* sH = d.s; 
		const deviceParameterSet<deviceFP, deviceComplex>* sD = d.dParamsDevice;

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
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					nonlinearPolarizationKernel{ sD });
				if((*sH).hasPlasma){
					d.deviceLaunch(
						(*sH).Nblock, 
						(*sH).Nthread, 
						plasmaCurrentKernel_twoStage_A{ sD });
					
					//CUDA and other platforms perform very 
					//differently for different versions of this kernel, use optimum per platform
					#ifdef __CUDACC__
					d.deviceLaunch(
						(unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 
						2 * minGridDimension, 
						plasmaCurrentKernel_twoStage_B{ sD });
					#else
					d.deviceLaunch(
						(unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 
						minGridDimension, 
						plasmaCurrentKernel_twoStage_B_simultaneous{ sD });
					#endif
					
					d.deviceLaunch(
						(*sH).Nblock, 
						(*sH).Nthread, 
						expandCylindricalBeam{ sD });
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFT::D2Z_Polarization);
					d.deviceLaunch(
						(*sH).Nblock / 2, 
						(*sH).Nthread, 
						updateKwithPlasmaKernelCylindric{ sD });
				}
				else{
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFT::D2Z_Polarization);
					d.deviceLaunch(
						(*sH).Nblock / 2, 
						(*sH).Nthread, 
						updateKwithPolarizationKernelCylindric{ sD });
				}
			}
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				radialLaplacianKernel{ sD });
			d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFT::D2Z);

			switch (stepNumber) {
			case 0:
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					rkKernel0Cylindric{ sD });
				break;
			case 1:
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					rkKernel1Cylindric{ sD });
				break;
			case 2:
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					rkKernel2Cylindric{ sD });
				break;
			case 3:
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					rkKernel3Cylindric{ sD });
				break;
			}
			return 0;
		}
		// 2D and 3D cartesian
		// Only one type of FFT
		// currently nonlinear polarization and plasma ffts are not batched, could give
		// minor speed boost by combining them, but requires additional memory, so
		// nonlinear polarization and plasma are fft-ed separately to accommodate larger
		// systems.
		else if ((*sH).isNonLinear) {
			//perform inverse FFT to get time-space electric field
			if((*sH).axesNumber==2){
				//undo delta rotation if biaxial
				d.deviceLaunch(
					(*sH).Nblock / 2,
					(*sH).Nthread,
					biaxialRotationKernel {sD,(*sH).workspace1,false}
				);
			}
			d.fft((*sH).workspace1, (*sH).gridETime1, deviceFFT::Z2D);
			//Plasma/multiphoton absorption
			if ((*sH).hasPlasma) {
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					plasmaCurrentKernel_twoStage_A{ sD });
				d.deviceLaunch(
					(unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 
					minGridDimension, 
					plasmaCurrentKernel_twoStage_B_simultaneous{ sD });
				d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFT::D2Z);
				d.deviceLaunch(
					(*sH).Nblock / 2, 
					(*sH).Nthread, 
					updateKwithPlasmaKernel{ sD });
			}
			//Nonlinear polarization
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				nonlinearPolarizationKernel{ sD });
			d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFT::D2Z);
			//rotate nonlinear polarization into the delta frame if biaxial
			if((*sH).axesNumber==2){
				d.deviceLaunch(
					(*sH).Nblock / 2,
					(*sH).Nthread,
					biaxialRotationKernel {sD,(*sH).workspace1,true}
				);
			}
		}

		//advance an RK4 step
		switch (stepNumber) {
		case 0:
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				rkKernel0{ sD });
			break;
		case 1:
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				rkKernel1{ sD });
			break;
		case 2:
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				rkKernel2{ sD });
			break;
		case 3:
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				rkKernel3{ sD });
			break;
		}
		return 0;
	}

	static unsigned long int solveNonlinearWaveEquationWithDevice(
		ActiveDevice& d, 
		simulationParameterSet* sCPU) {

		if (sCPU->isFDTD) {
			return solveFDTD(d, sCPU, 5, 0, 1e-6, 1e-6);
		}
		//prepare the propagation arrays
		preparePropagationGrids(d);
		prepareElectricFieldArrays(d);

		deviceFP* canaryPointer = 
			&d.deviceStruct.gridETime1[
				d.deviceStruct.Ntime / 2 
					+ d.deviceStruct.Ntime 
					* (d.deviceStruct.Nspace / 2 
						+ d.deviceStruct.Nspace 
						* (d.deviceStruct.Nspace2 / 2))];
		//Core propagation loop
		for (int64_t i = 0; i < d.deviceStruct.Nsteps; ++i) {

			//RK4
			runRK4Step(d, 0);
			runRK4Step(d, 1);
			runRK4Step(d, 2);
			runRK4Step(d, 3);

			//periodically check if the simulation diverged or was cancelled
			if ((*sCPU).cancellationCalled) break;
			if (i % 10 == 0 && d.isTheCanaryPixelNaN(canaryPointer)) break;
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
		}
		if ((*sCPU).isInFittingMode && !(*sCPU).isInSequence)(*(*sCPU).progressCounter)++;

		if(d.deviceStruct.axesNumber==2){
			d.deviceLaunch(
				d.s->Nblock / 2,
				d.s->Nthread,
				biaxialRotationKernel {d.dParamsDevice,d.deviceStruct.gridEFrequency1,false}
			);
		}
		//take final spectra and transfer the results to the CPU
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{ 
				d.deviceStruct.gridETime1, 
				static_cast<deviceFP>(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);
		getTotalSpectrum(d);

		return 13 * d.isTheCanaryPixelNaN(canaryPointer);
	}

	template<typename maxwellType>
	static unsigned int freeFDTD(
		ActiveDevice& d, 
		maxwellType& maxCalc) {
		unsigned int errorValue = 13 * d.isTheCanaryPixelNaN(&(maxCalc.Egrid[0].y));
		d.deviceFree(maxCalc.Egrid);
		d.deviceFree(maxCalc.EgridEstimate);
		d.deviceFree(maxCalc.EgridEstimate2);
		d.deviceFree(maxCalc.EgridNext);
		d.deviceFree(maxCalc.Hgrid);
		d.deviceFree(maxCalc.HgridEstimate);
		d.deviceFree(maxCalc.HgridEstimate2);
		d.deviceFree(maxCalc.HgridNext);
		d.deviceFree(maxCalc.materialGrid);
		d.deviceFree(maxCalc.materialGridNext);
		d.deviceFree(maxCalc.materialGridEstimate);
		d.deviceFree(maxCalc.materialGridEstimate2);
		d.deviceFree(maxCalc.deviceCopy);
		if(maxCalc.hasMaterialMap){
			d.deviceFree(maxCalc.materialMap);
			d.deviceFree(maxCalc.oscillatorIndexMap);
		}
		return errorValue;
	}


	template<typename maxwellType>
	static void calculateFDTDParameters(
		const simulationParameterSet* sCPU, 
		maxwellType& maxCalc,
		int mapValue = 0){

		if(mapValue == 0){
			double n0 = hostSellmeierFunc(
			0, twoPi<double>() * sCPU->pulse1.frequency, 
			(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			double nm1 = hostSellmeierFunc(
				0, twoPi<double>() * (-2e11 + sCPU->pulse1.frequency), 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			double np1 = hostSellmeierFunc(
				0, twoPi<double>() * (2e11 + sCPU->pulse1.frequency), 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			double nGroup = n0 + sCPU->pulse1.frequency * (np1 - nm1) / 4.0e11;
			maxCalc.waitFrames = 
				(maxCalc.frontBuffer + nGroup * maxCalc.crystalThickness + 10 * maxCalc.zStep) 
				/ (lightC<double>() * maxCalc.tStep);
			maxCalc.waitFrames = maxCalc.tGridFactor * (maxCalc.waitFrames / maxCalc.tGridFactor);
			maxCalc.Nt = maxCalc.waitFrames + (*sCPU).Ntime * maxCalc.tGridFactor;
		}
		

		//copy the crystal info
		if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].axisType == 0) {
			//isotropic
			for (int i = 0; i < 22; i++) {
				maxCalc.sellmeierEquations[i][mapValue].x = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
				maxCalc.sellmeierEquations[i][mapValue].y = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
				maxCalc.sellmeierEquations[i][mapValue].z = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
			}
		}
		else if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].axisType == 1) {
			//uniaxial
			for (int i = 0; i < 22; i++) {
				maxCalc.sellmeierEquations[i][mapValue].x = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
				maxCalc.sellmeierEquations[i][mapValue].y = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
				maxCalc.sellmeierEquations[i][mapValue].z = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i+22]);
			}
		}
		else {
			//biaxial
			for (int i = 0; i < 22; i++) {
				maxCalc.sellmeierEquations[i][mapValue] = maxwellPoint<deviceFP>{
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]),
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i+22]),
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i+44])
				};
			}
		}
		
		if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearSwitches[0]) 
			maxCalc.hasChi2[mapValue] = true;
		if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearSwitches[1] == 1) 
			maxCalc.hasFullChi3[mapValue] = true;
		if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearSwitches[1] == 2) 
			maxCalc.hasSingleChi3[mapValue] = true;

		//perform millers rule normalization on the nonlinear coefficients while copying the values
		double millersRuleFactorChi2 = 2e-12;
		double millersRuleFactorChi3 = 1.0;
		if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[0]) {
			double wRef = 
				twoPi<double>() * 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[0];
			double nRef = hostSellmeierFunc(
				0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi2 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[1];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi2 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[2];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi2 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[3];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi3 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[4];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi3 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[5];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi3 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[6];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi3 /= nRef * nRef - 1.0;
		}
		for (int i = 0; i < 6; i++) {
			maxCalc.chi2[i][mapValue] = maxwellPoint<deviceFP>{
				static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].d[3 * i] * millersRuleFactorChi2),
				static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].d[3 * i + 1] * millersRuleFactorChi2),
				static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].d[3 * i + 2] * millersRuleFactorChi2) };
			if (i > 2) maxCalc.chi2[i][mapValue] *= 2.0;
		}
		for (int i = 0; i < 27; i++) {
			maxCalc.chi3[i][mapValue].x = static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].chi3[i] * millersRuleFactorChi3);
			maxCalc.chi3[i][mapValue].y = static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].chi3[i+27] * millersRuleFactorChi3);
			maxCalc.chi3[i][mapValue].z = static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].chi3[i+54] * millersRuleFactorChi3);
		}

		//count the nonzero oscillators
		if(!maxCalc.hasMaterialMap){
			for (int i = 0; i < 7; i++) {
				if (maxCalc.sellmeierEquations[1 + i * 3][mapValue].x > 0.0) maxCalc.Noscillators++;
				else break;
			}

			//collect plasma properties
			if ((*sCPU).nonlinearAbsorptionStrength > 0.0) {
				maxCalc.Noscillators++;
				maxCalc.hasPlasma[mapValue] = true;
				maxCalc.kCarrierGeneration[mapValue] = 2.0 / ((*sCPU).bandGapElectronVolts);
				maxCalc.kDrude[mapValue] = - eps0<double>() * elCharge<double>() 
					/ ((*sCPU).effectiveMass * elMass<double>());
				maxCalc.gammaDrude[mapValue] = (*sCPU).drudeGamma;
				maxCalc.kNonlinearAbsorption[mapValue] = 0.5 * (*sCPU).nonlinearAbsorptionStrength;
				maxCalc.nonlinearAbsorptionOrder[mapValue] = static_cast<int>(
					std::ceil(eVtoHz<double>() * (*sCPU).bandGapElectronVolts 
						/ (*sCPU).pulse1.frequency)) - 1;
			}
			maxCalc.NMaterialGrid =
				(maxCalc.materialStop - maxCalc.materialStart)
				* maxCalc.Nx * maxCalc.Ny * maxCalc.Noscillators;
		}
	}

	static void prepareFDTD(
		ActiveDevice& d, 
		const simulationParameterSet* sCPU, 
		maxwell3D& maxCalc,
		const std::string& materialMapPath = "") {

		//Check if there is a material map and allocate/load it if necessary
		if (materialMapPath != "" && (*sCPU).runType != runTypes::counter && (*sCPU).runType != runTypes::cluster) {
			//throw std::runtime_error(std::string("path string:\n").append(materialMapPath));
			maxCalc.hasMaterialMap = false;
			int64_t NmaterialPoints = 0;
			std::vector<int8_t> materialMapCPU(maxCalc.Ngrid, 0);
			std::vector<int64_t> oscillatorIndexMapCPU(maxCalc.Ngrid, 0);
			std::ifstream fs(materialMapPath);

			if (fs.good()) {
				auto moveToColon = [&]() {
					char x = 0;
					while (x != ':' && fs.good()) {
						fs >> x;
					}
					return 0;
					};
				//First line: materials
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.materialKeys[i];
					//optimization: count oscillators
				}
				maxCalc.Noscillators = 7;
				//Second line: theta
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.materialTheta[i];
				}
				//Third line: phi
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.materialPhi[i];
				}
				//Fourth: Nonlinear absorption (NOTE: FIX PLASMA PARAMS)
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.kNonlinearAbsorption[i];
					maxCalc.hasPlasma[i] = maxCalc.kNonlinearAbsorption[i] > 0.0;
				}
				//Fifth: bandgap
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.kDrude[i];
				}
				//Sixth: gamma
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.gammaDrude[i];
				}
				//Seventh: effective mass
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.kCarrierGeneration[i];
				}
				//Eighth: start data
				int64_t gridCount{};
				while (fs.good() && gridCount < maxCalc.Ngrid) {
					int currentInt;
					fs >> currentInt;
					materialMapCPU[gridCount] = static_cast<int8_t>(currentInt);
					if (materialMapCPU[gridCount] > 0) {
						oscillatorIndexMapCPU[gridCount] = NmaterialPoints;
						NmaterialPoints++;
					}
					gridCount++;
				}

				d.deviceCalloc((void**)&(maxCalc.materialMap),
					maxCalc.Ngrid, sizeof(char));
				d.deviceMemcpy(
					(void*)maxCalc.materialMap,
					(void*)materialMapCPU.data(),
					maxCalc.Ngrid * sizeof(char),
					copyType::ToDevice);
				d.deviceCalloc((void**)&(maxCalc.oscillatorIndexMap),
					maxCalc.Ngrid, sizeof(int64_t));
				d.deviceMemcpy(
					(void*)maxCalc.oscillatorIndexMap,
					(void*)oscillatorIndexMapCPU.data(),
					maxCalc.Ngrid * sizeof(int64_t),
					copyType::ToDevice);
				maxCalc.hasMaterialMap = true;
				maxCalc.NMaterialGrid = NmaterialPoints * maxCalc.Noscillators;

				for(int i = 0; i<NmaterialMax; i++){
					calculateFDTDParameters(sCPU, maxCalc, i);
				}
				if(maxCalc.observationPoint==0){
					maxCalc.observationPoint = maxCalc.Nz - 10;
				}
				
			}
			else{
				std::runtime_error("Failed to load material map.\n");
			}

		}
		else{
			maxCalc.materialKeys[0] = (*sCPU).materialIndex;
			calculateFDTDParameters(sCPU, maxCalc);
		}

		
		//Make sure that the time grid is populated and do a 1D (time) FFT onto the frequency grid
		//make both grids available through the maxCalc class
		d.deviceMemcpy(
			d.s->gridETime1, 
			(*sCPU).ExtOut, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToDevice);
		maxCalc.inOutEx = d.s->gridETime1;
		maxCalc.inOutEy = d.s->gridETime2;
		d.fft(d.s->gridETime1, d.s->gridEFrequency1, deviceFFT::D2Z_1D);
		maxCalc.inputExFFT = reinterpret_cast<deviceFP*>(d.s->gridEFrequency1);
		maxCalc.inputEyFFT = reinterpret_cast<deviceFP*>(d.s->gridEFrequency2);
		d.deviceMemset(maxCalc.inOutEx, 0, 2*(*sCPU).Ngrid * sizeof(deviceFP));

		//allocate the new memory needed for the maxwell calculation
		d.deviceCalloc((void**)&(maxCalc.Egrid), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.EgridEstimate), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.EgridNext), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.EgridEstimate2), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.Hgrid), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.HgridEstimate), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.HgridNext), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.HgridEstimate2), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.materialGrid), 
			maxCalc.NMaterialGrid, sizeof(oscillator<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.materialGridEstimate), 
			maxCalc.NMaterialGrid, sizeof(oscillator<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.materialGridEstimate2), 
			maxCalc.NMaterialGrid, sizeof(oscillator<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.materialGridNext), 
			maxCalc.NMaterialGrid, sizeof(oscillator<deviceFP>));

		//make a device copy of the maxCalc class
		maxwell3D* maxCalcDevice{};
		d.deviceCalloc((void**)&maxCalcDevice, 1, sizeof(maxwell3D));
		d.deviceMemcpy(
			(void*)maxCalcDevice, 
			(void*)&maxCalc, 
			sizeof(maxwell3D), 
			copyType::ToDevice);
		maxCalc.deviceCopy = maxCalcDevice;
	}

	static unsigned long int solveFDTD(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		int64_t tFactor, 
		deviceFP dz, 
		deviceFP frontBuffer, 
		deviceFP backBuffer,
		deviceFP observationPoint,
		const std::string& materialMapPath,
		bool preserveNearField) {
		
		//initialize the grid if necessary
		if (!sCPU->isFollowerInSequence) {
			simulationParameterSet sCPUcopy = *sCPU;

			sCPUcopy.materialIndex = 0;
			sCPUcopy.crystalTheta = 0.0;
			sCPUcopy.crystalPhi = 0.0;
			sCPUcopy.crystalThickness = 0.0;
			sCPUcopy.propagationStep = 1e-9;

			sCPUcopy.nonlinearAbsorptionStrength = 0.0;
			sCPUcopy.chi2Tensor = sCPUcopy.crystalDatabase[0].d.data();
			sCPUcopy.chi3Tensor = sCPUcopy.crystalDatabase[0].chi3.data();
			sCPUcopy.nonlinearSwitches = 
				sCPUcopy.crystalDatabase[0].nonlinearSwitches.data();
			sCPUcopy.absorptionParameters = 
				sCPUcopy.crystalDatabase[0].absorptionParameters.data();
			sCPUcopy.sellmeierCoefficients = 
				sCPUcopy.crystalDatabase[0].sellmeierCoefficients.data();

			sCPUcopy.sellmeierType = sCPUcopy.crystalDatabase[0].sellmeierType;
			sCPUcopy.axesNumber = 0;
			sCPUcopy.isFDTD = false;
			d.reset(&sCPUcopy);
			solveNonlinearWaveEquationWithDevice(d, &sCPUcopy);
			d.reset(sCPU);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			(*sCPU).isFollowerInSequence = true;
		}
		if (dz == 0.0) dz = (*sCPU).propagationStep;
		//generate the FDTD data structure and prepare the device
		maxwell3D maxCalc = maxwell3D(sCPU, tFactor, dz, frontBuffer, backBuffer);
		if(observationPoint != 0.0){
			maxCalc.observationPoint = static_cast<int>(round(observationPoint/dz));
			if(maxCalc.observationPoint < 1 || maxCalc.observationPoint >= maxCalc.Nz){
				throw std::runtime_error("Invalid observation point in FDTD\n");
			}
		}
		prepareFDTD(d, sCPU, maxCalc, materialMapPath);
		
		//RK loop
		for (int64_t i = 0; i < maxCalc.Nt; i++) {
			d.deviceLaunch(
				maxCalc.Ngrid / 64, 
				64, 
				maxwellRKkernel0{ maxCalc.deviceCopy, i });
			d.deviceLaunch(
				maxCalc.Ngrid / 64, 
				64, 
				maxwellRKkernel1{ maxCalc.deviceCopy, i });
			d.deviceLaunch(
				maxCalc.Ngrid / 64,
				64, 
				maxwellRKkernel2{ maxCalc.deviceCopy, i });
			d.deviceLaunch(
				maxCalc.Ngrid / 64, 
				64, 
				maxwellRKkernel3{ maxCalc.deviceCopy, i });
			if (i % maxCalc.tGridFactor == 0 
				&& (i >= maxCalc.waitFrames)) {
				d.deviceLaunch(
					(maxCalc.Nx * maxCalc.Ny) / minGridDimension, 
					minGridDimension, 
					maxwellSampleGrid{ 
						maxCalc.deviceCopy, 
						(i - maxCalc.waitFrames) / maxCalc.tGridFactor });
			}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			if (i % 20 == 0 && d.isTheCanaryPixelNaN(&(maxCalc.Egrid[0].y))) break;
			if ((*sCPU).cancellationCalled) break;
		}
		

		
		
		//correct far-field amplitudes for vectorial effects
		if(preserveNearField){
			d.deviceMemcpy(
				(*sCPU).ExtOut, 
				maxCalc.inOutEx, 
				2*(*sCPU).Ngrid * sizeof(double), 
				copyType::ToHost);
			d.fft(maxCalc.inOutEx, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
			d.deviceMemcpy(
				(*sCPU).EkwOut,
				d.deviceStruct.gridEFrequency1,
				2 * d.deviceStruct.NgridC * sizeof(std::complex<double>),
				copyType::ToHost);
		}
		else {
			d.fft(maxCalc.inOutEx, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
			if((*sCPU).is3D){
				d.deviceLaunch(
					d.deviceStruct.Nblock / 2, 
					d.deviceStruct.Nthread,
					correctFDTDAmplitudesKernel{ d.dParamsDevice });
			}
			else{
				d.deviceLaunch(
					d.deviceStruct.Nblock / 2, 
					d.deviceStruct.Nthread,
					correctFDTDAmplitudesKernel2D{ d.dParamsDevice });
			}
			d.deviceMemcpy(
				(*sCPU).EkwOut,
				d.deviceStruct.gridEFrequency1,
				2 * d.deviceStruct.NgridC * sizeof(std::complex<double>),
				copyType::ToHost);
			d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
			//transfer result to CPU memory and take spectrum
			d.deviceMemcpy(
				(*sCPU).ExtOut, 
				d.deviceStruct.gridETime1, 
				2*(*sCPU).Ngrid * sizeof(double), 
				copyType::ToHost);
		}

		getTotalSpectrum(d);

		//free device memory		
		return freeFDTD(d, maxCalc);
	}

	static constexpr unsigned int functionID(const char* s) {
		unsigned int result = 0;
		for (int i = 0; s[i]; ++i) {
			result += s[i];
		}
		return result;
	}

	//Dispatcher of the sequence mode. 
	// New functions go here, and should have a unique ID (chances of a conflict are small, and 
	// will produce a compile-time error.
	// Functions cannot start with a number or the string "None".
	static int interpretCommand(
		const std::string& cc, 
		const double* iBlock, 
		double* vBlock, 
		ActiveDevice& d, 
		simulationParameterSet *sCPU) {
		if (cc.size() == 0) return 0;
		crystalEntry* db = (*sCPU).crystalDatabase;
		int error = 0;
		double parameters[32] = {};
		bool defaultMask[32] = {};
		auto functionNameEnd = cc.find_first_of('(');
		if (functionNameEnd ==std::string::npos) throw std::runtime_error(std::string("Didn't find opening parenthesis in\n").append(cc));
		switch (functionID(cc.substr(0, functionNameEnd).c_str())) {
		case functionID("rotate"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			rotateField(d, sCPU, deg2Rad<deviceFP>() * parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("set"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			if(parameters[0]<100 && parameters[0] >= 0){ 
				vBlock[static_cast<int>(parameters[0])] = parameters[1];
			}
			else throw std::runtime_error("set() index must be less than 100\n");
			break;
		case functionID("plasmaReinject"):
			(*sCPU).isReinjecting = true;
			[[fallthrough]];
		case functionID("plasma"):
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
		case functionID("nonlinear"):
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
		case functionID("fdtd"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			error = solveFDTD(
				d, 
				sCPU, 
				static_cast<int64_t>(parameters[0]), 
				parameters[1], 
				parameters[2], 
				parameters[3]);
			break;
		case functionID("fdtdGrid"):
			{std::string filepath = cc.substr(cc.find('"')+1, cc.find('"', cc.find('"') + 1) - cc.find('"') - 1);
			std::string newParameterString =cc.substr(cc.find('"', cc.find('"') + 1)+1,std::string::npos);
			newParameterString[0] = '(';
			interpretParameters(newParameterString, 2, iBlock, vBlock, parameters, defaultMask);}
			error = solveFDTD(d,
				sCPU,
				static_cast<int>(parameters[0]),
				(*sCPU).propagationStep,
				0.0,
				0.0,
				parameters[1],
				cc.substr(cc.find('"')+1, cc.find('"', cc.find('"') + 1) - cc.find('"') - 1));
			break;
		case functionID("fdtdGridNearField"):
			{std::string filepath = cc.substr(cc.find('"')+1, cc.find('"', cc.find('"') + 1) - cc.find('"') - 1);
			std::string newParameterString =cc.substr(cc.find('"', cc.find('"') + 1)+1,std::string::npos);
			newParameterString[0] = '(';
			interpretParameters(newParameterString, 2, iBlock, vBlock, parameters, defaultMask);}
			error = solveFDTD(d,
				sCPU,
				static_cast<int>(parameters[0]),
				(*sCPU).propagationStep,
				0.0,
				0.0,
				parameters[1],
				cc.substr(cc.find('"')+1, cc.find('"', cc.find('"') + 1) - cc.find('"') - 1),
				true);
			break;
		case functionID("default"):
			d.reset(sCPU);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU);
			break;
		case functionID("save"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			{
				int64_t saveLoc = (int64_t)parameters[0];
				if (saveLoc < (*sCPU).Nsims && saveLoc != 0 && (*sCPU).runType != runTypes::counter) {
					memcpy(
						&(*sCPU).ExtOut[saveLoc * (*sCPU).Ngrid * 2], 
						(*sCPU).ExtOut, 
						2 * (*sCPU).Ngrid * sizeof(double));
					memcpy(
						&(*sCPU).EkwOut[saveLoc * (*sCPU).NgridC * 2], 
						(*sCPU).EkwOut, 
						2 * (*sCPU).NgridC * sizeof(std::complex<double>));
					memcpy(
						&(*sCPU).totalSpectrum[saveLoc * 3 * (*sCPU).Nfreq], 
						(*sCPU).totalSpectrum, 
						3 * (*sCPU).Nfreq * sizeof(double));
				}
				else if ((*sCPU).runType != runTypes::counter) {
					throw std::runtime_error(
						std::string("Attempted out-of-bounds save() to index ")
						.append(std::to_string(saveLoc)).append("\n"));
				}
			}
			break;
		case functionID("savePlasma"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			{
				int64_t saveLoc = (int64_t)parameters[0];
				int64_t plasmaLoc = (int)parameters[1];
				if (saveLoc < (*sCPU).Nsims 
					&& saveLoc != 0 
					&& (*sCPU).runType != runTypes::counter) 
					savePlasma(d, saveLoc, plasmaLoc);
				else {
					throw std::runtime_error(
						std::string("Attempted out-of-bounds savePlasma() to index ")
						.append(std::to_string(saveLoc)).append("\n"));
				}
			}
			break;
		case functionID("init"):
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

		case functionID("linear"):
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
				}
				(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();
				(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
				(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
				d.reset(sCPU);
				applyLinearPropagation(d, sCPU, (*sCPU).materialIndex, (*sCPU).crystalThickness);
				if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			}

			break;
		case functionID("fresnelLoss"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
			d.reset(sCPU);
			applyFresnelLoss(d, sCPU, d.deviceStruct,
				(int)parameters[4],
				(int)parameters[5]);
			break;
		case functionID("sphericalMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applySphericalMirror(d, sCPU, d.deviceStruct, parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("parabolicMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyParabolicMirror(d, sCPU, d.deviceStruct, parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("aperture"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyAperature(d, sCPU,
				parameters[0],
				parameters[1]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("farFieldAperture"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			(*sCPU).materialIndex = 0;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();
			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			applyAperatureFarField(d, sCPU,
				parameters[0],
				parameters[1],
				parameters[2],
				parameters[3]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("farFieldInverseAperture"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			(*sCPU).materialIndex = 0;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();
			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			applyInverseAperatureFarField(d, sCPU,
				parameters[0],
				parameters[1],
				parameters[2],
				parameters[3]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("energy"):
			{
			if ((*sCPU).runType == runTypes::counter) break;
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			int targetVar = (int)parameters[0];
			int spectrumType = (int)parameters[1];
			double energy = 0.0;
			for(int i = 0; i < (*sCPU).Nfreq; i++){
				energy += (*sCPU).totalSpectrum[i + (*sCPU).Nfreq * spectrumType];
			}
			energy *= (*sCPU).fStep;
			if(targetVar<100) vBlock[targetVar] = energy;
			}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("filter"):
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
		case functionID("lorentzian"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyLorenzian(
				d, 
				sCPU, 
				parameters[0], 
				parameters[1], 
				parameters[2], 
				parameters[3], 
				parameters[4]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("applyOptic"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			if((*sCPU).optics.size() < static_cast<int>(parameters[0])+1){
				throw std::runtime_error("applyOptic: optic is not loaded");
			}
			applyOptic(
				d, 
				*sCPU, 
				static_cast<int>(parameters[0]),
				true,
				true);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("applyOpticX"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			if((*sCPU).optics.size() < static_cast<int>(parameters[0])+1){
				throw std::runtime_error("applyOpticX: optic is not loaded");
			}
			applyOptic(
				d, 
				*sCPU, 
				static_cast<int>(parameters[0]),
				true,
				false);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("applyOpticY"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			if((*sCPU).optics.size() < static_cast<int>(parameters[0])+1){
				throw std::runtime_error("applyOpticY: optic is not loaded");
			}
			applyOptic(
				d, 
				*sCPU, 
				static_cast<int>(parameters[0]),
				false,
				true);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("loadOptic"):{
				loadedInputData file(cc.substr(cc.find('"')+1, cc.find('"', cc.find('"') + 1) - cc.find('"') - 1));
				if(file.hasData) (*sCPU).optics.push_back(file);
				else throw std::runtime_error("loadOptic failed, check the file path. If you are using the Flatpak version, load the optic using the load optic button.");
			}
			break;
		case functionID("addPulse"):
			if ((*sCPU).runType == runTypes::counter) break;
		{
			interpretParameters(cc, 21, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			d.deviceMemcpy(
				d.deviceStruct.gridETime1, 
				(*sCPU).ExtOut, 
				2 * d.deviceStruct.Ngrid * sizeof(double), 
				copyType::ToDevice);
			d.deviceMemcpy(
				d.deviceStruct.gridEFrequency1, 
				(*sCPU).EkwOut, 
				2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
				copyType::ToDevice);

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
			d.deviceMemcpy(
				(*sCPU).EkwOut, 
				d.deviceStruct.gridEFrequency1, 
				2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
				copyType::ToHost);
			d.deviceMemcpy(
				(*sCPU).ExtOut, 
				d.deviceStruct.gridETime1, 
				2 * (*sCPU).Ngrid * sizeof(double), 
				copyType::ToHost);

			getTotalSpectrum(d);
		}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		
		case functionID("for"): {
				interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
				int counter = (int)parameters[0];
				int targetVar = (int)parameters[1];
				std::string currentString = cc.substr(cc.find_first_of('{') + 1, std::string::npos);
				std::string forStartString = currentString;
				vBlock[targetVar] = 0.0;
				for (int i = 0; i < counter; i++) {

					while (currentString.length() > 0 && currentString.at(0) != '}') {
						if (currentString.at(0) == '<') {
							currentString =
								currentString.substr(currentString.find_first_of('>'), std::string::npos);
							if (currentString.length() > 0) currentString =
								currentString.substr(1, std::string::npos);
						}
						if (currentString.at(0) == '{') {
							currentString = currentString.substr(1, std::string::npos);
							while (currentString.find_first_of('{') != std::string::npos
								&& currentString.find_first_of('{') < currentString.find_first_of('}')) {
								currentString =
									currentString.substr(currentString.find_first_of('}'), std::string::npos);
								currentString = currentString.substr(1, std::string::npos);
							}
							currentString =
								currentString.substr(currentString.find_first_of('}'), std::string::npos);
							if (currentString.length() < 5) break;
							currentString = currentString.substr(1, std::string::npos);
						}
						error = interpretCommand(currentString, iBlock, vBlock, d, sCPU);
						currentString = 
							currentString.substr(findParenthesesClosure(currentString)+1, std::string::npos);
						if (error || (*sCPU).cancellationCalled) break;
					}
					++vBlock[targetVar];
					currentString = forStartString;
					if (error || (*sCPU).cancellationCalled) break;
				}
			}
			break;
		default: throw std::runtime_error(std::string("Could not interpret:\n").append(cc));
		}
		return error;
	}
	
	static int solveSequenceWithDevice(
		ActiveDevice& d, 
		simulationParameterSet* sCPU) {

		int error = 0;
		//if it starts with 0, it's an old sequence; quit
		if ((*sCPU).sequenceString[0] == '0') {
			return 15;
		}

		//main text interpreter
		double iBlock[100] = { 0.0 };

		for (int k = 1; k < 38; k++) {
			iBlock[k] = (*sCPU).getByNumberWithMultiplier(k);
		}

		double vBlock[100] = { 0.0 };
		std::string currentString((*sCPU).sequenceString);
		simulationParameterSet backupSet = *sCPU;
		//shortest command is for(), if there's only 4 characters left, it can only
		//be whitespace or other trailing symbols, and even then something is wrong,
		//since for should be followed by a loop... next shortest is init(), and that
		//certainly doesn't belong at the end. So anything under six characters=bail
		size_t minLength = 5;
		while (currentString.length() > minLength) {
			//skip curly braces (for loops should have been handled by interpretCommand() already)
			if (currentString.at(0) == '{') {
				currentString = currentString.substr(1, std::string::npos);
				while (currentString.find_first_of('{') != std::string::npos
					&& currentString.find_first_of('{') < currentString.find_first_of('}')) {
					currentString = 
						currentString.substr(currentString.find_first_of('}'), std::string::npos);
					currentString = currentString.substr(1, std::string::npos);
				}
				currentString = 
					currentString.substr(currentString.find_first_of('}'), std::string::npos);
				if (currentString.length() < minLength) break;
				currentString = currentString.substr(1, std::string::npos);
			}
			//skip angle brackets (comments)
			while (currentString.at(0) == '<') {
				currentString = 
					currentString.substr(currentString.find_first_of('>'), std::string::npos);
				if (currentString.length() < minLength) break;
				currentString = currentString.substr(1, std::string::npos);
			}
			if (error || (*sCPU).cancellationCalled) break;
			error = interpretCommand(currentString, iBlock, vBlock, d, sCPU);
			if (error || (*sCPU).cancellationCalled) break;
			currentString = 
				currentString.substr(findParenthesesClosure(currentString), std::string::npos);
			if (currentString.length() < minLength) break;

			currentString = currentString.substr(1, std::string::npos);
			backupSet.isFollowerInSequence = (*sCPU).isFollowerInSequence;
			backupSet.optics = (*sCPU).optics;
			*sCPU = backupSet;
		}
		*sCPU = backupSet;
		return error;
	}

	// helper function for fitting mode, runs the simulation and returns difference from the desired outcome.
	static double getResidual(
		const dlib::matrix<double, 0, 1>& x){

		double result = 0.0;
		for (int i = 0; i < (*fittingSet).Nfitting; ++i) {
			(*fittingSet).setByNumberWithMultiplier((int64_t)(*fittingSet).fittingArray[3 * i], x(i));
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
				result += 
					(*fittingSet).totalSpectrum[(*fittingSet).fittingMode 
					* (*fittingSet).Nfreq + (*fittingSet).fittingROIstart + i];
			}
			return result;
		}

		//mode 3 & 4: match total spectrum to reference given in ascii file
		double a;
		double maxSim = 0.0;
		double maxRef = 0.0;
		double* simSpec = &(*fittingSet).totalSpectrum[
			2 * (*fittingSet).Nfreq + (*fittingSet).fittingROIstart];
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
unsigned long solveNonlinearWaveEquationX(simulationParameterSet* sCPU) {
	ActiveDevice d(sCPU);
	if (d.memoryStatus) return 1;
	return solveNonlinearWaveEquationWithDevice(d, sCPU);
}

// Main function for running a sequence
unsigned long solveNonlinearWaveEquationSequenceX(simulationParameterSet* sCPU) {
	if ((*sCPU).batchIndex == 36 && (*sCPU).batchLoc1 != 0) return 0;
	ActiveDevice d(sCPU);
	if (d.memoryStatus) return 1;
	return solveSequenceWithDevice(d, sCPU);
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
		parameters(i) = (*sCPU).getByNumber((int64_t)round((*sCPU).fittingArray[3 * i]));
		lowerBounds(i) = (*sCPU).fittingArray[3 * i + 1];
		upperBounds(i) = (*sCPU).fittingArray[3 * i + 2];
	}

	dlib::function_evaluation result;
	if ((*sCPU).fittingMode != 3) {
		result = dlib::find_max_global(
			getResidual, 
			lowerBounds, 
			upperBounds, 
			dlib::max_function_calls((*sCPU).fittingMaxIterations));
	}
	else {
		result = dlib::find_min_global(
			getResidual, 
			lowerBounds, 
			upperBounds, 
			dlib::max_function_calls((*sCPU).fittingMaxIterations));
	}

	for (int i = 0; i < (*sCPU).Nfitting; ++i) {
		(*sCPU).setByNumberWithMultiplier((int64_t)round((*sCPU).fittingArray[3 * i]), result.x(i));
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



