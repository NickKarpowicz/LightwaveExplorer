#pragma once
//ugly macro prefix because CUDA doesn't like my nice constexprs
#ifdef __CUDACC__
#define hostOrDevice __host__ __device__
#else
#define hostOrDevice
#endif
#include <algorithm>
//variadic template to constexpr the product of a bunch of values
//in a way that keeps Xe Graphics happy (no doubles)
//convention: if there are multiple types as inputs, 
//return type is the type of the first argument
template<typename T>
hostOrDevice static constexpr T constProd(T x) {
    return x;
}
template<typename T, typename... Args>
hostOrDevice static constexpr T constProd(T x, Args... args) {
    return (T)(x * constProd(args...));
}

template <typename T, typename U>
hostOrDevice static constexpr inline T maxN(T a, U b) {
    return (((a) > (b)) ? (a) : (b));
}

template <typename T, typename U>
hostOrDevice static constexpr inline T minN(T a, U b) {
    return (((a) < (b)) ? (a) : (b));
}
//pi to stupid digits
template <typename T>
hostOrDevice static constexpr T vPi() {
    return (T)3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384L;
}
template <typename T>
hostOrDevice static constexpr T twoPi() {
    return (T)(2.0 * vPi<T>());
}
template <typename T>
hostOrDevice static constexpr T angleTolerance() {
    return (T)1e-12;
}
template <typename T>
hostOrDevice static constexpr T invSqrtPi() {
    return (T)(1.0 / sqrt(vPi<T>()));
}
template <typename T>
hostOrDevice static constexpr T deg2Rad() {
    return (T)(vPi<T>() / 180.0);
}
template <typename T>
hostOrDevice static constexpr T rad2Deg() {
    return (T)(180.0 / vPi<T>());
}
template <typename T>
hostOrDevice static constexpr T lightC() {
    return (T)2.99792458e8;
}
template <typename T>
hostOrDevice static constexpr T eps0() {
    return (T)8.8541878128e-12;
}
template <typename T>
hostOrDevice static constexpr T sixth() {
    return (T)((T)1.0 / 6.0);
}
template <typename T>
hostOrDevice static constexpr T third() {
    return (T)((T)1.0 / 3.0);
}
template <typename T>
hostOrDevice static constexpr T kLorentzian() {
    return (T)3182.607353999257;
}

//give the Dawson function value at x
//used for Gaussian-based "Sellmeier" equation, as it is the Hilbert transform of the Gaussian function (susceptibility obeys Kramers-Kronig relations)
//based on Rybicki G.B., Computers in Physics, 3,85-87 (1989)
//this is the simplest implementation of the formula he provides, he also suggests speed improvements in case
//evaluation of this becomes a bottleneck
//macro for the two versions because I don't want to lose precision in the double version
//but Intel Xe graphics need the constants as floats, and CUDA doesn't work with the
//[[maybe_unused]] attribute... and I don't like compiler warnings.
[[maybe_unused]] hostOrDevice static float deviceDawson(const float x) {
	//parameters determining accuracy (higher n, smaller h -> more accurate but slower)
	int n = 15;
	float h = 0.3f;

	//series expansion for small x
	if (fabs(x) < 0.2f) {
		float x2 = x * x;
		float x4 = x2 * x2;
		return x * (1.0f - 2.0f * x2 / 3.0f + 4.0f * x4 / 15.0f - 8.0f * x2 * x4 / 105.0f + (16.0f / 945.0f) * x4 * x4 - (32.0f / 10395.0f) * x4 * x4 * x2);
	}

	int n0 = 2 * int(round(0.5f * x / h));
	float x0 = h * n0;
	float xp = x - x0;
	float d = 0.0f;
	for (int i = -n; i < n; i++) {
		if (i % 2 != 0) {
			d += expf(-(xp - i * h) * (xp - i * h)) / (i + n0);
		}
	}
	return invSqrtPi<float>() * d;
}

[[maybe_unused]] hostOrDevice static double deviceDawson(const double x) {
	//parameters determining accuracy (higher n, smaller h -> more accurate but slower)
	int n = 15;
	double h = 0.3;

	//series expansion for small x
	if (abs(x) < 0.2) {
		double x2 = x * x;
		double x4 = x2 * x2;
		return x * (1.0 - 2.0 * x2 / 3.0 + 4.0 * x4 / 15.0 - 8.0 * x2 * x4 / 105.0 + (16.0 / 945) * x4 * x4 - (32.0 / 10395) * x4 * x4 * x2);
	}

	int n0 = 2 * int(round(0.5 * x / h));
	double x0 = h * n0;
	double xp = x - x0;
	double d = 0.0;
	for (int i = -n; i < n; i++) {
		if (i % 2 != 0) {
			d += exp(-(xp - i * h) * (xp - i * h)) / (i + n0);
		}
	}
	return invSqrtPi<double>() * d;
}


//this is a job for std::erase, but when running the code on the cluster, everything
//is done with nvcc, with only an old version of cmake available. This means I can't
//use c++20 features there. So if it's compiled with c++17, use the more complicated
//function, otherwise just inline to std::erase.
#if __cplusplus<=201703L
inline void removeCharacterFromString(std::string& s, char removedChar) {
    s.erase(std::remove(s.begin(), s.end(), removedChar), s.end());
}
#else
inline void removeCharacterFromString(std::string& s, char removedChar) {
	std::erase(s,removedChar);
}
#endif