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
hostOrDevice static constexpr T sqrtTwo() {
    return (T)(1.4142135623730951);
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
hostOrDevice static constexpr T inverseEps0() {
    return (T)(1.0/8.8541878128e-12);
}
template <typename T>
hostOrDevice static constexpr T mu0() {
    return (T)1.2566370614e-6;
}
template <typename T>
hostOrDevice static constexpr T inverseMu0() {
    return (T)(1.0/1.2566370614e-6);
}
template <typename T>
hostOrDevice static constexpr T Zo() {
    return (T)376.730313668;
}
template <typename T>
hostOrDevice static constexpr T inverseZo() {
    return 0.0026544187287534486;
}
template <typename T>
hostOrDevice static constexpr T sixth() {
    return (T)((T)1.0 / 6.0);
}
template <typename T>
hostOrDevice static constexpr T sixtieth() {
    return (T)((T)1.0 / 60.0);
}
template <typename T>
hostOrDevice static constexpr T third() {
    return (T)((T)1.0 / 3.0);
}
template <typename T>
hostOrDevice static constexpr T kLorentzian() {
    return (T)3182.607353999257;
}
template <typename T>
hostOrDevice static constexpr T cZero() {
    return T{};
}
template <typename T>
hostOrDevice static constexpr T cOne() {
    return T(1.0, 0.0);
}
template <typename T>
hostOrDevice static constexpr T elCharge() {
    return (T)1.602176634e-19;
}
template <typename T>
hostOrDevice static constexpr T elMass() {
    return (T)9.1093837015e-31;
}
template <typename T>
hostOrDevice static constexpr T planckConstant() {
    return (T)6.62607015e-34;
}
template <typename T>
hostOrDevice static constexpr T eVtoHz() {
    return elCharge<T>() / planckConstant<T>();
}
template<typename T>
hostOrDevice static constexpr T firstDerivativeStencil(const int i){
    switch(i){
        case -3: return static_cast<T>(-1./60.);
        case -2: return static_cast<T>(3./20.);
        case -1: return static_cast<T>(-3./4.);
        case 1: return static_cast<T>(3./4.);
        case 2: return static_cast<T>(-3./20.);
        case 3: return static_cast<T>(1./60.);
        default: return static_cast<T>(0.0);
    }
    return static_cast<T>(0.0);
}

//this is a job for std::erase, but when running the code on the cluster, everything
//is done with nvcc, with only an old version of cmake available. This means I can't
//use c++20 features there. So if it's compiled with c++17, use the more complicated
//function, otherwise just inline to std::erase.
#if __cplusplus<=201703L
inline void removeCharacterFromString(std::string& s, const char removedChar) {
    s.erase(std::remove(s.begin(), s.end(), removedChar), s.end());
}
#else
inline void removeCharacterFromString(std::string& s, const char removedChar) {
	std::erase(s, removedChar);
}
#endif