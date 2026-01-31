#pragma once
//ugly macro prefix because CUDA doesn't like my nice constexprs
#ifdef __CUDACC__
#define hostOrDevice __host__ __device__
#else
#define hostOrDevice
#endif
#include <algorithm>
#include <vector>
#include <ranges>
#include <string_view>
#include <string>
//variadic template to constexpr the product of a bunch of values
//in a way that keeps Xe Graphics happy (no doubles)
//convention: if there are multiple types as inputs,
//return type is the type of the first argument
template<typename T>
hostOrDevice static constexpr T constProd(const T x) {
    return x;
}
template<typename T, typename... Args>
hostOrDevice static constexpr T constProd(const T x, const Args... args) {
    return static_cast<T>(x * constProd(args...));
}

template <typename T, typename U>
hostOrDevice static constexpr inline T maxN(const T a, const U b) {
    return (((a) > (b)) ? (a) : (b));
}

template <typename T, typename U>
hostOrDevice static constexpr inline T minN(const T a, const U b) {
    return (((a) < (b)) ? (a) : (b));
}
//pi to stupid digits
template <typename T>
hostOrDevice static consteval T vPi() {
    return static_cast<T>(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384L);
}
template <typename T>
hostOrDevice static constexpr T twoPi() {
    return (T)(2.0 * vPi<T>());
}
template <typename T>
hostOrDevice static consteval T sqrtTwo() {
    return (T)(1.4142135623730951);
}
template <typename T>
hostOrDevice static constexpr T invSqrtPi() {
    return (T)(0.56418958354775627928034964497783221304416656494141L);
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
hostOrDevice static consteval T lightC() {
    return (T)2.99792458e8;
}
template <typename T>
hostOrDevice static consteval T eps0() {
    return (T)8.8541878128e-12;
}
template <typename T>
hostOrDevice static consteval T inverseEps0() {
    return (T)(1.0/8.8541878128e-12);
}
template <typename T>
hostOrDevice static consteval T mu0() {
    return (T)1.2566370614e-6;
}
template <typename T>
hostOrDevice static constexpr T inverseMu0() {
    return (T)(1.0/1.2566370614e-6);
}
template <typename T>
hostOrDevice static consteval T Zo() {
    return (T)376.730313668;
}
template <typename T>
hostOrDevice static consteval T inverseZo() {
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
hostOrDevice static consteval T kLorentzian() {
    return (T)3182.607353999257;
}
template <typename T>
hostOrDevice static consteval T cZero() {
    return T{};
}
template <typename T>
hostOrDevice static constexpr T cOne() {
    return T(1.0, 0.0);
}
template <typename T>
hostOrDevice static consteval T elCharge() {
    return (T)1.602176634e-19;
}
template <typename T>
hostOrDevice static consteval T elMass() {
    return (T)9.1093837015e-31;
}
template <typename T>
hostOrDevice static consteval T planckConstant() {
    return (T)6.62607015e-34;
}
template <typename T>
hostOrDevice static consteval T eVtoHz() {
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

inline void removeCharacterFromString(std::string& s, const char removedChar) {
	std::erase(s, removedChar);
}


template <typename T>
std::vector<std::vector<double>> parse_string_to_vecs(
    const std::string& txt,
    char line_delimiter=';',
    char number_delimiter=' ')
{


    auto outer = txt | std::ranges::views::split(line_delimiter)
                 | std::ranges::views::transform([number_delimiter](auto&& row_range) {

                       std::string_view row_sv{
                           &*row_range.begin(),
                           static_cast<std::size_t>(row_range.end() - row_range.begin())
                       };
                       auto inner_view = row_sv
                                       | std::ranges::views::split(number_delimiter)
                                       | std::ranges::views::filter([](auto&& sub) { return !sub.empty(); })
                                       | std::ranges::views::transform([](auto&& token_range) {
                                           auto stod_or_zero = [](std::string_view sv){
                                               try{
                                                   return std::stod(std::string(sv));
                                               } catch (...) {
                                                   std::cout << "bad stod encountered!\n";
                                                   return 0.0;
                                               }
                                           };
                                             std::string_view token_sv{
                                                 &*token_range.begin(),
                                                 static_cast<std::size_t>(token_range.end() -
                                                                          token_range.begin())
                                             };
                                             return stod_or_zero(token_sv);
                                         });

                       return std::vector<T>{ std::ranges::begin(inner_view),
                                                    std::ranges::end(inner_view) };
                   });

    auto result = std::vector<std::vector<T>>{
        std::ranges::begin(outer), std::ranges::end(outer)
    };
    if(result.size() == 0){
        throw std::runtime_error("Unparsable or empty string.");
    }

    result.erase(std::remove_if(result.begin(), result.end(), [](auto&& c){ return c.empty(); }), result.end());

    return result;
}
