#include "LightwaveExplorerUtilities.h"

unsigned long   solveNonlinearWaveEquationSequence(void* lpParam);
unsigned long	solveNonlinearWaveEquation(void* lpParam);

#ifdef __CUDACC__
//In tests this mattered, since Thrust does math between complex and double up casting the double to a complex.
__device__ thrust::complex<double> operator/(double a, thrust::complex<double> b) {
    double divByDenominator = a / (b.real() * b.real() + b.imag() * b.imag());
    return thrust::complex<double>(b.real() * divByDenominator, -b.imag() * divByDenominator);
}
__device__ thrust::complex<double> operator/(thrust::complex<double> a, double b) { return thrust::complex<double>(a.real() / b, a.imag() / b); }

__device__ thrust::complex<double> operator*(double b, thrust::complex<double> a) { return thrust::complex<double>(a.real() * b, a.imag() * b); }
__device__ thrust::complex<double> operator*(thrust::complex<double> a, double b) { return thrust::complex<double>(a.real() * b, a.imag() * b); }

__device__ thrust::complex<double> operator+(double a, thrust::complex<double> b) { return thrust::complex<double>(b.real() + a, b.imag()); }
__device__ thrust::complex<double> operator+(thrust::complex<double> a, double b) { return thrust::complex<double>(a.real() + b, a.imag()); }

__device__ thrust::complex<double> operator-(double a, thrust::complex<double> b) { return thrust::complex<double>(a - b.real(), -b.imag()); }
__device__ thrust::complex<double> operator-(thrust::complex<double> a, double b) { return thrust::complex<double>(a.real() - b, a.imag()); }

#endif



