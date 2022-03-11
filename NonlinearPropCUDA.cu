#include "NonlinearPropCUDA.cuh"
#include "framework.h"
#include<complex>
#include<cstdlib>
#include<math.h>
#include "MPQ_Nonlinear_Propagation.h"
#include <cuComplex.h>
#include <cufft.h>
#include "qr_solve.hpp"
#include "MPQ_Nonlinear_Propagation.h"
//#include <complex.h>

#define THREADS_PER_BLOCK 64
#define MAX_LOADSTRING 1024
//overload the math operators for cuda complex numbers so this code fits inside the observable universe
__device__ cuDoubleComplex operator*(cuDoubleComplex a, cuDoubleComplex b) { return cuCmul(a, b); }
__device__ cuDoubleComplex operator+(cuDoubleComplex a, cuDoubleComplex b) { return cuCadd(a, b); }
__device__ cuDoubleComplex operator+(double a, cuDoubleComplex b) { return cuCadd(make_cuDoubleComplex(a, 0.0), b); }
__device__ cuDoubleComplex operator+(cuDoubleComplex a, double b) { return cuCadd(a, make_cuDoubleComplex(b, 0.0)); }
__device__ cuDoubleComplex operator-(cuDoubleComplex a, cuDoubleComplex b) { return cuCsub(a, b); }
__device__ cuDoubleComplex operator-(double a, cuDoubleComplex b) { return cuCsub(make_cuDoubleComplex(a, 0.0), b); }
__device__ cuDoubleComplex operator/(cuDoubleComplex a, double b) { return cuCdiv(a, make_cuDoubleComplex(b, 0.0)); }
__device__ cuDoubleComplex operator/(double b, cuDoubleComplex a) { return cuCdiv(make_cuDoubleComplex(b, 0.0), a); }
__device__ cuDoubleComplex operator*(cuDoubleComplex a, double b) { return cuCmul(a, make_cuDoubleComplex(b, 0.0)); }
__device__ cuDoubleComplex operator*(double b, cuDoubleComplex a) { return cuCmul(a, make_cuDoubleComplex(b, 0.0)); }


//complex exponential function for CUDA
__device__ __forceinline__ cuDoubleComplex cuCexpd(cuDoubleComplex z)
{
    cuDoubleComplex res;
    double t = exp(z.x);
    res.y = sin(z.y);
    res.x = cos(z.y);
    res.x *= t;
    res.y *= t;
    return res;
}

//copy and paste from
// https://forums.developer.nvidia.com/t/additional-cucomplex-functions-cucnorm-cucsqrt-cucexp-and-some-complex-double-functions/36892 
__device__ cuDoubleComplex cuCsqrt(cuDoubleComplex x)
{
    double radius = cuCabs(x);
    double cosA = x.x / radius;
    cuDoubleComplex out;
    out.x = sqrt(radius * (cosA + 1.0) / 2.0);
    out.y = sqrt(radius * (1.0 - cosA) / 2.0);
    // signbit should be false if x.y is negative
    if (signbit(x.y))
        out.y *= -1.0;

    return out;
}


__device__ cuDoubleComplex sellmeierCuda(cuDoubleComplex* ne, cuDoubleComplex* no, double* a, double f, double theta, double phi, int type, int eqn) {
    if (f == 0) return make_cuDoubleComplex(1.0,0.0); //exit immediately for f=0
    
    double c = 2.99792458e8; //speed of light
    double l = 1e6 * c / f; //wavelength in microns
    double ls = l * l;
    cuDoubleComplex ii = make_cuDoubleComplex(0.0, 1.0);
    double pi = 3.14159265358979323846264338327950288;
    double omega = 2 * pi * abs(f);
    double kL = 3183.9; //(e * e / (e_o *m_e)
    cuDoubleComplex one = make_cuDoubleComplex(1.0, 0);
    cuDoubleComplex na = one;
    cuDoubleComplex nb = one;
    //option 0: isotropic
    if (type == 0) {
        ne[0] = make_cuDoubleComplex(a[0]
            + (a[1] + a[2] * ls) / (ls + a[3]) + (a[4] + a[5] * ls) / (ls + a[6])
            + (a[7] + a[8] * ls) / (ls + a[9]) + (a[10] + a[11] * ls) / (ls + a[12])
            + a[13] * ls + a[14] * ls * ls + a[15] * ls * ls * ls, 0.0);
        if (cuCreal(ne[0]) < 1) {
            ne[0] = one;
        }
        ne[0] = ne[0] + kL * a[16] / (a[17] - omega * omega - ii * a[18] * omega)
            + kL * a[19] / (a[20] - omega * omega - ii * a[21] * omega);
        ne[0] = cuConj(cuCsqrt(ne[0]));
        if (isnan(cuCreal(ne[0]))) {
            ne[0] = one;
        }

        no[0] = ne[0];
        return ne[0];
    }
    //option 1: uniaxial
    else if (type == 1) {
        
        na = cuCsqrt(a[0]
            + (a[1] + a[2] * ls) / (ls + a[3]) + (a[4] + a[5] * ls) / (ls + a[6]) + (a[7]
                + a[8] * ls) / (ls + a[9]) + (a[10] + a[11] * ls) / (ls + a[12])
            + a[13] * ls + a[14] * ls * ls + a[15] * ls * ls * ls
            + kL * a[16] / (a[17] - omega * omega + ii * a[18] * omega)
            + kL * a[19] / (a[20] - omega * omega + ii * a[21] * omega));
        
        a = &a[22];
        nb = cuCsqrt(a[0]
            + (a[1] + a[2] * ls) / (ls + a[3]) + (a[4] + a[5] * ls) / (ls + a[6]) + (a[7]
                + a[8] * ls) / (ls + a[9]) + (a[10] + a[11] * ls) / (ls + a[12])
            + a[13] * ls + a[14] * ls * ls + a[15] * ls * ls * ls
            + kL * a[16] / (a[17] - omega * omega - ii * a[18] * omega)
            + kL * a[19] / (a[20] - omega * omega - ii * a[21] * omega));
        no[0] = na;
        ne[0] = 1.0 / cuCsqrt(cos(theta) * cos(theta) / (na * na) + sin(theta) * sin(theta) / (nb * nb));
        return ne[0];
    }
    else {
        //later, implement biaxial crystals, for now just return 1;
        return one;
    }
}
__global__ void radialLaplacianKernel(struct cudaLoop s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j = i / s.Ntime; //spatial coordinate
    long long k = i - j * s.Ntime; //temporal coordinate
    double rho = s.dx * j - (s.dx / 2) * s.Nspace;

    //zero at edges of grid and at origin
    if (abs(rho<1e-15) || j<3 || j>(s.Nspace-4)) {
        s.gridRadialLaplacian1[i] = make_cuDoubleComplex(0, 0);
        s.gridRadialLaplacian2[i] = make_cuDoubleComplex(0, 0);
    }
    else {
        s.gridRadialLaplacian1[i] = (s.firstDerivativeOperation[0] * s.gridETime[i - 3 * s.Ntime]
            + s.firstDerivativeOperation[1] * s.gridETime[i - 2 * s.Ntime]
            + s.firstDerivativeOperation[2] * s.gridETime[i - s.Ntime]
            + s.firstDerivativeOperation[3] * s.gridETime[i + s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime[i + 2 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime[i + 3 * s.Ntime])/rho;
        s.gridRadialLaplacian2[i] = (s.firstDerivativeOperation[0] * s.gridETime2[i - 3 * s.Ntime]
            + s.firstDerivativeOperation[1] * s.gridETime2[i - 2 * s.Ntime]
            + s.firstDerivativeOperation[2] * s.gridETime2[i - s.Ntime]
            + s.firstDerivativeOperation[3] * s.gridETime2[i + s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime2[i + 2 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime2[i + 3 * s.Ntime]) / rho;
    }
}

__global__ void prepareCartesianGridsKernel(double* theta, double* sellmeierCoefficients, struct cudaLoop s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j, k;
    long long Ntime = s.Ntime;
    long long Nspace = s.Nspace;
    int axesNumber = s.axesNumber;
    int sellmeierType = s.sellmeierType;
    double c = 2.99792458e8; //speed of light
    double pi = 3.14159265358979323846264338327950288; // pi to unneccessary precision
    cuDoubleComplex cuZero = make_cuDoubleComplex(0, 0);
    j = i / Ntime; //spatial coordinate
    k = i - j*Ntime; //temporal coordinate
    cuDoubleComplex ii = make_cuDoubleComplex(0, 1);
    double crystalTheta = sellmeierCoefficients[66];
    double crystalPhi = sellmeierCoefficients[67];
    double kStep = sellmeierCoefficients[70];
    double fStep = sellmeierCoefficients[71];
    double tol = sellmeierCoefficients[72];
    double dTheta = 0.1;
    double err, errPlus, errMinus;
    
    cuDoubleComplex ne, no, n0;
    double nePlus, neMinus;



    //frequency being resolved by current thread
    double f = k * fStep;
	if (k >= Ntime / 2) {
		f -= fStep * Ntime;
	}
	f *= -1;

    //transverse wavevector being resolved
	double dk = j * kStep - (j >= (Nspace / 2)) * (kStep * Nspace); //frequency grid in transverse direction


    //Find walkoff angle, starting from zero
    theta[i] = 0;
    double rhs = 2.99792458e8 * dk / (2 * 3.14159265358979323846264338327950288 * f);
    sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta + theta[i], crystalPhi, axesNumber, sellmeierType);    
    nePlus = cuCreal(ne);
    err = abs(nePlus * sin(theta[i]) - rhs);

    int iters = 0;
    errPlus = 2;
    errMinus = 2;
    while (err > tol && iters < 2048) {
        iters++;

        sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta + theta[i] + dTheta, crystalPhi, axesNumber, sellmeierType);
        nePlus = cuCreal(ne);
        errPlus = abs(nePlus * sin(theta[i] + dTheta) - rhs);

        sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta + theta[i] - dTheta, crystalPhi, axesNumber, sellmeierType);
        neMinus = cuCreal(ne);
        errMinus = abs(neMinus * sin(theta[i] - dTheta) - rhs);

        //Basic hill climbing algorithm
        //calculate the error at theta +/- dTheta
        // if theta + dTheta has lowest error, theta = theta+dTheta, err = errPlus
        // if theta - dTheta has lowest error, theta = theta-dTheta, err = errMinus
        // if theta has lowest error, step size is too large, dTheta /= 2;
        if (errPlus < err && errPlus < errMinus) {
            theta[i] += dTheta;
            err = errPlus;
        }
        else if (errMinus < err) {
            theta[i] -= dTheta;
            err = errMinus;
        }
        else {
            dTheta *= 0.5;
        }
        
    }


    //walkoff angle has been found, generate the rest of the grids
    f = k * fStep;
    if (k >= Ntime / 2) {
        f -= fStep * Ntime;
    }
    f *= -1;

    sellmeierCuda(&n0, &no, sellmeierCoefficients, abs(s.f0), crystalTheta, crystalPhi, axesNumber, sellmeierType);
    sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta + theta[i], crystalPhi, axesNumber, sellmeierType);
    if (isnan(cuCreal(ne)) || isnan(cuCreal(no))) {
        ne = make_cuDoubleComplex(1,0);
        no = make_cuDoubleComplex(1, 0);
    }
    s.ne[i] = ne;
    s.no[i] = no;

    cuDoubleComplex k0 = make_cuDoubleComplex(2 * pi * cuCreal(n0) * f / c,0);
    cuDoubleComplex ke = 2 * pi * ne * f / c;
    cuDoubleComplex ko = 2 * pi * no * f / c;

    if (cuCreal(ke) < 0 && cuCreal(ko) < 0) {
        s.gridPropagationFactor[i] = ii * (ke - k0 + dk * dk / (2. * cuCreal(ke))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor[i]))) {
            s.gridPropagationFactor[i] = cuZero;
        }

        s.gridPropagationFactor2[i] = ii * (ko - k0 + dk * dk / (2. * cuCreal(ko))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor2[i]))) {
            s.gridPropagationFactor2[i] = cuZero;
        }

        int posf = (int)(f < -20e12);
        s.gridPolarizationFactor[i] = ii * (posf * 2 * pi * f) / (2. * cuCreal(ne) * c) * s.h;
        s.gridPolarizationFactor2[i] = ii * (posf * 2 * pi * f) / (2. * cuCreal(no) * c) * s.h;
    }
    
    else {
        s.gridPropagationFactor[i] = cuZero;
        s.gridPropagationFactor2[i] = cuZero;
        s.gridPolarizationFactor[i] = cuZero;
        s.gridPolarizationFactor2[i] = cuZero;
    }
    
    
}

__global__ void prepareCylindricGridsKernel(double* sellmeierCoefficients, struct cudaLoop s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j, k;
    long long Ntime = s.Ntime;
    long long Nspace = s.Nspace;
    int axesNumber = s.axesNumber;
    int sellmeierType = s.sellmeierType;
    double c = 2.99792458e8; //speed of light
    double pi = 3.14159265358979323846264338327950288; // pi to unneccessary precision
    cuDoubleComplex cuZero = make_cuDoubleComplex(0, 0);
    j = i / Ntime; //spatial coordinate
    k = i - j * Ntime; //temporal coordinate
    cuDoubleComplex ii = make_cuDoubleComplex(0, 1);
    double crystalTheta = sellmeierCoefficients[66];
    double crystalPhi = sellmeierCoefficients[67];
    double kStep = sellmeierCoefficients[70];
    double fStep = sellmeierCoefficients[71];
    double tol = sellmeierCoefficients[72];
    double dTheta = 0.1;
    double err, errPlus, errMinus;

    cuDoubleComplex ne, no, n0;

    //frequency being resolved by current thread
    double f = k * fStep;
    if (k >= Ntime / 2) {
        f -= fStep * Ntime;
    }
    f *= -1;

    //transverse wavevector being resolved
    double dk = j * kStep - (j >= (Nspace / 2)) * (kStep * Nspace); //frequency grid in transverse direction



    sellmeierCuda(&n0, &no, sellmeierCoefficients, abs(s.f0), crystalTheta, crystalPhi, axesNumber, sellmeierType);
    sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), crystalTheta, crystalPhi, axesNumber, sellmeierType);
    if (isnan(cuCreal(ne)) || isnan(cuCreal(no))) {
        ne = make_cuDoubleComplex(1, 0);
        no = make_cuDoubleComplex(1, 0);
    }
    s.ne[i] = ne;
    s.no[i] = no;

    cuDoubleComplex k0 = make_cuDoubleComplex(2 * pi * cuCreal(n0) * f / c, 0);
    cuDoubleComplex ke = 2 * pi * ne * f / c;
    cuDoubleComplex ko = 2 * pi * no * f / c;

    if (cuCreal(ke) < 0 && cuCreal(ko) < 0) {
        s.gridPropagationFactor[i] = ii * (ke - k0 + dk * dk / (2. * cuCreal(ke))) * s.h;
        s.gridPropagationFactorRho1[i] = ii * (1 / (2. * cuCreal(ke))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor[i]))) {
            s.gridPropagationFactor[i] = cuZero;
            s.gridPropagationFactorRho1[i] = cuZero;
        }

        s.gridPropagationFactor2[i] = ii * (ko - k0 + dk * dk / (2. * cuCreal(ko))) * s.h;
        s.gridPropagationFactorRho2[i] = ii * (1 / (2. * cuCreal(ko))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor2[i]))) {
            s.gridPropagationFactor2[i] = cuZero;
            s.gridPropagationFactorRho2[i] = cuZero;
        }

        int posf = (int)(f < -20e12);
        s.gridPolarizationFactor[i] = ii * (posf * 2 * pi * f) / (2. * cuCreal(ne) * c) * s.h;
        s.gridPolarizationFactor2[i] = ii * (posf * 2 * pi * f) / (2. * cuCreal(no) * c) * s.h;
    }

    else {
        s.gridPropagationFactor[i] = cuZero;
        s.gridPropagationFactor2[i] = cuZero;
        s.gridPolarizationFactor[i] = cuZero;
        s.gridPolarizationFactor2[i] = cuZero;
        s.gridPropagationFactor[i] = cuZero;
        s.gridPropagationFactorRho2[i] = cuZero;
    }


}
//replaces E with its complex conjugate
__global__ void conjugateKernel(cuDoubleComplex* E) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    E[i] = cuConj(E[i]);
}

//replaces NaN values with 0
__global__ void fixnanKernel(cuDoubleComplex* E) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    if (isnan(cuCreal(E[i])) || isnan(cuCimag(E[i]))) {
        E[i] = make_cuDoubleComplex(0., 0.);
    }
}
__global__ void nonlinearpolarizationKernel(struct cudaLoop s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    double Ex = cuCreal(s.gridETime[i]) / s.propagationInts[0];
    double Ey = cuCreal(s.gridETime2[i]) / s.propagationInts[0];
    s.gridPolarizationTime[i] = 0.;
    s.gridPolarizationTime2[i] = 0.;

    //The d2eff tensor has the form
    // | d_xxx d_xyx d_yyx |
    // | d_xxy d_xyy d_yyy |
    if (s.nonlinearSwitches[0] == 1) {
        s.gridPolarizationTime[i] += s.chi2Tensor[0] * Ex * Ex + s.chi2Tensor[2] * Ex * Ey + s.chi2Tensor[4] * Ey * Ey;
        s.gridPolarizationTime2[i] += s.chi2Tensor[1] * Ex * Ex + s.chi2Tensor[3] * Ex * Ey + s.chi2Tensor[5] * Ey * Ey;
    }
    
    //to be implemented: full chi3 matrix on s.nonlinearSwitches[1]==1

    //using only one value of chi3, under assumption of centrosymmetry
    if (s.nonlinearSwitches[1] == 2) {
        s.gridPolarizationTime[i] += s.chi3Tensor[0] * (Ex * Ex * Ex + Ey * Ey * Ex / 3.);
        s.gridPolarizationTime2[i] += s.chi3Tensor[0] * (Ey * Ey * Ey + Ex * Ex * Ey / 3.);
    }

    //Nonlinear absorption
    if (s.nonlinearSwitches[2] == 1) {
        
        double Exi = cuCimag(s.gridETime[i]) / s.propagationInts[0];
        double Eyi = cuCimag(s.gridETime2[i]) / s.propagationInts[0];
        double fieldAmp2 = Exi * Exi + Eyi * Eyi;
        int j;
        for (j = 0; j < s.nonlinearSwitches[3]; j++) {
            Exi *= fieldAmp2;
            Eyi *= fieldAmp2;
        }
        s.gridPolarizationTime[i] += s.absorptionParameters[1] * Exi;
        s.gridPolarizationTime2[i] += s.absorptionParameters[1] * Eyi;
    }
}

//Plasma response with time-dependent carrier density
//This polarization needs a different factor in the nonlinear wave equation
//to account for the integration
//TO BE IMPLEMENTED!
__global__ void plasmaCurrentKernel(struct cudaLoop s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    int j,k,l;
    double N = 0;
    double integralx = 0;
    double integraly = 0;
    double t, w, Esquared, Ex, Ey;
    for (j = 0; j < s.Ntime; j++) {
        Ex = cuCreal(s.gridETime[l]);
        Ey = cuCreal(s.gridETime2[l]);
        l = j + i * s.Nspace;
        t = j * s.dt;
        Esquared = Ex * Ex + Ey * Ey;
        w = Esquared;
        for (k = 0; k < s.nonlinearSwitches[3]; k++) {
            w *= Esquared;
        }
        N += w;
        integralx += exp(-s.plasmaParameters[1] * t) * N * Ex;
        integraly += exp(-s.plasmaParameters[1] * t) * N * Ey;
        s.plasmaCurrent1[l] = s.plasmaParameters[0] * exp(-s.plasmaParameters[1] * t) * integralx;
        s.plasmaCurrent2[l] = s.plasmaParameters[0] * exp(-s.plasmaParameters[1] * t) * integraly;
    }
}


//Main kernel for RK4 propagation of the field
__global__ void rkKernel(struct cudaLoop s, int stepNumber) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j = i / s.Ntime; //spatial coordinate
    long long h = i - j * s.Ntime; //temporal coordinate


    if (s.isCylindric) {
        s.gridRadialLaplacian1[i] = s.gridPropagationFactorRho1[i] * s.k1[i];
        s.gridRadialLaplacian2[i] = s.gridPropagationFactorRho2[i] * s.k2[i];
    }

    //polarization is stored in a reduced format by cuFFT because the FFT is from real to complex, meaning if the output grid
    //were to be N_time x N_space, half of the points would be redundant. The extra steps below are to determine where in the grid the 
    //current point sits. Essentially, if in the negative frequency quadrants, reverse the frequency and take complex conjugate of the 
    //value
    if (h > s.propagationInts[3]) {
        h = s.Ntime - h;
        j = s.Nspace - j;
        h += j * s.propagationInts[3];
        s.k1[i] = s.gridPropagationFactor[i] * s.gridETemp[i] +s.gridPolarizationFactor[i] * cuConj(s.gridPolarizationFrequency[i]);
        s.k2[i] = s.gridPropagationFactor2[i] * s.gridETemp2[i] +s.gridPolarizationFactor2[i] * cuConj(s.gridPolarizationFrequency2[i]);
    }
    else {
        h += j * s.propagationInts[3];
        s.k1[i] = s.gridPropagationFactor[i] * s.gridETemp[i] +s.gridPolarizationFactor[i] * s.gridPolarizationFrequency[h];
        s.k2[i] = s.gridPropagationFactor2[i] * s.gridETemp2[i] +s.gridPolarizationFactor2[i] * s.gridPolarizationFrequency2[h];
    }
    if (s.isCylindric) {
        s.k1[i] = s.k1[i] + s.gridRadialLaplacian1[i];
        s.k2[i] = s.k2[i] + s.gridRadialLaplacian2[i];
    }

    //in the first substep, first construct the next intermediate field value
    //which will be used in the next substep. 
    if (stepNumber == 0) {
        s.gridETemp[i] = s.gridEFrequency[i] + 0.5 * s.k1[i];
        s.gridETemp2[i] = s.gridEFrequency2[i] + 0.5 * s.k2[i];
       
        s.gridEFrequencyNext1[i] = s.k1[i] / 6 + s.gridEFrequency[i];
        s.gridEFrequencyNext2[i] = s.k2[i] / 6 + s.gridEFrequency2[i];
    }

    //in the next substep, again construct the next intermediate field and add k/3 to solution
    else if (stepNumber == 1) {
        s.gridETemp[i] = s.gridEFrequency[i] + 0.5 * s.k1[i];
        s.gridETemp2[i] = s.gridEFrequency2[i] + 0.5 * s.k2[i];

        s.gridEFrequencyNext1[i] = s.gridEFrequencyNext1[i] + s.k1[i] / 3;
        s.gridEFrequencyNext2[i] = s.gridEFrequencyNext2[i] + s.k2[i] / 3;

    }

    //same action as previous substep, except the weight of k in the intermediate solution is 1 instead of 0.5
    else if (stepNumber == 2) {
        s.gridETemp[i] = s.gridEFrequency[i] + s.k1[i];
        s.gridETemp2[i] = s.gridEFrequency2[i] + s.k2[i];
        s.gridEFrequencyNext1[i] = s.gridEFrequencyNext1[i] + s.k1[i] / 3;
        s.gridEFrequencyNext2[i] = s.gridEFrequencyNext2[i] + s.k2[i] / 3;
    }

    //last substep. Solution is now complete and may be copied directly into the field arrays
    else {
        s.gridEFrequency[i] = s.gridEFrequencyNext1[i] + s.k1[i] / 6;
        s.gridEFrequency2[i] = s.gridEFrequencyNext2[i] + s.k2[i] / 6;
        s.gridETemp[i] = s.gridEFrequency[i];
        s.gridETemp2[i] = s.gridEFrequency2[i];
    }

}


//Take absolute value of complex array
__global__ void absKernel(double* absOut, cuDoubleComplex* complexIn) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    absOut[i] = cuCabs(complexIn[i]);
}

//Apply fft normalization
__global__ void fftNormalizeKernel(cuDoubleComplex* A, long long* fftSize) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    A[i] = A[i] / fftSize[0];
}

DWORD WINAPI propagationLoop(LPVOID lpParam) {

    //the struct s contains most of the simulation variables and pointers
    struct cudaLoop s;
    struct propthread* sCPU = (struct propthread*)lpParam;


    //initialize and take values from the struct handed over by the dispatcher
    long long i;
    s.Ntime = (*sCPU).Ntime;
    s.Nspace = (*sCPU).Nspace;
    s.dt = (*sCPU).tStep;
    s.dx = (*sCPU).rStep;
    s.fStep = (*sCPU).fStep;
    s.h = (*sCPU).propagationStep;
    s.Nsteps = (*sCPU).Npropagation;
    s.Ngrid = s.Ntime * s.Nspace;
    s.axesNumber = (*sCPU).axesNumber;
    s.sellmeierType = (*sCPU).sellmeierType;
    s.f0 = (*sCPU).frequency1;
    s.Nthread = THREADS_PER_BLOCK;
    s.Nblock = s.Ngrid / THREADS_PER_BLOCK;
    s.isCylindric =(*sCPU).isCylindric;
    s.isNonLinear = ((*sCPU).nonlinearSwitches[0] + (*sCPU).nonlinearSwitches[1] + (*sCPU).nonlinearSwitches[2]) > 0;
    (*sCPU).nonlinearSwitches[3] = (int)ceil((*sCPU).absorptionParameters[0] * 241.79893e12 / (*sCPU).frequency1) - 1;
    //CPU allocations
    std::complex<double>* gridPropagationFactorCPU = (std::complex<double>*)malloc(2 * s.Ngrid * sizeof(std::complex<double>));
    std::complex<double>* gridPolarizationFactorCPU = (std::complex<double>*)malloc(2 * s.Ngrid * sizeof(std::complex<double>));

    //GPU allocations
    int memErrors = 0;
    memErrors += cudaMalloc((void**)&s.gridETime, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridETime2, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridETemp, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridETemp2, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridEFrequency, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridEFrequency2, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactor, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPolarizationFactor, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactor2, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPolarizationFactor2, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactorRho1, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactorRho2, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridRadialLaplacian1, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridRadialLaplacian2, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridEFrequencyNext1, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridEFrequencyNext2, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.k1, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.k2, sizeof(cuDoubleComplex) * s.Ngrid);
    //the following two should have a size (s.Ntime / 2 + 1) * s.Nspace, but I get overruns during
    //the ffts if they're not larger. If I figure this out, it will save a complex grid worth of memory...
    memErrors += cudaMalloc((void**)&s.gridPolarizationFrequency, sizeof(cuDoubleComplex) * s.Ngrid); 
    memErrors += cudaMalloc((void**)&s.gridPolarizationFrequency2, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPolarizationTime, sizeof(double) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPolarizationTime2, sizeof(double) * s.Ngrid);

    memErrors += cudaMalloc((void**)&s.chi2Tensor, sizeof(double) * 9);
    memErrors += cudaMalloc((void**)&s.firstDerivativeOperation, sizeof(double) * 6);
    memErrors += cudaMalloc((void**)&s.chi3Tensor, sizeof(double) * 81);
    memErrors += cudaMalloc((void**)&s.nonlinearSwitches, sizeof(int) * 4);
    memErrors += cudaMalloc((void**)&s.absorptionParameters, sizeof(double) * 6);
    memErrors += cudaMalloc((void**)&s.propagationInts, sizeof(long long) * 4);
    (*sCPU).memoryError = memErrors;

    //prepare effective nonlinearity tensors and put them on the GPU
    long long propagationIntsCPU[4] = { s.Ngrid, s.Ntime, s.Nspace, (s.Ntime / 2 + 1) };
    double firstDerivativeOperation[6] = { -(1/s.dx) / 60, (1 / s.dx) * 3 / 20, (1 / s.dx) * -3 / 4, (1 / s.dx) * 3 / 4, (1 / s.dx) * -3 / 20, (1 / s.dx) / 60 };
    deff((*sCPU).deffTensor, (*sCPU).chi2Tensor, (*sCPU).crystalTheta, (*sCPU).crystalPhi);
    cudaMemcpy(s.chi2Tensor, (*sCPU).deffTensor, 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.nonlinearSwitches, (*sCPU).nonlinearSwitches, 4 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(s.propagationInts, propagationIntsCPU, 4 * sizeof(long long), cudaMemcpyHostToDevice);
    cudaMemcpy(s.chi3Tensor, (*sCPU).chi3Tensor, 27 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.absorptionParameters, (*sCPU).absorptionParameters, 6 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.firstDerivativeOperation, firstDerivativeOperation, 6 * sizeof(double), cudaMemcpyHostToDevice);

    //prepare FFT plans
    cufftPlan2d(&s.fftPlan, s.Nspace, s.Ntime, CUFFT_Z2Z);
    cufftPlan2d(&s.polfftPlan, s.Nspace, s.Ntime, CUFFT_D2Z);

    //prepare the propagation arrays
    if (s.isCylindric) {
        preparepropagation3Dcylindric(sCPU, s);
    }
    else {
        preparepropagation2Dcartesian(sCPU, s);
    }
    

    //generate the pulses, either through pulsegenerator() if this is the first in the series, or by copying
    //the output of the last simulation in the sequence
    if ((*sCPU).isFollowerInSequence) {
        cudaMemcpy(s.gridETime, (*sCPU).ExtOut, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        cudaMemcpy(s.gridETime2, &(*sCPU).ExtOut[(*sCPU).Ngrid], (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        cudaMemcpy(s.gridEFrequency, (*sCPU).EkwOut, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        cudaMemcpy(s.gridEFrequency2, &(*sCPU).EkwOut[(*sCPU).Ngrid], (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    }
    else {
        pulsegenerator(sCPU, &s);
    }
    
    //Copy the field into the temporary array
    cudaMemcpy(s.gridETemp, s.gridEFrequency, s.Nspace * s.Ntime * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    cudaMemcpy(s.gridETemp2, s.gridEFrequency2, s.Nspace * s.Ntime * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);

    //Core propagation loop
    for (i = 0; i < s.Nsteps; i++) {
        //calculate k1
        rkstep(s, 0);
        //calculate k2
        rkstep(s, 1);
        //calculate k3
        rkstep(s, 2);
        //calculate k4
        rkstep(s, 3);

        if ((*sCPU).imdone[0] == 2) {
            break;
        }
    }
    cudaDeviceSynchronize();

    //transform final result
    fixnanKernel<<<s.Nblock, s.Nthread>>>(s.gridEFrequency);
    fixnanKernel << <s.Nblock, s.Nthread >> > (s.gridEFrequency2);
    cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridEFrequency, (cufftDoubleComplex*)s.gridETime, CUFFT_INVERSE);
    cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridEFrequency2, (cufftDoubleComplex*)s.gridETime2, CUFFT_INVERSE);
    fftNormalizeKernel<<<s.Nblock, s.Nthread >>>(s.gridETime, s.propagationInts);
    fftNormalizeKernel<<<s.Nblock, s.Nthread >>>(s.gridETime2, s.propagationInts);
    cudaDeviceSynchronize();

    //copy the field arrays from the GPU to CPU memory
    cudaMemcpy((*sCPU).ExtOut, s.gridETime, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*sCPU).EkwOut, s.gridEFrequency, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*sCPU).ExtOut[s.Ngrid], s.gridETime2, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*sCPU).EkwOut[s.Ngrid], s.gridEFrequency2, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    
    //Free GPU memory
    cudaFree(s.propagationInts);
    cudaFree(s.nonlinearSwitches);
    cudaFree(s.absorptionParameters);
    cudaFree(s.gridETime); 
    cudaFree(s.gridETemp);
    cudaFree(s.gridPolarizationFrequency);
    cudaFree(s.gridEFrequency);
    cudaFree(s.gridPropagationFactor);
    cudaFree(s.gridPropagationFactorRho1);
    cudaFree(s.gridPropagationFactorRho2);
    cudaFree(s.gridRadialLaplacian1);
    cudaFree(s.gridRadialLaplacian2);
    cudaFree(s.firstDerivativeOperation);
    cudaFree(s.gridPolarizationFactor);
    cudaFree(s.gridEFrequencyNext1);
    cudaFree(s.k1);
    cudaFree(s.gridPolarizationTime);
    cudaFree(s.gridETime2);
    cudaFree(s.gridETemp2);
    cudaFree(s.gridPolarizationFrequency2);
    cudaFree(s.gridEFrequency2);
    cudaFree(s.gridPropagationFactor2);
    cudaFree(s.gridPolarizationFactor2);
    cudaFree(s.gridEFrequencyNext2);
    cudaFree(s.k2);
    cudaFree(s.gridPolarizationTime2);
    cudaFree(s.chi2Tensor);
    cudaFree(s.chi3Tensor);
    cufftDestroy(s.fftPlan);
    cufftDestroy(s.polfftPlan);

    //Free CPU memory
    free(gridPropagationFactorCPU);
    free(gridPolarizationFactorCPU);
    
    return 0;
}

//function to run a RK4 time step
//stepNumber is the sub-step index, from 0 to 3
int rkstep(struct cudaLoop s, int stepNumber) {

    //operations involving FFT
    if (s.isNonLinear || s.isCylindric) {
        //perform inverse FFT to get time-space electric field
        cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridETemp, (cufftDoubleComplex*)s.gridETime, CUFFT_INVERSE);
        cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridETemp2, (cufftDoubleComplex*)s.gridETime2, CUFFT_INVERSE);
        
        if (s.isNonLinear) {
            nonlinearpolarizationKernel << <s.Nblock, s.Nthread >> > (s);
            cufftExecD2Z(s.polfftPlan, s.gridPolarizationTime, (cufftDoubleComplex*)s.gridPolarizationFrequency);
            cufftExecD2Z(s.polfftPlan, s.gridPolarizationTime2, (cufftDoubleComplex*)s.gridPolarizationFrequency2);
        }

        if (s.isCylindric) {
            cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridRadialLaplacian1, (cufftDoubleComplex*)s.k1, CUFFT_FORWARD);
            cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridRadialLaplacian2, (cufftDoubleComplex*)s.k2, CUFFT_FORWARD);
            radialLaplacianKernel<<<s.Nblock, s.Nthread >>>(s);
        }
    }

    //calculate k
    rkKernel<<<s.Nblock, s.Nthread >>>(s, stepNumber);
    
    return 0;
}

int pulsegenerator(struct propthread* s, struct cudaLoop *sc) {
    long long i,j;
    double rB, zB, r, z; //r and z in the Beam and lab coordinates, respectively.
    double w0, wz, zR, Rz, phi; //Gaussian beam parameters
    double theta = 0; //rotation angle of the current beam
    double pulseSum = 0;
    std::complex<double> ne, no, n0; //active refractive index;
    double f, w; //active frequency;
    double pulseEnergySum;
    std::complex<double> ko, k0, specfac, specphase;
    double c = 2.99792458e8; //speed of light
    double eps0 = 8.8541878128e-12; //vacuum permittivity
    double pi = 3.14159265358979323846264338327950288; // pi to unneccessary precision
    std::complex<double> *pulse1, *pulse2, *pulse1f, *pulse2f;
    cufftHandle plan1;
    cufftHandle plan2;
    pulse1 = (std::complex<double>*)calloc((*s).Ngrid * 2, sizeof(std::complex<double>));
    pulse2 = (std::complex<double>*)calloc((*s).Ngrid * 2, sizeof(std::complex<double>));
    pulse1f = (std::complex<double>*)calloc((*s).Ngrid * 2, sizeof(std::complex<double>));
    pulse2f = (std::complex<double>*)calloc((*s).Ngrid * 2, sizeof(std::complex<double>));
    std::complex<double> Eb;
    std::complex<double> ii(0, 1);



    std::complex<double> polFactor1, polFactor2; //complex phase/amplitude factors for the polarization components
    sellmeier(&n0, &no, (*s).sellmeierCoefficients, (*s).frequency1, (*s).crystalTheta, (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);
    (*s).neref = real(n0);
    (*s).noref = imag(n0);


    //define pulse 1 in mixed space
    // Gaussian beam in x
    // Spectrum in frequency domain (supergaussian with phase terms)
    polFactor1 = cos((*s).polarizationAngle1) - ii * (*s).circularity1 * sin((*s).polarizationAngle1);
    polFactor2 = sin((*s).polarizationAngle1) + ii * (*s).circularity1 * cos((*s).polarizationAngle1);
    theta = (*s).propagationAngle1;
    zB = (*s).z01;
    w0 = (*s).beamwaist1;

    for (i = 1; i < (*s).Ntime; i++) {
        f = i * (*s).fStep;
        if (i >= (*s).Ntime / 2) {
            f -= (*s).fStep * (*s).Ntime;
        }
        f *= -1;
        w = 2 * pi * (f - (*s).frequency1);
        
        //supergaussian pulse spectrum, if no input pulse specified
        specfac = (f - (*s).frequency1)/(*s).bandwidth1;
        for (j = 0; j < (*s).sgOrder1; j++) {
            specfac *= specfac;
        }
        specphase = ii * ((*s).cephase1 + w * (*s).delay1 - (*s).gdd1 * w * w - (*s).tod1 * w * w * w);
        specfac = exp(-specfac - specphase);

        if ((*s).field1IsAllocated) {
            specfac = (*s).loadedField1[i] * exp(-specphase);
        }

        ne = (*s).refractiveIndex1[i + (*s).Ntime * j];
        no = (*s).refractiveIndex2[i + (*s).Ntime * j];
        ko = 2 * pi * no * f / c;
        k0 = 2 * pi * real(n0) * f / c;
        zR = pi * w0 * w0 * real(ne) * f / c;
        if (f == 0) {
            zR = 1e3;
        }

        for (j = 0; j < (*s).Nspace; j++) {
            rB = (*s).x01 + (*s).rStep * j - (*s).Nspace* (*s).rStep / 2.;
            r = rB * cos(theta) - zB * sin(theta);
            z = rB * sin(theta) + zB * cos(theta);
            
            wz = w0 * sqrt(1 + (z * z / (zR * zR)));
            Rz = z * (1. + (zR * zR / (z * z)));
            
            if (z == 0) {
                Rz = 1.0e15;
            }
            phi = atan(z / zR);
            //z = 0;
            Eb = (w0 / wz) * exp(-ii * (real(ko) * (z-zB) + real(ko) * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
            Eb *= specfac;
            if (isnan(cmodulussquared(Eb)) || f<=0) {
                Eb = 0;
            }
            
            pulse1[i + (*s).Ntime * j] = polFactor1 * Eb;
            pulse1[i + (*s).Ntime * j + (*s).Ngrid] = polFactor2 * Eb;
            pulseSum += abs(r)*(real(ne)*cmodulussquared(pulse1[i + (*s).Ntime * j]) + real(no)*cmodulussquared(pulse1[i + (*s).Ntime * j + (*s).Ngrid]));
        }
    }
    
    // copy the field and propagation grids to the GPU
    cudaMemcpy((*sc).gridETime, pulse1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridETime2, &pulse1[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    // fft along spatial dimention to get Fourier space beam
    // will take place in three steps:
    // 2D fft (x,f)->(k,t), temporary intermediate state (could be optimized out later)
    // 1D fft (k,t)->(k,f), copied to Fourier space beam
    // 2D fft (k,f)->(x,t), copied to real space beam

    cufftPlan1d(&plan1, (*sc).Ntime, CUFFT_Z2Z, (*sc).Nspace);
    cufftPlan2d(&plan2, (*sc).Nspace, (*sc).Ntime, CUFFT_Z2Z);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime, (cufftDoubleComplex*)(*sc).gridETemp, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).gridETemp, (cufftDoubleComplex*)(*sc).gridEFrequency, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridEFrequency, (cufftDoubleComplex*)(*sc).gridETime, CUFFT_INVERSE);

    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime2, (cufftDoubleComplex*)(*sc).gridETemp2, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).gridETemp2, (cufftDoubleComplex*)(*sc).gridEFrequency2, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridEFrequency2, (cufftDoubleComplex*)(*sc).gridETime2, CUFFT_INVERSE);

    //Take the conjugate of the field because me and cufft have different ideas of time
    conjugateKernel<<<(*sc).Nblock, (*sc).Nthread >>>((*sc).gridETime);
    conjugateKernel<<<(*sc).Nblock, (*sc).Nthread >>>((*sc).gridETime2);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime, (cufftDoubleComplex*)(*sc).gridEFrequency, CUFFT_INVERSE);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime2, (cufftDoubleComplex*)(*sc).gridEFrequency2, CUFFT_INVERSE);
    cudaDeviceSynchronize();

    //Copy the GPU grids to the CPU memory
    cudaMemcpy(pulse1, (*sc).gridETime, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pulse1[(*s).Ngrid], (*sc).gridETime2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(pulse1f, (*sc).gridEFrequency, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pulse1f[(*s).Ngrid], (*sc).gridEFrequency2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    //normalize the pulse energy and set it to the input value
    pulseSum *= c * eps0;
    pulseSum *= 59.958 * pi; //59.958 is emperical factor
    pulseSum *= (*s).rStep / (*s).fStep;
    pulseEnergySum = sqrt((*s).pulseEnergy1/pulseSum)/(*s).Ngrid;
    
    for (i = 0; i < (*s).Ngrid * 2; i++) {
        pulse1[i] = pulse1[i] * pulseEnergySum;
        pulse1f[i] = pulse1f[i] * pulseEnergySum;
    }
    cudaDeviceSynchronize();

    //do same for pulse 2 here
    pulseSum = 0;
    polFactor1 = cos((*s).polarizationAngle2) - ii * (*s).circularity2 * sin((*s).polarizationAngle2);
    polFactor2 = sin((*s).polarizationAngle2) + ii * (*s).circularity2 * cos((*s).polarizationAngle2);
    theta = (*s).propagationAngle2;
    zB = (*s).z02;
    w0 = (*s).beamwaist2;

    for (i = 1; i < (*s).Ntime; i++) {
        f = i * (*s).fStep;
        if (i >= (*s).Ntime / 2) {
            f -= (*s).fStep * (*s).Ntime;
        }
        f *= -1;
        w = 2 * pi * (f - (*s).frequency2);

        //supergaussian pulse spectrum, if no input pulse specified
        specfac = (f - (*s).frequency2) / (*s).bandwidth2;
        for (j = 0; j < (*s).sgOrder1; j++) {
            specfac *= specfac;
        }
        specphase = ii * ((*s).cephase2 + w * (*s).delay2 - (*s).gdd2 * w * w - (*s).tod2 * w * w * w);
        specfac = exp(-specfac - specphase);

        if ((*s).field2IsAllocated) {
            specfac = (*s).loadedField2[i] * exp(-specphase);
        }


        ne = (*s).refractiveIndex1[i + (*s).Ntime * j];
        no = (*s).refractiveIndex2[i + (*s).Ntime * j];
        ko = 2 * pi * no * f / c;
        k0 = 2 * pi * real(n0) * f / c;
        zR = pi * w0 * w0 * real(ne) * f / c;
        if (f == 0) {
            zR = 1e3;
        }

        for (j = 0; j < (*s).Nspace; j++) {

            rB = (*s).x01 + (*s).rStep * j - (*s).Nspace * (*s).rStep / 2.;
            r = rB * cos(theta) - zB * sin(theta);
            z = rB * sin(theta) + zB * cos(theta);

            wz = w0 * sqrt(1 + (z * z / (zR * zR)));
            Rz = z * (1. + (zR * zR / (z * z)));

            if (z == 0) {
                Rz = 1.0e15;
            }
            phi = atan(z / zR);
            Eb = (w0 / wz) * exp(-ii * (real(ko) * (z - zB) + real(ko) * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
            Eb *= specfac;
            if (isnan(cmodulussquared(Eb)) || f <= 0) {
                Eb = 0;
            }

            pulse2[i + (*s).Ntime * j] = polFactor1 * Eb;
            pulse2[i + (*s).Ntime * j + (*s).Ngrid] = polFactor2 * Eb;
            pulseSum += abs(r) * (real(ne) * cmodulussquared(pulse2[i + (*s).Ntime * j]) + real(no) * cmodulussquared(pulse2[i + (*s).Ntime * j + (*s).Ngrid]));
        }
    }

    // copy the field and propagation grids to the GPU
    cudaMemcpy((*sc).gridETime, pulse2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridETime2, &pulse2[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    // fft along spatial dimention to get Fourier space beam
    // will take place in three steps:
    // 2D fft (x,f)->(k,t), temporary intermediate state (could be optimized out later)
    // 1D fft (k,t)->(k,f), copied to Fourier space beam
    // 2D fft (k,f)->(x,t), copied to real space beam

    cufftPlan1d(&plan1, (*sc).Ntime, CUFFT_Z2Z, (*sc).Nspace);
    cufftPlan2d(&plan2, (*sc).Nspace, (*sc).Ntime, CUFFT_Z2Z);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime, (cufftDoubleComplex*)(*sc).gridETemp, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).gridETemp, (cufftDoubleComplex*)(*sc).gridEFrequency, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridEFrequency, (cufftDoubleComplex*)(*sc).gridETime, CUFFT_INVERSE);

    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime2, (cufftDoubleComplex*)(*sc).gridETemp2, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).gridETemp2, (cufftDoubleComplex*)(*sc).gridEFrequency2, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridEFrequency2, (cufftDoubleComplex*)(*sc).gridETime2, CUFFT_INVERSE);

    //Take the conjugate of the field because me and cufft have different ideas of time
    conjugateKernel << <(*sc).Nblock, (*sc).Nthread >> > ((*sc).gridETime);
    conjugateKernel << <(*sc).Nblock, (*sc).Nthread >> > ((*sc).gridETime2);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime, (cufftDoubleComplex*)(*sc).gridEFrequency, CUFFT_INVERSE);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime2, (cufftDoubleComplex*)(*sc).gridEFrequency2, CUFFT_INVERSE);
    cudaDeviceSynchronize();

    //Copy the GPU grids to the CPU memory
    cudaMemcpy(pulse2, (*sc).gridETime, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pulse2[(*s).Ngrid], (*sc).gridETime2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(pulse2f, (*sc).gridEFrequency, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pulse2f[(*s).Ngrid], (*sc).gridEFrequency2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    //normalize the pulse energy and set it to the input value
    pulseSum *= c * eps0;
    pulseSum *= 59.958 * pi; //59.958 is emperical factor
    pulseSum *= (*s).rStep / (*s).fStep;
    pulseEnergySum = sqrt((*s).pulseEnergy2 / pulseSum) / (*s).Ngrid;

    for (i = 0; i < (*s).Ngrid * 2; i++) {
        pulse2[i] = pulse2[i] * pulseEnergySum;
        pulse2f[i] = pulse2f[i] * pulseEnergySum;
    }
    cudaDeviceSynchronize();

    //make the combined fields
    for (i = 0; i < (*s).Ngrid * 2; i++) {
        (*s).Ext[i] = pulse1[i] + pulse2[i];
        (*s).Ekw[i] = pulse1f[i] + pulse2f[i];
    }
    //Copy the grids back to the GPU
    cudaMemcpy((*sc).gridETime, (*s).Ext, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridETime2, &(*s).Ext[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridEFrequency, (*s).Ekw, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridEFrequency2, &(*s).Ekw[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cufftDestroy(plan1);
    cufftDestroy(plan2);


    return 0;
}

int preparepropagation2Dcartesian(struct propthread* s, struct cudaLoop sc) {
    //recycle allocated device memory for the grids needed
    double* alphaGPU = (double*)sc.gridEFrequencyNext1;
    double* sellmeierCoefficients = (double*)sc.k1;
    sc.ne = sc.gridEFrequencyNext2;
    sc.no = sc.k2;

    //construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
    double* sellmeierCoefficientsAugmentedCPU = (double*)calloc(66 + 8, sizeof(double));
    memcpy(sellmeierCoefficientsAugmentedCPU, (*s).sellmeierCoefficients, 66 * (sizeof(double)));
    sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
    sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
    sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
    sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
    sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
    sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep; 
    sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
    cudaMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, (66+8) * sizeof(double), cudaMemcpyHostToDevice);

    //prepare the propagation grids
    prepareCartesianGridsKernel <<<sc.Nblock, sc.Nthread >>> (alphaGPU, sellmeierCoefficients, sc);
    cudaDeviceSynchronize();

    //copy the retrieved refractive indicies to the cpu
    cudaMemcpy((*s).refractiveIndex1, sc.ne, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*s).refractiveIndex2, sc.no, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    
    //clean up
    cudaMemset(sc.gridEFrequency, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.gridEFrequency2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.k1, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.k2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    free(sellmeierCoefficientsAugmentedCPU);
    return 0;
}

int preparepropagation3Dcylindric(struct propthread* s, struct cudaLoop sc) {
    //recycle allocated device memory for the grids needed
    double* alphaGPU = (double*)sc.gridEFrequencyNext1;
    double* sellmeierCoefficients = (double*)sc.k1;
    sc.ne = sc.gridEFrequencyNext2;
    sc.no = sc.k2;

    //construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
    double* sellmeierCoefficientsAugmentedCPU = (double*)calloc(66 + 8, sizeof(double));
    memcpy(sellmeierCoefficientsAugmentedCPU, (*s).sellmeierCoefficients, 66 * (sizeof(double)));
    sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
    sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
    sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
    sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
    sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
    sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
    sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
    cudaMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);

    //prepare the propagation grids
    prepareCylindricGridsKernel << <sc.Nblock, sc.Nthread >> > (sellmeierCoefficients, sc);
    cudaDeviceSynchronize();

    //copy the retrieved refractive indicies to the cpu
    cudaMemcpy((*s).refractiveIndex1, sc.ne, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*s).refractiveIndex2, sc.no, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);


    //clean up
    cudaMemset(sc.gridEFrequency, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.gridEFrequency2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.k1, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.k2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    free(sellmeierCoefficientsAugmentedCPU);
    return 0;
}

double thetasearch(struct propthread* s, double dk, double f, double tol) {
    double theta=0;
    double dTheta = 0.1;
    double err, errPlus, errMinus;
    double rhs = 2.99792458e8 * dk / (2 * 3.14159265358979323846264338327950288 * f);
    std::complex<double> ne, no;
    double nePlus, neMinus;
    f = abs(f);
    sellmeier(&ne, &no, (*s).sellmeierCoefficients, f, (*s).crystalTheta + theta, (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);
    nePlus = real(ne);
    err = abs(nePlus * sin(theta) - rhs);
    int iters = 0;
    while (err > tol && iters<65536) {
        iters++;

        sellmeier(&ne, &no, (*s).sellmeierCoefficients, f, (*s).crystalTheta + theta + dTheta, (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);
        nePlus = real(ne);
        errPlus = abs(nePlus * sin(theta+dTheta) - rhs);

        sellmeier(&ne, &no, (*s).sellmeierCoefficients, f, (*s).crystalTheta + theta - dTheta, (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);
        neMinus = real(ne);
        errMinus = abs(neMinus * sin(theta-dTheta) - rhs);

        //Basic hill climbing algorithm
        //calculate the error at theta +/- dTheta
        // if theta + dTheta has lowest error, theta = theta+dTheta, err = errPlus
        // if theta - dTheta has lowest error, theta = theta-dTheta, err = errMinus
        // if theta has lowest error, step size is too large, dTheta /= 2;
        if (errPlus < err && errPlus < errMinus) {
            theta += dTheta;
            err = errPlus;
        }
        else if (errMinus < err) {
            theta -= dTheta;
            err = errMinus;
        }
        else {
            dTheta *= 0.5;
        }
    }
    return theta;
}
//rearrange a matrix from row major order to column major (not used, maybe broken)
int swaprc(double* M, int dim1, int dim2) {
    double* Ms = (double*)malloc(dim1 * dim2 * sizeof(double));
    int i, j;
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim2; j++) {
            Ms[i + j * dim1] = M[j + i * dim2];
        }
    }
    free(Ms);
    return 0;
}

int deff(double* defftensor, double* dtensor, double theta, double phi) {
    double delta = 0.; //this angle is used for biaxial crystals, but I'm ignorning it for the moment
    int i, j, k;
    //Rotation matrix between the angles of the electric field and the crystal axes
    double R[] = { cos(theta) * cos(phi) * cos(delta) - sin(phi) * sin(delta), cos(theta) * sin(phi) * cos(delta) + cos(phi) * sin(delta),
        -sin(theta) * cos(delta), -cos(theta) * cos(phi) * sin(delta) - sin(phi) * cos(delta),
        -cos(theta) * sin(phi) * sin(delta) + cos(phi) * cos(delta), sin(theta) * sin(delta) };

    //Matrix to translate the mixed field matrix in the reduced notation into the crystalline frame
    double Ore[] = { R[0] * R[0], R[1] * R[1], R[2] * R[2], 2 * R[1] * R[2], 2 * R[0] * R[2], 2 * R[0] * R[1],
        2 * R[0] * R[3], 2 * R[1] * R[4], 2 * R[2] * R[5], 2 * (R[4] * R[2] + R[1] * R[5]), 2 * (R[3] * R[2] + R[0] * R[5]), 2 * (R[3] * R[1] + R[0] * R[4]),
        R[3] * R[3], R[4] * R[4], R[5] * R[5], 2 * R[4] * R[5], 2 * R[3] * R[5], 2 * R[3] * R[4]
};

    //The deff tensor is given by the equation R deff = d Ore, solve for deff, find d Ore first
    double dOre[9] = { 0 };
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 6; k++) {
                dOre[i + 3 * j] += dtensor[i + 3 * k] * Ore[k + 6 * j];
            }
        }
    }
    //Least squares solution to get the deff tensor
    qr_solve_mat(3, 2, 3, R, dOre, defftensor);

    //correct cross-terms
    for (i = 2; i < 4; i++) {
        defftensor[i] *= 0.5;
    }

    for (i = 0; i < 6; i++) {
        defftensor[i] *= 2e-12; //change from pm/V to m/V and multiply by 2 for chi(2) instead of d
    }
    return 0;
}
//c implementation of fftshift, working on complex double precision
//A is the input array, B is the output
//dim1: column length
//dim2: row length
int fftshiftZ(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2) {
    long long i, j;
    long long div1 = dim1 / 2;
    long long div2 = dim2 / 2;
    //Quadrant 1
    for (i = 0; i < div1; i++) {
        for (j = 0; j < div2; j++) {
            B[i + dim1 * j] = A[i + div1 + dim1 * (j + div2)];
        }
    }
    //Quadrant 2
    for (i = 0; i < div1; i++) {
        for (j = div2; j < dim2; j++) {
            B[i + dim1 * j] = A[i + div1 + dim1 * (j-div2)];
        }
    }
    //Quadrant 3
    for (i = div1; i < dim1; i++) {
        for (j = 0; j < div2; j++) {
            B[i + dim1 * j] = A[i - div1 + dim1 * (j + div2)];
        }
    }
    //Quadrant 4
    for (i = div1; i < dim1; i++) {
        for (j = div2; j < dim2; j++) {
            B[i + dim1 * j] = A[i - div1 + dim1 * (j - div2)];
        }
    }
    return 0;
}

//same as fftshiftZ, but flips the output array columns
int fftshiftZflip(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2) {
    long long i, j;
    long long div1 = dim1 / 2;
    long long div2 = dim2 / 2;
    //Quadrant 1
    for (i = 0; i < div1; i++) {
        for (j = 0; j < div2; j++) {
            B[(dim1-i-1) + dim1 * j] = A[i + div1 + dim1 * (j + div2)];
        }
    }
    //Quadrant 2
    for (i = 0; i < div1; i++) {
        for (j = div2; j < dim2; j++) {
            B[(dim1 - i-1) + dim1 * j] = A[i + div1 + dim1 * (j - div2)];
        }
    }
    //Quadrant 3
    for (i = div1; i < dim1; i++) {
        for (j = 0; j < div2; j++) {
            B[(dim1 - i-1) + dim1 * j] = A[i - div1 + dim1 * (j + div2)];
        }
    }
    //Quadrant 4
    for (i = div1; i < dim1; i++) {
        for (j = div2; j < dim2; j++) {
            B[(dim1 - i-1) + dim1 * j] = A[i - div1 + dim1 * (j - div2)];
        }
    }
    return 0;
}

//sellmeier equation
//outputs are pointers ne and no
//a is a 16-value array containing the coefficients
//f is frequency (Hz)
//theta is the crystal angle
//phi is the other crystal angle (currently unused because biaxials haven't been implemented)
//type is the kind of crystal (0: isotropic, 1: uniaxial, 2:biaxial) 
//eqn will switch to a different equation, in the future, currently not implemented
//current equation form:
//n^2 = a[0] //background (high freq) contribution
//      + (a[1] + a[2] * lambda^2) / (lambda^2 + a[3]) + (a[4] + a[5] * lambda^2)/ (lambda^2 + a[6]) //two resonances, purely real contribution
//      + (a[7] + a[8] * lambda^2) / (lambda^2 + a[9]) + (a[10] + a[11] * lambda^2) / (lambda^2 + a[12]) //two more resonances
//      + a[13] * lambda^2 + a[14] * lambda^4 + a[15] * lambda^6 //parametrized low-frequency correction
//      + 4*pi*e^2*a[16]/(a[17] - omega^2 + i * a[18] * omega) // complex-valued Lorenzian contribution (a[17] to zero for Drude)
//      + 4*pi*e^2*a[19]/(a[20] - omega^2 + i * a[21] * omega) // complex-valued Lorenzian contribution (a[21] to zero for Drude)
std::complex<double> sellmeier(std::complex<double>* ne, std::complex<double>* no, double* a, double f, double theta, double phi, int type, int eqn) {
    if (f == 0) return 1; //exit immediately for f=0

    double c = 2.99792458e8; //speed of light
    double l = 1e6*c / f; //wavelength in microns
    double ls = l * l;
    std::complex<double> ii(0, 1);
    double pi = 3.14159265358979323846264338327950288;
    double omega = 2*pi*abs(f);
    double kL = 3183.9; //(e * e / (e_o *m_e)
    //option 0: isotropic
    if (type == 0) {
        ne[0] = a[0]
            + (a[1] + a[2] * ls) / (ls + a[3]) + (a[4] + a[5] * ls) / (ls + a[6])
            + (a[7] + a[8] * ls) / (ls + a[9]) + (a[10] + a[11] * ls) / (ls + a[12])
            + a[13] * ls + a[14] * ls * ls + a[15] * ls * ls * ls;
        if (real(ne[0]) < 1) {
            ne[0] = 1.;
        }
        ne[0] += kL * a[16] / (a[17] - omega * omega - ii * a[18] * omega)
            + kL * a[19] / (a[20] - omega * omega - ii * a[21] * omega);
        ne[0] = conj(sqrt(ne[0]));
        if (isnan(real(ne[0]))) {
            ne[0] = 1;
        }
        no[0] = ne[0];
        return ne[0];
    }
    //option 1: uniaxial
    else if (type == 1) {
        std::complex<double> na = (sqrt(a[0]
            + (a[1] + a[2] * ls) / (ls + a[3]) + (a[4] + a[5] * ls) / (ls + a[6]) + (a[7]
                + a[8] * ls) / (ls + a[9]) + (a[10] + a[11] * ls) / (ls + a[12])
            + a[13] * ls + a[14] * ls * ls + a[15] * ls * ls * ls
            + kL * a[16] / (a[17] - omega * omega + ii * a[18] * omega)
            + kL * a[19] / (a[20] - omega * omega + ii * a[21] * omega)));
        a = &a[22];
        std::complex<double> nb = (sqrt(a[0]
            + (a[1] + a[2] * ls) / (ls + a[3]) + (a[4] + a[5] * ls) / (ls + a[6]) + (a[7]
                + a[8] * ls) / (ls + a[9]) + (a[10] + a[11] * ls) / (ls + a[12])
            + a[13] * ls + a[14] * ls * ls + a[15] * ls * ls * ls
            + kL * a[16] / (a[17] - omega * omega - ii * a[18] * omega)
            + kL * a[19] / (a[20] - omega * omega - ii * a[21] * omega)));
        if (isnan(real(na)) || isnan(real(nb))) {
            no[0] = 1;
            ne[0] = 1;
            return 1;
        }
        no[0] = na;
        ne[0] = 1.0 / sqrt(cos(theta) * cos(theta) / (na * na) + sin(theta) * sin(theta) / (nb * nb));
        return na;
    }
    else {
        //later, implement biaxial crystals, for now just return 1;
        return 1;
    }
}

int loadfrogspeck(char* frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel, int fieldIndex) {
    FILE* fp;
    int maxFileSize = 16384;
    double wavelength, R, phi, complexX, complexY, f, f0, f1, fmax;
    int i, k0, k1;
    double c = 1e9*2.99792458e8; //for conversion of wavelength in nm to frequency
    double df = 0;
    double fmin = 0;
    int currentRow = 0;
    std::complex<double>* E = (std::complex<double>*)calloc(maxFileSize, sizeof(std::complex<double>));

    //read the data
    fp = fopen(frogFilePath, "r");
    while (fscanf(fp, "%lf %lf %lf %lf %lf", &wavelength, &R, &phi, &complexX, &complexY) == 5 && currentRow < maxFileSize) {
        //get the complex field from the data
        E[currentRow].real(complexX);
        E[currentRow].imag(complexY);

        //keep track of the frequency step of the grid (running sum, divide by number of rows at end to get average)
        if (currentRow > 0) df += c / wavelength - fmax;

        //keep track of the highest frequency in the data
        fmax = c / wavelength;
        
        //store the lowest frequency in the data
        if (currentRow == 0) fmin = fmax;
        
        currentRow++;
    }
    fclose(fp);

    //return an error if nothing was loaded
    if (currentRow == 0) {
        free(E);
        return -1;
    }

    df /= currentRow; //average frequency step

    //interpolate the FROG data onto the simulation grid
    
    //fill the simulation grid based on the data
    for (i = 0; i < Ntime; i++) {

        //frequency grid used in the simulation
        f = i * fStep;
        if (i >= Ntime / 2) {
            f -= fStep * Ntime;
        }
        f *= -1;

        k0 = (int)floor((f - fmin) / df);
        k1 = (int)ceil((f - fmin) / df);
        if (k0 < 0 || k1 >= currentRow) {
            Egrid[i] = 0; //field is zero outside of data range
        }
        else {
            f0 = fmin + k0 * df;
            f1 = fmin + k1 * df;
            Egrid[i] = (E[k0] * (f1 - f) + E[k1] * (f - f0)) / df; //linear interpolation
            Egrid[i] *= (abs(Egrid[i]) > gateLevel);
        }
    }

    free(E);
    return currentRow;
}

