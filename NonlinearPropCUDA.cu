#include "NonlinearPropCUDA.cuh"
#include <complex>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <cuComplex.h>
#include <cufft.h>
#include <mkl.h>
#include <thread>

#define THREADS_PER_BLOCK 64
#define FALSE 0
#define TRUE 1
#define MAX_LOADSTRING 1024

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

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
__device__ __forceinline__ cuDoubleComplex cuCexpd(cuDoubleComplex z){
    cuDoubleComplex expZ;
    double r = exp(z.x);
    expZ.y = sin(z.y);
    expZ.x = cos(z.y);
    expZ.x *= r;
    expZ.y *= r;
    return expZ;
}

//check if the pixel (xi,yi) is in the box defined by points (x1,y1) and (x2,y2)
__device__ __forceinline__ bool checkIfInBox(int xi, int yi, int x1, int y1, int x2, int y2, int tolerance) {
    int xmax = max(x1, x2) + tolerance;
    int xmin = min(x1, x2) - tolerance;
    int ymax = max(y1, y2) + tolerance;
    int ymin = min(y1, y2) - tolerance;
    return xi >= xmin && xi <= xmax && yi >= ymin && yi <= ymax;
}

//sqrt for complex doubles on CUDA, copy and paste from
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

//Inner function for the Sellmeier equation to provide the refractive indicies
//current equation form:
//n^2 = a[0] //background (high freq) contribution
//      + four resonances, purely real contribution
//      + parametrized low-frequency correction
//      + 2 complex-valued Lorenzian contribution
//inputs:
//a: 22 component array of the coefficients
//ls: lamda^2 (microns^2)
//omega: frequency (rad/s)
//ii: sqrt(-1)
//kL: 3183.9 i.e. (e * e / (epsilon_o * m_e)
__device__ __forceinline__ cuDoubleComplex sellmeierSubfunctionCuda(double* a, double ls, double omega, cuDoubleComplex ii, double kL) {
    double realPart = a[0]
        + (a[1] + a[2] * ls) / (ls + a[3])
        + (a[4] + a[5] * ls) / (ls + a[6])
        + (a[7] + a[8] * ls) / (ls + a[9])
        + (a[10] + a[11] * ls) / (ls + a[12])
        + a[13] * ls
        + a[14] * ls * ls
        + a[15] * ls * ls * ls;

    //traditional sellmeier part is not allowed to give complex values because that almost always
    //means it's out of range and causes instability
    if (realPart < 0) realPart = 1;

    return cuCsqrt(realPart
        + kL * a[16] / (a[17] - omega * omega + ii * a[18] * omega)
        + kL * a[19] / (a[20] - omega * omega + ii * a[21] * omega));
}
//Sellmeier equation for refractive indicies
__device__ cuDoubleComplex sellmeierCuda(cuDoubleComplex* ne, cuDoubleComplex* no, double* a, double f, double theta, double phi, int type, int eqn) {
    if (f==0) return make_cuDoubleComplex(1.0,0.0); //exit immediately for f=0

    double ls = 2.99792458e14 / f; //wavelength in microns
    ls *= ls; //only wavelength^2 is ever used
    cuDoubleComplex ii = make_cuDoubleComplex(0.0, 1.0);
    double omega = 6.28318530718 * abs(f);
    double kL = 3183.9; //(e * e / (epsilon_o * m_e)


    //option 0: isotropic
    if (type == 0) {
        ne[0] = sellmeierSubfunctionCuda(a, ls, omega, ii, kL);
        no[0] = ne[0];
        return ne[0];
    }
    //option 1: uniaxial
    else if (type == 1) {
        
        cuDoubleComplex na = sellmeierSubfunctionCuda(a, ls, omega, ii, kL);
        cuDoubleComplex nb = sellmeierSubfunctionCuda(&a[22], ls, omega, ii, kL);
        no[0] = na;
        ne[0] = 1.0 / cuCsqrt(cos(theta) * cos(theta) / (na * na) + sin(theta) * sin(theta) / (nb * nb));
        return ne[0];
    }
    else {
        //type == 2: biaxial
        // X. Yin, S. Zhang and Z. Tian, Optics and Laser Technology 39 (2007) 510 - 513.
        // I am sorry if there is a bug and you're trying to find it, i did my best.
        cuDoubleComplex na = sellmeierSubfunctionCuda(a, ls, omega, ii, kL);
        cuDoubleComplex nb = sellmeierSubfunctionCuda(&a[22], ls, omega, ii, kL);
        cuDoubleComplex nc = sellmeierSubfunctionCuda(&a[44], ls, omega, ii, kL);
        double cosTheta = cos(theta);
        double cosTheta2 = cosTheta * cosTheta;
        double sinTheta = sin(theta);
        double sinTheta2 = sinTheta * sinTheta;
        double sinPhi = sin(phi);
        double sinPhi2 = sinPhi * sinPhi;
        double cosPhi = cos(phi);
        double cosPhi2 = cosPhi * cosPhi;
        double realna2 = cuCreal(na) * cuCreal(na);
        double realnb2 = cuCreal(nb) * cuCreal(nb);

        double delta = 0.5 * atan(-((1. / realna2 - 1. / realnb2) 
            * sin(2 * phi) * cosTheta) / ((cosPhi2 / realna2 + sinPhi2 / realnb2) 
            + ((sinPhi2 /realna2 + cosPhi2 / realnb2) 
                * cosTheta2 + sinTheta2 / (cuCreal(nc) * cuCreal(nc)))));

        ne[0] = 1.0/cuCsqrt(cos(delta) * cos(delta) * (cosTheta2 * (cosPhi2 / (na * na) 
            + sinPhi2 / (nb * nb)) + sinTheta2 / (nc*nc)) 
            + sin(delta) * sin(delta) * (sinPhi2 / (na*na) + cosPhi2 / (nb*nb)) 
            - 0.5*sin(2 * phi)*cosTheta*sin(2 * delta)*(1. / (na*na) - 1. / (nb*nb)));

        no[0] = 1.0/cuCsqrt(sin(delta) * sin(delta) * (cosTheta2 * (cosPhi2 / (na * na) 
            + sinPhi2 / (nb * nb)) + sinTheta2 / (nc*nc)) 
            + cos(delta)*cos(delta) * (sinPhi2 / (na*na) + cosPhi2 / (nb*nb)) 
            + 0.5 * sin(2 * phi)*cosTheta*sin(2 * delta)*(1. / (na*na) - 1. / (nb*nb)));
        return ne[0];
    }
}

//rotate the field around the propagation axis (basis change)
__global__ void rotateFieldKernel(cuDoubleComplex* Ein1, cuDoubleComplex* Ein2, cuDoubleComplex* Eout1, cuDoubleComplex* Eout2, double rotationAngle) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    Eout1[i] = cos(rotationAngle) * Ein1[i] - sin(rotationAngle) * Ein2[i];
    Eout2[i] = sin(rotationAngle) * Ein1[i] + cos(rotationAngle) * Ein2[i];
}

//Radial laplacian term of the nonlinear wave equation
// e.g., apply the 1/rho * d/drho operator
// this isn't really the optimal way to do this; revisit later
__global__ void radialLaplacianKernel(struct cudaParameterSet s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j = i / s.Ntime; //spatial coordinate
    double rho = s.dx * j - (s.dx / 2) * s.Nspace;

    //zero at edges of grid and at origin
    if ( j<3 || j>(s.Nspace-4)) {
        s.gridRadialLaplacian1[i] = make_cuDoubleComplex(0, 0);
        s.gridRadialLaplacian2[i] = make_cuDoubleComplex(0, 0);
    }
    else if (rho != 0.){
        rho = 1.0/rho;
        s.gridRadialLaplacian1[i] = rho*(s.firstDerivativeOperation[0] * s.gridETime1[i - 3 * s.Ntime]
            + s.firstDerivativeOperation[1] * s.gridETime1[i - 2 * s.Ntime]
            + s.firstDerivativeOperation[2] * s.gridETime1[i - s.Ntime]
            + s.firstDerivativeOperation[3] * s.gridETime1[i + s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime1[i + 2 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime1[i + 3 * s.Ntime]);
        s.gridRadialLaplacian2[i] = rho*(s.firstDerivativeOperation[0] * s.gridETime2[i - 3 * s.Ntime]
            + s.firstDerivativeOperation[1] * s.gridETime2[i - 2 * s.Ntime]
            + s.firstDerivativeOperation[2] * s.gridETime2[i - s.Ntime]
            + s.firstDerivativeOperation[3] * s.gridETime2[i + s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime2[i + 2 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime2[i + 3 * s.Ntime]);
    }
    else {
        //handle rho = 0 by Fornberg interpolation from surrounding points
        //r[0]= 1.5 r[1] - 0.6 r[2] + 0.1 r[3]
        s.gridRadialLaplacian1[i] = (1.5/(rho+s.dx)) * (s.firstDerivativeOperation[0] * s.gridETime1[i - 2 * s.Ntime]
            + s.firstDerivativeOperation[1] * s.gridETime1[i - 1 * s.Ntime]
            + s.firstDerivativeOperation[2] * s.gridETime1[i]
            + s.firstDerivativeOperation[3] * s.gridETime1[i + 2 *s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime1[i + 3 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime1[i + 4 * s.Ntime]);
        s.gridRadialLaplacian2[i] = (1.5 / (rho + s.dx)) * (s.firstDerivativeOperation[0] * s.gridETime2[i - 2 * s.Ntime]
            + s.firstDerivativeOperation[1] * s.gridETime2[i - 1 * s.Ntime]
            + s.firstDerivativeOperation[2] * s.gridETime2[i]
            + s.firstDerivativeOperation[3] * s.gridETime2[i + 2*s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime2[i + 3 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime2[i + 4 * s.Ntime]);

        s.gridRadialLaplacian1[i] = s.gridRadialLaplacian1[i] + (-0.6 / (rho + 2 * s.dx)) * (s.firstDerivativeOperation[0] * s.gridETime1[i - s.Ntime]
            + s.firstDerivativeOperation[1] * s.gridETime1[i]
            + s.firstDerivativeOperation[2] * s.gridETime1[i + s.Ntime]
            + s.firstDerivativeOperation[3] * s.gridETime1[i + 3 * s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime1[i + 4 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime1[i + 5 * s.Ntime]);
        s.gridRadialLaplacian2[i] = s.gridRadialLaplacian2[i] + (-0.6 / (rho + 2*s.dx)) * (s.firstDerivativeOperation[0] * s.gridETime2[i -  s.Ntime]
            + s.firstDerivativeOperation[1] * s.gridETime2[i]
            + s.firstDerivativeOperation[2] * s.gridETime2[i + s.Ntime]
            + s.firstDerivativeOperation[3] * s.gridETime2[i + 3 * s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime2[i + 4 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime2[i + 5 * s.Ntime]);

        s.gridRadialLaplacian1[i] = s.gridRadialLaplacian1[i] + (0.1 / (rho + 3*s.dx)) * (s.firstDerivativeOperation[0] * s.gridETime1[i]
            + s.firstDerivativeOperation[1] * s.gridETime1[i + 1 * s.Ntime]
            + s.firstDerivativeOperation[2] * s.gridETime1[i + 2 * s.Ntime]
            + s.firstDerivativeOperation[3] * s.gridETime1[i + 4 * s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime1[i + 5 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime1[i + 6 * s.Ntime]);
        s.gridRadialLaplacian2[i] = s.gridRadialLaplacian2[i] + (0.1 / (rho + 3*s.dx)) * (s.firstDerivativeOperation[0] * s.gridETime2[i]
            + s.firstDerivativeOperation[1] * s.gridETime2[i + 1 * s.Ntime]
            + s.firstDerivativeOperation[2] * s.gridETime2[i + 2 * s.Ntime]
            + s.firstDerivativeOperation[3] * s.gridETime2[i + 4 * s.Ntime]
            + s.firstDerivativeOperation[4] * s.gridETime2[i + 5 * s.Ntime]
            + s.firstDerivativeOperation[5] * s.gridETime2[i + 6 * s.Ntime]);
    }

}

//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
__global__ void prepareCartesianGridsKernel(double* theta, double* sellmeierCoefficients, struct cudaParameterSet s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j, k;
    long long Ntime = s.Ntime;
    long long Nspace = s.Nspace;
    int axesNumber = s.axesNumber;
    int sellmeierType = s.sellmeierType;
    double c = 2.99792458e8; //speed of light
    double twoPi = 2 * 3.14159265358979323846264338327950288;
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
    double rhs = 2.99792458e8 * dk / (twoPi * f);
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
        ne = make_cuDoubleComplex(1, 0);
        no = make_cuDoubleComplex(1, 0);
    }
    s.ne[i] = ne;
    s.no[i] = no;

    cuDoubleComplex k0 = make_cuDoubleComplex(twoPi * cuCreal(n0) * f / c, 0);
    cuDoubleComplex ke = twoPi * ne * f / c;
    cuDoubleComplex ko = twoPi * no * f / c;

    if (cuCreal(ke) < 0 && cuCreal(ko) < 0) {
        s.gridPropagationFactor1[i] = ii * (ke - k0 + dk * dk / (2. * cuCreal(ke))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor1[i]))) {
            s.gridPropagationFactor1[i] = cuZero;
        }

        s.gridPropagationFactor2[i] = ii * (ko - k0 + dk * dk / (2. * cuCreal(ko))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor2[i]))) {
            s.gridPropagationFactor2[i] = cuZero;
        }

        s.gridPolarizationFactor1[i] = ii * (twoPi * f) / (2. * cuCreal(ne) * c) * s.h;
        s.gridPolarizationFactor2[i] = ii * (twoPi * f) / (2. * cuCreal(no) * c) * s.h;
    }

    else {
        s.gridPropagationFactor1[i] = cuZero;
        s.gridPropagationFactor2[i] = cuZero;
        s.gridPolarizationFactor1[i] = cuZero;
        s.gridPolarizationFactor2[i] = cuZero;
    }
}
    
//prepare the propagation constants under the assumption of cylindrical symmetry of the beam
__global__ void prepareCylindricGridsKernel(double* sellmeierCoefficients, struct cudaParameterSet s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j, k;
    long long Ntime = s.Ntime;
    long long Nspace = s.Nspace;
    int axesNumber = s.axesNumber;
    int sellmeierType = s.sellmeierType;
    double c = 2.99792458e8; //speed of light
    double twoPi = 2 * 3.14159265358979323846264338327950288;
    cuDoubleComplex cuZero = make_cuDoubleComplex(0, 0);
    j = i / Ntime; //spatial coordinate
    k = i - j * Ntime; //temporal coordinate
    cuDoubleComplex ii = make_cuDoubleComplex(0, 1);
    double crystalTheta = sellmeierCoefficients[66];
    double crystalPhi = sellmeierCoefficients[67];
    double kStep = sellmeierCoefficients[70];
    double fStep = sellmeierCoefficients[71];

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

    cuDoubleComplex k0 = make_cuDoubleComplex(twoPi * cuCreal(n0) * f / c, 0);
    cuDoubleComplex ke = twoPi * ne * f / c;
    cuDoubleComplex ko = twoPi * no * f / c;

    if (cuCreal(ke) < 0 && cuCreal(ko) < 0) {
        s.gridPropagationFactor1[i] = ii * (ke - k0 + dk * dk / (2. * cuCreal(ke))) * s.h;
        s.gridPropagationFactor1Rho1[i] = ii * (1 / (2. * cuCreal(ke))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor1[i]))) {
            s.gridPropagationFactor1[i] = cuZero;
            s.gridPropagationFactor1Rho1[i] = cuZero;
        }

        s.gridPropagationFactor2[i] = ii * (ko - k0 + dk * dk / (2. * cuCreal(ko))) * s.h;
        s.gridPropagationFactor1Rho2[i] = ii * (1 / (2. * cuCreal(ko))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor2[i]))) {
            s.gridPropagationFactor2[i] = cuZero;
            s.gridPropagationFactor1Rho2[i] = cuZero;
        }

        s.gridPolarizationFactor1[i] = ii * (twoPi * f) / (2. * cuCreal(ne) * c) * s.h;
        s.gridPolarizationFactor2[i] = ii * (twoPi * f) / (2. * cuCreal(no) * c) * s.h;
    }

    else {
        s.gridPropagationFactor1[i] = cuZero;
        s.gridPropagationFactor2[i] = cuZero;
        s.gridPolarizationFactor1[i] = cuZero;
        s.gridPolarizationFactor2[i] = cuZero;
        s.gridPropagationFactor1[i] = cuZero;
        s.gridPropagationFactor1Rho2[i] = cuZero;
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

//calculate the nonlinear polarization, after FFT to get the field
//in the time domain
__global__ void nonlinearPolarizationKernel(struct cudaParameterSet s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    double Ex = 2*cuCreal(s.gridETime1[i]) / s.propagationInts[0];
    double Ey = 2*cuCreal(s.gridETime2[i]) / s.propagationInts[0];
    s.gridPolarizationTime1[i] = 0.;
    s.gridPolarizationTime2[i] = 0.;

    //The d2eff tensor has the form
    // | d_xxx d_xyx d_yyx |
    // | d_xxy d_xyy d_yyy |
    if (s.nonlinearSwitches[0] == 1) {
        s.gridPolarizationTime1[i] += s.chi2Tensor[0] * Ex * Ex + s.chi2Tensor[2] * Ex * Ey + s.chi2Tensor[4] * Ey * Ey;
        s.gridPolarizationTime2[i] += s.chi2Tensor[1] * Ex * Ex + s.chi2Tensor[3] * Ex * Ey + s.chi2Tensor[5] * Ey * Ey;
    }
    
    //to be implemented: full chi3 matrix on s.nonlinearSwitches[1]==1

    //using only one value of chi3, under assumption of centrosymmetry
    if (s.nonlinearSwitches[1] == 2) {
        s.gridPolarizationTime1[i] += s.chi3Tensor[0] * (Ex * Ex * Ex + (2. / 3.) * Ey * Ey * Ex);
        s.gridPolarizationTime2[i] += s.chi3Tensor[0] * (Ey * Ey * Ey + (2. / 3.) * Ex * Ex * Ey);
    }
}


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
__global__ void plasmaCurrentKernelPrep(struct cudaParameterSet s, double* workN, double* workEx) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;

    int k;
    double* workEy = &workEx[s.Ngrid];
    double w, Esquared, Ex, Ey;
    Ex = cuCreal(s.gridETime1[i]) / s.propagationInts[0];
    Ey = cuCreal(s.gridETime2[i]) / s.propagationInts[0];
    Esquared = Ex * Ex + Ey * Ey;
    //plasmaParameters[0] is the nonlinear absorption parameter
    w = s.plasmaParameters[0] * Esquared;
    //nonlinearSwitches[3] is Nphotons-2
    for (k = 0; k < s.nonlinearSwitches[3]; k++) {
        w *= Esquared;
    }
    //absorption currents
    s.gridPlasmaCurrent1[i] = w * Ex;
    s.gridPlasmaCurrent2[i] = w * Ey;

    //plasmaParameters[2] is the 1/photon energy, translating the loss of power
    //from the field to the number of free carriers
    //extra factor of (dt^2e^2/(m*photon energy*eo) included as it is needed for the amplitude
    //of the plasma current
    workN[i] = s.plasmaParameters[2] * (s.gridPlasmaCurrent1[i] * Ex + s.gridPlasmaCurrent2[i] * Ey);
    workEx[i] = Ex;
    workEy[i] = Ey;

}
__global__ void plasmaCurrentKernel2(struct cudaParameterSet s, double* workN, double* workEx) {
    long long j = threadIdx.x + blockIdx.x * blockDim.x;
    double N = 0;
    double integralx = 0;
    double integraly = 0;
    double* workEy = &workEx[s.Ngrid];
    double* expMinusGammaT = &s.expGammaT[s.Ntime];

    long long k, l;
    j *= s.Ntime;
    for (k = 0; k < s.Ntime; k++) {

        l = j + k;
        N += workN[l];

        integralx += s.expGammaT[k] * N * workEx[l];
        integraly += s.expGammaT[k] * N * workEy[l];


        s.gridPlasmaCurrent1[l] += expMinusGammaT[k] * integralx;
        s.gridPlasmaCurrent2[l] += expMinusGammaT[k] * integraly;
    }
}


//Main kernel for RK4 propagation of the field
__global__ void rkKernel(struct cudaParameterSet s, int stepNumber) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j = i / s.Ntime; //spatial coordinate
    long long h = i - j * s.Ntime; //temporal coordinate
    cuDoubleComplex plasmaJ1 = make_cuDoubleComplex(0, 0);
    cuDoubleComplex plasmaJ2 = make_cuDoubleComplex(0, 0);

    //note that the FFT of the radial laplacian is stored in k1 and k2
    //so that memory shouldn't be used for anything else
    if (s.isCylindric) {
        s.gridRadialLaplacian1[i] = s.gridPropagationFactor1Rho1[i] * s.k1[i];
        s.gridRadialLaplacian2[i] = s.gridPropagationFactor1Rho2[i] * s.k2[i];
    }

    if (s.hasPlasma) {
        double f = h * s.fStep;
        long long hp = h;
        if (h >= s.Ntime / 2) {
            f -= s.fStep * s.Ntime;
        }
        f *= -6.28318530718;
        cuDoubleComplex jfac = make_cuDoubleComplex(0, 1.0 / f);
        if (h > s.propagationInts[3]) {
            hp = s.Ntime - hp;
            j = s.Nspace - j;
            hp += j * s.propagationInts[3];

            if (f != 0) {
                plasmaJ1 = jfac * s.gridPolarizationFactor1[i] * cuConj(s.gridPlasmaCurrentFrequency1[hp]);
                plasmaJ2 = jfac * s.gridPolarizationFactor2[i] * cuConj(s.gridPlasmaCurrentFrequency2[hp]);
            }
        }
        else {
            hp += j * s.propagationInts[3];
            if (f != 0) {
                plasmaJ1 = jfac * s.gridPolarizationFactor1[i] * s.gridPlasmaCurrentFrequency1[hp];
                plasmaJ2 = jfac * s.gridPolarizationFactor2[i] * s.gridPlasmaCurrentFrequency2[hp];
            }
        }

    }


    //polarization is stored in a reduced format by cuFFT because the FFT is from real to complex, meaning if the output grid
    //were to be N_time x N_space, half of the points would be redundant. The extra steps below are to determine where in the grid the 
    //current point sits. Essentially, if in the negative frequency quadrants, reverse the frequency and take complex conjugate of the 
    //value
    if (h > s.propagationInts[3]) {
        h = s.Ntime - h;
        j = s.Nspace - j;
        h += j * s.propagationInts[3];

        s.k1[i] = s.gridPropagationFactor1[i] * s.gridETemp1[i] +s.gridPolarizationFactor1[i] * cuConj(s.gridPolarizationFrequency1[h]);
        s.k2[i] = s.gridPropagationFactor2[i] * s.gridETemp2[i] +s.gridPolarizationFactor2[i] * cuConj(s.gridPolarizationFrequency2[h]);
    }
    else {
        h += j * s.propagationInts[3];

        s.k1[i] = s.gridPropagationFactor1[i] * s.gridETemp1[i] +s.gridPolarizationFactor1[i] * s.gridPolarizationFrequency1[h];
        s.k2[i] = s.gridPropagationFactor2[i] * s.gridETemp2[i] +s.gridPolarizationFactor2[i] * s.gridPolarizationFrequency2[h];
    }
    if (s.isCylindric) {
        s.k1[i] = s.k1[i] + s.gridRadialLaplacian1[i];
        s.k2[i] = s.k2[i] + s.gridRadialLaplacian2[i];
    }

    if (s.hasPlasma) {
        s.k1[i] = s.k1[i] + plasmaJ1;
        s.k2[i] = s.k2[i] + plasmaJ2;
    }

    //in the first substep, first construct the next intermediate field value
    //which will be used in the next substep. 
    if (stepNumber == 0) {
        s.gridETemp1[i] = s.gridEFrequency1[i] + 0.5 * s.k1[i];
        s.gridETemp2[i] = s.gridEFrequency2[i] + 0.5 * s.k2[i];
       
        s.gridEFrequency1Next1[i] = s.k1[i] / 6 + s.gridEFrequency1[i];
        s.gridEFrequency1Next2[i] = s.k2[i] / 6 + s.gridEFrequency2[i];
    }

    //in the next substep, again construct the next intermediate field and add k/3 to solution
    else if (stepNumber == 1) {
        s.gridETemp1[i] = s.gridEFrequency1[i] + 0.5 * s.k1[i];
        s.gridETemp2[i] = s.gridEFrequency2[i] + 0.5 * s.k2[i];

        s.gridEFrequency1Next1[i] = s.gridEFrequency1Next1[i] + s.k1[i] / 3;
        s.gridEFrequency1Next2[i] = s.gridEFrequency1Next2[i] + s.k2[i] / 3;

    }

    //same action as previous substep, except the weight of k in the intermediate solution is 1 instead of 0.5
    else if (stepNumber == 2) {
        s.gridETemp1[i] = s.gridEFrequency1[i] + s.k1[i];
        s.gridETemp2[i] = s.gridEFrequency2[i] + s.k2[i];
        s.gridEFrequency1Next1[i] = s.gridEFrequency1Next1[i] + s.k1[i] / 3;
        s.gridEFrequency1Next2[i] = s.gridEFrequency1Next2[i] + s.k2[i] / 3;
    }

    //last substep. Solution is now complete and may be copied directly into the field arrays
    else {
        s.gridEFrequency1[i] = s.gridEFrequency1Next1[i] + s.k1[i] / 6;
        s.gridEFrequency2[i] = s.gridEFrequency1Next2[i] + s.k2[i] / 6;
        s.gridETemp1[i] = s.gridEFrequency1[i];
        s.gridETemp2[i] = s.gridEFrequency2[i];
    }

}

//Line plot calculation kernel
//runs through all the points and figures out the color of the current pixel
__global__ void plotDataKernel(double* dataX, double* dataY, double lineWidth, double markerWidth, int sizeData, double* plotGrid, double normFactorX, double normFactorY, double offsetX, double offsetY, int sizeX, int sizeY, double* xTicks, int NxTicks, double* yTicks, int NyTicks) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = i / sizeX;
    int k = i - j * sizeX;
    double x1, x2, y1, y2;
    int c;
    int tickLength = 12;
    double positionNormalization;
    double lineDistance;
    double pointDistance;
    plotGrid[i] = 0;

    //draw y ticks
    for (c = 0; c < NyTicks; c++) {
        x1 = 0;
        x2 = tickLength;
        y1 = normFactorY * (yTicks[c] + offsetY);
        y2 = y1;
        if (checkIfInBox(k, j, (int)x1, (int)y1, (int)x2, (int)y2, (int)(6*lineWidth))) {
            if (k < x1 || k > x2) {
                lineDistance = min(sqrt((x1 - k) * (x1 - k) + (y1 - j) * (y1 - j)), sqrt((x2 - k) * (x2 - k) + (y2 - j) * (y2 - j)));
            }
            else {
                positionNormalization = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                lineDistance = abs((x2 - x1) * (y1 - j) - (x1 - k) * (y2 - y1)) / positionNormalization;
            }
            lineDistance /= lineWidth;
            plotGrid[i] = max(plotGrid[i], exp(-lineDistance * lineDistance));
        }
    }

    //draw x ticks
    for (c = 0; c < NxTicks; c++) {
        x1 = normFactorX * (xTicks[c] + offsetX);
        x2 = x1;
        y1 = 0;
        y2 = tickLength;
        if (checkIfInBox(k, j, (int)x1, (int)y1, (int)x2, (int)y2, (int)(6 * lineWidth))) {
            if (k < (x1-1) || k > (x2+1)) {
                lineDistance = min(sqrt((x1 - k) * (x1 - k) + (y1 - j) * (y1 - j)), sqrt((x2 - k) * (x2 - k) + (y2 - j) * (y2 - j)));
            }
            else {
                positionNormalization = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                lineDistance = abs((x2 - x1) * (y1 - j) - (x1 - k) * (y2 - y1)) / positionNormalization;
            }
            lineDistance /= lineWidth;
            plotGrid[i] = max(plotGrid[i], exp(-lineDistance * lineDistance));
        }
    }

    for (c = 0; c < sizeData-1; c++) {
        x1 = normFactorX*(dataX[c] + offsetX);
        x2 = normFactorX*(dataX[c + 1] + offsetX);
        y1 = normFactorY*(dataY[c] + offsetY);
        y2 = normFactorY*(dataY[c + 1] + offsetY);
        pointDistance = sqrt((x1 - k) * (x1 - k) + (y1 - j) * (y1 - j));
        if (checkIfInBox(k, j, (int)x1, (int)y1, (int)x2, (int)y2, (int)(6 * lineWidth))) {
            
            //if (k < x1 || k > x2 || ((j > (y1 + 1)) && (j > (y2 + 1))) || ((j < (y1 - 1)) && (j < (y2 - 1)))) {
            if ((k < x1-1) || (k > x2+1) || ((j > (y1 + 1)) && (j > (y2 + 1))) || ((j < (y1 - 1)) && (j < (y2 - 1)))){
                lineDistance = min(sqrt((x1 - k) * (x1 - k) + (y1 - j) * (y1 - j)), sqrt((x2 - k) * (x2 - k) + (y2 - j) * (y2 - j)));
            }
            else {
                positionNormalization = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                lineDistance = abs((x2 - x1) * (y1 - j) - (x1 - k) * (y2 - y1)) / positionNormalization;
            }
            

            //end add
            lineDistance /= lineWidth;
            pointDistance /= markerWidth;
            plotGrid[i] = max(plotGrid[i], exp(-lineDistance * lineDistance));
            plotGrid[i] = max(plotGrid[i], exp(-pointDistance * pointDistance * pointDistance * pointDistance));
        }
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

//main function for running on CLI
//to be filled in!
int main(int argc, char *argv[]) {
    int i, j;
    int CUDAdevice;
    int CUDAdeviceCount = 0;
    cudaGetDeviceCount(&CUDAdeviceCount);
    cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
    struct cudaDeviceProp activeCUDADeviceProp;
    if (cuErr == cudaSuccess) {
        printf("Found %i GPU(s): \n", CUDAdeviceCount);
        for (i = 0; i < CUDAdeviceCount; i++) {
            cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
            printf("%s\r\n", activeCUDADeviceProp.name);
            printf(" Memory: %lli MB; Multiprocessors: %i\n", activeCUDADeviceProp.totalGlobalMem / (1024 * 1024), activeCUDADeviceProp.multiProcessorCount);
        }
    }
    else {
        printf("No GPU found.\n");
        return 1;
    }
    
    if (argc < 2) {
        printf("no input file specified.\n");
        return 2;
    }

    // allocate databases, main structs
    struct simulationParameterSet* sCPU = (struct simulationParameterSet*)calloc(512, sizeof(struct simulationParameterSet));
    struct crystalEntry* crystalDatabasePtr = (struct crystalEntry*)calloc(512, sizeof(struct crystalEntry));
    (*sCPU).sequenceString = (char*)calloc(256 * MAX_LOADSTRING, sizeof(char));
    (*sCPU).outputBasePath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    (*sCPU).field1FilePath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    (*sCPU).field2FilePath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    (*sCPU).crystalDatabase = crystalDatabasePtr;

    // read crystal database
    if (readCrystalDatabase(crystalDatabasePtr) == -2) {
        return 11;
    }
    if ((*crystalDatabasePtr).numberOfEntries == 0) {
        printf("Could not read crystal database.\n");
        free((*sCPU).sequenceString);
        free((*sCPU).outputBasePath);
        free(sCPU);
        free(crystalDatabasePtr);
        return 12;
    }
    printf("Read %i crystal database entries:\n", (*crystalDatabasePtr).numberOfEntries);
    for (j = 0; j < (*crystalDatabasePtr).numberOfEntries; j++) {
        printf("Material %i name: %ls", j, crystalDatabasePtr[j].crystalNameW);
    }
    


    // read from settings file
    if (readInputParametersFile(sCPU, crystalDatabasePtr, argv[1]) == 1) {
        printf("Could not read input file.\n");

        free((*sCPU).sequenceString);
        free((*sCPU).outputBasePath);
        free(sCPU);
        free(crystalDatabasePtr);
        return 13;
    }

    allocateGrids(sCPU);
    if (loadPulseFiles(sCPU) == 1) {
        printf("Could not read pulse file.\n");
        free((*sCPU).sequenceString);
        free((*sCPU).sequenceArray);
        free((*sCPU).refractiveIndex1);
        free((*sCPU).refractiveIndex2);
        free((*sCPU).imdone);
        free((*sCPU).deffTensor);
        free((*sCPU).loadedField1);
        free((*sCPU).loadedField2);
        free((*sCPU).outputBasePath);
        free((*sCPU).field1FilePath);
        free((*sCPU).field2FilePath);
        free(sCPU);
        free(crystalDatabasePtr);
        return 14;
    }

    readSequenceString(sCPU);
    configureBatchMode(sCPU);

    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();
    // run simulations
    std::thread *threadBlock = (std::thread*)calloc((*sCPU).Nsims, sizeof(std::thread));
    size_t maxThreads = min(CUDAdeviceCount, (*sCPU).Nsims);
    for (j = 0; j < (*sCPU).Nsims; j++) {

        sCPU[j].assignedGPU = j % CUDAdeviceCount;
        if (j >= maxThreads) {
            if (threadBlock[j - maxThreads].joinable()) {
                threadBlock[j - maxThreads].join();
            }
        }

        if ((*sCPU).isInSequence) {
            threadBlock[j] = std::thread(solveNonlinearWaveEquationSequence, &sCPU[j]);
        }
        else {
            threadBlock[j] = std::thread(solveNonlinearWaveEquation, &sCPU[j]);
        }
    }
    
	for (i = 0; i < (*sCPU).Nsims; i++) {
        if (sCPU[i].memoryError > 0) {
            printf("Warning: device memory error (%i).\n", sCPU[i].memoryError);
        }
		if (threadBlock[i].joinable()) {
			threadBlock[i].join();
		}
	}
    
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    printf("Finished after %8.4lf s. \n", 1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));


    saveDataSet(sCPU, crystalDatabasePtr, (*sCPU).outputBasePath);
    //free
    free(threadBlock);
    free((*sCPU).sequenceString);
    free((*sCPU).sequenceArray);
    free((*sCPU).refractiveIndex1);
    free((*sCPU).refractiveIndex2);
    free((*sCPU).imdone);
    free((*sCPU).deffTensor);
    free((*sCPU).loadedField1);
    free((*sCPU).loadedField2);
    free((*sCPU).outputBasePath);
    free((*sCPU).field1FilePath);
    free((*sCPU).field2FilePath);
    free(sCPU);
    free(crystalDatabasePtr);
    return 0;
}

unsigned long solveNonlinearWaveEquationSequence(void* lpParam) {
    struct simulationParameterSet* sCPU = (struct simulationParameterSet*)lpParam;
    const double pi = 3.1415926535897932384626433832795;
    int k;
    double rotationAngle;
    for (k = 0; k < (*sCPU).Nsequence; k++) {
        resolveSequence(k, sCPU, (*sCPU).crystalDatabase);
        rotationAngle = (pi / 180) * (*sCPU).sequenceArray[5 + 6 * k];
        solveNonlinearWaveEquation(sCPU);
        if (rotationAngle != 0.0) {
            rotateField(sCPU, rotationAngle);
        }

        if ((*sCPU).memoryError > 0) {
            printf("Warning: device memory error (%i).\n", (*sCPU).memoryError);
        }
    }
    return 0;
}
//main thread of the nonlinear wave equation implemented on CUDA
unsigned long solveNonlinearWaveEquation(void* lpParam) {

    //the struct s contains most of the simulation variables and pointers
    struct cudaParameterSet s;
    struct simulationParameterSet* sCPU = (struct simulationParameterSet*)lpParam;
    cudaSetDevice((*sCPU).assignedGPU);
    cudaStream_t CUDAStream;
    cudaStreamCreate(&CUDAStream);
    s.CUDAStream = CUDAStream;

    //initialize and take values from the struct handed over by the dispatcher
    unsigned long long i;
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
    s.Nblock = (int)(s.Ngrid / THREADS_PER_BLOCK);
    s.isCylindric =(*sCPU).isCylindric;
    s.isNonLinear = ((*sCPU).nonlinearSwitches[0] + (*sCPU).nonlinearSwitches[1] + (*sCPU).nonlinearSwitches[2]) > 0;
    



    //CPU allocations
    std::complex<double>* gridPropagationFactor1CPU = (std::complex<double>*)malloc(2 * s.Ngrid * sizeof(std::complex<double>));
    std::complex<double>* gridPolarizationFactor1CPU = (std::complex<double>*)malloc(2 * s.Ngrid * sizeof(std::complex<double>));
    
    //GPU allocations
    //I shouldn't need all these memsets but, they make me feel better
    int memErrors = 0;
    memErrors += cudaMalloc((void**)&s.gridETime1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridETime1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridETime2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridETime2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridETemp1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridETemp1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridETemp2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridETemp2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridEFrequency1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridEFrequency1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridEFrequency2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridEFrequency2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactor1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridPropagationFactor1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPolarizationFactor1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridPolarizationFactor1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactor2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridPropagationFactor2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPolarizationFactor2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridPolarizationFactor2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactor1Rho1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridPropagationFactor1Rho1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactor1Rho2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridPropagationFactor1Rho2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridRadialLaplacian1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridRadialLaplacian1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridRadialLaplacian2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridRadialLaplacian2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridEFrequency1Next1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridEFrequency1Next1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridEFrequency1Next2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridEFrequency1Next2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.k1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.k1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.k2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.k2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    //the following two should have a size (s.Ntime / 2 + 1) * s.Nspace, but I get overruns during
    //the ffts if they're not larger. If I figure this out, it will save a complex grid worth of memory...
    memErrors += cudaMalloc((void**)&s.gridPolarizationFrequency1, sizeof(cuDoubleComplex) * s.Ngrid); 
    cudaMemset(s.gridPolarizationFrequency1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPolarizationFrequency2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridPolarizationFrequency2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPlasmaCurrentFrequency1, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridPlasmaCurrentFrequency1, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPlasmaCurrentFrequency2, sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridPlasmaCurrentFrequency2, 0, sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPolarizationTime1, sizeof(double) * s.Ngrid);
    cudaMemset(s.gridPolarizationTime1, 0, sizeof(double) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPolarizationTime2, sizeof(double) * s.Ngrid);
    cudaMemset(s.gridPolarizationTime2, 0, sizeof(double) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPlasmaCurrent1, sizeof(double) * s.Ngrid);
    cudaMemset(s.gridPlasmaCurrent1, 0, sizeof(double) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPlasmaCurrent2, sizeof(double) * s.Ngrid);
    cudaMemset(s.gridPlasmaCurrent2, 0, sizeof(double) * s.Ngrid);

    memErrors += cudaMalloc((void**)&s.expGammaT, 2 * sizeof(double) * s.Ntime);
    double* expGammaTCPU = (double*)malloc(2 * sizeof(double) * s.Ntime);
    for (i = 0; i < s.Ntime; i++) {
        expGammaTCPU[i] = exp(s.dt * i * (*sCPU).drudeGamma);
        expGammaTCPU[i + s.Ntime] = exp(-s.dt * i * (*sCPU).drudeGamma);
    }
    cudaMemcpy(s.expGammaT, expGammaTCPU, 2 * sizeof(double) * s.Ntime, cudaMemcpyHostToDevice);
    free(expGammaTCPU);

    memErrors += cudaMalloc((void**)&s.chi2Tensor, sizeof(double) * 9);
    memErrors += cudaMalloc((void**)&s.firstDerivativeOperation, sizeof(double) * 6);
    memErrors += cudaMalloc((void**)&s.chi3Tensor, sizeof(double) * 81);
    memErrors += cudaMalloc((void**)&s.nonlinearSwitches, sizeof(int) * 4);
    memErrors += cudaMalloc((void**)&s.absorptionParameters, sizeof(double) * 6);
    memErrors += cudaMalloc((void**)&s.plasmaParameters, sizeof(double) * 6);
    memErrors += cudaMalloc((void**)&s.propagationInts, sizeof(long long) * 4);
    (*sCPU).memoryError = memErrors;
    if (memErrors > 0) {
        return memErrors;
    }

    //prepare effective nonlinearity tensors and put them on the GPU
    size_t propagationIntsCPU[4] = { s.Ngrid, s.Ntime, s.Nspace, (s.Ntime / 2 + 1) };
    double firstDerivativeOperation[6] = { -1. / 60.,  3. / 20., -3. / 4.,  3. / 4.,  -3. / 20., 1. / 60. };
    for (i = 0; i < 6; i++) {
        firstDerivativeOperation[i] *= (-1.0/(s.Ngrid * s.dx));
    }

    //set nonlinearSwitches[3] to the number of photons needed to overcome bandgap
    (*sCPU).nonlinearSwitches[3] = (int)ceil((*sCPU).bandGapElectronVolts * 241.79893e12 / (*sCPU).frequency1) - 2;

    double plasmaParametersCPU[6] = { 0 };
    
    if ((*sCPU).nonlinearAbsorptionStrength != 0.) {
        s.hasPlasma = TRUE;
    }
    else {
        s.hasPlasma = FALSE;
    }
    
    plasmaParametersCPU[0] = (*sCPU).nonlinearAbsorptionStrength; //nonlinear absorption strength parameter
    plasmaParametersCPU[1] = (*sCPU).drudeGamma; //gamma
    plasmaParametersCPU[2] = (1. / 8.8541878128e-12)*(*sCPU).tStep * (*sCPU).tStep * 2.817832e-08 / (1.6022e-19 * (*sCPU).bandGapElectronVolts * (*sCPU).effectiveMass); // (dt^2)*e* e / (m * band gap));

    calcEffectiveChi2Tensor((*sCPU).deffTensor, (*sCPU).chi2Tensor, (*sCPU).crystalTheta, (*sCPU).crystalPhi);
    cudaMemcpy(s.chi2Tensor, (*sCPU).deffTensor, 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.nonlinearSwitches, (*sCPU).nonlinearSwitches, 4 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(s.propagationInts, propagationIntsCPU, 4 * sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMemcpy(s.chi3Tensor, (*sCPU).chi3Tensor, 27 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.absorptionParameters, (*sCPU).absorptionParameters, 6 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.plasmaParameters, plasmaParametersCPU, 6 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.firstDerivativeOperation, firstDerivativeOperation, 6 * sizeof(double), cudaMemcpyHostToDevice);

    //prepare FFT plans
    cufftPlan2d(&s.fftPlan, (int)s.Nspace, (int)s.Ntime, CUFFT_Z2Z);
    cufftPlan2d(&s.polfftPlan, (int)s.Nspace, (int)s.Ntime, CUFFT_D2Z);
    cufftSetStream(s.fftPlan, CUDAStream);
    cufftSetStream(s.polfftPlan, CUDAStream);

    //prepare the propagation arrays
    if (s.isCylindric) {
        preparePropagation3DCylindric(sCPU, s);
    }
    else {
        preparePropagation2DCartesian(sCPU, s);
    }
    
    //generate the pulses, either through prepareElectricFieldArrays() if this is the first in the series, or by copying
    //the output of the last simulation in the sequence
    if ((*sCPU).isFollowerInSequence) {
        cudaMemcpy(s.gridETime1, (*sCPU).ExtOut, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        cudaMemcpy(s.gridETime2, &(*sCPU).ExtOut[(*sCPU).Ngrid], (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        cudaMemcpy(s.gridEFrequency1, (*sCPU).EkwOut, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        cudaMemcpy(s.gridEFrequency2, &(*sCPU).EkwOut[(*sCPU).Ngrid], (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    }
    else {
        prepareElectricFieldArrays(sCPU, &s);
    }
    
    //Copy the field into the temporary array
    cudaMemcpy(s.gridETemp1, s.gridEFrequency1, s.Nspace * s.Ntime * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    cudaMemcpy(s.gridETemp2, s.gridEFrequency2, s.Nspace * s.Ntime * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);

    //Core propagation loop
    for (i = 0; i < s.Nsteps; i++) {
        //calculate k1
        runRK4Step(s, 0);
        //calculate k2
        runRK4Step(s, 1);
        //calculate k3
        runRK4Step(s, 2);
        //calculate k4
        runRK4Step(s, 3);

        if ((*sCPU).imdone[0] == 2) {
            break;
        }

        if ((*sCPU).imdone[0] == 3) {
            //copy the field arrays from the GPU to CPU memory
            cudaMemcpy((*sCPU).ExtOut, s.gridETime1, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy((*sCPU).EkwOut, s.gridEFrequency1, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy(&(*sCPU).ExtOut[s.Ngrid], s.gridETime2, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy(&(*sCPU).EkwOut[s.Ngrid], s.gridEFrequency2, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
            (*sCPU).imdone[0] = 0;
        }
    }

    
    //transform final result
    fixnanKernel<<<s.Nblock, s.Nthread, 0, CUDAStream>>>(s.gridEFrequency1);
    fixnanKernel << <s.Nblock, s.Nthread, 0, CUDAStream >>> (s.gridEFrequency2);
    cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridEFrequency1, (cufftDoubleComplex*)s.gridETime1, CUFFT_INVERSE);
    cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridEFrequency2, (cufftDoubleComplex*)s.gridETime2, CUFFT_INVERSE);
    fftNormalizeKernel << <s.Nblock, s.Nthread, 0, CUDAStream>>> (s.gridETime1, s.propagationInts);
    fftNormalizeKernel<<<s.Nblock, s.Nthread, 0, CUDAStream >>> (s.gridETime2, s.propagationInts);


    //copy the field arrays from the GPU to CPU memory
    cudaMemcpy((*sCPU).ExtOut, s.gridETime1, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*sCPU).EkwOut, s.gridEFrequency1, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*sCPU).ExtOut[s.Ngrid], s.gridETime2, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*sCPU).EkwOut[s.Ngrid], s.gridEFrequency2, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);


    
    //Free GPU memory
    cudaFree(s.propagationInts);
    cudaFree(s.nonlinearSwitches);
    cudaFree(s.absorptionParameters);
    cudaFree(s.gridETime1); 
    cudaFree(s.gridETemp1);
    cudaFree(s.gridPolarizationFrequency1);
    cudaFree(s.gridEFrequency1);
    cudaFree(s.gridPropagationFactor1);
    cudaFree(s.gridPropagationFactor1Rho1);
    cudaFree(s.gridPropagationFactor1Rho2);
    cudaFree(s.gridRadialLaplacian1);
    cudaFree(s.gridRadialLaplacian2);
    cudaFree(s.firstDerivativeOperation);
    cudaFree(s.gridPolarizationFactor1);
    cudaFree(s.gridEFrequency1Next1);
    cudaFree(s.k1);
    cudaFree(s.gridPolarizationTime1);
    cudaFree(s.gridETime2);
    cudaFree(s.gridETemp2);
    cudaFree(s.gridPolarizationFrequency2);
    cudaFree(s.gridEFrequency2);
    cudaFree(s.gridPropagationFactor2);
    cudaFree(s.gridPolarizationFactor2);
    cudaFree(s.gridEFrequency1Next2);
    cudaFree(s.k2);
    cudaFree(s.gridPolarizationTime2);
    cudaFree(s.chi2Tensor);
    cudaFree(s.chi3Tensor);
    cudaFree(s.expGammaT);
    cufftDestroy(s.fftPlan);
    cufftDestroy(s.polfftPlan);
    cudaFree(s.plasmaParameters);
    cudaFree(s.gridPlasmaCurrent1);
    cudaFree(s.gridPlasmaCurrent2);
    cudaFree(s.gridPlasmaCurrentFrequency1);
    cudaFree(s.gridPlasmaCurrentFrequency2);
    
    cudaStreamDestroy(CUDAStream);

    //Free CPU memory
    free(gridPropagationFactor1CPU);
    free(gridPolarizationFactor1CPU);
    (*sCPU).imdone[0] = 1;
    return 0;
}

//function to run a RK4 time step
//stepNumber is the sub-step index, from 0 to 3
int runRK4Step(struct cudaParameterSet s, int stepNumber) {

    //operations involving FFT
    if (s.isNonLinear || s.isCylindric) {
        //perform inverse FFT to get time-space electric field
        cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridETemp1, (cufftDoubleComplex*)s.gridETime1, CUFFT_INVERSE);
        cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridETemp2, (cufftDoubleComplex*)s.gridETime2, CUFFT_INVERSE);
        
        if (s.isNonLinear) {
            nonlinearPolarizationKernel << <s.Nblock, s.Nthread, 0, s.CUDAStream >> > (s);
            cufftExecD2Z(s.polfftPlan, s.gridPolarizationTime1, (cufftDoubleComplex*)s.gridPolarizationFrequency1);
            cufftExecD2Z(s.polfftPlan, s.gridPolarizationTime2, (cufftDoubleComplex*)s.gridPolarizationFrequency2);
        }

        if (s.isCylindric) {
            radialLaplacianKernel << <s.Nblock, s.Nthread, 0, s.CUDAStream >> > (s);
            cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridRadialLaplacian1, (cufftDoubleComplex*)s.k1, CUFFT_FORWARD);
            cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.gridRadialLaplacian2, (cufftDoubleComplex*)s.k2, CUFFT_FORWARD);
            
        }

        if (s.hasPlasma) {
            plasmaCurrentKernelPrep << <s.Nblock, s.Nthread, 0, s.CUDAStream >> > (s, (double*)s.gridPlasmaCurrentFrequency1, (double*)s.gridPlasmaCurrentFrequency2);
            plasmaCurrentKernel2 << <(unsigned int)s.Nspace, 1, 0, s.CUDAStream >> > (s, (double*)s.gridPlasmaCurrentFrequency1, (double*)s.gridPlasmaCurrentFrequency2);
            //plasmaCurrentKernel << <(unsigned int)s.Nspace, 1 >> > (s);
            cufftExecD2Z(s.polfftPlan, s.gridPlasmaCurrent1, (cufftDoubleComplex*)s.gridPlasmaCurrentFrequency1);
            cufftExecD2Z(s.polfftPlan, s.gridPlasmaCurrent2, (cufftDoubleComplex*)s.gridPlasmaCurrentFrequency2);
            
        }
    }

    //calculate k
    rkKernel<<<s.Nblock, s.Nthread, 0, s.CUDAStream >>>(s, stepNumber);
    
    return 0;
}

int prepareElectricFieldArrays(struct simulationParameterSet* s, struct cudaParameterSet *sc) {
    size_t i,j;
    double rB, zB, r, z; //r and z in the Beam and lab coordinates, respectively.
    double w0, wz, zR, Rz, phi; //Gaussian beam parameters
    double theta = 0; //rotation angle of the current beam
    double pulseSum = 0;
    std::complex<double> ne, no; //active refractive index;
    double f, w; //active frequency;
    double pulseEnergySum;
    std::complex<double> ko, specfac, specphase;
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
        specphase = ii * ((*s).cephase1 + 2*pi*f * (*s).delay1 - (*s).gdd1 * w * w - (*s).tod1 * w * w * w);
        specfac = exp(-specfac - specphase);

        if ((*s).field1IsAllocated) {
            specfac = (*s).loadedField1[i] * exp(-specphase);
        }

        ne = (*s).refractiveIndex1[i + (*s).Ntime * j];
        no = (*s).refractiveIndex2[i + (*s).Ntime * j];
        ko = 2 * pi * no * f / c;
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
            Eb = (w0 / wz) * exp(-ii * (real(ko) * (z-zB) + real(ko) * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
            Eb *= specfac;
            if (isnan(cModulusSquared(Eb)) || f<=0) {
                Eb = 0;
            }
            
            pulse1[i + (*s).Ntime * j] = polFactor1 * Eb;
            pulse1[i + (*s).Ntime * j + (*s).Ngrid] = polFactor2 * Eb;
            pulseSum += abs(r)*(real(ne)*cModulusSquared(pulse1[i + (*s).Ntime * j]) + real(no)*cModulusSquared(pulse1[i + (*s).Ntime * j + (*s).Ngrid]));
        }
    }
    
    // copy the field and propagation grids to the GPU
    cudaMemcpy((*sc).gridETime1, pulse1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridETime2, &pulse1[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    // fft along spatial dimention to get Fourier space beam
    // will take place in three steps:
    // 2D fft (x,f)->(k,t), temporary intermediate state (could be optimized out later)
    // 1D fft (k,t)->(k,f), copied to Fourier space beam
    // 2D fft (k,f)->(x,t), copied to real space beam

    cufftPlan1d(&plan1, (int)(*sc).Ntime, CUFFT_Z2Z, (int)(*sc).Nspace);
    cufftSetStream(plan1, (*sc).CUDAStream);
    cufftPlan2d(&plan2, (int)(*sc).Nspace, (int)(*sc).Ntime, CUFFT_Z2Z);
    cufftSetStream(plan2, (*sc).CUDAStream);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime1, (cufftDoubleComplex*)(*sc).gridETemp1, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).gridETemp1, (cufftDoubleComplex*)(*sc).gridEFrequency1, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridEFrequency1, (cufftDoubleComplex*)(*sc).gridETime1, CUFFT_INVERSE);

    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime2, (cufftDoubleComplex*)(*sc).gridETemp2, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).gridETemp2, (cufftDoubleComplex*)(*sc).gridEFrequency2, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridEFrequency2, (cufftDoubleComplex*)(*sc).gridETime2, CUFFT_INVERSE);

    //Take the conjugate of the field because me and cufft have different ideas of time
    conjugateKernel<<<(*sc).Nblock, (*sc).Nthread, 0, (*sc).CUDAStream >>>((*sc).gridETime1);
    conjugateKernel<<<(*sc).Nblock, (*sc).Nthread, 0, (*sc).CUDAStream >>>((*sc).gridETime2);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime1, (cufftDoubleComplex*)(*sc).gridEFrequency1, CUFFT_INVERSE);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime2, (cufftDoubleComplex*)(*sc).gridEFrequency2, CUFFT_INVERSE);
    cudaDeviceSynchronize();

    //Copy the GPU grids to the CPU memory
    cudaMemcpy(pulse1, (*sc).gridETime1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pulse1[(*s).Ngrid], (*sc).gridETime2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(pulse1f, (*sc).gridEFrequency1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pulse1f[(*s).Ngrid], (*sc).gridEFrequency2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    //normalize the pulse energy and set it to the input value
    pulseSum *= c * eps0;
    pulseSum *= pi; //59.958 is emperical factor
    pulseSum *= (*s).rStep / (*s).fStep;
    pulseEnergySum = sqrt((*s).pulseEnergy1/pulseSum)/(*s).Ngrid;
    (*s).pulse1measEnergy = pulseEnergySum;
    for (i = 0; i < (*s).Ngrid * 2; i++) {
        pulse1[i] = pulse1[i] * pulseEnergySum;
        pulse1f[i] = pulse1f[i] * pulseEnergySum;
    }


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
        specphase = ii * ((*s).cephase2 + 2*pi*f * (*s).delay2 - (*s).gdd2 * w * w - (*s).tod2 * w * w * w);
        specfac = exp(-specfac - specphase);

        if ((*s).field2IsAllocated) {
            specfac = (*s).loadedField2[i] * exp(-specphase);
        }


        ne = (*s).refractiveIndex1[i + (*s).Ntime * j];
        no = (*s).refractiveIndex2[i + (*s).Ntime * j];
        ko = 2 * pi * no * f / c;
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
            if (isnan(cModulusSquared(Eb)) || f <= 0) {
                Eb = 0;
            }

            pulse2[i + (*s).Ntime * j] = polFactor1 * Eb;
            pulse2[i + (*s).Ntime * j + (*s).Ngrid] = polFactor2 * Eb;
            pulseSum += abs(r) * (real(ne) * cModulusSquared(pulse2[i + (*s).Ntime * j]) + real(no) * cModulusSquared(pulse2[i + (*s).Ntime * j + (*s).Ngrid]));
        }
    }

    // copy the field and propagation grids to the GPU
    cudaMemcpy((*sc).gridETime1, pulse2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridETime2, &pulse2[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    // fft along spatial dimention to get Fourier space beam
    // will take place in three steps:
    // 2D fft (x,f)->(k,t), temporary intermediate state (could be optimized out later)
    // 1D fft (k,t)->(k,f), copied to Fourier space beam
    // 2D fft (k,f)->(x,t), copied to real space beam

    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime1, (cufftDoubleComplex*)(*sc).gridETemp1, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).gridETemp1, (cufftDoubleComplex*)(*sc).gridEFrequency1, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridEFrequency1, (cufftDoubleComplex*)(*sc).gridETime1, CUFFT_INVERSE);

    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime2, (cufftDoubleComplex*)(*sc).gridETemp2, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).gridETemp2, (cufftDoubleComplex*)(*sc).gridEFrequency2, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridEFrequency2, (cufftDoubleComplex*)(*sc).gridETime2, CUFFT_INVERSE);

    //Take the conjugate of the field because me and cufft have different ideas of time
    conjugateKernel << <(*sc).Nblock, (*sc).Nthread, 0, (*sc).CUDAStream >> > ((*sc).gridETime1);
    conjugateKernel << <(*sc).Nblock, (*sc).Nthread, 0, (*sc).CUDAStream >> > ((*sc).gridETime2);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime1, (cufftDoubleComplex*)(*sc).gridEFrequency1, CUFFT_INVERSE);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).gridETime2, (cufftDoubleComplex*)(*sc).gridEFrequency2, CUFFT_INVERSE);
    cudaDeviceSynchronize();

    //Copy the GPU grids to the CPU memory
    cudaMemcpy(pulse2, (*sc).gridETime1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pulse2[(*s).Ngrid], (*sc).gridETime2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(pulse2f, (*sc).gridEFrequency1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pulse2f[(*s).Ngrid], (*sc).gridEFrequency2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    //normalize the pulse energy and set it to the input value
    pulseSum *= c * eps0;
    pulseSum *= pi; //59.958 is emperical factor
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
    cudaMemcpy((*sc).gridETime1, (*s).Ext, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridETime2, &(*s).Ext[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridEFrequency1, (*s).Ekw, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).gridEFrequency2, &(*s).Ekw[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cufftDestroy(plan1);
    cufftDestroy(plan2);

    free(pulse1);
    free(pulse2);
    free(pulse1f);
    free(pulse2f);

    return 0;
}

int preparePropagation2DCartesian(struct simulationParameterSet* s, struct cudaParameterSet sc) {
    //recycle allocated device memory for the grids needed
    double* alphaGPU = (double*)sc.gridEFrequency1Next1;
    double* sellmeierCoefficients = (double*)sc.k1;
    sc.ne = sc.gridEFrequency1Next2;
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
    prepareCartesianGridsKernel <<<sc.Nblock, sc.Nthread, 0, sc.CUDAStream >>> (alphaGPU, sellmeierCoefficients, sc);
    cudaDeviceSynchronize();

    //copy the retrieved refractive indicies to the cpu
    cudaMemcpy((*s).refractiveIndex1, sc.ne, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*s).refractiveIndex2, sc.no, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    
    //clean up
    cudaMemset(sc.gridEFrequency1Next1, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.gridEFrequency1Next2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.k1, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.k2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    free(sellmeierCoefficientsAugmentedCPU);
    return 0;
}

int preparePropagation3DCylindric(struct simulationParameterSet* s, struct cudaParameterSet sc) {
    //recycle allocated device memory for the grids needed
    double* sellmeierCoefficients = (double*)sc.k1;
    sc.ne = sc.gridEFrequency1Next2;
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
    prepareCylindricGridsKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sellmeierCoefficients, sc);
    cudaDeviceSynchronize();

    //copy the retrieved refractive indicies to the cpu
    cudaMemcpy((*s).refractiveIndex1, sc.ne, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*s).refractiveIndex2, sc.no, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    //clean up
    cudaMemset(sc.gridEFrequency1, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.gridEFrequency2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.k1, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.k2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    free(sellmeierCoefficientsAugmentedCPU);
    return 0;
}

double findWalkoffAngles(struct simulationParameterSet* s, double dk, double f, double tol) {
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

int calcEffectiveChi2Tensor(double* defftensor, double* dtensor, double theta, double phi) {
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
    double* work = (double*)malloc(128 * sizeof(double));
    int dgelsInfo;
    int dgelsParams[6] = { 3,2,3,3,3,64 };
    dgels("N", &dgelsParams[0], &dgelsParams[1], &dgelsParams[2], R, &dgelsParams[3], dOre, &dgelsParams[4], work, &dgelsParams[5], &dgelsInfo);
    defftensor[0] = dOre[0];
    defftensor[1] = dOre[1];
    defftensor[2] = dOre[3];
    defftensor[3] = dOre[4];
    defftensor[4] = dOre[6];
    defftensor[5] = dOre[7];
    free(work);

    //correct cross-terms
    for (i = 2; i < 4; i++) {
        defftensor[i] *= 0.5;
    }

    for (i = 0; i < 6; i++) {
        defftensor[i] *= 2e-12; //change from pm/V to m/V and multiply by 2 for chi(2) instead of d
    }
    return dgelsInfo;
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
int fftshiftAndFilp(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2) {
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
//REPLACED BY CUDA VERSION, DELETE LATER
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

int loadFrogSpeck(char* frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel, int fieldIndex) {
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
    if (fp == NULL) {
        free(E);
        return -1;
    }
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

int plotDataXY(double* X, double* Y, double minX, double maxX, double minY, double maxY, int N, int plotSizeX, int plotSizeY, double lineWidth, double markerWidth, double* plotGrid, double* xTicks, int NxTicks, double* yTicks, int NyTicks) {
    double normFactorX = plotSizeX/(maxX - minX);
    double normFactorY = plotSizeY/(maxY - minY);
    double offsetY = -minY;
    double offsetX = -minX;

    double* plotGridGPU;
    double* dataX;
    double* dataY;
    double* xTicksGPU;
    double* yTicksGPU;

    cudaMalloc((void**)&plotGridGPU, plotSizeX * plotSizeY * sizeof(double));
    cudaMalloc((void**)&dataX, N * sizeof(double));
    cudaMalloc((void**)&dataY, N * sizeof(double));
    cudaMalloc((void**)&xTicksGPU, NxTicks * sizeof(double));
    cudaMalloc((void**)&yTicksGPU, NyTicks * sizeof(double));
    cudaMemcpy(xTicksGPU, xTicks, NxTicks * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(yTicksGPU, yTicks, NyTicks * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(dataX, X, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dataY, Y, N * sizeof(double), cudaMemcpyHostToDevice);
    plotDataKernel<<<(plotSizeX*plotSizeY)/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(dataX, dataY, lineWidth, markerWidth, N, plotGridGPU, normFactorX, normFactorY, offsetX, offsetY, plotSizeX, plotSizeY, xTicksGPU, NxTicks, yTicksGPU, NyTicks);
    cudaMemcpy(plotGrid, plotGridGPU, (plotSizeX * plotSizeY) * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(plotGridGPU);
    cudaFree(yTicksGPU);
    cudaFree(xTicksGPU);
    cudaFree(dataX);
    cudaFree(dataY);

    return 0;
}

//Rotate the field on the GPU
//Allocates memory and copies from CPU, then copies back to CPU and deallocates
// - inefficient but the general principle is that only the CPU memory is preserved
// after simulations finish... and this only runs at the end of the simulation
int rotateField(struct simulationParameterSet *s, double rotationAngle) {
    cuDoubleComplex* Ein1, * Eout1, * Ein2, * Eout2;
    cudaMalloc((void**)&Ein1,  (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc((void**)&Ein2, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc((void**)&Eout1, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc((void**)&Eout2, (*s).Ngrid * sizeof(cuDoubleComplex));
    unsigned int Nthread = THREADS_PER_BLOCK;
    unsigned int Nblock = (unsigned int)((*s).Ngrid / THREADS_PER_BLOCK);

    cudaMemcpy(Ein1, (*s).EkwOut, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(Ein2, &(*s).EkwOut[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    rotateFieldKernel<<<Nblock, Nthread>>>(Ein1, Ein2, Eout1, Eout2, rotationAngle);

    cudaMemcpy((*s).EkwOut, Eout1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*s).EkwOut[(*s).Ngrid], Eout2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    cudaFree(Ein1);
    cudaFree(Ein2);
    cudaFree(Eout1);
    cudaFree(Eout2);
    return 0;
}

//calculates the squard modulus of a complex number, under the assumption that the
//machine's complex number format is interleaved doubles.
//c forced to run in c++ for nostalgia reasons
double cModulusSquared(std::complex<double>complexNumber) {
    double* xy = (double*)&complexNumber;
    return xy[0] * xy[0] + xy[1] * xy[1];
}

int allocateGrids(struct simulationParameterSet* sCPU) {
    (*sCPU).loadedField1 = (std::complex<double>*)calloc((*sCPU).Ntime, sizeof(std::complex<double>));
    (*sCPU).loadedField2 = (std::complex<double>*)calloc((*sCPU).Ntime, sizeof(std::complex<double>));

    (*sCPU).sequenceArray = (double*)calloc(256 * MAX_LOADSTRING, sizeof(double));

    (*sCPU).Ext = (std::complex<double>*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims, sizeof(std::complex<double>));
    (*sCPU).Ekw = (std::complex<double>*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims, sizeof(std::complex<double>));

    (*sCPU).ExtOut = (std::complex<double>*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims, sizeof(std::complex<double>));
    (*sCPU).EkwOut = (std::complex<double>*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims, sizeof(std::complex<double>));

    (*sCPU).refractiveIndex1 = (std::complex<double>*)calloc((*sCPU).Ngrid * (*sCPU).Nsims, sizeof(std::complex<double>));
    (*sCPU).refractiveIndex2 = (std::complex<double>*)calloc((*sCPU).Ngrid * (*sCPU).Nsims, sizeof(std::complex<double>));
    (*sCPU).deffTensor = (double*)calloc(9 * (*sCPU).Nsims, sizeof(double));
    (*sCPU).imdone = (int*)calloc((*sCPU).Nsims, sizeof(int));
    return 0;
}


int readCrystalDatabase(struct crystalEntry* db) {
    int i = 0;
    double* fd;
    FILE* fp;
    fp = fopen("CrystalDatabase.txt", "r");
    if (fp == NULL) {
        return -2;
    }

    //read the entries line
    int readErrors = 0;

    while (readErrors == 0 && !feof(fp) && i < MAX_LOADSTRING) {
        readErrors += 0 != fwscanf(fp, L"Name:\n");
        fgetws(db[i].crystalNameW, 256, fp);
        readErrors += 1 != fwscanf(fp, L"Type:\n%d\n", &db[i].axisType);
        readErrors += 1 != fwscanf(fp, L"Sellmeier equation:\n%d\n", &db[i].sellmeierType);
        fd = &db[i].sellmeierCoefficients[0];
        readErrors += 22 != fwscanf(fp, L"1st axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
        fd = &db[i].sellmeierCoefficients[22];
        readErrors += 22 != fwscanf(fp, L"2nd axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
        fd = &db[i].sellmeierCoefficients[44];
        readErrors += 22 != fwscanf(fp, L"3rd axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
        readErrors += 0 != fwscanf(fp, L"Sellmeier reference:\n");
        fgetws(db[i].sellmeierReference, 512, fp);
        readErrors += 1 != fwscanf(fp, L"chi2 type:\n%d\n", &db[i].nonlinearSwitches[0]);
        fd = &db[i].d[0];
        readErrors += 18 != fwscanf(fp, L"d:\n%lf %lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[3], &fd[6], &fd[9], &fd[12], &fd[15], &fd[1], &fd[4], &fd[7], &fd[10], &fd[13], &fd[16], &fd[2], &fd[5], &fd[8], &fd[11], &fd[14], &fd[17]);
        readErrors += 0 != fwscanf(fp, L"d reference:\n");
        fgetws(db[i].dReference, 512, fp);
        readErrors += 1 != fwscanf(fp, L"chi3 type:\n%d\n", &db[i].nonlinearSwitches[1]);
        fd = &db[i].chi3[0];
        readErrors += 9 != fwscanf(fp, L"chi3:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8]);
        fd = &db[i].chi3[9];
        readErrors += 9 != fwscanf(fp, L"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8]);
        fd = &db[i].chi3[18];
        readErrors += 9 != fwscanf(fp, L"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8]);
        readErrors += 0 != fwscanf(fp, L"chi3 reference:\n");
        fgetws(db[i].chi3Reference, 512, fp);
        readErrors += 0 != fwscanf(fp, L"Spectral file:\n");
        fgetws(db[i].spectralFile, 512, fp);
        readErrors += 0 != fwscanf(fp, L"~~~crystal end~~~\n");
        if (readErrors == 0) i++;
    }
    db[0].numberOfEntries = i;
    fclose(fp);

    return i;
}

int readSequenceString(struct simulationParameterSet* sCPU) {
    //read the sequence string (if there is one), convert it into an array if it exists
    char sequenceString[MAX_LOADSTRING];
    strcpy(sequenceString, (*sCPU).sequenceString);
    char* tokToken = strtok(sequenceString, ";");
    int sequenceCount = sscanf(sequenceString, "%lf %lf %lf %lf %lf %lf", &(*sCPU).sequenceArray[0], &(*sCPU).sequenceArray[1], &(*sCPU).sequenceArray[2], &(*sCPU).sequenceArray[3], &(*sCPU).sequenceArray[4], &(*sCPU).sequenceArray[5]);

    tokToken = strtok(NULL, ";");
    int lastread = sequenceCount;
    while (tokToken != NULL && lastread == 6) {
        lastread = sscanf(tokToken, "%lf %lf %lf %lf %lf %lf", &(*sCPU).sequenceArray[sequenceCount], &(*sCPU).sequenceArray[sequenceCount + 1], &(*sCPU).sequenceArray[sequenceCount + 2], &(*sCPU).sequenceArray[sequenceCount + 3], &(*sCPU).sequenceArray[sequenceCount + 4], &(*sCPU).sequenceArray[sequenceCount + 5]);
        sequenceCount += lastread;
        tokToken = strtok(NULL, ";");
    }

    (*sCPU).Nsequence = sequenceCount / 6;
    (*sCPU).isInSequence = ((*sCPU).Nsequence > 0);

    if (!(*sCPU).isInSequence) {
        char nopeString[] = "None.";
        strcpy((*sCPU).sequenceString, nopeString);
    }
    return 0;
}

int configureBatchMode(struct simulationParameterSet* sCPU) {
    int j;
    const double pi = 3.1415926535897932384626433832795;

    //Configure the struct array if in a batch
    for (j = 0; j < (*sCPU).Nsims; j++) {
        if (j > 0) {
            memcpy(&sCPU[j], sCPU, sizeof(struct simulationParameterSet));
        }
        
        if ((*sCPU).deffTensor != NULL) {
            sCPU[j].deffTensor = &(*sCPU).deffTensor[9 * j];;
        }

        sCPU[j].Ext = &(*sCPU).Ext[j * (*sCPU).Ngrid * 2];
        sCPU[j].Ekw = &(*sCPU).Ekw[j * (*sCPU).Ngrid * 2];
        sCPU[j].ExtOut = &(*sCPU).ExtOut[j * (*sCPU).Ngrid * 2];
        sCPU[j].EkwOut = &(*sCPU).EkwOut[j * (*sCPU).Ngrid * 2];
        

        sCPU[j].isFollowerInSequence = FALSE;
        
        if ((*sCPU).batchIndex == 1) {
            sCPU[j].delay2 += j * ((-1e-15 * (*sCPU).batchDestination) - ((*sCPU).delay2 - (*sCPU).timeSpan / 2)) / ((*sCPU).Nsims - 1.);
        }
        if ((*sCPU).batchIndex == 2) {
            sCPU[j].pulseEnergy1 += j * ((*sCPU).batchDestination - (*sCPU).pulseEnergy1) / ((*sCPU).Nsims - 1.);
        }
        if ((*sCPU).batchIndex == 3) {
            sCPU[j].cephase1 += j * (pi * (*sCPU).batchDestination - (*sCPU).cephase1) / ((*sCPU).Nsims - 1.);
        }
        if ((*sCPU).batchIndex == 5) {
            sCPU[j].crystalTheta += j * ((pi / 180) * (*sCPU).batchDestination - (*sCPU).crystalTheta) / ((*sCPU).Nsims - 1.);
        }
        if ((*sCPU).batchIndex == 6) {
            sCPU[j].gdd1 += j * (1e-30 * (*sCPU).batchDestination - (*sCPU).gdd1) / ((*sCPU).Nsims - 1.);
        }

        if ((*sCPU).batchIndex == 7) {
            sCPU[j].z01 += j * (1e-6 * (*sCPU).batchDestination - (*sCPU).z01) / ((*sCPU).Nsims - 1.);
        }

        if ((*sCPU).batchIndex == 8) {
            sCPU[j].drudeGamma += j * (1e12 * (*sCPU).batchDestination - (*sCPU).drudeGamma) / ((*sCPU).Nsims - 1.);
        }

        if ((*sCPU).batchIndex == 9) {
            sCPU[j].nonlinearAbsorptionStrength += j * ((*sCPU).batchDestination - (*sCPU).nonlinearAbsorptionStrength) / ((*sCPU).Nsims - 1.);
        }
        if ((*sCPU).batchIndex == 10) {
            sCPU[j].beamwaist1 += j * (1e-6 * (*sCPU).batchDestination - (*sCPU).beamwaist1) / ((*sCPU).Nsims - 1.);
        }
        
    }
    return 0;
}
int readInputParametersFile(struct simulationParameterSet* sCPU, struct crystalEntry* crystalDatabasePtr, char* filePath) {
    FILE* textfile;
    double pi = 3.1415926535897932384626433832795;
    textfile = fopen(filePath, "r");
    if (textfile == NULL) {
        return 1;
    }
    //read parameters using fscanf:
    //recipe for programming: copy/paste the block of fprintf statements in the saveDataSet() function,
    //then find/replace:
    // fprintf->fscanf
    // (*CPU). -> &(*CPU).
    // %e -> %lf
    // &(*sCPU).sequenceString -> (*sCPU).sequenceString
    // &(*sCPU).outputBasePath -> (*sCPU).outputBasePath
    fscanf(textfile, "Pulse energy 1 (J): %lf\nPulse energy 2(J): %lf\nFrequency 1 (Hz): %lf\nFrequency 2 (Hz): %lf\nBandwidth 1 (Hz): %lf\nBandwidth 2 (Hz): %lf\n", &(*sCPU).pulseEnergy1, &(*sCPU).pulseEnergy2, &(*sCPU).frequency1, &(*sCPU).frequency2, &(*sCPU).bandwidth1, &(*sCPU).bandwidth2);
    fscanf(textfile, "SG order: %i\nCEP 1 (rad): %lf\nCEP 2 (rad): %lf\nDelay 1 (s): %lf\nDelay 2 (s): %lf\nGDD 1 (s^-2): %lf\nGDD 2 (s^-2): %lf\nTOD 1 (s^-3): %lf\nTOD 2(s^-3): %lf\n", &(*sCPU).sgOrder1, &(*sCPU).cephase1, &(*sCPU).cephase2, &(*sCPU).delay1, &(*sCPU).delay2, &(*sCPU).gdd1, &(*sCPU).gdd2, &(*sCPU).tod1, &(*sCPU).tod2);
    fscanf(textfile, "Beamwaist 1 (m): %lf\nBeamwaist 2 (m): %lf\nx offset 1 (m): %lf\nx offset 2 (m): %lf\nz offset 1 (m): %lf\nz offset 2 (m): %lf\nNC angle 1 (rad): %lf\nNC angle 2 (rad): %lf\n", &(*sCPU).beamwaist1, &(*sCPU).beamwaist2, &(*sCPU).x01, &(*sCPU).x02, &(*sCPU).z01, &(*sCPU).z02, &(*sCPU).propagationAngle1, &(*sCPU).propagationAngle2);
    fscanf(textfile, "Polarization 1 (rad): %lf\nPolarization 2 (rad): %lf\nCircularity 1: %lf\nCircularity 2: %lf\n", &(*sCPU).polarizationAngle1, &(*sCPU).polarizationAngle2, &(*sCPU).circularity1, &(*sCPU).circularity2);
    fscanf(textfile, "Material index: %i\n", &(*sCPU).materialIndex);
    fscanf(textfile, "Crystal theta (rad): %lf\nCrystal phi (rad): %lf\nGrid width (m): %lf\ndx (m): %lf\nTime span (s): %lf\ndt (s): %lf\nThickness (m): %lf\ndz (m): %lf\n", &(*sCPU).crystalTheta, &(*sCPU).crystalPhi, &(*sCPU).spatialWidth, &(*sCPU).rStep, &(*sCPU).timeSpan, &(*sCPU).tStep, &(*sCPU).crystalThickness, &(*sCPU).propagationStep);
    fscanf(textfile, "Nonlinear absorption parameter: %lf\nBand gap (eV): %lf\nEffective mass (relative): %lf\nDrude gamma (Hz): %lf\n", &(*sCPU).nonlinearAbsorptionStrength, &(*sCPU).bandGapElectronVolts, &(*sCPU).effectiveMass, &(*sCPU).drudeGamma);
    fscanf(textfile, "Propagation mode: %i\n", &(*sCPU).symmetryType);
    fscanf(textfile, "Batch mode: %i\nBatch destination: %lf\nBatch steps: %lli\n", &(*sCPU).batchIndex, &(*sCPU).batchDestination, &(*sCPU).Nsims);
    fscanf(textfile, "Sequence: %s\n", (*sCPU).sequenceString);
    fscanf(textfile, "Output base path: %s\n", (*sCPU).outputBasePath);
    fscanf(textfile, "Field 1 from file type: %i\nField 2 from file type: %i\n", &(*sCPU).pulse1FileType, &(*sCPU).pulse2FileType);
    fscanf(textfile, "Field 1 file path: %s\n", (*sCPU).field1FilePath);
    fscanf(textfile, "Field 2 file path: %s\n", (*sCPU).field2FilePath);

    //derived parameters and cleanup:
    (*sCPU).sellmeierType = 0;
    (*sCPU).axesNumber = 0;
    (*sCPU).sgOrder2 = (*sCPU).sgOrder1;
    (*sCPU).Ntime = (size_t)round((*sCPU).timeSpan / (*sCPU).tStep);
    (*sCPU).Nspace = (size_t)round((*sCPU).spatialWidth / (*sCPU).rStep);
    (*sCPU).Ngrid = (*sCPU).Ntime * (*sCPU).Nspace;
    (*sCPU).kStep = 2 * pi / ((*sCPU).Nspace * (*sCPU).rStep);
    (*sCPU).fStep = 1.0 / ((*sCPU).Ntime * (*sCPU).tStep);
    (*sCPU).Npropagation = (size_t)round((*sCPU).crystalThickness / (*sCPU).propagationStep);

    (*sCPU).isCylindric = (*sCPU).symmetryType == 1;
    if ((*sCPU).isCylindric) {
        (*sCPU).x01 = 0;
        (*sCPU).x02 = 0;
        (*sCPU).propagationAngle1 = 0;
        (*sCPU).propagationAngle2 = 0;
    }

    if ((*sCPU).batchIndex == 0 || (*sCPU).batchIndex == 4 || (*sCPU).Nsims < 1) {
        (*sCPU).Nsims = 1;
    }

    (*sCPU).field1IsAllocated = FALSE;
    (*sCPU).field2IsAllocated = FALSE;

    //crystal from database (database must be loaded!)
    (*sCPU).chi2Tensor = crystalDatabasePtr[(*sCPU).materialIndex].d;
    (*sCPU).chi3Tensor = crystalDatabasePtr[(*sCPU).materialIndex].chi3;
    (*sCPU).nonlinearSwitches = crystalDatabasePtr[(*sCPU).materialIndex].nonlinearSwitches;
    (*sCPU).absorptionParameters = crystalDatabasePtr[(*sCPU).materialIndex].absorptionParameters;
    (*sCPU).sellmeierCoefficients = crystalDatabasePtr[(*sCPU).materialIndex].sellmeierCoefficients;
    (*sCPU).sellmeierType = crystalDatabasePtr[(*sCPU).materialIndex].sellmeierType;
    (*sCPU).axesNumber = crystalDatabasePtr[(*sCPU).materialIndex].axisType;


    fclose(textfile);
    return 0;
}



int saveSlurmScript(struct simulationParameterSet* sCPU, int gpuType, int gpuCount) {
    FILE* textfile;
    char* stringConversionBuffer = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    wchar_t* wideStringConversionBuffer = (wchar_t*)calloc(MAX_LOADSTRING, sizeof(char));
    char* outputpath = (char*)calloc(MAX_LOADSTRING, sizeof(char));

    char* fileName = (*sCPU).outputBasePath;
    while (strchr(fileName, '\\') != NULL) {
        fileName = strchr(fileName, '\\');
        fileName++;
    }
    char LF = '\x0A';
    strcpy(outputpath, (*sCPU).outputBasePath);
    strcat(outputpath, ".slurmScript");
    textfile = fopen(outputpath, "wb");
    fprintf(textfile, "#!/bin/bash -l"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "#SBATCH -o ./tjob.out.%%j"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "#SBATCH -e ./tjob.err.%%j"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "#SBATCH -D ./"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "#SBATCH -J lightwave"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "#SBATCH --constraint=\"gpu\""); fwrite(&LF, sizeof(char), 1, textfile);
    if (gpuType == 0) {
        fprintf(textfile, "#SBATCH --gres=gpu:rtx5000:%i", min(gpuCount,2)); fwrite(&LF, sizeof(char), 1, textfile);
    }
    if (gpuType == 1) {
        fprintf(textfile, "#SBATCH --gres=gpu:v100:%i", min(gpuCount, 2)); fwrite(&LF, sizeof(char), 1, textfile);
    }
    if (gpuType == 2) {
        fprintf(textfile, "#SBATCH --gres=gpu:a100:%i", min(gpuCount, 4)); fwrite(&LF, sizeof(char), 1, textfile);
        fprintf(textfile, "#SBATCH --cpus-per-task=%i", 2*min(gpuCount, 4)); fwrite(&LF, sizeof(char), 1, textfile);
    }
    fprintf(textfile, "#SBATCH --mem=%lliM",1024+(18 * sizeof(double) * (*sCPU).Ngrid * max(1,(*sCPU).Nsims))/1048576);
    fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "#SBATCH --nodes=1"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "#SBATCH --ntasks-per-node=1"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "#SBATCH --time=24:00:00"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "module purge"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "module load cuda/11.2"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "module load mkl/2022.0"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "module load gcc/9"); fwrite(&LF, sizeof(char), 1, textfile);
    fprintf(textfile, "export LD_LIBRARY_PATH=$MKL_HOME/lib/intel64:$LD_LIBRARY_PATH"); fwrite(&LF, sizeof(char), 1, textfile);
    if (gpuType == 0) {
        fprintf(textfile, "srun ./nnp75 %s.input > prog.out", fileName); fwrite(&LF, sizeof(char), 1, textfile);
    }
    if (gpuType == 1) {
        fprintf(textfile, "srun ./nnp70 %s.input > prog.out", fileName); fwrite(&LF, sizeof(char), 1, textfile);
    }
    if (gpuType == 2) {
        fprintf(textfile, "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}"); fwrite(&LF, sizeof(char), 1, textfile);
        fprintf(textfile, "srun ./nnp80 %s.input > $prog.out", fileName); fwrite(&LF, sizeof(char), 1, textfile);
    }
    fclose(textfile);
    free(outputpath);
    free(wideStringConversionBuffer);
    free(stringConversionBuffer);
    return 0;
}

int saveSettingsFile(struct simulationParameterSet* sCPU, struct crystalEntry* crystalDatabasePtr) {
    int j, k;
    FILE* textfile;
    char* stringConversionBuffer = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    wchar_t* wideStringConversionBuffer = (wchar_t*)calloc(MAX_LOADSTRING, sizeof(char));
    char* outputpath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    strcpy(outputpath, (*sCPU).outputBasePath);
    if ((*sCPU).runType > 0) {
        strcat(outputpath, ".input");
    }
    else {
        strcat(outputpath, ".txt");
    }
    
    textfile = fopen(outputpath, "w");
    fwprintf(textfile, L"Pulse energy 1 (J): %e\nPulse energy 2(J): %e\nFrequency 1 (Hz): %e\nFrequency 2 (Hz): %e\nBandwidth 1 (Hz): %e\nBandwidth 2 (Hz): %e\n", (*sCPU).pulseEnergy1, (*sCPU).pulseEnergy2, (*sCPU).frequency1, (*sCPU).frequency2, (*sCPU).bandwidth1, (*sCPU).bandwidth2);
    fwprintf(textfile, L"SG order: %i\nCEP 1 (rad): %e\nCEP 2 (rad): %e\nDelay 1 (s): %e\nDelay 2 (s): %e\nGDD 1 (s^-2): %e\nGDD 2 (s^-2): %e\nTOD 1 (s^-3): %e\nTOD 2(s^-3): %e\n", (*sCPU).sgOrder1, (*sCPU).cephase1, (*sCPU).cephase2, (*sCPU).delay1, (*sCPU).delay2, (*sCPU).gdd1, (*sCPU).gdd2, (*sCPU).tod1, (*sCPU).tod2);
    fwprintf(textfile, L"Beamwaist 1 (m): %e\nBeamwaist 2 (m): %e\nx offset 1 (m): %e\nx offset 2 (m): %e\nz offset 1 (m): %e\nz offset 2 (m): %e\nNC angle 1 (rad): %e\nNC angle 2 (rad): %e\n", (*sCPU).beamwaist1, (*sCPU).beamwaist2, (*sCPU).x01, (*sCPU).x02, (*sCPU).z01, (*sCPU).z02, (*sCPU).propagationAngle1, (*sCPU).propagationAngle2);
    fwprintf(textfile, L"Polarization 1 (rad): %e\nPolarization 2 (rad): %e\nCircularity 1: %e\nCircularity 2: %e\n", (*sCPU).polarizationAngle1, (*sCPU).polarizationAngle2, (*sCPU).circularity1, (*sCPU).circularity2);
    fwprintf(textfile, L"Material index: %i\n", (*sCPU).materialIndex);
    fwprintf(textfile, L"Crystal theta (rad): %e\nCrystal phi (rad): %e\nGrid width (m): %e\ndx (m): %e\nTime span (s): %e\ndt (s): %e\nThickness (m): %e\ndz (m): %e\n", (*sCPU).crystalTheta, (*sCPU).crystalPhi, (*sCPU).spatialWidth, (*sCPU).rStep, (*sCPU).timeSpan, (*sCPU).tStep, (*sCPU).crystalThickness, (*sCPU).propagationStep);
    fwprintf(textfile, L"Nonlinear absorption parameter: %e\nBand gap (eV): %e\nEffective mass (relative): %e\nDrude gamma (Hz): %e\n", (*sCPU).nonlinearAbsorptionStrength, (*sCPU).bandGapElectronVolts, (*sCPU).effectiveMass, (*sCPU).drudeGamma);
    fwprintf(textfile, L"Propagation mode: %i\n", (*sCPU).symmetryType);
    fwprintf(textfile, L"Batch mode: %i\nBatch destination: %e\nBatch steps: %lli\n", (*sCPU).batchIndex, (*sCPU).batchDestination, (*sCPU).Nsims);
    mbstowcs(wideStringConversionBuffer, (*sCPU).sequenceString, MAX_LOADSTRING);
    fwprintf(textfile, L"Sequence: %ls\n", wideStringConversionBuffer);


    if ((*sCPU).runType > 0) {
        char* fileName = (*sCPU).outputBasePath;
        while (strchr(fileName, '\\') != NULL) {
            fileName = strchr(fileName, '\\');
            fileName++;
        }
        mbstowcs(wideStringConversionBuffer, fileName, strlen(fileName));
        fwprintf(textfile, L"Output base path: %ls\n", wideStringConversionBuffer);
    }
    else {
        mbstowcs(wideStringConversionBuffer, (*sCPU).outputBasePath, MAX_LOADSTRING);
        fwprintf(textfile, L"Output base path: %ls\n", wideStringConversionBuffer);
    }

    fwprintf(textfile, L"Field 1 from file type: %i\nField 2 from file type: %i\n", (*sCPU).pulse1FileType, (*sCPU).pulse2FileType);
    mbstowcs(wideStringConversionBuffer, (*sCPU).field1FilePath, MAX_LOADSTRING);
    fwprintf(textfile, L"Field 1 file path: %ls\n", wideStringConversionBuffer);
    mbstowcs(wideStringConversionBuffer, (*sCPU).field2FilePath, MAX_LOADSTRING);
    fwprintf(textfile, L"Field 2 file path: %ls\n", wideStringConversionBuffer);

    fwprintf(textfile, L"Material name: %ls\nSellmeier reference: %ls\nChi2 reference: %ls\nChi3 reference: %ls\n", crystalDatabasePtr[(*sCPU).materialIndex].crystalNameW, crystalDatabasePtr[(*sCPU).materialIndex].sellmeierReference, crystalDatabasePtr[(*sCPU).materialIndex].dReference, crystalDatabasePtr[(*sCPU).materialIndex].chi3Reference);
    fwprintf(textfile, L"Sellmeier coefficients: \n");
    for (j = 0; j < 3; j++) {
        for (k = 0; k < 22; k++) {
            fwprintf(textfile, L"%e ", crystalDatabasePtr[(*sCPU).materialIndex].sellmeierCoefficients[j * 22 + k]);
        }
        fwprintf(textfile, L"\n");
    }
    fwprintf(textfile, L"Code version: 0.12 April 7, 2022\n");

    fclose(textfile);
    free(outputpath);
    free(wideStringConversionBuffer);
    free(stringConversionBuffer);
    return 0;
}

int saveDataSet(struct simulationParameterSet* sCPU, struct crystalEntry* crystalDatabasePtr, char* outputbase) {
    int j;

    saveSettingsFile(sCPU, crystalDatabasePtr);

    //Save the results as double instead of complex
    double* saveEout = (double*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims, sizeof(double));
    for (j = 0; j < ((*sCPU).Ngrid * (*sCPU).Nsims * 2); j++) {
        saveEout[j] = real((*sCPU).ExtOut[j]);
    }

    char* stringConversionBuffer = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    wchar_t* wideStringConversionBuffer = (wchar_t*)calloc(MAX_LOADSTRING, sizeof(char));
    char* outputpath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    char* outputbaseVar = strrchr(outputbase, '\\');
    if (!outputbaseVar) {
        outputbaseVar = outputbase;
    }
    else {
        outputbaseVar++;
    }
    double* matlabpadding = (double*)calloc(1024, sizeof(double));

    
    
    
    //write fields as binary
    for (j = 0; j < ((*sCPU).Ngrid * (*sCPU).Nsims * 2); j++) {
        saveEout[j] = real((*sCPU).ExtOut[j]);
    }
    FILE* ExtOutFile;
    size_t writeSize = 2 * ((*sCPU).Ngrid * (*sCPU).Nsims);
    strcpy(outputpath, outputbase);
    strcat(outputpath, "_ExtOut.dat");
    ExtOutFile = fopen(outputpath, "wb");
    fwrite(saveEout, sizeof(double), writeSize, ExtOutFile);
    fwrite(matlabpadding, sizeof(double), 1024, ExtOutFile);
    fclose(ExtOutFile);

    for (j = 0; j < ((*sCPU).Ngrid * (*sCPU).Nsims * 2); j++) {
        saveEout[j] = real((*sCPU).Ext[j]);
    }
    FILE* ExtInFile;
    strcpy(outputpath, outputbase);
    strcat(outputpath, "_ExtIn.dat");
    ExtInFile = fopen(outputpath, "wb");
    fwrite(saveEout, sizeof(double), writeSize, ExtInFile);
    fwrite(matlabpadding, sizeof(double), 1024, ExtInFile);
    fclose(ExtInFile);


    FILE* matlabfile;
    strcpy(outputpath, outputbase);
    strcat(outputpath, ".m");
    matlabfile = fopen(outputpath, "w");
    fprintf(matlabfile, "fid = fopen('%s_ExtIn.dat','rb'); \n", outputbaseVar);
    fprintf(matlabfile, "%s_ExtIn = fread(fid, %lli, 'double'); \n", outputbaseVar, 2 * (*sCPU).Ngrid * (*sCPU).Nsims);
    fprintf(matlabfile, "%s_ExtIn = reshape(%s_ExtIn,[%lli %lli %lli]); \n", outputbaseVar, outputbaseVar, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims);
    fprintf(matlabfile, "fclose(fid); \n");

    fprintf(matlabfile, "fid = fopen('%s_ExtOut.dat','rb'); \n", outputbaseVar);
    fprintf(matlabfile, "%s_ExtOut = fread(fid, %lli, 'double'); \n", outputbaseVar, 2 * (*sCPU).Ngrid * (*sCPU).Nsims);
    fprintf(matlabfile, "%s_ExtOut = reshape(%s_ExtOut,[%lli %lli %lli]); \n", outputbaseVar, outputbaseVar, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims);
    fprintf(matlabfile, "fclose(fid); \n");
    fclose(matlabfile);
    
    //write a python script for loading the output fields in a proper shape
    char scriptfilename[MAX_LOADSTRING];
    strcpy(scriptfilename, outputbase);
    strcat(scriptfilename, ".py");
    FILE* scriptfile;
    scriptfile = fopen(scriptfilename, "w");
    fprintf(scriptfile, "#!/usr/bin/python\nimport numpy as np\n");
    fprintf(scriptfile, "dt = %e\ndz = %e\ndx = %e\n", (*sCPU).tStep, (*sCPU).propagationStep, (*sCPU).rStep);
    fprintf(scriptfile, "%s_ExtIn = np.reshape(np.fromfile(\"", outputbaseVar);
    fprintf(scriptfile, "%s_ExtIn.dat", outputbaseVar);
    fprintf(scriptfile, "\",dtype=np.double)[0:%lli],(%lli,%lli,%lli),order='F')\n", 2 * (*sCPU).Ngrid * (*sCPU).Nsims, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims);
    fprintf(scriptfile, "%s_ExtOut = np.reshape(np.fromfile(\"", outputbaseVar);
    fprintf(scriptfile, "%s_ExtOut.dat", outputbaseVar);
    fprintf(scriptfile, "\",dtype=np.double)[0:%lli],(%lli,%lli,%lli),order='F')\n", 2 * (*sCPU).Ngrid * (*sCPU).Nsims, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims);
    fclose(scriptfile);
    
    free(saveEout);
    free(matlabpadding);
    free(stringConversionBuffer);
    free(wideStringConversionBuffer);
    return 0;
}

int resolveSequence(int currentIndex, struct simulationParameterSet* s, struct crystalEntry* db) {
    double pi = 3.1415926535897932384626433832795;

    //sequence format
    //material index, theta, phi, crystal length, propagation step, rotation angle
    int materialIndex = (int)(*s).sequenceArray[0 + 6 * currentIndex];
    double crystalTheta = (pi / 180) * (*s).sequenceArray[1 + 6 * currentIndex];
    double crystalPhi = (pi / 180) * (*s).sequenceArray[2 + 6 * currentIndex];
    double propagationStep = 1e-9 * (*s).sequenceArray[4 + 6 * currentIndex];
    size_t Npropagation = (size_t)(1e-6 * (*s).sequenceArray[3 + 6 * currentIndex] / propagationStep);
    double rotationAngle = (pi / 180) * (*s).sequenceArray[5 + 6 * currentIndex];

    if (currentIndex > 0) {
        (*s).isFollowerInSequence = TRUE;
    }

    (*s).propagationStep = propagationStep;
    (*s).Npropagation = Npropagation;


    (*s).materialIndex = materialIndex;
    (*s).crystalTheta = crystalTheta;
    (*s).crystalPhi = crystalPhi;
    (*s).chi2Tensor = db[materialIndex].d;
    (*s).chi3Tensor = db[materialIndex].chi3;
    (*s).nonlinearSwitches = db[materialIndex].nonlinearSwitches;
    (*s).absorptionParameters = db[materialIndex].absorptionParameters;
    (*s).sellmeierCoefficients = db[materialIndex].sellmeierCoefficients;

    (*s).sellmeierType = db[materialIndex].sellmeierType;
    (*s).axesNumber = db[materialIndex].axisType;
    return 0;
}

int loadPulseFiles(struct simulationParameterSet* sCPU) {

    //pulse type specifies if something has to be loaded to describe the pulses, or if they should be
    //synthesized later. 1: FROG .speck format; 2: EOS (not implemented yet)
    int frogLines = 0;
    if ((*sCPU).pulse1FileType == 1) {
        frogLines = loadFrogSpeck((*sCPU).field1FilePath, (*sCPU).loadedField1, (*sCPU).Ntime, (*sCPU).fStep, 0.0, 1);
        if (frogLines > 0) {
            (*sCPU).field1IsAllocated = TRUE;
        }
        else {
            return 1;
        }
    }

    if ((*sCPU).pulse2FileType == 1) {
        frogLines = loadFrogSpeck((*sCPU).field2FilePath, (*sCPU).loadedField2, (*sCPU).Ntime, (*sCPU).fStep, 0.0, 1);
        if (frogLines > 0) {
            (*sCPU).field2IsAllocated = TRUE;
        }
        else {
            return 1;
        }
    }
    return 0;
}

int loadSavedFields(struct simulationParameterSet* sCPU, char* outputBase) {
    char* outputpath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    size_t writeSize = 2 * ((*sCPU).Ngrid * (*sCPU).Nsims);
    double* loadE = (double*)malloc(writeSize * sizeof(double));
    size_t j;

    //read fields as binary
    FILE* ExtOutFile;
    
    strcpy(outputpath, outputBase);
    strcat(outputpath, "_ExtOut.dat");
    ExtOutFile = fopen(outputpath, "rb");
    if (ExtOutFile == NULL) {
        return 1;
    }
    fread(loadE, sizeof(double), writeSize, ExtOutFile);
    fclose(ExtOutFile);
    for (j = 0; j < writeSize; j++) {
        (*sCPU).ExtOut[j] = loadE[j];
    }

    cufftHandle fftPlan;
    cufftPlan2d(&fftPlan, (int)(*sCPU).Nspace, (int)(*sCPU).Ntime, CUFFT_Z2Z);

    cuDoubleComplex* fieldGridkw;
    cuDoubleComplex* fieldGridxt;
    cudaMalloc((void**)&fieldGridkw, sizeof(cuDoubleComplex) * (*sCPU).Ngrid);
    cudaMalloc((void**)&fieldGridxt, sizeof(cuDoubleComplex) * (*sCPU).Ngrid);

    for (j = 0; j < 2*(*sCPU).Nsims; j++) {
        cudaMemcpy(fieldGridxt, &(*sCPU).ExtOut[j*(*sCPU).Ngrid], (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        cufftExecZ2Z(fftPlan, fieldGridxt, fieldGridkw, CUFFT_FORWARD);
        cudaMemcpy(&(*sCPU).EkwOut[j*(*sCPU).Ngrid], fieldGridkw, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    }


    cudaFree(fieldGridkw);
    cudaFree(fieldGridxt);
    cufftDestroy(fftPlan);
    free(loadE);
    free(outputpath);
    return 0;
}