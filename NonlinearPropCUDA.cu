#include "cuda_runtime.h"
#include "device_launch_parameters.h"
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

//fitting parameter set as global variable
simulationParameterSet* fittingSet;
simulationParameterSet* fittingReferenceSet;

#define THREADS_PER_BLOCK 32
#define FALSE 0
#define TRUE 1
#define MAX_LOADSTRING 1024
#define TWOPI 6.2831853071795862319959269370884
#define PI 3.14159265358979323846264338327950288
#define DEG2RAD 0.017453292519943295474371680597869
#define LIGHTC 2.99792458e8
#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

//overload the math operators for cuda complex numbers so this code fits inside the observable universe
__device__ __forceinline__ cuDoubleComplex operator*(cuDoubleComplex a, cuDoubleComplex b) { return cuCmul(a, b); }
__device__ __forceinline__ cuDoubleComplex operator*(cuDoubleComplex a, double b) { return make_cuDoubleComplex(a.x * b, a.y * b); }
__device__ __forceinline__ cuDoubleComplex operator*(double b, cuDoubleComplex a) { return make_cuDoubleComplex(a.x * b, a.y * b); }
__device__ __forceinline__ cuDoubleComplex operator+(cuDoubleComplex a, cuDoubleComplex b) { return cuCadd(a, b); }
__device__ __forceinline__ cuDoubleComplex operator+(double a, cuDoubleComplex b) { return make_cuDoubleComplex(b.x + a, b.y); }
__device__ __forceinline__ cuDoubleComplex operator+(cuDoubleComplex a, double b) { return make_cuDoubleComplex(a.x + b, a.y); }
__device__ __forceinline__ cuDoubleComplex operator-(cuDoubleComplex a, cuDoubleComplex b) { return cuCsub(a, b); }
__device__ __forceinline__ cuDoubleComplex operator-(double a, cuDoubleComplex b) { return make_cuDoubleComplex(a-b.x, -b.y); }
__device__ __forceinline__ cuDoubleComplex operator-(cuDoubleComplex a, double b) { return make_cuDoubleComplex(a.x-b,a.y); }
__device__ __forceinline__ cuDoubleComplex operator/(cuDoubleComplex b, cuDoubleComplex a) { return cuCdiv(b, a); }
__device__ __forceinline__ cuDoubleComplex operator/(cuDoubleComplex a, double b) { return make_cuDoubleComplex(a.x / b, a.y / b); }
__device__  cuDoubleComplex operator/(double b, cuDoubleComplex a) {
    double divbByDenominator = b/(a.x * a.x + a.y * a.y);
    return make_cuDoubleComplex(a.x * divbByDenominator, -a.y * divbByDenominator);
}

//complex exponential function for CUDA
__device__ cuDoubleComplex cuCexpd(cuDoubleComplex z){
    double r = exp(z.x);
    return make_cuDoubleComplex(r*cos(z.y),r*sin(z.y));
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
__device__ cuDoubleComplex sellmeierSubfunctionCuda(
    double* a, double ls, double omega, cuDoubleComplex ii, double kL) {
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
__device__ cuDoubleComplex sellmeierCuda(
    cuDoubleComplex* ne, cuDoubleComplex* no, double* a, double f, double theta, double phi, int type, int eqn) {
    if (f == 0) return make_cuDoubleComplex(1.0, 0.0); //exit immediately for f=0

    double ls = 2.99792458e14 / f; //wavelength in microns
    ls *= ls; //only wavelength^2 is ever used
    cuDoubleComplex ii = make_cuDoubleComplex(0.0, 1.0);
    double omega = TWOPI * abs(f);
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
                + ((sinPhi2 / realna2 + cosPhi2 / realnb2)
                    * cosTheta2 + sinTheta2 / (cuCreal(nc) * cuCreal(nc)))));

        ne[0] = 1.0 / cuCsqrt(cos(delta) * cos(delta) * (cosTheta2 * (cosPhi2 / (na * na)
            + sinPhi2 / (nb * nb)) + sinTheta2 / (nc * nc))
            + sin(delta) * sin(delta) * (sinPhi2 / (na * na) + cosPhi2 / (nb * nb))
            - 0.5 * sin(2 * phi) * cosTheta * sin(2 * delta) * (1. / (na * na) - 1. / (nb * nb)));

        no[0] = 1.0 / cuCsqrt(sin(delta) * sin(delta) * (cosTheta2 * (cosPhi2 / (na * na)
            + sinPhi2 / (nb * nb)) + sinTheta2 / (nc * nc))
            + cos(delta) * cos(delta) * (sinPhi2 / (na * na) + cosPhi2 / (nb * nb))
            + 0.5 * sin(2 * phi) * cosTheta * sin(2 * delta) * (1. / (na * na) - 1. / (nb * nb)));
        return ne[0];
    }
}
__global__ void millersRuleNormalizationKernel(cudaParameterSet s, double* sellmeierCoefficients, double* referenceFrequencies) {
    if (!s.isUsingMillersRule) {
        return;
    }
    size_t i;
    double chi11[7];
    double chi12[7];
    cuDoubleComplex ne, no;
    for (i = 0; i < 7; i++) {
        if (referenceFrequencies[i] == 0) {
            chi11[i] = 100000.0;
            chi12[i] = 100000.0;
        }
        else {
            sellmeierCuda(&ne, &no, sellmeierCoefficients, referenceFrequencies[i], sellmeierCoefficients[66], sellmeierCoefficients[67], (int)sellmeierCoefficients[69], 0);
            chi11[i] = cuCreal(ne) * cuCreal(ne) - 1;
            chi12[i] = cuCreal(no) * cuCreal(no) - 1;
        }

    }

    //normalize chi2 tensor values
    s.chi2Tensor[0] /= chi11[0] * chi11[1] * chi11[2];
    s.chi2Tensor[1] /= chi11[0] * chi11[1] * chi12[2];
    s.chi2Tensor[2] /= chi11[0] * chi12[1] * chi11[2];
    s.chi2Tensor[3] /= chi11[0] * chi12[1] * chi12[2];
    s.chi2Tensor[4] /= chi12[0] * chi12[1] * chi11[2];
    s.chi2Tensor[5] /= chi12[0] * chi12[1] * chi12[2];

    //normalize chi3 tensor values
    // note that currently full chi3 isn't implemented so
    // this only applies to the first element, chi3_1111 under
    // the assumption of centrosymmetry
    s.chi3Tensor[0] /= chi11[3] * chi11[4] * chi11[5] * chi11[6];

}

__device__ __forceinline__ double cuCModSquared(cuDoubleComplex a) {
    return a.x * a.x + a.y * a.y;
}
__global__ void totalSpectrumKernel(cuDoubleComplex* fieldGrid1, cuDoubleComplex* fieldGrid2, double gridStep, size_t Ntime, size_t Nspace, double* spectrum) {
    size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    size_t j;
    double beamCenter1 = 0.;
    double beamCenter2 = 0.;
    double beamTotal1 = 0.;
    double beamTotal2 = 0.;
    double a, x;

    //find beam centers
    for (j = 0; j < Nspace; j++) {
        x = gridStep * j;
        a = cuCModSquared(fieldGrid1[i + j * Ntime]);
        beamTotal1 += a;
        beamCenter1 += x * a;
        a = cuCModSquared(fieldGrid2[i + j * Ntime]);
        beamTotal2 += a;
        beamCenter2 += x * a;
    }
    if (beamTotal1 > 0) {
        beamCenter1 /= beamTotal1;
    }
    if (beamTotal2 > 0) {
        beamCenter2 /= beamTotal2;
    }
    

    //Integrate total beam power, assuming radially-symmetric beam around
    //the center
    beamTotal1 = 0.;
    beamTotal2 = 0.;
    for (j = 0; j < Nspace; j++) {
        x = gridStep * j;
        beamTotal1 += abs(x - beamCenter1) * cuCModSquared(fieldGrid1[i + j * Ntime]);
        beamTotal2 += abs(x - beamCenter2) * cuCModSquared(fieldGrid2[i + j * Ntime]);
    }
    beamTotal1 *= gridStep/Ntime;
    beamTotal2 *= gridStep/Ntime;

    //put the values into the output spectrum
    spectrum[i] = beamTotal1;
    spectrum[i + Ntime] = beamTotal2;
    spectrum[i + 2 * Ntime] = beamTotal1 + beamTotal2;
}

//rotate the field around the propagation axis (basis change)
__global__ void rotateFieldKernel(
    cuDoubleComplex* Ein1, cuDoubleComplex* Ein2, cuDoubleComplex* Eout1, 
    cuDoubleComplex* Eout2, double rotationAngle) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    Eout1[i] = cos(rotationAngle) * Ein1[i] - sin(rotationAngle) * Ein2[i];
    Eout2[i] = sin(rotationAngle) * Ein1[i] + cos(rotationAngle) * Ein2[i];
}

//provide a list of nearest-3 neighbors for taking spatial derivatives
// exploiting the fact that the radial grid is offset by 1/4 step from 0
// this means that midpoints are available on the other side of the origin.
// returns rho at the given index j
__device__ __forceinline__ double resolveNeighborsInOffsetRadialSymmetry(
    long long* neighbors, long long N, int j, double dr, long long Ntime, long long h) {
	if (j < N / 2) {
		neighbors[0] = (N - j - 2) * Ntime + h;
		neighbors[1] = (j + 1) * Ntime + h;
		neighbors[2] = (N - j - 1) * Ntime + h;
		neighbors[3] = (N - j) * Ntime + h;
		neighbors[4] = (j - 1) * Ntime + h;
		neighbors[5] = (N - j + 1) * Ntime + h;
		return -(dr * (j - N / 2) + 0.25 * dr);
	}
	else {
		neighbors[0] = (N - j + 1) * Ntime + h;
		neighbors[1] = (j - 1) * Ntime + h;
		neighbors[2] = (N - j) * Ntime + h;
		neighbors[3] = (N - j - 1) * Ntime + h;
		neighbors[4] = (j + 1) * Ntime + h;
		neighbors[5] = (N - j - 2) * Ntime + h;
		return dr * (j - N / 2) + 0.25 * dr;
	}
}

__global__ void radialLaplacianKernel(cudaParameterSet* s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j = i / (*s).Ntime; //spatial coordinate
    long long h = i % (*s).Ntime; //temporal coordinate
    long long neighbors[6];

    //zero at edges of grid
    if (j<3 || j>((*s).Nspace - 4)) {
        (*s).gridRadialLaplacian1[i] = 0.;
        (*s).gridRadialLaplacian2[i] = 0.;
    }
    else {
        double rho = resolveNeighborsInOffsetRadialSymmetry(neighbors, (*s).Nspace, j, (*s).dx, (*s).Ntime, h);
        rho = -1.0 / rho;
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

}
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
__global__ void expandCylindricalBeam(cudaParameterSet* s, double* polarization1, double* polarization2) {
    size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    size_t j = i / (*s).Ntime; //spatial coordinate
    size_t k = i % (*s).Ntime; //temporal coordinate

    //positions on the expanded grid corresponding the the current index
    size_t pos1 = 2 * ((*s).Nspace - j - 1) * (*s).Ntime + k;
    size_t pos2 = (2 * j + 1) * (*s).Ntime + k;

    //reuse memory allocated for the radial Laplacian, casting complex double
    //to a 2x larger double real grid
    double* expandedBeam1 = (double*)(*s).gridRadialLaplacian1;
    double* expandedBeam2 = expandedBeam1 + 2 * (*s).Ngrid;

    expandedBeam1[pos1] = polarization1[i];
    expandedBeam1[pos2] = polarization1[i];
    expandedBeam2[pos1] = polarization2[i];
    expandedBeam2[pos2] = polarization2[i];
}
__device__ void findBirefingentCrystalAngle(double* alphaE, double* alphaO, long long j, double f, double* sellmeierCoefficients, cudaParameterSet s) {
    //Find walkoff angle, starting from zero
    // in the case of an extraordinary axis, the angle of propagation is related to the transverse
    // momentum in a complicated way:
    // sin(theta) * n(theta) = delta k * c/omega
    // theta depends on the refractive index, and the refractive index depends on theta
    // so we solve numerically
    double dAlpha = 0.1;
    double nePlus, neMinus;
    double err, errPlus, errMinus;
    cuDoubleComplex ne, no;

    
    cuDoubleComplex ii = make_cuDoubleComplex(0, 1);
    double crystalTheta = sellmeierCoefficients[66];
    double crystalPhi = sellmeierCoefficients[67];
    double kStep = sellmeierCoefficients[70];
    double tol = sellmeierCoefficients[72];
    double dk = j * kStep - (j >= (s.Nspace / 2)) * (kStep * s.Nspace); //frequency grid in transverse direction
    double rhs = LIGHTC * dk / (TWOPI * f);

    //if not biaxial, the o-axis can be solved analytically.
    sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
        crystalTheta, crystalPhi, s.axesNumber, s.sellmeierType);
    *alphaO = asin(rhs / cuCreal(no));
    if(s.axesNumber == 2){
        sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
            crystalTheta + *alphaO, crystalPhi, s.axesNumber, s.sellmeierType);
        nePlus = cuCreal(no);
        err = abs(nePlus * sin(*alphaO) - rhs);

        int iters = 0;
        errPlus = 2;
        errMinus = 2;
        while (err > tol && iters < 2048) {
            iters++;

            sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
                crystalTheta + *alphaO + dAlpha, crystalPhi, s.axesNumber, s.sellmeierType);
            nePlus = cuCreal(no);
            errPlus = abs(nePlus * sin(*alphaO + dAlpha) - rhs);

            sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
                crystalTheta + *alphaO - dAlpha, crystalPhi, s.axesNumber, s.sellmeierType);
            neMinus = cuCreal(no);
            errMinus = abs(neMinus * sin(*alphaO - dAlpha) - rhs);

            //Basic hill climbing algorithm
            //calculate the error at theta +/- dTheta
            // if theta + dTheta has lowest error, theta = theta+dTheta, err = errPlus
            // if theta - dTheta has lowest error, theta = theta-dTheta, err = errMinus
            // if theta has lowest error, step size is too large, dTheta /= 2;
            if (errPlus < err && errPlus < errMinus) {
                *alphaO += dAlpha;
                err = errPlus;
            }
            else if (errMinus < err) {
                *alphaO -= dAlpha;
                err = errMinus;
            }
            else {
                dAlpha *= 0.5;
            }

        }
    }

    //find the extraordinary angle if the crystal isn't isotropic
    *alphaE = *alphaO;
    if (s.axesNumber > 0) {
        sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
            crystalTheta + *alphaE, crystalPhi, s.axesNumber, s.sellmeierType);
        nePlus = cuCreal(ne);
        err = abs(nePlus * sin(*alphaE) - rhs);

        int iters = 0;
        errPlus = 2;
        errMinus = 2;
        dAlpha = 0.1;
        while (err > tol && iters < 2048) {
            iters++;

            sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
                crystalTheta + *alphaE + dAlpha, crystalPhi, s.axesNumber, s.sellmeierType);
            nePlus = cuCreal(ne);
            errPlus = abs(nePlus * sin(*alphaE + dAlpha) - rhs);

            sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
                crystalTheta + *alphaE - dAlpha, crystalPhi, s.axesNumber, s.sellmeierType);
            neMinus = cuCreal(ne);
            errMinus = abs(neMinus * sin(*alphaE - dAlpha) - rhs);

            //Basic hill climbing algorithm
            //calculate the error at theta +/- dTheta
            // if theta + dTheta has lowest error, theta = theta+dTheta, err = errPlus
            // if theta - dTheta has lowest error, theta = theta-dTheta, err = errMinus
            // if theta has lowest error, step size is too large, dTheta /= 2;
            if (errPlus < err && errPlus < errMinus) {
                *alphaE += dAlpha;
                err = errPlus;
            }
            else if (errMinus < err) {
                *alphaE -= dAlpha;
                err = errMinus;
            }
            else {
                dAlpha *= 0.5;
            }

        }
    }
    
    
}


//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
//note that the sellmeier coefficients have extra values appended to the end
//to give info about the current simulation
__global__ void applyFresnelLossKernel(double* sellmeierCoefficients1, double* sellmeierCoefficients2, cudaParameterSet s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    double alpha1, alpha2, alphaO1, alphaO2;
    long long j, k;
    long long Ntime = s.Ntime;
    int axesNumber = s.axesNumber;
    int sellmeierType = s.sellmeierType;
    cuDoubleComplex ne1, no1, ne2, no2, n0;
    cuDoubleComplex cuZero = make_cuDoubleComplex(0, 0);
    j = i / Ntime; //spatial coordinate
    k = i % Ntime; //temporal coordinate
    cuDoubleComplex ii = make_cuDoubleComplex(0, 1);
    double crystalTheta = sellmeierCoefficients1[66];
    double crystalPhi = sellmeierCoefficients1[67];
    double fStep = sellmeierCoefficients1[71];

    //frequency being resolved by current thread
    double f = k * fStep;
    if (k >= Ntime / 2) {
        f -= fStep * Ntime;
    }
    f *= -1;

    findBirefingentCrystalAngle(&alpha1, &alphaO1, j, f, sellmeierCoefficients1, s);
    findBirefingentCrystalAngle(&alpha2, &alphaO2, j, f, sellmeierCoefficients2, s);
    //walkoff angle has been found, generate the rest of the grids


    sellmeierCuda(&ne1, &no1, sellmeierCoefficients1, abs(f),
        crystalTheta + alpha1, crystalPhi, axesNumber, sellmeierType);
    sellmeierCuda(&n0, &no1, sellmeierCoefficients1, abs(f),
        crystalTheta + alphaO1, crystalPhi, axesNumber, sellmeierType);
    if (isnan(cuCreal(ne1)) || isnan(cuCreal(no1))) {
        ne1 = make_cuDoubleComplex(1, 0);
        no1 = make_cuDoubleComplex(1, 0);
    }


    sellmeierCuda(&ne2, &no2, sellmeierCoefficients2, abs(f),
        crystalTheta + alpha2, crystalPhi, axesNumber, sellmeierType);
    sellmeierCuda(&n0, &no2, sellmeierCoefficients2, abs(f),
        crystalTheta + alphaO2, crystalPhi, axesNumber, sellmeierType);
    if (isnan(cuCreal(ne2)) || isnan(cuCreal(no2))) {
        ne2 = make_cuDoubleComplex(1, 0);
        no2 = make_cuDoubleComplex(1, 0);
    }

    cuDoubleComplex ts = 2 * ne1 * cos(alpha1) / (ne1 * cos(alpha1) + ne2 * cos(alpha2));
    cuDoubleComplex tp = 2 * ne1 * cos(alpha1) / (ne2 * cos(alpha1) + ne1 * cos(alpha2));
    if (isnan(ts.x) || isnan(ts.y)) ts = make_cuDoubleComplex(0, 0);
    if (isnan(tp.x) || isnan(tp.y)) ts = make_cuDoubleComplex(0, 0);
    s.gridEFrequency1[i] = ts * s.gridEFrequency1[i];
    s.gridEFrequency2[i] = tp * s.gridEFrequency2[i];
}


__global__ void applyLinearPropagationKernel(double* sellmeierCoefficients, double thickness, cudaParameterSet s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    double alpha1, alphaO1;
    long long j, k;
    long long Ntime = s.Ntime;
    int axesNumber = s.axesNumber;
    int sellmeierType = s.sellmeierType;
    cuDoubleComplex ne, no, n0, n0o;
    cuDoubleComplex cuZero = make_cuDoubleComplex(0, 0);
    j = i / Ntime; //spatial coordinate
    k = i % Ntime; //temporal coordinate
    cuDoubleComplex ii = make_cuDoubleComplex(0, 1);
    double crystalTheta = sellmeierCoefficients[66];
    double crystalPhi = sellmeierCoefficients[67];
    double fStep = sellmeierCoefficients[71];

    //frequency being resolved by current thread
    double f = k * fStep;
    if (k >= Ntime / 2) {
        f -= fStep * Ntime;
    }
    f *= -1;
    double omega = TWOPI * f;
    findBirefingentCrystalAngle(&alpha1, &alphaO1, j, f, sellmeierCoefficients, s);
    //walkoff angle has been found, generate the rest of the grids


    sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f),
        crystalTheta + alpha1, crystalPhi, axesNumber, sellmeierType);
    sellmeierCuda(&n0, &n0o, sellmeierCoefficients, abs(f),
        crystalTheta + alphaO1, crystalPhi, axesNumber, sellmeierType);
    if (isnan(cuCreal(ne)) || isnan(cuCreal(no))) {
        ne = make_cuDoubleComplex(1, 0);
        no = make_cuDoubleComplex(1, 0);
    }



    cuDoubleComplex ts = cuCexpd(ii * ((ne-n0) * omega / LIGHTC) * thickness);
    cuDoubleComplex tp = cuCexpd(ii * ((no-n0) * omega / LIGHTC) * thickness);
    if (isnan(ts.x) || isnan(ts.y)) ts = make_cuDoubleComplex(0, 0);
    if (isnan(tp.x) || isnan(tp.y)) tp = make_cuDoubleComplex(0, 0);
    s.gridEFrequency1[i] = ts * s.gridEFrequency1[i];
    s.gridEFrequency2[i] = tp * s.gridEFrequency2[i];
}


//prepare propagation constants for the simulation, when it is taking place on a Cartesian grid
//note that the sellmeier coefficients have extra values appended to the end
//to give info about the current simulation
__global__ void prepareCartesianGridsKernel(double* sellmeierCoefficients, cudaParameterSet s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    double alpha, alphaO;
    long long j, k;
    long long Ntime = s.Ntime;
    int axesNumber = s.axesNumber;
    int sellmeierType = s.sellmeierType;
    cuDoubleComplex ne, no, n0;
    cuDoubleComplex cuZero = make_cuDoubleComplex(0, 0);
    j = i / Ntime; //spatial coordinate
    k = i % Ntime; //temporal coordinate
    if (k > s.propagationInts[3]) {
        return;
    }
    long long iC = s.propagationInts[3] * j + k;
    cuDoubleComplex ii = make_cuDoubleComplex(0, 1);
    double crystalTheta = sellmeierCoefficients[66];
    double crystalPhi = sellmeierCoefficients[67];
    double kStep = sellmeierCoefficients[70];
    double fStep = sellmeierCoefficients[71];

    //frequency being resolved by current thread
    double f = k * fStep;
    if (k >= Ntime / 2) {
        f -= fStep * Ntime;
    }
    f *= -1;

    //transverse wavevector being resolved
    double dk = j * kStep - (j >= (s.Nspace / 2)) * (kStep * s.Nspace); //frequency grid in transverse direction
    
    findBirefingentCrystalAngle(&alpha, &alphaO, j, f, sellmeierCoefficients, s);

    //walkoff angle has been found, generate the rest of the grids

    sellmeierCuda(&n0, &no, sellmeierCoefficients, abs(s.f0), 
        crystalTheta, crystalPhi, axesNumber, sellmeierType);
    sellmeierCuda(&ne, &no, sellmeierCoefficients, abs(f), 
        crystalTheta + alpha, crystalPhi, axesNumber, sellmeierType);
    if (s.axesNumber == 2) {
        cuDoubleComplex nDummy;
        sellmeierCuda(&nDummy, &no, sellmeierCoefficients, abs(f),
            crystalTheta + alphaO, crystalPhi, axesNumber, sellmeierType);
    }
    if (isnan(cuCreal(ne)) || isnan(cuCreal(no))) {
        ne = make_cuDoubleComplex(1, 0);
        no = make_cuDoubleComplex(1, 0);
    }
    s.ne[i] = ne;
    s.no[i] = no;

    cuDoubleComplex k0 = make_cuDoubleComplex(TWOPI * cuCreal(n0) * f / LIGHTC, 0);
    cuDoubleComplex ke = TWOPI * ne * f / LIGHTC;
    cuDoubleComplex ko = TWOPI * no * f / LIGHTC;

    if (s.isUsingMillersRule) {
        s.chiLinear1[iC] = -1. + ne * ne;
        s.chiLinear2[iC] = -1. + no * no;
        if ((cuCreal(s.chiLinear1[iC]) == 0) || (cuCreal(s.chiLinear2[iC]) == 0) || isnan(cuCreal(s.chiLinear1[iC])) || isnan(cuCreal(s.chiLinear2[iC]))) {
            s.chiLinear1[iC] = make_cuDoubleComplex(1, 0);
            s.chiLinear2[iC] = make_cuDoubleComplex(1, 0);
        }
    }
    else {
        s.chiLinear1[iC] = make_cuDoubleComplex(1, 0);
        s.chiLinear2[iC] = make_cuDoubleComplex(1, 0);
    }

    if (cuCreal(ke) < 0 && cuCreal(ko) < 0) {
        s.gridPropagationFactor1[iC] = ii * (ke - k0 - dk * dk / (2. * cuCreal(ke))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor1[iC]))) {
            s.gridPropagationFactor1[iC] = cuZero;
        }

        s.gridPropagationFactor2[iC] = ii * (ko - k0 - dk * dk / (2. * cuCreal(ko))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor2[iC]))) {
            s.gridPropagationFactor2[iC] = cuZero;
        }

        s.gridPolarizationFactor1[iC] = ii * s.chiLinear1[iC] * (TWOPI * f) / (2. * cuCreal(ne) * LIGHTC) * s.h;
        s.gridPolarizationFactor2[iC] = ii * s.chiLinear2[iC] * (TWOPI * f) / (2. * cuCreal(no) * LIGHTC) * s.h;
    }

    else {
        s.gridPropagationFactor1[iC] = cuZero;
        s.gridPropagationFactor2[iC] = cuZero;
        s.gridPolarizationFactor1[iC] = cuZero;
        s.gridPolarizationFactor2[iC] = cuZero;
    }
}
    
//prepare the propagation constants under the assumption of cylindrical symmetry of the beam
__global__ void prepareCylindricGridsKernel(double* sellmeierCoefficients, cudaParameterSet s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long j, k;
    long long Ntime = s.Ntime;
    long long Nspace = s.Nspace;
    int axesNumber = s.axesNumber;
    int sellmeierType = s.sellmeierType;
    cuDoubleComplex cuZero = make_cuDoubleComplex(0, 0);
    j = i / Ntime; //spatial coordinate
    k = i % Ntime; //temporal coordinate
    if (k > s.propagationInts[3]) {
        return;
    }
    long long iC = s.propagationInts[3] * j + k;

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

    cuDoubleComplex k0 = make_cuDoubleComplex(TWOPI * cuCreal(n0) * f / LIGHTC, 0);
    cuDoubleComplex ke = TWOPI * ne * f / LIGHTC;
    cuDoubleComplex ko = TWOPI * no * f / LIGHTC;

    if (s.isUsingMillersRule) {
        s.chiLinear1[iC] = -1. + ne * ne;
        s.chiLinear2[iC] = -1. + no * no;
        if ((cuCreal(s.chiLinear1[iC]) == 0) || (cuCreal(s.chiLinear2[iC]) == 0) || isnan(cuCreal(s.chiLinear1[iC])) || isnan(cuCreal(s.chiLinear2[iC]))) {
            s.chiLinear1[iC] = make_cuDoubleComplex(1, 0);
            s.chiLinear2[iC] = make_cuDoubleComplex(1, 0);
        }
    }
    else {
        s.chiLinear1[iC] = make_cuDoubleComplex(1, 0);
        s.chiLinear2[iC] = make_cuDoubleComplex(1, 0);
    }

    if (cuCreal(ke) < 0 && cuCreal(ko) < 0 && abs(dk) < cuCabs(ke)) {
        s.gridPropagationFactor1[iC] = ii * (ke - k0 - dk * dk / (2. * cuCreal(ke))) * s.h;
        s.gridPropagationFactor1Rho1[iC] = ii * (1 / (s.chiLinear1[iC] *2. * cuCreal(ke))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor1[iC]))) {
            s.gridPropagationFactor1[iC] = cuZero;
            s.gridPropagationFactor1Rho1[iC] = cuZero;
        }

        s.gridPropagationFactor2[iC] = ii * (ko - k0 - dk * dk / (2. * cuCreal(ko))) * s.h;
        s.gridPropagationFactor1Rho2[iC] = ii * (1 / (s.chiLinear2[iC] * 2. * cuCreal(ko))) * s.h;
        if (isnan(cuCreal(s.gridPropagationFactor2[iC]))) {
            s.gridPropagationFactor2[iC] = cuZero;
            s.gridPropagationFactor1Rho2[iC] = cuZero;
        }
        //factor of 0.5 comes from doubled grid size in cylindrical symmetry mode after expanding the beam
        s.gridPolarizationFactor1[iC] = 0.5 * s.chiLinear1[iC] * ii * (TWOPI * f) / (2. * cuCreal(ne) * LIGHTC) * s.h;
        s.gridPolarizationFactor2[iC] = 0.5 * s.chiLinear2[iC] * ii * (TWOPI * f) / (2. * cuCreal(no) * LIGHTC) * s.h;


    }

    else {
        s.gridPropagationFactor1[iC] = cuZero;
        s.gridPropagationFactor2[iC] = cuZero;
        s.gridPolarizationFactor1[iC] = cuZero;
        s.gridPolarizationFactor2[iC] = cuZero;
        s.gridPropagationFactor1[iC] = cuZero;
        s.gridPropagationFactor1Rho2[iC] = cuZero;
    }


}

//replaces E with its complex conjugate
__global__ void conjugateKernel(cuDoubleComplex* E) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    E[i] = cuConj(E[i]);
}

__global__ void realToComplexKernel(double* in, cuDoubleComplex* out) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    out[i] = make_cuDoubleComplex(in[i], 0.0);
}

__global__ void complexToRealKernel(cuDoubleComplex* in, double* out) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    out[i] = cuCreal(in[i]);
}

__global__ void materialPhaseKernel(double df, size_t Ntime, double* a, double f01, double f02, double thickness1,  double thickness2, double* phase1, double* phase2) {
    size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    //frequency being resolved by current thread
    double f = i * df;
    if (i >= Ntime / 2) {
        f -= df * Ntime;
    }
    
    //give phase shift relative to group velocity (approximated 
    // with low-order finite difference) so the pulse doesn't move
    cuDoubleComplex ne, no, no0, n0p, n0m;
    sellmeierCuda(&ne, &no, a, abs(f), 0, 0, 0, 0);
    f *= -TWOPI;
    sellmeierCuda(&ne, &no0, a, f01, 0, 0, 0, 0);
    sellmeierCuda(&ne, &n0p, a, f01+1e11, 0, 0, 0, 0);
    sellmeierCuda(&ne, &n0m, a, f01-1e11, 0, 0, 0, 0);
    no0 = no0 + f01 * (n0p - n0m) / 2e11;
    phase1[i] = thickness1 * f * cuCreal(no - no0) / LIGHTC;
    sellmeierCuda(&ne, &no0, a, f02, 0, 0, 0, 0);
    sellmeierCuda(&ne, &n0p, a, f02 + 1e11, 0, 0, 0, 0);
    sellmeierCuda(&ne, &n0m, a, f02 - 1e11, 0, 0, 0, 0);
    no0 = no0 + f02 * (n0p - n0m) / 2e11;
    phase2[i] = thickness2 * f * cuCreal(no - no0) / LIGHTC;
    
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
__global__ void nonlinearPolarizationKernel(cudaParameterSet* s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    double Ex = (*s).gridETime1[i] / (*s).propagationInts[0];
    double Ey = (*s).gridETime2[i] / (*s).propagationInts[0];

    double Ex2 = Ex * Ex;
    double Ey2 = Ey * Ey;
    (*s).gridPolarizationTime1[i] = 0.;
    (*s).gridPolarizationTime2[i] = 0.;

    //The d2eff tensor has the form
    // | d_xxx d_xyx d_yyx |
    // | d_xxy d_xyy d_yyy |
    if ((*s).nonlinearSwitches[0] == 1) {
        (*s).gridPolarizationTime1[i] += (*s).chi2Tensor[0] * Ex2 + (*s).chi2Tensor[2] * Ex * Ey + (*s).chi2Tensor[4] * Ey2;
        (*s).gridPolarizationTime2[i] += (*s).chi2Tensor[1] * Ex2 + (*s).chi2Tensor[3] * Ex * Ey + (*s).chi2Tensor[5] * Ey2;
    }
    
    //resolve the full chi3 matrix when (*s).nonlinearSwitches[1]==1
    if ((*s).nonlinearSwitches[1] == 1) {

        //rotate field into crystal frame
        double E3[3] = { (*s).rotationForward[0] * Ex + (*s).rotationForward[1] * Ey, 
            (*s).rotationForward[3] * Ex + (*s).rotationForward[4] * Ey,
            (*s).rotationForward[6] * Ex + (*s).rotationForward[7] * Ey };
        
        //loop over tensor element X_abcd
        //i hope the compiler unrolls this, but no way am I writing that out by hand
        unsigned char a, b, c, d;
        double P3[3] = { 0 };
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++) {
                for (c = 0; c < 3; c++) {
                    for (d = 0; d < 3; d++) {
                        P3[d] += (*s).chi3Tensor[a + 3 * b + 9 * c + 27 * d] * E3[a] * E3[b] * E3[c];
                    }
                }
            }
        }

        //rotate back into simulation frame
        (*s).gridPolarizationTime1[i] += (*s).rotationBackward[0] * P3[0] + (*s).rotationBackward[1] * P3[1] + (*s).rotationBackward[2] * P3[2];
        (*s).gridPolarizationTime2[i] += (*s).rotationBackward[3] * P3[0] + (*s).rotationBackward[4] * P3[1] + (*s).rotationBackward[5] * P3[2];
    }
    //using only one value of chi3, under assumption of centrosymmetry
    if ((*s).nonlinearSwitches[1] == 2) {
        double Esquared = (*s).chi3Tensor[0] * (Ex2 + Ey2);
        (*s).gridPolarizationTime1[i] += Ex * Esquared;
        (*s).gridPolarizationTime2[i] += Ey * Esquared;
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
__global__ void plasmaCurrentKernelPrep(cudaParameterSet* s, double* workN, double* workEx) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;

    int k;
    double* workEy = &workEx[(*s).Ngrid];
    double w, Esquared, Ex, Ey;
    Ex = 0.5 * (*s).gridETime1[i] / (*s).propagationInts[0];
    Ey = 0.5 * (*s).gridETime2[i] / (*s).propagationInts[0];
    Esquared = Ex * Ex + Ey * Ey;
    //plasmaParameters[0] is the nonlinear absorption parameter
    w = (*s).plasmaParameters[0] * Esquared;
    //nonlinearSwitches[3] is Nphotons-2
    for (k = 0; k < (*s).nonlinearSwitches[3]; k++) {
        w *= Esquared;
    }
    //absorption currents
    (*s).gridPlasmaCurrent1[i] = w * Ex;
    (*s).gridPlasmaCurrent2[i] = w * Ey;

    //plasmaParameters[2] is the 1/photon energy, translating the loss of power
    //from the field to the number of free carriers
    //extra factor of (dt^2e^2/(m*photon energy*eo) included as it is needed for the amplitude
    //of the plasma current
    workN[i] = (*s).plasmaParameters[2] * ((*s).gridPlasmaCurrent1[i] * Ex + (*s).gridPlasmaCurrent2[i] * Ey);
    workEx[i] = Ex;
    workEy[i] = Ey;

}
__global__ void plasmaCurrentKernel2(cudaParameterSet* s, double* workN, double* workEx) {
    long long j = threadIdx.x + blockIdx.x * blockDim.x;
    double N = 0;
    double integralx = 0;
    double integraly = 0;
    double* workEy = &workEx[(*s).Ngrid];
    double* expMinusGammaT = &(*s).expGammaT[(*s).Ntime];

    long long k;
    j *= (*s).Ntime;
    double a;
    for (k = 0; k < (*s).Ntime; k++) {
        N += workN[j];
        a = N * (*s).expGammaT[k];
        integralx += a * workEx[j];
        integraly += a * workEy[j];
        (*s).gridPlasmaCurrent1[j] += expMinusGammaT[k] * integralx;
        (*s).gridPlasmaCurrent2[j] += expMinusGammaT[k] * integraly;
        j++;
    }
}
//size_t propagationIntsCPU[4] = { s.Ngrid, s.Ntime, s.Nspace, (s.Ntime / 2 + 1) };
__device__ cuDoubleComplex findFFTValueD2Z_2D(cuDoubleComplex* Arr, long long* propagationInts, long long i) {
    long long h = i % propagationInts[1];
    if (h > propagationInts[3]) {
        h = propagationInts[1] - h + (propagationInts[2] - i / propagationInts[1]) * propagationInts[3];
        return cuConj(Arr[h]);
    }
    else {
        h += (i / propagationInts[1]) * propagationInts[3];
        return Arr[h];
    }
}

//Main kernel for RK4 propagation of the field
__global__ void rkKernel(cudaParameterSet *sP, int stepNumber) {
    long long iO = threadIdx.x + blockIdx.x * blockDim.x;
    long long i = (iO % (*sP).propagationInts[3]) + (iO / (*sP).propagationInts[3]) * (*sP).Ntime;
    long long h = i % ((*sP).Ntime); //temporal coordinate
    long long j = i / ((*sP).Ntime); //spatial coordinate

    long long iC = (*sP).propagationInts[3] * j + h;
    //long long iC = ((*sP).Ntime / 2) * j + h;
    cuDoubleComplex plasmaJ1 = make_cuDoubleComplex(0, 0);
    cuDoubleComplex plasmaJ2 = make_cuDoubleComplex(0, 0);

    if ((*sP).isCylindric) {
        //correct for double grid size after expandCylindricalBeam()
        if (j > (*sP).Nspace / 2) {
            j += (*sP).Nspace;
        }
    }
    //polarization is stored in a reduced format by cuFFT because the FFT is from real to complex, meaning if the output grid
    //were to be N_time x N_space, half of the points would be redundant. The extra steps below are to determine where in the grid the 
    //current point sit(*sP). Essentially, if in the negative frequency quadrants, reverse the frequency and take complex conjugate of the 
    //value
    h += j * (*sP).propagationInts[3];

    (*sP).k1[iC] = (*sP).gridPropagationFactor1[iC] * (*sP).gridETemp1[iC] + (*sP).gridPolarizationFactor1[iC] * (*sP).gridPolarizationFrequency1[h];
    (*sP).k2[iC] = (*sP).gridPropagationFactor2[iC] * (*sP).gridETemp2[iC] + (*sP).gridPolarizationFactor2[iC] * (*sP).gridPolarizationFrequency2[h];

    if ((*sP).hasPlasma) {
        cuDoubleComplex jfac = make_cuDoubleComplex(0, -1.0 / (TWOPI * ((i % ((*sP).Ntime/2)) * (*sP).fStep)));
        if ((i % ((*sP).Ntime/2)) != 0) {
            (*sP).k1[iC] = (*sP).k1[iC] + jfac * (*sP).gridPolarizationFactor1[iC] * (*sP).gridPlasmaCurrentFrequency1[h] / (*sP).chiLinear1[iC];
            (*sP).k2[iC] = (*sP).k2[iC] + jfac * (*sP).gridPolarizationFactor2[iC] * (*sP).gridPlasmaCurrentFrequency2[h] / (*sP).chiLinear2[iC];
        }
    }

    if ((*sP).isCylindric) {
        (*sP).k1[iC] = (*sP).k1[iC] + (*sP).gridPropagationFactor1Rho1[iC] * (*sP).workspace1C[iC];
        (*sP).k2[iC] = (*sP).k2[iC] + (*sP).gridPropagationFactor1Rho2[iC] * (*sP).workspace2C[iC];
    }

    //in the first substep, first construct the next intermediate field value
    //which will be used in the next substep. 
    if (stepNumber == 0) {
        (*sP).gridETemp1[iC] = (*sP).gridEFrequency1[iC] + 0.5 * (*sP).k1[iC];
        (*sP).gridETemp2[iC] = (*sP).gridEFrequency2[iC] + 0.5 * (*sP).k2[iC];
        (*sP).gridEFrequency1Next1[iC] = (*sP).k1[iC] / 6 + (*sP).gridEFrequency1[iC];
        (*sP).gridEFrequency1Next2[iC] = (*sP).k2[iC] / 6 + (*sP).gridEFrequency2[iC];
        (*sP).workspace1C[iC] = (*sP).chiLinear1[iC] * (*sP).gridETemp1[iC];
        (*sP).workspace2C[iC] = (*sP).chiLinear2[iC] * (*sP).gridETemp2[iC];
    }

    //in the next substep, again construct the next intermediate field and add k/3 to solution
    else if (stepNumber == 1) {
        (*sP).gridETemp1[iC] = (*sP).gridEFrequency1[iC] + 0.5 * (*sP).k1[iC];
        (*sP).gridETemp2[iC] = (*sP).gridEFrequency2[iC] + 0.5 * (*sP).k2[iC];
        (*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).k1[iC] / 3;
        (*sP).gridEFrequency1Next2[iC] = (*sP).gridEFrequency1Next2[iC] + (*sP).k2[iC] / 3;
        (*sP).workspace1C[iC] = (*sP).chiLinear1[iC] * (*sP).gridETemp1[iC];
        (*sP).workspace2C[iC] = (*sP).chiLinear2[iC] * (*sP).gridETemp2[iC];
    }

    //same action as previous substep, except the weight of k in the intermediate solution is 1 instead of 0.5
    else if (stepNumber == 2) {
        (*sP).gridETemp1[iC] = (*sP).gridEFrequency1[iC] + (*sP).k1[iC];
        (*sP).gridETemp2[iC] = (*sP).gridEFrequency2[iC] + (*sP).k2[iC];
        (*sP).gridEFrequency1Next1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).k1[iC] / 3;
        (*sP).gridEFrequency1Next2[iC] = (*sP).gridEFrequency1Next2[iC] + (*sP).k2[iC] / 3;
        (*sP).workspace1C[iC] = (*sP).chiLinear1[iC] * (*sP).gridETemp1[iC];
        (*sP).workspace2C[iC] = (*sP).chiLinear2[iC] * (*sP).gridETemp2[iC];
    }

    //last substep. Solution is now complete and may be copied directly into the field arrays
    else {
        (*sP).gridEFrequency1[iC] = (*sP).gridEFrequency1Next1[iC] + (*sP).k1[iC] / 6;
        (*sP).gridEFrequency2[iC] = (*sP).gridEFrequency1Next2[iC] + (*sP).k2[iC] / 6;
        (*sP).gridETemp1[iC] = (*sP).gridEFrequency1[iC];
        (*sP).gridETemp2[iC] = (*sP).gridEFrequency2[iC];
        (*sP).workspace1C[iC] = (*sP).chiLinear1[iC] * (*sP).gridETemp1[iC];
        (*sP).workspace2C[iC] = (*sP).chiLinear2[iC] * (*sP).gridETemp2[iC];
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


//Apply fft normalization
__global__ void multiplyByConstantKernel(cuDoubleComplex* A, double val) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    A[i] = val * A[i];
}

//Apply fft normalization
__global__ void multiplyByConstantKernelD(double* A, double val) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    A[i] = val * A[i];
}


//element-wise B*A = C;
__global__ void multiplicationKernel(cuDoubleComplex* A, cuDoubleComplex* B, cuDoubleComplex* C) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    C[i] = B[i] * A[i];
}

__global__ void multiplicationKernelCompact(cuDoubleComplex* A, cuDoubleComplex* B, cuDoubleComplex* C, cudaParameterSet s) {
    long long i = threadIdx.x + blockIdx.x * blockDim.x;
    long long h = i % s.Ntime; //temporal coordinate
    if (h > s.propagationInts[3]) {
        //C[i] = make_cuDoubleComplex(0., 0.);
        return;
    }
    long long j = i / s.Ntime; //spatial coordinate
    long long iC = s.propagationInts[3] * j + h;
    C[iC] = A[iC] * B[iC];
}


//main function for running on CLI
int main(int argc, char *argv[]) {
    int i, j;
    int CUDAdevice;
    int CUDAdeviceCount = 0;
    size_t progressCounter = 0;
    cudaGetDeviceCount(&CUDAdeviceCount);
    cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
    struct cudaDeviceProp activeCUDADeviceProp;
    if (cuErr == cudaSuccess) {
        printf("Found %i GPU(s): \n", CUDAdeviceCount);
        for (i = 0; i < CUDAdeviceCount; i++) {
            cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
            printf("%s\r\n", activeCUDADeviceProp.name);
            printf(" Memory: %lli MB; Multiprocessors: %i\n", 
                activeCUDADeviceProp.totalGlobalMem / (1024 * 1024), activeCUDADeviceProp.multiProcessorCount);
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
    simulationParameterSet* sCPU = (simulationParameterSet*)calloc(512, sizeof(simulationParameterSet));
    crystalEntry* crystalDatabasePtr = (crystalEntry*)calloc(512, sizeof(crystalEntry));
    (*sCPU).crystalDatabase = crystalDatabasePtr;
    (*sCPU).progressCounter = &progressCounter;
    // read crystal database
    if (readCrystalDatabase(crystalDatabasePtr) == -2) {
        return 11;
    }
    if ((*crystalDatabasePtr).numberOfEntries == 0) {
        printf("Could not read crystal database.\n");
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
        free(sCPU);
        free(crystalDatabasePtr);
        return 13;
    }

    allocateGrids(sCPU);
    if (loadPulseFiles(sCPU) == 1) {
        printf("Could not read pulse file.\n");
        free((*sCPU).refractiveIndex1);
        free((*sCPU).refractiveIndex2);
        free((*sCPU).imdone);
        free((*sCPU).deffTensor);
        free((*sCPU).loadedField1);
        free((*sCPU).loadedField2);
        free(sCPU);
        free(crystalDatabasePtr);
        return 14;
    }

    readSequenceString(sCPU);
    printf("Found %i steps in sequence\n", (*sCPU).Nsequence);
    readFittingString(sCPU);
    configureBatchMode(sCPU);

    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

    // run simulations
    if ((*sCPU).isInFittingMode) {
        if ((*sCPU).fittingMode == 3) {
            if (loadReferenceSpectrum((*sCPU).fittingPath, sCPU)) {
                printf("Could not load reference spectrum!\n");
                free((*sCPU).refractiveIndex1);
                free((*sCPU).refractiveIndex2);
                free((*sCPU).imdone);
                free((*sCPU).deffTensor);
                free((*sCPU).loadedField1);
                free((*sCPU).loadedField2);
                free((*sCPU).Ext);
                free((*sCPU).Ekw);
                free((*sCPU).ExtOut);
                free((*sCPU).EkwOut);
                free((*sCPU).totalSpectrum);
                free(sCPU);
                free(crystalDatabasePtr);
                return 10;
            }
        }
        printf("Running in fitting mode -- I don't know how long this will take!\n");
        runFitting(sCPU);
        
        auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
        printf("Finished after %8.4lf s. \n",
            1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
        
        saveDataSet(sCPU, crystalDatabasePtr, (*sCPU).outputBasePath, FALSE);
        //free
        free((*sCPU).refractiveIndex1);
        free((*sCPU).refractiveIndex2);
        free((*sCPU).imdone);
        free((*sCPU).deffTensor);
        free((*sCPU).loadedField1);
        free((*sCPU).loadedField2);
        free((*sCPU).Ext);
        free((*sCPU).Ekw);
        free((*sCPU).ExtOut);
        free((*sCPU).EkwOut);
        free((*sCPU).totalSpectrum);
        free(sCPU);
        free(crystalDatabasePtr);

        return 0;
    }
    std::thread *threadBlock = (std::thread*)calloc((*sCPU).Nsims * (*sCPU).Nsims2, sizeof(std::thread));
    size_t maxThreads = min(CUDAdeviceCount, (*sCPU).Nsims * (*sCPU).Nsims2);
    for (j = 0; j < (*sCPU).Nsims * (*sCPU).Nsims2; j++) {

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
    
	for (i = 0; i < (*sCPU).Nsims * (*sCPU).Nsims2; i++) {
        if (sCPU[i].memoryError > 0) {
            printf("Warning: device memory error (%i).\n", sCPU[i].memoryError);
        }
		if (threadBlock[i].joinable()) {
			threadBlock[i].join();
		}
	}
    
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    printf("Finished after %8.4lf s. \n", 
        1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));


    saveDataSet(sCPU, crystalDatabasePtr, (*sCPU).outputBasePath, FALSE);
    //free
    free(threadBlock);
    free((*sCPU).refractiveIndex1);
    free((*sCPU).refractiveIndex2);
    free((*sCPU).imdone);
    free((*sCPU).deffTensor);
    free((*sCPU).loadedField1);
    free((*sCPU).loadedField2);
    free((*sCPU).Ext);
    free((*sCPU).Ekw);
    free((*sCPU).ExtOut);
    free((*sCPU).EkwOut);
    free((*sCPU).totalSpectrum);
    free(sCPU);
    free(crystalDatabasePtr);
    return 0;
}

unsigned long solveNonlinearWaveEquationSequence(void* lpParam) {
    simulationParameterSet* sCPU = (simulationParameterSet*)lpParam;
    simulationParameterSet* sCPUbackup = (simulationParameterSet*)calloc(1, sizeof(simulationParameterSet));
    memcpy(sCPUbackup, sCPU, sizeof(simulationParameterSet));
    int k;
    for (k = 0; k < (*sCPU).Nsequence; k++) {
        resolveSequence(k, sCPU, (*sCPU).crystalDatabase);
        memcpy(sCPU, sCPUbackup, sizeof(simulationParameterSet));
    }
    free(sCPUbackup);
    return 0;
}
//main thread of the nonlinear wave equation implemented on CUDA
unsigned long solveNonlinearWaveEquation(void* lpParam) {

    //the struct s contains most of the simulation variables and pointers
    cudaParameterSet s;
    simulationParameterSet* sCPU = (simulationParameterSet*)lpParam;
    cudaSetDevice((*sCPU).assignedGPU);
    cudaStreamCreate(&s.CUDAStream);
    fillRotationMatricies(sCPU, &s);

    //initialize and take values from the struct handed over by the dispatcher
    unsigned long long i;
    s.Ntime = (*sCPU).Ntime;
    s.Nspace = (*sCPU).Nspace;
    s.dt = (*sCPU).tStep;
    s.dx = (*sCPU).rStep;
    s.fStep = (*sCPU).fStep;
    s.Nsteps = (size_t)round((*sCPU).crystalThickness / (*sCPU).propagationStep);
    s.h = (*sCPU).crystalThickness / (s.Nsteps); //adjust step size so that thickness can be varied continuously by fitting
    s.Ngrid = s.Ntime * s.Nspace;
    s.axesNumber = (*sCPU).axesNumber;
    s.sellmeierType = (*sCPU).sellmeierType;
    s.f0 = (*sCPU).frequency1;
    s.Nthread = THREADS_PER_BLOCK;
    s.Nblock = (int)(s.Ngrid / THREADS_PER_BLOCK);
    s.NblockC = (int)((s.Nspace * (s.Ntime / 2 + 1)) / THREADS_PER_BLOCK);
    s.isCylindric =(*sCPU).isCylindric;
    s.isNonLinear = ((*sCPU).nonlinearSwitches[0] + (*sCPU).nonlinearSwitches[1]) > 0;
    s.isUsingMillersRule = ((*sCPU).crystalDatabase[(*sCPU).materialIndex].nonlinearReferenceFrequencies[0]) != 0;
    size_t NgridC = (s.Ntime / 2 + 1) * s.Nspace; //size of the positive frequency side of the grid
    size_t beamExpansionFactor = 1;
    if (s.isCylindric) {
        beamExpansionFactor++;
    }
    //CPU allocations
    std::complex<double>* gridPropagationFactor1CPU = (std::complex<double>*)malloc(2 * s.Ngrid * sizeof(std::complex<double>));
    std::complex<double>* gridPolarizationFactor1CPU = (std::complex<double>*)malloc(2 * s.Ngrid * sizeof(std::complex<double>));
    
    //GPU allocations
    //I shouldn't need all these memsets but, they make me feel better
    int memErrors = 0;
    memErrors += cudaMalloc((void**)&s.gridETime1, 2 * sizeof(double) * s.Ngrid);
    cudaMemset(s.gridETime1, 0, 2 * sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.workspace1, 2 * sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.workspace1, 0, 2 * sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridETemp1, 2 * sizeof(cuDoubleComplex) * NgridC);
    cudaMemset(s.gridETemp1, 0, 2 * sizeof(cuDoubleComplex) * NgridC);
    memErrors += cudaMalloc((void**)&s.gridEFrequency1, 2 * sizeof(cuDoubleComplex) * NgridC);
    cudaMemset(s.gridEFrequency1, 0, 2 * sizeof(cuDoubleComplex) * NgridC);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactor1, 2 * sizeof(cuDoubleComplex) * NgridC);
    cudaMemset(s.gridPropagationFactor1, 0, 2 * sizeof(cuDoubleComplex) * NgridC);
    memErrors += cudaMalloc((void**)&s.gridPolarizationFactor1, 2 * sizeof(cuDoubleComplex) * NgridC);
    cudaMemset(s.gridPolarizationFactor1, 0, 2 * sizeof(cuDoubleComplex) * NgridC);
    memErrors += cudaMalloc((void**)&s.gridPropagationFactor1Rho1, 2 * sizeof(cuDoubleComplex) * NgridC);
    cudaMemset(s.gridPropagationFactor1Rho1, 0, 2 * sizeof(cuDoubleComplex) * NgridC);
    memErrors += cudaMalloc((void**)&s.gridRadialLaplacian1, 2 * sizeof(cuDoubleComplex) * s.Ngrid);
    cudaMemset(s.gridRadialLaplacian1, 0, 2 * sizeof(cuDoubleComplex) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridEFrequency1Next1, 2 * sizeof(cuDoubleComplex) * NgridC);
    cudaMemset(s.gridEFrequency1Next1, 0, 2 * sizeof(cuDoubleComplex) * NgridC);
    memErrors += cudaMalloc((void**)&s.chiLinear1, 2 * sizeof(cuDoubleComplex) * NgridC);
    cudaMemset(s.chiLinear1, 0, 2 * sizeof(cuDoubleComplex) * NgridC);
    memErrors += cudaMalloc((void**)&s.k1, 2 * sizeof(cuDoubleComplex) * NgridC);
    cudaMemset(s.k1, 0, 2 * sizeof(cuDoubleComplex) * NgridC);

    memErrors += cudaMalloc((void**)&s.gridPolarizationFrequency1, sizeof(cuDoubleComplex) * 2 * beamExpansionFactor * NgridC);
    cudaMemset(s.gridPolarizationFrequency1, 0, sizeof(cuDoubleComplex) * 2 * beamExpansionFactor * NgridC);
    memErrors += cudaMalloc((void**)&s.gridPlasmaCurrentFrequency1, sizeof(cuDoubleComplex) * 2 * beamExpansionFactor *  NgridC);
    cudaMemset(s.gridPlasmaCurrentFrequency1, 0, sizeof(cuDoubleComplex) * 2 * beamExpansionFactor * NgridC);
    memErrors += cudaMalloc((void**)&s.gridPolarizationTime1, 2 * sizeof(double) * s.Ngrid);
    cudaMemset(s.gridPolarizationTime1, 0, 2 * sizeof(double) * s.Ngrid);
    memErrors += cudaMalloc((void**)&s.gridPlasmaCurrent1, 2 * sizeof(double) * s.Ngrid);
    cudaMemset(s.gridPlasmaCurrent1, 0, 2 * sizeof(double) * s.Ngrid);
    

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
    memErrors += cudaMalloc((void**)&s.propagationInts, sizeof(long long) * 5);
    (*sCPU).memoryError = memErrors;
    if (memErrors > 0) {
        return memErrors;
    }

    //second polarization grids are to pointers within the first polarization
    //to have contiguous memory
    s.gridETime2 = s.gridETime1 + s.Ngrid;
    s.workspace2 = s.workspace1 + s.Ngrid;
    
    s.gridPlasmaCurrent2 = s.gridPlasmaCurrent1 + s.Ngrid;
    s.gridPolarizationTime2 = s.gridPolarizationTime1 + s.Ngrid;
    s.gridPlasmaCurrentFrequency2 = s.gridPlasmaCurrentFrequency1 + beamExpansionFactor * NgridC;
    s.gridPolarizationFrequency2 = s.gridPolarizationFrequency1 + beamExpansionFactor * NgridC;
    s.k2 = s.k1 + NgridC;
    s.chiLinear2 = s.chiLinear1 + NgridC;
    s.gridRadialLaplacian2 = s.gridRadialLaplacian1 + s.Ngrid;
    s.gridPropagationFactor1Rho2 = s.gridPropagationFactor1Rho1 + NgridC;
    s.gridPolarizationFactor2 = s.gridPolarizationFactor1 + NgridC;
    s.gridEFrequency1Next2 = s.gridEFrequency1Next1 + NgridC;
    s.gridPropagationFactor2 = s.gridPropagationFactor1 + NgridC;
    s.gridEFrequency2 = s.gridEFrequency1 + NgridC;
    s.gridETemp2 = s.gridETemp1 + NgridC;

    //prepare effective nonlinearity tensors and put them on the GPU
    size_t propagationIntsCPU[5] = { s.Ngrid, s.Ntime, s.Nspace, (s.Ntime / 2 + 1), NgridC };
    double firstDerivativeOperation[6] = { -1. / 60.,  3. / 20., -3. / 4.,  3. / 4.,  -3. / 20., 1. / 60. };
    for (i = 0; i < 6; i++) {
        firstDerivativeOperation[i] *= (-2.0/(s.Ngrid * s.dx));
    }

    //set nonlinearSwitches[3] to the number of photons needed to overcome bandgap
    (*sCPU).nonlinearSwitches[3] = (int)ceil((*sCPU).bandGapElectronVolts * 241.79893e12 / (*sCPU).frequency1) - 2;

    double plasmaParametersCPU[6] = { 0 };
    
    if ((*sCPU).nonlinearAbsorptionStrength > 0.) {
        s.hasPlasma = TRUE;
        s.isNonLinear = TRUE;
    }
    else {
        s.hasPlasma = FALSE;
    }
    
    plasmaParametersCPU[0] = (*sCPU).nonlinearAbsorptionStrength; //nonlinear absorption strength parameter
    plasmaParametersCPU[1] = (*sCPU).drudeGamma; //gamma
    if ((*sCPU).nonlinearAbsorptionStrength > 0.) {
        //plasmaParametersCPU[2] = (1. / 8.8541878128e-12) * (*sCPU).tStep * (*sCPU).tStep
            //* 2.817832e-08 / (1.6022e-19 * (*sCPU).bandGapElectronVolts * (*sCPU).effectiveMass); // (dt^2)*e* e / (m * band gap));
        plasmaParametersCPU[2] = (*sCPU).tStep * (*sCPU).tStep
            * 2.817832e-08 / (1.6022e-19 * (*sCPU).bandGapElectronVolts * (*sCPU).effectiveMass); // (dt^2)*e* e / (m * band gap));
    }
    else {
        plasmaParametersCPU[2] = 0;
    }
    

    calcEffectiveChi2Tensor((*sCPU).deffTensor, (*sCPU).chi2Tensor, (*sCPU).crystalTheta, (*sCPU).crystalPhi);
    cudaMemcpy(s.chi2Tensor, (*sCPU).deffTensor, 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.nonlinearSwitches, (*sCPU).nonlinearSwitches, 4 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(s.propagationInts, propagationIntsCPU, 5 * sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMemcpy(s.chi3Tensor, (*sCPU).chi3Tensor, 81 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.absorptionParameters, (*sCPU).absorptionParameters, 6 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.plasmaParameters, plasmaParametersCPU, 6 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(s.firstDerivativeOperation, firstDerivativeOperation, 6 * sizeof(double), cudaMemcpyHostToDevice);

    //prepare FFT plans
    size_t workSize;
    int cufftSizes1[] = { (int)s.Nspace, (int)s.Ntime };
    cufftCreate(&s.fftPlan);
    cufftGetSizeMany(s.fftPlan, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_Z2Z, 2, &workSize);
    cufftMakePlanMany(s.fftPlan, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_Z2Z, 2, &workSize);

    cufftCreate(&s.fftPlanD2Z);
    cufftGetSizeMany(s.fftPlanD2Z, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_D2Z, 2, &workSize);
    cufftMakePlanMany(s.fftPlanD2Z, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_D2Z, 2, &workSize);

    cufftCreate(&s.fftPlanZ2D);
    cufftGetSizeMany(s.fftPlanZ2D, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_Z2D, 2, &workSize);
    cufftMakePlanMany(s.fftPlanZ2D, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_Z2D, 2, &workSize);

    //the inputs are written out explicitly here even though 6x NULL would also do the same thing
    //this is in the hope that someday I figure out a way to do a truncated transform that doesn't
    //calculate the >f_nyquist/2 elements that I ignore
    int cufftSizes2[]{ 2*(int)s.Nspace, (int)s.Ntime };
    int istride = 1;
    int ostride = 1;
    int idist = (int)(2 * s.Ngrid);
    int odist = (int)(2 * s.Nspace * (s.Ntime / 2 + 1));
    int inembed[] = { (int)(2*s.Nspace), (int)s.Ntime };
    int onembed[] = { (int)(2*s.Nspace), (int)s.Ntime / 2 + 1 };

    cufftCreate(&s.doublePolfftPlan);
    cufftGetSizeMany(s.doublePolfftPlan, 2, cufftSizes2, inembed, istride, idist, onembed, ostride, odist, CUFFT_D2Z, 2, &workSize);
    cufftMakePlanMany(s.doublePolfftPlan, 2, cufftSizes2, inembed, istride, idist, onembed, ostride, odist, CUFFT_D2Z, 2, &workSize);

    cufftCreate(&s.polfftPlan);
    cufftGetSizeMany(s.polfftPlan, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_D2Z, 2, &workSize);
    cufftMakePlanMany(s.polfftPlan, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_D2Z, 2, &workSize);

    cufftSetStream(s.fftPlan, s.CUDAStream);
    cufftSetStream(s.polfftPlan, s.CUDAStream);
    cufftSetStream(s.doublePolfftPlan, s.CUDAStream);

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
        cudaMemcpy(s.workspace1, (*sCPU).ExtOut, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        cudaMemcpy(s.workspace2, &(*sCPU).ExtOut[(*sCPU).Ngrid], (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        complexToRealKernel << <s.Nblock, s.Nthread, 0, s.CUDAStream >> > (s.workspace1, s.gridETime1);
        complexToRealKernel << <s.Nblock, s.Nthread, 0, s.CUDAStream >> > (s.workspace2, s.gridETime2);
        cufftExecD2Z(s.fftPlanD2Z, s.gridETime1, s.gridEFrequency1);
        multiplyByConstantKernel << <(unsigned int)((s.Ntime / 2 + 1) * (s.Nspace/16)), 32, 0, s.CUDAStream >> > (s.gridEFrequency1, 2.0);

    }
    else {
        prepareElectricFieldArrays(sCPU, &s);
    }
    s.gridEFrequency2 = s.gridEFrequency1 + NgridC;
    s.gridETemp2 = s.gridETemp1 + NgridC;
    s.workspace1C = s.workspace1;
    s.workspace2C = s.workspace1C + NgridC;
    //Copy the field into the temporary array
    cudaMemcpy(s.gridETemp1, s.gridEFrequency1, NgridC * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    cudaMemcpy(s.gridETemp2, s.gridEFrequency2, NgridC * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    multiplicationKernelCompact<<<s.Nblock, s.Nthread, 0, s.CUDAStream>>>(s.chiLinear1, s.gridETemp1, s.workspace1C, s);
    multiplicationKernelCompact<<<s.Nblock, s.Nthread, 0, s.CUDAStream>>>(s.chiLinear2, s.gridETemp2, s.workspace2C, s);

    cudaParameterSet* sP;
    cudaMalloc(&sP, sizeof(cudaParameterSet));
    cudaMemcpy(sP, &s, sizeof(cudaParameterSet), cudaMemcpyHostToDevice);
    //Core propagation loop
    for (i = 0; i < s.Nsteps; i++) {
        
        //RK4
        runRK4Step(&s, sP, 0);
        runRK4Step(&s, sP, 1);
        runRK4Step(&s, sP, 2);
        runRK4Step(&s, sP, 3);

        if ((*sCPU).imdone[0] == 2) {
            break;
        }

        if ((*sCPU).imdone[0] == 3) {
            //copy the field arrays from the GPU to CPU memory
            cudaMemcpy((*sCPU).ExtOut, s.gridETime1, 2 * (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy((*sCPU).EkwOut, s.gridEFrequency1, 2 * (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

            (*sCPU).imdone[0] = 0;
        }
        (*(*sCPU).progressCounter)++;
    }

    
    //transform final result
    cufftExecZ2D(s.fftPlanZ2D, (cufftDoubleComplex*)s.gridEFrequency1, s.gridETime1);
    multiplyByConstantKernelD<<<(int)(s.Ngrid/16), 32, 0, s.CUDAStream >> >(s.gridETime1, 0.5);
    realToComplexKernel << <2*s.Nblock, s.Nthread, 0, s.CUDAStream >> > (s.gridETime1, s.workspace1);
    fftNormalizeKernel << <2*s.Nblock, s.Nthread, 0, s.CUDAStream >> > (s.workspace1, &s.propagationInts[0]);

    //copy the field arrays from the GPU to CPU memory
    cudaMemcpy((*sCPU).ExtOut, s.workspace1, 2*(*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    multiplyByConstantKernel << <(int)(s.Ngrid/16), 32, 0, s.CUDAStream >> > (s.workspace1, 2);
    getTotalSpectrum(sCPU, &s);
    realToComplexKernel << <2*s.Nblock, s.Nthread, 0, s.CUDAStream >> > (s.gridETime1, s.workspace1);
    fftNormalizeKernel << <2*s.Nblock, s.Nthread, 0, s.CUDAStream >> > (s.workspace1, &s.propagationInts[0]);
    multiplyByConstantKernel << <(int)(s.Ngrid/16), 32, 0, s.CUDAStream >> > (s.workspace1, 2);
    cufftExecZ2Z(s.fftPlan, (cufftDoubleComplex*)s.workspace1, (cufftDoubleComplex*)s.workspace1, CUFFT_FORWARD);
    cudaMemcpy((*sCPU).EkwOut, s.workspace1, 2*(*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    

    //Free GPU memory
    cudaFree(s.propagationInts);
    cudaFree(s.nonlinearSwitches);
    cudaFree(s.absorptionParameters);
    cudaFree(s.gridETime1); 
    cudaFree(s.workspace1);
    cudaFree(s.gridETemp1);
    cudaFree(s.gridPolarizationFrequency1);
    cudaFree(s.gridEFrequency1);
    cudaFree(s.gridPropagationFactor1);
    cudaFree(s.gridPropagationFactor1Rho1);
    cudaFree(s.gridRadialLaplacian1);
    cudaFree(s.firstDerivativeOperation);
    cudaFree(s.gridPolarizationFactor1);
    cudaFree(s.gridEFrequency1Next1);
    cudaFree(s.k1);
    cudaFree(s.gridPolarizationTime1);
    cudaFree(s.chi2Tensor);
    cudaFree(s.chi3Tensor);
    cudaFree(s.expGammaT);
    cudaFree(s.chiLinear1);
    cufftDestroy(s.fftPlan);
    cufftDestroy(s.fftPlanD2Z);
    cufftDestroy(s.fftPlanZ2D);
    cufftDestroy(s.polfftPlan);
    cufftDestroy(s.doublePolfftPlan);
    cudaFree(s.plasmaParameters);
    cudaFree(s.gridPlasmaCurrent1);
    cudaFree(s.gridPlasmaCurrentFrequency1);
    cudaFree(sP);
    cudaStreamDestroy(s.CUDAStream);


    //Free CPU memory
    free(gridPropagationFactor1CPU);
    free(gridPolarizationFactor1CPU);
    (*sCPU).imdone[0] = 1;
    return 0;
}

//function to run a RK4 time step
//stepNumber is the sub-step index, from 0 to 3
int runRK4Step(cudaParameterSet* sH, cudaParameterSet* sD, int stepNumber) {

    //operations involving FFT
    if ((*sH).isNonLinear || (*sH).isCylindric) {
        //perform inverse FFT to get time-space electric field
        cufftExecZ2D((*sH).fftPlanZ2D, (cufftDoubleComplex*)(*sH).workspace1C, (*sH).gridETime1);
        if ((*sH).isNonLinear) {
            nonlinearPolarizationKernel<<<(*sH).Nblock, (*sH).Nthread, 0, (*sH).CUDAStream>>>(sD);

            if ((*sH).isCylindric) {
                expandCylindricalBeam <<< (*sH).Nblock, (*sH).Nthread, 0, (*sH).CUDAStream >>> 
                    (sD, (*sH).gridPolarizationTime1, (*sH).gridPolarizationTime2);
                cufftExecD2Z((*sH).doublePolfftPlan, (double*)(*sH).gridRadialLaplacian1, (cufftDoubleComplex*)(*sH).gridPolarizationFrequency1);
            }
            else {
                cufftExecD2Z((*sH).polfftPlan, (*sH).gridPolarizationTime1, (cufftDoubleComplex*)(*sH).gridPolarizationFrequency1);
            }
        }
        if ((*sH).hasPlasma) {
            plasmaCurrentKernelPrep <<<(*sH).Nblock, (*sH).Nthread, 0, (*sH).CUDAStream >>> 
                (sD, (double*)(*sH).gridPlasmaCurrentFrequency1, (double*)(*sH).gridPlasmaCurrentFrequency2);
            plasmaCurrentKernel2 <<<(unsigned int)(*sH).Nspace/16, 16, 0, (*sH).CUDAStream >>> 
                (sD, (double*)(*sH).gridPlasmaCurrentFrequency1, (double*)(*sH).gridPlasmaCurrentFrequency2);
            
            if ((*sH).isCylindric) {
                expandCylindricalBeam <<< (*sH).Nblock, (*sH).Nthread, 0, (*sH).CUDAStream >>> 
                    (sD, (*sH).gridPlasmaCurrent1, (*sH).gridPlasmaCurrent2);
                cufftExecD2Z((*sH).doublePolfftPlan, (double*)(*sH).gridRadialLaplacian1, (cufftDoubleComplex*)(*sH).gridPlasmaCurrentFrequency1);
            }
            else {
                cufftExecD2Z((*sH).polfftPlan, (*sH).gridPlasmaCurrent1, (cufftDoubleComplex*)(*sH).gridPlasmaCurrentFrequency1);
            }

        }
        
        if ((*sH).isCylindric) {
            radialLaplacianKernel <<<(*sH).Nblock, (*sH).Nthread, 0, (*sH).CUDAStream >>> (sD);
            cufftExecD2Z((*sH).fftPlanD2Z, (*sH).gridRadialLaplacian1, (cufftDoubleComplex*)(*sH).workspace1C);
        }
    }

    //calculate k
    rkKernel <<<(*sH).NblockC, (*sH).Nthread, 0, (*sH).CUDAStream >>> (sD, stepNumber);
    return 0;
}

int prepareElectricFieldArrays(simulationParameterSet* s, cudaParameterSet *sc) {
    size_t i,j;
    double rB, zB, r, z; //r and z in the Beam and lab coordinates, respectively.
    double w0, wz, zR, Rz, phi; //Gaussian beam parameters
    double theta = 0; //rotation angle of the current beam
    double pulseSum = 0;
    std::complex<double> ne, no; //active refractive index;
    double f, w; //active frequency;
    double pulseEnergySum;
    std::complex<double> ko, specfac, specphase;
    double eps0 = 8.8541878128e-12; //vacuum permittivity
    
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


    //get the material phase
    double* materialCoefficientsCUDA;
    //NOTE TO SELF: add second phase material
    double* materialPhase1CUDA;
    double* materialPhase2CUDA;
    cudaMalloc((void**)&materialCoefficientsCUDA, 66 * sizeof(double));
    cudaMalloc((void**)&materialPhase1CUDA, (*s).Ntime * sizeof(double));
    cudaMalloc((void**)&materialPhase2CUDA, (*s).Ntime * sizeof(double));
    cudaMemcpy(materialCoefficientsCUDA, (*s).crystalDatabase[(*s).phaseMaterialIndex1].sellmeierCoefficients, 66 * sizeof(double), cudaMemcpyHostToDevice);
    materialPhaseKernel<<<(unsigned int)(*s).Ntime, 1, 0, (*sc).CUDAStream >> >((*s).fStep, (*s).Ntime, materialCoefficientsCUDA, (*s).frequency1, (*s).frequency2, (*s).phaseMaterialThickness1, (*s).phaseMaterialThickness2, materialPhase1CUDA, materialPhase2CUDA);
    double* materialPhase1 = (double*)calloc((*s).Ntime, sizeof(double));
    double* materialPhase2 = (double*)calloc((*s).Ntime, sizeof(double));
    cudaMemcpy(materialPhase1, materialPhase1CUDA, (*s).Ntime * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(materialPhase2, materialPhase2CUDA, (*s).Ntime * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(materialPhase2CUDA);
    cudaFree(materialPhase1CUDA);
    cudaFree(materialCoefficientsCUDA);



    for (i = 1; i < (*s).Ntime; i++) {
        f = i * (*s).fStep;
        if (i >= (*s).Ntime / 2) {
            f -= (*s).fStep * (*s).Ntime;
        }
        //f *= -1;
        w = TWOPI * (f - (*s).frequency1);
        
        //supergaussian pulse spectrum, if no input pulse specified
        specfac = pow((f - (*s).frequency1)/(*s).bandwidth1,(*s).sgOrder1);

        specphase = ii * ((*s).cephase1 + TWOPI*f * ((*s).delay1 - 0.5*(*s).tStep*(*s).Ntime) + 0.5 * (*s).gdd1 * w * w + (*s).tod1 * w * w * w/6.0 + materialPhase1[i]);
        specfac = exp(-specfac - specphase);

        if ((*s).field1IsAllocated) {
            specfac = (*s).loadedField1[i] * exp(-specphase);
        }



        for (j = 0; j < (*s).Nspace; j++) {
            ne = (*s).refractiveIndex1[i + (*s).Ntime * j];
            no = (*s).refractiveIndex2[i + (*s).Ntime * j];
            ko = TWOPI * no * f / LIGHTC;
            zR = PI * w0 * w0 * real(ne) * f / LIGHTC;
            if (f == 0) {
                zR = 1e3;
            }
            rB = ((*s).x01 - (*s).rStep * (j - (*s).Nspace / 2.0) + 0.25*(*s).rStep);
            r = rB * cos(theta) - zB * sin(theta);
            z = rB * sin(theta) + zB * cos(theta);
            
            wz = w0 * sqrt(1 + (z * z / (zR * zR)));
            Rz = z * (1. + (zR * zR / (z * z)));
            
            if (z == 0) {
                Rz = 1.0e15;
            }
            phi = atan(z / zR);
            Eb = (w0 / wz) * exp(ii * (real(ko) * (z-zB) + real(ko) * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
            Eb *= specfac;
            if (isnan(cModulusSquared(Eb)) || f<=0) {
                Eb = 0;
            }
            
            pulse1[i + (*s).Ntime * j] = polFactor1 * Eb;
            pulse1[i + (*s).Ntime * j + (*s).Ngrid] = polFactor2 * Eb;
            pulseSum += abs(r)*(real(ne)*cModulusSquared(pulse1[i + (*s).Ntime * j]) 
                + real(no)*cModulusSquared(pulse1[i + (*s).Ntime * j + (*s).Ngrid]));
        }
    }
    
    // copy the field and propagation grids to the GPU
    cudaMemcpy((*sc).workspace1, pulse1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).workspace2, &pulse1[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    // fft along spatial dimention to get Fourier space beam
    // will take place in three steps:
    // 2D fft (x,f)->(k,t), temporary intermediate state (could be optimized out later)
    // 1D fft (k,t)->(k,f), copied to Fourier space beam
    // 2D fft (k,f)->(x,t), copied to real space beam
    cufftPlan1d(&plan1, (int)(*sc).Ntime, CUFFT_Z2Z, (int)(*sc).Nspace);
    cufftSetStream(plan1, (*sc).CUDAStream);
    cufftPlan2d(&plan2, (int)(*sc).Nspace, (int)(*sc).Ntime, CUFFT_Z2Z);
    cufftSetStream(plan2, (*sc).CUDAStream);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_INVERSE);

    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace2, (cufftDoubleComplex*)(*sc).workspace2, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).workspace2, (cufftDoubleComplex*)(*sc).workspace2, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace2, (cufftDoubleComplex*)(*sc).workspace2, CUFFT_INVERSE);

    cudaMemcpy(pulse1, (*sc).workspace1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_INVERSE);
    cudaMemcpy(pulse1f, (*sc).workspace1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    cudaMemcpy(&pulse1[(*s).Ngrid], (*sc).workspace2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace2, (cufftDoubleComplex*)(*sc).workspace2, CUFFT_INVERSE);
    cudaMemcpy(&pulse1f[(*s).Ngrid], (*sc).workspace2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    //Copy the GPU grids to the CPU memory
    
    
    
    

    //normalize the pulse energy and set it to the input value
    pulseSum *= PI * LIGHTC * eps0;
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
        //f *= -1;
        w = TWOPI * (f - (*s).frequency2);

        //supergaussian pulse spectrum, if no input pulse specified
        specfac = pow((f - (*s).frequency2) / (*s).bandwidth2,(*s).sgOrder2);

        specphase = ii * ((*s).cephase2 + TWOPI*f * ((*s).delay2 - 0.5*(*s).tStep*(*s).Ntime) + 0.5*(*s).gdd2 * w * w + (*s).tod2 * w * w * w/6.0 + materialPhase2[i]);
        specfac = exp(-specfac - specphase);

        if ((*s).field2IsAllocated) {
            specfac = (*s).loadedField2[i] * exp(-specphase);
        }




        for (j = 0; j < (*s).Nspace; j++) {
            ne = (*s).refractiveIndex1[i + (*s).Ntime * j];
            no = (*s).refractiveIndex2[i + (*s).Ntime * j];
            ko = TWOPI * no * f / LIGHTC;
            zR = PI * w0 * w0 * real(ne) * f / LIGHTC;
            if (f == 0) {
                zR = 1e3;
            }
            rB = ((*s).x01 - (*s).rStep * (j - (*s).Nspace / 2.0) + 0.25 * (*s).rStep);
            r = rB * cos(theta) - zB * sin(theta);
            z = rB * sin(theta) + zB * cos(theta);

            wz = w0 * sqrt(1 + (z * z / (zR * zR)));
            Rz = z * (1. + (zR * zR / (z * z)));

            if (z == 0) {
                Rz = 1.0e15;
            }
            phi = atan(z / zR);
            Eb = (w0 / wz) * exp(ii * (real(ko) * (z - zB) + real(ko) * r * r / (2 * Rz) - phi) - r * r / (wz * wz));
            Eb *= specfac;
            if (isnan(cModulusSquared(Eb)) || f <= 0) {
                Eb = 0;
            }

            pulse2[i + (*s).Ntime * j] = polFactor1 * Eb;
            pulse2[i + (*s).Ntime * j + (*s).Ngrid] = polFactor2 * Eb;
            pulseSum += abs(r) * (real(ne) * cModulusSquared(pulse2[i + (*s).Ntime * j]) 
                + real(no) * cModulusSquared(pulse2[i + (*s).Ntime * j + (*s).Ngrid]));
        }
    }

    // copy the field and propagation grids to the GPU
    cudaMemcpy((*sc).workspace1, pulse2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy((*sc).workspace2, &pulse2[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    // fft along spatial dimention to get Fourier space beam
    // will take place in three steps:
    // 2D fft (x,f)->(k,t), temporary intermediate state (could be optimized out later)
    // 1D fft (k,t)->(k,f), copied to Fourier space beam
    // 2D fft (k,f)->(x,t), copied to real space beam

    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_INVERSE);

    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace2, (cufftDoubleComplex*)(*sc).workspace2, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).workspace2, (cufftDoubleComplex*)(*sc).workspace2, CUFFT_FORWARD);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace2, (cufftDoubleComplex*)(*sc).workspace2, CUFFT_INVERSE);

    cudaMemcpy(pulse2, (*sc).workspace1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_INVERSE);
    cudaMemcpy(pulse2f, (*sc).workspace1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    cudaMemcpy(&pulse2[(*s).Ngrid], (*sc).workspace2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cufftExecZ2Z(plan2, (cufftDoubleComplex*)(*sc).workspace2, (cufftDoubleComplex*)(*sc).workspace2, CUFFT_INVERSE);
    cudaMemcpy(&pulse2f[(*s).Ngrid], (*sc).workspace2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();


    //normalize the pulse energy and set it to the input value
    pulseSum *= PI * LIGHTC * eps0;
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
    cudaMemcpy((*sc).workspace1, (*s).Ekw, 2*(*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);


    //Adapt to the compact grid based on double-to-complex FFTs
    cufftExecZ2Z((*sc).fftPlan, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_INVERSE);
    fftNormalizeKernel << <2*(*sc).Nblock, (*sc).Nthread, 0, (*sc).CUDAStream >> > ((*sc).workspace1, (*sc).propagationInts);
    complexToRealKernel << <2*(*sc).Nblock, (*sc).Nthread, 0, (*sc).CUDAStream >> > ((*sc).workspace1, (*sc).gridETime1);
    cufftExecD2Z((*sc).fftPlanD2Z, (*sc).gridETime1, (*sc).gridEFrequency1);
    multiplyByConstantKernel<<<(unsigned int)(((*sc).Ntime/2 + 1) * ((*sc).Nspace/16)), 32, 0, (*sc).CUDAStream >> >((*sc).gridEFrequency1, 2.0);
    cudaDeviceSynchronize();
    cufftDestroy(plan1);
    cufftDestroy(plan2);

    free(pulse1);
    free(pulse2);
    free(pulse1f);
    free(pulse2f);
    free(materialPhase1);
    free(materialPhase2);
    return 0;
}
int applyFresnelLoss(simulationParameterSet* s, int materialIndex1, int materialIndex2) {
    cudaParameterSet sc;
    sc.Ntime = (*s).Ntime;
    sc.Nspace = (*s).Nspace;
    sc.dt = (*s).tStep;
    sc.dx = (*s).rStep;
    sc.fStep = (*s).fStep;
    sc.h = (*s).propagationStep;
    sc.Nsteps = (*s).Npropagation;
    sc.Ngrid = (*s).Ntime * (*s).Nspace;
    sc.axesNumber = (*s).axesNumber;
    sc.sellmeierType = (*s).sellmeierType;
    sc.f0 = (*s).frequency1;
    sc.Nthread = THREADS_PER_BLOCK;
    sc.Nblock = (int)((*s).Ngrid / THREADS_PER_BLOCK);
    sc.isCylindric = (*s).isCylindric;
    sc.isNonLinear = ((*s).nonlinearSwitches[0] + (*s).nonlinearSwitches[1]) > 0;
    cudaStreamCreate(&sc.CUDAStream);
    cufftPlan2d(&sc.fftPlan, (int)sc.Nspace, (int)sc.Ntime, CUFFT_Z2Z);
    double* sellmeierCoefficientsAugmentedCPU = (double*)calloc(66 + 8, sizeof(double));
    memcpy(sellmeierCoefficientsAugmentedCPU, (*s).crystalDatabase[materialIndex1].sellmeierCoefficients, 66 * (sizeof(double)));
    sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
    sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
    sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
    sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
    sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
    sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
    sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
    double* sellmeierCoefficients1;
    double* sellmeierCoefficients2;
    cudaMalloc(&sellmeierCoefficients1, 74 * sizeof(double));
    cudaMalloc(&sellmeierCoefficients2, 74 * sizeof(double));
    cudaMemcpy(sellmeierCoefficients1, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    memcpy(sellmeierCoefficientsAugmentedCPU, (*s).crystalDatabase[materialIndex2].sellmeierCoefficients, 66 * (sizeof(double)));
    sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
    sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
    sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
    sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
    sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
    sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
    sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
    cudaMemcpy(sellmeierCoefficients2, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);

    size_t propagationIntsCPU[4] = { sc.Ngrid, sc.Ntime, sc.Nspace, (sc.Ntime / 2 + 1) };
    cudaMalloc((void**)&sc.propagationInts, sizeof(size_t) * 4);
    cudaMemcpy(sc.propagationInts, propagationIntsCPU, 4 * sizeof(size_t), cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();
    cudaMalloc(&sc.gridEFrequency1, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc(&sc.gridEFrequency2, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc(&sc.gridETime1, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc(&sc.gridETime2, (*s).Ngrid * sizeof(cuDoubleComplex));

    cudaMemcpy(sc.gridEFrequency1, (*s).EkwOut, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(sc.gridEFrequency2, &(*s).EkwOut[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    applyFresnelLossKernel<<<sc.Nblock, sc.Nthread, 0, sc.CUDAStream>>>(sellmeierCoefficients1, sellmeierCoefficients2, sc);

    //transform final result
    fixnanKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sc.gridEFrequency1);
    fixnanKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sc.gridEFrequency2);
    cufftExecZ2Z(sc.fftPlan, (cufftDoubleComplex*)sc.gridEFrequency1, (cufftDoubleComplex*)sc.gridETime1, CUFFT_INVERSE);
    cufftExecZ2Z(sc.fftPlan, (cufftDoubleComplex*)sc.gridEFrequency2, (cufftDoubleComplex*)sc.gridETime2, CUFFT_INVERSE);
    //fftNormalizeKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sc.gridETime1, sc.propagationInts);
    //fftNormalizeKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sc.gridETime2, sc.propagationInts);

    //copy the field arrays from the GPU to CPU memory
    cudaMemcpy((*s).ExtOut, sc.gridETime1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*s).EkwOut, sc.gridEFrequency1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*s).ExtOut[sc.Ngrid], sc.gridETime2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*s).EkwOut[sc.Ngrid], sc.gridEFrequency2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);


    free(sellmeierCoefficientsAugmentedCPU);
    
    cufftDestroy(sc.fftPlan);
    cudaFree(sc.gridEFrequency1);
    cudaFree(sc.gridEFrequency2);
    cudaFree(sc.gridETime1);
    cudaFree(sc.gridETime2);
    cudaFree(sellmeierCoefficients1);
    cudaFree(sellmeierCoefficients2);
    cudaStreamDestroy(sc.CUDAStream);
    return 0;
}


int applyLinearPropagation(simulationParameterSet* s, int materialIndex, double thickness) {
    cudaParameterSet sc;
    sc.Ntime = (*s).Ntime;
    sc.Nspace = (*s).Nspace;
    sc.dt = (*s).tStep;
    sc.dx = (*s).rStep;
    sc.fStep = (*s).fStep;
    sc.h = (*s).propagationStep;
    sc.Nsteps = (*s).Npropagation;
    sc.Ngrid = (*s).Ntime * (*s).Nspace;
    sc.axesNumber = (*s).axesNumber;
    sc.sellmeierType = (*s).sellmeierType;
    sc.f0 = (*s).frequency1;
    sc.Nthread = THREADS_PER_BLOCK;
    sc.Nblock = (int)((*s).Ngrid / THREADS_PER_BLOCK);
    sc.isCylindric = (*s).isCylindric;
    sc.isNonLinear = ((*s).nonlinearSwitches[0] + (*s).nonlinearSwitches[1]) > 0;
    cudaStreamCreate(&sc.CUDAStream);
    cufftPlan2d(&sc.fftPlan, (int)sc.Nspace, (int)sc.Ntime, CUFFT_Z2Z);
    double* sellmeierCoefficientsAugmentedCPU = (double*)calloc(66 + 8, sizeof(double));
    memcpy(sellmeierCoefficientsAugmentedCPU, (*s).crystalDatabase[materialIndex].sellmeierCoefficients, 66 * (sizeof(double)));
    sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
    sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
    sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
    sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
    sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
    sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
    sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
    double* sellmeierCoefficients;

    cudaMalloc(&sellmeierCoefficients, 74 * sizeof(double));

    cudaMemcpy(sellmeierCoefficients, sellmeierCoefficientsAugmentedCPU, (66 + 8) * sizeof(double), cudaMemcpyHostToDevice);


    size_t propagationIntsCPU[4] = { sc.Ngrid, sc.Ntime, sc.Nspace, (sc.Ntime / 2 + 1) };
    cudaMalloc((void**)&sc.propagationInts, sizeof(size_t) * 4);
    cudaMemcpy(sc.propagationInts, propagationIntsCPU, 4 * sizeof(size_t), cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();
    cudaMalloc(&sc.gridEFrequency1, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc(&sc.gridEFrequency2, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc(&sc.gridETime1, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc(&sc.gridETime2, (*s).Ngrid * sizeof(cuDoubleComplex));

    cudaMemcpy(sc.gridEFrequency1, (*s).EkwOut, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(sc.gridEFrequency2, &(*s).EkwOut[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    applyLinearPropagationKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sellmeierCoefficients, thickness, sc);

    //transform final result
    fixnanKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sc.gridEFrequency1);
    fixnanKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sc.gridEFrequency2);
    cufftExecZ2Z(sc.fftPlan, (cufftDoubleComplex*)sc.gridEFrequency1, (cufftDoubleComplex*)sc.gridETime1, CUFFT_INVERSE);
    cufftExecZ2Z(sc.fftPlan, (cufftDoubleComplex*)sc.gridEFrequency2, (cufftDoubleComplex*)sc.gridETime2, CUFFT_INVERSE);

    //copy the field arrays from the GPU to CPU memory
    cudaMemcpy((*s).ExtOut, sc.gridETime1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*s).EkwOut, sc.gridEFrequency1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*s).ExtOut[sc.Ngrid], sc.gridETime2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*s).EkwOut[sc.Ngrid], sc.gridEFrequency2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);


    free(sellmeierCoefficientsAugmentedCPU);

    cufftDestroy(sc.fftPlan);
    cudaFree(sc.gridEFrequency1);
    cudaFree(sc.gridEFrequency2);
    cudaFree(sc.gridETime1);
    cudaFree(sc.gridETime2);
    cudaFree(sellmeierCoefficients);
    cudaStreamDestroy(sc.CUDAStream);
    return 0;
}

int preparePropagation2DCartesian(simulationParameterSet* s, cudaParameterSet sc) {
    //recycle allocated device memory for the grids needed
    double* sellmeierCoefficients = (double*)sc.gridEFrequency1Next1;
    sc.ne = sc.workspace1;
    sc.no = sc.workspace2;

    double* referenceFrequencies;
    cudaMalloc(&referenceFrequencies, 7 * sizeof(double));
    cudaMemcpy(referenceFrequencies, (*s).crystalDatabase[(*s).materialIndex].nonlinearReferenceFrequencies, 7 * sizeof(double), cudaMemcpyHostToDevice);

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
    prepareCartesianGridsKernel <<<sc.Nblock, sc.Nthread, 0, sc.CUDAStream >>> (sellmeierCoefficients, sc);
    millersRuleNormalizationKernel << <1, 1, 0, sc.CUDAStream >> > (sc, sellmeierCoefficients, referenceFrequencies);
    cudaDeviceSynchronize();

    //copy the retrieved refractive indicies to the cpu
    cudaMemcpy((*s).refractiveIndex1, sc.ne, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*s).refractiveIndex2, sc.no, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    
    //clean up
    cudaMemset(sc.gridEFrequency1Next1, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.gridEFrequency1Next2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    free(sellmeierCoefficientsAugmentedCPU);
    cudaFree(referenceFrequencies);
    return 0;
}

int preparePropagation3DCylindric(simulationParameterSet* s, cudaParameterSet sc) {
    //recycle allocated device memory for the grids needed
    double* sellmeierCoefficients = (double*)sc.gridEFrequency1Next1;
    sc.ne = sc.workspace1;
    sc.no = sc.workspace2;
    double* referenceFrequencies;
    cudaMalloc(&referenceFrequencies, 7*sizeof(double));
    cudaMemcpy(referenceFrequencies, (*s).crystalDatabase[(*s).materialIndex].nonlinearReferenceFrequencies, 7 * sizeof(double), cudaMemcpyHostToDevice);

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
    cudaDeviceSynchronize();
    //prepare the propagation grids
    
    prepareCylindricGridsKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sellmeierCoefficients, sc);
    millersRuleNormalizationKernel<<<1, 1, 0, sc.CUDAStream>>>(sc, sellmeierCoefficients, referenceFrequencies);
    cudaDeviceSynchronize();

    //copy the retrieved refractive indicies to the cpu
    cudaMemcpy((*s).refractiveIndex1, sc.ne, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy((*s).refractiveIndex2, sc.no, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    //clean up
    cudaMemset(sc.gridEFrequency1Next1, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMemset(sc.gridEFrequency1Next2, 0, (*s).Ngrid * sizeof(cuDoubleComplex));
    cudaFree(referenceFrequencies);
    free(sellmeierCoefficientsAugmentedCPU);
    return 0;
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
    double l = 1e6* LIGHTC / f; //wavelength in microns
    double ls = l * l;
    std::complex<double> ii(0, 1);
    double omega = TWOPI*abs(f);
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
int loadReferenceSpectrum(char* spectrumPath, simulationParameterSet* sCPU) {
    FILE* fp = fopen(spectrumPath, "r");
    if (fp == NULL) {
        printf("Could not read reference file\r\n");
        return 1;
    }
    size_t maxFileSize = 16384;
    size_t currentRow = 0;
    double c = 1e9 * LIGHTC;
    double* loadedWavelengths = (double*)calloc(8192, sizeof(double));
    double* loadedFrequencies = (double*)calloc(8192, sizeof(double));
    double* loadedIntensities = (double*)calloc(8192, sizeof(double));
    double maxWavelength = 0;
    double minWavelength = 0;
    if (fp == NULL) {
        free(loadedWavelengths);
        free(loadedIntensities);
        free(loadedFrequencies);
        return -1;
    }

    while (fscanf(fp, "%lf %lf", &loadedWavelengths[currentRow], &loadedIntensities[currentRow]) == 2 && currentRow < maxFileSize) {
        if (currentRow == 0) {
            maxWavelength = loadedWavelengths[currentRow];
            minWavelength = loadedWavelengths[currentRow];
        }
        else {
            maxWavelength = max(maxWavelength, loadedWavelengths[currentRow]);
            minWavelength = min(minWavelength, loadedWavelengths[currentRow]);
        }
        //rescale to frequency spacing
        loadedIntensities[currentRow] *= loadedWavelengths[currentRow] * loadedWavelengths[currentRow];
        loadedFrequencies[currentRow] = c / loadedWavelengths[currentRow];
        currentRow++;
    }
    size_t sizeData = currentRow-1;
    size_t i, j;

    double maxFrequency = c / minWavelength;
    double minFrequency = c / maxWavelength;
    double currentFrequency = 0;
    double df;
    
    for(i = 1; i < (*sCPU).Ntime; i++) {
        currentFrequency = i * (*sCPU).fStep;
        if ((currentFrequency > minFrequency) && (currentFrequency < maxFrequency)) {
            //find the first frequency greater than the current value
            j = sizeData-1;
            while((loadedFrequencies[j] <= currentFrequency) && (j > 2)) {
                j--;
            }
            df = loadedFrequencies[j] - loadedFrequencies[j - 1];
            (*sCPU).fittingArray[i] = 
                (loadedIntensities[j-1] * (loadedFrequencies[j] - currentFrequency) 
                    + loadedIntensities[j] * (currentFrequency - loadedFrequencies[j-1])) / df; //linear interpolation
        }
    }
    
    fclose(fp);
    free(loadedWavelengths);
    free(loadedIntensities);
    free(loadedFrequencies);
    return 0;
}
int loadFrogSpeck(char* frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel, int fieldIndex) {
    FILE* fp;
    int maxFileSize = 16384;
    double wavelength, R, phi, complexX, complexY, f, f0, f1, fmax;
    int i, k0, k1;
    double c = 1e9 * LIGHTC; //for conversion of wavelength in nm to frequency
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
            f -= fStep * (Ntime-1);
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


//Rotate the field on the GPU
//Allocates memory and copies from CPU, then copies back to CPU and deallocates
// - inefficient but the general principle is that only the CPU memory is preserved
// after simulations finish... and this only runs at the end of the simulation
int rotateField(simulationParameterSet *s, double rotationAngle) {
    cudaParameterSet sc;
    cuDoubleComplex* Ein1, * Eout1, * Ein2, * Eout2;
    double *spectrum;
    cudaMalloc((void**)&Ein1, 2*(*s).Ngrid * sizeof(cuDoubleComplex));
    cudaMalloc((void**)&Eout1, 2*(*s).Ngrid * sizeof(cuDoubleComplex));
    Ein2 = &Ein1[(*s).Ngrid];
    Eout2 = &Eout1[(*s).Ngrid];

    cudaMalloc((void**)&spectrum, (*s).Ntime * 3 * sizeof(double));
    cudaMalloc((void**)&sc.propagationInts, 4 * sizeof(size_t));
    unsigned int Nthread = THREADS_PER_BLOCK;
    unsigned int Nblock = (unsigned int)((*s).Ngrid / THREADS_PER_BLOCK);

    cudaMemcpy(Ein1, (*s).EkwOut, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(Ein2, &(*s).EkwOut[(*s).Ngrid], (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    rotateFieldKernel<<<Nblock, Nthread>>>(Ein1, Ein2, Eout1, Eout2, rotationAngle);

    cudaMemcpy((*s).EkwOut, Eout1, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*s).EkwOut[(*s).Ngrid], Eout2, (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    //recalculate spectra
    //first populate only the necessary parts of the struct and make the FFTplan
    sc.gridETemp1 = Eout1;
    sc.gridETemp2 = Eout2;
    sc.workspace1 = Ein1;
    sc.workspace2 = Ein2;
    sc.Nthread = THREADS_PER_BLOCK;
    sc.Nblock = (int)((*s).Ngrid / THREADS_PER_BLOCK);
    sc.gridPolarizationTime1 = spectrum;
    size_t propagationIntsCPU[4] = { (*s).Ngrid, (*s).Ntime, (*s).Nspace, ((*s).Ntime / 2 + 1) };
    cudaMemcpy(sc.propagationInts, propagationIntsCPU, 4 * sizeof(size_t), cudaMemcpyHostToDevice);
    cudaSetDevice((*s).assignedGPU);
    cudaStreamCreate(&sc.CUDAStream);
    size_t workSize;
    int cufftSizes1[] = { (int)(*s).Nspace, (int)(*s).Ntime };
    cufftCreate(&sc.fftPlan);
    cufftGetSizeMany(sc.fftPlan, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_Z2Z, 2, &workSize);
    cufftMakePlanMany(sc.fftPlan, 2, cufftSizes1, NULL, NULL, NULL, NULL, NULL, NULL, CUFFT_Z2Z, 2, &workSize);
    cufftSetStream(sc.fftPlan, sc.CUDAStream);

    //fft the rotated fields and normalize
    cufftExecZ2Z(sc.fftPlan, (cufftDoubleComplex*)Eout1, (cufftDoubleComplex*)Ein1, CUFFT_INVERSE);
    fftNormalizeKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sc.workspace1, sc.propagationInts);
    fftNormalizeKernel << <sc.Nblock, sc.Nthread, 0, sc.CUDAStream >> > (sc.workspace2, sc.propagationInts);
    cufftDestroy(sc.fftPlan);

    //update the CPU fields
    cudaDeviceSynchronize();
    cudaMemcpy((*s).ExtOut, Ein1, 2 * (*s).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);

    //new spectrum
    getTotalSpectrum(s, &sc);

    cudaFree(sc.propagationInts);
    cudaFree(Ein1);
    cudaFree(Eout1);
    cudaFree(spectrum);
    cudaStreamDestroy(sc.CUDAStream);
    return 0;
}

//calculates the squard modulus of a complex number, under the assumption that the
//machine's complex number format is interleaved doubles.
//c forced to run in c++ for nostalgia reasons
double cModulusSquared(std::complex<double>complexNumber) {
    double* xy = (double*)&complexNumber;
    return xy[0] * xy[0] + xy[1] * xy[1];
}

int allocateGrids(simulationParameterSet* sCPU) {
    (*sCPU).loadedField1 = (std::complex<double>*)calloc((*sCPU).Ntime, sizeof(std::complex<double>));
    (*sCPU).loadedField2 = (std::complex<double>*)calloc((*sCPU).Ntime, sizeof(std::complex<double>));

    (*sCPU).Ext = (std::complex<double>*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims * (*sCPU).Nsims2, sizeof(std::complex<double>));
    (*sCPU).Ekw = (std::complex<double>*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims * (*sCPU).Nsims2, sizeof(std::complex<double>));

    (*sCPU).ExtOut = (std::complex<double>*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims * (*sCPU).Nsims2, sizeof(std::complex<double>));
    (*sCPU).EkwOut = (std::complex<double>*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims * (*sCPU).Nsims2, sizeof(std::complex<double>));

    (*sCPU).refractiveIndex1 = (std::complex<double>*)calloc((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2, sizeof(std::complex<double>));
    (*sCPU).refractiveIndex2 = (std::complex<double>*)calloc((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2, sizeof(std::complex<double>));
    (*sCPU).deffTensor = (double*)calloc(9 * (*sCPU).Nsims * (*sCPU).Nsims2, sizeof(double));
    (*sCPU).totalSpectrum = (double*)calloc((*sCPU).Nsims * (*sCPU).Nsims2 * (*sCPU).Ntime * 3, sizeof(double));
    (*sCPU).imdone = (int*)calloc((*sCPU).Nsims * (*sCPU).Nsims2, sizeof(int));
    return 0;
}


int readCrystalDatabase(crystalEntry* db) {
    int i = 0;
    double* fd;
    FILE* fp;
    fp = fopen("CrystalDatabase.txt", "r");
    if (fp == NULL) {
        return -2;
    }
    wchar_t lineBuffer[MAX_LOADSTRING] = { 0 };

    //read the entries line
    int readErrors = 0;

    while (readErrors == 0 && !feof(fp) && i < MAX_LOADSTRING) {
        readErrors += 0 != fwscanf(fp, L"Name:\n");
        fgetws(db[i].crystalNameW, 256, fp);
        readErrors += 1 != fwscanf(fp, L"Type:\n%d\n", &db[i].axisType);
        readErrors += 1 != fwscanf(fp, L"Sellmeier equation:\n%d\n", &db[i].sellmeierType);
        fd = &db[i].sellmeierCoefficients[0];
        readErrors += 22 != fwscanf(fp, L"1st axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
            &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
        fd = &db[i].sellmeierCoefficients[22];
        readErrors += 22 != fwscanf(fp, L"2nd axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
            &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
        fd = &db[i].sellmeierCoefficients[44];
        readErrors += 22 != fwscanf(fp, L"3rd axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
            &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
        readErrors += 0 != fwscanf(fp, L"Sellmeier reference:\n");
        fgetws(db[i].sellmeierReference, 512, fp);
        readErrors += 1 != fwscanf(fp, L"chi2 type:\n%d\n", &db[i].nonlinearSwitches[0]);
        fd = &db[i].d[0];
        readErrors += 18 != fwscanf(fp, L"d:\n%lf %lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf %lf\n",
            &fd[0], &fd[3], &fd[6], &fd[9], &fd[12], &fd[15],
            &fd[1], &fd[4], &fd[7], &fd[10], &fd[13], &fd[16],
            &fd[2], &fd[5], &fd[8], &fd[11], &fd[14], &fd[17]);
        readErrors += 0 != fwscanf(fp, L"d reference:\n");
        fgetws(db[i].dReference, 512, fp);
        readErrors += 1 != fwscanf(fp, L"chi3 type:\n%d\nchi3:\n", &db[i].nonlinearSwitches[1]);
        //handle chi3 in a flexible way to avoid making the user write 81 zeroes when not needed
        fd = &db[i].chi3[0];
        memset(fd, 0, 81 * sizeof(double));
        if (db[i].nonlinearSwitches[1] != 1 && db[i].nonlinearSwitches[1] != 2) {
            fgetws(lineBuffer, MAX_LOADSTRING, fp);
            fgetws(lineBuffer, MAX_LOADSTRING, fp);
            fgetws(lineBuffer, MAX_LOADSTRING, fp);
        }
        else if (db[i].nonlinearSwitches[1] == 1) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 27; k++) {
                    readErrors += 1 != fwscanf(fp, L"%lf", &fd[j + 3 * k]);
                }
                fwscanf(fp, L"\n");
            }
        }
        else if (db[i].nonlinearSwitches[1] == 2) {
            readErrors += 1 != fwscanf(fp, L"%lf", &fd[0]);
            fgetws(lineBuffer, MAX_LOADSTRING, fp);
            fgetws(lineBuffer, MAX_LOADSTRING, fp);
            fgetws(lineBuffer, MAX_LOADSTRING, fp);
        }
        readErrors += 0 != fwscanf(fp, L"chi3 reference:\n");
        fgetws(db[i].chi3Reference, 512, fp);
        readErrors += 0 != fwscanf(fp, L"Spectral file:\n");
        fgetws(db[i].spectralFile, 512, fp);
        fd = db[i].nonlinearReferenceFrequencies;
        readErrors += 0 != fwscanf(fp, L"Nonlinear reference frequencies:\n");
        readErrors += 7 != fwscanf(fp, L"%lf %lf %lf %lf %lf %lf %lf\n",
            &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6]);
        readErrors += 0 != fwscanf(fp, L"~~~crystal end~~~\n");
        if (readErrors == 0) i++;
    }
    db[0].numberOfEntries = i;
    fclose(fp);

    return i;
}

int readSequenceString(simulationParameterSet* sCPU) {
    //read the sequence string (if there is one), convert it into an array if it exists
    char sequenceString[MAX_LOADSTRING];
    strcpy(sequenceString, (*sCPU).sequenceString);
char* tokToken = strtok(sequenceString, ";");
int sequenceCount = sscanf(sequenceString, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
    &(*sCPU).sequenceArray[0], &(*sCPU).sequenceArray[1], &(*sCPU).sequenceArray[2],
    &(*sCPU).sequenceArray[3], &(*sCPU).sequenceArray[4], &(*sCPU).sequenceArray[5],
    &(*sCPU).sequenceArray[6], &(*sCPU).sequenceArray[7], &(*sCPU).sequenceArray[8],
    &(*sCPU).sequenceArray[9], &(*sCPU).sequenceArray[10]);

tokToken = strtok(NULL, ";");
int lastread = sequenceCount;
while (tokToken != NULL && lastread == 11) {
    lastread = sscanf(tokToken, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &(*sCPU).sequenceArray[sequenceCount], &(*sCPU).sequenceArray[sequenceCount + 1],
        &(*sCPU).sequenceArray[sequenceCount + 2], &(*sCPU).sequenceArray[sequenceCount + 3],
        &(*sCPU).sequenceArray[sequenceCount + 4], &(*sCPU).sequenceArray[sequenceCount + 5],
        &(*sCPU).sequenceArray[sequenceCount + 6], &(*sCPU).sequenceArray[sequenceCount + 7],
        &(*sCPU).sequenceArray[sequenceCount + 8], &(*sCPU).sequenceArray[sequenceCount + 9],
        &(*sCPU).sequenceArray[sequenceCount + 10]);
    if (lastread > 0) {
        sequenceCount += lastread;
    }
    tokToken = strtok(NULL, ";");
}
(*sCPU).Nsequence = sequenceCount / 11;
(*sCPU).isInSequence = ((*sCPU).Nsequence > 0);

if (!(*sCPU).isInSequence) {
    char nopeString[] = "None.";
    strcpy((*sCPU).sequenceString, nopeString);
}
return 0;
}

int readFittingString(simulationParameterSet* sCPU) {
    //read the fitting string (if there is one), convert it into an array if it exists
    char fittingString[MAX_LOADSTRING];
    double ROIbegin;
    double ROIend;
    strcpy(fittingString, (*sCPU).fittingString);
    char* tokToken = strtok(fittingString, ";");
    bool paramsRead = (4 == sscanf(fittingString, "%lf %lf %lf %d",
        &ROIbegin, &ROIend, &(*sCPU).fittingPrecision, &(*sCPU).fittingMaxIterations));
    (*sCPU).fittingROIstart = (size_t)(ROIbegin / (*sCPU).fStep);
    (*sCPU).fittingROIstop = (size_t)min(ROIend / (*sCPU).fStep, (*sCPU).Ntime / 2);
    (*sCPU).fittingROIsize = min(max(1, (*sCPU).fittingROIstop - (*sCPU).fittingROIstart), (*sCPU).Ntime / 2);
    int fittingCount = 0;
    tokToken = strtok(NULL, ";");
    int lastread = 3;
    while (tokToken != NULL && lastread == 3) {
        lastread = sscanf(tokToken, "%lf %lf %lf",
            &(*sCPU).fittingArray[fittingCount], &(*sCPU).fittingArray[fittingCount + 1],
            &(*sCPU).fittingArray[fittingCount + 2]);
        if (lastread > 0) {
            fittingCount += lastread;
        }
        tokToken = strtok(NULL, ";");
    }
    (*sCPU).Nfitting = fittingCount / 3;
    (*sCPU).isInFittingMode = (((*sCPU).Nfitting) > 0 && paramsRead);

    if (!(*sCPU).isInFittingMode) {
        char nopeString[] = "None.";
        strcpy((*sCPU).fittingString, nopeString);
    }
    return 0;
}

int configureBatchMode(simulationParameterSet* sCPU) {
    size_t i,j, currentRow;
    if ((*sCPU).batchIndex == 0 || (*sCPU).Nsims == 1) {
        return 0;
    }

    //pointers to values that can be scanned in batch mode
    double* targets[36] = { 0,
        &(*sCPU).pulseEnergy1, &(*sCPU).pulseEnergy2, &(*sCPU).frequency1, &(*sCPU).frequency2,
        &(*sCPU).bandwidth1, &(*sCPU).bandwidth2, &(*sCPU).cephase1, &(*sCPU).cephase2,
        &(*sCPU).delay1, &(*sCPU).delay2, &(*sCPU).gdd1, &(*sCPU).gdd2,
        &(*sCPU).tod1, &(*sCPU).tod2, &(*sCPU).phaseMaterialThickness1, &(*sCPU).phaseMaterialThickness2,
        &(*sCPU).beamwaist1, &(*sCPU).beamwaist2,
        &(*sCPU).x01, &(*sCPU).x02, &(*sCPU).z01, &(*sCPU).z02,
        &(*sCPU).propagationAngle1, &(*sCPU).propagationAngle2, &(*sCPU).polarizationAngle1, &(*sCPU).polarizationAngle2,
        &(*sCPU).circularity1, &(*sCPU).circularity2, &(*sCPU).crystalTheta, &(*sCPU).crystalPhi,
        &(*sCPU).nonlinearAbsorptionStrength, &(*sCPU).drudeGamma, &(*sCPU).effectiveMass, &(*sCPU).crystalThickness,
        &(*sCPU).propagationStep };

    //multipliers to the Batch end value from the interface
    // (e.g. frequency in THz requires 1e12 multiplier)
    double multipliers[36] = { 0,
        1, 1, 1e12, 1e12,
        1e12, 1e12, PI, PI,
        1e-15, 1e-15, 1e-30, 1e-30,
        1e-45, 1e-45, 1e-6, 1e-6,
        1e-6, 1e-6,
        1e-6, 1e-6, 1e-6, 1e-6,
        DEG2RAD, DEG2RAD, DEG2RAD, DEG2RAD,
        1, 1, DEG2RAD, DEG2RAD,
        1, 1e12, 1, 1e-6,
        1e-9 };

    //Configure the struct array if in a batch
    for (i = 0; i < (*sCPU).Nsims2; i++) {
        currentRow = i * (*sCPU).Nsims;
        
        if (currentRow > 0) {
            memcpy(&sCPU[currentRow], sCPU, sizeof(simulationParameterSet));
        }
        if ((*sCPU).Nsims2 > 1) {
            *((double*)((simulationParameterSet*)targets[(*sCPU).batchIndex2] + currentRow)) +=
                i * (multipliers[(*sCPU).batchIndex2] * (*sCPU).batchDestination2 - *targets[(*sCPU).batchIndex2])
                / ((*sCPU).Nsims2 - 1);
        }
        
        for (j = 0; j < (*sCPU).Nsims; j++) {
            if (j > 0) {
                memcpy(&sCPU[j + currentRow], &sCPU[currentRow], sizeof(simulationParameterSet));
            }

            if ((*sCPU).deffTensor != NULL) {
                sCPU[j + currentRow].deffTensor = &(*sCPU).deffTensor[9 * (j + currentRow)];;
            }

            sCPU[j + currentRow].Ext = &(*sCPU).Ext[(j + currentRow) * (*sCPU).Ngrid * 2];
            sCPU[j + currentRow].Ekw = &(*sCPU).Ekw[(j + currentRow) * (*sCPU).Ngrid * 2];
            sCPU[j + currentRow].ExtOut = &(*sCPU).ExtOut[(j + currentRow) * (*sCPU).Ngrid * 2];
            sCPU[j + currentRow].EkwOut = &(*sCPU).EkwOut[(j + currentRow) * (*sCPU).Ngrid * 2];
            sCPU[j + currentRow].totalSpectrum = &(*sCPU).totalSpectrum[(j + currentRow) * (*sCPU).Ntime * 3];
            sCPU[j + currentRow].isFollowerInSequence = FALSE;

            // To add new modes, append values to the two arrays above, and to the combobox in the UI.
            // Cast the pointer to the original value to a pointer to a struct, 
            // increment, recast to a pointer to double and resolve then add j times the scan step size.
            *((double*)((simulationParameterSet*)targets[(*sCPU).batchIndex] + (j + currentRow))) +=
                j * (multipliers[(*sCPU).batchIndex] * (*sCPU).batchDestination - *targets[(*sCPU).batchIndex])
                / ((*sCPU).Nsims - 1);
        }
    }
    
    return 0;
}
int skipFileUntilCharacter(FILE* fstream, char target) {
    int c = 0;
    while (c != EOF && c != target) {
        c = fgetc(fstream);
    }
    return 0;
}
int readInputParametersFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, char* filePath) {
    FILE* textfile;
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
    int readValueCount = 0;
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).pulseEnergy1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).pulseEnergy2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).frequency1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).frequency2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).bandwidth1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).bandwidth2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%i", &(*sCPU).sgOrder1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%i", &(*sCPU).sgOrder2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).cephase1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).cephase2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).delay1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).delay2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).gdd1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).gdd2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).tod1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).tod2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%d", &(*sCPU).phaseMaterialIndex1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%d", &(*sCPU).phaseMaterialIndex2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).phaseMaterialThickness1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).phaseMaterialThickness2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).beamwaist1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).beamwaist2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).x01);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).x02);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).z01);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).z02);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).propagationAngle1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).propagationAngle2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).polarizationAngle1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).polarizationAngle2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).circularity1);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).circularity2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%i", &(*sCPU).materialIndex);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%i", &(*sCPU).materialIndexAlternate);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).crystalTheta);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).crystalPhi);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).spatialWidth);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).spatialHeight);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).rStep);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).timeSpan);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).tStep);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).crystalThickness);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).propagationStep);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).nonlinearAbsorptionStrength);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).bandGapElectronVolts);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).effectiveMass);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).drudeGamma);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%i", &(*sCPU).symmetryType);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%i", &(*sCPU).batchIndex);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).batchDestination);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lli", &(*sCPU).Nsims);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%i", &(*sCPU).batchIndex2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lf", &(*sCPU).batchDestination2);
    skipFileUntilCharacter(textfile, ':');
    readValueCount += fscanf(textfile, "%lli", &(*sCPU).Nsims2);
    readValueCount += fscanf(textfile, "\nSequence: ");
    fgets((*sCPU).sequenceString, MAX_LOADSTRING, textfile);
    readValueCount += fscanf(textfile, "Fitting: ");
    fgets((*sCPU).fittingString, MAX_LOADSTRING, textfile);
    readValueCount += fscanf(textfile, "Fitting mode : %i\n", &(*sCPU).fittingMode);
    readValueCount += fscanf(textfile, "Output base path: ");
    fgets((*sCPU).outputBasePath, MAX_LOADSTRING, textfile);
    readValueCount += fscanf(textfile, "Field 1 from file type: %i\nField 2 from file type: %i\n",
        &(*sCPU).pulse1FileType, &(*sCPU).pulse2FileType);
    readValueCount += fscanf(textfile, "Field 1 file path: ");
    fgets((*sCPU).field1FilePath, MAX_LOADSTRING, textfile);
    readValueCount += fscanf(textfile, "Field 2 file path: ");
    fgets((*sCPU).field2FilePath, MAX_LOADSTRING, textfile);
    readValueCount += fscanf(textfile, "Fitting reference file path: ");
    fgets((*sCPU).fittingPath, MAX_LOADSTRING, textfile);
    
    removeCharacterFromString((*sCPU).field1FilePath, MAX_LOADSTRING, '\r');
    removeCharacterFromString((*sCPU).field1FilePath, MAX_LOADSTRING, '\n');
    removeCharacterFromString((*sCPU).field2FilePath, MAX_LOADSTRING, '\r');
    removeCharacterFromString((*sCPU).field2FilePath, MAX_LOADSTRING, '\n');
    removeCharacterFromString((*sCPU).fittingPath, MAX_LOADSTRING, '\r');
    removeCharacterFromString((*sCPU).fittingPath, MAX_LOADSTRING, '\n');
    removeCharacterFromString((*sCPU).fittingString, MAX_LOADSTRING, '\r');
    removeCharacterFromString((*sCPU).fittingString, MAX_LOADSTRING, '\n');
    removeCharacterFromString((*sCPU).sequenceString, MAX_LOADSTRING, '\r');
    removeCharacterFromString((*sCPU).sequenceString, MAX_LOADSTRING, '\n');

    //derived parameters and cleanup:
    (*sCPU).sellmeierType = 0;
    (*sCPU).axesNumber = 0;
    (*sCPU).Ntime = (size_t)round((*sCPU).timeSpan / (*sCPU).tStep);
    (*sCPU).Nspace = (size_t)round((*sCPU).spatialWidth / (*sCPU).rStep);
    (*sCPU).Ngrid = (*sCPU).Ntime * (*sCPU).Nspace;
    (*sCPU).kStep = TWOPI / ((*sCPU).Nspace * (*sCPU).rStep);
    (*sCPU).fStep = 1.0 / ((*sCPU).Ntime * (*sCPU).tStep);
    (*sCPU).Npropagation = (size_t)round((*sCPU).crystalThickness / (*sCPU).propagationStep);

    (*sCPU).isCylindric = (*sCPU).symmetryType == 1;
    if ((*sCPU).isCylindric) {
        (*sCPU).x01 = 0;
        (*sCPU).x02 = 0;
        (*sCPU).propagationAngle1 = 0;
        (*sCPU).propagationAngle2 = 0;
    }

    if ((*sCPU).batchIndex == 0 || (*sCPU).Nsims < 1) {
        (*sCPU).Nsims = 1;
    }
    if ((*sCPU).batchIndex2 == 0 || (*sCPU).Nsims2 < 1) {
        (*sCPU).Nsims2 = 1;
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
    return readValueCount;
}

//print a linefeed without a carriage return so that linux systems don't complain
//about impure scripts from DOS machines
//fopen() should be called with "wb"
void unixNewLine(FILE* iostream) {
    char LF = '\x0A';
    fwrite(&LF, sizeof(char), 1, iostream);
}


int saveSlurmScript(simulationParameterSet* sCPU, int gpuType, int gpuCount) {
    FILE* textfile;
    char* stringConversionBuffer = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    wchar_t* wideStringConversionBuffer = (wchar_t*)calloc(MAX_LOADSTRING, sizeof(char));
    char* outputpath = (char*)calloc(MAX_LOADSTRING, sizeof(char));

    char* fileName = (*sCPU).outputBasePath;
    while (strchr(fileName, '\\') != NULL) {
        fileName = strchr(fileName, '\\');
        fileName++;
    }
    strcpy(outputpath, (*sCPU).outputBasePath);
    strcat(outputpath, ".slurmScript");
    textfile = fopen(outputpath, "wb");
    fprintf(textfile, "#!/bin/bash -l"); unixNewLine(textfile);
    fprintf(textfile, "#SBATCH -o ./tjob.out.%%j"); unixNewLine(textfile);
    fprintf(textfile, "#SBATCH -e ./tjob.err.%%j"); unixNewLine(textfile);
    fprintf(textfile, "#SBATCH -D ./"); unixNewLine(textfile);
    fprintf(textfile, "#SBATCH -J lightwave");  unixNewLine(textfile);
    fprintf(textfile, "#SBATCH --constraint=\"gpu\""); unixNewLine(textfile);
    if (gpuType == 0) {
        fprintf(textfile, "#SBATCH --gres=gpu:rtx5000:%i", min(gpuCount,2)); unixNewLine(textfile);
    }
    if (gpuType == 1) {
        fprintf(textfile, "#SBATCH --gres=gpu:v100:%i", min(gpuCount, 2)); unixNewLine(textfile);
    }
    if (gpuType == 2) {
        fprintf(textfile, "#SBATCH --gres=gpu:a100:%i", min(gpuCount, 4)); unixNewLine(textfile);
        fprintf(textfile, "#SBATCH --cpus-per-task=%i", 2*min(gpuCount, 4)); unixNewLine(textfile);
    }
    fprintf(textfile, "#SBATCH --mem=%lliM",1024+(18 * sizeof(double) * (*sCPU).Ngrid * max(1,(*sCPU).Nsims))/1048576);
    unixNewLine(textfile);
    fprintf(textfile, "#SBATCH --nodes=1"); unixNewLine(textfile);
    fprintf(textfile, "#SBATCH --ntasks-per-node=1"); unixNewLine(textfile);
    fprintf(textfile, "#SBATCH --time=24:00:00"); unixNewLine(textfile);
    fprintf(textfile, "module purge"); unixNewLine(textfile);
    fprintf(textfile, "module load cuda/11.2"); unixNewLine(textfile);
    fprintf(textfile, "module load mkl/2022.0"); unixNewLine(textfile);
    fprintf(textfile, "module load gcc/9"); unixNewLine(textfile);
    fprintf(textfile, "export LD_LIBRARY_PATH=$MKL_HOME/lib/intel64:$LD_LIBRARY_PATH"); unixNewLine(textfile);
    if (gpuType == 0 || gpuType == 1) {
        fprintf(textfile, "srun ./nnp %s.input > %s.out", fileName, fileName); unixNewLine(textfile);
    }
    if (gpuType == 2) {
        fprintf(textfile, "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}"); unixNewLine(textfile);
        fprintf(textfile, "srun ./nnp %s.input > %s.out", fileName, fileName); unixNewLine(textfile);
    }
    fclose(textfile);
    free(outputpath);
    free(wideStringConversionBuffer);
    free(stringConversionBuffer);
    return 0;
}

int saveSettingsFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr) {
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
    fwprintf(textfile, L"Pulse energy 1 (J): %14.14e\nPulse energy 2 (J): %14.14e\nFrequency 1 (Hz): %14.14e\n", 
        (*sCPU).pulseEnergy1, (*sCPU).pulseEnergy2, (*sCPU).frequency1);
    fwprintf(textfile, L"Frequency 2 (Hz): %14.14e\nBandwidth 1 (Hz): %14.14e\nBandwidth 2 (Hz): %14.14e\n", 
        (*sCPU).frequency2, (*sCPU).bandwidth1, (*sCPU).bandwidth2);
    fwprintf(textfile, L"SG order 1: %i\nSG order 2: %i\nCEP 1 (rad): %14.14e\nCEP 2 (rad): %14.14e\nDelay 1 (s): %14.14e\nDelay 2 (s): %14.14e\nGDD 1 (s^-2): %14.14e\nGDD 2 (s^-2): %14.14e\nTOD 1 (s^-3): %14.14e\nTOD 2 (s^-3): %14.14e\n", 
        (*sCPU).sgOrder1, (*sCPU).sgOrder2, (*sCPU).cephase1, (*sCPU).cephase2, (*sCPU).delay1, (*sCPU).delay2, (*sCPU).gdd1, (*sCPU).gdd2, (*sCPU).tod1, (*sCPU).tod2);
    fwprintf(textfile, L"Phase material 1 index: %i\nPhase material 2 index: %i\nPhase material thickness 1 (mcr.): %14.14e\nPhase material thickness 2 (mcr.): %14.14e\n", 
        (*sCPU).phaseMaterialIndex1, (*sCPU).phaseMaterialIndex2, (*sCPU).phaseMaterialThickness1, (*sCPU).phaseMaterialThickness2);
    fwprintf(textfile, L"Beamwaist 1 (m): %14.14e\nBeamwaist 2 (m): %14.14e\nx offset 1 (m): %14.14e\nx offset 2 (m): %14.14e\nz offset 1 (m): %14.14e\nz offset 2 (m): %14.14e\nNC angle 1 (rad): %14.14e\nNC angle 2 (rad): %14.14e\n", 
        (*sCPU).beamwaist1, (*sCPU).beamwaist2, (*sCPU).x01, (*sCPU).x02, (*sCPU).z01, (*sCPU).z02, (*sCPU).propagationAngle1, (*sCPU).propagationAngle2);
    fwprintf(textfile, L"Polarization 1 (rad): %14.14e\nPolarization 2 (rad): %14.14e\nCircularity 1: %14.14e\nCircularity 2: %14.14e\n", 
        (*sCPU).polarizationAngle1, (*sCPU).polarizationAngle2, (*sCPU).circularity1, (*sCPU).circularity2);
    fwprintf(textfile, L"Material index: %i\nAlternate material index: %i\n", 
        (*sCPU).materialIndex, (*sCPU).materialIndexAlternate);
    fwprintf(textfile, L"Crystal theta (rad): %14.14e\nCrystal phi (rad): %14.14e\nGrid width (m): %14.14e\nGrid height (m): %14.14e\ndx (m): %14.14e\nTime span (s): %14.14e\ndt (s): %14.14e\nThickness (m): %14.14e\ndz (m): %14.14e\n", 
        (*sCPU).crystalTheta, (*sCPU).crystalPhi, (*sCPU).spatialWidth, (*sCPU).spatialHeight, (*sCPU).rStep, (*sCPU).timeSpan, (*sCPU).tStep, (*sCPU).crystalThickness, (*sCPU).propagationStep);
    fwprintf(textfile, L"Nonlinear absorption parameter: %14.14e\nBand gap (eV): %14.14e\nEffective mass (relative): %14.14e\nDrude gamma (Hz): %14.14e\n", 
        (*sCPU).nonlinearAbsorptionStrength, (*sCPU).bandGapElectronVolts, (*sCPU).effectiveMass, (*sCPU).drudeGamma);
    fwprintf(textfile, L"Propagation mode: %i\n", 
        (*sCPU).symmetryType);
    fwprintf(textfile, L"Batch mode: %i\nBatch destination: %14.14e\nBatch steps: %lli\n", 
        (*sCPU).batchIndex, (*sCPU).batchDestination, (*sCPU).Nsims);
    fwprintf(textfile, L"Batch mode 2: %i\nBatch destination 2: %14.14e\nBatch steps 2: %lli\n",
        (*sCPU).batchIndex2, (*sCPU).batchDestination2, (*sCPU).Nsims2);
    mbstowcs(wideStringConversionBuffer, (*sCPU).sequenceString, MAX_LOADSTRING);
    fwprintf(textfile, L"Sequence: %ls\n", wideStringConversionBuffer);
    mbstowcs(wideStringConversionBuffer, (*sCPU).fittingString, MAX_LOADSTRING);
    fwprintf(textfile, L"Fitting: %ls\n", wideStringConversionBuffer);
    fwprintf(textfile, L"Fitting mode: %i\n", (*sCPU).fittingMode);

    if ((*sCPU).runType > 0) {
        char* fileName = (*sCPU).outputBasePath;
        while (strchr(fileName, '\\') != NULL) {
            fileName = strchr(fileName, '\\');
            fileName++;
        }
        mbstowcs(wideStringConversionBuffer, fileName, strlen(fileName));
        wideStringConversionBuffer[strlen(fileName)] = L'\0';
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
    mbstowcs(wideStringConversionBuffer, (*sCPU).fittingPath, MAX_LOADSTRING);
    fwprintf(textfile, L"Fitting reference file path: %ls\n", wideStringConversionBuffer);

    fwprintf(textfile, L"Material name: %ls\nSellmeier reference: %ls\nChi2 reference: %ls\nChi3 reference: %ls\n", crystalDatabasePtr[(*sCPU).materialIndex].crystalNameW, crystalDatabasePtr[(*sCPU).materialIndex].sellmeierReference, crystalDatabasePtr[(*sCPU).materialIndex].dReference, crystalDatabasePtr[(*sCPU).materialIndex].chi3Reference);
    fwprintf(textfile, L"Sellmeier coefficients: \n");
    for (j = 0; j < 3; j++) {
        for (k = 0; k < 22; k++) {
            fwprintf(textfile, L"%14.14e ", crystalDatabasePtr[(*sCPU).materialIndex].sellmeierCoefficients[j * 22 + k]);
        }
        fwprintf(textfile, L"\n");
    }
    fwprintf(textfile, L"Code version: 0.31 August 15, 2022\n");

    fclose(textfile);
    free(outputpath);
    free(wideStringConversionBuffer);
    free(stringConversionBuffer);
    return 0;
}

int removeCharacterFromString(char* cString, size_t N, char removedChar) {
    size_t i = 0;
    size_t r = 0;
    while (i < N-1) {
        if (cString[i] == removedChar) {
            memmove(&cString[i], &cString[i + 1], N - i - r - 1);
            cString[N - r - 1] = 0;
            r++;
        }
        else {
            i++;
        }
    }
    if (cString[N - 1] == removedChar) {
        cString[N - 1] = '\0';
    }
    return 0;
}

int saveDataSet(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, char* outputbase, bool saveInputs) {
    int j;

    saveSettingsFile(sCPU, crystalDatabasePtr);

    //Save the results as double instead of complex
    double* saveEout = (double*)calloc((*sCPU).Ngrid * 2 * (*sCPU).Nsims * (*sCPU).Nsims2, sizeof(double));
    for (j = 0; j < ((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2 * 2); j++) {
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
    for (j = 0; j < ((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2 * 2); j++) {
        saveEout[j] = real((*sCPU).ExtOut[j]);
    }
    FILE* ExtOutFile;
    size_t writeSize = 2 * ((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2);
    strcpy(outputpath, outputbase);
    strcat(outputpath, "_Ext.dat");
    ExtOutFile = fopen(outputpath, "wb");
    fwrite(saveEout, sizeof(double), writeSize, ExtOutFile);
    fwrite(matlabpadding, sizeof(double), 1024, ExtOutFile);
    fclose(ExtOutFile);
    if (saveInputs) {
        for (j = 0; j < ((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2 * 2); j++) {
            saveEout[j] = real((*sCPU).Ext[j]);
        }
        FILE* ExtInFile;
        strcpy(outputpath, outputbase);
        strcat(outputpath, "_ExtIn.dat");
        ExtInFile = fopen(outputpath, "wb");
        fwrite(saveEout, sizeof(double), writeSize, ExtInFile);
        fwrite(matlabpadding, sizeof(double), 1024, ExtInFile);
        fclose(ExtInFile);
    }

    //Save the spectrum
    FILE* totalSpectrumFile;
    strcpy(outputpath, outputbase);
    strcat(outputpath, "_spectrum.dat");
    totalSpectrumFile = fopen(outputpath, "wb");
    fwrite((*sCPU).totalSpectrum, sizeof(double), 3 * (*sCPU).Ntime * (*sCPU).Nsims * (*sCPU).Nsims2, totalSpectrumFile);
    fwrite(matlabpadding, sizeof(double), 1024, totalSpectrumFile);
    fclose(totalSpectrumFile);

    FILE* matlabfile;
    strcpy(outputpath, outputbase);
    strcat(outputpath, ".m");
    matlabfile = fopen(outputpath, "w");
    
    if (saveInputs) {
        fprintf(matlabfile, "fid = fopen('%s_ExtIn.dat','rb'); \n", outputbaseVar);
        fprintf(matlabfile, "%s_ExtIn = fread(fid, %lli, 'double'); \n", outputbaseVar, 2 * (*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2);
        fprintf(matlabfile, "%s_ExtIn = reshape(%s_ExtIn,[%lli %lli %lli]); \n", outputbaseVar, outputbaseVar, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims * (*sCPU).Nsims2);
        fprintf(matlabfile, "fclose(fid); \n");
    }
    
    fprintf(matlabfile, "fid = fopen('%s_Ext.dat','rb'); \n", outputbaseVar);
    fprintf(matlabfile, "%s_Ext = fread(fid, %lli, 'double'); \n", outputbaseVar, 2 * (*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2);
    fprintf(matlabfile, "%s_Ext = reshape(%s_Ext,[%lli %lli %lli]); \n", outputbaseVar, outputbaseVar, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims * (*sCPU).Nsims2);
    fprintf(matlabfile, "fclose(fid); \n");
    fprintf(matlabfile, "fid = fopen('%s_spectrum.dat','rb'); \n", outputbaseVar);
    fprintf(matlabfile, "%s_spectrum = fread(fid, %lli, 'double'); \n", outputbaseVar, 3 * (*sCPU).Ntime * (*sCPU).Nsims * (*sCPU).Nsims2);
    fprintf(matlabfile, "%s_spectrum = reshape(%s_spectrum,[%lli %i %zi]); \n", outputbaseVar, outputbaseVar, (*sCPU).Ntime, 3, (*sCPU).Nsims * (*sCPU).Nsims2);
    fprintf(matlabfile, "fclose(fid); \n");
    fprintf(matlabfile, "dt = %e;\ndz = %e;\ndx = %e;\ndf = %e;\n", (*sCPU).tStep, (*sCPU).propagationStep, (*sCPU).rStep, (*sCPU).fStep);
    fclose(matlabfile);
    
    //write a python script for loading the output fields in a proper shape
    char scriptfilename[MAX_LOADSTRING];
    strcpy(scriptfilename, outputbase);
    strcat(scriptfilename, ".py");
    FILE* scriptfile;
    scriptfile = fopen(scriptfilename, "w");
    fprintf(scriptfile, "#!/usr/bin/python\nimport numpy as np\n");
    fprintf(scriptfile, "dt = %e\ndz = %e\ndx = %e\ndf = %e\n", (*sCPU).tStep, (*sCPU).propagationStep, (*sCPU).rStep, (*sCPU).fStep);
    if (saveInputs) {
        fprintf(scriptfile, "%s_ExtIn = np.reshape(np.fromfile(\"", outputbaseVar);
        fprintf(scriptfile, "%s_ExtIn.dat", outputbaseVar);
        fprintf(scriptfile, "\",dtype=np.double)[0:%lli],(%lli,%lli,%lli),order='F')\n", 2 * (*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims * (*sCPU).Nsims2);
    }
    fprintf(scriptfile, "%s_Ext = np.reshape(np.fromfile(\"", outputbaseVar);
    fprintf(scriptfile, "%s_Ext.dat", outputbaseVar);
    fprintf(scriptfile, "\",dtype=np.double)[0:%lli],(%lli,%lli,%lli),order='F')\n", 2 * (*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims * (*sCPU).Nsims2);
    fprintf(scriptfile, "%s_spectrum = np.reshape(np.fromfile(\"", outputbaseVar);
    fprintf(scriptfile, "%s_spectrum.dat", outputbaseVar);
    fprintf(scriptfile, "\",dtype=np.double)[0:%lli],(%lli,%i,%zi),order='F')\n", 3 * (*sCPU).Ntime * (*sCPU).Nsims * (*sCPU).Nsims2, (*sCPU).Ntime, 3, (*sCPU).Nsims * (*sCPU).Nsims2);
    fclose(scriptfile);
    
    free(saveEout);
    free(matlabpadding);
    free(stringConversionBuffer);
    free(wideStringConversionBuffer);
    return 0;
}

int resolveSequence(int currentIndex, simulationParameterSet* s, crystalEntry* db) {
    double* offsetArray = &(*s).sequenceArray[11 * currentIndex];
    
    //sequence format
    //0: step type
    int stepType = (int)offsetArray[0];
    // 
    // if stepType == 0, normal propagation
    //1: material index
    //2: theta,
    //3: phi, 
    //4: NL absorption
    //5: Band gap
    //6: Drude relaxation
    //7: Effective mass
    //8: Crystal thickness
    //9: Propagation step size
    //10: rotation angle
    //
    // if stepType == 1, linear propagation
    // same parameters as 0, but only 1,2,3,8, and 10 matter
    //
    // if stepType == 2, fresnel loss
    // 1: incidence material index
    // 2: transmission material index
    // other parameters don't matter
    // 
    // if stepType == 3, spherical mirror
    // 1: ROC (m)

    
    switch (stepType) {
    case 0:
        if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
        if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
        if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
        if ((int)offsetArray[4] != -1) (*s).nonlinearAbsorptionStrength = offsetArray[4];
        if ((int)offsetArray[5] != -1) (*s).bandGapElectronVolts = offsetArray[5];
        if ((int)offsetArray[6] != -1) (*s).drudeGamma = offsetArray[6];
        if ((int)offsetArray[7] != -1) (*s).effectiveMass = offsetArray[7];
        if ((int)offsetArray[8] != -1) (*s).crystalThickness = 1e-6 * offsetArray[8];
        if ((int)offsetArray[9] != -1) (*s).propagationStep = 1e-9 * offsetArray[9];
        if ((int)offsetArray[8] != -1 && (int)offsetArray[8] != -1) (*s).Npropagation 
            = (size_t)(1e-6 * offsetArray[8] / (*s).propagationStep);
        if (currentIndex > 0) {
            (*s).isFollowerInSequence = TRUE;
        }
        (*s).chi2Tensor = db[(*s).materialIndex].d;
        (*s).chi3Tensor = db[(*s).materialIndex].chi3;
        (*s).nonlinearSwitches = db[(*s).materialIndex].nonlinearSwitches;
        (*s).absorptionParameters = db[(*s).materialIndex].absorptionParameters;
        (*s).sellmeierCoefficients = db[(*s).materialIndex].sellmeierCoefficients;

        (*s).sellmeierType = db[(*s).materialIndex].sellmeierType;
        (*s).axesNumber = db[(*s).materialIndex].axisType;


        solveNonlinearWaveEquation(s);
        if (offsetArray[10] != 0.0) {
            rotateField(s, DEG2RAD * offsetArray[10]);
        }

        if ((*s).memoryError > 0) {
            printf("Warning: device memory error (%i).\n", (*s).memoryError);
        }
        return 0;
    case 1:
        if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
        if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
        if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
        applyLinearPropagation(s, (int)offsetArray[1], 1.0e6*offsetArray[8]);
        if (offsetArray[10] != 0.0) {
            rotateField(s, DEG2RAD * offsetArray[10]);
        }

        return 0;
    case 2:
        if ((int)offsetArray[1] != -1) (*s).materialIndex = (int)offsetArray[1];
        if ((int)offsetArray[2] != -1) (*s).crystalTheta = DEG2RAD * offsetArray[2];
        if ((int)offsetArray[3] != -1) (*s).crystalPhi = DEG2RAD * offsetArray[3];
        applyFresnelLoss(s, (int)offsetArray[4], (int)offsetArray[5]);
        return 0;
    }
    


    return 1;
}

int loadPulseFiles(simulationParameterSet* sCPU) {

    //pulse type specifies if something has to be loaded to describe the pulses, or if they should be
    //synthesized later. 1: FROG .speck format; 2: EOS (not implemented yet)
    int frogLines = 0;
    int errCount = 0;
    if ((*sCPU).pulse1FileType == 1) {
        frogLines = loadFrogSpeck((*sCPU).field1FilePath, (*sCPU).loadedField1, (*sCPU).Ntime, (*sCPU).fStep, 0.0, 1);
        if (frogLines > 1) {
            (*sCPU).field1IsAllocated = TRUE;
        }
        else {
            (*sCPU).field1IsAllocated = FALSE;
            errCount++;
        }
    }

    if ((*sCPU).pulse2FileType == 1) {
        frogLines = loadFrogSpeck((*sCPU).field2FilePath, (*sCPU).loadedField2, (*sCPU).Ntime, (*sCPU).fStep, 0.0, 1);
        if (frogLines > 1) {
            (*sCPU).field2IsAllocated = TRUE;
        }
        else {
            (*sCPU).field2IsAllocated = FALSE;
            errCount++;
        }
    }
    return errCount;
}

int loadSavedFields(simulationParameterSet* sCPU, char* outputBase, bool GPUisPresent) {
    char* outputpath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    size_t writeSize = 2 * ((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2);
    double* loadE = (double*)malloc(writeSize * sizeof(double));
    size_t j;

    //read fields as binary
    FILE* ExtOutFile;
    
    strcpy(outputpath, outputBase);
    strcat(outputpath, "_Ext.dat");
    ExtOutFile = fopen(outputpath, "rb");
    if (ExtOutFile == NULL) {
        return 1;
    }
    fread(loadE, sizeof(double), writeSize, ExtOutFile);
    fclose(ExtOutFile);
    for (j = 0; j < writeSize; j++) {
        (*sCPU).ExtOut[j] = loadE[j];
    }
    free(loadE);

    FILE* spectrumFile;
    strcpy(outputpath, outputBase);
    strcat(outputpath, "_spectrum.dat");
    spectrumFile = fopen(outputpath, "rb");
    fread((*sCPU).totalSpectrum, sizeof(double), (*sCPU).Nsims * (*sCPU).Nsims2 * 3 * (*sCPU).Ntime, spectrumFile);
    fclose(spectrumFile);
    free(outputpath);
    if (GPUisPresent) {
        cufftHandle fftPlan;
        cufftPlan2d(&fftPlan, (int)(*sCPU).Nspace, (int)(*sCPU).Ntime, CUFFT_Z2Z);

        cuDoubleComplex* fieldGridkw;
        cuDoubleComplex* fieldGridxt;
        cudaMalloc((void**)&fieldGridkw, sizeof(cuDoubleComplex) * (*sCPU).Ngrid);
        cudaMalloc((void**)&fieldGridxt, sizeof(cuDoubleComplex) * (*sCPU).Ngrid);

        for (j = 0; j < 2 * (*sCPU).Nsims * (*sCPU).Nsims2; j++) {
            cudaMemcpy(fieldGridxt, &(*sCPU).ExtOut[j * (*sCPU).Ngrid], (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
            cufftExecZ2Z(fftPlan, fieldGridxt, fieldGridkw, CUFFT_FORWARD);
            cudaMemcpy(&(*sCPU).EkwOut[j * (*sCPU).Ngrid], fieldGridkw, (*sCPU).Ngrid * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
        }
        cudaFree(fieldGridkw);
        cudaFree(fieldGridxt);
        cufftDestroy(fftPlan);
    }
    else {
        DFTI_DESCRIPTOR_HANDLE dftiHandle = NULL;
        MKL_LONG fftDimensions[2] = { (long)(*sCPU).Nspace , (long)(*sCPU).Ntime };
        MKL_LONG mklError = 0;
        mklError = DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 2, fftDimensions);
        DftiSetValue(dftiHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
        if (mklError != DFTI_NO_ERROR) return 1;
        mklError = DftiCommitDescriptor(dftiHandle);
        if (mklError != DFTI_NO_ERROR) return 2;
        for (j = 0; j < (2 * (*sCPU).Nsims * (*sCPU).Nsims2); j++) {
            mklError = DftiComputeForward(dftiHandle, &(*sCPU).ExtOut[j*(*sCPU).Ngrid], &(*sCPU).EkwOut[j*(*sCPU).Ngrid]);
            if (mklError != DFTI_NO_ERROR) return 3;
        }

        DftiFreeDescriptor(&dftiHandle);
    }


    return 0;
}


int getTotalSpectrum(simulationParameterSet* sCPU, cudaParameterSet* sc) {
    cufftHandle plan1;
    cufftPlan1d(&plan1, (int)(*sCPU).Ntime, CUFFT_Z2Z, (int)(*sCPU).Nspace);
    cufftSetStream(plan1, (*sc).CUDAStream);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).workspace1, (cufftDoubleComplex*)(*sc).workspace1, CUFFT_FORWARD);
    cufftExecZ2Z(plan1, (cufftDoubleComplex*)(*sc).workspace2, (cufftDoubleComplex*)(*sc).workspace2, CUFFT_FORWARD);

    totalSpectrumKernel<<<(unsigned int)(*sCPU).Ntime, 1, 0, (*sc).CUDAStream>>>((*sc).workspace1, (*sc).workspace2, (*sCPU).rStep, (*sCPU).Ntime, (*sCPU).Nspace, (*sc).gridPolarizationTime1);
    cudaDeviceSynchronize();
    cudaMemcpy((*sCPU).totalSpectrum, (*sc).gridPolarizationTime1, 3 * (*sCPU).Ntime * sizeof(double), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cufftDestroy(plan1);
    return 0;
}

unsigned long runFitting(simulationParameterSet* sCPU) {
    int n = (int)(*sCPU).Nfitting;
    int m = (int)(*sCPU).fittingROIsize;
    fittingReferenceSet = sCPU;
    fittingSet = (simulationParameterSet*)malloc((*sCPU).Nsims * sizeof(simulationParameterSet));
    memcpy(fittingSet, sCPU, (*sCPU).Nsims * sizeof(simulationParameterSet));
    
    double commonPrecision = (*sCPU).fittingPrecision;
    const double eps[6] = { commonPrecision,commonPrecision,commonPrecision,commonPrecision,commonPrecision,commonPrecision }; /* set precisions for stop-criteria */
    double jacobianPrecision = commonPrecision;
    double* x = (double*)mkl_malloc(sizeof(double) * n, 64);
    double* fittingValues = (double*)mkl_malloc(sizeof(double) * m, 64);
    double* fjac = (double*)mkl_malloc(sizeof(double) * m * n, 64);
    double* lowerBounds = (double*)mkl_malloc(sizeof(double) * n, 64);
    double* upperBounds = (double*)mkl_malloc(sizeof(double) * n, 64);
    const int maxIterations = max((*sCPU).fittingMaxIterations,2);
    const int maxTrialIterations = max(maxIterations/10, 2);
    /* initial step bound */
    double rs = 0.0;
    int RCI_Request;
    int successful;

    int iter;
    int stopCriterion;
    double inputResiduals = 0.0, outputResiduals = 0.0;
    _TRNSPBC_HANDLE_t handle;
    int i;
    int error = 0;

    //initial guess and bounds
    for (i = 0; i < n; i++) {
        x[i] = 1.;
        upperBounds[i] = (*fittingSet).fittingArray[3 * i + 2];
        lowerBounds[i] = (*fittingSet).fittingArray[3 * i + 1];
    }

    //initialize fitting function and jacobian
    for (i = 0; i < m; i++) {
        fittingValues[i] = 0.0;
    }
    for (i = 0; i < m * n; i++) {
        fjac[i] = 0.0;
    }

    error += dtrnlspbc_init(&handle, &n, &m, x, lowerBounds, upperBounds, eps, &maxIterations, &maxTrialIterations, &rs) != TR_SUCCESS;
    size_t currentIteration = 0;
    if (error == 0) {
        RCI_Request = 0;
        successful = 0;
        while (successful == 0 && (*sCPU).imdone[0] != 2 && currentIteration < maxIterations)
        {
            currentIteration++;
            if (dtrnlspbc_solve(&handle, fittingValues, fjac, &RCI_Request) != TR_SUCCESS)
            {
                successful = -1;
            }

            //check convergence
            if (RCI_Request > -7 && RCI_Request < -1) successful = 1;

            //recalculate
            if (RCI_Request == 1)
            {
                runFittingIteration(&m, &n, x, fittingValues);
            }

            //make jacobian
            if (RCI_Request == 2)
            {
                djacobi(runFittingIteration, &n, &m, fjac, x, &jacobianPrecision);
            }
        }
    }

    /* get solution statuses
       handle            in:        TR solver handle
       iter              out:       number of iterations
       stopCriterion             out:       number of stop criterion
       inputResiduals                out:       initial residuals
       outputResiduals                out:       final residuals */
    dtrnlspbc_get(&handle, &iter, &stopCriterion, &inputResiduals, &outputResiduals);
    memcpy(sCPU, fittingSet, (*fittingSet).Nsims * sizeof(simulationParameterSet));
    //free memory
    dtrnlspbc_delete(&handle);
    mkl_free(upperBounds);
    mkl_free(lowerBounds);
    mkl_free(fjac);
    mkl_free(fittingValues);
    mkl_free(x);
    MKL_Free_Buffers();
    free(fittingSet);
    return 0;
}


void runFittingIteration(int* m, int* n, double* fittingValues, double* fittingFunction) {
    int i;
    int fitLocation;
    double referenceValue;
    //pointers to values that can be scanned in batch mode
    double* targets[36] = { 0,
        &(*fittingSet).pulseEnergy1, &(*fittingSet).pulseEnergy2, &(*fittingSet).frequency1, &(*fittingSet).frequency2,
        &(*fittingSet).bandwidth1, &(*fittingSet).bandwidth2, &(*fittingSet).cephase1, &(*fittingSet).cephase2,
        &(*fittingSet).delay1, &(*fittingSet).delay2, &(*fittingSet).gdd1, &(*fittingSet).gdd2,
        &(*fittingSet).tod1, &(*fittingSet).tod2, &(*fittingSet).phaseMaterialThickness1, &(*fittingSet).phaseMaterialThickness2,
        &(*fittingSet).beamwaist1, &(*fittingSet).beamwaist2,
        &(*fittingSet).x01, &(*fittingSet).x02, &(*fittingSet).z01, &(*fittingSet).z02,
        &(*fittingSet).propagationAngle1, &(*fittingSet).propagationAngle2, &(*fittingSet).polarizationAngle1, &(*fittingSet).polarizationAngle2,
        &(*fittingSet).circularity1, &(*fittingSet).circularity2, &(*fittingSet).crystalTheta, &(*fittingSet).crystalPhi,
        &(*fittingSet).nonlinearAbsorptionStrength, &(*fittingSet).drudeGamma, &(*fittingSet).effectiveMass, &(*fittingSet).crystalThickness,
        &(*fittingSet).propagationStep };

    double* references[36] = { 0,
    &(*fittingReferenceSet).pulseEnergy1, &(*fittingReferenceSet).pulseEnergy2, &(*fittingReferenceSet).frequency1, &(*fittingReferenceSet).frequency2,
    &(*fittingReferenceSet).bandwidth1, &(*fittingReferenceSet).bandwidth2, &(*fittingReferenceSet).cephase1, &(*fittingReferenceSet).cephase2,
    &(*fittingReferenceSet).delay1, &(*fittingReferenceSet).delay2, &(*fittingReferenceSet).gdd1, &(*fittingReferenceSet).gdd2,
    &(*fittingReferenceSet).tod1, &(*fittingReferenceSet).tod2, &(*fittingReferenceSet).phaseMaterialThickness1, &(*fittingReferenceSet).phaseMaterialThickness2,
    &(*fittingReferenceSet).beamwaist1, &(*fittingReferenceSet).beamwaist2,
    &(*fittingReferenceSet).x01, &(*fittingReferenceSet).x02, &(*fittingReferenceSet).z01, &(*fittingReferenceSet).z02,
    &(*fittingReferenceSet).propagationAngle1, &(*fittingReferenceSet).propagationAngle2, &(*fittingReferenceSet).polarizationAngle1, &(*fittingReferenceSet).polarizationAngle2,
    &(*fittingReferenceSet).circularity1, &(*fittingReferenceSet).circularity2, &(*fittingReferenceSet).crystalTheta, &(*fittingReferenceSet).crystalPhi,
    &(*fittingReferenceSet).nonlinearAbsorptionStrength, &(*fittingReferenceSet).drudeGamma, &(*fittingReferenceSet).effectiveMass, &(*fittingReferenceSet).crystalThickness,
    &(*fittingReferenceSet).propagationStep };

    
    for (i = 0; i < *n; i++) {
        fitLocation = (int)round((*fittingSet).fittingArray[3 * i]);
        referenceValue = *references[fitLocation];
        if (referenceValue == 0.0) {
            referenceValue = 1.;
        }
        *targets[fitLocation] = fittingValues[i] * referenceValue;
    }
    if ((*fittingSet).isInSequence) {
        solveNonlinearWaveEquationSequence(fittingSet);
        (*fittingSet).isFollowerInSequence = FALSE;
    }
    else {
        solveNonlinearWaveEquation(fittingSet);
    }
    

    //mode 0: maximize total spectrum in ROI
    if ((*fittingSet).fittingMode == 0) {
        for (i = 0; i < *m; i++) {
            fittingFunction[i] = (1.0e8 / ((*fittingSet).totalSpectrum[2 * (*fittingSet).Ntime + (*fittingSet).fittingROIstart + i]));
        }
    }
    //mode 1: maximize s-polarized spectrum in ROI
    if ((*fittingSet).fittingMode == 1) {
        for (i = 0; i < *m; i++) {
            fittingFunction[i] = (1.0e8 / ((*fittingSet).totalSpectrum[(*fittingSet).fittingROIstart + i]));
        }
    }
    //mode 2: maximize p-polarized spectrum in ROI
    if ((*fittingSet).fittingMode == 2) {
        for (i = 0; i < *m; i++) {
            fittingFunction[i] = (1.0e8 / ((*fittingSet).totalSpectrum[(*fittingSet).Ntime + (*fittingSet).fittingROIstart + i]));
        }
    }
    //mode 3: match total spectrum to reference given in ascii file
    if ((*fittingSet).fittingMode == 3) {
        double maxSim = 0;
        double maxRef = 0;
        double* simSpec = &(*fittingSet).totalSpectrum[2 * (*fittingSet).Ntime + (*fittingSet).fittingROIstart];
        double* refSpec = &(*fittingSet).fittingArray[(*fittingSet).fittingROIstart];
        for (i = 0; i < *m; i++) {
            maxSim = max(maxSim, simSpec[i]);
            maxRef = max(maxRef, refSpec[i]);
        }
        if (maxSim == 0) {
            maxSim = 1;
        }
        if (maxRef == 0) {
            maxRef = 1;
        }
        for (i = 0; i < *m; i++) {
            fittingFunction[i] = log10(1e5*refSpec[i] / maxRef) - log10(1e5*simSpec[i] / maxSim);
            //fittingFunction[i] = 1.0e8 / ((*fittingSet).totalSpectrum[(*fittingSet).Ntime + (*fittingSet).fittingROIstart + i]);
        }
    }


    return;
}

//generate the rotation matricies for translating between the beam coordinates and
//crystal coordinates
int fillRotationMatricies(simulationParameterSet* sCPU, cudaParameterSet* s) {
    double cosT = cos((*sCPU).crystalTheta);
    double sinT = sin((*sCPU).crystalTheta);
    double cosP = cos((*sCPU).crystalPhi);
    double sinP = sin((*sCPU).crystalPhi);
    double forward[9] = 
        { cosP, sinP, 0, -cosT * sinP, cosT * cosP, sinT, sinT * sinP, -sinT * cosP, cosT };

    //reverse direction (same array contents)
    sinT *= -1;
    sinP *= -1;
    double backward[9] = 
        { cosP, sinP, 0, -cosT * sinP, cosT * cosP, sinT, sinT * sinP, -sinT * cosP, cosT };

    memcpy((*s).rotationForward, forward, 9 * sizeof(double));
    memcpy((*s).rotationBackward, backward, 9 * sizeof(double));
    return 0;
}