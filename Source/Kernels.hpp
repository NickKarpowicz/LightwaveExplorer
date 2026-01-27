#include "DeviceFunctions.hpp"
#include "LightwaveExplorerTrilingual.h"
#include "MaxwellDeviceFunctions.hpp"

using namespace deviceFunctions;
namespace kernelNamespace {
// Kernels are written in the form of functors, which can be passed to the three
// main
//  platforms (CUDA, SYCL, c++) without using macros or much in the way of
//  wrappers. All members are public so that I can do aggregate initliazation
//  without writing constructors for all of them. (Could therefore be slightly
//  simplified by using structs instead of classes, but I think a class makes
//  more sense in this context). A macro is used to name the namespace between
//  different modes to avoid ODR violations.

// calculate the total energy spectrum of the beam for the 2D modes. Note that
// for the cartesian one, it will be treated as a round beam instead of an
// infinite plane wave in the transverse direction. Thus, the 2D Cartesian
// spectra are approximations.
class totalSpectrumKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t i) const {
    deviceFP beamCenter1{};
    deviceFP beamCenter2{};
    deviceFP beamTotal1{};
    deviceFP beamTotal2{};

    // find beam centers
    if ((*s).isCylindric) {
      beamCenter1 = ((*s).Nspace / 2.0f * (*s).dx) + 0.25f * (*s).dx;
      beamCenter2 = beamCenter1;
    } else {
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

    // Integrate total beam power, assuming radially-symmetric beam around
    // the center
    beamTotal1 = 0.0f;
    beamTotal2 = 0.0f;
    for (int64_t j{}; j < (*s).Nspace; ++j) {
      deviceFP x = (*s).dx * j;
      beamTotal1 += deviceFPLib::abs(x - beamCenter1) *
                    modulusSquared((*s).workspace1[i + j * (*s).Nfreq]);
      beamTotal2 += deviceFPLib::abs(x - beamCenter2) *
                    modulusSquared((*s).workspace2[i + j * (*s).Nfreq]);
    }

    // put the values into the output spectrum
    (*s).gridPolarizationTime1[i] = beamTotal1;
    (*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
    (*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
  }
};

// Calculate the energy spectrum after a 2D propagation assuming that the beam
// height in the non-resolved direction is == the grid width (i.e. square grid)
// More quantitative than the mapping to round beams, but rather specific
//  DISABLED IN FAVOR OF ROUND-BEAM APPROXIMATION
//	int64_t j;
//	deviceFP beamTotal1 = 0.0f;
//	deviceFP beamTotal2 = 0.0f;
//	//Integrate total beam power
//	beamTotal1 = 0.0f;
//	beamTotal2 = 0.0f;
//	for (j = 0; j < (*s).Nspace; ++j) {
//		beamTotal1 += modulusSquared((*s).workspace1[i + j *
//(*s).Nfreq]); 		beamTotal2 += modulusSquared((*s).workspace2[i + j *
//(*s).Nfreq]);
//	}
//	beamTotal1 *= 2.0f * LIGHTC * eps0<deviceFP>() * (*s).dx * (*s).dx *
//(*s).Nspace * (*s).dt * (*s).dt; 	beamTotal2 *= 2.0f * LIGHTC *
//eps0<deviceFP>() * (*s).dx * (*s).dx * (*s).Nspace * (*s).dt * (*s).dt;
//	//put the values into the output spectrum
//	(*s).gridPolarizationTime1[i] = beamTotal1;
//	(*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
//	(*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 +
//beamTotal2;
// };

// Calculate the energy spectrum after a 3D propagation
class totalSpectrum3DKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(int64_t i) const {
    deviceFP beamTotal1{};
    deviceFP beamTotal2{};
    // Integrate total beam power
    for (int64_t j{}; j < (*s).Nspace * (*s).Nspace2; ++j) {
      beamTotal1 += modulusSquared((*s).workspace1[i + j * (*s).Nfreq]);
      beamTotal2 += modulusSquared((*s).workspace2[i + j * (*s).Nfreq]);
    }

    // put the values into the output spectrum
    (*s).gridPolarizationTime1[i] = beamTotal1;
    (*s).gridPolarizationTime1[i + (*s).Nfreq] = beamTotal2;
    (*s).gridPolarizationTime1[i + 2 * (*s).Nfreq] = beamTotal1 + beamTotal2;
  }
};

// perform a Hankel transform by direct quadrature
// the offset radial grid allows the sum to be done with a midpoint method
// with no numerical effort and the rho=0 point is excluded from the grid
// this function is slow and order N^2 as it is not used in the core loop.
// the numerical accuracy of Hankel transforms that I've seen is relatively
// low due to Gibbs phenomena and I find the FFT-based propagation implemented
// below better for nonlinear phenomena. I might later use this for linear
// propagation in sequences however.
class hankelKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceFP *in;
  deviceFP *out;
  deviceFunction void operator()(const int64_t i) const {
    int64_t col = i / (*s).Ntime; // spatial coordinate
    deviceFP dk = constProd((deviceFP)(1.0 / vPi<deviceFP>()), 2) /
                  ((*s).dx * (*s).Nspace);
    out[i] = {};
    out[i + (*s).Ngrid] = {};
    deviceFP r0;
    deviceFP J0 = 1.0f;
    deviceFP k0 = col * dk;
    for (int64_t r{}; r < (*s).Nspace; ++r) {
      r0 = rhoInRadialSymmetry((*s).Nspace, r, (*s).dx);
      J0 = r0 * j0Device(r0 * k0);
      out[i] += J0 * in[r * (*s).Ntime + i % (*s).Ntime];
      out[i + (*s).Ngrid] +=
          J0 * in[r * (*s).Ntime + (*s).Ngrid + i % (*s).Ntime];
    }
    out[i] *= (*s).dx;
    out[i + (*s).Ngrid] *= (*s).dx;
  }
};

// inverse Hankel transform from the k-space back to the offset spatial grid
class inverseHankelKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceFP *in;
  deviceFP *out;
  deviceFunction void operator()(const int64_t i) const {
    int64_t col = i / (*s).Ntime; // spatial coordinate
    deviceFP dk = constProd((deviceFP)(1.0 / vPi<deviceFP>()), 2.0) /
                  ((*s).dx * (*s).Nspace);
    out[i] = {};
    out[i + (*s).Ngrid] = {};
    deviceFP r0 = rhoInRadialSymmetry((*s).Nspace, col, (*s).dx);
    deviceFP J0 = 1.0f;
    deviceFP k0 = col * dk;
    for (int64_t k{}; k < (*s).Nspace; ++k) {
      k0 = k * dk;
      J0 = k0 * j0Device(r0 * k0);
      out[i] += J0 * in[k * (*s).Ntime + i % (*s).Ntime];
      out[i + (*s).Ngrid] +=
          J0 * in[k * (*s).Ntime + (*s).Ngrid + i % (*s).Ntime];
    }
    out[i] *= 0.5f * dk / ((*s).Ntime);
    out[i + (*s).Ngrid] *= 0.5f * dk / ((*s).Ntime);
  }
};

// rotate the field around the propagation axis (basis change)
class rotateFieldKernel {
public:
  const deviceComplex *Ein1;
  const deviceComplex *Ein2;
  deviceComplex *Eout1;
  deviceComplex *Eout2;
  const deviceFP rotationAngle;
  deviceFunction void operator()(const int64_t i) const {
    Eout1[i] = deviceFPLib::cos(rotationAngle) * Ein1[i] -
               deviceFPLib::sin(rotationAngle) * Ein2[i];
    Eout2[i] = deviceFPLib::sin(rotationAngle) * Ein1[i] +
               deviceFPLib::cos(rotationAngle) * Ein2[i];
  }
};

class rotateField90Kernel {
public:
  const deviceComplex *Ein1;
  const deviceComplex *Ein2;
  deviceComplex *Eout1;
  deviceComplex *Eout2;
  deviceFunction void operator()(const int64_t i) const {
    Eout1[i] = -Ein2[i];
    Eout2[i] = Ein1[i];
  }
};

class rotateField180Kernel {
public:
  const deviceComplex *Ein1;
  const deviceComplex *Ein2;
  deviceComplex *Eout1;
  deviceComplex *Eout2;
  deviceFunction void operator()(const int64_t i) const {
    Eout1[i] = -Ein1[i];
    Eout2[i] = -Ein2[i];
  }
};

// apply a perfect polarizer to the field (same pattern as rotate with negative
// angle, and one rotated component is rotated back to orgin)
class polarizerKernel {
public:
  const deviceComplex *Ein1;
  const deviceComplex *Ein2;
  deviceComplex *Eout1;
  deviceComplex *Eout2;
  const deviceFP rotationAngle;
  deviceFunction void operator()(const int64_t i) const {
    deviceComplex projectedField = deviceFPLib::cos(rotationAngle) * Ein1[i] +
                                   deviceFPLib::sin(rotationAngle) * Ein2[i];
    Eout1[i] = deviceFPLib::cos(rotationAngle) * projectedField;
    Eout2[i] = deviceFPLib::sin(rotationAngle) * projectedField;
  }
};

// calculate the extra term in the Laplacian encountered in cylindrical
// coordinates (1/rho d/drho)
class radialLaplacianKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t j = i / (*s).Ntime;

    // zero at edges of grid
    [[unlikely]] if (j < 3 || j > ((*s).Nspace - 4)) {
      (*s).gridRadialLaplacian1[i] = {};
      (*s).gridRadialLaplacian2[i] = {};
    } else {
      const int64_t h = i - j * (*s).Ntime;
      const deviceFP *E1 = (*s).gridETime1 + h;
      const deviceFP *E2 = (*s).gridETime2 + h;
      if (j < (*s).Nspace / 2) {
        const deviceFP rhoFac =
            2.0f / ((*s).dx * (*s).dx *
                    (static_cast<deviceFP>((*s).Nspace / 2 - j) - 0.25f));
        const int64_t neighbors[6]{((*s).Nspace - j - 2) * (*s).Ntime,
                                   (j + 1) * (*s).Ntime,
                                   ((*s).Nspace - j - 1) * (*s).Ntime,
                                   ((*s).Nspace - j) * (*s).Ntime,
                                   (j - 1) * (*s).Ntime,
                                   ((*s).Nspace - j + 1) * (*s).Ntime};
        (*s).gridRadialLaplacian1[i] =
            rhoFac * (firstDerivativeStencil<deviceFP>(-3) * E1[neighbors[0]] +
                      firstDerivativeStencil<deviceFP>(-2) * E1[neighbors[1]] +
                      firstDerivativeStencil<deviceFP>(-1) * E1[neighbors[2]] +
                      firstDerivativeStencil<deviceFP>(1) * E1[neighbors[3]] +
                      firstDerivativeStencil<deviceFP>(2) * E1[neighbors[4]] +
                      firstDerivativeStencil<deviceFP>(3) * E1[neighbors[5]]);
        (*s).gridRadialLaplacian2[i] =
            rhoFac * (firstDerivativeStencil<deviceFP>(-3) * E2[neighbors[0]] +
                      firstDerivativeStencil<deviceFP>(-2) * E2[neighbors[1]] +
                      firstDerivativeStencil<deviceFP>(-1) * E2[neighbors[2]] +
                      firstDerivativeStencil<deviceFP>(1) * E2[neighbors[3]] +
                      firstDerivativeStencil<deviceFP>(2) * E2[neighbors[4]] +
                      firstDerivativeStencil<deviceFP>(3) * E2[neighbors[5]]);
      } else {
        const deviceFP rhoFac =
            2.0f / ((*s).dx * (*s).dx *
                    (static_cast<deviceFP>(j - (*s).Nspace / 2) + 0.25f));
        const int64_t neighbors[6]{((*s).Nspace - j + 1) * (*s).Ntime,
                                   (j - 1) * (*s).Ntime,
                                   ((*s).Nspace - j) * (*s).Ntime,
                                   ((*s).Nspace - j - 1) * (*s).Ntime,
                                   (j + 1) * (*s).Ntime,
                                   ((*s).Nspace - j - 2) * (*s).Ntime};
        (*s).gridRadialLaplacian1[i] =
            rhoFac * (firstDerivativeStencil<deviceFP>(-3) * E1[neighbors[0]] +
                      firstDerivativeStencil<deviceFP>(-2) * E1[neighbors[1]] +
                      firstDerivativeStencil<deviceFP>(-1) * E1[neighbors[2]] +
                      firstDerivativeStencil<deviceFP>(1) * E1[neighbors[3]] +
                      firstDerivativeStencil<deviceFP>(2) * E1[neighbors[4]] +
                      firstDerivativeStencil<deviceFP>(3) * E1[neighbors[5]]);
        (*s).gridRadialLaplacian2[i] =
            rhoFac * (firstDerivativeStencil<deviceFP>(-3) * E2[neighbors[0]] +
                      firstDerivativeStencil<deviceFP>(-2) * E2[neighbors[1]] +
                      firstDerivativeStencil<deviceFP>(-1) * E2[neighbors[2]] +
                      firstDerivativeStencil<deviceFP>(1) * E2[neighbors[3]] +
                      firstDerivativeStencil<deviceFP>(2) * E2[neighbors[4]] +
                      firstDerivativeStencil<deviceFP>(3) * E2[neighbors[5]]);
      }
    }
  }
};

class apertureFarFieldKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceFP radius;
  const deviceFP activationParameter;
  const deviceFP xOffset;
  const deviceFP yOffset;
  deviceFunction void operator()(int64_t i) const {
    const int64_t col = i / ((*s).Nfreq - 1);   // spatial coordinate
    const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
    i = j + col * (*s).Nfreq;
    const int64_t l = col / (*s).Nspace;
    const int64_t k = col - l * (*s).Nspace;

    // magnitude of k vector
    const deviceFP ko =
        constProd(twoPi<deviceFP>(), 1.0 / lightC<double>()) * j * (*s).fStep;

    // transverse wavevector being resolved
    const deviceFP dk1 =
        k * (*s).dk1 -
        (k >= ((int64_t)(*s).Nspace / 2)) *
            ((*s).dk1 * (int64_t)(*s).Nspace); // frequency grid in x direction
    deviceFP dk2 = {};
    if ((*s).is3D)
      dk2 = l * (*s).dk2 -
            (l >= ((int64_t)(*s).Nspace2 / 2)) *
                ((*s).dk2 *
                 (int64_t)(*s).Nspace2); // frequency grid in y direction

    // light that won't go the the farfield is immediately zero
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
    const deviceFP a =
        1.0f -
        (1.0f / (1.0f + deviceFPLib::exp(-activationParameter * (r - radius))));
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
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceFP radius;
  const deviceFP activationParameter;
  const deviceFP xOffset;
  const deviceFP yOffset;
  deviceFunction void operator()(int64_t i) const {
    const int64_t col = i / ((*s).Nfreq - 1);   // spatial coordinate
    const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
    i = j + col * (*s).Nfreq;
    const int64_t l = col / (*s).Nspace;
    const int64_t k = col - l * (*s).Nspace;

    // magnitude of k vector
    const deviceFP ko =
        constProd(twoPi<deviceFP>(), 1.0 / lightC<double>()) * j * (*s).fStep;

    // transverse wavevector being resolved
    const deviceFP dk1 =
        k * (*s).dk1 -
        (k >= ((int64_t)(*s).Nspace / 2)) *
            ((*s).dk1 * (int64_t)(*s).Nspace); // frequency grid in x direction
    deviceFP dk2 = {};
    if ((*s).is3D)
      dk2 = l * (*s).dk2 -
            (l >= ((int64_t)(*s).Nspace2 / 2)) *
                ((*s).dk2 *
                 (int64_t)(*s).Nspace2); // frequency grid in y direction

    // light that won't go the the farfield is immediately zero
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
    const deviceFP a =
        (1.0f / (1.0f + deviceFPLib::exp(-activationParameter * (r - radius))));
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
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceFP radius;
  const deviceFP activationParameter;
  const deviceFP xOffset;
  const deviceFP yOffset;
  deviceFunction void operator()(int64_t i) const {
    const int64_t col = i / ((*s).Nfreq - 1);   // spatial coordinate
    const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
    i = j + col * (*s).Nfreq;
    const int64_t k = col % (*s).Nspace;

    // magnitude of k vector
    deviceFP ko =
        constProd(twoPi<deviceFP>(), 1.0 / lightC<double>()) * j * (*s).fStep;

    // transverse wavevector being resolved
    deviceFP dk1 = constProd((deviceFP)2.0, 1.0 / vPi<double>()) * k /
                   ((*s).dx * (*s).Nspace);
    ; // frequency grid in x direction

    // light that won't go the the farfield is immediately zero
    if (dk1 * dk1 > ko * ko) {
      (*s).gridEFrequency1[i] = {};
      (*s).gridEFrequency2[i] = {};
      return;
    }

    const deviceFP theta1 = deviceFPLib::asin(dk1 / ko);
    const deviceFP a =
        1.0f -
        (1.0f / (1.0f + deviceFPLib::exp(-activationParameter *
                                         (deviceFPLib::abs(theta1) - radius))));
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
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceFP radius;
  const deviceFP activationParameter;
  const deviceFP xOffset;
  const deviceFP yOffset;
  deviceFunction void operator()(int64_t i) const {
    const int64_t col = i / ((*s).Nfreq - 1);   // spatial coordinate
    const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
    i = j + col * (*s).Nfreq;
    const int64_t k = col % (*s).Nspace;

    // magnitude of k vector
    deviceFP ko =
        constProd(twoPi<deviceFP>(), 1.0 / lightC<double>()) * j * (*s).fStep;

    // transverse wavevector being resolved
    deviceFP dk1 = constProd((deviceFP)2.0, 1.0 / vPi<double>()) * k /
                   ((*s).dx * (*s).Nspace);
    ; // frequency grid in x direction

    // light that won't go the the farfield is immediately zero
    if (dk1 * dk1 > ko * ko) {
      (*s).gridEFrequency1[i] = {};
      (*s).gridEFrequency2[i] = {};
      return;
    }

    const deviceFP theta1 = deviceFPLib::asin(dk1 / ko);
    const deviceFP a =
        (1.0f / (1.0f + deviceFPLib::exp(-activationParameter *
                                         (deviceFPLib::abs(theta1) - radius))));
    (*s).gridEFrequency1[i] *= a;
    (*s).gridEFrequency2[i] *= a;
    if (j == 1) {
      (*s).gridEFrequency1[i - 1] = deviceComplex{};
      (*s).gridEFrequency2[i - 1] = deviceComplex{};
    }
  }
};

// apply a spectral filter to the beam (full details in docs)
class filterKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFP f0;
  const deviceFP bandwidth;
  const int order;
  const deviceFP inBandAmplitude;
  const deviceFP outOfBandAmplitude;
  deviceFunction void operator()(int64_t i) const {
    const int64_t col = i / ((*s).Nfreq - 1);   // spatial coordinate
    const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
    i = j + col * (*s).Nfreq;

    deviceFP f = ((*s).fStep * j - f0) / bandwidth;
    for (int p = 1; p < order; p++) {
      f *= ((*s).fStep * j - f0) / bandwidth;
    }
    const deviceFP filterFunction =
        outOfBandAmplitude + inBandAmplitude * deviceFPLib::exp(-0.5f * f);
    (*s).gridEFrequency1[i] *= filterFunction;
    (*s).gridEFrequency2[i] *= filterFunction;
  }
};

class applyOpticKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceComplex *complexReflectivity;
  const bool applyX;
  const bool applyY;
  deviceFunction void operator()(int64_t i) const {
    const int64_t col = i / ((*s).Nfreq - 1);   // spatial coordinate
    const int64_t j = 1 + i % ((*s).Nfreq - 1); // frequency coordinate
    i = j + col * (*s).Nfreq;
    if (applyX)
      (*s).gridEFrequency1[i] *= complexReflectivity[j] / (*s).Ntime;
    else
      (*s).gridEFrequency1[i] /= (*s).Ntime;
    if (applyY)
      (*s).gridEFrequency2[i] *= complexReflectivity[j] / (*s).Ntime;
    else
      (*s).gridEFrequency2[i] /= (*s).Ntime;
    if (j == 1) {
      (*s).gridEFrequency1[i - 1] = {};
      (*s).gridEFrequency2[i - 1] = {};
    }
  }
};

// apply a lorentzian gain or loss in a certain cross-section of the beam.
//  amplitude - strength of the copy of the pulse applied
//  f0 - resonance frequency of the Lorentzian (THz)
//  gamma - linewidth (radHz)
//  radius - radius of the spot (m)
//  order - supergaussian order of the spot shape
class lorentzianSpotKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
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
    } else {
      j = i / ((*s).Nfreq - 1);
      i = h + j * ((*s).Nfreq);
      f = h * (*s).fStep;
      r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) +
                           0.25f * (*s).dx);
    }

    deviceFP w0 = twoPi<deviceFP>() * f0;
    deviceFP w = twoPi<deviceFP>() * f;
    deviceComplex lorentzian =
        gamma * w0 * amplitude /
        (w0 * w0 - w * w + deviceComplex(0.0f, gamma * w));
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

// Apply a (soft, possibly) aperture
class ApertureKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceFP radius;
  const deviceFP activationParameter;
  const deviceFP x_offset;
  const deviceFP y_offset;
  deviceFunction void operator()(int64_t i) const {
    const int64_t col = i / (*s).Ntime;
    const int64_t j = col % (*s).Nspace;
    const int64_t k = col / (*s).Nspace;
    deviceFP r;
    if ((*s).is3D) {
      const deviceFP x = ((*s).dx * (j - (*s).Nspace / 2.0f)) - x_offset;
      const deviceFP y = ((*s).dx * (k - (*s).Nspace2 / 2.0f)) - y_offset;
      r = deviceFPLib::hypot(x, y);
    } else if (s->isCylindric){
      r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) +
                           0.25f * (*s).dx);
    }
    else{
        r = deviceFPLib::abs( (*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) +
                             0.25f * (*s).dx - x_offset);
    }

    const deviceFP a =
        1.0f - (1.0f / (1.0f + deviceFPLib::exp(-activationParameter *
                                                (r - radius) / (*s).dx)));
    (*s).gridETime1[i] *= a;
    (*s).gridETime2[i] *= a;
  }
};

// apply a spatial phase corresponding to a parabolic mirror (on-axis)
class parabolicMirrorKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
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
    } else {
      r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) +
                           0.25f * (*s).dx);
    }

    const deviceComplex u = deviceLib::exp(
        deviceComplex(0.0f, w * r * r * (0.5f / focus) / lightC<deviceFP>()));

    (*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
    (*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
  }
};

// apply a spatial phase corresponding to a spherical mirror (on axis)
class sphericalMirrorKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
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
    } else {
      r = deviceFPLib::abs((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f) +
                           0.25f * (*s).dx);
    }
    const bool isNegative = ROC_in < 0.0f;
    const deviceFP ROC = deviceFPLib::abs(ROC_in);
    deviceComplex u = deviceComplex{};
    if (r >= ROC) {
      u = deviceComplex{};
    } else if (r > 0.5f * ROC) {
      u = deviceLib::exp(deviceComplex(
          0.0f, 2.0f * deviceFPLib::pow(-1.0f, isNegative) * w * ROC *
                    ((deviceFPLib::sqrt(1.0f - r * r / (ROC * ROC))) - 1.0f) /
                    lightC<deviceFP>()));
    } else {
      deviceFP ratio = r / ROC;
      ratio *= ratio;
      u = deviceLib::exp(deviceComplex(
          0.0f, 2.0f * deviceFPLib::pow(-1.0f, isNegative) * w * ROC *
                    (-0.5f * ratio - 0.125f * ratio * ratio -
                     0.0625f * ratio * ratio * ratio) /
                    lightC<deviceFP>()));
    }

    (*s).gridEFrequency1[i] = u * (*s).gridEFrequency1[i];
    (*s).gridEFrequency2[i] = u * (*s).gridEFrequency2[i];
    if (isComplexNaN((*s).gridEFrequency1[i]) ||
        isComplexNaN((*s).gridEFrequency2[i])) {
      (*s).gridEFrequency1[i] = deviceComplex{};
      (*s).gridEFrequency2[i] = deviceComplex{};
    }
  }
};

// apply linear propagation through a given medium to the fields
class correctFDTDAmplitudesKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t localIndex) const {
    int64_t i = localIndex;
    const int64_t h = 1 + i % ((*s).Nfreq - 1);
    const int64_t col = i / ((*s).Nfreq - 1);
    i = h + col * ((*s).Nfreq);
    const int64_t k = col / (*s).Nspace;
    const int64_t j = col - k * (*s).Nspace;

    const deviceFP kMagnitude =
        (h * (*s).fStep * twoPi<deviceFP>()) / lightC<deviceFP>();
    const deviceFP dk1 =
        j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace);
    const deviceFP dk2 =
        k * (*s).dk2 - (k >= ((*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2);
    if (dk2 * dk2 + dk1 * dk1 < kMagnitude * kMagnitude) {
      s->gridEFrequency1[i] *=
          (*s).fftNorm * kMagnitude /
          deviceFPLib::sqrt((kMagnitude - dk1) * (kMagnitude + dk1));
      s->gridEFrequency2[i] *=
          (*s).fftNorm * kMagnitude /
          deviceFPLib::sqrt((kMagnitude - dk2) * (kMagnitude + dk2));
    } else {
      s->gridEFrequency1[i] = {};
      s->gridEFrequency2[i] = {};
    }
    if (h == 1) {
      s->gridEFrequency1[i - 1] = {};
      s->gridEFrequency2[i - 1] = {};
    }
  }
};

class correctFDTDAmplitudesKernel2D {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t localIndex) const {
    int64_t i = localIndex;
    const int64_t h = 1 + i % ((*s).Nfreq - 1);
    const int64_t col = i / ((*s).Nfreq - 1);
    i = h + col * ((*s).Nfreq);
    const int64_t j = col % (*s).Nspace;

    const deviceFP kMagnitude =
        (h * (*s).fStep * twoPi<deviceFP>()) / lightC<deviceFP>();
    const deviceFP dk1 =
        j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace);
    if (dk1 * dk1 < kMagnitude * kMagnitude && h > 2) {
      s->gridEFrequency1[i] *=
          (*s).fftNorm * kMagnitude /
          deviceFPLib::sqrt((kMagnitude - dk1) * (kMagnitude + dk1));
      s->gridEFrequency2[i] *= (*s).fftNorm;
    } else {
      s->gridEFrequency1[i] = {};
      s->gridEFrequency2[i] = {};
    }
    if (h == 1) {
      s->gridEFrequency1[i - 1] = {};
      s->gridEFrequency2[i - 1] = {};
    }
  }
};

// apply linear propagation through a given medium to the fields
class applyLinearPropagationKernel {
public:
  const deviceFP *sellmeierCoefficients;
  const deviceFP thickness;
  const deviceParameterSet<deviceFP, deviceComplex> *s;
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

    // frequency being resolved by current thread
    const deviceFP f = h * (*s).fStep;
    const deviceFP omega = twoPi<deviceFP>() * f;
    deviceFP delta = findBirefringentCrystalIndex(s, sellmeierCoefficients,
                                                  localIndex, &ne, &no);
    if (s->axesNumber == 2 && i < s->NgridC)
      s->gridBiaxialDelta[i] = delta;

    const deviceFP dk1 =
        j * (*s).dk1 - (j >= ((*s).Nspace / 2)) * ((*s).dk1 * (*s).Nspace);
    deviceFP dk2 =
        k * (*s).dk2 - (k >= ((*s).Nspace2 / 2)) * ((*s).dk2 * (*s).Nspace2);
    if (!(*s).is3D)
      dk2 = 0.0f;
    sellmeierDevice(&n0, &n0o, sellmeierCoefficients, (*s).f0, crystalTheta,
                    crystalPhi, axesNumber, sellmeierType);
    if (isComplexNaN(ne) || isComplexNaN(no)) {
      ne = cOne<deviceComplex>();
      no = cOne<deviceComplex>();
    }

    const deviceComplex ke = ne * omega / lightC<deviceFP>();
    const deviceComplex ko = no * omega / lightC<deviceFP>();
    const deviceComplex k0 = n0 * omega / lightC<deviceFP>();

    deviceComplex ts = fourierPropagator(ke, dk1, dk2, k0.real(), thickness);
    deviceComplex tp = fourierPropagator(ko, dk1, dk2, k0.real(), thickness);
    if (((dk1 * dk1 + dk2 * dk2 > modulusSquared(ke)) ||
         (dk1 * dk1 + dk2 * dk2 > modulusSquared(ko))) &&
        thickness < 0.0f) {
      ts = deviceComplex{};
      tp = deviceComplex{};
    }
    if (isComplexNaN(ts))
      ts = deviceComplex{};
    if (isComplexNaN(tp))
      tp = deviceComplex{};
    (*s).gridEFrequency1[i] = ts * (*s).gridEFrequency1[i];
    (*s).gridEFrequency2[i] = tp * (*s).gridEFrequency2[i];
    if (isComplexNaN((*s).gridEFrequency1[i]) ||
        isComplexNaN((*s).gridEFrequency2[i])) {
      (*s).gridEFrequency1[i] = deviceComplex{};
      (*s).gridEFrequency2[i] = deviceComplex{};
    }
    if (h == 1) {
      (*s).gridEFrequency1[i - 1] = deviceComplex{};
      (*s).gridEFrequency2[i - 1] = deviceComplex{};
    }
  }
};

// prepare propagation constants for the simulation, when it is taking place on
// a Cartesian grid note that the sellmeier coefficients have extra values
// appended to the end to give info about the current simulation
class prepareCartesianGridsKernel {
public:
  const deviceFP *sellmeierCoefficients;
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t localIndex) const {
    const int64_t j = localIndex / ((*s).Nfreq - 1);       // spatial coordinate
    const int64_t k = 1 + (localIndex % ((*s).Nfreq - 1)); // temporal
                                                           // coordinate
    const int64_t i = k + j * (*s).Nfreq;

    // frequency being resolved by current thread
    const deviceFP f = k * (*s).fStep;

    // transverse wavevector being resolved
    const deviceFP dk =
        j * (*s).dk1 -
        (j >= ((*s).Nspace / 2)) *
            ((*s).dk1 * (*s).Nspace); // frequency grid in transverse direction
    deviceComplex ne, no;
    deviceFP d = findBirefringentCrystalIndex(s, sellmeierCoefficients,
                                              localIndex, &ne, &no);
    if (s->axesNumber == 2 && i < s->NgridC)
      s->gridBiaxialDelta[i] = d;
    // if the refractive index was returned weird,
    // then the index isn't valid, so set the propagator to zero for that
    // frequency
    if (minN(ne.real(), no.real()) < 0.9f || isComplexNaN(ne) ||
        isComplexNaN(no)) {
      (*s).gridPropagationFactor1[i] = {};
      (*s).gridPropagationFactor2[i] = {};
      (*s).gridPolarizationFactor1[i] = {};
      (*s).gridPolarizationFactor2[i] = {};
      return;
    }

    const deviceComplex k0 = deviceComplex(
        twoPi<deviceFP>() * (*s).n0.real() * f / lightC<deviceFP>(), 0.0f);
    const deviceComplex ke = twoPi<deviceFP>() * ne * f / lightC<deviceFP>();
    const deviceComplex ko = twoPi<deviceFP>() * no * f / lightC<deviceFP>();

    // chi11 factor, also multiplied by -sqrt(-1)!
    const deviceComplex chi11 = ((*s).isUsingMillersRule)
                                    ? deviceComplex((*s).chiLinear1[k].imag(),
                                                    -(*s).chiLinear1[k].real())
                                    : cOne<deviceComplex>();
    const deviceComplex chi12 = ((*s).isUsingMillersRule)
                                    ? deviceComplex((*s).chiLinear2[k].imag(),
                                                    -(*s).chiLinear2[k].real())
                                    : cOne<deviceComplex>();

    const deviceComplex kz1 =
        deviceLib::sqrt(ke - dk) * deviceLib::sqrt(ke + dk);
    const deviceComplex kz2 =
        deviceLib::sqrt(ko - dk) * deviceLib::sqrt(ko + dk);

    if (kz1.real() > 0.0f && kz2.real() > 0.0f && !isComplexNaN(k0) &&
        !isComplexNaN(kz1) && !isComplexNaN(kz2)) {
      (*s).gridPropagationFactor1[i] =
          fourierPropagator(ke, dk, deviceFP{}, k0.real(), 0.5f * (*s).h);
      (*s).gridPropagationFactor2[i] =
          fourierPropagator(ko, dk, deviceFP{}, k0.real(), 0.5f * (*s).h);
      (*s).gridPolarizationFactor1[i] =
          deviceLib::pow((*s).chiLinear1[k] + 1.0f, (deviceFP)0.25f) * chi11 *
          (twoPi<deviceFP>() * twoPi<deviceFP>() * f * f) /
          ((2.0f * lightC<deviceFP>() * lightC<deviceFP>() * kz1)) * (*s).h;
      (*s).gridPolarizationFactor2[i] =
          deviceLib::pow((*s).chiLinear2[k] + 1.0f, (deviceFP)0.25f) * chi12 *
          (twoPi<deviceFP>() * twoPi<deviceFP>() * f * f) /
          ((2.0f * lightC<deviceFP>() * lightC<deviceFP>() * kz2)) * (*s).h;
    } else {
      (*s).gridPropagationFactor1[i] = {};
      (*s).gridPropagationFactor2[i] = {};
      (*s).gridPolarizationFactor1[i] = {};
      (*s).gridPolarizationFactor2[i] = {};
    }

    if (isComplexNaN((*s).gridPropagationFactor1[i]) ||
        isComplexNaN((*s).gridPropagationFactor2[i]) ||
        isComplexNaN((*s).gridPolarizationFactor1[i]) ||
        isComplexNaN((*s).gridPolarizationFactor2[i])) {
      (*s).gridPropagationFactor1[i] = {};
      (*s).gridPropagationFactor2[i] = {};
      (*s).gridPolarizationFactor1[i] = {};
      (*s).gridPolarizationFactor2[i] = {};
    }
  }
};

// prepare propagation constants for the simulation, when it is taking place on
// a Cartesian grid note that the sellmeier coefficients have extra values
// appended to the end to give info about the current simulation
class prepare3DGridsKernel {
public:
  const deviceFP *sellmeierCoefficients;
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t localIndex) const {

    const int64_t col = localIndex / ((*s).Nfreq - 1);   // spatial coordinate
    const int64_t j = 1 + localIndex % ((*s).Nfreq - 1); // frequency coordinate
    const int64_t i = j + col * (*s).Nfreq;
    const int64_t k = col % (*s).Nspace;
    const int64_t l = col / (*s).Nspace;

    // frequency being resolved by current thread
    const deviceFP f = j * (*s).fStep;

    // transverse wavevector being resolved
    const deviceFP dk1 =
        k * (*s).dk1 -
        (k >= ((*s).Nspace / 2)) *
            ((*s).dk1 * (*s).Nspace); // frequency grid in x direction
    const deviceFP dk2 =
        l * (*s).dk2 -
        (l >= ((*s).Nspace2 / 2)) *
            ((*s).dk2 * (*s).Nspace2); // frequency grid in y direction

    deviceComplex ne, no;
    deviceFP d = findBirefringentCrystalIndex(s, sellmeierCoefficients,
                                              localIndex, &ne, &no);
    if (s->axesNumber == 2 && i < s->NgridC)
      s->gridBiaxialDelta[i] = d;
    if (minN(ne.real(), no.real()) < 0.9f || isComplexNaN(ne) ||
        isComplexNaN(no)) {
      (*s).gridPropagationFactor1[i] = {};
      (*s).gridPropagationFactor2[i] = {};
      (*s).gridPolarizationFactor1[i] = {};
      (*s).gridPolarizationFactor2[i] = {};
      return;
    }

    deviceComplex k0 = twoPi<deviceFP>() * (*s).n0 * f / lightC<deviceFP>();
    deviceComplex ke = twoPi<deviceFP>() * ne * f / lightC<deviceFP>();
    deviceComplex ko =
        (deviceFP)twoPi<deviceFP>() * no * f / lightC<deviceFP>();

    // chi11 factor, also multiplied by -sqrt(-1)!
    const deviceComplex chi11 = ((*s).isUsingMillersRule)
                                    ? deviceComplex((*s).chiLinear1[j].imag(),
                                                    -(*s).chiLinear1[j].real())
                                    : cOne<deviceComplex>();
    const deviceComplex chi12 = ((*s).isUsingMillersRule)
                                    ? deviceComplex((*s).chiLinear2[j].imag(),
                                                    -(*s).chiLinear2[j].real())
                                    : cOne<deviceComplex>();

    deviceComplex kz1 = deviceLib::sqrt(ke * ke - dk1 * dk1 - dk2 * dk2);
    deviceComplex kz2 = deviceLib::sqrt(ko * ko - dk1 * dk1 - dk2 * dk2);
    if (kz1.real() > 0.0f && kz2.real() > 0.0f) {
      (*s).gridPropagationFactor1[i] =
          fourierPropagator(ke, dk1, dk2, k0.real(), 0.5f * (*s).h);
      if (isnan(((*s).gridPropagationFactor1[i].real()))) {
        (*s).gridPropagationFactor1[i] = {};
      }

      (*s).gridPropagationFactor2[i] =
          fourierPropagator(ko, dk1, dk2, k0.real(), 0.5f * (*s).h);
      if (isnan(((*s).gridPropagationFactor2[i].real()))) {
        (*s).gridPropagationFactor2[i] = {};
      }

      (*s).gridPolarizationFactor1[i] =
          deviceLib::pow((*s).chiLinear1[j] + 1.0f, 0.25f) * chi11 *
          (twoPi<deviceFP>() * twoPi<deviceFP>() * f * f) /
          (2.0f * lightC<deviceFP>() * lightC<deviceFP>() * kz1) * (*s).h;
      (*s).gridPolarizationFactor2[i] =
          deviceLib::pow((deviceComplex)(*s).chiLinear2[j] + 1.0f, 0.25f) *
          chi12 * (twoPi<deviceFP>() * twoPi<deviceFP>() * f * f) /
          (2.0f * lightC<deviceFP>() * lightC<deviceFP>() * kz2) * (*s).h;
    } else {
      (*s).gridPropagationFactor1[i] = {};
      (*s).gridPropagationFactor2[i] = {};
      (*s).gridPolarizationFactor1[i] = {};
      (*s).gridPolarizationFactor2[i] = {};
    }

    if (isComplexNaN((*s).gridPropagationFactor1[i]) ||
        isComplexNaN((*s).gridPropagationFactor2[i]) ||
        isComplexNaN((*s).gridPolarizationFactor1[i]) ||
        isComplexNaN((*s).gridPolarizationFactor2[i])) {
      (*s).gridPropagationFactor1[i] = {};
      (*s).gridPropagationFactor2[i] = {};
      (*s).gridPolarizationFactor1[i] = {};
      (*s).gridPolarizationFactor2[i] = {};
    }
  }
};

// prepare the chi(1) arrays that will be needed in the simulation
class getChiLinearKernel {
public:
  deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceFP *sellmeierCoefficients;
  deviceFunction void operator()(const int64_t i) const {
    // frequency being resolved by current thread
    const deviceFP f = static_cast<deviceFP>(i) * (*s).fStep;
    deviceComplex ne, no;
    if (i < (*s).Ntime / 2) {
      sellmeierDevice(&ne, &no, sellmeierCoefficients, f, (*s).crystalTheta,
                      (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType,
                      false); // note: false to applySqrt, so it returns n^2

      (*s).chiLinear1[i] = ne - 1.0f;
      (*s).chiLinear2[i] = no - 1.0f;
      if (i > 0 && i < ((*s).Nfreq - 1) && (*s).chiLinear1[i].real() != 0.0f &&
          (*s).chiLinear2[i].real() != 0.0f) {
        (*s).inverseChiLinear1[i] = 1.0f / (*s).chiLinear1[i].real();
        (*s).inverseChiLinear2[i] = 1.0f / (*s).chiLinear2[i].real();
      } else {
        (*s).inverseChiLinear1[i] = {};
        (*s).inverseChiLinear2[i] = {};
      }

      (*s).fieldFactor1[i] =
          deviceFPLib::pow((*s).chiLinear1[i].real() + 1.0f, -0.25f);
      // account for the effective field strength in the medium (1/n)
      (*s).fieldFactor2[i] =
          deviceFPLib::pow((*s).chiLinear2[i].real() + 1.0f, -0.25f);
      if ((*s).isUsingMillersRule) {
        (*s).fieldFactor1[i] *= (*s).chiLinear1[i].real();
        (*s).fieldFactor2[i] *= (*s).chiLinear2[i].real();
      }
      (*s).fieldFactor1[i] *= (*s).fftNorm;
      (*s).fieldFactor2[i] *= (*s).fftNorm;
      if (isComplexNaN(ne) || isComplexNaN(no) ||
          isComplexNaN((*s).chiLinear1[i]) ||
          isComplexNaN((*s).chiLinear2[i]) || isnan((*s).fieldFactor1[i]) ||
          isnan((*s).fieldFactor2[i]) || isnan((*s).inverseChiLinear1[i]) ||
          isnan((*s).inverseChiLinear2[i]) || ne.real() < 0.9f ||
          no.real() < 0.9f || ne.imag() > 0.0f || no.imag() > 0.0f) {
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
      sellmeierDevice(&n0, &no, sellmeierCoefficients,
                      deviceFPLib::abs((*s).f0), (*s).crystalTheta,
                      (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);
      (*s).n0 = no;
      (*s).chiLinear1[(*s).Ntime / 2] = cOne<deviceComplex>();
      (*s).chiLinear2[(*s).Ntime / 2] = cOne<deviceComplex>();
      (*s).fieldFactor1[(*s).Ntime / 2] = {};
      (*s).fieldFactor2[(*s).Ntime / 2] = {};
      (*s).inverseChiLinear1[(*s).Ntime / 2] =
          1.0f / (*s).chiLinear1[(*s).Ntime / 2].real();
      (*s).inverseChiLinear2[(*s).Ntime / 2] =
          1.0f / (*s).chiLinear2[(*s).Ntime / 2].real();
    }

    // apply Miller's rule to nonlinear coefficients
    if (!(*s).isUsingMillersRule || i > 80) {
      return;
    }
    const deviceFP *referenceFrequencies = &sellmeierCoefficients[72];
    deviceFP chi11[7];

    for (int im = (i > 17) * 3; im < 7; ++im) {
      if (referenceFrequencies[im] == 0.0f) {
        chi11[im] = 100000.0f;
      } else {
        sellmeierDevice(&ne, &no, sellmeierCoefficients,
                        referenceFrequencies[im], (*s).crystalTheta,
                        (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType,
                        false); // note: false to applySqrt, so it returns n^2
        chi11[im] = no.real() - 1.0f;
      }
    }

    // normalize chi2 tensor values
    if (i < 18) {
      (*s).chi2Tensor[i] /= chi11[0] * chi11[1] * chi11[2];
      if (isnan((*s).chi2Tensor[i]))
        (*s).chi2Tensor[i] = {};
    }

    // normalize chi3 tensor values
    (*s).chi3Tensor[i] /= chi11[3] * chi11[4] * chi11[5] * chi11[6];
    if (isnan((*s).chi3Tensor[i]))
      (*s).chi3Tensor[i] = {};
  }
};

// prepare the propagation constants under the assumption of cylindrical
// symmetry of the beam
class prepareCylindricGridsKernel {
public:
  deviceFP *sellmeierCoefficients;
  deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(int64_t i) const {
    const int64_t j = i / ((*s).Nfreq - 1);     // spatial coordinate
    const int64_t k = 1 + i % ((*s).Nfreq - 1); // temporal coordinate
    i = static_cast<deviceFP>(k) + static_cast<deviceFP>(j) * (*s).Nfreq;
    const deviceComplex ii = deviceComplex(0.0f, 1.0f);
    // frequency being resolved by current thread
    const deviceFP f = -(static_cast<deviceFP>(k) * (*s).fStep);
    // transverse wavevector being resolved
    const deviceFP dk =
        static_cast<deviceFP>(j) * (*s).dk1 -
        static_cast<deviceFP>(j >= ((*s).Nspace / 2)) *
            ((*s).dk1 * (*s).Nspace); // frequency grid in transverse direction
    deviceComplex ne, no;
    deviceFP delta = sellmeierDevice(
        &ne, &no, sellmeierCoefficients, (*s).fStep * k, (*s).crystalTheta,
        (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);
    if ((*s).axesNumber == 2)
      (*s).gridBiaxialDelta[i] = delta;
    // if the refractive index was returned weird, then the index isn't valid,
    // so set the propagator to zero for that frequency
    if (minN(ne.real(), no.real()) < 0.95f || ne.real() > 6.0f ||
        no.real() > 6.0f || isComplexNaN(ne) || isComplexNaN(no)) {
      (*s).gridPropagationFactor1[i] = {};
      (*s).gridPropagationFactor2[i] = {};
      (*s).gridPolarizationFactor1[i] = {};
      (*s).gridPolarizationFactor2[i] = {};
      (*s).gridPropagationFactor1Rho1[i] = {};
      (*s).gridPropagationFactor1Rho2[i] = {};
      return;
    }

    const deviceComplex k0 = deviceComplex(
        twoPi<deviceFP>() * (*s).n0.real() * f / lightC<deviceFP>(), 0.0f);
    const deviceComplex ke = twoPi<deviceFP>() * ne * f / lightC<deviceFP>();
    const deviceComplex ko = twoPi<deviceFP>() * no * f / lightC<deviceFP>();

    const deviceComplex chi11 =
        ((*s).isUsingMillersRule) ? (*s).chiLinear1[k] : cOne<deviceComplex>();
    const deviceComplex chi12 =
        ((*s).isUsingMillersRule) ? (*s).chiLinear2[k] : cOne<deviceComplex>();

    if ((dk * dk < minN(ke.real() * ke.real() + ke.imag() * ke.imag(),
                        ko.real() * ko.real() + ko.imag() * ko.imag())) &&
        (*s).fieldFactor1[k] > 0.0f && (*s).fieldFactor2[k] > 0.0f) {
      (*s).gridPropagationFactor1[i] =
          0.5f * ii * (ke - k0 - dk * dk / (2.0f * ke.real())) * (*s).h;
      (*s).gridPropagationFactor1[i] =
          isComplexNaN((*s).gridPropagationFactor1[i])
              ? deviceComplex{}
              : deviceLib::exp((*s).gridPropagationFactor1[i]);
      (*s).gridPropagationFactor1Rho1[i] =
          (*s).fftNorm * ii * (*s).h / ((*s).fieldFactor1[k] * 2.0f * ke);
      if (isnan((deviceLib::abs((*s).gridPropagationFactor1Rho1[i] +
                                (*s).gridPropagationFactor1[i])))) {
        (*s).gridPropagationFactor1[i] = {};
        (*s).gridPropagationFactor1Rho1[i] = {};
      }

      (*s).gridPropagationFactor2[i] =
          0.5f * ii * (ko - k0 - dk * dk / (2.0f * ko.real())) * (*s).h;
      (*s).gridPropagationFactor2[i] =
          isComplexNaN((*s).gridPropagationFactor2[i])
              ? deviceComplex{}
              : deviceLib::exp((*s).gridPropagationFactor2[i]);
      (*s).gridPropagationFactor1Rho2[i] =
          (*s).fftNorm * ii * (*s).h / ((*s).fieldFactor2[k] * 2.0f * ko);
      if (isnan((deviceLib::abs((*s).gridPropagationFactor1Rho2[i] +
                                (*s).gridPropagationFactor2[i])))) {
        (*s).gridPropagationFactor2[i] = {};
        (*s).gridPropagationFactor1Rho2[i] = {};
      }

      // factor of 0.5 comes from deviceFPd grid size in cylindrical symmetry
      // mode after expanding the beam
      (*s).gridPolarizationFactor1[i] =
          0.5f *
          deviceLib::pow((deviceComplex)(*s).chiLinear1[k] + 1.0f, 0.25f) *
          chi11 * ii * (twoPi<deviceFP>() * f) /
          (2.0f * ne.real() * lightC<deviceFP>()) * (*s).h;
      (*s).gridPolarizationFactor2[i] =
          0.5f *
          deviceLib::pow((deviceComplex)(*s).chiLinear2[k] + 1.0f, 0.25f) *
          chi12 * ii * (twoPi<deviceFP>() * f) /
          (2.0f * no.real() * lightC<deviceFP>()) * (*s).h;
      if (isComplexNaN((*s).gridPolarizationFactor1[i]) ||
          isComplexNaN((*s).gridPolarizationFactor2[i])) {
        (*s).gridPolarizationFactor1[i] = {};
        (*s).gridPolarizationFactor2[i] = {};
      }
    } else {
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
  const deviceFP *a;
  const deviceFP f0;
  const deviceFP thickness;
  const int sellmeierEquation;
  deviceFP *phase1;
  deviceFunction void operator()(const int64_t i) const {
    // frequency being resolved by current thread
    const deviceFP f = i * df;
    // give phase shift relative to group velocity (approximated
    //  with low-order finite difference) so the pulse doesn't move
    deviceComplex ne, no, no0, n0p, n0m;
    sellmeierDevice<deviceFP, deviceComplex>(&ne, &no, a, f, {}, {}, 0,
                                             sellmeierEquation);
    sellmeierDevice<deviceFP, deviceComplex>(&ne, &no0, a, f0, {}, {}, 0,
                                             sellmeierEquation);
    sellmeierDevice<deviceFP, deviceComplex>(&ne, &n0p, a, f0 + 1.0e11f, {}, {},
                                             0, sellmeierEquation);
    sellmeierDevice<deviceFP, deviceComplex>(&ne, &n0m, a, f0 - 1.0e11f, {}, {},
                                             0, sellmeierEquation);
    no0 = no0 + f0 * (n0p - n0m) / 2.0e11f;
    phase1[i] = thickness * twoPi<deviceFP>() * f * (no.real() - no0.real()) /
                lightC<deviceFP>();
  }
};

// calculate the nonlinear polarization, after FFT to get the field
// in the time domain
class nonlinearPolarizationKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t i) const {
    const deviceFP Ex = (*s).gridETime1[i];
    const deviceFP Ey = (*s).gridETime2[i];

    if ((*s).nonlinearSwitches.hasChi2 || (*s).nonlinearSwitches.hasChi3) {

      deviceFP P[3]{};

      // rotate field into crystal frame
      const deviceFP E3[3] = {
          (*s).rotationForward[0] * Ex + (*s).rotationForward[1] * Ey,
          (*s).rotationForward[3] * Ex + (*s).rotationForward[4] * Ey,
          (*s).rotationForward[6] * Ex + (*s).rotationForward[7] * Ey};

      // second order nonlinearity, element-by-element in the reduced tensor
      // note that for historical reasons, chi2 is column-order and chi3 is
      // row-order...
      if ((*s).nonlinearSwitches.hasChi2) {
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
      if ((*s).nonlinearSwitches.hasChi3 &&
          !((*s).nonlinearSwitches.assumeCentrosymmetric)) {
        // loop over tensor element X_abcd
        // i hope the compiler unrolls this, but no way am I writing that out by
        // hand
        for (auto a = 0; a < 3; ++a) {
          for (auto b = 0; b < 3; ++b) {
            for (auto c = 0; c < 3; ++c) {
              for (auto d = 0; d < 3; ++d) {
                P[d] += (*s).chi3Tensor[a + 3 * b + 9 * c + 27 * d] * E3[a] *
                        E3[b] * E3[c];
              }
            }
          }
        }
      }

      (*s).gridPolarizationTime1[i] = (*s).rotationBackward[0] * P[0] +
                                      (*s).rotationBackward[1] * P[1] +
                                      (*s).rotationBackward[2] * P[2];
      (*s).gridPolarizationTime2[i] = (*s).rotationBackward[3] * P[0] +
                                      (*s).rotationBackward[4] * P[1] +
                                      (*s).rotationBackward[5] * P[2];

      // using only one value of chi3, under assumption of centrosymmetry when
      // nonlinearSwitches[1]==2
      if ((*s).nonlinearSwitches.assumeCentrosymmetric) {
        const deviceFP Esquared = (*s).chi3Tensor[0] * (Ex * Ex + Ey * Ey);
        (*s).gridPolarizationTime1[i] += Ex * Esquared;
        (*s).gridPolarizationTime2[i] += Ey * Esquared;
      }
    } else {
      const deviceFP Esquared = (*s).chi3Tensor[0] * (Ex * Ex + Ey * Ey);
      (*s).gridPolarizationTime1[i] = Ex * Esquared;
      (*s).gridPolarizationTime2[i] = Ey * Esquared;
    }
    if ((*s).isCylindric)
      expandCylindricalBeamDevice(
          s, i, (*s).gridRadialLaplacian1 + (*s).Ngrid * 4 * (*s).hasPlasma,
          (*s).gridPolarizationTime1, (*s).gridPolarizationTime2);
  }
};

// Plasma response with time-dependent carrier density
// This polarization needs a different factor in the nonlinear wave equation
// to account for the integration
// equation for the plasma current:
// J_drude(t) = (e/m)*exp(-gamma*t)*\int_-infty^t dt' exp(gamma*t)*N(t)*E(t)
// J_absorption(t) = (beta*E^2)^(Nphot-1)*E
// The energy factor sets the amplitude based on the bandgap of the crystal,
// such that the integral of E*J_absorption is carrier density applied in 3
// parts due to the different patterns of calculating ionization rate (A) and
// integrating over trajectories (B). Added to RK4 propagation array afterwards
// (C).
class plasmaCurrentKernel_twoStage_A {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t i) const {
    const int pMax = (*s).plasmaParameters.fieldExponent;

    // save values in workspaces, casting to deviceFP
    deviceFP *dN = (deviceFP *)(*s).workspace1;
    const deviceFP Esquared = (*s).plasmaParameters.nonlinearAbsorption *
                              ((*s).gridETime1[i] * (*s).gridETime1[i] +
                               (*s).gridETime2[i] * (*s).gridETime2[i]);
    deviceFP EtoThe2N = Esquared;
    for (int p = 0; p < pMax; ++p) {
      EtoThe2N *= Esquared;
    }
    (*s).gridPolarizationTime1[i] = EtoThe2N * (*s).gridETime1[i];
    (*s).gridPolarizationTime2[i] = EtoThe2N * (*s).gridETime2[i];
    dN[i] = (*s).plasmaParameters.energyFactor *
            ((*s).gridPolarizationTime1[i] * (*s).gridETime1[i] +
             (*s).gridPolarizationTime2[i] * (*s).gridETime2[i]);
  }
};

class plasmaCurrentKernel_twoStage_B {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t j = (i) * (*s).Ntime;
    const deviceFP *expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
    const deviceFP *dN =
        (j % (*s).Ngrid) + reinterpret_cast<deviceFP *>((*s).workspace1);
    const deviceFP *E = &(*s).gridETime1[j];
    deviceFP *P = &(*s).gridPolarizationTime1[j];
    SimpsonIntegrator<deviceFP> N(s->plasmaParameters.initialDensity);
    SimpsonIntegrator<deviceFP> integral;
    for (int64_t k{}; k < (*s).Ntime; ++k) {
      N.step(dN[k]);
      deviceFP plasmaFactor = N.current_integral() * (*s).expGammaT[k] *
                              (*s).plasmaParameters.integrationFactor;
      integral.step(plasmaFactor * E[k]);
      P[k] += expMinusGammaT[k] * integral.current_integral();
    }
  }
};

class plasmaCurrentKernel_twoStage_B_simultaneous {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t j = (i) * (*s).Ntime;
    const deviceFP *expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
    const deviceFP *dN = j + reinterpret_cast<deviceFP *>((*s).workspace1);
    const deviceFP *E = &(*s).gridETime1[j];
    const deviceFP *Ey = E + (*s).Ngrid;
    deviceFP *P = &(*s).gridPolarizationTime1[j];
    deviceFP *P2 = P + (*s).Ngrid;

    SimpsonIntegrator<deviceFP> N(s->plasmaParameters.initialDensity);
    SimpsonIntegrator<deviceFP> integral_x;
    SimpsonIntegrator<deviceFP> integral_y;
    for (int64_t k{}; k < (*s).Ntime; ++k) {
      N.step(dN[k]);
      deviceFP plasmaFactor = N.current_integral() * (*s).expGammaT[k] *
                              (*s).plasmaParameters.integrationFactor;
      integral_x.step(plasmaFactor * E[k]);
      integral_y.step(plasmaFactor * Ey[k]);
      P[k] += expMinusGammaT[k] * integral_x.current_integral();
      P2[k] += expMinusGammaT[k] * integral_y.current_integral();
    }
  }
};

class plasmaCurrentKernel_SaveOutput {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t j = (i) * (*s).Ntime;
    const deviceFP *expMinusGammaT = &(*s).expGammaT[(*s).Ntime];
    const deviceFP *dN = j + reinterpret_cast<deviceFP *>((*s).workspace1);
    const deviceFP *E = &(*s).gridETime1[j];
    const deviceFP *Ey = E + (*s).Ngrid;
    deviceFP *P = &(*s).gridPolarizationTime1[j];
    deviceFP *P2 = P + (*s).Ngrid;
    deviceFP *Ndest = reinterpret_cast<deviceFP *>((*s).gridEFrequency1);
    Ndest += j;

    SimpsonIntegrator<deviceFP> N(s->plasmaParameters.initialDensity);
    SimpsonIntegrator<deviceFP> integral_x;
    SimpsonIntegrator<deviceFP> integral_y;
    for (int64_t k{}; k < (*s).Ntime; ++k) {
      N.step(dN[k]);
      deviceFP plasmaFactor = N.current_integral() * (*s).expGammaT[k] *
                              (*s).plasmaParameters.integrationFactor;
      integral_x.step(plasmaFactor * E[k]);
      integral_y.step(plasmaFactor * Ey[k]);
      P[k] += expMinusGammaT[k] * integral_x.current_integral();
      P2[k] += expMinusGammaT[k] * integral_y.current_integral();
      Ndest[k] = N.current_integral();
    }
  }
};

class updateKwithPolarizationKernelCylindric {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
    const int64_t radIndex = i / ((*sP).Nfreq - 1);
    const int64_t gridIndex = freqIndex + radIndex * ((*sP).Nfreq);
    const int64_t polIndex =
        freqIndex +
        (radIndex + ((radIndex > ((*sP).Nspace / 2))) * (*sP).Nspace) *
            (*sP).Nfreq;
    (*sP).k1[gridIndex] +=
        (*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[polIndex];
    (*sP).k2[gridIndex] +=
        (*sP).gridPolarizationFactor2[gridIndex] * (*sP).workspace2P[polIndex];
  }
};

class updateKwithPlasmaKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    int64_t h = 1 + i % ((*sP).Nfreq - 1);   // temporal coordinate
    const int64_t j = i / ((*sP).Nfreq - 1); // spatial coordinate
    const deviceFP jfac = -1.0f / (twoPi<deviceFP>() * h * (*sP).fStep);
    h += j * (*sP).Nfreq;
    (*sP).k1[h] +=
        deviceComplex{-jfac * (*sP).gridPolarizationFactor1[h].imag(),
                      jfac * (*sP).gridPolarizationFactor1[h].real()} *
        (*sP).workspace1[h] * (*sP).inverseChiLinear1[h % ((*sP).Nfreq)];
    (*sP).k2[h] +=
        deviceComplex{-jfac * (*sP).gridPolarizationFactor2[h].imag(),
                      jfac * (*sP).gridPolarizationFactor2[h].real()} *
        (*sP).workspace2P[h] * (*sP).inverseChiLinear2[h % ((*sP).Nfreq)];
  }
};

class updateKwithPlasmaKernelCylindric {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t spaceIndex = i / ((*sP).Nfreq - 1);
    int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
    const deviceFP jfac = -1.0f / (twoPi<deviceFP>() * freqIndex * (*sP).fStep);
    const int64_t gridIndex = freqIndex + spaceIndex * ((*sP).Nfreq);
    int64_t fftIndex =
        freqIndex +
        (spaceIndex + ((spaceIndex > ((*sP).Nspace / 2))) * (*sP).Nspace) *
            (*sP).Nfreq;
    freqIndex = gridIndex % (*sP).Nfreq;
    (*sP).k1[gridIndex] +=
        deviceComplex{-jfac * (*sP).gridPolarizationFactor1[gridIndex].imag(),
                      jfac * (*sP).gridPolarizationFactor1[gridIndex].real()} *
        (*sP).workspace1[fftIndex] * (*sP).inverseChiLinear1[freqIndex];
    (*sP).k2[gridIndex] +=
        deviceComplex{-jfac * (*sP).gridPolarizationFactor2[gridIndex].imag(),
                      jfac * (*sP).gridPolarizationFactor2[gridIndex].real()} *
        (*sP).workspace2P[fftIndex] * (*sP).inverseChiLinear2[freqIndex];

    fftIndex += 4 * (*sP).NgridC;
    (*sP).k1[gridIndex] +=
        (*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[fftIndex];
    (*sP).k2[gridIndex] +=
        (*sP).gridPolarizationFactor2[gridIndex] * (*sP).workspace2P[fftIndex];
  }
};

class biaxialRotationKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceComplex *x;
  bool backwards;
  deviceFunction void operator()(const int64_t gridIndex) const {
    deviceComplex *y = x + sP->NgridC;
    const int64_t freqIndex =
        1 + gridIndex % ((*sP).Nfreq - 1); // frequency coordinate
    const int64_t i =
        freqIndex + (gridIndex / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);

    deviceFP cosDelta = deviceFPLib::cos(sP->gridBiaxialDelta[i]);
    deviceFP sinDelta = deviceFPLib::sin(sP->gridBiaxialDelta[i]);
    if (backwards)
      sinDelta = -sinDelta;

    deviceComplex newX = cosDelta * x[i] - sinDelta * y[i];
    y[i] = sinDelta * x[i] + cosDelta * y[i];
    x[i] = newX;
  }
};

// Slightly different kernels for the four stages of RK4.
// They used to be one big kernel with a switch case
// but this has slightly better utilization.
class rkKernel0 {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(int64_t gridIndex) const {
    const int64_t freqIndex =
        1 + gridIndex % ((*sP).Nfreq - 1); // frequency coordinate
    gridIndex = freqIndex + (gridIndex / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
    (*sP).k1[gridIndex] +=
        (*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[gridIndex];
    (*sP).gridEFrequency1Next1[gridIndex] =
        (*sP).gridPropagationFactor1[gridIndex] *
        (*sP).gridPropagationFactor1[gridIndex] *
        (sixth<deviceFP>() * (*sP).k1[gridIndex] +
         (*sP).gridEFrequency1[gridIndex]);
    const deviceFP ff = (gridIndex > (*sP).NgridC)
                            ? (*sP).fieldFactor2[freqIndex]
                            : (*sP).fieldFactor1[freqIndex];
    (*sP).workspace1[gridIndex] =
        ff * (*sP).gridPropagationFactor1[gridIndex] *
        ((*sP).gridEFrequency1[gridIndex] + 0.5f * (*sP).k1[gridIndex]);
    (*sP).k1[gridIndex] = {};
    [[unlikely]] if (freqIndex == 1)
      (*sP).workspace1[gridIndex - 1] = {};
  }
};

class rkKernel1 {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1); // frequency coordinate
    const int64_t gridIndex =
        freqIndex + (i / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
    (*sP).k1[gridIndex] +=
        (*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[gridIndex];
    (*sP).gridEFrequency1Next1[gridIndex] =
        (*sP).gridEFrequency1Next1[gridIndex] +
        (*sP).gridPropagationFactor1[gridIndex] * (deviceFP)third<deviceFP>() *
            (*sP).k1[gridIndex];
    const deviceFP ff = (gridIndex > (*sP).NgridC)
                            ? (*sP).fieldFactor2[freqIndex]
                            : (*sP).fieldFactor1[freqIndex];
    (*sP).workspace1[gridIndex] =
        ff * ((*sP).gridPropagationFactor1[gridIndex] *
                  (*sP).gridEFrequency1[gridIndex] +
              0.5f * (*sP).k1[gridIndex]);
    (*sP).k1[gridIndex] = {};
    [[unlikely]] if (freqIndex == 1)
      (*sP).workspace1[gridIndex - 1] = {};
  }
};

class rkKernel2 {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1); // frequency coordinate
    const int64_t gridIndex =
        freqIndex + (i / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
    (*sP).k1[gridIndex] +=
        (*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[gridIndex];
    (*sP).gridEFrequency1Next1[gridIndex] =
        (*sP).gridEFrequency1Next1[gridIndex] +
        (*sP).gridPropagationFactor1[gridIndex] * (deviceFP)third<deviceFP>() *
            (*sP).k1[gridIndex];
    const deviceFP ff = (gridIndex > (*sP).NgridC)
                            ? (*sP).fieldFactor2[freqIndex]
                            : (*sP).fieldFactor1[freqIndex];
    (*sP).workspace1[gridIndex] =
        ff * ((*sP).gridPropagationFactor1[gridIndex] *
              ((*sP).gridPropagationFactor1[gridIndex] *
                   (*sP).gridEFrequency1[gridIndex] +
               (*sP).k1[gridIndex]));
    (*sP).k1[gridIndex] = {};
    [[unlikely]] if (freqIndex == 1)
      (*sP).workspace1[gridIndex - 1] = {};
  }
};

class rkKernel3 {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1); // frequency coordinate
    const int64_t gridIndex =
        freqIndex + (i / ((*sP).Nfreq - 1)) * ((*sP).Nfreq);
    (*sP).k1[gridIndex] +=
        (*sP).gridPolarizationFactor1[gridIndex] * (*sP).workspace1[gridIndex];
    (*sP).gridEFrequency1[gridIndex] = (*sP).gridEFrequency1Next1[gridIndex] +
                                       sixth<deviceFP>() * (*sP).k1[gridIndex];
    const deviceFP ff = (gridIndex > (*sP).NgridC)
                            ? (*sP).fieldFactor2[freqIndex]
                            : (*sP).fieldFactor1[freqIndex];
    (*sP).workspace1[gridIndex] = ff * (*sP).gridEFrequency1[gridIndex];
    (*sP).k1[gridIndex] = {};
    [[unlikely]] if (freqIndex == 1)
      (*sP).workspace1[gridIndex - 1] = {};
  }
};

// Kernels for symmetry around z axis use a different form, adding the radial
// Laplacian instead of the nonlinear polarization
class rkKernel0Cylindric {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
    const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * (*sP).Nfreq;
    (*sP).k1[gridIndex] += (*sP).gridPropagationFactor1Rho1[gridIndex] *
                           (*sP).workspace1[gridIndex];
    (*sP).gridEFrequency1Next1[gridIndex] =
        (*sP).gridPropagationFactor1[gridIndex] *
        (*sP).gridPropagationFactor1[gridIndex] *
        (sixth<deviceFP>() * (*sP).k1[gridIndex] +
         (*sP).gridEFrequency1[gridIndex]);
    const deviceFP ff = (gridIndex > (*sP).NgridC)
                            ? (*sP).fieldFactor2[freqIndex]
                            : (*sP).fieldFactor1[freqIndex];
    (*sP).workspace1[gridIndex] =
        ff * ((*sP).gridPropagationFactor1[gridIndex] *
              ((*sP).gridEFrequency1[gridIndex] + 0.5f * (*sP).k1[gridIndex]));
    (*sP).k1[gridIndex] = {};
    [[unlikely]] if (freqIndex == 1)
      (*sP).workspace1[gridIndex - 1] = {};
  }
};

class rkKernel1Cylindric {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
    const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * (*sP).Nfreq;
    (*sP).k1[gridIndex] += (*sP).gridPropagationFactor1Rho1[gridIndex] *
                           (*sP).workspace1[gridIndex];
    (*sP).gridEFrequency1Next1[gridIndex] =
        (*sP).gridEFrequency1Next1[gridIndex] +
        (*sP).gridPropagationFactor1[gridIndex] * third<deviceFP>() *
            (*sP).k1[gridIndex];
    const deviceFP ff = (gridIndex > (*sP).NgridC)
                            ? (*sP).fieldFactor2[freqIndex]
                            : (*sP).fieldFactor1[freqIndex];
    (*sP).workspace1[gridIndex] =
        ff * ((*sP).gridPropagationFactor1[gridIndex] *
                  (*sP).gridEFrequency1[gridIndex] +
              0.5f * (*sP).k1[gridIndex]);
    (*sP).k1[gridIndex] = {};
    [[unlikely]] if (freqIndex == 1)
      (*sP).workspace1[gridIndex - 1] = {};
  }
};

class rkKernel2Cylindric {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
    const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * (*sP).Nfreq;
    (*sP).k1[gridIndex] += (*sP).gridPropagationFactor1Rho1[gridIndex] *
                           (*sP).workspace1[gridIndex];
    (*sP).gridEFrequency1Next1[gridIndex] =
        (*sP).gridEFrequency1Next1[gridIndex] +
        (*sP).gridPropagationFactor1[gridIndex] * third<deviceFP>() *
            (*sP).k1[gridIndex];
    const deviceFP ff = (gridIndex > (*sP).NgridC)
                            ? (*sP).fieldFactor2[freqIndex]
                            : (*sP).fieldFactor1[freqIndex];
    (*sP).workspace1[gridIndex] =
        ff * ((*sP).gridPropagationFactor1[gridIndex] *
              ((*sP).gridPropagationFactor1[gridIndex] *
                   (*sP).gridEFrequency1[gridIndex] +
               (*sP).k1[gridIndex]));
    (*sP).k1[gridIndex] = {};
    [[unlikely]] if (freqIndex == 1)
      (*sP).workspace1[gridIndex - 1] = {};
  }
};

class rkKernel3Cylindric {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *sP;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t freqIndex = 1 + i % ((*sP).Nfreq - 1);
    const int64_t gridIndex = freqIndex + (i / ((*sP).Nfreq - 1)) * (*sP).Nfreq;
    (*sP).k1[gridIndex] += (*sP).gridPropagationFactor1Rho1[gridIndex] *
                           (*sP).workspace1[gridIndex];
    (*sP).gridEFrequency1[gridIndex] = (*sP).gridEFrequency1Next1[gridIndex] +
                                       sixth<deviceFP>() * (*sP).k1[gridIndex];
    const deviceFP ff = (gridIndex > (*sP).NgridC)
                            ? (*sP).fieldFactor2[freqIndex]
                            : (*sP).fieldFactor1[freqIndex];
    (*sP).workspace1[gridIndex] = ff * (*sP).gridEFrequency1[gridIndex];
    (*sP).k1[gridIndex] = {};
    [[unlikely]] if (freqIndex == 1)
      (*sP).workspace1[gridIndex - 1] = {};
  }
};

class maxwellRKkernel0 {
public:
  const maxwell3D *s;
  const int64_t t;
  deviceFunction void operator()(const int64_t i) const {
    maxwellKPoint<deviceFP> k =
        maxwellDerivativeTerms(s, i, s->Egrid, s->Hgrid);
    maxwellCurrentTerms(s, i, t, false, s->Egrid, s->materialGrid, k, 0);
    s->EgridEstimate[i] = s->Egrid[i] + k.kE * (s->tStep * 0.5f);
    s->EgridNext[i] = s->Egrid[i] + k.kE * (sixth<deviceFP>() * s->tStep);
    s->HgridEstimate[i] = s->Hgrid[i] + k.kH * (s->tStep * 0.5f);
    s->HgridNext[i] = s->Hgrid[i] + k.kH * (sixth<deviceFP>() * s->tStep);
  }
};
class maxwellRKkernel1 {
public:
  const maxwell3D *s;
  const int64_t t;
  deviceFunction void operator()(const int64_t i) const {
    maxwellKPoint<deviceFP> k =
        maxwellDerivativeTerms(s, i, s->EgridEstimate, s->HgridEstimate);
    maxwellCurrentTerms(s, i, t, true, s->EgridEstimate,
                        s->materialGridEstimate, k, 1);
    s->EgridEstimate2[i] = s->Egrid[i] + k.kE * (s->tStep * 0.5f);
    s->EgridNext[i] += k.kE * (third<deviceFP>() * s->tStep);
    s->HgridEstimate2[i] = s->Hgrid[i] + k.kH * (s->tStep * 0.5f);
    s->HgridNext[i] += k.kH * (third<deviceFP>() * s->tStep);
  }
};
class maxwellRKkernel2 {
public:
  const maxwell3D *s;
  const int64_t t;
  deviceFunction void operator()(const int64_t i) const {
    maxwellKPoint<deviceFP> k =
        maxwellDerivativeTerms(s, i, s->EgridEstimate2, s->HgridEstimate2);
    maxwellCurrentTerms(s, i, t, true, s->EgridEstimate2,
                        s->materialGridEstimate2, k, 2);
    s->EgridEstimate[i] = s->Egrid[i] + k.kE * s->tStep;
    s->EgridNext[i] += k.kE * (third<deviceFP>() * s->tStep);
    s->HgridEstimate[i] = s->Hgrid[i] + k.kH * s->tStep;
    s->HgridNext[i] += k.kH * (third<deviceFP>() * s->tStep);
  }
};
class maxwellRKkernel3 {
public:
  const maxwell3D *s;
  const int64_t t;
  deviceFunction void operator()(const int64_t i) const {
    maxwellKPoint<deviceFP> k =
        maxwellDerivativeTerms(s, i, s->EgridEstimate, s->HgridEstimate);
    maxwellCurrentTerms(s, i, t + 1, false, s->EgridEstimate,
                        s->materialGridEstimate, k, 3);
    s->Egrid[i] = s->EgridNext[i] + k.kE * (sixth<deviceFP>() * s->tStep);
    s->Hgrid[i] = s->HgridNext[i] + k.kH * (sixth<deviceFP>() * s->tStep);
  }
};

// store the field in the observation plane in the in/out Ex and Ey arrays
class maxwellSampleGrid {
public:
  maxwell3D *s;
  int64_t time;
  deviceFunction void operator()(int64_t i) const {
    int64_t gridIndex = i * s->Nz + s->observationPoint;
    s->inOutEx[i * s->NtIO + time] = s->Egrid[gridIndex].x;
    s->inOutEy[i * s->NtIO + time] = s->Egrid[gridIndex].y;
  }
};

class maxwellSetInitialCarrierDensity {
public:
  maxwell3D *s;
  deviceFunction void operator()(int64_t i) const {
    const int64_t xIndex = i / s->Nz;
    const int64_t zIndex = i - s->Nz * xIndex;
    bool solveMaterialEquations =
        s->hasMaterialMap
            ? s->materialMap[i] > 0
            : zIndex >= s->materialStart && zIndex < s->materialStop;

    if (solveMaterialEquations) {
      const int64_t oscillatorIndex =
          s->hasMaterialMap
              ? s->oscillatorIndexMap[i] * s->Noscillators
              : (zIndex - s->materialStart) * s->Noscillators +
                    xIndex * (s->materialStop - s->materialStart) *
                        s->Noscillators;
      const int oscillatorType = s->hasMaterialMap ? s->materialMap[i] - 1 : 0;
      if (s->hasPlasma[oscillatorType]) {
        const deviceFP rho =
            elCharge<deviceFP>() * s->startingCarriers[oscillatorType];
        s->materialGrid[oscillatorIndex + s->Noscillators - 1].P.x = rho;
        s->materialGridEstimate[oscillatorIndex + s->Noscillators - 1].P.x =
            rho;
        s->materialGridEstimate2[oscillatorIndex + s->Noscillators - 1].P.x =
            rho;
        s->materialGridNext[oscillatorIndex + s->Noscillators - 1].P.x = rho;
      }
    }
  }
};

class beamNormalizeKernel {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  const deviceFP *rawSum;
  deviceFP *field;
  const deviceFP pulseEnergy;
  deviceFunction void operator()(const int64_t i) const {
    field[i] *=
        deviceFPLib::sqrt(pulseEnergy / ((deviceFP)(*s).Ntime * (*rawSum)));
  }
};

class addDoubleArraysKernel {
public:
  deviceFP *A;
  deviceFP *B;
  deviceFunction void operator()(const int64_t i) const { A[i] += B[i]; }
};

// crease a pulse on the grid for the 2D modes.
//  Note that normalization of the 2D mode assumes radial symmetry (i.e. that
//  it's
// a gaussian beam, not an infinite plane wave, which would have zero amplitude
// for finite energy).
class beamGenerationKernel2D {
public:
  deviceComplex *field;
  const Pulse<deviceFP> *p;
  deviceFP *pulseSum;
  deviceParameterSet<deviceFP, deviceComplex> *s;
  const bool hasLoadedField;
  const deviceComplex *loadedField;
  const deviceFP *materialPhase;
  const deviceFP *sellmeierCoefficients;
  deviceFunction void operator()(int64_t i) const {
    const int64_t h = 1 + i % ((*s).Nfreq - 1);
    const int64_t j = i / ((*s).Nfreq - 1);
    i = h + j * ((*s).Nfreq);
    const deviceFP f = h * (*s).fStep;
    const deviceFP w = twoPi<deviceFP>() * (f - (*p).frequency);

    // supergaussian pulse spectrum, if no input pulse specified
    deviceComplex specfac = deviceComplex(
        -deviceFPLib::pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder),
        0.0f);

    if (isComplexNaN(materialPhase[h])) {
      field[i] = deviceComplex{};
      field[i + (*s).NgridC] = deviceComplex{};
      return;
    }
    deviceComplex specphase = deviceComplex(
        0.0f,
        -((*p).cep +
          twoPi<deviceFP>() * f * ((*p).delay - 0.5f * (*s).dt * (*s).Ntime) +
          0.5f * (*p).gdd * w * w + sixth<deviceFP>() * (*p).tod * w * w * w +
          materialPhase[h]));
    specfac = deviceLib::exp(specfac + specphase);

    if (hasLoadedField) {
      specfac = loadedField[h] * deviceLib::exp(specphase);
    }
    deviceComplex ne, no;
    sellmeierDevice(&ne, &no, sellmeierCoefficients, f, (*s).crystalTheta,
                    (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);

    if (isComplexNaN(ne) || isComplexNaN(no)) {
      field[i] = deviceComplex{};
      field[i + (*s).NgridC] = deviceComplex{};
      return;
    }

    const deviceFP ko = twoPi<deviceFP>() * no.real() * f / lightC<deviceFP>();
    const deviceFP zR = vPi<deviceFP>() * (*p).beamwaist * (*p).beamwaist *
                        no.real() * f / lightC<deviceFP>();

    const deviceFP rB =
        ((*p).x0 - (*s).dx * (-0.5f * (*s).Nspace + (deviceFP)j) -
         0.25f * (*s).dx);
    const deviceFP r = rB * deviceFPLib::cos((*p).beamAngle) -
                       (*p).z0 * deviceFPLib::sin((*p).beamAngle);
    const deviceFP z = rB * deviceFPLib::sin((*p).beamAngle) +
                       (*p).z0 * deviceFPLib::cos((*p).beamAngle);

    const deviceFP wz =
        (*p).beamwaist * deviceFPLib::sqrt(1.0f + (z * z / (zR * zR)));
    const deviceFP Rz =
        (z != 0.0f) ? z * (1.0f + (zR * zR / (z * z))) : 1.0e15f;
    const deviceFP phi = deviceFPLib::atan(z / zR);
    deviceComplex Eb =
        deviceComplex(0.0f, 1.0f) *
            (ko * (z - (*p).z0) + ko * r * r / (2.0f * Rz) - phi) -
        r * r / (wz * wz);
    Eb = isComplexNaN(Eb) ? deviceComplex{}
                          : ((*p).beamwaist / wz) * deviceLib::exp(Eb);
    Eb = Eb * specfac;

    field[i] = deviceComplex(deviceFPLib::cos((*p).polarizationAngle),
                             -(*p).circularity *
                                 deviceFPLib::sin((*p).polarizationAngle)) *
               Eb;
    field[i + (*s).NgridC] =
        deviceComplex(deviceFPLib::sin((*p).polarizationAngle),
                      (*p).circularity *
                          deviceFPLib::cos((*p).polarizationAngle)) *
        Eb;
    if (isComplexNaN(field[i]))
      field[i] = deviceComplex{};
    if (isComplexNaN(field[i + (*s).NgridC]))
      field[i + (*s).NgridC] = deviceComplex{};
    deviceFP pointEnergy =
        deviceFPLib::abs(r) *
        (modulusSquared(field[i]) + modulusSquared(field[i + (*s).NgridC]));
    pointEnergy *= 2.0f * vPi<deviceFP>() * lightC<deviceFP>() *
                   eps0<deviceFP>() * (*s).dx * (*s).dt;
    // two factors of two cancel here - there should be one for the missing
    // frequency plane, but the sum is over x instead of r
    // accordingly we already multiplied by two
    atomicAdd(pulseSum, pointEnergy);
  }
};

// Generate a beam in full 3D mode
class beamGenerationKernel3D {
public:
  deviceComplex *field;
  const Pulse<deviceFP> *p;
  deviceFP *pulseSum;
  deviceParameterSet<deviceFP, deviceComplex> *s;
  const bool hasLoadedField;
  const deviceComplex *loadedField;
  const deviceFP *materialPhase;
  const deviceFP *sellmeierCoefficients;
  deviceFunction void operator()(int64_t i) const {
    const int64_t h = 1 + i % ((*s).Nfreq - 1);
    const int64_t col = i / ((*s).Nfreq - 1);
    i = h + col * ((*s).Nfreq);
    const int64_t j = col % (*s).Nspace;
    const int64_t k = col / (*s).Nspace;
    const deviceFP f = h * (*s).fStep;
    const deviceFP w = twoPi<deviceFP>() * (f - (*p).frequency);
    if (isComplexNaN(materialPhase[h])) {
      field[i] = deviceComplex{};
      field[i + (*s).NgridC] = deviceComplex{};
      return;
    }
    // supergaussian pulse spectrum, if no input pulse specified
    deviceComplex specfac = deviceComplex(
        -deviceFPLib::pow((f - (*p).frequency) / (*p).bandwidth, (*p).sgOrder),
        0.0f);

    const deviceComplex specphase =
        deviceComplex(0.0f, -((*p).cep +
                              twoPi<deviceFP>() * f *
                                  ((*p).delay - 0.5f * (*s).dt * (*s).Ntime) +
                              0.5f * (*p).gdd * w * w +
                              (*p).tod * w * w * w / 6.0f + materialPhase[h]));
    specfac = isComplexNaN(specphase) ? deviceComplex{}
                                      : deviceLib::exp(specfac + specphase);

    if (hasLoadedField) {
      specfac = loadedField[h] * deviceLib::exp(specphase);
    }
    deviceComplex ne, no;
    sellmeierDevice(&ne, &no, sellmeierCoefficients, f, (*s).crystalTheta,
                    (*s).crystalPhi, (*s).axesNumber, (*s).sellmeierType);
    if (isComplexNaN(ne) || isComplexNaN(no)) {
      field[i] = deviceComplex{};
      field[i + (*s).NgridC] = deviceComplex{};
      return;
    }
    const deviceFP ko = twoPi<deviceFP>() * no.real() * f / lightC<deviceFP>();
    const deviceFP zR = vPi<deviceFP>() * (*p).beamwaist * (*p).beamwaist *
                        no.real() * f / lightC<deviceFP>();

    const deviceFP xo =
        ((*s).dx * ((deviceFP)j - (*s).Nspace / 2.0f)) - (*p).x0;
    const deviceFP yo =
        ((*s).dx * ((deviceFP)k - (*s).Nspace2 / 2.0f)) - (*p).y0;
    const deviceFP zo = (*p).z0;
    const deviceFP cB = deviceFPLib::cos((*p).beamAngle);
    const deviceFP cA = deviceFPLib::cos((*p).beamAnglePhi);
    const deviceFP sB = deviceFPLib::sin((*p).beamAngle);
    const deviceFP sA = deviceFPLib::sin((*p).beamAnglePhi);
    const deviceFP x = cB * xo + sA * sB * yo + sA * sB * zo;
    const deviceFP y = cA * yo - sA * zo;
    const deviceFP z = -sB * xo + sA * cB * yo + cA * cB * zo;
    const deviceFP r = deviceFPLib::hypot(x, y);

    const deviceFP wz =
        (*p).beamwaist * deviceFPLib::sqrt(1.0f + (z * z / (zR * zR)));
    const deviceFP Rz =
        (z != 0.0f) ? z * (1.0f + (zR * zR / (z * z))) : 1.0e15f;
    const deviceFP phi = deviceFPLib::atan(z / zR);

    deviceComplex Eb = ((*p).beamwaist / wz) *
                       deviceLib::exp(deviceComplex(0.0f, 1.0f) *
                                          (ko * (z - (*p).z0) +
                                           ko * r * r / (2.0f * Rz) - phi) -
                                      r * r / (wz * wz));
    Eb = Eb * specfac;
    if (isComplexNaN(Eb) || f <= 0.0f) {
      Eb = deviceComplex{};
    }

    field[i] = deviceComplex(deviceFPLib::cos((*p).polarizationAngle),
                             -(*p).circularity *
                                 deviceFPLib::sin((*p).polarizationAngle)) *
               Eb;
    field[i + (*s).NgridC] =
        deviceComplex(deviceFPLib::sin((*p).polarizationAngle),
                      (*p).circularity *
                          deviceFPLib::cos((*p).polarizationAngle)) *
        Eb;
    deviceFP pointEnergy =
        (modulusSquared(field[i]) + modulusSquared(field[i + (*s).NgridC]));
    pointEnergy *= constProd(lightC<deviceFP>() * eps0<deviceFP>(), 2) *
                   (*s).dx * (*s).dx * (*s).dt;

    // factor 2 accounts for the missing negative frequency plane
    atomicAdd(pulseSum, pointEnergy);
  }
};

class multiplyByConstantKernelD {
public:
  deviceFP *A;
  const deviceFP val;
  deviceFunction void operator()(const int64_t i) const { A[i] *= val; }
};

class multiplicationKernelCompactDoubleVector {
public:
  const deviceFP *A;
  const deviceComplex *B;
  deviceComplex *C;
  const deviceParameterSet<deviceFP, deviceComplex> *s;
  deviceFunction void operator()(const int64_t i) const {
    const int64_t h = i % (*s).Nfreq; // temporal coordinate
    C[i] = A[h] * B[i];
  }
};

// Expand the information contained in the radially-symmetric beam in the offset
// grid
//  representation.
//  The grid is offset from the origin; rather than ...-2 -1 0 1 2... etc, which
//  would contain redundant information (the symmetry means that -1 and -1 are
//  equivalent) the grid is at the points -1.75 -0.75 0.25 1.25 2.25, etc. the
//  grid spacing is the same, but now the two sides of the origin contain
//  different information. This has effectively doubled the resolution of the
//  nonlinear polarization. We make use of this by expanding into the
//  full-resolution beam on the grid -2.25 -1.75 -1.25 -0.75 -0.25 0.25
//  0.75 1.25 1.75 2.25... after FFT, we can discard the high frequencies. Thus
//  we have downsampled in such a way as to avoid aliasing, which inside the
//  simulation is most likely the appear (and cause instability) in the
//  nonlinear terms.
class expandCylindricalBeam {
public:
  const deviceParameterSet<deviceFP, deviceComplex> *s;

  deviceFunction void operator()(const int64_t i) const {
    const int64_t j = i / (*s).Ntime; // spatial coordinate
    const int64_t k = i % (*s).Ntime; // temporal coordinate

    // positions on the expanded grid corresponding the the current index
    const int64_t pos1 = 2 * ((*s).Nspace - j - 1) * (*s).Ntime + k;
    const int64_t pos2 = (2 * j + 1) * (*s).Ntime + k;

    // reuse memory allocated for the radial Laplacian, casting complex double
    // to a 2x larger double real grid
    deviceFP *expandedBeam1 = (*s).gridRadialLaplacian1;
    deviceFP *expandedBeam2 = expandedBeam1 + 2 * (*s).Ngrid;

    expandedBeam1[pos1] = (*s).gridPolarizationTime1[i];
    expandedBeam1[pos2] = (*s).gridPolarizationTime1[i];
    expandedBeam2[pos1] = (*s).gridPolarizationTime2[i];
    expandedBeam2[pos2] = (*s).gridPolarizationTime2[i];
  }
};
} // namespace kernelNamespace
