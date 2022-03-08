#pragma once

#include "resource.h"
#include<complex>
//Simulation parameter struct to pass to the simulations running in threads
struct propthread {
    double rStep;
    double tStep;
    double fStep;
    double kStep;
    double propagationStep;
    long long Npropagation;
    long long Ntime;
    long long Nspace;
    long long Ngrid;
    int NSims;

    int materialIndex;
    double crystalTheta;
    double crystalPhi;
    double crystalThickness;
    double* chi2Tensor;
    double* deffTensor;
    double* chi3Tensor;
    double* sellmeierCoefficients;
    std::complex<double>* refractiveIndex1;
    std::complex<double>* refractiveIndex2;
    double* absorptionParameters;
    int sellmeierType;
    int axesNumber;
    double neref;
    double noref;
    int* nonlinearSwitches;

    //spatial beam properties;
    double beamwaist1;
    double beamwaist2;
    double z01;
    double z02;
    double x01;
    double x02;
    double propagationAngle1;
    double propagationAngle2;
    int isCylindric;

    //spectral/temporal field properties
    double pulseEnergy1;
    double pulseEnergy2;
    double frequency1;
    double frequency2;
    int sgOrder1;
    int sgOrder2;
    double bandwidth1;
    double bandwidth2;
    double cephase1;
    double cephase2;
    double delay1;
    double delay2;
    double gdd1;
    double gdd2;
    double tod1;
    double tod2;

    //loaded FROG/EOS fields
    std::complex<double>* loadedField1;
    std::complex<double>* loadedField2;
    bool field1IsAllocated = 0;
    bool field2IsAllocated = 0;



    //polarization properties
    double polarizationAngle1;
    double polarizationAngle2;
    double circularity1;
    double circularity2;


    int pulsetype;
    double* InfoVec;
    std::complex<double>* Ext;
    std::complex<double>* Ekw;
    std::complex<double>* ExtOut;
    std::complex<double>* EkwOut;
    int* imdone;
    
    struct crystalentry* crystalDatabase;

};

struct crystalentry {
    wchar_t crystalNameW[256];
    int axisType;
    int sellmeierType;
    int nonlinearSwitches[4];
    double sellmeierCoefficients[66];
    wchar_t sellmeierReference[512];
    double d[18];
    wchar_t dReference[512];
    double chi3[48];
    wchar_t chi3Reference[512];
    double absorptionParameters[6];
    wchar_t spectralFile[512];
    double spectralData[2048];
};

double              cmodulussquared(std::complex<double>complexNumber);