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

    //sequence
    bool isInSequence;
    bool isFollowerInSequence;
    double* sequenceArray;
    int Nsequence;

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

struct guiStruct {
    HWND mainWindow;

    HWND tbMaterialIndex;
    HWND tbCrystalTheta;
    HWND tbCrystalPhi;
    HWND tbCrystalThickness;

    HWND tbGridXdim;
    HWND tbRadialStepSize;
    HWND tbTimeStepSize;

    HWND tbTimeSpan;
    HWND tbXstep;
    HWND tbBatchMode;
    HWND tbNumberSims;
    HWND tbBatchDestination;

    HWND tbPulse1Delay;
    HWND tbPulse2Delay;
    HWND tbFieldStrength1;
    HWND tbFieldStrength2;
    HWND tbBandwidth1;
    HWND tbBandwidth2;
    HWND tbFrequency1;
    HWND tbFrequency2;
    HWND tbCEPhase1;
    HWND tbCEPhase2;
    HWND tbPulseType;


    HWND tbGDD1;
    HWND tbGDD2;
    HWND tbTOD1;
    HWND tbTOD2;

    HWND tbXoffset1;
    HWND tbXoffset2;
    HWND tbZoffset1;
    HWND tbZoffset2;
    HWND tbBeamwaist1;
    HWND tbBeamwaist2;
    HWND tbPropagationAngle1;
    HWND tbPropagationAngle2;
    HWND tbPolarizationAngle1;
    HWND tbPolarizationAngle2;
    HWND tbCircularity1;
    HWND tbCircularity2;


    HWND pdPropagationMode;
    HWND pdBatchMode;
    HWND pdPulseType;
    HWND pdRStep;

    HWND tbSequence;

    HWND cbSavePsi;
    HWND tbFileNameBase;
    HWND tbPulse1Path;
    HWND pdPulse1Type;
    HWND tbPulse2Path;
    HWND pdPulse2Type;
    HWND buttonPulse1Path;
    HWND buttonPulse2Path;
    HWND buttonRun;
    HWND buttonFile;
    HWND buttonPlot;
    HWND buttonRefreshDB;
    HWND buttonStop;
    HWND tbWhichSimToPlot;
    HWND textboxSims;
};

double              vmaxa(double* v, int vlength);
DWORD WINAPI        mainsimthread(LPVOID lpParam);
int                 DrawLabels(HDC hdc);
int                 LabelTextBox(HDC hdc, HWND parentWindow, HWND targetTextBox, const wchar_t* labelText, int xOffset, int yOffset);
int                 HWNDToString(HWND inputA, char* outputString, int bufferSize);
double              HWNDToDouble(HWND inputA);
int                 AppendTextToWindow(HWND inputA, wchar_t* messagebuffer, int buffersize);
int                 getFileNameBaseFromDlg(HWND hWnd, HWND outputTextbox);
int                 getFileNameBaseFromDlgDat(HWND hWnd, HWND outputTextbox);
int                 DrawArrayAsBitmap(HDC hdc, INT64 Nx, INT64 Ny, INT64 x, INT64 y, INT64 height, INT64 width, double* data, int cm);
double              cmodulussquared(std::complex<double>complexNumber);
int                 drawsimplots(int simIndex);
int                 linearremap(double* A, int nax, int nay, double* B, int nbx, int nby, int modeInterp);
int                 readcrystaldatabase(struct crystalentry* db, bool isVerbose);
int                 resolvesequence(int currentIndex, struct propthread* s);