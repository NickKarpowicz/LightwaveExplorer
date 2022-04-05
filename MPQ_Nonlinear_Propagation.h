#pragma once

#include "resource.h"
#include<complex>


struct guiStruct {
    HWND mainWindow;

    HWND tbMaterialIndex;
    HWND tbCrystalTheta;
    HWND tbCrystalPhi;
    HWND tbCrystalThickness;
    HWND tbNonlinearAbsortion;
    HWND tbDrudeGamma;
    HWND tbBandGap;
    HWND tbEffectiveMass;

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
    HWND tbPulseEnergy1;
    HWND tbPulseEnergy2;
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
    HWND buttonRunOnCluster;
    HWND tbWhichSimToPlot;
    HWND textboxSims;
};

DWORD WINAPI        mainSimThread(LPVOID lpParam);
DWORD WINAPI        runOnCluster(LPVOID lpParam);
int                 drawLabels(HDC hdc);
int                 labelTextBox(HDC hdc, HWND parentWindow, HWND targetTextBox, const wchar_t* labelText, int xOffset, int yOffset);
int                 getStringFromHWND(HWND inputA, char* outputString, int bufferSize);
double              getDoubleFromHWND(HWND inputA);
int                 appendTextToWindow(HWND inputA, wchar_t* messageString, int buffersize);
int                 getFileNameBaseFromDlg(HWND hWnd, HWND outputTextbox);
int                 getFileNameBaseFromDlgDat(HWND hWnd, HWND outputTextbox);
int                 drawArrayAsBitmap(HDC hdc, INT64 Nx, INT64 Ny, INT64 x, INT64 y, INT64 height, INT64 width, double* data, int cm);
DWORD WINAPI        drawSimPlots(LPVOID lpParam);
int                 linearRemap(double* A, int nax, int nay, double* B, int nbx, int nby, int modeInterp);
int                 drawLabeledXYPlot(HDC hdc, int N, double* Y, double xStep, int posX, int posY, int pixelsWide, int pixelsTall, int forceYOrigin, double YOrigin, double yDiv);
int                 readParametersFromInterface();
int                 freeSemipermanentGrids();
template<typename... Args> void printToConsole(HWND console, const wchar_t* format, Args... args);
