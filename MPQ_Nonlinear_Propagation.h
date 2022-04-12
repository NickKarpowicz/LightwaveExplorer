#pragma once

#include "resource.h"
#include<complex>
#include<d2d1.h>
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
    HWND buttonLoad;
    HWND tbPlotNumber;
    HWND pdClusterSelector;
    HWND buttonRunOnCluster;
    HWND tbWhichSimToPlot;
    HWND textboxSims;
    HWND plotBox1;
    HWND plotBox2;
    HWND plotBox3;
    HWND plotBox4;
    HWND plotBox5;
    HWND plotBox6;
    HWND plotBox7;
    HWND plotBox8;
    ID2D1Factory* pFactory = NULL;
    int xOffsetRow1 = 160;
    int xOffsetRow2 = 480;
    int xOffsetRow3 = 640;
    int vs = 26;
    int radioButtonOffset = 349;
    int btnwidth = 120;
    int btnoffset = 160;
    int btnoffset0 = 5;
    int btnoffset2 = 510;
    int btnoffset2a = 160 + 150 + 5;
    int rbsize = 18;
    int consoleSize = 630;
    int textboxwidth = 150;
};


DWORD WINAPI        mainSimThread(LPVOID lpParam);
DWORD WINAPI        createRunFile(LPVOID lpParam);
int                 drawLabels(HDC hdc);
int                 labelTextBox(HDC hdc, HWND parentWindow, HWND targetTextBox, const wchar_t* labelText, int xOffset, int yOffset);
int                 getStringFromHWND(HWND inputA, char* outputString, int bufferSize);
double              getDoubleFromHWND(HWND inputA);
int                 appendTextToWindow(HWND inputA, wchar_t* messageString, int buffersize);
int                 getFileNameBaseFromDlg(HWND hWnd, HWND outputTextbox);
int                 getFileNameBaseFromDlgDat(HWND hWnd, HWND outputTextbox);
int                 drawArrayAsBitmap(HDC hdc, INT64 Nx, INT64 Ny, INT64 x, INT64 y, INT64 height, INT64 width, float* data, int cm);
DWORD WINAPI        drawSimPlots(LPVOID lpParam);
int                 linearRemap(float* A, int nax, int nay, float* B, int nbx, int nby);
int                 readParametersFromInterface();
int                 freeSemipermanentGrids();
template<typename... Args> void printToConsole(HWND console, const wchar_t* format, Args... args);
int                 floatyText(HDC hdc, HWND parentWindow, const wchar_t* labelText, int xOffset, int yOffset);
int                 openDialogBoxAndLoad(HWND hWnd);
void                setTitleBarDark(HWND hWnd);
void                plotXYDirect2d(HWND targetWindow, float dX, float* Y, size_t Npts, float unitY, bool forceminY, float forcedminY);
