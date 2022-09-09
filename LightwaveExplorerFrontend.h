#include "resource.h"
#include<complex>
#include<d2d1.h>
typedef struct guiStruct {
    HWND mainWindow = NULL;
    
    HWND tbMaterialIndex = NULL;
    HWND tbMaterialIndexAlternate = NULL;
    HWND tbCrystalTheta = NULL;
    HWND tbCrystalPhi = NULL;
    HWND tbCrystalThickness = NULL;
    HWND tbNonlinearAbsortion = NULL;
    HWND tbDrudeGamma = NULL;
    HWND tbBandGap = NULL;
    HWND tbEffectiveMass = NULL;

    HWND tbGridXdim = NULL;
    HWND tbRadialStepSize = NULL;
    HWND tbTimeStepSize = NULL;

    HWND tbTimeSpan = NULL;
    HWND tbXstep = NULL;
    HWND tbBatchMode = NULL;
    HWND tbNumberSims = NULL;
    HWND tbNumberSims2 = NULL;
    HWND tbBatchDestination = NULL;
    HWND tbBatchDestination2 = NULL;

    HWND tbPulse1Delay = NULL;
    HWND tbPulse2Delay = NULL;
    HWND tbPulseEnergy1 = NULL;
    HWND tbPulseEnergy2 = NULL;
    HWND tbBandwidth1 = NULL;
    HWND tbBandwidth2 = NULL;
    HWND tbFrequency1 = NULL;
    HWND tbFrequency2 = NULL;
    HWND tbCEPhase1 = NULL;
    HWND tbCEPhase2 = NULL;
    HWND tbPulseType1 = NULL;
    HWND tbPulseType2 = NULL;


    HWND tbGDD1 = NULL;
    HWND tbGDD2 = NULL;
    HWND tbTOD1 = NULL;
    HWND tbTOD2 = NULL;
    HWND tbPhaseMaterialIndex1 = NULL;
    HWND tbPhaseMaterialIndex2 = NULL;
    HWND tbPhaseMaterialThickness1 = NULL;
    HWND tbPhaseMaterialThickness2 = NULL;
    HWND tbXoffset1 = NULL;
    HWND tbXoffset2 = NULL;
    HWND tbZoffset1 = NULL;
    HWND tbZoffset2 = NULL;
    HWND tbBeamwaist1 = NULL;
    HWND tbBeamwaist2 = NULL;
    HWND tbPropagationAngle1 = NULL;
    HWND tbPropagationAngle2 = NULL;
    HWND tbPolarizationAngle1 = NULL;
    HWND tbPolarizationAngle2 = NULL;
    HWND tbCircularity1 = NULL;
    HWND tbCircularity2 = NULL;


    HWND pdPropagationMode = NULL;
    HWND pdBatchMode = NULL;
    HWND pdBatchMode2 = NULL;
    HWND pdRStep = NULL;

    HWND tbSequence = NULL;
    HWND tbFitting = NULL;

    HWND pbProgress = NULL;
    HWND pbProgressB = NULL;
    HWND cbSavePsi = NULL;
    HWND tbFileNameBase = NULL;
    HWND tbPulse1Path = NULL;
    HWND pdPulse1Type = NULL;
    HWND tbPulse2Path = NULL;
    HWND pdPulse2Type = NULL;
    HWND buttonPulse1Path = NULL;
    HWND buttonPulse2Path = NULL;
    HWND tbFittingReferencePath = NULL;
    HWND pdFittingType = NULL;
    HWND buttonFittingReference = NULL;
    HWND buttonRun = NULL;
    HWND buttonFile = NULL;
    HWND buttonPlot = NULL;
    HWND buttonRefreshDB = NULL;
    HWND buttonStop = NULL;
    HWND buttonLoad = NULL;
    HWND buttonAddCrystalToSequence = NULL;
    HWND buttonAddEchoSequence = NULL;
    HWND tbPlotNumber = NULL;
    HWND cbLogPlot = NULL;
    HWND pdClusterSelector = NULL;
    HWND buttonRunOnCluster = NULL;
    HWND tbWhichSimToPlot = NULL;
    HWND buttonFit = NULL;
    HWND textboxSims = NULL;
    HWND tbGPUStatus = NULL;
    HWND plotBox1 = NULL;
    HWND plotBox2 = NULL;
    HWND plotBox3 = NULL;
    HWND plotBox4 = NULL;
    HWND plotBox5 = NULL;
    HWND plotBox6 = NULL;
    HWND plotBox7 = NULL;
    HWND plotBox8 = NULL;
    ID2D1Factory* pFactory = NULL;
    int xOffsetRow1 = 160;
    int xOffsetRow2 = 480;
    int xOffsetRow3 = 640;
    int vs = 26;
    int radioButtonOffset = 349;
    int btnwidth = 100;
    int btnHeight = 26;
    int comboBoxHeight = 18;
    int btnoffset = 160;
    int btnoffset0 = 5;
    int btnoffset2 = 530;
    int btnoffset2a = 160 + 150 + 5;
    int rbsize = 18;
    int consoleSize = 630;
    int textboxwidth = 150;
    int plotSpacerX = 64;
    int plotSpacerY = 40;

} guiStruct;

typedef struct imagePlotStruct {
    HWND plotBox = NULL;
    double* data = NULL;
    std::complex<double>* complexData = NULL;
    int colorMap = 4;
    bool logScale = FALSE;
    double logMin;
    int dataType;
} imagePlotStruct;

typedef struct plotStruct {
    HWND plotBox = NULL;
    double* data = NULL;
    std::complex<double>* complexData = NULL;
    bool logScale = FALSE;
    double logMin = 0;
    int dataType = 0;
    double dx = 1.0;
    float x0 = 0.0;
    size_t Npts = 0;
    float unitY = 1.0;
    bool forceYmin = FALSE;
    float forcedYmin = 0.0;
    D2D1_COLOR_F color = D2D1::ColorF(1, 1, 1, 1);
} plotStruct;


DWORD WINAPI        mainSimThread(LPVOID lpParam);
DWORD WINAPI        createRunFile(LPVOID lpParam);
int                 drawLabels(HDC hdc);
int                 labelTextBox(HDC hdc, HWND parentWindow, HWND targetTextBox, const wchar_t* labelText, int xOffset, int yOffset);
int                 getStringFromHWND(HWND inputA, char* outputString, int bufferSize);
double              getDoubleFromHWND(HWND inputA);
int                 appendTextToWindow(HWND inputA, wchar_t* messageString, int buffersize);
int                 getFileNameBaseFromDlg(HWND hWnd, HWND outputTextbox);
int                 getFileNameBaseFromDlgDat(HWND hWnd, HWND outputTextbox);
int                 drawArrayAsBitmap(HWND plotBox, INT64 Nx, INT64 Ny, INT64 x, INT64 y, INT64 height, INT64 width, float* data, int cm);
DWORD WINAPI        drawSimPlots(LPVOID lpParam);
int                 linearRemap(float* A, int nax, int nay, float* B, int nbx, int nby);
int                 readParametersFromInterface();
int                 freeSemipermanentGrids();
template<typename... Args> void printToConsole(HWND console, const wchar_t* format, Args... args);
int                 floatyText(HDC hdc, HWND parentWindow, const wchar_t* labelText, int xOffset, int yOffset);
int                 openDialogBoxAndLoad(HWND hWnd);
void                setTitleBarDark(HWND hWnd);
DWORD WINAPI        plotXYDirect2d(LPVOID inputStruct);
int                 setInterfaceValuesToActiveValues();
DWORD WINAPI        fittingThread(LPVOID lpParam);
int                 insertLineBreaksAfterSemicolons(char* cString, size_t N);
DWORD WINAPI        statusMonitorThread(LPVOID lpParam);
int                 setWindowTextToDoubleExp(HWND win, double in);
int                 setWindowTextToDouble(HWND win, double in);
std::complex<double> getDoubleDoublesfromHWND(HWND inputA);
int                 linearRemapDoubleToFloat(double* A, int nax, int nay, float* B, int nbx, int nby);
int                 linearRemapZToLogFloat(std::complex<double>* A, int nax, int nay, float* B, int nbx, int nby, double logMin);
