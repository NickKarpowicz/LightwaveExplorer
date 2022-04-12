#define _CRT_SECURE_NO_DEPRECATE
#include "framework.h"
#include "MPQ_Nonlinear_Propagation.h"
#include "NonlinearPropCUDA.cuh"
#include<cstdio>
#include<complex>
#include<math.h>
#include<Commdlg.h>
#include<stdlib.h>
#include<chrono>
#include<Windows.h>
#include<Uxtheme.h>
#include<dwmapi.h>
#include<d2d1.h>

#define MAX_LOADSTRING 1024
#define ID_BTNRUN 11110
#define ID_BTNPLOT 11111
#define ID_BTNGETFILENAME 11112
#define ID_BTNREFRESHDB 11113
#define ID_BTNSTOP 11114
#define ID_BTNPULSE1 11115
#define ID_BTNPULSE2 11116
#define ID_BTNRUNONCLUSTER 11117
#define ID_BTNLOAD 11118

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name
struct guiStruct maingui;                       // struct containing all of the GUI elements
struct simulationParameterSet* activeSetPtr;    // Main structure containing simulation parameters and pointers
struct crystalEntry* crystalDatabasePtr;        // Crystal info database
WCHAR programDirectory[MAX_LOADSTRING];         // Program working directory (useful if the crystal database has to be reloaded)
bool isRunning = FALSE;
bool isPlotting = FALSE;
bool isGridAllocated = FALSE;
bool cancellationCalled = FALSE;

COLORREF uiWhite = RGB(255, 255, 255);
COLORREF uiGrey = RGB(216, 216, 216);
COLORREF uiDarkGrey = RGB(32, 32, 32);
COLORREF uiBlack = RGB(16, 16, 16);
HBRUSH blackBrush = CreateSolidBrush(uiBlack);

// Forward declarations of (Microsoft) functions included in this code module:
ATOM                MyRegisterClass(HINSTANCE hInstance);
bool                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

DWORD WINAPI mainSimThread(LPVOID lpParam) {
    cancellationCalled = FALSE;
    int j;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();
    HANDLE plotThread;
    DWORD hplotThread;
    readParametersFromInterface();
    (*activeSetPtr).runType = 0;
    if (isGridAllocated) {
        freeSemipermanentGrids();
    }

    allocateGrids(activeSetPtr);
    isGridAllocated = TRUE;
    (*activeSetPtr).crystalDatabase = crystalDatabasePtr;
    loadPulseFiles(activeSetPtr);
    readSequenceString(activeSetPtr);
    configureBatchMode(activeSetPtr);

    //run the simulations
    for (j = 0; j < (*activeSetPtr).Nsims; j++) {
        if ((*activeSetPtr).isInSequence) {
            solveNonlinearWaveEquationSequence(&activeSetPtr[j]);
            (*activeSetPtr).plotSim = j;
            plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);
        }
        else {
            solveNonlinearWaveEquation(&activeSetPtr[j]);
            //printToConsole(maingui.textboxSims, _T("Sellmeier check: f: %lf n1: %lf n2: %lf.\r\n"), 1e-12 * activeSetPtr[j].fStep * 64, real(activeSetPtr[j].refractiveIndex1[64]), real(activeSetPtr[j].refractiveIndex2[64]));
            
            if (activeSetPtr[j].memoryError > 0) {
                printToConsole(maingui.textboxSims, _T("Warning: device memory error (%i).\r\n"), activeSetPtr[j].memoryError);
            }
        }

        if (cancellationCalled) {
            printToConsole(maingui.textboxSims, _T("Warning: series cancelled, stopping after %i simulations.\r\n"), j + 1);
            break;
        }

        (*activeSetPtr).plotSim = j;
        plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);

    }

    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    printToConsole(maingui.textboxSims, _T("Finished after %8.4lf s. \r\n"), 1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
    saveDataSet(activeSetPtr, crystalDatabasePtr, (*activeSetPtr).outputBasePath);

    free((*activeSetPtr).sequenceArray);
    free((*activeSetPtr).refractiveIndex1);
    free((*activeSetPtr).refractiveIndex2);
    free((*activeSetPtr).imdone);
    free((*activeSetPtr).deffTensor);
    free((*activeSetPtr).loadedField1);
    free((*activeSetPtr).loadedField2);
    
    isRunning = FALSE;
    return 0;
}


//Instead of running locally, make the files needed to run it on
//one of the GPU-equipped clusters. Specifically, Raven and Cobra
//of the MPCDF
DWORD WINAPI createRunFile(LPVOID lpParam) {

    readParametersFromInterface();

    if (isGridAllocated) {
        freeSemipermanentGrids();
    }

    allocateGrids(activeSetPtr);
    isGridAllocated = TRUE;
    (*activeSetPtr).crystalDatabase = crystalDatabasePtr;
    loadPulseFiles(activeSetPtr);
    readSequenceString(activeSetPtr);
    configureBatchMode(activeSetPtr);


    free((*activeSetPtr).sequenceArray);
    free((*activeSetPtr).refractiveIndex1);
    free((*activeSetPtr).refractiveIndex2);
    free((*activeSetPtr).imdone);
    free((*activeSetPtr).deffTensor);
    free((*activeSetPtr).loadedField1);
    free((*activeSetPtr).loadedField2);
    char* fileName = (*activeSetPtr).outputBasePath;
    wchar_t wideBuffer[MAX_LOADSTRING];
    while (strchr(fileName, '\\') != NULL) {
        fileName = strchr(fileName, '\\');
        fileName++;
    }
    
    mbstowcs(wideBuffer, fileName, strlen(fileName)+1);
    printToConsole(maingui.textboxSims, 
        L"Run on cluster with:\r\nsbatch %ls.slurmScript\r\n",
        wideBuffer);
    
    //create command line settings file
    saveSettingsFile(activeSetPtr, crystalDatabasePtr);

    //create SLURM script
    int gpuType = 0;
    int gpuCount = 1;
	switch ((*activeSetPtr).runType) {
	case 1:
		gpuType = 0;
		gpuCount = 1;
		break;
	case 2:
		gpuType = 0;
		gpuCount = 2;
		break;
	case 3:
		gpuType = 1;
		gpuCount = 1;
		break;
	case 4:
		gpuType = 1;
		gpuCount = 2;
		break;
	case 5:
		gpuType = 2;
		gpuCount = 1;
		break;
	case 6:
		gpuType = 2;
		gpuCount = 2;
		break;
	case 7:
		gpuType = 2;
		gpuCount = 4;
		break;

	}
    saveSlurmScript(activeSetPtr, gpuType, gpuCount);

    isRunning = FALSE;
    return 0;
}

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
    _In_opt_ HINSTANCE hPrevInstance,
    _In_ LPWSTR    lpCmdLine,
    _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // Initialize global strings
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_MPQNONLINEARPROPAGATION, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Perform application initialization:
    if (!InitInstance(hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_MPQNONLINEARPROPAGATION));

    MSG msg;

    // Main message loop:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        //this is literally the only change I made from the default template - dialog message handling
        //allows the tab key to switch between text fields
        if (!IsDialogMessage(maingui.mainWindow, &msg))
        {
            if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
            {
                TranslateMessage(&msg);
                DispatchMessage(&msg);
            }
        }

    }
    return (int)msg.wParam;
}



//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc = WndProc;
    wcex.cbClsExtra = 0;
    wcex.cbWndExtra = 0;
    wcex.hInstance = hInstance;
    wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_MPQNONLINEARPROPAGATION));
    wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground = CreateSolidBrush(uiDarkGrey);
    wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_MPQNONLINEARPROPAGATION);
    wcex.lpszClassName = szWindowClass;
    wcex.hIconSm = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

// InitInstance is called by Windows when the program starts, this creates the main application window and
// populates it with elements, and does initial preparations:
// - dynamic allocations of crystal database and main struct pointer
// - loading crystal database
// - identify GPU
bool InitInstance(HINSTANCE hInstance, int nCmdShow)
{
    int k = 0;
    
    hInst = hInstance; // Store instance handle in our global variable
    int xOffsetRow1 = maingui.xOffsetRow1;
    int xOffsetRow2 = maingui.xOffsetRow2;
    int xOffsetRow3 = maingui.xOffsetRow3;
    int vs = maingui.vs;
    int radioButtonOffset = maingui.radioButtonOffset;
    int btnwidth = maingui.btnwidth;
    int btnoffset = maingui.btnoffset;
    int btnoffset0 = maingui.btnoffset0;
    int btnoffset2 = maingui.btnoffset2;
    int btnoffset2a = maingui.btnoffset2a;
    int rbsize = maingui.rbsize;
    int consoleSize = maingui.consoleSize;
    int textboxwidth = maingui.textboxwidth;
    maingui.mainWindow = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW | WS_EX_CONTROLPARENT,
        CW_USEDEFAULT, CW_USEDEFAULT, 2200, 33 * vs + consoleSize + 60, nullptr, nullptr, hInstance, nullptr);
    SetMenu(maingui.mainWindow, NULL);
    SetWindowTextA(maingui.mainWindow, "Nick's nonlinear propagator");

    //text boxes for input parameters
    maingui.tbPulseEnergy1 = CreateWindow(WC_EDIT, TEXT("24e-9"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 0 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulseEnergy2 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 1 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFrequency1 = CreateWindow(WC_EDIT, TEXT("130"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 2 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFrequency2 = CreateWindow(WC_EDIT, TEXT("200"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 3 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBandwidth1 = CreateWindow(WC_EDIT, TEXT("40"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 4 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBandwidth2 = CreateWindow(WC_EDIT, TEXT("40"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 5 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulseType = CreateWindow(WC_EDIT, TEXT("2"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 6 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCEPhase1 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 7 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCEPhase2 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 8 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulse1Delay = CreateWindow(WC_EDIT, TEXT("-90"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 9 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulse2Delay = CreateWindowW(WC_EDIT, TEXT("-20"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 10 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbGDD1 = CreateWindow(WC_EDIT, TEXT("-27.262"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 11 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbGDD2 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 12 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTOD1 = CreateWindow(WC_EDIT, TEXT("230.57"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 13 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTOD2 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 14 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.tbBeamwaist1 = CreateWindow(WC_EDIT, TEXT("2.8"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 15 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBeamwaist2 = CreateWindow(WC_EDIT, TEXT("50"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 16 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbXoffset1 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 17 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbXoffset2 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 18 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbZoffset1 = CreateWindow(WC_EDIT, TEXT("-100"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 19 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbZoffset2 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 20 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPropagationAngle1 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 21 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPropagationAngle2 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 22 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPolarizationAngle1 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 23 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPolarizationAngle2 = CreateWindow(WC_EDIT, TEXT("90"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 24 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCircularity1 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 25 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCircularity2 = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 26 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFileNameBase = CreateWindow(WC_EDIT, TEXT("TestFile"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow3, 1 * vs, 775, 20, maingui.mainWindow, NULL, hInstance, NULL);
    
    maingui.plotBox1 = CreateWindow(WC_STATIC, NULL, WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox2 = CreateWindow(WC_STATIC, NULL, WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox3 = CreateWindow(WC_STATIC, NULL, WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox4 = CreateWindow(WC_STATIC, NULL, WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox5 = CreateWindow(WC_STATIC, NULL, WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox6 = CreateWindow(WC_STATIC, NULL, WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox7 = CreateWindow(WC_STATIC, NULL, WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox8 = CreateWindow(WC_STATIC, NULL, WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    D2D1CreateFactory(D2D1_FACTORY_TYPE_MULTI_THREADED, &maingui.pFactory);

    maingui.tbMaterialIndex = CreateWindow(WC_EDIT, TEXT("3"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 0 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCrystalTheta = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 1 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCrystalPhi = CreateWindow(WC_EDIT, TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 2 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.tbNonlinearAbsortion = CreateWindow(WC_EDIT, TEXT("2e-87"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 3 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBandGap = CreateWindow(WC_EDIT, TEXT("3"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 4 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbDrudeGamma = CreateWindow(WC_EDIT, TEXT("10"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 5 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbEffectiveMass = CreateWindow(WC_EDIT, TEXT("0.081"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 6 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.tbGridXdim = CreateWindow(WC_EDIT, TEXT("198.4"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 7 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbRadialStepSize = CreateWindow(WC_EDIT, TEXT("0.62"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 8 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTimeSpan = CreateWindow(WC_EDIT, TEXT("358.4"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 9 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTimeStepSize = CreateWindow(WC_EDIT, TEXT("0.7"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 10 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCrystalThickness = CreateWindow(WC_EDIT, TEXT("500"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 11 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbXstep = CreateWindow(WC_EDIT, TEXT("25"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 12 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    
    maingui.tbBatchDestination = CreateWindow(WC_EDIT, TEXT("20"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 16 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbNumberSims = CreateWindow(WC_EDIT, TEXT("1"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 17 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
  
    maingui.buttonRun = CreateWindowW(WC_BUTTON, TEXT("Run"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | BS_CENTER | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2, 22 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNRUN, hInstance, NULL);
    maingui.buttonStop = CreateWindow(WC_BUTTON, TEXT("Stop"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2, 23 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNSTOP, hInstance, NULL);
    maingui.buttonRefreshDB = CreateWindow(WC_BUTTON, TEXT("Refresh DB"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2, 24 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNREFRESHDB, hInstance, NULL);
    maingui.buttonRunOnCluster = CreateWindow(WC_BUTTON, TEXT("Cluster script"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2a, 25 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNRUNONCLUSTER, hInstance, NULL);
    maingui.pdClusterSelector = CreateWindow(WC_COMBOBOX, TEXT(""), CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, btnoffset2a + btnwidth + 5, 25 * vs -2, 190, 8 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    
    TCHAR A[64];
    memset(&A, 0, sizeof(A));
    TCHAR clusterNames[7][64] = { L"Cobra 1xRTX 5000", L"Cobra 2xRTX 5000", L"Cobra 1xV100", L"Cobra 2xV100", L"Raven A100", L"Raven 2xA100", L"Raven 4xA100"};
    for (k = 0; k < 7; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)clusterNames[k]);
        SendMessage(maingui.pdClusterSelector, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdClusterSelector, CB_SETCURSEL, (WPARAM)0, 0);
    maingui.buttonPlot = CreateWindow(WC_BUTTON, TEXT("Plot"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2a, 22 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNPLOT, hInstance, NULL);
    maingui.tbWhichSimToPlot = CreateWindow(WC_EDIT, TEXT("1"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2a + btnwidth + 5, 22 * vs, 60, 20, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.buttonLoad = CreateWindow(WC_BUTTON, TEXT("Load"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2a, 23 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNLOAD, hInstance, NULL);


    maingui.tbSequence = CreateWindow(WC_EDIT, TEXT(""), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL, xOffsetRow1 + textboxwidth + 4, 19 * vs-2, xOffsetRow2-xOffsetRow1, 66, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonFile = CreateWindow(WC_BUTTON, TEXT("Set Path"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow3, 0 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNGETFILENAME, hInstance, NULL);

    
    
    maingui.pdPropagationMode = CreateWindow(WC_COMBOBOX, TEXT(""), CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, xOffsetRow2, 13 * vs, textboxwidth, 5 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR energyModeNames[2][64] = {
        TEXT("2D Cartesian"), TEXT("3D radial symmetry")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 2; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)energyModeNames[k]);
        SendMessage(maingui.pdPropagationMode, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdPropagationMode, CB_SETCURSEL, (WPARAM)1, 0);

    maingui.pdBatchMode = CreateWindow(WC_COMBOBOX, TEXT(""), CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, xOffsetRow2, 15 * vs - 4, textboxwidth, 23 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR batchModeNames[11][64] = {
        TEXT("none"), TEXT("Pulse 2 delay"), TEXT("Pulse 1 Energy"), TEXT("CEP"), TEXT("Propagation"), TEXT("Theta"), TEXT("Pulse 1 GDD"), TEXT("Pulse 1 z"), TEXT("Gamma"), TEXT("NL absorption"), TEXT("Beamwaist 1")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 11; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)batchModeNames[k]);
        SendMessage(maingui.pdBatchMode, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdBatchMode, CB_SETCURSEL, (WPARAM)0, 0);

    maingui.pdPulse1Type = CreateWindow(WC_COMBOBOX, TEXT(""), CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, xOffsetRow1, 27 * vs - 4, textboxwidth, 9 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR pdPulse1Names[3][64] = {
        TEXT("Synthetic"), TEXT("load FROG"), TEXT("load EOS")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 3; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)pdPulse1Names[k]);
        SendMessage(maingui.pdPulse1Type, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdPulse1Type, CB_SETCURSEL, (WPARAM)0, 0);
    maingui.tbPulse1Path = CreateWindow(WC_EDIT, TEXT("pulse1.speck.dat"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL, 0, 28 * vs, xOffsetRow2 + 150, 46, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonPulse1Path = CreateWindow(WC_BUTTON, TEXT("Set path 1"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2a, 27 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNPULSE1, hInstance, NULL);

    maingui.pdPulse2Type = CreateWindow(WC_COMBOBOX, TEXT(""), CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, xOffsetRow1, 30 * vs - 4, textboxwidth, 9 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR pdPulse2Names[3][64] = {
        TEXT("Synthetic"), TEXT("load FROG"), TEXT("load EOS")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 3; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)pdPulse2Names[k]);
        SendMessage(maingui.pdPulse2Type, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdPulse2Type, CB_SETCURSEL, (WPARAM)0, 0);
    maingui.tbPulse2Path = CreateWindow(WC_EDIT, TEXT("pulse2.speck.dat"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL, 0, 31 * vs, xOffsetRow2 + 150, 46, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonPulse2Path = CreateWindow(WC_BUTTON, TEXT("Set path 2"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1 + textboxwidth + 5, 30 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNPULSE2, hInstance, NULL);


    //Text message window
    maingui.textboxSims = CreateWindow(WC_EDIT, TEXT(""), WS_CHILD | WS_VISIBLE | WS_BORDER | ES_LEFT | ES_MULTILINE | WS_VSCROLL | WS_HSCROLL, 0, 33 * vs, consoleSize, consoleSize, maingui.mainWindow, NULL, hInstance, NULL);

    if (!maingui.mainWindow)
    {
        return FALSE;
    }


    SetWindowTheme(maingui.mainWindow, L"DarkMode_Explorer", NULL);
    ShowWindow(maingui.mainWindow, nCmdShow);

    //make the active set pointer
    activeSetPtr = (struct simulationParameterSet*)calloc(2048, sizeof(struct simulationParameterSet));
    (*activeSetPtr).outputBasePath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    (*activeSetPtr).sequenceString = (char*)calloc(MAX_LOADSTRING * 256, sizeof(char));
    (*activeSetPtr).field1FilePath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    (*activeSetPtr).field2FilePath = (char*)calloc(MAX_LOADSTRING, sizeof(char));

    //Find, count, and name the GPUs
    int CUDAdevice, i;
    int CUDAdeviceCount = 0;
    cudaGetDeviceCount(&CUDAdeviceCount);    
    cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
    struct cudaDeviceProp activeCUDADeviceProp;
    wchar_t wcstring[514];
    size_t convertedChars = 0;
    if (cuErr == cudaSuccess) {
        printToConsole(maingui.textboxSims, _T("Found %i GPU(s): \r\n"), CUDAdeviceCount);
        for (i = 0; i < CUDAdeviceCount; i++) {
            cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
            memset(wcstring, 0, sizeof(wchar_t));
            mbstowcs_s(&convertedChars, wcstring, 256, activeCUDADeviceProp.name, _TRUNCATE);
            printToConsole(maingui.textboxSims, _T("%ls\r\n"), wcstring);
            printToConsole(maingui.textboxSims, _T(" Memory: %lli MB; Multiprocessors: %i\r\n"), activeCUDADeviceProp.totalGlobalMem/(1024*1024), activeCUDADeviceProp.multiProcessorCount);
        }

    }
    else {
        printToConsole(maingui.textboxSims, L"No compatible GPU found.\r\n");
    }
    
    //read the crystal database
    crystalDatabasePtr = (struct crystalEntry*)calloc(MAX_LOADSTRING, sizeof(struct crystalEntry));
    if (crystalDatabasePtr != NULL) {
        GetCurrentDirectory(MAX_LOADSTRING - 1, programDirectory);
        readCrystalDatabase(crystalDatabasePtr);
        printToConsole(maingui.textboxSims, _T("Read %i entries:\r\n"), (*crystalDatabasePtr).numberOfEntries);
        for (i = 0; i < (*crystalDatabasePtr).numberOfEntries; i++) {
            printToConsole(maingui.textboxSims, _T("Material %i name: %s\r\n"), i, crystalDatabasePtr[i].crystalNameW);
        }
    }


    return TRUE;
}

// This function handles button presses and other interaction with the application window
// to add a button, just give it an ID number and add it to the case structure
// then put whatever code should run when the button is pressed in that case
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    int plotSim;
    HANDLE mainthread;
    DWORD hMainThread;
    HANDLE plotThread;
    DWORD hplotThread;

    switch (message)
    {
    case WM_COMMAND:
    {
        int wmId = LOWORD(wParam);

        // Parse the menu selections:
        switch (wmId)
        {
        case IDM_ABOUT:
            DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
            break;
        case IDM_EXIT:
            if (isGridAllocated) {
                free((*activeSetPtr).ExtOut);
                free((*activeSetPtr).EkwOut);
                free((*activeSetPtr).Ext);
                free((*activeSetPtr).Ekw);
            }
            free((*activeSetPtr).outputBasePath);
            free((*activeSetPtr).field1FilePath);
            free((*activeSetPtr).field2FilePath);
            free(activeSetPtr);
            free(crystalDatabasePtr);
            DestroyWindow(hWnd);
            break;
        case ID_BTNRUN:
            if (!isRunning) {
                isRunning = TRUE;
                mainthread = CreateThread(NULL, 0, mainSimThread, activeSetPtr, 0, &hMainThread);
            }
            break;
        case ID_BTNRUNONCLUSTER:
            if (!isRunning) {
                isRunning = TRUE;
                mainthread = CreateThread(NULL, 0, createRunFile, activeSetPtr, 0, &hMainThread);
            }
            break;
        case ID_BTNSTOP:
            if (isRunning) {
                cancellationCalled = TRUE;
                for (int i = 0; i < (*activeSetPtr).Nsims; i++) {
                    (*activeSetPtr).imdone[i] = 2;
                }
            }
            break;
        case ID_BTNGETFILENAME:
            getFileNameBaseFromDlg(hWnd, maingui.tbFileNameBase);
            break;
        case ID_BTNPULSE1:
            getFileNameBaseFromDlgDat(hWnd, maingui.tbPulse1Path);
            break;
        case ID_BTNPULSE2:
            getFileNameBaseFromDlgDat(hWnd, maingui.tbPulse2Path);
            break;
        case ID_BTNREFRESHDB:
            SetCurrentDirectory(programDirectory);
            memset(crystalDatabasePtr, 0, 512 * sizeof(struct crystalEntry));
            readCrystalDatabase(crystalDatabasePtr);
            printToConsole(maingui.textboxSims, _T("Read %i entries:\r\n"), (*crystalDatabasePtr).numberOfEntries);
            for (int i = 0; i < (*crystalDatabasePtr).numberOfEntries; i++) {
                printToConsole(maingui.textboxSims, _T("Material %i name: %s\r\n"), i, crystalDatabasePtr[i].crystalNameW);
            }
            break;
        case ID_BTNPLOT:
            if (isGridAllocated && !isPlotting) {
                plotSim = (int)getDoubleFromHWND(maingui.tbWhichSimToPlot);
                plotSim = min(plotSim, (int)(*activeSetPtr).Nsims);
                plotSim--;
                plotSim = max(plotSim, 0);
                if (isRunning && (*activeSetPtr).imdone[plotSim] == 0) {
                    (*activeSetPtr).imdone[plotSim] = 3;
                    while ((*activeSetPtr).imdone[plotSim] == 3) {
                        Sleep(30);
                    }
                }

                (*activeSetPtr).plotSim = plotSim;
                plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);
            }
            break;
        case ID_BTNLOAD:
            openDialogBoxAndLoad(hWnd);
            
            plotSim = (int)getDoubleFromHWND(maingui.tbWhichSimToPlot);
            plotSim = min(plotSim, (int)(*activeSetPtr).Nsims);
            plotSim--;
            plotSim = max(plotSim, 0);

            (*activeSetPtr).plotSim = plotSim;
            plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);

            break;
        default:
            return DefWindowProc(hWnd, message, wParam, lParam);
        }
    }
    case WM_CTLCOLORBTN:
    {
        HDC hdc = (HDC)wParam;
        SetBkColor(hdc, uiBlack);
        SetTextColor(hdc, uiGrey);
        setTitleBarDark(maingui.buttonRun);
        return (LRESULT)blackBrush;
    }
    case WM_CTLCOLORSCROLLBAR:
    {
        HDC hdc = (HDC)wParam;
        //SetBkColor(hdc, uiBlack);
        return (LRESULT)blackBrush;
    }
    case WM_CTLCOLORLISTBOX:
    {
        HDC hdc = (HDC)wParam;
        SetBkColor(hdc, uiBlack);
        SetTextColor(hdc, uiGrey);
        return (LRESULT)blackBrush;
    }
    case WM_CTLCOLORSTATIC:
    {
        HDC hdc = (HDC)wParam;
        SetBkColor(hdc, uiBlack);
        SetTextColor(hdc, uiGrey);
        return (LRESULT)blackBrush;
    }
    case WM_CTLCOLOREDIT:
    {
        HDC hdc = (HDC)wParam;
        SetBkColor(hdc, uiBlack);
        SetTextColor(hdc, uiGrey);
        return (LRESULT)blackBrush;
    }
    case WM_PAINT:
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hWnd, &ps);
        SetTextColor(hdc, uiGrey);
        SetBkColor(hdc, uiDarkGrey);
        drawLabels(hdc);
        EndPaint(hWnd, &ps);
       
        break;
    }
    case WM_SIZE:
    {
        RECT mainRect;
        GetWindowRect(maingui.mainWindow, &mainRect);
        SetWindowPos(maingui.textboxSims, HWND_TOP, 0, 33*maingui.vs, maingui.consoleSize, mainRect.bottom - mainRect.top - 35*maingui.vs -4, NULL);
        SetWindowPos(maingui.tbFileNameBase, HWND_TOP, maingui.xOffsetRow3, maingui.vs, mainRect.right - mainRect.left - maingui.xOffsetRow3- 30, 20, NULL);

        int spacerX = 50;
        int spacerY = 40;
        int x0 = maingui.consoleSize + spacerX + 10;
        int y0 = 90 + spacerY;
        int imagePanelSizeX = mainRect.right - mainRect.left - x0 - 2 * spacerX;
        int imagePanelSizeY = mainRect.bottom - mainRect.top - y0 - 5 * spacerY;

        size_t simIndex = (*activeSetPtr).plotSim;
        int x = x0;
        int y = y0;
        int dx = imagePanelSizeX / 2;
        int dy = imagePanelSizeY / 4;
        int plotMargin = 0;
        int xCorrection = 10;
        int yCorrection = 45;
        if (!isGridAllocated) {
            SetWindowPos(maingui.plotBox1, HWND_BOTTOM, x - xCorrection, y - yCorrection, dx, dy, NULL);
            SetWindowPos(maingui.plotBox2, HWND_BOTTOM, x - xCorrection, y + 1 * dy + 1 * spacerY - yCorrection, dx, dy, NULL);
            SetWindowPos(maingui.plotBox3, HWND_BOTTOM, x - xCorrection, y + 2 * dy + 2 * spacerY - yCorrection, dx, dy, NULL);
            SetWindowPos(maingui.plotBox4, HWND_BOTTOM, x - xCorrection, y + 3 * dy + 3 * spacerY - yCorrection, dx, dy, NULL);
            SetWindowPos(maingui.plotBox5, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 0 * dy + 0 * spacerY - yCorrection, dx, dy, NULL);
            SetWindowPos(maingui.plotBox6, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 1 * dy + 1 * spacerY - yCorrection, dx, dy, NULL);
            SetWindowPos(maingui.plotBox7, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 2 * dy + 2 * spacerY - yCorrection, dx, dy, NULL);
            SetWindowPos(maingui.plotBox8, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 3 * dy + 3 * spacerY - yCorrection, dx, dy, NULL);
        }


        if (isGridAllocated && !isRunning) {
            drawSimPlots(activeSetPtr);
        }
        break;
    }
    case WM_THEMECHANGED:
		setTitleBarDark(hWnd);
		SetWindowTheme(maingui.textboxSims, L"DarkMode_Explorer", NULL);
		SetWindowTheme(maingui.tbPulse1Path, L"DarkMode_Explorer", NULL);
		SetWindowTheme(maingui.tbPulse2Path, L"DarkMode_Explorer", NULL);
		SetWindowTheme(maingui.tbSequence, L"DarkMode_Explorer", NULL);
		SetWindowTheme(maingui.tbFileNameBase, L"DarkMode_Explorer", NULL);
		SetWindowTheme(maingui.pdBatchMode, L"DarkMode_Explorer", NULL);
		UpdateWindow(hWnd);
		break;
    case WM_DESTROY:
        free(crystalDatabasePtr);
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}


// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
int freeSemipermanentGrids() {
    isGridAllocated = FALSE;
    free((*activeSetPtr).ExtOut);
    free((*activeSetPtr).EkwOut);
    free((*activeSetPtr).Ext);
    free((*activeSetPtr).Ekw);
    return 0;
}

int readParametersFromInterface() {
    const double pi = 3.1415926535897932384626433832795;
    (*activeSetPtr).pulseEnergy1 = getDoubleFromHWND(maingui.tbPulseEnergy1);
    (*activeSetPtr).pulseEnergy2 = getDoubleFromHWND(maingui.tbPulseEnergy2);
    (*activeSetPtr).frequency1 = 1e12 * getDoubleFromHWND(maingui.tbFrequency1);
    (*activeSetPtr).frequency2 = 1e12 * getDoubleFromHWND(maingui.tbFrequency2);
    (*activeSetPtr).bandwidth1 = 1e12 * getDoubleFromHWND(maingui.tbBandwidth1);
    (*activeSetPtr).bandwidth2 = 1e12 * getDoubleFromHWND(maingui.tbBandwidth2);
    (*activeSetPtr).sgOrder1 = 2 * ((int)ceil(getDoubleFromHWND(maingui.tbPulseType) / 2));
    if ((*activeSetPtr).sgOrder1 < 2) {
        (*activeSetPtr).sgOrder1 = 2;
    }
    (*activeSetPtr).cephase1 = getDoubleFromHWND(maingui.tbCEPhase1);
    (*activeSetPtr).cephase2 = getDoubleFromHWND(maingui.tbCEPhase2);
    (*activeSetPtr).delay1 = -1e-15 * getDoubleFromHWND(maingui.tbPulse1Delay); 
    (*activeSetPtr).delay2 = -1e-15 * getDoubleFromHWND(maingui.tbPulse2Delay);
    (*activeSetPtr).gdd1 = 1e-30 * getDoubleFromHWND(maingui.tbGDD1);
    (*activeSetPtr).gdd2 = 1e-30 * getDoubleFromHWND(maingui.tbGDD2);
    (*activeSetPtr).tod1 = 1e-45 * getDoubleFromHWND(maingui.tbTOD1);
    (*activeSetPtr).tod2 = 1e-45 * getDoubleFromHWND(maingui.tbTOD2);

    (*activeSetPtr).beamwaist1 = 1e-6 * getDoubleFromHWND(maingui.tbBeamwaist1);
    (*activeSetPtr).beamwaist2 = 1e-6 * getDoubleFromHWND(maingui.tbBeamwaist2);
    (*activeSetPtr).x01 = 1e-6 * getDoubleFromHWND(maingui.tbXoffset1);
    (*activeSetPtr).x02 = 1e-6 * getDoubleFromHWND(maingui.tbXoffset2);
    (*activeSetPtr).z01 = 1e-6 * getDoubleFromHWND(maingui.tbZoffset1);
    (*activeSetPtr).z02 = 1e-6 * getDoubleFromHWND(maingui.tbZoffset2);
    (*activeSetPtr).propagationAngle1 = (pi / 180) * getDoubleFromHWND(maingui.tbPropagationAngle1);
    (*activeSetPtr).propagationAngle2 = (pi / 180) * getDoubleFromHWND(maingui.tbPropagationAngle2);
    (*activeSetPtr).polarizationAngle1 = (pi / 180) * getDoubleFromHWND(maingui.tbPolarizationAngle1);
    (*activeSetPtr).polarizationAngle2 = (pi / 180) * getDoubleFromHWND(maingui.tbPolarizationAngle2);
    (*activeSetPtr).circularity1 = getDoubleFromHWND(maingui.tbCircularity1);
    (*activeSetPtr).circularity2 = getDoubleFromHWND(maingui.tbCircularity2);

    (*activeSetPtr).materialIndex = (int)getDoubleFromHWND(maingui.tbMaterialIndex);;
    (*activeSetPtr).crystalTheta = (pi / 180) * getDoubleFromHWND(maingui.tbCrystalTheta);
    (*activeSetPtr).crystalPhi = (pi / 180) * getDoubleFromHWND(maingui.tbCrystalPhi);

    (*activeSetPtr).spatialWidth = 1e-6 * getDoubleFromHWND(maingui.tbGridXdim);
    (*activeSetPtr).rStep = 1e-6 * getDoubleFromHWND(maingui.tbRadialStepSize);
    (*activeSetPtr).timeSpan = 1e-15 * getDoubleFromHWND(maingui.tbTimeSpan);
    (*activeSetPtr).tStep = 1e-15 * getDoubleFromHWND(maingui.tbTimeStepSize);

    //fix the grid dimensions so that each is divisble by 32
    (*activeSetPtr).spatialWidth = (*activeSetPtr).rStep * (32 * ceil((*activeSetPtr).spatialWidth / ((*activeSetPtr).rStep * 32)));
    (*activeSetPtr).timeSpan = (*activeSetPtr).tStep * (32 * ceil((*activeSetPtr).timeSpan / ((*activeSetPtr).tStep * 32)));

    (*activeSetPtr).crystalThickness = 1e-6 * getDoubleFromHWND(maingui.tbCrystalThickness);
    (*activeSetPtr).propagationStep = 1e-9 * getDoubleFromHWND(maingui.tbXstep);

    (*activeSetPtr).nonlinearAbsorptionStrength = getDoubleFromHWND(maingui.tbNonlinearAbsortion);
    (*activeSetPtr).bandGapElectronVolts = getDoubleFromHWND(maingui.tbBandGap);
    (*activeSetPtr).effectiveMass = getDoubleFromHWND(maingui.tbEffectiveMass);
    (*activeSetPtr).drudeGamma = 1e12 * getDoubleFromHWND(maingui.tbDrudeGamma);

    (*activeSetPtr).runType = 1+(int)SendMessage(maingui.pdClusterSelector, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).batchIndex = (int)SendMessage(maingui.pdBatchMode, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).symmetryType = (int)SendMessage(maingui.pdPropagationMode, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);

    (*activeSetPtr).pulse1FileType = (int)SendMessage(maingui.pdPulse1Type, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).pulse2FileType = (int)SendMessage(maingui.pdPulse2Type, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);

    char noneString[] = "None";

    memset((*activeSetPtr).sequenceString, 0, 256 * MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbSequence, (*activeSetPtr).sequenceString, MAX_LOADSTRING * 256);
    if (strnlen_s((*activeSetPtr).sequenceString, 256 * MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).sequenceString, noneString);
    }

    memset((*activeSetPtr).outputBasePath, 0, MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbFileNameBase, (*activeSetPtr).outputBasePath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).outputBasePath, 256 * MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).outputBasePath, noneString);
    }

    memset((*activeSetPtr).field1FilePath, 0, MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbPulse1Path, (*activeSetPtr).field1FilePath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).field1FilePath, 256 * MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).field1FilePath, noneString);
    }

    memset((*activeSetPtr).field2FilePath, 0, MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbPulse2Path, (*activeSetPtr).field2FilePath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).field2FilePath, 256 * MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).field2FilePath, noneString);
    }

    (*activeSetPtr).batchDestination = getDoubleFromHWND(maingui.tbBatchDestination);
    (*activeSetPtr).Nsims = (size_t)getDoubleFromHWND(maingui.tbNumberSims);

    //derived parameters and cleanup:
    (*activeSetPtr).delay1 += (*activeSetPtr).timeSpan / 2;
    (*activeSetPtr).delay2 += (*activeSetPtr).timeSpan / 2;
    (*activeSetPtr).sellmeierType = 0;
    (*activeSetPtr).axesNumber = 0;
    (*activeSetPtr).sgOrder2 = (*activeSetPtr).sgOrder1;
    (*activeSetPtr).Ntime = (size_t)round((*activeSetPtr).timeSpan / (*activeSetPtr).tStep);
    (*activeSetPtr).Nspace = (size_t)round((*activeSetPtr).spatialWidth / (*activeSetPtr).rStep);
    (*activeSetPtr).Ngrid = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace;
    (*activeSetPtr).kStep = 2 * pi / ((*activeSetPtr).Nspace * (*activeSetPtr).rStep);
    (*activeSetPtr).fStep = 1.0 / ((*activeSetPtr).Ntime * (*activeSetPtr).tStep);
    (*activeSetPtr).Npropagation = (size_t)round((*activeSetPtr).crystalThickness / (*activeSetPtr).propagationStep);

    (*activeSetPtr).isCylindric = (*activeSetPtr).symmetryType == 1;
    if ((*activeSetPtr).isCylindric) {
        (*activeSetPtr).x01 = 0;
        (*activeSetPtr).x02 = 0;
        (*activeSetPtr).propagationAngle1 = 0;
        (*activeSetPtr).propagationAngle2 = 0;
    }

    if ((*activeSetPtr).batchIndex == 0 || (*activeSetPtr).batchIndex == 4 || (*activeSetPtr).Nsims < 1) {
        (*activeSetPtr).Nsims = 1;
    }

    (*activeSetPtr).field1IsAllocated = FALSE;
    (*activeSetPtr).field2IsAllocated = FALSE;

    //crystal from database (database must be loaded!)
    (*activeSetPtr).chi2Tensor = crystalDatabasePtr[(*activeSetPtr).materialIndex].d;
    (*activeSetPtr).chi3Tensor = crystalDatabasePtr[(*activeSetPtr).materialIndex].chi3;
    (*activeSetPtr).nonlinearSwitches = crystalDatabasePtr[(*activeSetPtr).materialIndex].nonlinearSwitches;
    (*activeSetPtr).absorptionParameters = crystalDatabasePtr[(*activeSetPtr).materialIndex].absorptionParameters;
    (*activeSetPtr).sellmeierCoefficients = crystalDatabasePtr[(*activeSetPtr).materialIndex].sellmeierCoefficients;
    (*activeSetPtr).sellmeierType = crystalDatabasePtr[(*activeSetPtr).materialIndex].sellmeierType;
    (*activeSetPtr).axesNumber = crystalDatabasePtr[(*activeSetPtr).materialIndex].axisType;



    return 0;
}

//quality of life function - put a text label on a text box window, relative to its position
int labelTextBox(HDC hdc, HWND parentWindow, HWND targetTextBox, const wchar_t* labelText, int xOffset, int yOffset) {
    RECT rectTextBox;
    POINT positionTextBox;
    GetWindowRect(targetTextBox, &rectTextBox);
    positionTextBox.x = rectTextBox.left;
    positionTextBox.y = rectTextBox.top;
    ScreenToClient(parentWindow, &positionTextBox);
    TextOutW(hdc, positionTextBox.x + xOffset, positionTextBox.y + yOffset, labelText, (int)_tcslen(labelText));
    return 0;
}

int floatyText(HDC hdc, HWND parentWindow, const wchar_t* labelText, int xOffset, int yOffset) {
    TextOutW(hdc, xOffset, yOffset, labelText, (int)_tcslen(labelText));
    return 0;
}

//reads the content of a text box and returns a double containing its numerical value
double getDoubleFromHWND(HWND inputA)
{
    double sdata;
    int len = GetWindowTextLength(inputA);
    TCHAR s[100];
    if (len > 0)
    {
        len = GetWindowText(inputA, s, 100);
        sdata = _wtof(s);
        return sdata;
    }
    else {
        return 0.;
    }
}

//Add a text string contained in messageBuffer to the text box inputA
int appendTextToWindow(HWND inputA, wchar_t* messageString, int buffersize) {
    int len = (int)GetWindowTextLength(inputA);
    wchar_t* newbuffer = (wchar_t*)calloc(2 * ((size_t)(len) + buffersize), sizeof(wchar_t));
    if (newbuffer != NULL) {
        if (len > 0) {
            len = GetWindowText(inputA, newbuffer, len + 1);
        }
        memcpy(newbuffer + len, messageString, buffersize * sizeof(wchar_t));
        SetWindowText(inputA, newbuffer);
        SendMessage(inputA, EM_LINESCROLL, 0, 99999);
        free(newbuffer);
    }
    return 0;
}


//template function that works as a wrapper for swprintf_s, for writing to a text control working as a console
//don't give it a format string approaching the size of MAX_LOADSTRING, but come on, that's over a thousand characters
template<typename... Args> void printToConsole(HWND console, const wchar_t* format, Args... args) {
    wchar_t newBuffer[MAX_LOADSTRING] = { 0 };
    swprintf_s(newBuffer, MAX_LOADSTRING, format, args...);
    appendTextToWindow(console, newBuffer, MAX_LOADSTRING);
}


//returns a string containing the text in a text box
int getStringFromHWND(HWND inputA, char* outputString, int bufferSize)
{
    int len = GetWindowTextLength(inputA);
    if (len > 0)
    {
        len = GetWindowTextA(inputA, outputString, bufferSize);
    }
    return 0;
}

int drawLabels(HDC hdc) {
    int labos = -160;
    int x0plots = 380;
    int dxplots = 404;
    int dyplots = 214;
    int vs = 26;

    labelTextBox(hdc, maingui.mainWindow, maingui.tbMaterialIndex, _T("Material index"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCrystalTheta, _T("Crystal theta (deg)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCrystalPhi, _T("Crystal phi (deg)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbNonlinearAbsortion, _T("NL absorption"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbBandGap, _T("Band gap (eV)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbDrudeGamma, _T("Gamma (THz)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbEffectiveMass, _T("Eff. Mass (rel.)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbBeamwaist1, _T("Beamwaist 1 (mcr.)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbBeamwaist2, _T("Beamwaist 2 (mcr.)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbGridXdim, _T("Grid width (mcr.)"), labos, 0);

    labelTextBox(hdc, maingui.mainWindow, maingui.tbTimeStepSize, _T("dt (fs)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbRadialStepSize, _T("dx or dr (mcr.)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbTimeSpan, _T("Time span (fs)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbXstep, _T("dz (nm)"), labos, 0); 
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCrystalThickness, _T("Thickness (mcr.)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdBatchMode, _T("Batch mode"), labos, 4);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbNumberSims, _T("Batch steps"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbBatchDestination, _T("Batch end"), labos, 0);

    labelTextBox(hdc, maingui.mainWindow, maingui.tbPulse1Delay, _T("Delay 1 (fs)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPulse2Delay, _T("Delay 2 (fs)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPulseEnergy1, _T("Energy 1 (J)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPulseEnergy2, _T("Energy 2 (J)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbBandwidth1, _T("Bandwidth 1 (THz)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbBandwidth2, _T("Bandwidth 2 (THz)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbFrequency1, _T("Frequency 1 (THz)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbFrequency2, _T("Frequency 2 (THz)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCEPhase1, _T("CEP/pi 1"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCEPhase2, _T("CEP/pi 2"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPulseType, _T("SG order"), labos, 0);
    
    labelTextBox(hdc, maingui.mainWindow, maingui.tbGDD1, _T("GDD 1 (fs^2)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbGDD2, _T("GDD 2 (fs^2)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbTOD1, _T("TOD 1 (fs^3)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbTOD2, _T("TOD 2 (fs^3)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbXoffset1, _T("x offset 1 (mcr.)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbXoffset2, _T("x offset 2 (mcr.)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbZoffset1, _T("z offset 1 (mcr.)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbZoffset2, _T("z offset 2 (mcr.)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPropagationAngle1, _T("NC angle 1 (deg)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPropagationAngle2, _T("NC angle 2 (deg)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPolarizationAngle1, _T("Polarization 1 (deg)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPolarizationAngle2, _T("Polarization 2 (deg)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCircularity1, _T("Circularity 1"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCircularity2, _T("Circularity 2"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdPulse1Type, _T("Pulse 1 type:"), labos, 4);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdPulse2Type, _T("Pulse 2 type:"), labos, 4);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdPropagationMode, _T("Propagation mode"), labos, 4);

    labelTextBox(hdc, maingui.mainWindow, maingui.tbSequence, _T("Crystal sequence:"), 4, -24);

    //plot labels
    RECT mainRect;
    GetWindowRect(maingui.mainWindow, &mainRect);
    int x0 = maingui.consoleSize + 10;
    int y0 = 2 * maingui.vs;
    int textHeight = 21;
    int imagePanelSizeX = mainRect.right - mainRect.left - x0;
    int imagePanelSizeY = mainRect.bottom - mainRect.top - y0 -5*textHeight;

    int column1X = x0;
    int column2X = x0 + imagePanelSizeX/2;
    int unitLabelRow1X = x0 + imagePanelSizeX / 4;
    int unitLabelRow2X = x0 + (3 * imagePanelSizeX) / 4;
    int row1Y = y0;
    int row2Y = y0 + imagePanelSizeY/4;
    int row3Y = y0 + imagePanelSizeY/2;
    int row4Y = y0 + (3 * imagePanelSizeY) / 4;
    int row5Y = y0 + imagePanelSizeY + 15;
    int dy = 256;
    int plotMargin = 0;
    int spacerX = 50;
    int spacerY = 40;
    floatyText(hdc, maingui.mainWindow, L"s-polarization, space/time:", column1X, row1Y);
    floatyText(hdc, maingui.mainWindow, L"p-polarization, space/time:", column1X, row2Y);
    floatyText(hdc, maingui.mainWindow, L"s-polarization waveform (GV/m):", column1X, row3Y);
    floatyText(hdc, maingui.mainWindow, L"p-polarization waveform (GV/m):", column1X, row4Y);
    floatyText(hdc, maingui.mainWindow, L"Time (fs)", unitLabelRow1X-40, row5Y);
    floatyText(hdc, maingui.mainWindow, L"s-polarization, Fourier, Log:", column2X, row1Y);
    floatyText(hdc, maingui.mainWindow, L"p-polarization, Fourier, Log:", column2X, row2Y);
    floatyText(hdc, maingui.mainWindow, L"s-polarization spectrum, log-scale:", column2X, row3Y);
    floatyText(hdc, maingui.mainWindow, L"p-polarization spectrum, log-scale:", column2X, row4Y);
    floatyText(hdc, maingui.mainWindow, L"Frequency (THz)", unitLabelRow2X-90, row5Y);
    return 0;
}

int getFileNameBaseFromDlg(HWND hWnd, HWND outputTextbox) {

    //create the dialog box and get the file path
    OPENFILENAME ofn;
    TCHAR szFileName[MAX_PATH];
    TCHAR szFileNameNoExt[MAX_PATH];
    ZeroMemory(&ofn, sizeof(ofn));
    WORD fbasedirend;
    WORD fbaseloc = 0;
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = hWnd;
    ofn.lpstrFilter = TEXT("Text Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0");
    ofn.lpstrFile = szFileName;
    ofn.nMaxFile = MAX_PATH;
    ofn.Flags = OFN_EXPLORER;
    ofn.lpstrDefExt = TEXT("txt");
    ofn.lpstrFile[0] = '\0';
    ofn.nFileExtension = 0;
    ofn.nFileOffset = 0;
    //only save if there is memory allocated and the simulation is complete
    if (GetSaveFileNameW(&ofn)) {

        //get the base of the file name, so that different files can be made with different extensions based on that
        if (ofn.nFileExtension > 0) {
            fbaseloc = ofn.nFileExtension - 1;
            fbasedirend = ofn.nFileOffset;
        }
        _tcsncpy_s(szFileNameNoExt, szFileName, fbaseloc);
        szFileNameNoExt[MAX_PATH - 1] = 0;
        SetWindowText(outputTextbox, szFileNameNoExt);
    }

    return 0;
}


int getFileNameBaseFromDlgDat(HWND hWnd, HWND outputTextbox) {

    //create the dialog box and get the file path
    OPENFILENAME ofn;
    TCHAR szFileName[MAX_PATH];
    ZeroMemory(&ofn, sizeof(ofn));
    WORD fbaseloc = 0;
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = hWnd;
    ofn.lpstrFilter = TEXT("Dat Files (*.dat)\0*.dat\0All Files (*.*)\0*.*\0");
    ofn.lpstrFile = szFileName;
    ofn.nMaxFile = MAX_PATH;
    ofn.Flags = OFN_EXPLORER;
    ofn.lpstrDefExt = TEXT("dat");
    ofn.lpstrFile[0] = '\0';
    ofn.nFileExtension = 0;
    ofn.nFileOffset = 0;

    if (GetSaveFileNameW(&ofn)) {

        SetWindowText(outputTextbox, szFileName);
    }

    return 0;
}

int openDialogBoxAndLoad(HWND hWnd) {
    WORD fbasedirend;
    OPENFILENAME ofn;
    TCHAR szFileName[MAX_LOADSTRING];
    TCHAR szFileNameNoExt[MAX_LOADSTRING];
    ZeroMemory(&ofn, sizeof(ofn));
    WORD fbaseloc = 0;
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = hWnd;
    ofn.lpstrFilter = TEXT("Text Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0");
    ofn.lpstrFile = szFileName;
    ofn.nMaxFile = MAX_LOADSTRING;
    ofn.Flags = OFN_EXPLORER;
    ofn.lpstrDefExt = TEXT("dat");
    ofn.lpstrFile[0] = '\0';
    ofn.nFileExtension = 0;
    ofn.nFileOffset = 0;
    char fileNameString[MAX_LOADSTRING];
    wcstombs(fileNameString, szFileName, MAX_LOADSTRING);
    if (GetSaveFileNameW(&ofn)) {
        wcstombs(fileNameString, szFileName, MAX_LOADSTRING);
        if (isGridAllocated) {
            freeSemipermanentGrids();
            isGridAllocated = FALSE;
        }
        readInputParametersFile(activeSetPtr, crystalDatabasePtr, fileNameString);

        //get the base of the file name, so that different files can be made with different extensions based on that
        if (ofn.nFileExtension > 0) {
            fbaseloc = ofn.nFileExtension - 1;
            fbasedirend = ofn.nFileOffset;
        }
        _tcsncpy_s(szFileNameNoExt, szFileName, fbaseloc);
        szFileNameNoExt[MAX_LOADSTRING - 1] = 0;
        wcstombs(fileNameString, szFileNameNoExt, MAX_LOADSTRING);


        allocateGrids(activeSetPtr);
        isGridAllocated = TRUE;
        loadSavedFields(activeSetPtr, fileNameString, FALSE);

    }
    
    return 0;
}

//Take an array of doubles and make an image in the window
//Color maps:
//  cm = 1: grayscale
//  cm = 2: similar to matlab's jet
//  cm = 3: similar to Colorcet L07
//  cn = 4: vaporwave (symmetric amplitude)
int drawArrayAsBitmap(HDC hdc, INT64 Nx, INT64 Ny, INT64 x, INT64 y, INT64 height, INT64 width, float* data, int cm) {

    // creating input
    unsigned char* pixels = (unsigned char*)calloc(4 * Nx * Ny, sizeof(unsigned char));
    if (pixels == NULL) {
        return 1;
    }
    INT64 i;
    INT64 Ntot = Nx * Ny;
    float nval;
    int stride = 4;
    //Find the image maximum and minimum
    float imin = data[0];
    float imax = data[0];
    for (i = 1; i < Ntot; i++) {
        if (data[i] > imax) imax = data[i];
        if (data[i] < imin) imin = data[i];
    }

    if (imin != imax) {
        if (cm == 0) {
            for (i = 0; i < Ntot; i++) {
                nval = (data[i] - imin) / (imax - imin);
                pixels[stride * i + 0] = (unsigned char)(255 * nval); //blue channel
                pixels[stride * i + 1] = (unsigned char)(255 * nval); //green channel
                pixels[stride * i + 2] = (unsigned char)(255 * nval); //red channel
            }
        }
        if (cm == 1) {
            for (i = 0; i < Ntot; i++) {
                nval = (data[i] - imin) / (imax - imin);
                pixels[stride * i + 0] = (unsigned char)(255 * cos(3.1415 * nval / 2.)); //blue channel
                pixels[stride * i + 1] = (unsigned char)(255 * cos(3.1415 * (nval - 0.5))); //green channel
                pixels[stride * i + 2] = (unsigned char)(255 * sin(3.1415 * nval / 2.)); //red channel
            }
        }
        if (cm == 2) {
            for (i = 0; i < Ntot; i++) {
                nval = (data[i] - imin) / (imax - imin);
                pixels[stride * i + 0] = (unsigned char)(255 * cos(3.1415 * nval / 2.)); //blue channel
                pixels[stride * i + 1] = (unsigned char)(255 * cos(3.1415 * (nval - 0.5))); //green channel
                pixels[stride * i + 2] = (unsigned char)(255 * sin(3.1415 * nval / 2.)); //red channel
                if (nval < 0.02) {
                    pixels[stride * i + 0] = 255;
                    pixels[stride * i + 1] = 128;
                    pixels[stride * i + 2] = 128;
                }
                if (nval < 0.01) {
                    pixels[stride * i + 0] = 255;
                    pixels[stride * i + 1] = 255;
                    pixels[stride * i + 2] = 255;
                }
            }
        }
        if (cm == 3) {
            for (i = 0; i < Ntot; i++) {
                nval = 255*(data[i] - imin) / (imax - imin);
                
                pixels[stride * i + 0] = (unsigned char)(255 *
                    (0.998*exp(-pow(7.7469e-03 * (nval - 160),6))
                    + 0.22*exp(-pow(0.016818 * (nval - 305), 4)))); //blue channel
                
                
                pixels[stride * i + 1] = (unsigned char)(255 *
                    (0.022 * exp(-pow(0.042045*(nval - 25), 4))
                    + 0.11 * exp(-pow(0.015289*(nval - 120), 4))
                    + 1 * exp(-pow(4.6889e-03*(nval - 400), 6)))); //green channel
                pixels[stride * i + 2] = (unsigned char)(255 *
                    (exp(-pow(3.1101e-03*(nval - 415), 10)))); //red channel
                    
            }
        }
        if (cm == 4) {
            imax = max(imax, -imin);
            imin = min(imin, -imax);
            for (i = 0; i < Ntot; i++) {
                
                nval = (data[i] - imin) / (imax - imin);
                pixels[stride * i + 0] = (unsigned char)(255 * (1.00 * exp(-pow(4*(nval-0.05),2))
                    + 1 * exp(-pow(4 * (nval - 1.05), 2)))); //blue channel
                pixels[stride * i + 1] = (unsigned char)(255 * (1.02 * exp(-pow(3.5*(nval - 1.05), 2)))); //green channel
                pixels[stride * i + 2] = (unsigned char)(255 * (0.8 * exp(-pow(4*(nval - 0.05), 2)))); //red channel
            }
        }
    }
    BITMAPINFOHEADER bmih;

    bmih.biSize = sizeof(BITMAPINFOHEADER);
    bmih.biWidth = (long)Nx;
    bmih.biHeight = (long)Ny;
    bmih.biPlanes = 1;
    bmih.biBitCount = 32;
    bmih.biCompression = BI_RGB;
    bmih.biSizeImage = 0;
    bmih.biXPelsPerMeter = 0;
    bmih.biYPelsPerMeter = 0;
    bmih.biClrUsed = 0;
    bmih.biClrImportant = 0;

    BITMAPINFO dbmi;
    ZeroMemory(&dbmi, sizeof(dbmi));
    dbmi.bmiHeader = bmih;
    dbmi.bmiColors->rgbBlue = 0;
    dbmi.bmiColors->rgbGreen = 0;
    dbmi.bmiColors->rgbRed = 0;
    dbmi.bmiColors->rgbReserved = 0;

    HBITMAP hbmp;
    ZeroMemory(&hbmp, sizeof(HBITMAP));
    hbmp = CreateDIBitmap(hdc, &bmih, CBM_INIT, pixels, &dbmi, DIB_RGB_COLORS);
    BITMAP bm;
    ZeroMemory(&bm, sizeof(BITMAP));
    HDC hdcMem = CreateCompatibleDC(hdc);

    SelectObject(hdcMem, hbmp);
    GetObject(hbmp, sizeof(bm), &bm);

    //BitBlt(hdc, x, y, bm.bmWidth, bm.bmHeight, hdcMem, 0, 0, SRCCOPY);
    StretchBlt(hdc, (int)x, (int)y, (int)width, (int)height, hdcMem, 0, 0, bm.bmWidth, bm.bmHeight, SRCCOPY);
    free(pixels);
    DeleteDC(hdcMem);
    return 0;
}

DWORD WINAPI drawSimPlots(LPVOID lpParam) {
    if (isGridAllocated) {
        isPlotting = TRUE;
        RECT mainRect;
        GetWindowRect(maingui.mainWindow, &mainRect);

        int spacerX = 50;
        int spacerY = 40;
        int x0 = maingui.consoleSize + spacerX + 10;
        int y0 = 90 + spacerY;
        int imagePanelSizeX = mainRect.right - mainRect.left - x0 - 2*spacerX;
        int imagePanelSizeY = mainRect.bottom - mainRect.top - y0 - 5*spacerY;

        size_t simIndex = (*activeSetPtr).plotSim;
        int x = x0;
        int y = y0;
        int dx = imagePanelSizeX/2;
        int dy = imagePanelSizeY/4;
        int plotMargin = 0;
        int xCorrection = 10;
        int yCorrection = 45;

        ShowWindow(maingui.plotBox1, 0);
        ShowWindow(maingui.plotBox2, 0);
        ShowWindow(maingui.plotBox5, 0);
        ShowWindow(maingui.plotBox6, 0);
        SetWindowPos(maingui.plotBox3, HWND_BOTTOM, x - xCorrection, y + 2 * dy + 2 * spacerY - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox4, HWND_BOTTOM, x - xCorrection, y + 3 * dy + 3 * spacerY - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox7, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 2 * dy + 2 * spacerY - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox8, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 3 * dy + 3 * spacerY - yCorrection, dx, dy, NULL);


        float logPlotOffset = 1e11f;
        int i,j;
        HDC hdc;

        if (simIndex < 0) {
            simIndex = 0;
        }
        if (simIndex > (*activeSetPtr).Nsims) {
            simIndex = (*activeSetPtr).Nsims - 1;
        }

        hdc = GetWindowDC(maingui.mainWindow);
        float* plotarr = (float*)calloc((*activeSetPtr).Ngrid, sizeof(float));
        float* plotarrC = (float*)calloc((*activeSetPtr).Ntime * (*activeSetPtr).Nspace, sizeof(float));
        float* plotarr2 = (float*)calloc((size_t)(dx) * (size_t)(dy), sizeof(float));
        
        std::complex<double>* shiftedFFT = (std::complex<double>*)calloc((*activeSetPtr).Ngrid, sizeof(std::complex<double>));

        //Plot Time Domain, s-polarization
        for (i = 0; i < (*activeSetPtr).Ngrid; i++) {
            plotarr[i] = (float)(real((*activeSetPtr).ExtOut[i + simIndex * (*activeSetPtr).Ngrid * 2]));
        }
        linearRemap(plotarr, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime, plotarr2, (int)dy, (int)dx);
        drawArrayAsBitmap(hdc, dx, dy, x, y, dy, dx, plotarr2, 4);

        plotXYDirect2d(maingui.plotBox3, (float)(*activeSetPtr).tStep / 1e-15f, &plotarr[(*activeSetPtr).Ngrid / 2], (*activeSetPtr).Ntime, 1e9, FALSE, 0);

        //Plot Time Domain, p-polarization
        for (i = 0; i < (*activeSetPtr).Ngrid; i++) {
            plotarr[i] = (float)(real((*activeSetPtr).ExtOut[i + (*activeSetPtr).Ngrid + simIndex * (*activeSetPtr).Ngrid * 2]));
        }

        linearRemap(plotarr, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime, plotarr2, (int)dy, (int)dx);
        drawArrayAsBitmap(hdc, dx, dy, x, (size_t)(y) + (size_t)(dy) + spacerY, dy, dx, plotarr2, 4);
        plotXYDirect2d(maingui.plotBox4, (float)(*activeSetPtr).tStep / 1e-15f, &plotarr[(*activeSetPtr).Ngrid / 2], (*activeSetPtr).Ntime, 1e9, FALSE, 0);

        //Plot Fourier Domain, s-polarization
        fftshiftZ(&(*activeSetPtr).EkwOut[simIndex * (*activeSetPtr).Ngrid * 2], shiftedFFT, (*activeSetPtr).Ntime, (*activeSetPtr).Nspace);
        for (i = 0; i < (*activeSetPtr).Ngrid; i++) {
            plotarr[i] = log10((float)cModulusSquared(shiftedFFT[i]) + logPlotOffset);
        }

        for (i = 0; i < ((*activeSetPtr).Ntime / 2); i++) {
            for (j = 0; j < (*activeSetPtr).Nspace; j++) {
                plotarrC[i + ((*activeSetPtr).Ntime / 2) * j] = plotarr[i + (*activeSetPtr).Ntime * j + (*activeSetPtr).Ntime / 2];
            }
        }

        linearRemap(plotarrC, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime / 2, plotarr2, (int)dy, (int)dx);
        drawArrayAsBitmap(hdc, dx, dy, (size_t)(x) + (size_t)(dx) + (size_t)(spacerX), y, dy, dx, plotarr2, 3);
        plotXYDirect2d(maingui.plotBox7, (float)(*activeSetPtr).fStep / 1e12f, &plotarrC[(*activeSetPtr).Ngrid / 4], (*activeSetPtr).Ntime / 2, 1, TRUE, 12);
        
        //Plot Fourier Domain, p-polarization
        fftshiftZ(&(*activeSetPtr).EkwOut[simIndex * (*activeSetPtr).Ngrid * 2 + (*activeSetPtr).Ngrid], shiftedFFT, (*activeSetPtr).Ntime, (*activeSetPtr).Nspace);

        for (i = 0; i < (*activeSetPtr).Ngrid; i++) {
            plotarr[i] = log10((float)cModulusSquared(shiftedFFT[i]) + logPlotOffset);
        }
        for (i = 0; i < ((*activeSetPtr).Ntime/2); i++) {
            for (j = 0; j < (*activeSetPtr).Nspace; j++) {
                plotarrC[i + ((*activeSetPtr).Ntime / 2) * j] = plotarr[i + (*activeSetPtr).Ntime * j + (*activeSetPtr).Ntime/2];
            }
        }
        
        linearRemap(plotarrC, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime/2, plotarr2, (int)dy, (int)dx);
        drawArrayAsBitmap(hdc, dx, dy, (size_t)(x) + (size_t)(dx) + (size_t)(spacerX), (size_t)(y) + (size_t)(dy) + (size_t)(spacerY), dy, dx, plotarr2, 3);

        plotXYDirect2d(maingui.plotBox8, (float)(*activeSetPtr).fStep / 1e12f, &plotarrC[(*activeSetPtr).Ngrid / 4], (*activeSetPtr).Ntime/2, 1, TRUE, 12);

        free(shiftedFFT);
        free(plotarr);
        free(plotarr2);
        free(plotarrC);
        ReleaseDC(maingui.mainWindow, hdc);
        isPlotting = FALSE;
    }
    return 0;
}



//resize matrix A to the size of matrix B
//B is overwritten with the resized matrix
int linearRemap(float* A, int nax, int nay, float* B, int nbx, int nby) {
    int i, j;
    float A00;
    float f;

    int nx0, ny0;
    int Ni, Nj;
	for (i = 0; i < nbx; i++) {
		f = i*(nax / (float)nbx);
		Ni = (int)f;
		nx0 = nay * min(Ni, nax);
		for (j = 0; j < nby; j++) {
			f = (j*(nay / (float)nby));
			Nj = (int)f;
			ny0 = min(nay, Nj);
			A00 = A[ny0 + nx0];
			B[i * nby + j] = A00;
		}
	}
    return 0;
}

void setTitleBarDark(HWND hWnd)
{
    BOOL isDarkMode = TRUE;
    DwmSetWindowAttribute(hWnd, 20, &isDarkMode, sizeof(isDarkMode));
}

void plotXYDirect2d(HWND targetWindow, float dX, float* Y, size_t Npts, float unitY, bool forceminY, float forcedminY) {
    size_t i;
    D2D1_POINT_2F p1;
    D2D1_POINT_2F p2;
    D2D1_ELLIPSE marker;
    D2D1_SIZE_F sizeF;
    ZeroMemory(&sizeF, sizeof(D2D1_SIZE_F));
    ID2D1HwndRenderTarget *pRenderTarget;
    ID2D1SolidColorBrush* pBrush;
    float markerSize = 1.5;
    float lineWidth = 1.25;

    //get limits of Y
    float maxY = 0;
    float minY = 0;
    for (i = 0; i < Npts; i++) {
        maxY = max(Y[i], maxY);
        minY = min(Y[i], minY);
    }
    if (minY == maxY) {
        minY = -1;
        maxY = 1;
    }
    if (forceminY) {
        minY = forcedminY;
    }


    RECT targetRectangle;
    GetClientRect(targetWindow, &targetRectangle);
    D2D1_SIZE_U size = D2D1::SizeU(targetRectangle.right, targetRectangle.bottom);

    HRESULT hr = maingui.pFactory->CreateHwndRenderTarget(
        D2D1::RenderTargetProperties(),
        D2D1::HwndRenderTargetProperties(targetWindow, size),
        &pRenderTarget);
    if (pRenderTarget != NULL) {
        sizeF = pRenderTarget->GetSize();
    }
    
    float scaleX = sizeF.width / (Npts * (float)dX);
    float scaleY = sizeF.height / ((float)(maxY - minY));

    if (SUCCEEDED(hr))
    {
        const D2D1_COLOR_F color = D2D1::ColorF(1, 1, 1, 1);
        hr = pRenderTarget->CreateSolidColorBrush(color, &pBrush);

        if (SUCCEEDED(hr))
        {
            PAINTSTRUCT ps;
            BeginPaint(targetWindow, &ps);
            pRenderTarget->BeginDraw();

            pRenderTarget->Clear(D2D1::ColorF(D2D1::ColorF::Black));
            for (i = 0; i < Npts - 1; i++) {
                p1.x = scaleX * (i * (float)dX);
                p1.y = sizeF.height - scaleY * ((float)Y[i] - (float)minY);
                p2.x = scaleX * ((i + 1) * (float)dX);
                p2.y = sizeF.height - scaleY * ((float)Y[i + 1] - (float)minY);
                pRenderTarget->DrawLine(p1, p2, pBrush, lineWidth, 0);
                marker.point = p1;
                marker.radiusX = markerSize;
                marker.radiusY = markerSize;
                pRenderTarget->FillEllipse(&marker, pBrush);
            }

            hr = pRenderTarget->EndDraw();
            EndPaint(targetWindow, &ps);
        }
        pBrush->Release();
        pRenderTarget->Release();

    }

    //Draw the labels
    int NyTicks = 3;
    int NxTicks = 3;
    wchar_t messageBuffer[MAX_LOADSTRING];
    double yTicks1[3] = { maxY, 0.5 * (maxY + minY), minY };
    double xTicks1[3] = { 0.25 * dX * Npts, 0.5 * dX * Npts, 0.75 * dX * Npts };
    HDC hdc = GetWindowDC(maingui.mainWindow);
    SetTextColor(hdc, uiGrey);
    SetBkColor(hdc, uiDarkGrey);
    RECT windowRect;
    GetWindowRect(targetWindow, &windowRect);

    int posX = windowRect.left;
    int posY = windowRect.top;
    int pixelsTall = windowRect.bottom - windowRect.top;
    int pixelsWide = windowRect.right - windowRect.left;

    GetWindowRect(maingui.mainWindow, &windowRect);
    posX -= windowRect.left;
    posY -= windowRect.top;
    for (i = 0; i < NyTicks; i++) {
        memset(messageBuffer, 0, MAX_LOADSTRING * sizeof(wchar_t));
        swprintf_s(messageBuffer, MAX_LOADSTRING,
            _T("%1.1f"), yTicks1[i] / unitY);
        TextOutW(hdc, posX - 32, posY + (int)(i * 0.96 * pixelsTall / 2), messageBuffer, (int)_tcslen(messageBuffer));
    }
    for (i = 0; i < 3; i++) {
        memset(messageBuffer, 0, MAX_LOADSTRING * sizeof(wchar_t));
        swprintf_s(messageBuffer, MAX_LOADSTRING,
            _T("%3.0f"), xTicks1[i]);
        TextOutW(hdc, posX + (int)(0.25 * pixelsWide * ((size_t)(i)+1) - 12), posY + pixelsTall, messageBuffer, (int)_tcslen(messageBuffer));
    }
    ReleaseDC(maingui.mainWindow, hdc);

}