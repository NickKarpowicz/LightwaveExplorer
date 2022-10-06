#pragma comment(linker,"\"/manifestdependency:type='win32' \
name='Microsoft.Windows.Common-Controls' version='6.0.0.0' \
processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

#include "framework.h"
#include "LightwaveExplorerFrontend.h"
#include "LightwaveExplorerCore.cuh"
#include "LightwaveExplorerCoreCPU.h"
#include "LightwaveExplorerUtilities.h"
#include "LightwaveExplorerSYCL/LightwaveExplorerSYCL.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <Commdlg.h>
#include <chrono>
#include <Windows.h>
#include <Uxtheme.h>
#include <dwmapi.h>
#include <d2d1.h>
#include <nvml.h>
#include <thread>

#define ID_BTNRUN 11110
#define ID_BTNPLOT 11111
#define ID_BTNGETFILENAME 11112
#define ID_BTNREFRESHDB 11113
#define ID_BTNSTOP 11114
#define ID_BTNPULSE1 11115
#define ID_BTNPULSE2 11116
#define ID_BTNRUNONCLUSTER 11117
#define ID_BTNLOAD 11118
#define ID_BTNFIT 11119
#define ID_BTNFITREFERENCE 11120
#define ID_BTNADDECHO 11121
#define ID_BTNADDCRYSTAL 11122
#define ID_CBLOGPLOT 12000
#define ID_CBFORCECPU 12001
#define ID_PLOTSCRUBBER 13001
#define MAX_SIMULATIONS 16192
#define MIN_GRIDDIM 8
#define TWOPI 6.2831853071795862
#define PI 3.1415926535897931
#define DEG2RAD 1.7453292519943295e-02
#define RAD2DEG 57.2957795130823229
#define LIGHTC 2.99792458e8
#define EPS0 8.8541878128e-12
#define SIXTH 0.1666666666666667
#define THIRD 0.3333333333333333

// Global Variables:
HINSTANCE hInst;                            // current instance
WCHAR szTitle[MAX_LOADSTRING];              // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];        // the main window class name
guiStruct maingui;                          // struct containing all of the GUI elements
simulationParameterSet* activeSetPtr;       // Main structure containing simulation parameters and pointers
crystalEntry* crystalDatabasePtr;           // Crystal info database
WCHAR programDirectory[MAX_LOADSTRING];     // Program working directory (useful if the crystal database has to be reloaded)
bool isRunning = FALSE;
bool isPlotting = FALSE;
bool isGridAllocated = FALSE;
bool cancellationCalled = FALSE;
bool hasGPU = FALSE;
int cudaCount = 0;
size_t progressCounter = 0;
//UI colors
COLORREF uiWhite = RGB(255, 255, 255);
COLORREF uiGrey = RGB(216, 216, 216);
COLORREF uiDarkGrey = RGB(32, 32, 32);
COLORREF uiBlack = RGB(16, 16, 16);
HBRUSH blackBrush = CreateSolidBrush(uiBlack);
HBRUSH greyBrush = CreateSolidBrush(uiDarkGrey);
// Forward declarations of (Microsoft) functions
ATOM                MyRegisterClass(HINSTANCE hInstance);
bool                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);


DWORD WINAPI mainSimThread(LPVOID lpParam) {
    cancellationCalled = FALSE;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();
    HANDLE plotThread = NULL;
    DWORD hplotThread;
    HANDLE cpuThread = NULL;
    DWORD hCpuThread;

    //bool forcingCPU = IsDlgButtonChecked(maingui.mainWindow, ID_CBFORCECPU) == BST_CHECKED;
    bool forcingCPU = 0;
    if (forcingCPU) {
        printToConsole(maingui.textboxSims, L"Forcing to run on CPU!\r\n");
    }
    if (isGridAllocated) {
        freeSemipermanentGrids();
    }
    memset(activeSetPtr, 0, sizeof(simulationParameterSet));
    readParametersFromInterface();
    if ((*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 > MAX_SIMULATIONS) {
        printToConsole(maingui.textboxSims, L"Too many simulations in batch mode. Must be under %i total.\r\n", MAX_SIMULATIONS);
    }
    (*activeSetPtr).runType = 0;
    allocateGrids(activeSetPtr);
    isGridAllocated = TRUE;
    setTrackbarLimitsToActiveSet();
    (*activeSetPtr).isFollowerInSequence = FALSE;
    (*activeSetPtr).crystalDatabase = crystalDatabasePtr;
    loadPulseFiles(activeSetPtr);
    readSequenceString(activeSetPtr);
    configureBatchMode(activeSetPtr);
    int error = 0;
    //run the simulations
    isRunning = TRUE;
    progressCounter = 0;
    cpuThread = CreateThread(NULL, 0, offloadToCPU, &activeSetPtr[(*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 - (*activeSetPtr).NsimsCPU], 0, &hCpuThread);

    int pulldownSelection = (int)SendMessage(maingui.pdPrimaryQueue, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    auto sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
    auto normalFunction = &solveNonlinearWaveEquationCPU;
    int assignedGPU = 0;

    if ((pulldownSelection - cudaCount) == 0) {
        sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
        normalFunction = &solveNonlinearWaveEquationSYCL;
    }
    else if ((pulldownSelection - cudaCount) == 1) {
        sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
        normalFunction = &solveNonlinearWaveEquationCPU;
    }
    else if ((pulldownSelection - cudaCount) < 0) {
        sequenceFunction = &solveNonlinearWaveEquationSequence;
        normalFunction = &solveNonlinearWaveEquation;
        assignedGPU = pulldownSelection;
    }
 
    for (int j = 0; j < ((*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 - (*activeSetPtr).NsimsCPU); j++) {
        if ((*activeSetPtr).isInSequence) {
            error = sequenceFunction(&activeSetPtr[j]);

            if (activeSetPtr[j].memoryError != 0) {
                if (activeSetPtr[j].memoryError == -1) {
                    printToConsole(maingui.textboxSims, _T("Not enough free GPU memory, sorry.\r\n"), activeSetPtr[j].memoryError);
                }
                else {
                    printToConsole(maingui.textboxSims, _T("Warning: device memory error (%i).\r\n"), activeSetPtr[j].memoryError);
                } 
            }
            if (error) break;
            if (!isPlotting) {
                (*activeSetPtr).plotSim = j;
                plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);
                if(plotThread != 0)CloseHandle(plotThread);
            }

        }
        else {
            error = normalFunction(&activeSetPtr[j]); 
            if (activeSetPtr[j].memoryError != 0) {
                if (activeSetPtr[j].memoryError == -1) {
                    printToConsole(maingui.textboxSims, _T("Not enough free GPU memory, sorry.\r\n"), activeSetPtr[j].memoryError);
                }
                else {
                    printToConsole(maingui.textboxSims, _T("Warning: device memory error (%i).\r\n"), activeSetPtr[j].memoryError);
                }
            }

            if (error) break;
        }


        if (cancellationCalled) {
            printToConsole(maingui.textboxSims, _T("Warning: series cancelled, stopping after %i simulations.\r\n"), j + 1);
            break;
        }

        if (!isPlotting) {
            (*activeSetPtr).plotSim = j;
            plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);
            if(plotThread != NULL)CloseHandle(plotThread);
        }

    }

    if ((*activeSetPtr).NsimsCPU != 0 && cpuThread != 0) {
        WaitForSingleObject(cpuThread, INFINITE);
        CloseHandle(cpuThread);
    }
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    if (error==13) {
        printToConsole(maingui.textboxSims, 
            L"NaN detected in grid!\r\nTry using a larger spatial/temporal step\r\nor smaller propagation step.\r\nSimulation was cancelled.\r\n");
    }
    else {
        printToConsole(maingui.textboxSims, _T("Finished after %8.4lf s. \r\n"), 1e-6 *
            (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
    }

    saveDataSet(activeSetPtr, crystalDatabasePtr, (*activeSetPtr).outputBasePath, FALSE);
    deallocateGrids(activeSetPtr, FALSE);
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
        L"Run %ls on cluster with:\r\nsbatch %ls.slurmScript\r\n",
        wideBuffer, wideBuffer);
    
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
    ShowWindow(maingui.plotBox1, nCmdShow);
    ShowWindow(maingui.plotBox2, nCmdShow);
    ShowWindow(maingui.plotBox3, nCmdShow);
    ShowWindow(maingui.plotBox4, nCmdShow);
    ShowWindow(maingui.plotBox5, nCmdShow);
    ShowWindow(maingui.plotBox6, nCmdShow);
    ShowWindow(maingui.plotBox7, nCmdShow);
    ShowWindow(maingui.plotBox8, nCmdShow);
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

//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;
    memset(&wcex, 0, sizeof(WNDCLASSEXW));
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
    int textboxwidth = maingui.textboxwidth;
    int halfBox = textboxwidth / 2 - 4;
    int xOffsetRow1 = maingui.xOffsetRow1;
    int xOffsetRow1b = xOffsetRow1 + halfBox + 8;
    int xOffsetRow2 = maingui.xOffsetRow2;
    int xOffsetRow2b = xOffsetRow2 + halfBox + 8;
    int xOffsetRow3 = maingui.xOffsetRow3;
    int vs = maingui.vs;
    int btnwidth = maingui.btnwidth;
    int btnoffset2 = maingui.btnoffset2;
    int btnoffset2a = maingui.btnoffset2a;
    int consoleSize = maingui.consoleSize;
    
    int btnHeight = maingui.btnHeight;
    maingui.mainWindow = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW | WS_EX_CONTROLPARENT,
        CW_USEDEFAULT, CW_USEDEFAULT, 1900, 1000, nullptr, nullptr, hInstance, nullptr);
    SetMenu(maingui.mainWindow, NULL);
    SetWindowTextA(maingui.mainWindow, "Lightwave Explorer");

    //text boxes for input parameters
    maingui.tbPulseEnergy1 = CreateWindow(WC_EDIT, TEXT("24e-9"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 0 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulseEnergy2 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 0 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFrequency1 = CreateWindow(WC_EDIT, TEXT("130"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 1 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFrequency2 = CreateWindow(WC_EDIT, TEXT("200"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 1 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBandwidth1 = CreateWindow(WC_EDIT, TEXT("40"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 2 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBandwidth2 = CreateWindow(WC_EDIT, TEXT("40"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 2 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulseType1 = CreateWindow(WC_EDIT, TEXT("2"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 3 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulseType2 = CreateWindow(WC_EDIT, TEXT("2"),
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow1b, 3 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCEPhase1 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 4 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCEPhase2 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 4 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulse1Delay = CreateWindow(WC_EDIT, TEXT("-90"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 5 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulse2Delay = CreateWindowW(WC_EDIT, TEXT("-20"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 5 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbGDD1 = CreateWindow(WC_EDIT, TEXT("-40.262"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 6 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbGDD2 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 6 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTOD1 = CreateWindow(WC_EDIT, TEXT("150"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 7 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTOD2 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 7 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.tbPhaseMaterialIndex1 = CreateWindow(WC_EDIT, L"0",
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow1, 8 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPhaseMaterialIndex2 = CreateWindow(WC_EDIT, L"0",
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow1b, 8 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPhaseMaterialThickness1 = CreateWindow(WC_EDIT, L"0",
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow1, 9 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPhaseMaterialThickness2 = CreateWindow(WC_EDIT, L"0",
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow1b, 9 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.tbBeamwaist1 = CreateWindow(WC_EDIT, TEXT("3.2"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 10 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBeamwaist2 = CreateWindow(WC_EDIT, TEXT("50"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 10 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbXoffset1 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 11 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbXoffset2 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 11 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbZoffset1 = CreateWindow(WC_EDIT, TEXT("-100"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 12 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbZoffset2 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 12 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPropagationAngle1 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 13 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPropagationAngle2 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 13 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPolarizationAngle1 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 14 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPolarizationAngle2 = CreateWindow(WC_EDIT, TEXT("90"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 14 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCircularity1 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1, 15 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCircularity2 = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1b, 15 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFileNameBase = CreateWindow(WC_EDIT, TEXT("TestFile"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 1 * vs+5, 775, 20, maingui.mainWindow, NULL, hInstance, NULL);
    
    maingui.plotBox1 = CreateWindow(WC_STATIC, NULL, 
        WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox2 = CreateWindow(WC_STATIC, NULL, 
        WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox3 = CreateWindow(WC_STATIC, NULL, 
        WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox4 = CreateWindow(WC_STATIC, NULL, 
        WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox5 = CreateWindow(WC_STATIC, NULL, 
        WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox6 = CreateWindow(WC_STATIC, NULL, 
        WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox7 = CreateWindow(WC_STATIC, NULL, 
        WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.plotBox8 = CreateWindow(WC_STATIC, NULL, 
        WS_CHILD | WS_VISIBLE | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 3 * vs, 50, 50, maingui.mainWindow, NULL, hInstance, NULL);
    D2D1CreateFactory(D2D1_FACTORY_TYPE_SINGLE_THREADED, &maingui.pFactory);

    maingui.tbMaterialIndex = CreateWindow(WC_EDIT, TEXT("3"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2, 0 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbMaterialIndexAlternate = CreateWindow(WC_EDIT, TEXT("3"),
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow2b, 0 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCrystalTheta = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2, 1 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCrystalPhi = CreateWindow(WC_EDIT, TEXT("0"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2b, 1 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.tbNonlinearAbsortion = CreateWindow(WC_EDIT, TEXT("24e-87"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2, 2 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBandGap = CreateWindow(WC_EDIT, TEXT("3"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2b, 2 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbDrudeGamma = CreateWindow(WC_EDIT, TEXT("10"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2, 3 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbEffectiveMass = CreateWindow(WC_EDIT, TEXT("0.081"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2b, 3 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.tbGridXdim = CreateWindow(WC_EDIT, TEXT("210"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_AUTOHSCROLL,
        xOffsetRow2, 4 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbRadialStepSize = CreateWindow(WC_EDIT, TEXT("0.62"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2b, 4 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTimeSpan = CreateWindow(WC_EDIT, TEXT("512"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2, 5 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTimeStepSize = CreateWindow(WC_EDIT, TEXT("0.52"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2b, 5 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCrystalThickness = CreateWindow(WC_EDIT, TEXT("500"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2, 6 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbXstep = CreateWindow(WC_EDIT, TEXT("20"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2b, 6 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    

    maingui.tbBatchDestination = CreateWindow(WC_EDIT, TEXT("20"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2, 10 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBatchDestination2 = CreateWindow(WC_EDIT, TEXT("0"),
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow2+textboxwidth / 2 + 4, 10 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbNumberSims = CreateWindow(WC_EDIT, TEXT("1"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow2, 11 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbNumberSims2 = CreateWindow(WC_EDIT, TEXT("1"),
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow2 + textboxwidth / 2 + 4, 11 * vs, halfBox, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbGPUStatus = CreateWindow(WC_EDIT, TEXT("1"),
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow2 + textboxwidth / 2 + 32, 24 * vs, halfBox-28, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonRun = CreateWindowW(WC_BUTTON, TEXT("Run"), 
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | BS_CENTER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        btnoffset2, 17 * vs, btnwidth, btnHeight, maingui.mainWindow, (HMENU)ID_BTNRUN, hInstance, NULL);
    maingui.buttonStop = CreateWindow(WC_BUTTON, TEXT("Stop"), 
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        btnoffset2, 19 * vs, btnwidth, btnHeight, maingui.mainWindow, (HMENU)ID_BTNSTOP, hInstance, NULL);
    maingui.buttonRefreshDB = CreateWindow(WC_BUTTON, TEXT("Reload DB"), 
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        btnoffset2-btnwidth-8, 19 * vs, btnwidth, btnHeight, maingui.mainWindow, (HMENU)ID_BTNREFRESHDB, hInstance, NULL);
    maingui.buttonRunOnCluster = CreateWindow(WC_BUTTON, TEXT("Cluster file"), 
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        btnoffset2, 20 * vs, btnwidth, btnHeight, maingui.mainWindow, (HMENU)ID_BTNRUNONCLUSTER, hInstance, NULL);
    maingui.pdClusterSelector = CreateWindow(WC_COMBOBOX, TEXT(""), 
        CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, 
        btnoffset2a+1, 20 * vs+1, 205, 9 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonFit = CreateWindow(WC_BUTTON, TEXT("Fit"),
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT,
        btnoffset2, 21 * vs, btnwidth, btnHeight, maingui.mainWindow, (HMENU)ID_BTNFIT, hInstance, NULL);
    
    maingui.pbProgress = CreateWindowEx(0, PROGRESS_CLASS, (LPTSTR)NULL, WS_CHILD | WS_VISIBLE, 
        xOffsetRow2 - 160, 24 * vs+5, 180, 10, maingui.mainWindow, NULL, hInstance, NULL);
    SendMessage(maingui.pbProgress, PBM_SETRANGE, 0, MAKELPARAM(0, 100));
    maingui.pbProgressB = CreateWindowEx(0, PROGRESS_CLASS, (LPTSTR)NULL, WS_CHILD | WS_VISIBLE,
        xOffsetRow2 - 160, 24 * vs+8, 180, 4, maingui.mainWindow, NULL, hInstance, NULL);
    SendMessage(maingui.pbProgressB, PBM_SETRANGE, 0, MAKELPARAM(0, 100));

    TCHAR A[64];
    memset(&A, 0, sizeof(A));
    TCHAR clusterNames[7][64] = { 
        L"Cobra 1xRTX 5000", 
        L"Cobra 2xRTX 5000", 
        L"Cobra 1xV100", 
        L"Cobra 2xV100", 
        L"Raven A100", 
        L"Raven 2xA100", 
        L"Raven 4xA100"};
    for (k = 0; k < 7; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)clusterNames[k]);
        SendMessage(maingui.pdClusterSelector, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    
    SendMessage(maingui.pdClusterSelector, CB_SETCURSEL, (WPARAM)0, 0);
    maingui.buttonPlot = CreateWindow(WC_BUTTON, TEXT("Plot"), 
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        btnoffset2, 18 * vs, btnwidth, btnHeight, maingui.mainWindow, (HMENU)ID_BTNPLOT, hInstance, NULL);
    maingui.tbWhichSimToPlot = CreateWindow(WC_EDIT, TEXT("1"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        btnoffset2 - btnwidth - 8, 18 * vs + 2, 40, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.cbLogPlot = CreateWindow(WC_BUTTON, TEXT(""), WS_CHILD | WS_VISIBLE | BS_AUTOCHECKBOX, 
        btnoffset2a + btnwidth+60, 18 * vs+5, 12, 12, maingui.mainWindow, (HMENU)ID_CBLOGPLOT, hInstance, NULL);
    SendMessage(maingui.cbLogPlot, BM_SETCHECK, BST_CHECKED, 0);
    maingui.trackbarPlot = CreateWindowEx(0, TRACKBAR_CLASS, L"Plot scrubber", WS_CHILD | WS_VISIBLE | TBS_AUTOTICKS | TBS_ENABLESELRANGE,
        btnoffset2a, 18 * vs  + 2, btnwidth, 20, maingui.mainWindow, (HMENU)ID_PLOTSCRUBBER, hInstance, NULL);
    SendMessageW(maingui.trackbarPlot, TBM_SETRANGE, (WPARAM)TRUE, (LPARAM)MAKELONG(0, 1));



    maingui.buttonLoad = CreateWindow(WC_BUTTON, TEXT("Load"), 
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        btnoffset2a, 19 * vs, btnwidth, btnHeight, maingui.mainWindow, (HMENU)ID_BTNLOAD, hInstance, NULL);


    maingui.buttonAddEchoSequence = CreateWindow(WC_BUTTON, TEXT("\x2B83"),
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow2, 12 * vs, btnHeight, btnHeight, maingui.mainWindow, (HMENU)ID_BTNADDECHO, hInstance, NULL);
    maingui.buttonAddCrystalToSequence = CreateWindow(WC_BUTTON, TEXT("\x21AF"),
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow2 + btnHeight + 4, 12 * vs, btnHeight, btnHeight, maingui.mainWindow, (HMENU)ID_BTNADDCRYSTAL, hInstance, NULL);
    maingui.tbSequence = CreateWindow(WC_EDIT, TEXT(""), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL | ES_WANTRETURN, 
        xOffsetRow1 + textboxwidth + 4, 13 * vs, xOffsetRow2-xOffsetRow1, 4*vs-6, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.tbFitting = CreateWindow(WC_EDIT, TEXT(""),
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL | ES_WANTRETURN,
        xOffsetRow1 + textboxwidth + 4, 22 * vs + 2, xOffsetRow2 - xOffsetRow1, 71-vs, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.buttonFile = CreateWindow(WC_BUTTON, TEXT("Set Path"), 
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow3, 0 * vs, btnwidth, btnHeight, maingui.mainWindow, (HMENU)ID_BTNGETFILENAME, hInstance, NULL);

    
    
    maingui.pdPropagationMode = CreateWindow(WC_COMBOBOX, TEXT(""), 
        CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, 
        xOffsetRow2, 7 * vs, textboxwidth, 8 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR energyModeNames[3][64] = {
        TEXT("2D Cartesian"), TEXT("3D radial symmetry"), TEXT("3D")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 3; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)energyModeNames[k]);
        SendMessage(maingui.pdPropagationMode, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdPropagationMode, CB_SETCURSEL, (WPARAM)1, 0);

    maingui.pdBatchMode = CreateWindow(WC_COMBOBOX, TEXT(""), 
        CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, 
        xOffsetRow2, 8 * vs, textboxwidth, 38 * 24, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.pdBatchMode2 = CreateWindow(WC_COMBOBOX, TEXT(""),
        CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
        xOffsetRow2, 9 * vs, textboxwidth, 38 * 24, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR batchModeNames[36][64] = {
        L"none", 
        L"01: Energy 1", 
        L"02: Energy 2", 
        L"03: Frequency 1", 
        L"04: Frequency 2", 
        L"05: Bandwidth 1", 
        L"06: Bandwidth 2", 
        L"07: CEP 1", 
        L"08: CEP 2", 
        L"09: Delay 1", 
        L"10: Delay 2",
        L"11: GDD 1",
        L"12: GDD 2",
        L"13: TOD 1",
        L"14: TOD 2",
        L"15: Thickness 1",
        L"16: Thickness 2",
        L"17: Beamwaist 1",
        L"18: Beamwaist 2",
        L"19: x offset 1",
        L"20: x offset 2",
        L"21: z offset 1",
        L"22: z offset 2",
        L"23: NC angle 1",
        L"24: NC angle 2",
        L"25: Polarization 1",
        L"26: Polarization 2",
        L"27: Circularity 1",
        L"28: Circularity 2",
        L"29: Crystal Theta",
        L"30: Crystal Phi",
        L"31: NL absorption",
        L"32: Gamma",
        L"33: Eff. mass",
        L"34: Thickness",
        L"35: dz",
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 36; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)batchModeNames[k]);
        SendMessage(maingui.pdBatchMode, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 36; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)batchModeNames[k]);
        SendMessage(maingui.pdBatchMode2, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdBatchMode, CB_SETCURSEL, (WPARAM)0, 0);
    SendMessage(maingui.pdBatchMode2, CB_SETCURSEL, (WPARAM)0, 0);
    maingui.pdPulse1Type = CreateWindow(WC_COMBOBOX, TEXT(""), 
        CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, 
        xOffsetRow1, 16 * vs+1, textboxwidth, 9 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR pdPulse1Names[3][64] = {
        TEXT("Synthetic"), 
        TEXT("load FROG"), 
        TEXT("load EOS")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 3; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)pdPulse1Names[k]);
        SendMessage(maingui.pdPulse1Type, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdPulse1Type, CB_SETCURSEL, (WPARAM)0, 0);
    maingui.tbPulse1Path = CreateWindow(WC_EDIT, TEXT("pulse1.speck.dat"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL, 
        0, 17 * vs, xOffsetRow1b+halfBox, 2*vs, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonPulse1Path = CreateWindow(WC_BUTTON, TEXT("Set path"), 
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1 - btnwidth + 8, 16 * vs, btnwidth-8, btnHeight, maingui.mainWindow, (HMENU)ID_BTNPULSE1, hInstance, NULL);

    maingui.pdPulse2Type = CreateWindow(WC_COMBOBOX, TEXT(""), 
        CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, 
        xOffsetRow1, 19 * vs+1, textboxwidth, 9 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR pdPulse2Names[3][64] = {
        TEXT("Synthetic"), 
        TEXT("load FROG"), 
        TEXT("load EOS")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 3; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)pdPulse2Names[k]);
        SendMessage(maingui.pdPulse2Type, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }

    SendMessage(maingui.pdPulse2Type, CB_SETCURSEL, (WPARAM)0, 0);
    maingui.tbPulse2Path = CreateWindow(WC_EDIT, TEXT("pulse2.speck.dat"), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL, 
        0, 20 * vs, xOffsetRow1b + halfBox, 2 * vs, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonPulse2Path = CreateWindow(WC_BUTTON, TEXT("Set path"), 
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, 
        xOffsetRow1 - btnwidth + 8, 19 * vs, btnwidth-8, btnHeight, maingui.mainWindow, (HMENU)ID_BTNPULSE2, hInstance, NULL);

    maingui.pdFittingType = CreateWindow(WC_COMBOBOX, TEXT(""),
        CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
        xOffsetRow1, 22 * vs+1, textboxwidth, 9 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR pdFittingNames[4][64] = {
        TEXT("Maximize x"),
        TEXT("Maximize y"),
        TEXT("Maximize total"),
        TEXT("Match spectrum")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 4; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)pdFittingNames[k]);
        SendMessage(maingui.pdFittingType, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdFittingType, CB_SETCURSEL, (WPARAM)0, 0);
    maingui.tbFittingReferencePath = CreateWindow(WC_EDIT, TEXT("ReferenceFile.txt"),
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL,
        0, 23 * vs, xOffsetRow1b + halfBox, 2 * vs, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonFittingReference = CreateWindow(WC_BUTTON, TEXT("Set path"),
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT,
        xOffsetRow1 - btnwidth + 8, 22 * vs, btnwidth-8, btnHeight, maingui.mainWindow, (HMENU)ID_BTNFITREFERENCE, hInstance, NULL);


    //Text message window
    maingui.textboxSims = CreateWindow(WC_EDIT, TEXT(""), 
        WS_CHILD | WS_VISIBLE | WS_BORDER | ES_LEFT | ES_MULTILINE | WS_VSCROLL | WS_HSCROLL, 
        0, 25 * vs, consoleSize, consoleSize, maingui.mainWindow, NULL, hInstance, NULL);

    if (!maingui.mainWindow)
    {
        return FALSE;
    }

    maingui.tbCPUsims = CreateWindow(WC_EDIT, TEXT("0"),
        WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT,
        btnoffset2a + 84 + 80 + 6, 17 * vs + 2, 40, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.pdPrimaryQueue = CreateWindow(WC_COMBOBOX, TEXT(""),
        CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
        btnoffset2a, 17 * vs, 80, 9 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.pdSecondaryQueue = CreateWindow(WC_COMBOBOX, TEXT(""),
        CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
        btnoffset2a + 84, 17 * vs, 80, 9 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    SetWindowTheme(maingui.mainWindow, L"DarkMode_Explorer", NULL);
    ShowWindow(maingui.mainWindow, nCmdShow);
    SetWindowTheme(maingui.mainWindow, L"DarkMode_Explorer", NULL);
    //make the active set pointer
    activeSetPtr = (simulationParameterSet*)calloc(MAX_SIMULATIONS, sizeof(simulationParameterSet));

    //Find, count, and name the GPUs

    int CUDAdevice, i;
    int CUDAdeviceCount = 0;
    cudaGetDeviceCount(&CUDAdeviceCount);    
    cudaCount = CUDAdeviceCount;
    cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
    struct cudaDeviceProp activeCUDADeviceProp;
    wchar_t wcstring[514];
    size_t convertedChars = 0;
    if (cuErr == cudaSuccess) {
        if (CUDAdeviceCount > 0) {
            hasGPU = TRUE;
        }
        
        printToConsole(maingui.textboxSims, _T("CUDA found %i GPU(s): \r\n"), CUDAdeviceCount);
        for (i = 0; i < CUDAdeviceCount; i++) {
            cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
            memset(wcstring, 0, sizeof(wchar_t));
            mbstowcs_s(&convertedChars, wcstring, 256, activeCUDADeviceProp.name, _TRUNCATE);
            printToConsole(maingui.textboxSims, _T("%ls\r\n"), wcstring);
            printToConsole(maingui.textboxSims, _T(" Memory: %lli MB; Multiprocessors: %i\r\n"), 
                activeCUDADeviceProp.totalGlobalMem/1048576, activeCUDADeviceProp.multiProcessorCount);
        }

    }
    else {
        printToConsole(maingui.textboxSims, L"No CUDA-compatible GPU found.\r\n");
        hasGPU = FALSE;
    }

    //read SYCL devices
    wchar_t syclDeviceList[MAX_LOADSTRING] = { 0 };
    wchar_t syclDefault[MAX_LOADSTRING] = { 0 };
    int syclCount = readSYCLDevices(syclDeviceList, syclDefault);
    //printToConsole(maingui.textboxSims, syclDeviceList);
    printToConsole(maingui.textboxSims, syclDefault);

    if (CUDAdeviceCount > 0) {
        swprintf_s(A, MAX_LOADSTRING, L"CUDA");
        SendMessage(maingui.pdPrimaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
        SendMessage(maingui.pdSecondaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
        memset(&A, 0, sizeof(A));
        for (i = 1; i < CUDAdeviceCount; i++) {
            swprintf_s(A, MAX_LOADSTRING, L"CUDA %i", i);
            SendMessage(maingui.pdPrimaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
            SendMessage(maingui.pdSecondaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
            memset(&A, 0, sizeof(A));
        }
    }
    swprintf_s(A, MAX_LOADSTRING, L"SYCL");
    SendMessage(maingui.pdPrimaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    SendMessage(maingui.pdSecondaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    memset(&A, 0, sizeof(A));
    swprintf_s(A, MAX_LOADSTRING, L"OpenMP");
    SendMessage(maingui.pdPrimaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    SendMessage(maingui.pdSecondaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    memset(&A, 0, sizeof(A));
    SendMessage(maingui.pdPrimaryQueue, (UINT)CB_SETCURSEL, (WPARAM)0, (LPARAM)0);
    SendMessage(maingui.pdSecondaryQueue, (UINT)CB_SETCURSEL, (WPARAM)1, (LPARAM)1);
    //read the crystal database
    crystalDatabasePtr = (crystalEntry*)calloc(MAX_LOADSTRING, sizeof(crystalEntry));
    if (crystalDatabasePtr != NULL) {
        GetCurrentDirectory(MAX_LOADSTRING - 1, programDirectory);
        readCrystalDatabase(crystalDatabasePtr);
        printToConsole(maingui.textboxSims, _T("Read %i entries:\r\n"), (*crystalDatabasePtr).numberOfEntries);
        for (i = 0; i < (*crystalDatabasePtr).numberOfEntries; i++) {
            printToConsole(maingui.textboxSims, _T("Material %i name: %s\r\n"), i, crystalDatabasePtr[i].crystalNameW);
        }
    }



    char defaultFilename[] = "DefaultValues.ini";
    readInputParametersFile(activeSetPtr, crystalDatabasePtr, defaultFilename);
    setInterfaceValuesToActiveValues();
    DWORD hMonitorThread;
    CreateThread(NULL, 0, statusMonitorThread, activeSetPtr, 0, &hMonitorThread);
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
                freeSemipermanentGrids();
            }
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
        case ID_BTNFIT:
            if (!isRunning) {
                isRunning = TRUE;
                mainthread = CreateThread(NULL, 0, fittingThread, activeSetPtr, 0, &hMainThread);
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
        case ID_BTNFITREFERENCE:
            getFileNameBaseFromDlgDat(hWnd, maingui.tbFittingReferencePath);
            break;
        case ID_BTNREFRESHDB:
            SetCurrentDirectory(programDirectory);
            memset(crystalDatabasePtr, 0, 512 * sizeof(crystalEntry));
            readCrystalDatabase(crystalDatabasePtr);
            printToConsole(maingui.textboxSims, _T("Read %i entries:\r\n"), (*crystalDatabasePtr).numberOfEntries);
            for (int i = 0; i < (*crystalDatabasePtr).numberOfEntries; i++) {
                printToConsole(maingui.textboxSims, _T("Material %i name: %s\r\n"), i, crystalDatabasePtr[i].crystalNameW);
            }
            break;
        case ID_BTNPLOT:
            if (isGridAllocated && !isPlotting) {
                plotSim = (int)getDoubleFromHWND(maingui.tbWhichSimToPlot);
                plotSim--;
                plotSim = max(plotSim, 0);
                if (isRunning && (*activeSetPtr).imdone[plotSim] == 0 && !(*activeSetPtr).isInFittingMode) {
                    (*activeSetPtr).imdone[plotSim] = 3;
                    int failCtr = 0;
                    while ((*activeSetPtr).imdone[plotSim] == 3 && failCtr<1000) {
                        failCtr++;
                        Sleep(30);
                    }
                }

                (*activeSetPtr).plotSim = plotSim;
                plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);
            }
            break;
        case ID_BTNLOAD:
            if (openDialogBoxAndLoad(hWnd)) {
                (*activeSetPtr).plotSim = 0;
                
                setInterfaceValuesToActiveValues();
                //plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);
                //Sleep(2000);
                if (isGridAllocated) {
                    setTrackbarLimitsToActiveSet();
                    drawSimPlots(activeSetPtr);
                }
                
            }
            else {
                printToConsole(maingui.textboxSims, L"Read failure.\r\n");
            }
            break;
        case ID_BTNADDCRYSTAL:
            if (GetWindowTextLength(maingui.tbSequence) < 11) {
                SetWindowText(maingui.tbSequence, L"");
            }
            printToConsole(maingui.tbSequence, L"%i %i %2.1lf %2.1lf %2.1e %2.1lf %2.1lf %2.1lf %2.1lf %2.1lf %i;\r\n",
                0, (int)getDoubleFromHWND(maingui.tbMaterialIndex), getDoubleFromHWND(maingui.tbCrystalTheta),
                getDoubleFromHWND(maingui.tbCrystalPhi), getDoubleFromHWND(maingui.tbNonlinearAbsortion),
                getDoubleFromHWND(maingui.tbBandGap),getDoubleFromHWND(maingui.tbDrudeGamma),
                getDoubleFromHWND(maingui.tbEffectiveMass), getDoubleFromHWND(maingui.tbCrystalThickness),
                getDoubleFromHWND(maingui.tbXstep), 0);
            break;
        case ID_BTNADDECHO:
            if (GetWindowTextLength(maingui.tbSequence) < 11) {
                SetWindowText(maingui.tbSequence, L"");
            }
            printToConsole(maingui.tbSequence, L"%i %i %i %i %i %i %i %i %i %i %i;\r\n", 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0);
            break;
        default:
            return DefWindowProc(hWnd, message, wParam, lParam);
        }
    }
    case WM_CTLCOLORBTN:
    {
        HDC hdc = (HDC)wParam;
        SetBkColor(hdc, uiGrey);
        SetTextColor(hdc, uiGrey);
        return (LRESULT)greyBrush;
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
        SetBkColor(hdc, uiGrey);
        SetTextColor(hdc, uiGrey);
        if (HWND(lParam) == maingui.trackbarPlot) {
            return (LRESULT)greyBrush;
        }
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
        GetClientRect(maingui.mainWindow, &mainRect);
        SetWindowPos(maingui.textboxSims, HWND_TOP, 0, 25*maingui.vs, maingui.consoleSize, mainRect.bottom - mainRect.top - 25*maingui.vs, NULL);
        SetWindowPos(maingui.tbFileNameBase, HWND_TOP, maingui.xOffsetRow3, maingui.vs+5, mainRect.right - mainRect.left - maingui.xOffsetRow3- 10, 20, NULL);
 
        int spacerX = maingui.plotSpacerX;
        int spacerY = maingui.plotSpacerY;
        int x0 = maingui.consoleSize + spacerX + 12;
        int y0 = 90 + spacerY;
        int imagePanelSizeX = mainRect.right - mainRect.left - x0 - spacerX;
        int imagePanelSizeY = mainRect.bottom - mainRect.top - y0 - 3 * spacerY;

        int x = x0;
        int y = y0;
        int dx = imagePanelSizeX / 2;
        int dy = imagePanelSizeY / 4;
        int xCorrection = 10;
        int yCorrection = 45;

        SetWindowPos(maingui.plotBox1, HWND_BOTTOM, x - xCorrection, y - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox2, HWND_BOTTOM, x - xCorrection, y + 1 * dy + 1 * spacerY - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox3, HWND_BOTTOM, x - xCorrection, y + 2 * dy + 2 * spacerY - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox4, HWND_BOTTOM, x - xCorrection, y + 3 * dy + 3 * spacerY - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox5, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 0 * dy + 0 * spacerY - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox6, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 1 * dy + 1 * spacerY - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox7, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 2 * dy + 2 * spacerY - yCorrection, dx, dy, NULL);
        SetWindowPos(maingui.plotBox8, HWND_BOTTOM, x + dx + spacerX - xCorrection, y + 3 * dy + 3 * spacerY - yCorrection, dx, dy, NULL);
        
        if (isGridAllocated) {
            ShowWindow(maingui.plotBox1, 0);
            ShowWindow(maingui.plotBox2, 0);
            ShowWindow(maingui.plotBox5, 0);
            ShowWindow(maingui.plotBox6, 0);
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
        SetWindowTheme(maingui.tbFitting, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.tbFittingReferencePath, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonRun, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonFit, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.pdBatchMode, L"DarkMode_CFD", NULL);
        SetWindowTheme(maingui.pdBatchMode2, L"DarkMode_CFD", NULL);
        SetWindowTheme(maingui.pdClusterSelector, L"DarkMode_CFD", NULL);
        SetWindowTheme(maingui.pdFittingType, L"DarkMode_CFD", NULL);
        SetWindowTheme(maingui.pdPropagationMode, L"DarkMode_CFD", NULL);
        SetWindowTheme(maingui.pdPrimaryQueue, L"DarkMode_CFD", NULL);
        SetWindowTheme(maingui.pdSecondaryQueue, L"DarkMode_CFD", NULL);
        SetWindowTheme(maingui.pdPulse1Type, L"DarkMode_CFD", NULL);
        SetWindowTheme(maingui.pdPulse2Type, L"DarkMode_CFD", NULL);
        SetWindowTheme(maingui.buttonFile, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonFittingReference, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonLoad, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonPulse1Path, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonPulse2Path, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonRefreshDB, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonStop, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonRunOnCluster, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonPlot, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonAddCrystalToSequence, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.buttonAddEchoSequence, L"DarkMode_Explorer", NULL);
        SetWindowTheme(maingui.cbLogPlot, L"DarkMode_Explorer", L"wstr");
        //SetWindowTheme(maingui.cbForceCPU, L"DarkMode_Explorer", L"wstr");
        SetWindowTheme(maingui.pbProgress, L"", L"");
        SetWindowTheme(maingui.pbProgressB, L"", L"");
        SendMessage(maingui.pbProgress, PBM_SETBKCOLOR, 0, RGB(0,0,0));
        SendMessage(maingui.pbProgress, PBM_SETBARCOLOR, 0, RGB(87, 254, 255));
        SendMessage(maingui.pbProgressB, PBM_SETBKCOLOR, 0, RGB(0, 0, 0));
        SendMessage(maingui.pbProgressB, PBM_SETBARCOLOR, 0, RGB(255, 84, 255));
        SendMessage(maingui.pdClusterSelector, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
        SendMessage(maingui.pdBatchMode, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
        SendMessage(maingui.pdBatchMode2, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
        SendMessage(maingui.pdFittingType, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
        SendMessage(maingui.pdPropagationMode, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
        SendMessage(maingui.pdPulse1Type, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
        SendMessage(maingui.pdPulse2Type, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
        SendMessage(maingui.pdRStep, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
        SendMessage(maingui.pdPrimaryQueue, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
        SendMessage(maingui.pdSecondaryQueue, CB_SETITEMHEIGHT, -1, (WPARAM)maingui.comboBoxHeight);
		UpdateWindow(hWnd);
		break;
    case WM_HSCROLL:
        setWindowTextToInt(maingui.tbWhichSimToPlot, (int)SendMessage(maingui.trackbarPlot, TBM_GETPOS, 0, 0));
        (*activeSetPtr).plotSim = (int)SendMessage(maingui.trackbarPlot, TBM_GETPOS, 0, 0) - 1;
        plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);
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
    delete[] (*activeSetPtr).ExtOut;
    delete[] (*activeSetPtr).EkwOut;
    delete[] (*activeSetPtr).totalSpectrum;
    return 0;
}

int readParametersFromInterface() {
    std::complex<double> tmp;

    (*activeSetPtr).pulseEnergy1 = getDoubleFromHWND(maingui.tbPulseEnergy1);
    (*activeSetPtr).pulseEnergy2 = getDoubleFromHWND(maingui.tbPulseEnergy2);
    (*activeSetPtr).frequency1 = 1e12 * getDoubleFromHWND(maingui.tbFrequency1);
    (*activeSetPtr).frequency2 = 1e12 * getDoubleFromHWND(maingui.tbFrequency2);
    (*activeSetPtr).bandwidth1 = 1e12 * getDoubleFromHWND(maingui.tbBandwidth1);
    (*activeSetPtr).bandwidth2 = 1e12 * getDoubleFromHWND(maingui.tbBandwidth2);
    (*activeSetPtr).sgOrder1 = 2 * ((int)ceil(getDoubleFromHWND(maingui.tbPulseType1) / 2));
    if ((*activeSetPtr).sgOrder1 < 2) {
        (*activeSetPtr).sgOrder1 = 2;
    }
    (*activeSetPtr).sgOrder2 = 2 * ((int)ceil(getDoubleFromHWND(maingui.tbPulseType2) / 2));
    if ((*activeSetPtr).sgOrder2 < 2) {
        (*activeSetPtr).sgOrder2 = 2;
    }
    (*activeSetPtr).cephase1 = PI * getDoubleFromHWND(maingui.tbCEPhase1);
    (*activeSetPtr).cephase2 = PI * getDoubleFromHWND(maingui.tbCEPhase2);
    (*activeSetPtr).delay1 = 1e-15 * getDoubleFromHWND(maingui.tbPulse1Delay); 
    (*activeSetPtr).delay2 = 1e-15 * getDoubleFromHWND(maingui.tbPulse2Delay);
    (*activeSetPtr).gdd1 = 1e-30 * getDoubleFromHWND(maingui.tbGDD1);
    (*activeSetPtr).gdd2 = 1e-30 * getDoubleFromHWND(maingui.tbGDD2);
    (*activeSetPtr).tod1 = 1e-45 * getDoubleFromHWND(maingui.tbTOD1);
    (*activeSetPtr).tod2 = 1e-45 * getDoubleFromHWND(maingui.tbTOD2);
    (*activeSetPtr).phaseMaterialIndex1 = (int)getDoubleFromHWND(maingui.tbPhaseMaterialIndex1);
    (*activeSetPtr).phaseMaterialIndex2 = (int)getDoubleFromHWND(maingui.tbPhaseMaterialIndex2);
    (*activeSetPtr).phaseMaterialThickness1 = 1e-6*getDoubleFromHWND(maingui.tbPhaseMaterialThickness1);
    (*activeSetPtr).phaseMaterialThickness2 = 1e-6*getDoubleFromHWND(maingui.tbPhaseMaterialThickness2);
    (*activeSetPtr).beamwaist1 = 1e-6 * getDoubleFromHWND(maingui.tbBeamwaist1);
    (*activeSetPtr).beamwaist2 = 1e-6 * getDoubleFromHWND(maingui.tbBeamwaist2);
    tmp = 1e-6 * getDoubleDoublesfromHWND(maingui.tbXoffset1);
    (*activeSetPtr).x01 = real(tmp);
    (*activeSetPtr).y01 = imag(tmp);
    tmp = 1e-6 * getDoubleDoublesfromHWND(maingui.tbXoffset2);
    (*activeSetPtr).x02 = real(tmp);
    (*activeSetPtr).y02 = imag(tmp);
    (*activeSetPtr).z01 = 1e-6 * getDoubleFromHWND(maingui.tbZoffset1);
    (*activeSetPtr).z02 = 1e-6 * getDoubleFromHWND(maingui.tbZoffset2);
    tmp = DEG2RAD * getDoubleDoublesfromHWND(maingui.tbPropagationAngle1);
    (*activeSetPtr).propagationAngle1 = real(tmp);
    (*activeSetPtr).propagationAnglePhi1 = imag(tmp);
    tmp = DEG2RAD * getDoubleDoublesfromHWND(maingui.tbPropagationAngle2);
    (*activeSetPtr).propagationAngle2 = real(tmp);
    (*activeSetPtr).propagationAnglePhi2 = imag(tmp);
    (*activeSetPtr).polarizationAngle1 = DEG2RAD * getDoubleFromHWND(maingui.tbPolarizationAngle1);
    (*activeSetPtr).polarizationAngle2 = DEG2RAD * getDoubleFromHWND(maingui.tbPolarizationAngle2);
    (*activeSetPtr).circularity1 = getDoubleFromHWND(maingui.tbCircularity1);
    (*activeSetPtr).circularity2 = getDoubleFromHWND(maingui.tbCircularity2);

    (*activeSetPtr).materialIndex = (int)getDoubleFromHWND(maingui.tbMaterialIndex);
    (*activeSetPtr).materialIndexAlternate = (int)getDoubleFromHWND(maingui.tbMaterialIndexAlternate);
    (*activeSetPtr).crystalTheta = DEG2RAD * getDoubleFromHWND(maingui.tbCrystalTheta);
    (*activeSetPtr).crystalPhi = DEG2RAD * getDoubleFromHWND(maingui.tbCrystalPhi);

    //(*activeSetPtr).spatialWidth = 1e-6 * getDoubleFromHWND(maingui.tbGridXdim);
    tmp = 1e-6 * getDoubleDoublesfromHWND(maingui.tbGridXdim);
    (*activeSetPtr).spatialWidth = real(tmp);
    (*activeSetPtr).spatialHeight = imag(tmp);
    (*activeSetPtr).rStep = 1e-6 * getDoubleFromHWND(maingui.tbRadialStepSize);
    (*activeSetPtr).timeSpan = 1e-15 * getDoubleFromHWND(maingui.tbTimeSpan);
    (*activeSetPtr).tStep = 1e-15 * getDoubleFromHWND(maingui.tbTimeStepSize);

    (*activeSetPtr).timeSpan = (*activeSetPtr).tStep * (16 * ceil((*activeSetPtr).timeSpan / ((*activeSetPtr).tStep * 16)));

    (*activeSetPtr).crystalThickness = 1e-6 * getDoubleFromHWND(maingui.tbCrystalThickness);
    (*activeSetPtr).propagationStep = 1e-9 * getDoubleFromHWND(maingui.tbXstep);

    (*activeSetPtr).nonlinearAbsorptionStrength = getDoubleFromHWND(maingui.tbNonlinearAbsortion);
    (*activeSetPtr).bandGapElectronVolts = getDoubleFromHWND(maingui.tbBandGap);
    (*activeSetPtr).effectiveMass = getDoubleFromHWND(maingui.tbEffectiveMass);
    (*activeSetPtr).drudeGamma = 1e12 * getDoubleFromHWND(maingui.tbDrudeGamma);

    (*activeSetPtr).runType = 1+(int)SendMessage(maingui.pdClusterSelector, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).batchIndex = (int)SendMessage(maingui.pdBatchMode, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).batchIndex2 = (int)SendMessage(maingui.pdBatchMode2, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).symmetryType = (int)SendMessage(maingui.pdPropagationMode, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);

    (*activeSetPtr).pulse1FileType = (int)SendMessage(maingui.pdPulse1Type, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).pulse2FileType = (int)SendMessage(maingui.pdPulse2Type, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).fittingMode = (int)SendMessage(maingui.pdFittingType, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);

    char noneString[] = "None";

    memset((*activeSetPtr).sequenceString, 0, MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbSequence, (*activeSetPtr).sequenceString, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).sequenceString, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).sequenceString, noneString);
    }
    else {
        removeCharacterFromString((*activeSetPtr).sequenceString, strnlen_s((*activeSetPtr).sequenceString, MAX_LOADSTRING), '\r');
        removeCharacterFromString((*activeSetPtr).sequenceString, strnlen_s((*activeSetPtr).sequenceString, MAX_LOADSTRING), '\n');
    }

    memset((*activeSetPtr).fittingString, 0, MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbFitting, (*activeSetPtr).fittingString, 1024);
    if (strnlen_s((*activeSetPtr).fittingString, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).fittingString, noneString);
    }
    else {
        removeCharacterFromString((*activeSetPtr).fittingString, strnlen_s((*activeSetPtr).fittingString, 1024), '\r');
        removeCharacterFromString((*activeSetPtr).fittingString, strnlen_s((*activeSetPtr).fittingString, 1024), '\n');
    }

    memset((*activeSetPtr).outputBasePath, 0, MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbFileNameBase, (*activeSetPtr).outputBasePath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).outputBasePath, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).outputBasePath, noneString);
    }

    memset((*activeSetPtr).field1FilePath, 0, MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbPulse1Path, (*activeSetPtr).field1FilePath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).field1FilePath, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).field1FilePath, noneString);
    }
    else {
        removeCharacterFromString((*activeSetPtr).field1FilePath, strnlen_s((*activeSetPtr).field1FilePath, MAX_LOADSTRING), '\r');
        removeCharacterFromString((*activeSetPtr).field1FilePath, strnlen_s((*activeSetPtr).field1FilePath, MAX_LOADSTRING), '\n');
    }

    memset((*activeSetPtr).field2FilePath, 0, MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbPulse2Path, (*activeSetPtr).field2FilePath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).field2FilePath, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).field2FilePath, noneString);
    }
    else {
        removeCharacterFromString((*activeSetPtr).field2FilePath, strnlen_s((*activeSetPtr).field2FilePath, MAX_LOADSTRING), '\r');
        removeCharacterFromString((*activeSetPtr).field2FilePath, strnlen_s((*activeSetPtr).field2FilePath, MAX_LOADSTRING), '\n');
    }

    memset((*activeSetPtr).fittingPath, 0, MAX_LOADSTRING * sizeof(char));
    getStringFromHWND(maingui.tbFittingReferencePath, (*activeSetPtr).fittingPath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).fittingPath, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).fittingPath, noneString);
    }
    else {
        removeCharacterFromString((*activeSetPtr).fittingPath, strnlen_s((*activeSetPtr).fittingPath, MAX_LOADSTRING), '\r');
        removeCharacterFromString((*activeSetPtr).fittingPath, strnlen_s((*activeSetPtr).fittingPath, MAX_LOADSTRING), '\n');
    }

    (*activeSetPtr).batchDestination = getDoubleFromHWND(maingui.tbBatchDestination);
    (*activeSetPtr).batchDestination2 = getDoubleFromHWND(maingui.tbBatchDestination2);
    (*activeSetPtr).Nsims = (size_t)getDoubleFromHWND(maingui.tbNumberSims);
    (*activeSetPtr).Nsims2 = (size_t)getDoubleFromHWND(maingui.tbNumberSims2);
    (*activeSetPtr).NsimsCPU = (size_t)getDoubleFromHWND(maingui.tbCPUsims);
    //derived parameters and cleanup:



    (*activeSetPtr).sellmeierType = 0;
    (*activeSetPtr).axesNumber = 0;
    (*activeSetPtr).Ntime = (size_t)(MIN_GRIDDIM*ceil((*activeSetPtr).timeSpan / (MIN_GRIDDIM*(*activeSetPtr).tStep)));
   
    if ((*activeSetPtr).symmetryType == 2) {
        (*activeSetPtr).is3D = TRUE;
        (*activeSetPtr).spatialWidth = (*activeSetPtr).rStep * (MIN_GRIDDIM * ceil((*activeSetPtr).spatialWidth / ((*activeSetPtr).rStep * MIN_GRIDDIM)));
        (*activeSetPtr).Nspace = (size_t)round((*activeSetPtr).spatialWidth / (*activeSetPtr).rStep);
        if ((*activeSetPtr).spatialHeight > 0) {
            (*activeSetPtr).spatialHeight = (*activeSetPtr).rStep * (MIN_GRIDDIM * ceil((*activeSetPtr).spatialHeight / ((*activeSetPtr).rStep * MIN_GRIDDIM)));
        }
        else {
            (*activeSetPtr).spatialHeight = (*activeSetPtr).spatialWidth;
        }
        (*activeSetPtr).Nspace2 = (size_t)round((*activeSetPtr).spatialHeight / (*activeSetPtr).rStep);
    }
    else {
        (*activeSetPtr).is3D = FALSE;
        (*activeSetPtr).Nspace2 = 1;
        (*activeSetPtr).spatialHeight = 0;
        (*activeSetPtr).spatialWidth = (*activeSetPtr).rStep * (MIN_GRIDDIM * ceil((*activeSetPtr).spatialWidth / ((*activeSetPtr).rStep * MIN_GRIDDIM)));
        (*activeSetPtr).Nspace = (size_t)round((*activeSetPtr).spatialWidth / (*activeSetPtr).rStep);
    }
    //printToConsole(maingui.textboxSims, L"(x,y) = %lli, %lli\r\n", (*activeSetPtr).Nspace, (*activeSetPtr).Nspace2);
    (*activeSetPtr).Nfreq = (*activeSetPtr).Ntime / 2 + 1;
    (*activeSetPtr).NgridC = (*activeSetPtr).Nfreq * (*activeSetPtr).Nspace * (*activeSetPtr).Nspace2;
    (*activeSetPtr).Ngrid = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * (*activeSetPtr).Nspace2;
    (*activeSetPtr).kStep = TWOPI / ((*activeSetPtr).Nspace * (*activeSetPtr).rStep);
    (*activeSetPtr).fStep = 1.0 / ((*activeSetPtr).Ntime * (*activeSetPtr).tStep);
    (*activeSetPtr).Npropagation = (size_t)round((*activeSetPtr).crystalThickness / (*activeSetPtr).propagationStep);

    (*activeSetPtr).isCylindric = (*activeSetPtr).symmetryType == 1;
    if ((*activeSetPtr).isCylindric) {
        (*activeSetPtr).x01 = 0;
        (*activeSetPtr).x02 = 0;
        (*activeSetPtr).propagationAngle1 = 0;
        (*activeSetPtr).propagationAngle2 = 0;
    }

    if ((*activeSetPtr).batchIndex == 0 || (*activeSetPtr).Nsims < 1) {
        (*activeSetPtr).Nsims = 1;
    }
    if ((*activeSetPtr).batchIndex2 == 0 || (*activeSetPtr).Nsims2 < 1) {
        (*activeSetPtr).Nsims2 = 1;
    }
    (*activeSetPtr).NsimsCPU = min((*activeSetPtr).NsimsCPU, (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2);

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
    (*activeSetPtr).progressCounter = &progressCounter;
    return 0;
}

//quality of life function - put a text label on a text box window, relative to its position
int labelTextBox(HDC hdc, HWND parentWindow, HWND targetTextBox, const wchar_t* labelText, int xOffset, int yOffset) {
    RECT rectTextBox;
    POINT positionTextBox{};
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

//reads the content of a text box and returns a double containing its numerical value
std::complex<double> getDoubleDoublesfromHWND(HWND inputA)
{
    std::complex<double> sdata = (0,0);
    double* xloc = (double*)&sdata;
    double* yloc = xloc + 1;
    int len = GetWindowTextLength(inputA);
    TCHAR s[256];
    if (len > 0)
    {
        len = GetWindowText(inputA, s, 256);
        swscanf_s(s, L"%lf; %lf", xloc, yloc);
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
int setWindowTextToInt(HWND win, int in) {
    wchar_t textBuffer[128];
    swprintf_s(textBuffer, 128, L"%i", in);
    SetWindowTextW(win, textBuffer);
    return 0;
}

int getNumberOfDecimalsToDisplay(double in, bool isExponential) {
    if (in == 0) return 0;
    in = abs(in);
    int digits = -1;
    int logValue = (int)floor(log10(in));
    in /= pow((double)10.0, logValue);
    while (digits < 15 && in > 1e-3 && in < 9.999) {
        in -= (int)in;
        in *= 10;
        digits++;
    }
    if (isExponential) {
        return max(0, digits);
    }
    return max(0,digits-logValue);
}
int setWindowTextToDouble(HWND win, double in) {
    wchar_t textBuffer[128];
    if (abs(in) < 1e-3) {
        setWindowTextToDoubleExp(win, in);
        return 0;
    }
    int digits = getNumberOfDecimalsToDisplay(in, FALSE);
    if (digits == 0) {
        swprintf_s(textBuffer, 128, L"%i", (int)round(in));
    }
    else if (digits == 1) {
        swprintf_s(textBuffer, 128, L"%2.1lf", in);
    }
    else if (digits == 2) {
        swprintf_s(textBuffer, 128, L"%3.2lf", in);
    }
    else if (digits == 3) {
        swprintf_s(textBuffer, 128, L"%4.3lf", in);
    }
    else {
        swprintf_s(textBuffer, 128, L"%5.4lf", in);
    }
    SetWindowTextW(win, textBuffer);
    return 0;
}

int setWindowTextToDoubleExp(HWND win, double in) {
    wchar_t textBuffer[128];
    int digits = getNumberOfDecimalsToDisplay(in, TRUE);
    if (in == 0) {
        swprintf_s(textBuffer, 128, L"0");
    }
    else if (digits == 0) {
        swprintf_s(textBuffer, 128, L"%4.0e", in);
    }
    else if (digits == 1) {
        swprintf_s(textBuffer, 128, L"%4.1e", in);
    }
    else if (digits == 2) {
        swprintf_s(textBuffer, 128, L"%4.2e", in);
    }
    else if (digits == 3) {
        swprintf_s(textBuffer, 128, L"%4.3e", in);
    }
    else if (digits == 4) {
        swprintf_s(textBuffer, 128, L"%4.4e", in);
    }
    else {
        swprintf_s(textBuffer, 128, L"%e", in);
    }

    SetWindowText(win, textBuffer);
    return 0;
}
int setInterfaceValuesToActiveValues() {

    setWindowTextToDoubleExp(maingui.tbPulseEnergy1, (*activeSetPtr).pulseEnergy1);
    setWindowTextToDoubleExp(maingui.tbPulseEnergy2, (*activeSetPtr).pulseEnergy2);
    setWindowTextToDouble(maingui.tbFrequency1, 1e-12*(*activeSetPtr).frequency1);
    setWindowTextToDouble(maingui.tbFrequency2, 1e-12*(*activeSetPtr).frequency2);
    setWindowTextToDouble(maingui.tbBandwidth1, 1e-12 * (*activeSetPtr).bandwidth1);
    setWindowTextToDouble(maingui.tbBandwidth2, 1e-12 * (*activeSetPtr).bandwidth2);
    setWindowTextToInt(maingui.tbPulseType1, (*activeSetPtr).sgOrder1);
    setWindowTextToInt(maingui.tbPulseType2, (*activeSetPtr).sgOrder2);
    setWindowTextToDouble(maingui.tbCEPhase1, PI * (*activeSetPtr).cephase1);
    setWindowTextToDouble(maingui.tbCEPhase2, PI * (*activeSetPtr).cephase2);
    setWindowTextToDouble(maingui.tbPulse1Delay, 1e15 * (*activeSetPtr).delay1);
    setWindowTextToDouble(maingui.tbPulse2Delay, 1e15 * (*activeSetPtr).delay2);
    setWindowTextToDouble(maingui.tbGDD1, 1e30*(*activeSetPtr).gdd1);
    setWindowTextToDouble(maingui.tbGDD2, 1e30*(*activeSetPtr).gdd2);
    setWindowTextToDouble(maingui.tbTOD1, 1e45*(*activeSetPtr).tod1);
    setWindowTextToDouble(maingui.tbTOD2, 1e45*(*activeSetPtr).tod2);
    setWindowTextToInt(maingui.tbPhaseMaterialIndex1, (*activeSetPtr).phaseMaterialIndex1);
    setWindowTextToInt(maingui.tbPhaseMaterialIndex2, (*activeSetPtr).phaseMaterialIndex2);
    setWindowTextToDouble(maingui.tbPhaseMaterialThickness1, 1e6 * (*activeSetPtr).phaseMaterialThickness1);
    setWindowTextToDouble(maingui.tbPhaseMaterialThickness2, 1e6 * (*activeSetPtr).phaseMaterialThickness2);
    setWindowTextToDouble(maingui.tbBeamwaist1, 1e6 * (*activeSetPtr).beamwaist1);
    setWindowTextToDouble(maingui.tbBeamwaist2, 1e6 * (*activeSetPtr).beamwaist2);
    setWindowTextToDouble(maingui.tbXoffset1, 1e6 * (*activeSetPtr).x01);
    setWindowTextToDouble(maingui.tbXoffset1, 1e6 * (*activeSetPtr).x02);
    setWindowTextToDouble(maingui.tbZoffset1, 1e6 * (*activeSetPtr).z01);
    setWindowTextToDouble(maingui.tbZoffset2, 1e6 * (*activeSetPtr).z02);
    setWindowTextToDouble(maingui.tbPropagationAngle1, RAD2DEG * (*activeSetPtr).propagationAngle1);
    setWindowTextToDouble(maingui.tbPropagationAngle2, RAD2DEG * (*activeSetPtr).propagationAngle2);
    setWindowTextToDouble(maingui.tbPolarizationAngle1, 0.001*round(1000*RAD2DEG * (*activeSetPtr).polarizationAngle1));
    setWindowTextToDouble(maingui.tbPolarizationAngle2, 0.001*round(1000*RAD2DEG * (*activeSetPtr).polarizationAngle2));
    setWindowTextToDouble(maingui.tbCircularity1, (*activeSetPtr).circularity1);
    setWindowTextToDouble(maingui.tbCircularity2, (*activeSetPtr).circularity2);
    SendMessage(maingui.pdPulse1Type, CB_SETCURSEL, (WPARAM)(*activeSetPtr).pulse1FileType, 0);
    SendMessage(maingui.pdPulse2Type, CB_SETCURSEL, (WPARAM)(*activeSetPtr).pulse2FileType, 0);
    setWindowTextToInt(maingui.tbMaterialIndex, (*activeSetPtr).materialIndex);
    setWindowTextToInt(maingui.tbMaterialIndexAlternate, (*activeSetPtr).materialIndexAlternate);
    setWindowTextToDouble(maingui.tbCrystalTheta, RAD2DEG * asin(sin((*activeSetPtr).crystalTheta)));
    setWindowTextToDouble(maingui.tbCrystalPhi, RAD2DEG * asin(sin((*activeSetPtr).crystalPhi)));
    setWindowTextToDoubleExp(maingui.tbNonlinearAbsortion, (*activeSetPtr).nonlinearAbsorptionStrength);
    setWindowTextToDouble(maingui.tbBandGap, (*activeSetPtr).bandGapElectronVolts);
    setWindowTextToDouble(maingui.tbDrudeGamma, 1e-12 * (*activeSetPtr).drudeGamma);

    if (!(*activeSetPtr).is3D) {
        setWindowTextToDouble(maingui.tbGridXdim, 1e6 * (*activeSetPtr).spatialWidth);
    }
    else if((*activeSetPtr).spatialHeight == (*activeSetPtr).spatialWidth
        || (*activeSetPtr).spatialHeight == 1 || (*activeSetPtr).spatialHeight == 0){
        setWindowTextToDouble(maingui.tbGridXdim, 1e6 * (*activeSetPtr).spatialWidth);
    }
    else {
        wchar_t buffer[128];
        swprintf_s(buffer, 128, L"%i;%i", (int)(1e6 * (*activeSetPtr).spatialWidth), (int)(1e6 * (*activeSetPtr).spatialHeight));
        SetWindowTextW(maingui.tbGridXdim, buffer);
    }
    
    setWindowTextToDouble(maingui.tbRadialStepSize, 1e6 * (*activeSetPtr).rStep);
    setWindowTextToDouble(maingui.tbTimeSpan, 1e15 * (*activeSetPtr).timeSpan);
    setWindowTextToDouble(maingui.tbTimeStepSize, 1e15 * (*activeSetPtr).tStep);
    setWindowTextToDouble(maingui.tbCrystalThickness, 1e6 * (*activeSetPtr).crystalThickness);
    setWindowTextToDouble(maingui.tbXstep, 1e9 * (*activeSetPtr).propagationStep);
    SendMessage(maingui.pdPropagationMode, CB_SETCURSEL, (WPARAM)(*activeSetPtr).symmetryType, 0);
    SendMessage(maingui.pdBatchMode, CB_SETCURSEL, (WPARAM)(*activeSetPtr).batchIndex, 0);
    setWindowTextToDouble(maingui.tbBatchDestination, (*activeSetPtr).batchDestination);
    setWindowTextToInt(maingui.tbNumberSims, (int)(*activeSetPtr).Nsims);
    SendMessage(maingui.pdBatchMode2, CB_SETCURSEL, (WPARAM)(*activeSetPtr).batchIndex2, 0);
    setWindowTextToDouble(maingui.tbBatchDestination2, (*activeSetPtr).batchDestination2);
    setWindowTextToInt(maingui.tbNumberSims2, (int)(*activeSetPtr).Nsims2);

    if (strlen((*activeSetPtr).sequenceString)>12) {
        insertLineBreaksAfterSemicolons((*activeSetPtr).sequenceString, MAX_LOADSTRING);
        SetWindowTextA(maingui.tbSequence, (*activeSetPtr).sequenceString);
        removeCharacterFromString((*activeSetPtr).sequenceString, MAX_LOADSTRING, '\r');
        removeCharacterFromString((*activeSetPtr).sequenceString, MAX_LOADSTRING, '\n');
    }
    else {
        SetWindowText(maingui.tbSequence, L"");
    }


        

    if (strcmp((*activeSetPtr).field1FilePath, "None") == 0) {
        SetWindowText(maingui.tbPulse1Path, L"");
    }
    else {
        SetWindowTextA(maingui.tbPulse1Path, (*activeSetPtr).field1FilePath);
    }

    if (strcmp((*activeSetPtr).field2FilePath, "None") == 0) {
        SetWindowText(maingui.tbPulse2Path, L"");
    }
    else {
        SetWindowTextA(maingui.tbPulse2Path, (*activeSetPtr).field2FilePath);
    }

    if (strcmp((*activeSetPtr).fittingPath, "None") == 0) {
        SetWindowText(maingui.tbFittingReferencePath, L"");
    }
    else {
        SetWindowTextA(maingui.tbFittingReferencePath, (*activeSetPtr).fittingPath);
    }
    
    if ((*activeSetPtr).fittingString[0] == 'N') {
        SetWindowText(maingui.tbFitting, L"");
    }
    else {
        insertLineBreaksAfterSemicolons((*activeSetPtr).fittingString, MAX_LOADSTRING);
        SetWindowTextA(maingui.tbFitting, (*activeSetPtr).fittingString);
        removeCharacterFromString((*activeSetPtr).fittingString, MAX_LOADSTRING, '\r');
        removeCharacterFromString((*activeSetPtr).fittingString, MAX_LOADSTRING, '\n');
    }
    SendMessage(maingui.pdFittingType, CB_SETCURSEL, (WPARAM)(*activeSetPtr).fittingMode, 0);

    return 0;
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
    
    labelTextBox(hdc, maingui.mainWindow, maingui.tbMaterialIndex, _T("Material index"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCrystalTheta, _T("Theta, phi (\x00B0)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbNonlinearAbsortion, _T("NL absorption"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbDrudeGamma, _T("Drude: gamma, m"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbBeamwaist1, _T("Beamwaist (\x00B5m)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbGridXdim, _T("Grid width, dx (\x00B5m)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbTimeSpan, _T("Time span, dt (fs)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCrystalThickness, _T("Length, dz (\x00B5m, nm)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdBatchMode, _T("Batch mode"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdBatchMode2, _T("Batch mode 2"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbNumberSims, _T("Batch steps"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbBatchDestination, _T("Batch end"), labos, 0);

    labelTextBox(hdc, maingui.mainWindow, maingui.tbPulse1Delay, _T("Delay (fs)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPulseEnergy1, _T("Energy (J)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbBandwidth1, _T("Bandwidth (THz)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbFrequency1, _T("Frequency (THz)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCEPhase1, _T("CEP/pi"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPulseType1, _T("SG order"), labos, 0);
    
    labelTextBox(hdc, maingui.mainWindow, maingui.tbGDD1, _T("GDD 1 (fs^2)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbTOD1, _T("TOD 1 (fs^3)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPhaseMaterialIndex1, _T("Phase material"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPhaseMaterialThickness1, L"Thickness (\x00B5m)", labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbXoffset1, _T("x offset (\x00B5m)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbZoffset1, _T("z offset (\x00B5m)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPropagationAngle1, _T("NC angle (\x00B0)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbPolarizationAngle1, _T("Polarization (\x00B0)"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbCircularity1, _T("Circularity"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbGPUStatus, _T("GPU (W)"), -72, -1);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdPulse1Type, _T("Pulse 1:"), labos, 4);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdPulse2Type, _T("Pulse 2:"), labos, 4);

    //labelTextBox(hdc, maingui.mainWindow, maingui.tbCPUsims, _T("Offloads:"), -70, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdPropagationMode, _T("Propagation mode"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.pdFittingType, _T("Fit type:"), labos, 0);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbSequence, _T("Crystal sequence:"), 4, -22);
    labelTextBox(hdc, maingui.mainWindow, maingui.tbFitting, _T("Fitting command:"), 4, -22);
    labelTextBox(hdc, maingui.mainWindow, maingui.cbLogPlot, _T("Log"), 16, -1);
    //labelTextBox(hdc, maingui.mainWindow, maingui.cbForceCPU, _T("CPU"), 16, -1);
    //plot labels
    RECT mainRect;
    GetWindowRect(maingui.mainWindow, &mainRect);

    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox1, _T("x-polarization, space/time:"), -maingui.plotSpacerX +8, -22);
    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox2, _T("y-polarization, space/time:"),  -maingui.plotSpacerX +8, -22);
    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox3, _T("x-polarization waveform (GV/m):"),  -maingui.plotSpacerX +8, -22);
    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox4, _T("y-polarization waveform (GV/m):"),  -maingui.plotSpacerX +8, -22);
    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox5, _T("x-polarization, Fourier, Log:"),  -maingui.plotSpacerX +8, -22);
    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox6, _T("y-polarization, Fourier, Log:"),  -maingui.plotSpacerX +8, -22);
    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox7, _T("x-polarization spectrum:"),  -maingui.plotSpacerX +8, -22);
    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox8, _T("y-polarization spectrum:"),  -maingui.plotSpacerX +8, -22);


    GetWindowRect(maingui.plotBox8, &mainRect);
    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox8, _T("Frequency (THz)"), (mainRect.right-mainRect.left)/2 - 15*4, (mainRect.bottom-mainRect.top) + 20);
    GetWindowRect(maingui.plotBox4, &mainRect);
    labelTextBox(hdc, maingui.mainWindow, maingui.plotBox4, _T("Time (fs)"), (mainRect.right - mainRect.left) / 2 - 9 * 4, (mainRect.bottom - mainRect.top) + 20);
    return 0;
}

int getFileNameBaseFromDlg(HWND hWnd, HWND outputTextbox) {
    //create the dialog box and get the file path
    OPENFILENAME ofn;
    TCHAR szFileName[MAX_PATH]{};
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
    if (GetOpenFileNameW(&ofn)) {
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
    TCHAR szFileName[MAX_PATH]{};
    ZeroMemory(&ofn, sizeof(ofn));
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

    if (GetOpenFileNameW(&ofn)) {

        SetWindowText(outputTextbox, szFileName);
    }

    return 0;
}

int openDialogBoxAndLoad(HWND hWnd) {
    WORD fbasedirend;
    OPENFILENAME ofn;
    TCHAR szFileName[MAX_LOADSTRING]{};
    TCHAR szFileNameNoExt[MAX_LOADSTRING]{};
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
    int readParameters;
    wcstombs(fileNameString, szFileName, MAX_LOADSTRING);
    if (GetOpenFileNameW(&ofn)) {
        wcstombs(fileNameString, szFileName, MAX_LOADSTRING);
        if (isGridAllocated) {
            freeSemipermanentGrids();
            isGridAllocated = FALSE;
        }
        readParameters = readInputParametersFile(activeSetPtr, crystalDatabasePtr, fileNameString);
        allocateGrids(activeSetPtr);
        isGridAllocated = TRUE;
        //There should be 50 parameters, remember to update this if adding new ones!
        if (readParameters == 61) {
            //get the base of the file name, so that different files can be made with different extensions based on that
            if (ofn.nFileExtension > 0) {
                fbaseloc = ofn.nFileExtension - 1;
                fbasedirend = ofn.nFileOffset;
            }
            _tcsncpy_s(szFileNameNoExt, szFileName, fbaseloc);
            szFileNameNoExt[MAX_LOADSTRING - 1] = 0;
            wcstombs(fileNameString, szFileNameNoExt, MAX_LOADSTRING);


            int res;
            res = loadSavedFields(activeSetPtr, fileNameString);
            //printToConsole(maingui.textboxSims,L"loaded with %i\r\n", res);
            return TRUE;
        }
        else {
            printToConsole(maingui.textboxSims, L"Read %i\r\n", readParameters);
        }


    }
    
    return FALSE;
}

//Take an array of doubles and make an image in the window
//Color maps:
//  cm = 1: grayscale
//  cm = 2: similar to matlab's jet
//  cm = 3: similar to Colorcet L07
//  cn = 4: vaporwave (symmetric amplitude)
int drawArrayAsBitmap(HWND plotBox, INT64 Nx, INT64 Ny, INT64 x, INT64 y, INT64 height, INT64 width, float* data, int cm) {
    if (Nx * Ny == 0) {
        printToConsole(maingui.textboxSims, L"what are you doing");
    }
    // creating input
    unsigned char* pixels = (unsigned char*)calloc(4 * Nx * Ny, sizeof(unsigned char));
    if (pixels == NULL) {
        return 1;
    }



    size_t Ntot = Nx * Ny;
    float nval;
    int stride = 4;
    //Find the image maximum and minimum
    float imin = data[0];
    float imax = data[0];
    for (size_t i = 1; i < Ntot; i++) {
        if (data[i] > imax) imax = data[i];
        if (data[i] < imin) imin = data[i];
    }

    float oneOver255 = 1.0f/255;
    unsigned char colorMap[256][3];
\
    for (int j = 0; j < 256; j++) {
        switch (cm) {
        case 0:
            
            colorMap[j][0] = j;
            colorMap[j][1] = j;
            colorMap[j][2] = j;
            break;
        case 1:
            nval = j * oneOver255;
            colorMap[j][0] = (unsigned char)(255 * cos(PI * nval / 2.));
            colorMap[j][1] = (unsigned char)(255 * cos(PI * (nval - 0.5)));
            colorMap[j][2] = (unsigned char)(255 * sin(PI * nval / 2.));
            break;
        case 2:
            nval = j * oneOver255;
            colorMap[j][0] = (unsigned char)(255 * cos(PI * nval / 2.));
            colorMap[j][1] = (unsigned char)(255 * cos(PI * (nval - 0.5)));
            colorMap[j][2] = (unsigned char)(255 * sin(PI * nval / 2.));
            if (nval < 0.02) {
                colorMap[j][0] = 255;
                colorMap[j][1] = 128;
                colorMap[j][2] = 128;
            }
            if (nval < 0.01) {
                colorMap[j][0] = 255;
                colorMap[j][1] = 255;
                colorMap[j][2] = 255;
            }
            break;
        case 3:
            nval = 255 * (j*oneOver255);
            colorMap[j][0] = (unsigned char)(255 *
                (0.998 * exp(-pow(7.7469e-03 * (nval - 160), 6))
                    + 0.22 * exp(-pow(0.016818 * (nval - 305), 4))));
            colorMap[j][1] = (unsigned char)(255 *
                (0.022 * exp(-pow(0.042045 * (nval - 25), 4))
                    + 0.11 * exp(-pow(0.015289 * (nval - 120), 4))
                    + 1 * exp(-pow(4.6889e-03 * (nval - 400), 6))));
            colorMap[j][2] = (unsigned char)(255 *
                (exp(-pow(3.1101e-03 * (nval - 415), 10))));
            break;
        case 4:
            nval = j * oneOver255;
            colorMap[j][0] = (unsigned char)(255. * (1.0 * exp(-pow(4.5 * (nval - 0.05), 2))
                + 1.00 * exp(-pow(3.5 * (nval - 1.05), 2))));
            colorMap[j][1] = (unsigned char)(255. * (0.95 * exp(-pow(3.5 * (nval - 1.05), 2))));
            colorMap[j][2] = (unsigned char)(255. * (0.9 * exp(-pow(4.5 * (nval - 0.05), 2))+ 0.2 * exp(-pow(3.5 * (nval - 1.05), 2))));
        }


    }

    if (cm == 4) {
        imax = max(imax, -imin);
        imin = min(imin, -imax);
    }
    unsigned char currentValue;
    if (imin != imax) {
#pragma omp parallel for private(currentValue) num_threads(2)
        for (int p = 0; p < Ntot; p++) {
            currentValue = (unsigned char)(255 * (data[p] - imin) / (imax - imin));
            pixels[stride * p + 0] = colorMap[currentValue][0];
            pixels[stride * p + 1] = colorMap[currentValue][1];
            pixels[stride * p + 2] = colorMap[currentValue][2];
        }
        
    }
    BITMAPINFOHEADER bmih{};

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

    HDC hdc = GetWindowDC(plotBox);
    if (hdc == NULL) {
        printToConsole(maingui.textboxSims,L"no context\r\n");
        return 1;
    }
    HBITMAP hbmp;
    ZeroMemory(&hbmp, sizeof(HBITMAP));
    hbmp = CreateDIBitmap(hdc, &bmih, CBM_INIT, pixels, &dbmi, DIB_RGB_COLORS);
    BITMAP bm;
    ZeroMemory(&bm, sizeof(BITMAP));
    HDC hdcMem = CreateCompatibleDC(hdc);
    if (hdcMem == NULL) {
        printToConsole(maingui.textboxSims, L"no cdc\r\n");
        return 1;
    }
    SelectObject(hdcMem, hbmp);
    GetObject(hbmp, sizeof(bm), &bm);
    //StretchBlt(hdc, (int)x, (int)y, (int)width, (int)height, hdcMem, 0, 0, bm.bmWidth, bm.bmHeight, SRCCOPY);
    if (0 == BitBlt(hdc, (int)x, (int)y, (int)width, (int)height, hdcMem, 0, 0, SRCCOPY)) {
        printToConsole(maingui.textboxSims, L"draw fail\r\n");
    }
    free(pixels);
    DeleteDC(hdcMem);
    DeleteObject(hbmp);
    ReleaseDC(plotBox, hdc);
    return 0;
}

DWORD WINAPI imagePlotThread(LPVOID lpParam) {
    imagePlotStruct *s = (imagePlotStruct*)lpParam;
    RECT plotRect, mainRect;
    
    GetWindowRect((*s).plotBox, &plotRect);
    GetWindowRect(maingui.mainWindow, &mainRect);


    int dx = plotRect.right - plotRect.left;
    int dy = plotRect.bottom - plotRect.top;
    if (dx == 0 || dy == 0) {
        return 1;
    }
    size_t plotSize = (size_t)dx * (size_t)dy;
    float* plotarr2 = (float*)malloc(plotSize * sizeof(float));

    switch ((*s).dataType) {
    case 0:
        //if ((*s).data == (*activeSetPtr).ExtOut)printToConsole(maingui.textboxSims, L"Remapping\r\n");
        linearRemapDoubleToFloat((*s).data, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime, plotarr2, (int)dy, (int)dx);
        break;
    case 1:

        std::complex<double> *shiftedFFT = (std::complex<double>*)malloc((*activeSetPtr).Nspace * (*activeSetPtr).Nfreq * sizeof(std::complex<double>));
        fftshiftD2Z((*s).complexData, shiftedFFT, (*activeSetPtr).Nfreq, (*activeSetPtr).Nspace);
        linearRemapZToLogFloat(shiftedFFT, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Nfreq, plotarr2, (int)dy, (int)dx, (*s).logMin);
        free(shiftedFFT);
        break;
    }
    //if ((*s).data == (*activeSetPtr).ExtOut)printToConsole(maingui.textboxSims, L"Remapped\r\n");
    drawArrayAsBitmap(maingui.mainWindow,
        plotRect.right - plotRect.left, 
        plotRect.bottom - plotRect.top, 
        plotRect.left - mainRect.left, 
        plotRect.top - mainRect.top, 
        plotRect.bottom - plotRect.top, 
        plotRect.right - plotRect.left, 
        plotarr2, (*s).colorMap);
    //if ((*s).data == (*activeSetPtr).ExtOut)printToConsole(maingui.textboxSims, L"Drew\r\n");

    free(plotarr2);

    return 0;
}

DWORD WINAPI drawSimPlots(LPVOID lpParam) {
    if (!isGridAllocated || isPlotting) {
        return 0;
    }
    isPlotting = TRUE;

    imagePlotStruct sTime1, sTime2, sFreq1, sFreq2;
    HANDLE plotThreads[4];
    DWORD plotHandles[4];
	
	bool logPlot = TRUE;
	if (IsDlgButtonChecked(maingui.mainWindow, ID_CBLOGPLOT) != BST_CHECKED) {
		logPlot = FALSE;
	}
	size_t simIndex = (*activeSetPtr).plotSim;
    if (simIndex < 0) {
        simIndex = 0;
    }
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }


	float logPlotOffset = (float)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
	if ((*activeSetPtr).is3D) {
		logPlotOffset = (float)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
	}

	size_t cubeMiddle = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * ((*activeSetPtr).Nspace2 / 2);

    if ((*activeSetPtr).ExtOut == NULL) {
        printToConsole(maingui.textboxSims, L"WTF m8\r\n");
        return 1;
    }
	sTime1.data = 
        &(*activeSetPtr).ExtOut[simIndex * (*activeSetPtr).Ngrid * 2 + cubeMiddle];
	sTime1.plotBox = maingui.plotBox1;
	sTime1.dataType = 0;
	sTime1.colorMap = 4;
	plotThreads[0] = CreateThread(NULL, 0, imagePlotThread, &sTime1, 0, &plotHandles[0]);
    
	sTime2.data = 
        &(*activeSetPtr).ExtOut[(*activeSetPtr).Ngrid + simIndex * (*activeSetPtr).Ngrid * 2 + cubeMiddle];
	sTime2.plotBox = maingui.plotBox2;
	sTime2.dataType = 0;
	sTime2.colorMap = 4;
    plotThreads[1] = CreateThread(NULL, 0, imagePlotThread, &sTime2, 0, &plotHandles[1]);

	sFreq1.complexData = 
        &(*activeSetPtr).EkwOut[simIndex * (*activeSetPtr).NgridC * 2];
	sFreq1.plotBox = maingui.plotBox5;
	sFreq1.dataType = 1;
	sFreq1.logMin = logPlotOffset;
	sFreq1.colorMap = 3;
    plotThreads[2] = CreateThread(NULL, 0, imagePlotThread, &sFreq1, 0, &plotHandles[2]);

	sFreq2.complexData = 
        &(*activeSetPtr).EkwOut[simIndex * (*activeSetPtr).NgridC * 2 + (*activeSetPtr).NgridC];
	sFreq2.plotBox = maingui.plotBox6;
	sFreq2.dataType = 1;
	sFreq2.logMin = logPlotOffset;
	sFreq2.colorMap = 3;
    plotThreads[3] = CreateThread(NULL, 0, imagePlotThread, &sFreq2, 0, &plotHandles[3]);


    plotStruct sWave1, sWave2, sSpectrum1, sSpectrum2;
    sWave1.plotBox = maingui.plotBox3;
    sWave1.dx = (*activeSetPtr).tStep / 1e-15f;
    sWave1.x0 = -(float)(sWave1.dx * (*activeSetPtr).Ntime) / 2;
    sWave1.data = 
        &(*activeSetPtr).ExtOut[simIndex * (*activeSetPtr).Ngrid * 2 + cubeMiddle + (*activeSetPtr).Ntime * (*activeSetPtr).Nspace / 2];
    sWave1.Npts = (*activeSetPtr).Ntime;
    sWave1.unitY = 1e9;
    sWave1.color = D2D1::ColorF(0, 1, 1, 1);
    plotXYDirect2d(&sWave1);

    sWave2.plotBox = maingui.plotBox4;
    sWave2.dx = (*activeSetPtr).tStep / 1e-15f;
    sWave2.x0 = -(float)(sWave2.dx * (*activeSetPtr).Ntime) / 2;
    sWave2.data = 
        &(*activeSetPtr).ExtOut[(*activeSetPtr).Ngrid + simIndex * (*activeSetPtr).Ngrid * 2 + cubeMiddle + (*activeSetPtr).Ntime * (*activeSetPtr).Nspace / 2];
    sWave2.Npts = (*activeSetPtr).Ntime;
    sWave2.unitY = 1e9;
    sWave2.color = D2D1::ColorF(1, 0, 1, 1);
    plotXYDirect2d(&sWave2);

    sSpectrum1.plotBox = maingui.plotBox7;
    sSpectrum1.dx = (float)(*activeSetPtr).fStep / 1e12f;
    sSpectrum1.data = &(*activeSetPtr).totalSpectrum[simIndex * 3 * (*activeSetPtr).Nfreq];
    sSpectrum1.Npts = (*activeSetPtr).Nfreq;
    sSpectrum1.logScale = logPlot;
    sSpectrum1.forceYmin = logPlot;
    sSpectrum1.forcedYmin = -4.0;
    sSpectrum1.color = D2D1::ColorF(0.5, 0, 1, 1);
    plotXYDirect2d(&sSpectrum1);

    sSpectrum2.plotBox = maingui.plotBox8;
    sSpectrum2.dx = (float)(*activeSetPtr).fStep / 1e12f;
    sSpectrum2.data = &(*activeSetPtr).totalSpectrum[(1 + simIndex * 3) * (*activeSetPtr).Nfreq];
    sSpectrum2.Npts = (*activeSetPtr).Nfreq;
    sSpectrum2.logScale = logPlot;
    sSpectrum2.forceYmin = logPlot;
    sSpectrum2.forcedYmin = -4.0;
    sSpectrum2.color = D2D1::ColorF(1, 0, 0.5, 1);
    plotXYDirect2d(&sSpectrum2);


    if (WAIT_TIMEOUT == WaitForMultipleObjects(4, plotThreads, TRUE, 1000)) {
        printToConsole(maingui.textboxSims,L"Warning, an image thread timed out!\r\n");
    }
    for (unsigned int i = 0; i < 4; i++) {
        if (plotThreads[i] != 0)CloseHandle(plotThreads[i]);
    }

    isPlotting = FALSE;
    return 0;
}


int linearRemapZToLogFloat(std::complex<double>* A, int nax, int nay, float* B, int nbx, int nby, double logMin) {
    int i, j;
    float A00;
    float f;

    int nx0, ny0;
    int Ni, Nj;
    for (i = 0; i < nbx; i++) {
        f = i * (nax / (float)nbx);
        Ni = (int)f;
        nx0 = nay * min(Ni, nax);
        for (j = 0; j < nby; j++) {
            f = (j * (nay / (float)nby));
            Nj = (int)f;
            ny0 = min(nay, Nj);
            A00 = (float)log10(cModulusSquared(A[ny0 + nx0]) + logMin);
            B[i * nby + j] = A00;
        }
    }
    return 0;
}
//resize double matrix A to the size of float matrix B
//B is overwritten with the resized matrix
int linearRemapDoubleToFloat(double* A, int nax, int nay, float* B, int nbx, int nby) {
    int i, j;
    int nx0, ny0;
    int Ni, Nj;
    for (i = 0; i < nbx; i++) {
        Ni = (int)(i * (nax / (float)nbx));
        nx0 = nay * min(Ni, nax);
        for (j = 0; j < nby; j++) {
            Nj = (int)((j * (nay / (float)nby)));
            ny0 = min(nay, Nj);
            B[i * nby + j] = (float)A[ny0 + nx0];
        }
    }
    return 0;
}

//resize matrix A to the size of matrix B
//B is overwritten with the resized matrix
int linearRemap(float* A, int nax, int nay, float* B, int nbx, int nby) {
    int j;
    float A00;
    float f;

    int nx0, ny0;
    int Ni, Nj;
#pragma omp parallel for private(j,A00,f,nx0,ny0,Ni,Nj) num_threads(2)
	for (int i = 0; i < nbx; i++) {
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

void setTitleBarDark(HWND hWnd){
    BOOL isDarkMode = TRUE;
    DwmSetWindowAttribute(hWnd, 20, &isDarkMode, sizeof(BOOL));
}


DWORD WINAPI plotXYDirect2d(LPVOID inputStruct) {
    plotStruct* s = (plotStruct*)inputStruct;
    if ((*s).Npts == 0) return 1;
    size_t i;
    D2D1_POINT_2F p1{};
    D2D1_POINT_2F p2{};
    D2D1_ELLIPSE marker{};
    D2D1_SIZE_F sizeF;
    ZeroMemory(&sizeF, sizeof(D2D1_SIZE_F));
    ID2D1HwndRenderTarget *pRenderTarget;
    ID2D1SolidColorBrush* pBrush;
    float markerSize = 1.5;
    float lineWidth = 1.25;

    //get limits of Y
    float maxY = 0;
    float minY = 0;
    if ((*s).logScale) {
        for (i = 0; i < (*s).Npts; i++) {
            maxY = max((float)log10((*s).data[i]), maxY);
            minY = min((float)log10((*s).data[i]), minY);
        }
    }
    else {
        for (i = 0; i < (*s).Npts; i++) {
            maxY = max((float)(*s).data[i], maxY);
            minY = min((float)(*s).data[i], minY);
        }
    }

    if (minY == maxY) {
        minY = -1;
        maxY = 1;
    }
    if ((*s).forceYmin) {
        minY = (*s).forcedYmin;
    }

    //Draw the labels
    int NyTicks = 3;
    wchar_t messageBuffer[MAX_LOADSTRING];
    double yTicks1[3] = { maxY, 0.5 * (maxY + minY), minY };
    double xTicks1[3] = { 0.25 * (*s).dx*(*s).Npts + (*s).x0, 0.5 * (*s).dx * (*s).Npts + (*s).x0, 0.75 * (*s).dx * (*s).Npts + (*s).x0 };
    HDC hdc = GetWindowDC(maingui.mainWindow);
    SetTextColor(hdc, uiGrey);
    SetBkColor(hdc, uiDarkGrey);
    RECT windowRect;
    GetWindowRect((*s).plotBox, &windowRect);

    int posX = windowRect.left;
    int posY = windowRect.top;
    int pixelsTall = windowRect.bottom - windowRect.top;
    int pixelsWide = windowRect.right - windowRect.left;

    GetWindowRect(maingui.mainWindow, &windowRect);
    posX -= windowRect.left;
    posY -= windowRect.top;
    wchar_t blankOut[] = L"          ";
    for (i = 0; i < NyTicks; i++) {
        memset(messageBuffer, 0, MAX_LOADSTRING * sizeof(wchar_t));
        if (abs(yTicks1[i] / (*s).unitY) > 10.0 || abs(yTicks1[i] / (*s).unitY) < 0.01) {
            swprintf_s(messageBuffer, MAX_LOADSTRING,
                _T("%1.1e "), yTicks1[i] / (*s).unitY);
            TextOutW(hdc, posX - 59, posY + (int)(i * 0.96 * pixelsTall / 2), blankOut, (int)_tcslen(blankOut));
            TextOutW(hdc, posX - 59, posY + (int)(i * 0.96 * pixelsTall / 2), messageBuffer, (int)_tcslen(messageBuffer));
        }
        else {
            swprintf_s(messageBuffer, MAX_LOADSTRING,
                _T("%1.4f "), yTicks1[i] / (*s).unitY);
            TextOutW(hdc, posX - 59, posY + (int)(i * 0.96 * pixelsTall / 2), blankOut, (int)_tcslen(blankOut));
            TextOutW(hdc, posX - 54, posY + (int)(i * 0.96 * pixelsTall / 2), messageBuffer, (int)_tcslen(messageBuffer));
        }


    }
    for (i = 0; i < 3; i++) {
        memset(messageBuffer, 0, MAX_LOADSTRING * sizeof(wchar_t));
        swprintf_s(messageBuffer, MAX_LOADSTRING,
            _T("%i"), (int)round(xTicks1[i]));
        TextOutW(hdc, posX + (int)(0.25 * pixelsWide * ((size_t)(i)+1) - 12), posY + pixelsTall, blankOut, (int)_tcslen(blankOut));
        TextOutW(hdc, posX + (int)(0.25 * pixelsWide * ((size_t)(i)+1) - 12), posY + pixelsTall, messageBuffer, (int)_tcslen(messageBuffer));
    }
    ReleaseDC(maingui.mainWindow, hdc);

    RECT targetRectangle;
    GetClientRect((*s).plotBox, &targetRectangle);
    D2D1_SIZE_U size = D2D1::SizeU(targetRectangle.right, targetRectangle.bottom);

    HRESULT hr = maingui.pFactory->CreateHwndRenderTarget(
        D2D1::RenderTargetProperties(),
        D2D1::HwndRenderTargetProperties((*s).plotBox, size),
        &pRenderTarget);
    if (pRenderTarget != NULL) {
        sizeF = pRenderTarget->GetSize();
    }
    
    float scaleX = sizeF.width / ((*s).Npts * (float)(*s).dx);
    float scaleY = sizeF.height / ((float)(maxY - minY));

    if (SUCCEEDED(hr))
    {
        hr = pRenderTarget->CreateSolidColorBrush((*s).color, &pBrush);

        if (SUCCEEDED(hr))
        {
            PAINTSTRUCT ps;
            BeginPaint((*s).plotBox, &ps);
            pRenderTarget->BeginDraw();

            pRenderTarget->Clear(D2D1::ColorF(D2D1::ColorF::Black));
            for (i = 0; i < (*s).Npts - 1; i++) {
                p1.x = scaleX * (i * (float)(*s).dx);
                p2.x = scaleX * ((i + 1) * (float)(*s).dx);
                if ((*s).logScale) {
                    p1.y = sizeF.height - scaleY * ((float)log10((*s).data[i]) - (float)minY);
                    p2.y = sizeF.height - scaleY * ((float)log10((*s).data[i + 1]) - (float)minY);
                }
                else {
                    p1.y = sizeF.height - scaleY * ((float)(*s).data[i] - (float)minY);
                    p2.y = sizeF.height - scaleY * ((float)(*s).data[i + 1] - (float)minY);
                }

                pRenderTarget->DrawLine(p1, p2, pBrush, lineWidth, 0);
                marker.point = p1;
                marker.radiusX = markerSize;
                marker.radiusY = markerSize;
                pRenderTarget->FillEllipse(&marker, pBrush);
            }

            hr = pRenderTarget->EndDraw();
            EndPaint((*s).plotBox, &ps);
        }
        pBrush->Release();
        pRenderTarget->Release();

    }
    return 0;
}

DWORD WINAPI fittingThread(LPVOID lpParam) {
    cancellationCalled = FALSE;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

    readParametersFromInterface();
    (*activeSetPtr).runType = 0;
    if (isGridAllocated) {
        freeSemipermanentGrids();
    }

    allocateGrids(activeSetPtr);
    isGridAllocated = TRUE;
    (*activeSetPtr).isFollowerInSequence = FALSE;
    (*activeSetPtr).crystalDatabase = crystalDatabasePtr;
    loadPulseFiles(activeSetPtr);
    readSequenceString(activeSetPtr);
    configureBatchMode(activeSetPtr);
    readFittingString(activeSetPtr);
    progressCounter = 0;
    (*activeSetPtr).progressCounter = &progressCounter;
    if ((*activeSetPtr).fittingMode == 3) {
        if (loadReferenceSpectrum((*activeSetPtr).fittingPath, activeSetPtr)) {
            printToConsole(maingui.textboxSims, L"Could not read reference file!\r\n");
            free((*activeSetPtr).imdone);
            free((*activeSetPtr).deffTensor);
            free((*activeSetPtr).loadedField1);
            free((*activeSetPtr).loadedField2);
            return 1;
        }

    }

    printToConsole(maingui.textboxSims, L"Fitting %i values in mode %i.\r\nRegion of interest contains %lli elements\r\n", 
        (*activeSetPtr).Nfitting, (*activeSetPtr).fittingMode, (*activeSetPtr).fittingROIsize);

    (*activeSetPtr).runningOnCPU = (!hasGPU || IsDlgButtonChecked(maingui.mainWindow, ID_CBFORCECPU));

    //run the simulations
    if ((*activeSetPtr).runningOnCPU) {
        runDlibFittingCPU(activeSetPtr);
    }
    else {
        runDlibFitting(activeSetPtr);
    }
    

    (*activeSetPtr).plotSim = 0;
    drawSimPlots(activeSetPtr);
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    printToConsole(maingui.textboxSims, _T("Finished fitting after %8.4lf s. \r\n"), 1e-6 *
        (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
    saveDataSet(activeSetPtr, crystalDatabasePtr, (*activeSetPtr).outputBasePath, FALSE);
    setInterfaceValuesToActiveValues();
    printToConsole(maingui.textboxSims, L"Fitting result:\r\n (index, value)\r\n");
    for (int i = 0; i < (*activeSetPtr).Nfitting; i++) {
        printToConsole(maingui.textboxSims, L"%i,  %lf\r\n", i, (*activeSetPtr).fittingResult[i]);
    }
    free((*activeSetPtr).imdone);
    free((*activeSetPtr).deffTensor);
    free((*activeSetPtr).loadedField1);
    free((*activeSetPtr).loadedField2);
    isRunning = FALSE;
    return 0;
}

int insertLineBreaksAfterSemicolons(char* cString, size_t N) {
    size_t i = 0;
    while (i < N-1) {
        if (cString[i] == ';') {
            if (cString[i + 1] != ' ' && cString[i + 2] != ' ') {
                memmove(&cString[i + 3], &cString[i + 1], N - i - 3);
            }
            else if (cString[i + 1] != ' ') {
                memmove(&cString[i + 3], &cString[i + 1], N - i - 3);
            }
            else if (cString[i + 1] == ' ') {
                memmove(&cString[i + 3], &cString[i + 2], N - i - 2);
            }
            cString[i + 1] = '\r';
            cString[i + 2] = '\n';
            i += 2;

            //trim space following semicolon
            if (cString[i + 1] == ' ') {
                memmove(&cString[i + 1], &cString[i + 2], N - i - 1);
            }
        }
        i++;
    }
    return 0;
}

DWORD WINAPI statusMonitorThread(LPVOID lpParam) {

    nvmlDevice_t nvmlDevice = 0;

    unsigned int devicePower = 0;
    int i,j;
    nvmlReturn_t nvmlError;
    size_t lengthEstimate = 0;
    double length;
    double step;

    if (hasGPU) {
        nvmlInit_v2();
        nvmlDeviceGetHandleByIndex_v2(0, &nvmlDevice);
    }

    while (TRUE) {
        
        lengthEstimate = 0;
        if (!(*activeSetPtr).isInFittingMode) {
            //estimate the total steps in the current simulation;
            if ((*activeSetPtr).Nsequence > 0) {
                for (i = 0; i < (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2; i++) {
                    for (j = 0; j < (*activeSetPtr).Nsequence; j++) {
                        if ((int)(*activeSetPtr).sequenceArray[11 * j] == 0) {
                            length = activeSetPtr[i].crystalThickness;
                            step = activeSetPtr[i].propagationStep;
                            if ((int)(*activeSetPtr).sequenceArray[11 * j + 8] != -1) {
                                length = (*activeSetPtr).sequenceArray[11 * j + 8] * 1e-6;
                            }
                            if ((int)(*activeSetPtr).sequenceArray[11 * j + 9] != -1) {
                                step = (*activeSetPtr).sequenceArray[11 * j + 9] * 1e-9;
                            }
                            lengthEstimate += (size_t)round(length / step);
                        }
                    }
                }
            }
            else {
                for (i = 0; i < (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2; i++) {
                    lengthEstimate += (size_t)round(activeSetPtr[i].crystalThickness / activeSetPtr[i].propagationStep);
                }
            }
            if (lengthEstimate > 0) {
                SendMessage(maingui.pbProgress, PBM_SETPOS, (int)((100 * progressCounter) / lengthEstimate), 0);
                SendMessage(maingui.pbProgressB, PBM_SETPOS, (int)((100 * progressCounter) / lengthEstimate), 0);
            }
            
        }
        else {
            lengthEstimate = (*activeSetPtr).fittingMaxIterations;
            if (lengthEstimate > 0) {
                SendMessage(maingui.pbProgress, PBM_SETPOS, (int)((100 * progressCounter) / lengthEstimate), 0);
                SendMessage(maingui.pbProgressB, PBM_SETPOS, (int)((100 * progressCounter) / lengthEstimate), 0);
                //SendMessage(maingui.pbProgress, PBM_SETPOS, (int)(10*progressCounter), 0);
            }
        }

        if (hasGPU) {
            nvmlError = nvmlDeviceGetPowerUsage(nvmlDevice, &devicePower);
            setWindowTextToDouble(maingui.tbGPUStatus, round(devicePower / 1000));

        }
        
        Sleep(500);
    }
    if (hasGPU) {
        nvmlShutdown();
    }
    
    return 0;
}

int setTrackbarLimitsToActiveSet() {
    if (isGridAllocated) {
        SendMessage(maingui.trackbarPlot, TBM_SETRANGE, (WPARAM)TRUE, (LPARAM)MAKELONG(1, (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2));
    }
    return 0;
}

DWORD WINAPI offloadToCPU(LPVOID lpParam) {
    if ((*activeSetPtr).NsimsCPU < 1) return 0;
    simulationParameterSet* cpuSims = (simulationParameterSet*)lpParam;
    int pulldownSelection = (int)SendMessage(maingui.pdSecondaryQueue, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    auto sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
    auto normalFunction = &solveNonlinearWaveEquationCPU;
    int assignedGPU = 0;

    if ( (pulldownSelection - cudaCount) == 0) {
        sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
        normalFunction = &solveNonlinearWaveEquationSYCL;
        if (pulldownSelection == (int)SendMessage(maingui.pdPrimaryQueue, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0)) {
            printToConsole(maingui.textboxSims, L"Warning: can't launch two SYCL simulations simultaneously.\r\n The secondary queue will be done on OpenMP.\r\n");
            sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
            normalFunction = &solveNonlinearWaveEquationCPU;
        }
    }
    else if ((pulldownSelection - cudaCount) == 1) {
        sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
        normalFunction = &solveNonlinearWaveEquationCPU;
    }
    else if ((pulldownSelection - cudaCount) < 0) {
        sequenceFunction = &solveNonlinearWaveEquationSequence;
        normalFunction = &solveNonlinearWaveEquation;
        assignedGPU = pulldownSelection;
    }


    int error = 0;
    if ((*activeSetPtr).isInSequence) {
        for (unsigned int i = 0; i < (*activeSetPtr).NsimsCPU; i++) {
            cpuSims[i].assignedGPU = assignedGPU;
            error = sequenceFunction(&cpuSims[i]);
            if (error) break;
        }
    }
    else {
        for (unsigned int i = 0; i < (*activeSetPtr).NsimsCPU; i++) {
            cpuSims[i].assignedGPU = assignedGPU;
            error = normalFunction(&cpuSims[i]);
            if (error) break;
        }
    }

    if (error) {
        printToConsole(maingui.textboxSims, L"Encountered error %i in CPU run.\r\n", error);
        return 1;
    }

    return 0;
}