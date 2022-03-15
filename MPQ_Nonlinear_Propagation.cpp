// WinTDSE.cpp : Defines the entry point for the application.
//
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS
#include "framework.h"
#include "MPQ_Nonlinear_Propagation.h"
#include "NonlinearPropCUDA.cuh"
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<string.h>
#include<sstream>
#include<complex>
#include<math.h>
#include<Commdlg.h>
#include<stdlib.h>
#include<time.h>
#include<vector>
#include<Windows.h>
#include<CommCtrl.h>

#define MAX_LOADSTRING 1024
#define ID_BTNRUN 11110
#define ID_BTNPLOT 11111
#define ID_BTNGETFILENAME 11112
#define ID_BTNREFRESHDB 11114
#define ID_CBSAVEPSI 11113
#define ID_BTNSTOP 11115
#define ID_BTNPULSE1 11116
#define ID_BTNPULSE2 11117


// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name
struct guiStruct maingui;                       // struct containing all of the GUI elements
struct propthread* activeSetPtr;                // Main structure containing simulation parameters and pointers
struct crystalentry* crystalDatabasePtr;        // Crystal info database
bool isRunning = FALSE;
bool isGridAllocated = FALSE;
bool cancellationCalled = FALSE;


// Forward declarations of (Microsoft) functions included in this code module:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);


int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
    _In_opt_ HINSTANCE hPrevInstance,
    _In_ LPWSTR    lpCmdLine,
    _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Place code here.

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
    wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_MPQNONLINEARPROPAGATION);
    wcex.lpszClassName = szWindowClass;
    wcex.hIconSm = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
    hInst = hInstance; // Store instance handle in our global variable
    int xOffsetRow1 = 160;
    int xOffsetRow2 = 430 + 50;
    int xOffsetRow3 = 540 + 100;
    int vs = 26;
    int radioButtonOffset = xOffsetRow2 - 155 + 24;
    int btnwidth = 100;
    int btnoffset = xOffsetRow1;
    int btnoffset0 = 5;
    int btnoffset2 = xOffsetRow2 + 50;
    int btnoffset2a = btnoffset2 - btnwidth - 10;
    int rbsize = 18;
    int consoleSize = 420;
    int textboxwidth = 150;
    maingui.mainWindow = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW | WS_EX_CONTROLPARENT,
        CW_USEDEFAULT, CW_USEDEFAULT, 2200, 33 * vs + consoleSize + 60, nullptr, nullptr, hInstance, nullptr);
    SetMenu(maingui.mainWindow, NULL);
    SetWindowTextA(maingui.mainWindow, "Nick's nonlinear propagator");


    //text boxes for input parameters
    maingui.tbFieldStrength1 = CreateWindow(TEXT("Edit"), TEXT("1e-4"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 0 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFieldStrength2 = CreateWindow(TEXT("Edit"), TEXT("1e-4"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 1 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFrequency1 = CreateWindow(TEXT("Edit"), TEXT("300"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 2 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFrequency2 = CreateWindow(TEXT("Edit"), TEXT("200"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 3 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBandwidth1 = CreateWindow(TEXT("Edit"), TEXT("80"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 4 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBandwidth2 = CreateWindow(TEXT("Edit"), TEXT("40"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 5 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulseType = CreateWindow(TEXT("Edit"), TEXT("4"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 6 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCEPhase1 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 7 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCEPhase2 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 8 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulse1Delay = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 9 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPulse2Delay = CreateWindowW(TEXT("Edit"), TEXT("-20"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 10 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbGDD1 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 11 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbGDD2 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 12 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTOD1 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 13 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTOD2 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 14 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    
    maingui.tbBeamwaist1 = CreateWindow(TEXT("Edit"), TEXT("400"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 15 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBeamwaist2 = CreateWindow(TEXT("Edit"), TEXT("50"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 16 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbXoffset1 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 17 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbXoffset2 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 18 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbZoffset1 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 19 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbZoffset2 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 20 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPropagationAngle1 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 21 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPropagationAngle2 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 22 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPolarizationAngle1 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 23 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbPolarizationAngle2 = CreateWindow(TEXT("Edit"), TEXT("90"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 24 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCircularity1 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 25 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCircularity2 = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1, 26 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbFileNameBase = CreateWindow(TEXT("Edit"), TEXT("TestFile"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow3, 1 * vs, 775, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbMaterialIndex = CreateWindow(TEXT("Edit"), TEXT("0"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 0 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCrystalTheta = CreateWindow(TEXT("Edit"), TEXT("20.1"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 1 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCrystalPhi = CreateWindow(TEXT("Edit"), TEXT("30"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 2 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbGridXdim = CreateWindow(TEXT("Edit"), TEXT("4096"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 3 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbRadialStepSize = CreateWindow(TEXT("Edit"), TEXT("16"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 4 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTimeSpan = CreateWindow(TEXT("Edit"), TEXT("128"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 5 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbTimeStepSize = CreateWindow(TEXT("Edit"), TEXT("0.5"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 6 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbCrystalThickness = CreateWindow(TEXT("Edit"), TEXT("200"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 7 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbXstep = CreateWindow(TEXT("Edit"), TEXT("250"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 8 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbBatchDestination = CreateWindow(TEXT("Edit"), TEXT("20"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 12 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.tbNumberSims = CreateWindow(TEXT("Edit"), TEXT("41"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow2, 13 * vs, textboxwidth, 20, maingui.mainWindow, NULL, hInstance, NULL);

    maingui.buttonRun = CreateWindow(TEXT("button"), TEXT("Run"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2, 20 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNRUN, hInstance, NULL);
    maingui.buttonStop = CreateWindow(TEXT("button"), TEXT("Stop"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2, 21 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNSTOP, hInstance, NULL);
    maingui.buttonRefreshDB = CreateWindow(TEXT("button"), TEXT("Refresh DB"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, btnoffset2, 22 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNREFRESHDB, hInstance, NULL);
    maingui.tbSequence = CreateWindow(TEXT("Edit"), TEXT(""), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL, xOffsetRow1 + textboxwidth + 4, 15 * vs-2, xOffsetRow2-xOffsetRow1, 66, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonFile = CreateWindow(TEXT("button"), TEXT("Set Path"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow3, 0 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNGETFILENAME, hInstance, NULL);

    int k = 0;
    TCHAR A[64];
    maingui.pdPropagationMode = CreateWindow(WC_COMBOBOX, TEXT(""), CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, xOffsetRow2, 9 * vs, textboxwidth, 5 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR energyModeNames[2][64] = {
        TEXT("2D Cartesian"), TEXT("3D radial symmetry")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 2; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)energyModeNames[k]);
        SendMessage(maingui.pdPropagationMode, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdPropagationMode, CB_SETCURSEL, (WPARAM)0, 0);

    maingui.pdBatchMode = CreateWindow(WC_COMBOBOX, TEXT(""), CBS_DROPDOWNLIST | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE, xOffsetRow2, 11 * vs - 4, textboxwidth, 9 * 20, maingui.mainWindow, NULL, hInstance, NULL);
    TCHAR batchModeNames[7][64] = {
        TEXT("none"), TEXT("Pulse 2 delay"), TEXT("Pulse 1 Energy"), TEXT("CEP"), TEXT("Propagation"), TEXT("Theta"), TEXT("Pulse 1 GDD")
    };
    memset(&A, 0, sizeof(A));
    for (k = 0; k < 7; k++) {
        wcscpy_s(A, sizeof(A) / sizeof(TCHAR), (TCHAR*)batchModeNames[k]);
        SendMessage(maingui.pdBatchMode, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
    }
    SendMessage(maingui.pdBatchMode, CB_SETCURSEL, (WPARAM)1, 0);

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
    maingui.tbPulse1Path = CreateWindow(TEXT("Edit"), TEXT("pulse1.speck.dat"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL, 0, 28 * vs, xOffsetRow2 + 150, 46, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonPulse1Path = CreateWindow(TEXT("button"), TEXT("Set path 1"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1 + textboxwidth +5, 27 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNPULSE1, hInstance, NULL);

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
    maingui.tbPulse2Path = CreateWindow(TEXT("Edit"), TEXT("pulse2.speck.dat"), WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP | WS_EX_CONTROLPARENT | ES_MULTILINE | WS_VSCROLL, 0, 31 * vs, xOffsetRow2 + 150, 46, maingui.mainWindow, NULL, hInstance, NULL);
    maingui.buttonPulse2Path = CreateWindow(TEXT("button"), TEXT("Set path 2"), WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | WS_TABSTOP | WS_EX_CONTROLPARENT, xOffsetRow1 + textboxwidth + 5, 30 * vs, btnwidth, 20, maingui.mainWindow, (HMENU)ID_BTNPULSE2, hInstance, NULL);


    //Text message window
    maingui.textboxSims = CreateWindow(TEXT("Edit"), TEXT(""), WS_CHILD | WS_VISIBLE | WS_BORDER | ES_LEFT | ES_MULTILINE | WS_VSCROLL | WS_HSCROLL, 0, 33 * vs, xOffsetRow2 + 150, consoleSize, maingui.mainWindow, NULL, hInstance, NULL);

    if (!maingui.mainWindow)
    {
        return FALSE;
    }
    ShowWindow(maingui.mainWindow, nCmdShow);

    UpdateWindow(maingui.mainWindow);
    
    
    int CUDAdevice;
    cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
    struct cudaDeviceProp activeCUDADeviceProp;
    cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
    if (cuErr == cudaSuccess) {
        wchar_t* messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
        swprintf_s(messagebuffer, 1024, TEXT("Found GPU: "));
        AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
        size_t origsize = 256 + 1;
        const size_t newsize = 256;
        size_t convertedChars = 0;
        wchar_t wcstring[newsize];
        mbstowcs_s(&convertedChars, wcstring, origsize, activeCUDADeviceProp.name, _TRUNCATE);
        AppendTextToWindow(maingui.textboxSims, wcstring, 256);
        
        swprintf_s(messagebuffer, 1024, TEXT("\r\n\r\n"));
        AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
        free(messagebuffer);
    }
    
    //read the crystal database
    crystalDatabasePtr = (struct crystalentry*)calloc(512, sizeof(struct crystalentry));
    readcrystaldatabase(crystalDatabasePtr, TRUE);
    return TRUE;
}

//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE: Processes messages for the main window.
//
//  WM_COMMAND  - process the application menu
//  WM_PAINT    - Paint the main window
//  WM_DESTROY  - post a quit message and return
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    HANDLE mainthread;
    DWORD hMainThread;

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
                free(activeSetPtr);
            }
            DestroyWindow(hWnd);
            break;
        case ID_BTNRUN:
            if (!isRunning) {
                isRunning = TRUE;
                mainthread = CreateThread(NULL, 0, mainsimthread, activeSetPtr, 0, &hMainThread);
            }
            break;
        case ID_BTNSTOP:
            if (isRunning) {
                cancellationCalled = TRUE;
                for (int i = 0; i < (*activeSetPtr).NSims; i++) {
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
            memset(crystalDatabasePtr, 0, 512 * sizeof(struct crystalentry));
            readcrystaldatabase(crystalDatabasePtr, TRUE);
            break;

        default:
            return DefWindowProc(hWnd, message, wParam, lParam);
        }
    }
    break;
    case WM_PAINT:
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hWnd, &ps);
        DrawLabels(hdc);
        EndPaint(hWnd, &ps);
    }
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

DWORD WINAPI mainsimthread(LPVOID lpParam) {
    cancellationCalled = FALSE;
    wchar_t* messagebuffer;
    int j, k;
    time_t tstart, tthreadmid;
    time(&tstart);
    HANDLE plotThread;
    DWORD hplotThread;

    activeSetPtr = (struct propthread*)malloc(1024 * sizeof(struct propthread));

    readParametersFromInterfaceAndAllocate();
    
    //run the simulations
    for (j = 0; j < (*activeSetPtr).NSims; j++) {
        if ((*activeSetPtr).isInSequence) {
            for (k = 0; k < (*activeSetPtr).Nsequence; k++) {
                resolvesequence(k, &activeSetPtr[j]);
                propagationLoop(&activeSetPtr[j]);
                (*activeSetPtr).plotSim = j;
                plotThread = CreateThread(NULL, 0, drawsimplots, activeSetPtr, 0, &hplotThread);
                if (activeSetPtr[j].memoryError > 0) {
                    messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
                    swprintf_s(messagebuffer, 1024,
                        _T("Warning: device memory error (%i).\r\n"), activeSetPtr[j].memoryError);
                    AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
                    free(messagebuffer);
                }
            }
        }
        else {
            propagationLoop(&activeSetPtr[j]);
            (*activeSetPtr).plotSim = j;
            plotThread = CreateThread(NULL, 0, drawsimplots, activeSetPtr, 0, &hplotThread);
        }

        if (cancellationCalled) {
            messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
            swprintf_s(messagebuffer, 1024,
                _T("Warning: series cancelled, stopping after %i simulations.\r\n"), j+1);
            AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
            free(messagebuffer);
            break;
        }
    }
    time(&tthreadmid);

    messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
    swprintf_s(messagebuffer, 1024,
        _T("Finished after %i s. \r\n"), (int)(tthreadmid - tstart));
    AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
    free(messagebuffer);

    saveDataSet();

    free((*activeSetPtr).sequenceString);
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
int readParametersFromInterfaceAndAllocate() {
    wchar_t* messagebuffer;
    if (isGridAllocated) {
        isGridAllocated = FALSE;
        free((*activeSetPtr).ExtOut);
        free((*activeSetPtr).EkwOut);
        free((*activeSetPtr).Ext);
        free((*activeSetPtr).Ekw);
        free(activeSetPtr);
    }
    int j;
    const double pi = 3.1415926535897932384626433832795;

    (*activeSetPtr).materialIndex = (int)HWNDToDouble(maingui.tbMaterialIndex);;
    (*activeSetPtr).crystalTheta = (pi / 180) * HWNDToDouble(maingui.tbCrystalTheta);
    (*activeSetPtr).crystalPhi = (pi / 180) * HWNDToDouble(maingui.tbCrystalPhi);
    (*activeSetPtr).crystalThickness = 1e-6 * HWNDToDouble(maingui.tbCrystalThickness);
    (*activeSetPtr).sellmeierType = 0;
    (*activeSetPtr).axesNumber = 0;

    (*activeSetPtr).beamwaist1 = 1e-6 * HWNDToDouble(maingui.tbBeamwaist1);
    (*activeSetPtr).beamwaist2 = 1e-6 * HWNDToDouble(maingui.tbBeamwaist2);
    (*activeSetPtr).x01 = 1e-6 * HWNDToDouble(maingui.tbXoffset1);
    (*activeSetPtr).x02 = 1e-6 * HWNDToDouble(maingui.tbXoffset2);
    (*activeSetPtr).z01 = 1e-6 * HWNDToDouble(maingui.tbZoffset1);
    (*activeSetPtr).z02 = 1e-6 * HWNDToDouble(maingui.tbZoffset2);
    (*activeSetPtr).propagationAngle1 = (pi / 180) * HWNDToDouble(maingui.tbPropagationAngle1);
    (*activeSetPtr).propagationAngle2 = (pi / 180) * HWNDToDouble(maingui.tbPropagationAngle2);


    (*activeSetPtr).spatialWidth = 1e-6 * HWNDToDouble(maingui.tbGridXdim);
    (*activeSetPtr).rStep = 1e-6 * HWNDToDouble(maingui.tbRadialStepSize);
    (*activeSetPtr).Nspace = (long long)round((*activeSetPtr).spatialWidth / (*activeSetPtr).rStep);
    (*activeSetPtr).timeSpan = 1e-15 * HWNDToDouble(maingui.tbTimeSpan);
    (*activeSetPtr).tStep = 1e-15 * HWNDToDouble(maingui.tbTimeStepSize);
    (*activeSetPtr).Ntime = (long long)round((*activeSetPtr).timeSpan / (*activeSetPtr).tStep);
    (*activeSetPtr).Ngrid = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace;
    (*activeSetPtr).kStep = 2 * pi / ((*activeSetPtr).Nspace * (*activeSetPtr).rStep);
    (*activeSetPtr).fStep = 1.0 / ((*activeSetPtr).Ntime * (*activeSetPtr).tStep);
    (*activeSetPtr).propagationStep = 1e-9 * HWNDToDouble(maingui.tbXstep);
    (*activeSetPtr).Npropagation = (long long)round((*activeSetPtr).crystalThickness / (*activeSetPtr).propagationStep);
    (*activeSetPtr).pulseEnergy1 = HWNDToDouble(maingui.tbFieldStrength1);
    (*activeSetPtr).pulseEnergy2 = HWNDToDouble(maingui.tbFieldStrength2);
    (*activeSetPtr).frequency1 = 1e12 * HWNDToDouble(maingui.tbFrequency1);
    (*activeSetPtr).frequency2 = 1e12 * HWNDToDouble(maingui.tbFrequency2);
    (*activeSetPtr).bandwidth1 = 1e12 * HWNDToDouble(maingui.tbBandwidth1);
    (*activeSetPtr).bandwidth2 = 1e12 * HWNDToDouble(maingui.tbBandwidth2);
    (*activeSetPtr).cephase1 = HWNDToDouble(maingui.tbCEPhase1);
    (*activeSetPtr).cephase2 = HWNDToDouble(maingui.tbCEPhase2);
    (*activeSetPtr).delay1 = -1e-15 * HWNDToDouble(maingui.tbPulse1Delay) + (*activeSetPtr).timeSpan / 2;
    (*activeSetPtr).delay2 = -1e-15 * HWNDToDouble(maingui.tbPulse2Delay) + (*activeSetPtr).timeSpan / 2;
    (*activeSetPtr).gdd1 = 1e-30 * HWNDToDouble(maingui.tbGDD1);
    (*activeSetPtr).gdd2 = 1e-30 * HWNDToDouble(maingui.tbGDD2);
    (*activeSetPtr).tod1 = 1e-45 * HWNDToDouble(maingui.tbTOD1);
    (*activeSetPtr).tod2 = 1e-45 * HWNDToDouble(maingui.tbTOD2);
    (*activeSetPtr).sgOrder1 = 2 * ((int)ceil(HWNDToDouble(maingui.tbPulseType) / 2));
    if ((*activeSetPtr).sgOrder1 < 2) {
        (*activeSetPtr).sgOrder1 = 2;
    }
    (*activeSetPtr).sgOrder2 = (*activeSetPtr).sgOrder1;
    (*activeSetPtr).polarizationAngle1 = (pi / 180) * HWNDToDouble(maingui.tbPolarizationAngle1);
    (*activeSetPtr).polarizationAngle2 = (pi / 180) * HWNDToDouble(maingui.tbPolarizationAngle2);
    (*activeSetPtr).circularity1 = HWNDToDouble(maingui.tbCircularity1);
    (*activeSetPtr).circularity2 = HWNDToDouble(maingui.tbCircularity2);
    (*activeSetPtr).batchIndex = (int)SendMessage(maingui.pdBatchMode, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).symmetryType = (int)SendMessage(maingui.pdPropagationMode, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    (*activeSetPtr).isCylindric = (*activeSetPtr).symmetryType == 1;
    if ((*activeSetPtr).isCylindric) {
        (*activeSetPtr).x01 = 0;
        (*activeSetPtr).x02 = 0;
        (*activeSetPtr).propagationAngle1 = 0;
        (*activeSetPtr).propagationAngle2 = 0;
    }
    (*activeSetPtr).batchDestination = HWNDToDouble(maingui.tbBatchDestination);

    (*activeSetPtr).NSims = (int)HWNDToDouble(maingui.tbNumberSims);
    if ((*activeSetPtr).batchIndex == 0 || (*activeSetPtr).batchIndex == 4 || (*activeSetPtr).NSims < 1) {
        (*activeSetPtr).NSims = 1;
    }


    (*activeSetPtr).loadedField1 = (std::complex<double>*)calloc((*activeSetPtr).Ntime, sizeof(std::complex<double>));
    (*activeSetPtr).loadedField2 = (std::complex<double>*)calloc((*activeSetPtr).Ntime, sizeof(std::complex<double>));
    int pulse1FileType = (int)SendMessage(maingui.pdPulse1Type, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
    int pulse2FileType = (int)SendMessage(maingui.pdPulse2Type, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);

    (*activeSetPtr).field1IsAllocated = FALSE;
    (*activeSetPtr).field2IsAllocated = FALSE;

    int frogLines = 0;
    if (pulse1FileType == 1) {
        char pulse1Path[MAX_LOADSTRING];
        HWNDToString(maingui.tbPulse1Path, pulse1Path, MAX_LOADSTRING);


        frogLines = loadfrogspeck(pulse1Path, (*activeSetPtr).loadedField1, (*activeSetPtr).Ntime, (*activeSetPtr).fStep, 0.0, 1);
        if (frogLines > 0) (*activeSetPtr).field1IsAllocated = TRUE;
        messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
        swprintf_s(messagebuffer, 1024,
            _T("loaded FROG file 1 (%i lines, %i).\r\n"), frogLines, (*activeSetPtr).field1IsAllocated);
        AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
        free(messagebuffer);

    }
    if (pulse2FileType == 1) {
        char pulse2Path[MAX_LOADSTRING];
        HWNDToString(maingui.tbPulse2Path, pulse2Path, MAX_LOADSTRING);

        frogLines = loadfrogspeck(pulse2Path, (*activeSetPtr).loadedField2, (*activeSetPtr).Ntime, (*activeSetPtr).fStep, 0.0, 1);
        if (frogLines > 0) (*activeSetPtr).field2IsAllocated = TRUE;
        messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
        swprintf_s(messagebuffer, 1024,
            _T("loaded FROG file 2 (%i lines).\r\n"), frogLines);
        AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
        free(messagebuffer);
    }

    (*activeSetPtr).sequenceString = (char*)calloc(256 * MAX_LOADSTRING, sizeof(char));
    (*activeSetPtr).sequenceArray = (double*)calloc(256 * MAX_LOADSTRING, sizeof(double));
    HWNDToString(maingui.tbSequence, (*activeSetPtr).sequenceString, MAX_LOADSTRING * 256);

    char* tokToken = strtok((*activeSetPtr).sequenceString, ";");
    int sequenceCount = sscanf((*activeSetPtr).sequenceString, "%lf %lf %lf %lf %lf %lf", (*activeSetPtr).sequenceArray, &(*activeSetPtr).sequenceArray[1], &(*activeSetPtr).sequenceArray[2], &(*activeSetPtr).sequenceArray[3], &(*activeSetPtr).sequenceArray[4], &(*activeSetPtr).sequenceArray[5]);

    tokToken = strtok(NULL, ";");
    int lastread = sequenceCount;
    while (tokToken != NULL && lastread == 6) {
        lastread = sscanf(tokToken, "%lf %lf %lf %lf %lf %lf", &(*activeSetPtr).sequenceArray[sequenceCount], &(*activeSetPtr).sequenceArray[sequenceCount + 1], &(*activeSetPtr).sequenceArray[sequenceCount + 2], &(*activeSetPtr).sequenceArray[sequenceCount + 3], &(*activeSetPtr).sequenceArray[sequenceCount + 4], &(*activeSetPtr).sequenceArray[sequenceCount + 5]);
        sequenceCount += lastread;
        tokToken = strtok(NULL, ";");
    }

    (*activeSetPtr).Nsequence = sequenceCount / 6;
    (*activeSetPtr).isInSequence = ((*activeSetPtr).Nsequence > 0);


    (*activeSetPtr).Ext = (std::complex<double>*)calloc((*activeSetPtr).Ngrid * 2 * (*activeSetPtr).NSims, sizeof(std::complex<double>));
    (*activeSetPtr).Ekw = (std::complex<double>*)calloc((*activeSetPtr).Ngrid * 2 * (*activeSetPtr).NSims, sizeof(std::complex<double>));

    (*activeSetPtr).ExtOut = (std::complex<double>*)calloc((*activeSetPtr).Ngrid * 2 * (*activeSetPtr).NSims, sizeof(std::complex<double>));
    (*activeSetPtr).EkwOut = (std::complex<double>*)calloc((*activeSetPtr).Ngrid * 2 * (*activeSetPtr).NSims, sizeof(std::complex<double>));

    isGridAllocated = TRUE;
    (*activeSetPtr).refractiveIndex1 = (std::complex<double>*)calloc((*activeSetPtr).Ngrid * (*activeSetPtr).NSims, sizeof(std::complex<double>));
    (*activeSetPtr).refractiveIndex2 = (std::complex<double>*)calloc((*activeSetPtr).Ngrid * (*activeSetPtr).NSims, sizeof(std::complex<double>));
    (*activeSetPtr).deffTensor = (double*)calloc(9 * (*activeSetPtr).NSims, sizeof(double));
    (*activeSetPtr).imdone = (int*)calloc((*activeSetPtr).NSims, sizeof(int));


    (*activeSetPtr).chi2Tensor = crystalDatabasePtr[(*activeSetPtr).materialIndex].d;
    (*activeSetPtr).chi3Tensor = crystalDatabasePtr[(*activeSetPtr).materialIndex].chi3;
    (*activeSetPtr).nonlinearSwitches = crystalDatabasePtr[(*activeSetPtr).materialIndex].nonlinearSwitches;
    (*activeSetPtr).absorptionParameters = crystalDatabasePtr[(*activeSetPtr).materialIndex].absorptionParameters;
    (*activeSetPtr).sellmeierCoefficients = crystalDatabasePtr[(*activeSetPtr).materialIndex].sellmeierCoefficients;
    (*activeSetPtr).sellmeierType = crystalDatabasePtr[(*activeSetPtr).materialIndex].sellmeierType;
    (*activeSetPtr).axesNumber = crystalDatabasePtr[(*activeSetPtr).materialIndex].axisType;


    double batchstart = 0;
    if ((*activeSetPtr).batchIndex == 1) {
        batchstart = (*activeSetPtr).delay2 - (*activeSetPtr).timeSpan / 2;
    }
    if ((*activeSetPtr).batchIndex == 2) {
        batchstart = (*activeSetPtr).pulseEnergy1;
    }
    if ((*activeSetPtr).batchIndex == 3) {
        batchstart = (*activeSetPtr).cephase1;
    }
    if ((*activeSetPtr).batchIndex == 5) {
        batchstart = (*activeSetPtr).crystalTheta;
    }
    if ((*activeSetPtr).batchIndex == 6) {
        batchstart = (*activeSetPtr).gdd1;
    }

    for (j = 0; j < (*activeSetPtr).NSims; j++) {
        if (j > 0) {
            memcpy(&activeSetPtr[j], activeSetPtr, sizeof(struct propthread));
        }


        activeSetPtr[j].deffTensor = &(*activeSetPtr).deffTensor[9 * j];
        activeSetPtr[j].Ext = &(*activeSetPtr).Ext[j * (*activeSetPtr).Ngrid * 2];
        activeSetPtr[j].Ekw = &(*activeSetPtr).Ekw[j * (*activeSetPtr).Ngrid * 2];
        activeSetPtr[j].ExtOut = &(*activeSetPtr).ExtOut[j * (*activeSetPtr).Ngrid * 2];
        activeSetPtr[j].EkwOut = &(*activeSetPtr).EkwOut[j * (*activeSetPtr).Ngrid * 2];


        activeSetPtr[j].isFollowerInSequence = FALSE;

        if ((*activeSetPtr).batchIndex == 1) {
            activeSetPtr[j].delay2 += j * ((-1e-15 * (*activeSetPtr).batchDestination) - batchstart) / ((*activeSetPtr).NSims - 1);
        }
        if ((*activeSetPtr).batchIndex == 2) {
            activeSetPtr[j].pulseEnergy1 += j * ((*activeSetPtr).batchDestination - batchstart) / ((*activeSetPtr).NSims - 1);
        }
        if ((*activeSetPtr).batchIndex == 3) {
            activeSetPtr[j].cephase1 += j * (pi * (*activeSetPtr).batchDestination - batchstart) / ((*activeSetPtr).NSims - 1);
        }
        if ((*activeSetPtr).batchIndex == 5) {
            activeSetPtr[j].crystalTheta += j * ((pi / 180) * (*activeSetPtr).batchDestination - batchstart) / ((*activeSetPtr).NSims - 1);
        }
        if ((*activeSetPtr).batchIndex == 6) {
            activeSetPtr[j].gdd1 += j * (1e-30 * (*activeSetPtr).batchDestination - batchstart) / ((*activeSetPtr).NSims - 1);
        }
    }
    return 0;
}
int saveDataSet() {
    int j, k;

    //Save the results as double instead of complex
    double* saveEout = (double*)&(*activeSetPtr).refractiveIndex1[0];
    double* saveEin = (double*)&(*activeSetPtr).refractiveIndex2[0];
    for (j = 0; j < ((*activeSetPtr).Ngrid * (*activeSetPtr).NSims * 2); j++) {
        saveEout[j] = real((*activeSetPtr).ExtOut[j]);
        saveEin[j] = real((*activeSetPtr).Ext[j]);
    }

    char outputbase[MAX_LOADSTRING];
    HWNDToString(maingui.tbFileNameBase, outputbase, MAX_LOADSTRING);

    FILE* textfile;
    char* outputpath = (char*)calloc(MAX_LOADSTRING, sizeof(char));
    strcpy(outputpath, outputbase);
    strcat(outputpath, ".txt");
    textfile = fopen(outputpath, "w");
    fprintf(textfile, "Amplitude 1: %e\nAmplitude 2: %e\nFrequency 1: %e\nFrequency 2: %e\nBandwidth 1: %e\nBandwidth 2: %e\n", (*activeSetPtr).pulseEnergy1, (*activeSetPtr).pulseEnergy2, (*activeSetPtr).frequency1, (*activeSetPtr).frequency2, (*activeSetPtr).bandwidth1, (*activeSetPtr).bandwidth2);
    fprintf(textfile, "SG order: %i\nCEP 1: %e\nCEP 2: %e\nDelay 1: %e\nDelay 2: %e\nGDD 1: %e\nGDD 2: %e\nTOD 1: %e\nTOD 2: %e\n", (*activeSetPtr).sgOrder1, (*activeSetPtr).cephase1, (*activeSetPtr).cephase2, (*activeSetPtr).delay1, (*activeSetPtr).delay2, (*activeSetPtr).gdd1, (*activeSetPtr).gdd2, (*activeSetPtr).tod1, (*activeSetPtr).tod2);
    fprintf(textfile, "Beamwaist 1: %e\nBeamwaist 2: %e\nx offset 1: %e\nx offset 2: %e\nz offset 1: %e\nz offset 2: %e\nNC angle 1: %e\nNC angle 2: %e\n", (*activeSetPtr).beamwaist1, (*activeSetPtr).beamwaist2, (*activeSetPtr).x01, (*activeSetPtr).x02, (*activeSetPtr).z01, (*activeSetPtr).z02, (*activeSetPtr).propagationAngle1, (*activeSetPtr).propagationAngle2);
    fprintf(textfile, "Polarization 1: %e\nPolarization 2: %e\nCircularity 1: %e\nCircularity 2: %e\n", (*activeSetPtr).polarizationAngle1, (*activeSetPtr).polarizationAngle2, (*activeSetPtr).circularity1, (*activeSetPtr).circularity2);
    fprintf(textfile, "Material index: %i\n", (*activeSetPtr).materialIndex);
    fwprintf(textfile, _T("Material name: %s\nSellmeier reference: %s\nChi2 reference: %s\nChi3 reference: %s\n"), crystalDatabasePtr[(*activeSetPtr).materialIndex].crystalNameW, crystalDatabasePtr[(*activeSetPtr).materialIndex].sellmeierReference, crystalDatabasePtr[(*activeSetPtr).materialIndex].dReference, crystalDatabasePtr[(*activeSetPtr).materialIndex].chi3Reference);
    fprintf(textfile, "Sellmeier coefficients: \n");
    for (j = 0; j < 3; j++) {
        for (k = 0; k < 22; k++) {
            fprintf(textfile, "%e ", crystalDatabasePtr[(*activeSetPtr).materialIndex].sellmeierCoefficients[j * 22 + k]);
        }
        fprintf(textfile, "\n");
    }
    fprintf(textfile, "Crystal theta: %e\nCrystal phi: %e\nGrid width: %e\ndx: %e\nTime span: %e\ndt: %e\nThickness: %e\ndz: %e\n", (*activeSetPtr).crystalTheta, (*activeSetPtr).crystalPhi, (*activeSetPtr).spatialWidth, (*activeSetPtr).rStep, (*activeSetPtr).timeSpan, (*activeSetPtr).tStep, (*activeSetPtr).crystalThickness, (*activeSetPtr).propagationStep);
    fprintf(textfile, "Propagation mode: %i\n", (*activeSetPtr).symmetryType);
    fprintf(textfile, "Batch mode: %i\nBatch destination: %e\nBatch steps: %i\n", (*activeSetPtr).batchIndex, (*activeSetPtr).batchDestination, (*activeSetPtr).NSims);
    if ((*activeSetPtr).isInSequence) {
        HWNDToString(maingui.tbSequence, (*activeSetPtr).sequenceString, MAX_LOADSTRING * 256);
        fprintf(textfile, "Sequence: %s\n", (*activeSetPtr).sequenceString);
    }

    fprintf(textfile, "Code version: 0.00 Feb. 15, 2022\n");

    fclose(textfile);

    FILE* ExtOutFile;
    strcpy(outputpath, outputbase);
    strcat(outputpath, "_ExtOut.dat");
    ExtOutFile = fopen(outputpath, "wb");
    fwrite(saveEout, sizeof(double), 2 * ((*activeSetPtr).Ngrid * (*activeSetPtr).NSims) + 1024, ExtOutFile);
    fclose(ExtOutFile);

    FILE* ExtInFile;
    strcpy(outputpath, outputbase);
    strcat(outputpath, "_ExtIn.dat");
    ExtInFile = fopen(outputpath, "wb");
    fwrite(saveEin, sizeof(double), 2 * ((*activeSetPtr).Ngrid * (*activeSetPtr).NSims) + 1024, ExtInFile);
    fclose(ExtInFile);


    FILE* matlabfile;
    strcpy(outputpath, outputbase);
    strcat(outputpath, ".m");

    char* outputbaseVar = strrchr(outputbase, '\\');
    if (!outputbaseVar) {
        outputbaseVar = outputbase;
    }
    else {
        outputbaseVar++;
    }

    matlabfile = fopen(outputpath, "w");
    fprintf(matlabfile, "fid = fopen('%s_ExtIn.dat','rb'); \n", outputbaseVar);
    fprintf(matlabfile, "%s_ExtIn = fread(fid, %lli, 'double'); \n", outputbaseVar, 2 * (*activeSetPtr).Ngrid * (*activeSetPtr).NSims);
    fprintf(matlabfile, "%s_ExtIn = reshape(%s_ExtIn,[%lli %lli %i]); \n", outputbaseVar, outputbaseVar, (*activeSetPtr).Ntime, (*activeSetPtr).Nspace, 2 * (*activeSetPtr).NSims);
    fprintf(matlabfile, "fclose(fid); \n");

    fprintf(matlabfile, "fid = fopen('%s_ExtOut.dat','rb'); \n", outputbaseVar);
    fprintf(matlabfile, "%s_ExtOut = fread(fid, %lli, 'double'); \n", outputbaseVar, 2 * (*activeSetPtr).Ngrid * (*activeSetPtr).NSims);
    fprintf(matlabfile, "%s_ExtOut = reshape(%s_ExtOut,[%lli %lli %i]); \n", outputbaseVar, outputbaseVar, (*activeSetPtr).Ntime, (*activeSetPtr).Nspace, 2 * (*activeSetPtr).NSims);
    fprintf(matlabfile, "fclose(fid); \n");
    fclose(matlabfile);

    /*write a python script for loading the output fields in a proper shape*/

    char scriptfilename[MAX_PATH];
    strcpy(scriptfilename, outputbase);
    strcat(scriptfilename, ".py");
    FILE* scriptfile;
    scriptfile = fopen(scriptfilename, "w");
    fprintf(scriptfile, "#!/usr/bin/python\nimport numpy as np\n");
    fprintf(scriptfile, "dt = %e\ndz = %e\ndx = %e\n", (*activeSetPtr).tStep, (*activeSetPtr).propagationStep, (*activeSetPtr).rStep);
    fprintf(scriptfile, "%s_ExtIn = np.reshape(np.fromfile(\"", outputbaseVar);
    fprintf(scriptfile, "%s_ExtIn.dat", outputbaseVar);
    fprintf(scriptfile, "\",dtype=np.double)[0:%lli],(%lli,%lli,%i),order='F')\n", 2 * (*activeSetPtr).Ngrid * (*activeSetPtr).NSims, (*activeSetPtr).Ntime, (*activeSetPtr).Nspace, 2 * (*activeSetPtr).NSims);
    fprintf(scriptfile, "%s_ExtOut = np.reshape(np.fromfile(\"", outputbaseVar);
    fprintf(scriptfile, "%s_ExtOut.dat", outputbaseVar);
    fprintf(scriptfile, "\",dtype=np.double)[0:%lli],(%lli,%lli,%i),order='F')\n", 2 * (*activeSetPtr).Ngrid * (*activeSetPtr).NSims, (*activeSetPtr).Ntime, (*activeSetPtr).Nspace, 2 * (*activeSetPtr).NSims);
    fclose(scriptfile);
    return 0;
}

double vmaxa(double* v, int vlength) {
    double maxval = fabs(v[0]);
    int i;
    int imax = 0;
    for (i = 1; i < (vlength - 1); i++) {
        if (fabs(v[i]) > maxval) {
            maxval = fabs(v[i]);
            imax = i;
        }
    }

    if (imax > 0) {
        double no = fabs(v[imax]) / maxval;
        double np = fabs(v[imax + 1]) / maxval;
        double nm = fabs(v[imax - 1]) / maxval;
        maxval *= abs(2.0 * no * sqrt(no * no - nm * np) / sqrt(-(nm - 2.0 * no + np) * (nm + 2.0 * no + np)));
    }

    return maxval;
}

int resolvesequence(int currentIndex, struct propthread* s) {
    double pi = 3.1415926535897932384626433832795;
    
    //sequence format
    //material index, theta, phi, crystal length, propagation step, rotation angle
    int materialIndex = (int)(*s).sequenceArray[0 + 6*currentIndex];
    double crystalTheta = (pi/180) * (*s).sequenceArray[1 + 6 * currentIndex];
    double crystalPhi = (pi/180) * (*s).sequenceArray[2 + 6 * currentIndex];
    double propagationStep = 1e-9 * (*s).sequenceArray[4 + 6 * currentIndex];
    long long Npropagation = (long long)(1e-6*(*s).sequenceArray[3 + 6 * currentIndex]/propagationStep);
    double rotationAngle = (pi / 180) * (*s).sequenceArray[5 + 6 * currentIndex];

    if (currentIndex > 0) {
        (*s).isFollowerInSequence = TRUE;
    }

    (*s).propagationStep = propagationStep;
    (*s).Npropagation = Npropagation;


    (*s).materialIndex = materialIndex;
    (*s).crystalTheta = crystalTheta;
    (*s).crystalPhi = crystalPhi;
    (*s).chi2Tensor = crystalDatabasePtr[materialIndex].d;
    (*s).chi3Tensor = crystalDatabasePtr[materialIndex].chi3;
    (*s).nonlinearSwitches = crystalDatabasePtr[materialIndex].nonlinearSwitches;
    (*s).absorptionParameters = crystalDatabasePtr[materialIndex].absorptionParameters;
    (*s).sellmeierCoefficients = crystalDatabasePtr[materialIndex].sellmeierCoefficients;

    (*s).sellmeierType = crystalDatabasePtr[materialIndex].sellmeierType;
    (*s).axesNumber = crystalDatabasePtr[materialIndex].axisType;
    return 0;
}

//quality of life function - put a text label on a text box window, relative to its position
int LabelTextBox(HDC hdc, HWND parentWindow, HWND targetTextBox, const wchar_t* labelText, int xOffset, int yOffset) {
    RECT rectTextBox;
    POINT positionTextBox;
    GetWindowRect(targetTextBox, &rectTextBox);
    positionTextBox.x = rectTextBox.left;
    positionTextBox.y = rectTextBox.top;
    ScreenToClient(parentWindow, &positionTextBox);
    TextOutW(hdc, positionTextBox.x + xOffset, positionTextBox.y + yOffset, labelText, (int)_tcslen(labelText));
    return 0;
}

//reads the content of a text box and returns a double containing its numerical value
double HWNDToDouble(HWND inputA)
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

//Add a text string contained in messagebuffer to the text box inputA
int AppendTextToWindow(HWND inputA, wchar_t* messagebuffer, int buffersize) {
    int len = GetWindowTextLength(inputA);
    wchar_t* newbuffer = (wchar_t*)calloc(2 * len + 2 * buffersize, sizeof(wchar_t));
    if (len > 0) {
        len = GetWindowText(inputA, (LPWSTR)&newbuffer[0], len + 1);
    }
    for (int i = 0; i < buffersize; i++) {
        newbuffer[len + i] = messagebuffer[i];
    }
    SetWindowText(inputA, newbuffer);
    SendMessage(inputA, EM_LINESCROLL, 0, 99999);
    free(newbuffer);
    return 0;
}

int readcrystaldatabase(struct crystalentry* db, bool isVerbose) {
    int maxEntries = 64;
    int i;
    double* fd;
    wchar_t* messagebuffer;
    FILE* fp;
    fp = fopen("CrystalDatabase.txt","r");
    if (fp == NULL) {
        messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
        swprintf_s(messagebuffer, 1024,
            _T("Could not open database!\r\n"));
        AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
        free(messagebuffer);
        return 1;
    }
    if (isVerbose) {
        messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
        swprintf_s(messagebuffer, 1024,
            _T("Reading crystal database file in verbose mode\r\n"));
        AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
        free(messagebuffer);
    }
    //read the entries line
    fscanf(fp, "Total entries: %d\n", &maxEntries);

    for (i = 0; i < maxEntries; i++){
        fwscanf(fp, _T("Name:\n%[^\n]\n"), db[i].crystalNameW);
        fscanf(fp, "Type:\n%d\n", &db[i].axisType);
        fscanf(fp, "Sellmeier equation:\n%d\n", &db[i].sellmeierType);
        fd = &db[i].sellmeierCoefficients[0];
        fscanf(fp, "1st axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
        fd = &db[i].sellmeierCoefficients[22];
        fscanf(fp, "2nd axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
        fd = &db[i].sellmeierCoefficients[44];
        fscanf(fp, "3rd axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
        fwscanf(fp, _T("Sellmeier reference:\n%[^\n]\n"), db[i].sellmeierReference);
        fscanf(fp, "chi2 type:\n%d\n", &db[i].nonlinearSwitches[0]);
        fd = &db[i].d[0];
        fscanf(fp, "d:\n%lf %lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[3], &fd[6], &fd[9], &fd[12], &fd[15], &fd[1], &fd[4], &fd[7], &fd[10], &fd[13], &fd[16], &fd[2], &fd[5], &fd[8], &fd[11], &fd[14], &fd[17]);
        fwscanf(fp, _T("d reference:\n%[^\n]\n"), db[i].dReference);
        fscanf(fp, "chi3 type:\n%d\n", &db[i].nonlinearSwitches[1]);
        fd = &db[i].chi3[0];
        fscanf(fp, "chi3:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8]);
        fd = &db[i].chi3[9];
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8]);
        fd = &db[i].chi3[18];
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8]);
        fwscanf(fp, _T("chi3 reference:\n%[^\n]\n"), db[i].chi3Reference);
        fscanf(fp, "Nonlinear absorption type:\n%d\nAbsorption parameters:\n%lf %lf %lf %lf %lf %lf\n", &db[i].nonlinearSwitches[2], &db[i].absorptionParameters[0], &db[i].absorptionParameters[1], &db[i].absorptionParameters[2], &db[i].absorptionParameters[3], &db[i].absorptionParameters[4], &db[i].absorptionParameters[5]);
        fwscanf(fp, _T("Spectral file:\n%[^\n]\n"), db[i].spectralFile);
        fscanf(fp, "~~~crystal end~~~\n");
        if (isVerbose) {
            messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
            swprintf_s(messagebuffer, 1024,
                _T("Material %i name: %s\r\nSellmeier reference: %s\r\nChi2 reference: %s\r\nChi3 reference: %s\r\n\r\n"), i, db[i].crystalNameW, db[i].sellmeierReference, db[i].dReference, db[i].chi3Reference);
            AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
            free(messagebuffer);
        }


    }

    messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
    swprintf_s(messagebuffer, 1024,
        _T("Read %i entries\r\n"), (int)(maxEntries));
    AppendTextToWindow(maingui.textboxSims, messagebuffer, 1024);
    free(messagebuffer);

    fclose(fp);

    return 0;
}

//returns a string containing the text in a text box
int HWNDToString(HWND inputA, char* outputString, int bufferSize)
{
    int len = GetWindowTextLength(inputA);
    if (len > 0)
    {
        len = GetWindowTextA(inputA, outputString, bufferSize);
    }
    return 0;
}

int DrawLabels(HDC hdc) {
    int labos = -160;
    int x0plots = 380;
    int dxplots = 404;
    int dyplots = 214;
    int vs = 26;

    LabelTextBox(hdc, maingui.mainWindow, maingui.tbMaterialIndex, _T("Material index"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbCrystalTheta, _T("Crystal theta (deg)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbCrystalPhi, _T("Crystal phi (deg)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbBeamwaist1, _T("Beamwaist 1 (mcr.)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbBeamwaist2, _T("Beamwaist 2 (mcr.)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbGridXdim, _T("Grid width (mcr.)"), labos, 0);

    LabelTextBox(hdc, maingui.mainWindow, maingui.tbTimeStepSize, _T("dt (fs)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbRadialStepSize, _T("dx or dr (mcr.)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbTimeSpan, _T("Time span (fs)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbXstep, _T("dz (nm)"), labos, 0); 
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbCrystalThickness, _T("Thickness (mcr.)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.pdBatchMode, _T("Batch mode"), labos, 4);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbNumberSims, _T("Batch steps"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbBatchDestination, _T("Batch end"), labos, 0);

    LabelTextBox(hdc, maingui.mainWindow, maingui.tbPulse1Delay, _T("Delay 1 (fs)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbPulse2Delay, _T("Delay 2 (fs)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFieldStrength1, _T("Energy 1 (J)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFieldStrength2, _T("Energy 2 (J)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbBandwidth1, _T("Bandwidth 1 (THz)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbBandwidth2, _T("Bandwidth 2 (THz)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFrequency1, _T("Frequency 1 (THz)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFrequency2, _T("Frequency 2 (THz)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbCEPhase1, _T("CEP/pi 1"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbCEPhase2, _T("CEP/pi 2"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbPulseType, _T("SG order"), labos, 0);
    
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbGDD1, _T("GDD 1 (fs^2)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbGDD2, _T("GDD 2 (fs^2)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbTOD1, _T("TOD 1 (fs^3)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbTOD2, _T("TOD 2 (fs^3)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbXoffset1, _T("x offset 1 (mcr.)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbXoffset2, _T("x offset 2 (mcr.)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbZoffset1, _T("z offset 1 (mcr.)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbZoffset2, _T("z offset 2 (mcr.)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbPropagationAngle1, _T("NC angle 1 (deg)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbPropagationAngle2, _T("NC angle 2 (deg)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbPolarizationAngle1, _T("Polarization 1 (deg)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbPolarizationAngle2, _T("Polarization 2 (deg)"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbCircularity1, _T("Circularity 1"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbCircularity2, _T("Circularity 2"), labos, 0);
    LabelTextBox(hdc, maingui.mainWindow, maingui.pdPulse1Type, _T("Pulse 1 type:"), labos, 4);
    LabelTextBox(hdc, maingui.mainWindow, maingui.pdPulse2Type, _T("Pulse 2 type:"), labos, 4);
    LabelTextBox(hdc, maingui.mainWindow, maingui.pdPropagationMode, _T("Propagation mode"), labos, 4);

    LabelTextBox(hdc, maingui.mainWindow, maingui.tbSequence, _T("Crystal sequence:"), 4, -24);

    //plot labels
    int x = 690;
    int y = 125;
    int dx = 64 * 11;
    int dy = 256;
    int plotMargin = 0;
    int spacerX = 50;
    int spacerY = 40;
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("s-polarization, space/time:"), 0, 32);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("p-polarization, space/time:"), 0, 32+dy+spacerY);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("s-polarization waveform (GV/m):"), 0, 32 + 2*(dy + spacerY));
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("p-polarization waveform (GV/m):"), 0, 32 + 3*(dy + spacerY));
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("Time (fs)"), dx/2, 36 + 4 * (dy + spacerY));
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("s-polarization, Fourier, Log:"), 0+dx+spacerX, 32);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("p-polarization, Fourier, Log:"), 0+dx+spacerX, 32+dy+spacerY);
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("s-polarization spectrum, log-scale:"), 0 + dx + spacerX, 32 + 2*(dy + spacerY));
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("p-polarization spectrum, log-scale:"), 0 + dx + spacerX, 32 + 3*(dy + spacerY));
    LabelTextBox(hdc, maingui.mainWindow, maingui.tbFileNameBase, _T("Frequency (THz)"), dx / 2 + 0 + dx + spacerX, 36 + 4 * (dy + spacerY));
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

//Take an array of doubles and make an image in the window
//Color maps:
//  cm = 1: grayscale
//  cm = 2: similar to matlab's jet
int DrawArrayAsBitmap(HDC hdc, INT64 Nx, INT64 Ny, INT64 x, INT64 y, INT64 height, INT64 width, double* data, int cm) {

    // creating input
    unsigned char* pixels = (unsigned char*)calloc(4 * Nx * Ny, sizeof(unsigned char));
    INT64 i;
    INT64 Ntot = Nx * Ny;
    double nval;
    int stride = 4;
    //Find the image maximum and minimum
    double imin = data[0];
    double imax = data[0];
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

DWORD WINAPI drawsimplots(LPVOID lpParam) {
    bool wasRunning = isRunning;
    isRunning = TRUE; //this locks the grid memory so it doesn't get freed while plotting, set back to wasRunning at the end
    if (isGridAllocated) {
        int simIndex = (*activeSetPtr).plotSim;
        int x = 690;
        int y = 125;
        int dx = 64*11;
        int dy = 256;
        int plotMargin = 0;
        int spacerX = 50;
        int spacerY = 40;
        int i,j;
        HDC hdc;

        //int simIndex = (int)HWNDToDouble(maingui.tbWhichSimToPlot);
        //simIndex--;
        if (simIndex < 0) {
            simIndex = 0;
        }
        if (simIndex > (*activeSetPtr).NSims) {
            simIndex = (*activeSetPtr).NSims - 1;
        }

        hdc = GetWindowDC(maingui.mainWindow);
        double* plotarr = (double*)calloc((*activeSetPtr).Ngrid, sizeof(double));
        double* plotarrC = (double*)calloc((*activeSetPtr).Ngrid, sizeof(double));
        double* plotarr2 = (double*)calloc(dx * dy, sizeof(double));
        
        std::complex<double>* shiftedFFT = (std::complex<double>*)calloc((*activeSetPtr).Ngrid, sizeof(std::complex<double>));

        //Plot Time Domain, s-polarization
        for (i = 0; i < (*activeSetPtr).Ngrid; i++) {
            plotarr[i] = (real((*activeSetPtr).ExtOut[i + simIndex * (*activeSetPtr).Ngrid * 2]));
        }

        
        linearremap(plotarr, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime, plotarr2, (int)dy, (int)dx, 0);
        DrawArrayAsBitmap(hdc, dx, dy, x, y, dy, dx, plotarr2, 1);
        drawLabeledXYPlot(hdc, (int)(*activeSetPtr).Ntime, &plotarr[(*activeSetPtr).Ngrid / 2], (*activeSetPtr).tStep / 1e-15, x, y + 2 * dy + 2 * spacerY, (int)dx, (int)dy, 0, 0, 1e9);

        //Plot Time Domain, p-polarization
        for (i = 0; i < (*activeSetPtr).Ngrid; i++) {
            plotarr[i] = (real((*activeSetPtr).ExtOut[i + (*activeSetPtr).Ngrid + simIndex * (*activeSetPtr).Ngrid * 2]));
        }

        linearremap(plotarr, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime, plotarr2, (int)dy, (int)dx, 0);
        DrawArrayAsBitmap(hdc, dx, dy, x, y + dy + spacerY, dy, dx, plotarr2, 1);
        drawLabeledXYPlot(hdc, (int)(*activeSetPtr).Ntime, &plotarr[(*activeSetPtr).Ngrid / 2], (*activeSetPtr).tStep / 1e-15, x, y + 3 * dy + 3 * spacerY, (int)dx, (int)dy, 0, 0, 1e9);

        //Plot Fourier Domain, s-polarization
        fftshiftZ(&(*activeSetPtr).EkwOut[simIndex * (*activeSetPtr).Ngrid * 2], shiftedFFT, (*activeSetPtr).Ntime, (*activeSetPtr).Nspace);
        for (i = 0; i < (*activeSetPtr).Ngrid; i++) {
            plotarr[i] = log10(cmodulussquared(shiftedFFT[i]) + 1e-0);
        }

        for (i = 0; i < ((*activeSetPtr).Ntime / 2); i++) {
            for (j = 0; j < (*activeSetPtr).Nspace; j++) {
                plotarrC[i + ((*activeSetPtr).Ntime / 2) * j] = plotarr[i + (*activeSetPtr).Ntime * j + (*activeSetPtr).Ntime / 2];
            }
        }

        linearremap(plotarrC, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime / 2, plotarr2, (int)dy, (int)dx, 0);
        DrawArrayAsBitmap(hdc, dx, dy, x + dx + spacerX, y, dy, dx, plotarr2, 1);
        drawLabeledXYPlot(hdc, (int)(*activeSetPtr).Ntime / 2, &plotarrC[(*activeSetPtr).Ngrid / 4], (*activeSetPtr).fStep/1e12, x + dx + spacerX + plotMargin, y + 2 * dy + 2 * spacerY, dx, dy, 2, 8, 1);

        //Plot Fourier Domain, p-polarization
        fftshiftZ(&(*activeSetPtr).EkwOut[simIndex * (*activeSetPtr).Ngrid * 2 + (*activeSetPtr).Ngrid], shiftedFFT, (*activeSetPtr).Ntime, (*activeSetPtr).Nspace);
        for (i = 0; i < (*activeSetPtr).Ngrid; i++) {
            plotarr[i] = log10(cmodulussquared(shiftedFFT[i]) + 1e-0);
        }
        for (i = 0; i < ((*activeSetPtr).Ntime/2); i++) {
            for (j = 0; j < (*activeSetPtr).Nspace; j++) {
                plotarrC[i + ((*activeSetPtr).Ntime / 2) * j] = plotarr[i + (*activeSetPtr).Ntime * j + (*activeSetPtr).Ntime/2];
            }
        }
        
        linearremap(plotarrC, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime/2, plotarr2, (int)dy, (int)dx, 0);

        DrawArrayAsBitmap(hdc, dx, dy, x + dx + spacerX, y + dy + spacerY, dy, dx, plotarr2, 1);

        drawLabeledXYPlot(hdc, (int)(*activeSetPtr).Ntime/2, &plotarrC[(*activeSetPtr).Ngrid / 4], (*activeSetPtr).fStep/1e12, x + dx + spacerX + plotMargin, y + 3 * dy + 3 * spacerY, (int)dx, (int)dy, 2, 8, 1);
        
        free(shiftedFFT);
        free(plotarr);
        free(plotarr2);
        free(plotarrC);
        ReleaseDC(maingui.mainWindow, hdc);
    }
    isRunning = wasRunning;
    return 0;
}
int drawLabeledXYPlot(HDC hdc, int N, double* Y, double xStep, int posX, int posY, int pixelsWide, int pixelsTall, int forceYOrigin, double YOrigin, double yDiv) {
    double maxY = 0;
    double minY = 0;
    int i;
    double* X = (double*)calloc(N, sizeof(double));
    double* plotArray = (double*)calloc(pixelsWide * pixelsTall, sizeof(double));
    for (i = 0; i < N; i++) {
        X[i] = xStep * i;
        maxY = max(Y[i], maxY);
        minY = min(Y[i], minY);
    }

    int NyTicks = 3;
    double yTicks1[3] = { maxY, 0, minY };
    if (forceYOrigin == 2) {
        minY = maxY - YOrigin;
        NyTicks = 2;
        yTicks1[1] = minY;

    }
    else if (forceYOrigin == 1) {
        minY = YOrigin;
        NyTicks = 2;
        yTicks1[1] = minY;
    }


    double xTicks1[3] = { 0.25 * xStep * N, 0.5 * xStep * N, 0.75 * xStep * N };

    
    plotDataXY(X, Y, 0, xStep * N, minY, 1.02 * maxY, N, pixelsWide, pixelsTall, 1, 2.2, plotArray, xTicks1, 3, yTicks1, NyTicks);
    const wchar_t labelText[4] = _T("0.5");


    DrawArrayAsBitmap(hdc, pixelsWide, pixelsTall, posX, posY, pixelsTall, pixelsWide, plotArray, 0);
    wchar_t* messagebuffer;
    for (i = 0; i < NyTicks; i++) {
        messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
        swprintf_s(messagebuffer, 1024,
            _T("%1.1f"), yTicks1[i]/yDiv);
        TextOutW(hdc, posX - 32, posY + (int)(i * 0.96 * pixelsTall / 2), messagebuffer, (int)_tcslen(messagebuffer));
        free(messagebuffer);
    }
    for (i = 0; i < 3; i++) {
        messagebuffer = (wchar_t*)calloc(1024, sizeof(wchar_t));
        swprintf_s(messagebuffer, 1024,
            _T("%3.0f"), xTicks1[i]);
        TextOutW(hdc, posX + (int)(0.25 * pixelsWide * (i + 1) - 12), posY + pixelsTall, messagebuffer, (int)_tcslen(messagebuffer));
        free(messagebuffer);
    }
    free(X);
    free(plotArray);
    return 0;
}
//calculates the squard modulus of a complex number, under the assumption that the
//machine's complex number format is interleaved doubles.
//Orders of magnitude faster than the standard library abs()
double cmodulussquared(std::complex<double>complexNumber) {
    double* xy = (double*)&complexNumber;
    double modulusSquared = xy[0] * xy[0] + xy[1] * xy[1];
    return modulusSquared;
}


//use linear interpolation to resize matrix A to the size of matrix B
//B is overwritten with the resized matrix
int linearremap(double* A, int nax, int nay, double* B, int nbx, int nby, int modeInterp) {
    int i, j;
    double a, b, c, d;
    double d00, d01, d10, d11;
    double w00, w01, w10, w11;
    double A00, A01, A10, A11;
    double f, nf;

    int nx0, nx1, ny0, ny1;
    int Ni, Nj;

    if (modeInterp == 1) {
        for (i = 0; i < nbx; i++) {
            f = ((double)nax / (double)nbx) * (double)i;
            Ni = (int)f;
            a = f - Ni;
            b = 1. - a;
            nx0 = nay * min(Ni, nax - 1);
            nx1 = nay * min(Ni + 1, nax - 1);

            for (j = 0; j < nby; j++) {
                f = ((double)nay / (double)nby) * (double)j;
                Nj = (int)f;
                c = f - Nj;
                d = 1. - c;
                ny0 = min(nay - 1, Nj);
                ny1 = min(nay - 1, Nj + 1);

                A00 = A[ny0 + nx0];
                A01 = A[ny1 + nx0];
                A10 = A[ny0 + nx1];
                A11 = A[ny1 + nx1];

                d00 = sqrt(a * a + c * c);
                d10 = sqrt(b * b + c * c);
                d01 = sqrt(a * a + d * d);
                d11 = sqrt(b * b + d * d);

                w00 = 1 - d00;
                w00 *= w00;
                w01 = 1 - d01;
                w01 *= w01;
                w10 = 1 - d10;
                w10 *= w10;
                w11 = 1 - d11;
                w11 *= w11;
                nf = 1. / (w00 + w01 + w10 + w11);

                B[i * nby + j] = nf * (w00 * A00 + w10 * A10 + w01 * A01 + w11 * A11);
                //B[i * nby + j] = A00;
            }
        }
    }
    else {
        for (i = 0; i < nbx; i++) {
            f = ((double)nax / (double)nbx) * (double)i;
            Ni = (int)f;
            nx0 = nay * min(Ni, nax);
            for (j = 0; j < nby; j++) {
                f = ((double)nay / (double)nby) * (double)j;
                Nj = (int)f;
                ny0 = min(nay, Nj);
                A00 = A[ny0 + nx0];
                B[i * nby + j] = A00;
            }
        }
    }

    return 0;
}

