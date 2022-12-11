#include "LightwaveExplorerFrontendGTK.h"
#include "../LightwaveExplorerCore.cuh"
#include "../LightwaveExplorerCoreCPU.h"
#include "../LightwaveExplorerSYCL/LightwaveExplorerSYCL.h"

#include <chrono>
#include <Windows.h>
#include <delayimp.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <nvml.h>
#include <thread>
#define LABELWIDTH 6
#define MAX_LOADSTRING 1024
#define PI 3.1415926535897931
#define maxN(a,b)            (((a) > (b)) ? (a) : (b))
#define minN(a,b)            (((a) < (b)) ? (a) : (b))
bool isRunning = FALSE;
bool isPlotting = FALSE;
bool isGridAllocated = FALSE;
bool cancellationCalled = FALSE;
bool CUDAavailable = FALSE;
bool SYCLavailable = FALSE;
int cudaGPUCount = 0;
int syclGPUCount = 0;
size_t progressCounter = 0;

class mainGui {
    bool isActivated;
public:
    LweTextBox textBoxes[54];
    LweButton buttons[16];
    LweConsole console;
    LweConsole sequence;
    LweConsole fitCommand;
    LweTextBox filePaths[4];
    LwePulldown pulldowns[10];
    LweDrawBox drawBoxes[8];
    LweCheckBox checkBoxes[3];
    LweSlider plotSlider;
    LweWindow window;
    size_t pathTarget;
    bool isRunning = FALSE;
    bool isPlotting = FALSE;
    bool isGridAllocated = FALSE;
    bool cancellationCalled = FALSE;
    bool CUDAavailable = FALSE;
    bool SYCLavailable = FALSE;
    int cudaGPUCount = 0;
    int syclGPUCount = 0;
    size_t progressCounter = 0;
    mainGui() : isActivated(0), pathTarget(0), isRunning(0), isPlotting(0),
    isGridAllocated(0), cancellationCalled(0), CUDAavailable(0), SYCLavailable(0),
    cudaGPUCount(0), syclGPUCount(0){}
    ~mainGui() {}

    void activate(GtkApplication* app) {
        int buttonWidth = 4;
        int textWidth = 3;
        int labelWidth = 6;
        int plotWidth = 12;
        int plotHeight = 6;
        int colWidth = labelWidth + 2 * textWidth;
        int textCol1a = 6;
        int textCol2a = textCol1a + 2 * textWidth + labelWidth;
        int textCol1b = textCol1a + textWidth;
        int textCol2b = textCol2a + textWidth;
        int buttonCol1 = textCol2a - labelWidth;
        int buttonCol2 = buttonCol1 + buttonWidth;
        int buttonCol3 = buttonCol2 + buttonWidth;

        g_object_set(gtk_settings_get_default(),
            "gtk-application-prefer-dark-theme", TRUE,
            NULL);
        window.init(app, _T("Lightwave Explorer GTK"), 1080, 800);
        GtkWidget* parentHandle = window.parentHandle();
        for (int i = 0; i < 16; ++i) {
            textBoxes[i].init(parentHandle, textCol1a, i, textWidth, 1);
        }
        for (int i = 0; i < 16; ++i) {
            textBoxes[16 + i].init(parentHandle, textCol1b, i, textWidth, 1);
        }
        for (int i = 0; i < 6; ++i) {
            textBoxes[32 + 2 * i].init(parentHandle, textCol2a, i + 1, textWidth, 1);
            textBoxes[32 + 2 * i + 1].init(parentHandle, textCol2b, i + 1, textWidth, 1);
        }
        textBoxes[44].init(parentHandle, textCol2a, 10, textWidth, 1);
        textBoxes[45].init(parentHandle, textCol2b, 10, textWidth, 1);
        textBoxes[46].init(parentHandle, textCol2a, 11, textWidth, 1);
        textBoxes[47].init(parentHandle, textCol2b, 11, textWidth, 1);

        filePaths[0].init(parentHandle, 0, 17, colWidth, 1);
        filePaths[0].setMaxCharacters(36);
        pulldowns[0].addElement(_T("Synthetic"));
        pulldowns[0].addElement(_T("FROG"));
        pulldowns[0].addElement(_T("EOS"));
        pulldowns[0].init(parentHandle, labelWidth, 16, 2 * textWidth, 1);
        filePaths[0].setLabel(0, -1, _T("Data 1:"));

        filePaths[1].init(parentHandle, 0, 19, colWidth, 1);
        filePaths[1].setMaxCharacters(36);
        filePaths[1].setLabel(0, -1, _T("Data 2:"));
        pulldowns[1].addElement(_T("Synthetic"));
        pulldowns[1].addElement(_T("FROG"));
        pulldowns[1].addElement(_T("EOS"));
        pulldowns[1].init(parentHandle, labelWidth, 18, 2 * textWidth, 1);

        filePaths[2].init(parentHandle, 0, 21, colWidth, 1);
        filePaths[2].setMaxCharacters(36);
        filePaths[2].setLabel(0, -1, _T("Fit data:"));
        pulldowns[2].addElement(_T("Maximize x"));
        pulldowns[2].addElement(_T("Maximize y"));
        pulldowns[2].addElement(_T("Maximize Total"));
        pulldowns[2].addElement(_T("Fit spectrum"));
        pulldowns[2].init(parentHandle, labelWidth, 20, 2 * textWidth, 1);

        filePaths[3].init(parentHandle, 0, 23, colWidth, 1);
        filePaths[3].setLabel(0, -1, _T("Output:"));
        filePaths[3].setMaxCharacters(36);

        drawBoxes[0].init(window.parentHandle(2), 0, 0, plotWidth, plotHeight);
        drawBoxes[0].setDrawingFunction(drawTimeImage1);

        drawBoxes[1].init(window.parentHandle(2), 0, plotHeight, plotWidth, plotHeight);
        drawBoxes[1].setDrawingFunction(drawTimeImage2);

        drawBoxes[2].init(window.parentHandle(2), 0, 2 * plotHeight, plotWidth, plotHeight);
        drawBoxes[2].setDrawingFunction(drawField1Plot);

        drawBoxes[3].init(window.parentHandle(2), 0, 3 * plotHeight, plotWidth, plotHeight);
        drawBoxes[3].setDrawingFunction(drawField2Plot);

        drawBoxes[4].init(window.parentHandle(2), plotWidth, 0, plotWidth, plotHeight);
        drawBoxes[4].setDrawingFunction(drawField1Plot);

        drawBoxes[5].init(window.parentHandle(2), plotWidth, plotHeight, plotWidth, plotHeight);
        drawBoxes[5].setDrawingFunction(drawField1Plot);

        drawBoxes[6].init(window.parentHandle(2), plotWidth, 2 * plotHeight, plotWidth, plotHeight);
        drawBoxes[6].setDrawingFunction(drawSpectrum1Plot);
        drawBoxes[7].init(window.parentHandle(2), plotWidth, 3 * plotHeight, plotWidth, plotHeight);
        drawBoxes[7].setDrawingFunction(drawSpectrum2Plot);
        textBoxes[48].init(window.parentHandle(4), 2, 0, 2, 1);
        textBoxes[49].init(window.parentHandle(4), 4, 0, 2, 1);
        textBoxes[50].init(window.parentHandle(4), 8, 0, 2, 1);
        textBoxes[51].init(window.parentHandle(4), 10, 0, 2, 1);

        //textBoxes[52].init(window.parentHandle(3), 4, 0, 1, 1);
        plotSlider.init(window.parentHandle(3), 0, 0, 4, 1);
        plotSlider.setRange(0.0, 10.0);
        plotSlider.setDigits(0);

        checkBoxes[0].init(_T("Total"), window.parentHandle(4), 12, 0, 1, 1);
        checkBoxes[1].init(_T("Log"), window.parentHandle(4), 13, 0, 1, 1);
        pulldowns[3].addElement(_T("00: Vacuum"));
        pulldowns[3].addElement(_T("01: BBO"));
        pulldowns[3].init(parentHandle, textCol2a, 0, 2 * textWidth, 1);

        pulldowns[4].addElement(_T("2D Cartesian"));
        pulldowns[4].addElement(_T("3D radial symm."));
        pulldowns[4].addElement(_T("3D"));
        pulldowns[4].init(parentHandle, textCol2a, 7, 2 * textWidth, 1);

        pulldowns[5].addElement(_T("none"));
        pulldowns[5].addElement(_T("oh lord"));
        pulldowns[5].addElement(_T("write the fill"));
        pulldowns[5].init(parentHandle, textCol2a, 8, 2 * textWidth, 1);

        pulldowns[6].addElement(_T("none"));
        pulldowns[6].addElement(_T("oh lord"));
        pulldowns[6].addElement(_T("write the fill"));
        pulldowns[6].init(parentHandle, textCol2a, 9, 2 * textWidth, 1);

        pulldowns[7].addElement(_T("CUDA"));
        pulldowns[7].addElement(_T("SYCL"));
        pulldowns[7].addElement(_T("OpenMP"));
        pulldowns[7].init(window.parentHandle(5), labelWidth, 0, buttonWidth, 1);

        pulldowns[8].addElement(_T("CUDA"));
        pulldowns[8].addElement(_T("SYCL"));
        pulldowns[8].addElement(_T("OpenMP"));
        pulldowns[8].init(window.parentHandle(5), 2 * labelWidth + buttonWidth, 0, buttonWidth, 1);
        textBoxes[52].init(window.parentHandle(5), 2 * labelWidth + 2 * buttonWidth, 0, 1, 1);
        pulldowns[9].addElement(_T("Cobra 1x R5000"));
        pulldowns[9].addElement(_T("Cobra 2x R5000"));
        pulldowns[9].addElement(_T("Raven 1x A100"));
        pulldowns[9].init(parentHandle, 0, 24, 2 * buttonWidth, 1);

        sequence.init(parentHandle, buttonCol1, 13, colWidth, 6);
        fitCommand.init(parentHandle, buttonCol1, 21, colWidth, 4);

        buttons[0].init(_T("Run"), parentHandle, buttonCol3, 19, buttonWidth, 1, runButtonClick);
        buttons[1].init(_T("Stop"), parentHandle, buttonCol2, 19, buttonWidth, 1, runButtonClick);
        buttons[2].init(_T("Cluster"), parentHandle, 2 * buttonWidth, 24, buttonWidth, 1, runButtonClick);
        buttons[3].init(_T("Fit"), parentHandle, buttonCol3, 20, buttonWidth, 1, runButtonClick);
        buttons[4].init(_T("Load"), parentHandle, buttonCol1, 19, buttonWidth, 1, runButtonClick);
        //buttons[5].init(_T("Reload"), parentHandle, buttonCol2, 18, buttonWidth, 1, runButtonClick);
        buttons[6].init(_T("Path"), parentHandle, textWidth, 16, textWidth, 1, openFileDialogCallback, 0);
        buttons[7].init(_T("Path"), parentHandle, textWidth, 18, textWidth, 1, openFileDialogCallback, (gpointer)1);
        buttons[8].init(_T("Path"), parentHandle, textWidth, 20, textWidth, 1, openFileDialogCallback, (gpointer)2);
        buttons[9].init(_T("Path"), parentHandle, textWidth, 22, textWidth, 1, saveFileDialogCallback, (gpointer)3);
        buttons[10].init(_T("xlim"), window.parentHandle(4), 0, 0, 1, 1, runButtonClick);
        buttons[11].init(_T("ylim"), window.parentHandle(4), 6, 0, 1, 1, runButtonClick);
        buttons[12].init(_T("SVG"), window.parentHandle(3), 5, 0, 1, 1, runButtonClick);

        console.init(window.parentHandle(1), 0, 0, 1, 1);

        textBoxes[0].setLabel(-labelWidth, 0, _T("Pulse energy (J)"));
        textBoxes[1].setLabel(-labelWidth, 0, _T("Frequency (THz)"));
        textBoxes[2].setLabel(-labelWidth, 0, _T("Bandwidth (THz)"));
        textBoxes[3].setLabel(-labelWidth, 0, _T("SG order"));
        textBoxes[4].setLabel(-labelWidth, 0, _T("CEP/pi"));
        textBoxes[5].setLabel(-labelWidth, 0, _T("Delay (fs)"));
        textBoxes[6].setLabel(-labelWidth, 0, _T("GDD (fs^2)"));
        textBoxes[7].setLabel(-labelWidth, 0, _T("TOD (fs^2)"));
        textBoxes[8].setLabel(-labelWidth, 0, _T("Phase material"));
        textBoxes[9].setLabel(-labelWidth, 0, _T("Thickness (\xce\xbcm)"));
        textBoxes[10].setLabel(-labelWidth, 0, _T("Beamwaist (\xce\xbcm)"));
        textBoxes[11].setLabel(-labelWidth, 0, _T("x offset (\xce\xbcm)"));
        textBoxes[12].setLabel(-labelWidth, 0, _T("z offset (\xce\xbcm)"));
        textBoxes[13].setLabel(-labelWidth, 0, _T("NC angle (deg)"));
        textBoxes[14].setLabel(-labelWidth, 0, _T("Polarization (deg)"));
        textBoxes[15].setLabel(-labelWidth, 0, _T("Circularity"));
        textBoxes[32].setLabel(-labelWidth, 0, _T("Theta, phi (deg)"));
        textBoxes[34].setLabel(-labelWidth, 0, _T("NL absorption"));
        textBoxes[36].setLabel(-labelWidth, 0, _T("Drude: gamma, m"));
        textBoxes[38].setLabel(-labelWidth, 0, _T("Max x, dx (\xce\xbcm)"));
        textBoxes[40].setLabel(-labelWidth, 0, _T("Time span, dt (fs)"));
        textBoxes[42].setLabel(-labelWidth, 0, _T("Max z, dz (\xce\xbcm,nm)"));
        pulldowns[3].setLabel(-labelWidth, 0, _T("Material"));
        textBoxes[44].setLabel(-labelWidth, 0, _T("Batch end"));
        textBoxes[46].setLabel(-labelWidth, 0, _T("Batch steps"));
        pulldowns[4].setLabel(-labelWidth, 0, _T("Propagation"));
        pulldowns[5].setLabel(-labelWidth, 0, _T("Batch mode"));
        pulldowns[6].setLabel(-labelWidth, 0, _T("Batch mode 2"));
        pulldowns[7].setLabel(-labelWidth, 0, _T("Run type:"));
        pulldowns[8].setLabel(-labelWidth, 0, _T(" Offload type:"));
        fitCommand.setLabel(0, -1, _T("Fitting:"));
        sequence.setLabel(0, -1, _T("Sequence:"));

        window.present();

        checkLibraryAvailability();
    }
};
mainGui theGui;
simulationParameterSet* activeSetPtr;       // Main structure containing simulation parameters and pointers
crystalEntry* crystalDatabasePtr;           // Crystal info database
char programDirectory[MAX_LOADSTRING];     // Program working directory (useful if the crystal database has to be reloaded)

///////////////////.
//Definitions over
///////////////////.
void checkLibraryAvailability() {
    theGui.console.cPrint("\r\n");
    __try {
        HRESULT hr = __HrLoadAllImportsForDll("cufft64_10.dll");
        if (SUCCEEDED(hr)) {
            CUDAavailable = TRUE;
            theGui.console.cPrint("gofft.\r\n");
        }
    }
    __except (EXCEPTION_EXECUTE_HANDLER) {
        CUDAavailable = FALSE;
        theGui.console.cPrint("CUDA not available because cufft64_10.dll was not found.\r\n");
    }
    if (CUDAavailable) {
        CUDAavailable = FALSE;
        __try {
            HRESULT hr = __HrLoadAllImportsForDll("nvml.dll");
            if (SUCCEEDED(hr)) {
                CUDAavailable = TRUE;
                theGui.console.cPrint("gofnml.\r\n");
            }
        }
        __except (EXCEPTION_EXECUTE_HANDLER) {
            CUDAavailable = FALSE;
            theGui.console.cPrint("CUDA not available.\r\n");
        }
    }
    if (CUDAavailable) {
        CUDAavailable = FALSE;
        __try {
            HRESULT hr = __HrLoadAllImportsForDll("cudart64_110.dll");
            if (SUCCEEDED(hr)) {
                CUDAavailable = TRUE;
                theGui.console.cPrint("nort.\r\n");
            }
        }
        __except (EXCEPTION_EXECUTE_HANDLER) {
            CUDAavailable = FALSE;
            theGui.console.cPrint("CUDA not available because cudart64_110.dll was not found.\r\n\r\n");
        }
    }
    CUDAavailable = TRUE;
    if (CUDAavailable) {
        //Find, count, and name the GPUs
        int CUDAdevice, i;

        cudaGetDeviceCount(&cudaGPUCount);
        cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
        struct cudaDeviceProp activeCUDADeviceProp;
        wchar_t wcstring[514];
        size_t convertedChars = 0;
        //if (cuErr == cudaSuccess) {
        if(TRUE){
            if (cudaGPUCount > 0) {
                CUDAavailable = TRUE;
            }

            theGui.console.cPrint("CUDA found %i GPU(s): \r\n", cudaGPUCount);
            for (i = 0; i < cudaGPUCount; ++i) {
                cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
                memset(wcstring, 0, sizeof(wchar_t));
                mbstowcs_s(&convertedChars, wcstring, 256, activeCUDADeviceProp.name, _TRUNCATE);
                theGui.console.cPrint("%ls\r\n", wcstring);
                theGui.console.cPrint("    Memory: %i MB\r\n    Multiprocessors: %i\r\n",
                    (int)ceil(((float)activeCUDADeviceProp.totalGlobalMem) / 1048576), activeCUDADeviceProp.multiProcessorCount);
            }

        }
        else {
            theGui.console.cPrint("No CUDA-compatible GPU found.\r\n");
            CUDAavailable = FALSE;
        }
    }
    else {
        theGui.console.cPrint("No CUDA-compatible GPU found.\r\n");
        CUDAavailable = FALSE;
    }

    //read SYCL devices
    wchar_t syclDeviceList[MAX_LOADSTRING] = { 0 };
    wchar_t syclDefault[MAX_LOADSTRING] = { 0 };
    __try {
        HRESULT hr = __HrLoadAllImportsForDll("LightwaveExplorerSYCL.dll");
        if (SUCCEEDED(hr)) {
            SYCLavailable = TRUE;
        }
    }
    __except (EXCEPTION_EXECUTE_HANDLER) {
        SYCLavailable = FALSE;
        DWORD SYCLfile = GetFileAttributes("LightwaveExplorerSYCL.dll");
        if (SYCLfile != 0xFFFFFFFF) {
            theGui.console.cPrint("Couldn't run SYCL... \r\nHave you installed the DPC++ runtime from Intel?\r\n");
            theGui.console.cPrint("https://www.intel.com/content/www/us/en/developer/articles/tool/compilers-redistributable-libraries-by-version.html\r\n");
        }
        else {
            theGui.console.cPrint("No SYCL file...\r\n");
        }
    }

    if (SYCLavailable) {
        size_t syclDevices = 0;
        __try {
            syclDevices = readSYCLDevices(syclDeviceList, syclDefault);
        }
        __except (EXCEPTION_EXECUTE_HANDLER) {
            SYCLavailable = FALSE;
            theGui.console.cPrint("Couldn't run SYCL... \r\nHave you installed the DPC++ runtime from Intel?\r\n");
            theGui.console.cPrint("https://www.intel.com/content/www/us/en/developer/articles/tool/compilers-redistributable-libraries-by-version.html\r\n");
        }
        unsigned char* deviceCounts = (unsigned char*)&syclDevices;
        if (deviceCounts[0] == 0u) {
            theGui.console.cPrint("Something is wrong - SYCL doesn't think you have a CPU.\r\n");
            SYCLavailable = FALSE;
        }
        else {
            syclGPUCount = deviceCounts[1];
            //printToConsole(maingui.plotBox2, syclDeviceList);
            //if (syclGPUCount > 0) printToConsole(maingui.plotBox2, syclDefault);
        }
    }
}
int drawArrayAsBitmap(cairo_t* cr, int Nx, int Ny, float* data, int cm) {
    if (Nx * Ny == 0) return 1;
    
    // creating input
    unsigned char* pixels = (unsigned char*)calloc(4 * Nx * Ny, sizeof(unsigned char));
    if (pixels == NULL) return 1;


    size_t Ntot = Nx * Ny;
    float nval;
    int stride = 4;
    //Find the image maximum and minimum
    float imin = data[0];
    float imax = data[0];
    for (size_t i = 1; i < Ntot; ++i) {
        if (data[i] > imax) imax = data[i];
        if (data[i] < imin) imin = data[i];
    }

    float oneOver255 = 1.0f / 255;
    unsigned char colorMap[256][3];

    for (int j = 0; j < 256; ++j) {
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
            nval = 255 * (j * oneOver255);
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
            colorMap[j][2] = (unsigned char)(255. * (0.9 * exp(-pow(4.5 * (nval - 0.05), 2)) + 0.2 * exp(-pow(3.5 * (nval - 1.05), 2))));
        }
    }

    if (cm == 4) {
        imax = maxN(imax, -imin);
        imin = minN(imin, -imax);
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

    int caiStride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, Nx);
    cairo_surface_t* cSurface = cairo_image_surface_create_for_data(pixels, CAIRO_FORMAT_RGB24, Nx, Ny, caiStride);
    cairo_set_source_surface(cr, cSurface, 0, 0);
    cairo_paint(cr);
    cairo_surface_destroy(cSurface);
    return 0;
}

void pathFromDialogBox(GtkDialog* dialog, int response) {
    if (response == GTK_RESPONSE_ACCEPT) {
        GtkFileChooser* chooser = GTK_FILE_CHOOSER(dialog);
        GFile* file = gtk_file_chooser_get_file(chooser);
        char* path = g_file_get_path(file);
        theGui.filePaths[theGui.pathTarget].overwritePrint("%s", path);
    }
    g_object_unref(dialog);
}

void openFileDialog(size_t pathTarget) {
    theGui.pathTarget = pathTarget;
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new("Open File", theGui.window.windowHandle(), GTK_FILE_CHOOSER_ACTION_OPEN, "Ok", "Cancel");
    g_signal_connect(fileC, "response", G_CALLBACK(pathFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
}

void openFileDialogCallback(GtkWidget* widget, gpointer pathTarget) {
    theGui.pathTarget = (size_t)pathTarget;
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new("Open File", theGui.window.windowHandle(), GTK_FILE_CHOOSER_ACTION_OPEN, "Ok", "Cancel");
    g_signal_connect(fileC, "response", G_CALLBACK(pathFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
}

void saveFileDialogCallback(GtkWidget* widget, gpointer pathTarget) {
    theGui.pathTarget = (size_t)pathTarget;
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new("Save File", theGui.window.windowHandle(), GTK_FILE_CHOOSER_ACTION_SAVE, "Ok", "Cancel");
    GtkFileChooser* chooser = GTK_FILE_CHOOSER(fileC);
    gtk_file_chooser_set_file(chooser, NULL, NULL);

    g_signal_connect(fileC, "response", G_CALLBACK(pathFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
}

static void runButtonClick() {
    //printToConsole("Testing... %lf\n",theGui.textBoxes[0].valueDouble());
    theGui.console.cPrint("Testing...\nFirst field: %lf\nPulldown: %i\nCheckbox: %s\n", theGui.textBoxes[0].valueDouble(), theGui.pulldowns[0].getValue(), theGui.checkBoxes[0].isChecked() ? "true" : "false");
    theGui.textBoxes[0].setToDouble(4.0);
    theGui.textBoxes[1].setToDouble(2.1e-9);
}

static void activate(GtkApplication* app, gpointer user_data) {
    theGui.activate(app);
}

int LwePlot2d(plotStruct* inputStruct) {
    plotStruct* s = (plotStruct*)inputStruct;
    if ((*s).Npts == 0) return 1;
    size_t iMin = 0;
    size_t iMax = (*s).Npts;

    cairo_t* cr = (*s).cr;
    cairo_font_extents_t fe;
    cairo_text_extents_t te;

    double fontSize = 14.0;
    cairo_set_font_size(cr, fontSize);
    cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_font_extents(cr, &fe);
    double x1, y1, x2, y2;
    double width;
    double height;
    size_t strLen;

    double layoutTop = 0.0;
    double layoutBottom = 0.0;
    double layoutLeft = 0.0;
    double layoutRight = 0.0;

    //Fixed parameters affecting aesthetics
    double radius = 2;
    double lineWidth = 1.5;
    double axisSpaceX = 70.0;
    double axisSpaceY = 35.0;
    double axisLabelSpaceX = 16.0;

    //get limits and make the plotting arrays
    double maxY = 0.0;
    double minY = 0.0;
    double maxX = 0.0;
    double minX = 0.0;
    double currentY;
    double currentX;
    double* xValues = new double[(*s).Npts + 2]();
    for (int i = 0; i < (*s).Npts; ++i) {
        if ((*s).hasDataX) currentX = (double)(*s).dataX[i];
        else { currentX = (double)(i * (*s).dx + (*s).x0); }
        if (i == 0) {
            minX = currentX;
            maxX = currentX;
        }
        xValues[i] = currentX;
        if ((*s).forceXmin && (currentX < (*s).forcedXmin)) {
            iMin = i + 1;
        }
        if ((*s).forceXmax && (currentX > (*s).forcedXmax)) {
            iMax = i;
            break;
        }
        maxX = maxN(currentX, maxX);
        minX = minN(currentX, minX);
    }

    for (size_t i = iMin; i < iMax; ++i) {
        if ((*s).logScale) { currentY = (double)log10((*s).data[i]); }
        else { currentY = (double)(*s).data[i]; }
        maxY = maxN(currentY, maxY);
        minY = minN(currentY, minY);
        if ((*s).ExtraLines > 0) {
            if ((*s).logScale) { currentY = (double)log10((*s).data2[i]); }
            else { currentY = (double)(*s).data2[i]; }
            maxY = maxN(currentY, maxY);
            minY = minN(currentY, minY);
        }
    }

    if (minY == maxY) {
        minY = -1;
        maxY = 1;
    }
    if ((*s).forceYmin) {
        minY = (double)(*s).forcedYmin;
    }
    if ((*s).forceYmax) {
        maxY = (double)(*s).forcedYmax;
    }
    if ((*s).forceXmin) {
        minX = (double)(*s).forcedXmin;
    }
    if ((*s).forceXmax) {
        maxX = (double)(*s).forcedXmax;
    }

    //Tickmark labels
    int NyTicks = 3;
    char messageBuffer[4096] = { 0 };
    double yTicks1[3] = { maxY, 0.5 * (maxY + minY), minY };
    double xTicks1[3] = { minX + 0.25 * (maxX - minX), minX + 0.5 * (maxX - minX), minX + 0.75 * (maxX - minX) };
    height = (double)(*s).height;
    width = (double)(*s).width;

    //Start SVG file if building one
    auto SVGh = [&](double x) {
        return (int)(15 * x);
    };
    if ((*s).makeSVG) {
        snprintf((*s).svgString, (*s).svgBuffer, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
        snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "<svg width=\"%f\" height=\"%f\" viewBox=\"0 0 %f %f\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n",
            width, height, width, height);
        snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "<rect fill=\"#%X%X%X\" stroke=\"#000\" x=\"0\" y=\"0\" width=\"%f\" height=\"%f\"/>\n",
            SVGh(0.0f), SVGh(0.0f), SVGh(0.0f), width, height);
    }
    LweColor black(0, 0, 0, 0);
    cairo_rectangle(cr, 0, 0, width, height);
    black.setCairo(cr);
    cairo_fill(cr);
    width -= axisSpaceX;
    height -= axisSpaceY;
    double scaleX = width / ((double)(maxX - minX));
    double scaleY = height / ((double)(maxY - minY));
    LweColor currentColor = (*s).textColor;


    //make the paintbrush
    //hr = pRenderTarget->CreateSolidColorBrush((*s).textColor, &pBrush);
    currentColor = (*s).textColor;
    //lambdas for writing components of SVG file
    auto SVGstdline = [&]() {
        if ((*s).makeSVG)snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"#%X%X%X\" stroke-width=\"%f\"/>\n", x1, y1, x2, y2, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), lineWidth);
    };

    auto SVGstdcircle = [&]() {
        if ((*s).makeSVG)snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"none\" fill=\"#%X%X%X\" />\n", x1, y1, radius, currentColor.rHex(), currentColor.gHex(), currentColor.bHex());
    };

    auto SVGstartgroup = [&]() {
        if ((*s).makeSVG)snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "<g>\n");
    };

    auto SVGendgroup = [&]() {
        if ((*s).makeSVG)snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "</g>\n");
    };

    auto SVGcentertext = [&]() {
        if ((*s).makeSVG)snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "<text font-family=\"Arial\" font-size=\"%f\" fill=\"#%X%X%X\" x=\"%f\" y=\"%f\" text-anchor=\"middle\">\n%s\n</text>\n", fontSize, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 0.5 * (layoutLeft + layoutRight), layoutTop + fontSize, messageBuffer);
    };

    auto SVGlefttext = [&]() {
        if ((*s).makeSVG)snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "<text font-family=\"Arial\" font-size=\"%f\" fill=\"#%X%X%X\" x=\"%f\" y=\"%f\">\n%s\n</text>\n", fontSize, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), layoutLeft, layoutTop + fontSize, messageBuffer);
    };


    cairo_set_line_width(cr, lineWidth);
    auto cairoLine = [&]() {
        currentColor.setCairo(cr);
        cairo_move_to(cr, x1, y1);
        cairo_line_to(cr, x2, y2);
        cairo_stroke(cr);
    };
    auto cairoCircle = [&]() {
        currentColor.setCairo(cr);
        cairo_arc(cr, x1, y1, radius, 0, 6.2831853071795862);
        cairo_fill(cr);
    };
    auto cairoLeftText = [&]() {
        currentColor.setCairo(cr);
        cairo_text_extents(cr, messageBuffer, &te);
        cairo_move_to(cr, layoutLeft, 0.5 * (layoutBottom + layoutTop - te.height));
        cairo_show_text(cr, messageBuffer);
    };
    auto cairoCenterText = [&]() {
        currentColor.setCairo(cr);
        cairo_text_extents(cr, messageBuffer, &te);
        cairo_move_to(cr, 0.5 * (layoutLeft + layoutRight - te.width), 0.5 * (layoutBottom + layoutTop - te.height));
        cairo_show_text(cr, messageBuffer);
    };

    auto cairoVerticalText = [&]() {
        currentColor.setCairo(cr);
        cairo_text_extents(cr, messageBuffer, &te);
        cairo_move_to(cr, 0.0, height);
        cairo_rotate(cr, -3.1415926535897931 / 2);
        cairo_rel_move_to(cr, 0.5 * (layoutLeft + layoutRight - te.width), fontSize);
        cairo_show_text(cr, messageBuffer);
        cairo_rotate(cr, 3.1415926535897931 / 2);
    };

    currentColor = (*s).textColor;
    //y-tick text labels
    for (int i = 0; i < NyTicks; ++i) {
        memset(messageBuffer, 0, MAX_LOADSTRING * sizeof(wchar_t));
        if (abs(yTicks1[i] / (*s).unitY) > 10.0 || abs(yTicks1[i] / (*s).unitY) < 0.01) {
            snprintf(messageBuffer, MAX_LOADSTRING,
                _T("%1.1e "), yTicks1[i] / (*s).unitY);
        }
        else {
            snprintf(messageBuffer, MAX_LOADSTRING,
                _T("%1.3f "), yTicks1[i] / (*s).unitY);

        }
        layoutLeft = axisLabelSpaceX;
        layoutTop = (double)(i * (0.5 * (height)));
        if (i == 2) layoutTop -= 8.0f;
        if (i == 1) layoutTop -= 6.0f;
        layoutBottom = layoutTop + axisSpaceY;
        layoutRight = axisSpaceX;
        strLen = strlen(messageBuffer);

        cairoLeftText();
        SVGlefttext();
    }
    //y-axis name
    if ((*s).yLabel != NULL) {
        memset(messageBuffer, 0, MAX_LOADSTRING * sizeof(wchar_t));
        snprintf(messageBuffer, MAX_LOADSTRING,
            _T("%s"), (*s).yLabel);
        layoutLeft = 0;
        layoutTop = height;
        layoutBottom = height + axisSpaceY;
        layoutRight = height;
        strLen = strlen(messageBuffer);
        cairoVerticalText();
        if ((*s).makeSVG)snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "<text font-family=\"Arial\" font-size=\"%f\" fill=\"#%X%X%X\" x=\"%f\" y=\"%f\" text-anchor=\"middle\" transform=\"translate(%f, %f) rotate(-90)\">\n%s\n</text>\n", fontSize, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 0.5 * (layoutLeft + layoutRight), layoutTop + fontSize, -(layoutLeft + layoutRight), height, messageBuffer);
    }

    //x-axis name
    if ((*s).xLabel != NULL) {
        layoutLeft = axisSpaceX;
        layoutTop = height + 2 * fontSize;
        layoutBottom = height + axisSpaceY;
        layoutRight = axisSpaceX + width;
        strLen = strlen(messageBuffer);
        memset(messageBuffer, 0, MAX_LOADSTRING * sizeof(wchar_t));
        snprintf(messageBuffer, MAX_LOADSTRING,
            _T("%s"), (*s).xLabel);
        strLen = strlen(messageBuffer);
        cairoCenterText();
        SVGcentertext();

    }

    //x-axis tick labels
    for (int i = 0; i < 3; ++i) {
        memset(messageBuffer, 0, MAX_LOADSTRING * sizeof(wchar_t));
        snprintf(messageBuffer, MAX_LOADSTRING,
            _T("%i"), (int)round(xTicks1[i]));
        layoutLeft = (double)(axisSpaceX + 0.25 * width * ((size_t)(i)+1) - axisSpaceX / 2);
        layoutTop = height;
        layoutBottom = height + axisSpaceY;
        layoutRight = layoutLeft + axisSpaceX;
        strLen = strlen(messageBuffer);
        cairoCenterText();
        SVGcentertext();
    }

    //Draw axes and tickmarks
    SVGstartgroup();
    currentColor = (*s).axisColor;
    x1 = axisSpaceX;
    y1 = height;
    x2 = scaleX * (maxX - minX) + axisSpaceX;
    y2 = height;
    cairoLine();
    SVGstdline();
    x1 = axisSpaceX;
    y1 = height;
    x2 = x1;
    y2 = 0.0;
    cairoLine();
    SVGstdline();
    for (int i = 0; i < 2; ++i) {
        y1 = (double)(height - scaleY * (yTicks1[i] - minY));
        y2 = y1;
        x1 = axisSpaceX;
        x2 = x1 + 10;
        cairoLine();
        SVGstdline();
    }
    for (int i = 0; i < 3; ++i) {
        y1 = height;
        y2 = y1 - 10;
        x1 = (double)(axisSpaceX + scaleX * (xTicks1[i] - minX));
        x2 = x1;
        cairoLine();
        SVGstdline();
    }
    SVGendgroup();

    //Lambda for plotting a single line
    //note: SVG could be more compact (and maybe faster to load?) using a big polyline call rather
    //than individual lines like it is now.
    auto plotLine = [&](double* y) {
        SVGstartgroup();
        for (size_t i = iMin; i < iMax - 1; ++i) {
            x1 = scaleX * (xValues[i] - minX);
            x2 = scaleX * (xValues[i + 1] - minX);
            if ((*s).logScale) {
                y1 = height - scaleY * ((double)log10(y[i]) - (double)minY);
                y2 = height - scaleY * ((double)log10(y[i + 1]) - (double)minY);
            }
            else {
                y1 = height - scaleY * ((double)y[i] - (double)minY);
                y2 = height - scaleY * ((double)y[i + 1] - (double)minY);
            }
            x1 += axisSpaceX;
            x2 += axisSpaceX;

            if (y1 <= height) {
                if (y2 <= height) {
                    cairoLine();
                    SVGstdline();
                }
                else {
                    x2 = x1 + (height - y1) / ((y2 - y1) / (x2 - x1));
                    y2 = height;
                    cairoLine();
                    SVGstdline();
                }
                cairoCircle();
                SVGstdcircle();
            }
            else if (y2 <= height) {
                x1 = x1 + (height - y1) / ((y2 - y1) / (x2 - x1));
                y1 = height;
                cairoLine();
                SVGstdline();
            }
        }
        SVGendgroup();
    };

    //Plot the main line
    currentColor = (*s).color;
    plotLine((*s).data);

    //Optional overlay curves
    if ((*s).ExtraLines > 0) {
        currentColor = (*s).color2;
        plotLine((*s).data2);
    }
    if ((*s).ExtraLines > 1) {
        currentColor = (*s).color3;
        plotLine((*s).data3);
    }
    if ((*s).ExtraLines > 2) {
        currentColor = (*s).color2;
        plotLine((*s).data4);
    }
    delete[] xValues;
    if ((*s).makeSVG) {
        snprintf((*s).svgString + strlen((*s).svgString), (*s).svgBuffer, "</svg>");
    }
    return 0;
}

void drawField1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    GtkStyleContext* context;
    LweColor black(0, 0, 0, 0);
    LweColor white(0.8, 0.8, 0.8, 0);
    LweColor magenta(1, 0, 1, 0);

    double xData[10] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. };
    double yData[10] = { 1., 2., 3., 4., 5., -6., -7., -8., -9., -10. };
    plotStruct testPlot;
    testPlot.area = area;
    testPlot.cr = cr;
    testPlot.width = width;
    testPlot.height = height;
    testPlot.data = yData;
    testPlot.x0 = 1.0;
    testPlot.dx = 1.0;
    testPlot.Npts = 10;
    testPlot.axisColor = white;
    testPlot.color = magenta;
    testPlot.xLabel = "x label";
    testPlot.yLabel = "y label";
    testPlot.unitY = 1.0;
    testPlot.makeSVG = FALSE;
    LwePlot2d(&testPlot);
    context = gtk_widget_get_style_context(GTK_WIDGET(area));
}

void drawField2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    GtkStyleContext* context;
    LweColor black(0, 0, 0, 0);
    LweColor white(0.8, 0.8, 0.8, 0);
    LweColor magenta(1, 0, 1, 0);

    double xData[10] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. };
    double yData[10] = { 1., 2., 3., 4., 5., -6., -7., -8., -9., -10. };
    plotStruct testPlot;
    testPlot.area = area;
    testPlot.cr = cr;
    testPlot.width = width;
    testPlot.height = height;
    testPlot.data = yData;
    testPlot.x0 = 1.0;
    testPlot.dx = 1.0;
    testPlot.Npts = 10;
    testPlot.axisColor = white;
    testPlot.color = magenta;
    testPlot.xLabel = "x label";
    testPlot.yLabel = "y label";
    testPlot.unitY = 1.0;
    testPlot.makeSVG = FALSE;
    LwePlot2d(&testPlot);
    context = gtk_widget_get_style_context(GTK_WIDGET(area));
}

void drawSpectrum1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    GtkStyleContext* context;
    LweColor black(0, 0, 0, 0);
    LweColor white(0.8, 0.8, 0.8, 0);
    LweColor magenta(1, 0, 1, 0);

    double xData[10] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. };
    double yData[10] = { 1., 2., 3., 4., 5., -6., -7., -8., -9., -10. };
    plotStruct testPlot;
    testPlot.area = area;
    testPlot.cr = cr;
    testPlot.width = width;
    testPlot.height = height;
    testPlot.data = yData;
    testPlot.x0 = 1.0;
    testPlot.dx = 1.0;
    testPlot.Npts = 10;
    testPlot.axisColor = white;
    testPlot.color = magenta;
    testPlot.xLabel = "x label";
    testPlot.yLabel = "y label";
    testPlot.unitY = 1.0;
    testPlot.makeSVG = FALSE;
    LwePlot2d(&testPlot);
    context = gtk_widget_get_style_context(GTK_WIDGET(area));
}

void drawSpectrum2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    GtkStyleContext* context;
    LweColor black(0, 0, 0, 0);
    LweColor white(0.8, 0.8, 0.8, 0);
    LweColor magenta(1, 0, 1, 0);

    double xData[10] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. };
    double yData[10] = { 1., 2., 3., 4., 5., -6., -7., -8., -9., -10. };
    plotStruct testPlot;
    testPlot.area = area;
    testPlot.cr = cr;
    testPlot.width = width;
    testPlot.height = height;
    testPlot.data = yData;
    testPlot.x0 = 1.0;
    testPlot.dx = 1.0;
    testPlot.Npts = 10;
    testPlot.axisColor = white;
    testPlot.color = magenta;
    testPlot.xLabel = "x label";
    testPlot.yLabel = "y label";
    testPlot.unitY = 1.0;
    testPlot.makeSVG = FALSE;
    LwePlot2d(&testPlot);
    context = gtk_widget_get_style_context(GTK_WIDGET(area));
}

int linearRemapDoubleToFloat(double* A, int nax, int nay, float* B, int nbx, int nby) {
    int i, j;
    int nx0, ny0;
    int Ni, Nj;
    for (i = 0; i < nbx; ++i) {
        Ni = (int)(i * (nax / (float)nbx));
        nx0 = nay * minN(Ni, nax);
        for (j = 0; j < nby; ++j) {
            Nj = (int)((j * (nay / (float)nby)));
            ny0 = minN(nay, Nj);
            B[i * nby + j] = (float)A[ny0 + nx0];
        }
    }
    return 0;
}

void drawTimeImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    int Nx = 64;
    int Ny = 32;
    double dataC[64 * 64] = { 0 };
    for (int i = 0; i < 64; i++) {
        for (int j = 0; j < 64; j++) {
            dataC[i + 64 * j] = 0.1 * j;
        }
    }
    float* scaledData = new float[width * height]();
    linearRemapDoubleToFloat(dataC, Nx, Ny, scaledData, width, height);
    drawArrayAsBitmap(cr, width, height, scaledData, 3);
    delete[] scaledData;
}

void drawTimeImage2(GtkDrawingArea* area,
    cairo_t* cr,
    int             width,
    int             height,
    gpointer        data) {}

void drawFourierImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {}
void drawFourierImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {}


int main(int argc, char** argv) {
    GtkApplication* app;
    int status;

    app = gtk_application_new("lweTest.wtf", G_APPLICATION_FLAGS_NONE);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    status = g_application_run(G_APPLICATION(app), argc, argv);
    g_object_unref(app);

    return status;
}