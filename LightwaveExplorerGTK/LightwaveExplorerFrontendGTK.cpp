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
#define MAX_SIMULATIONS 4096
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
simulationParameterSet* activeSetPtr;       // Main structure containing simulation parameters and pointers
crystalEntry* crystalDatabasePtr;           // Crystal info database
size_t progressCounter = 0;

void updateDisplay();

class mainGui {
    bool isActivated;
    bool queueUpdate;
    std::thread threadPoolSim[3];
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
    bool saveSVG = FALSE;
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
    cudaGPUCount(0), syclGPUCount(0), saveSVG(0), queueUpdate(0){}
    ~mainGui() {}
    void requestPlotUpdate() {
        queueUpdate = TRUE;
    }
    void applyUpdate() {
        if (queueUpdate) {
            queueUpdate = FALSE;
            for (int i = 0; i < 8; ++i) {
                drawBoxes[i].queueDraw();
            }
        }
    }
    void activate(GtkApplication* app) {
        activeSetPtr = (simulationParameterSet*)calloc(MAX_SIMULATIONS, sizeof(simulationParameterSet));
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
        drawBoxes[4].setDrawingFunction(drawFourierImage1);

        drawBoxes[5].init(window.parentHandle(2), plotWidth, plotHeight, plotWidth, plotHeight);
        drawBoxes[5].setDrawingFunction(drawFourierImage2);

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


        pulldowns[4].addElement(_T("2D Cartesian"));
        pulldowns[4].addElement(_T("3D radial symm."));
        pulldowns[4].addElement(_T("3D"));
        pulldowns[4].init(parentHandle, textCol2a, 7, 2 * textWidth, 1);

        char batchModeNames[36][64] = {
        "none",
        "01: Energy 1",
        "02: Energy 2",
        "03: Frequency 1",
        "04: Frequency 2",
        "05: Bandwidth 1",
        "06: Bandwidth 2",
        "07: CEP 1",
        "08: CEP 2",
        "09: Delay 1",
        "10: Delay 2",
        "11: GDD 1",
        "12: GDD 2",
        "13: TOD 1",
        "14: TOD 2",
        "15: Thickness 1",
        "16: Thickness 2",
        "17: Beamwaist 1",
        "18: Beamwaist 2",
        "19: x offset 1",
        "20: x offset 2",
        "21: z offset 1",
        "22: z offset 2",
        "23: NC angle 1",
        "24: NC angle 2",
        "25: Polarization 1",
        "26: Polarization 2",
        "27: Circularity 1",
        "28: Circularity 2",
        "29: Crystal Theta",
        "30: Crystal Phi",
        "31: NL absorption",
        "32: Gamma",
        "33: Eff. mass",
        "34: Thickness",
        "35: dz",
        };

        for (int i = 0; i < 36; ++i) {
            pulldowns[5].addElement(batchModeNames[i]);
            pulldowns[6].addElement(batchModeNames[i]);
        }
        pulldowns[5].init(parentHandle, textCol2a, 8, 2 * textWidth, 1);
        pulldowns[6].init(parentHandle, textCol2a, 9, 2 * textWidth, 1);





        int openMPposition = 0;
        char A[128] = { 0 };
        checkLibraryAvailability();
        CUDAavailable = TRUE;
        if (CUDAavailable) {
            pulldowns[7].addElement("CUDA");
            pulldowns[8].addElement("CUDA");
            //swprintf_s(A, MAX_LOADSTRING, L"CUDA");
            //SendMessage(maingui.pdPrimaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
            //SendMessage(maingui.pdSecondaryQueue, (UINT)CB_ADDSTRING, (WPARAM)0, (LPARAM)A);
            openMPposition++;
            memset(&A, 0, sizeof(A));
            for (int i = 1; i < cudaGPUCount; ++i) {
                snprintf(A, 128, "CUDA %i", i);
                pulldowns[7].addElement(A);
                pulldowns[8].addElement(A);
                memset(&A, 0, sizeof(A));
                openMPposition++;
            }
        }
        if (SYCLavailable) {
            snprintf(A, 128, "SYCL");
            pulldowns[7].addElement(A);
            pulldowns[8].addElement(A);
            memset(&A, 0, sizeof(A));
            openMPposition++;
            if (syclGPUCount > 0) {
                snprintf(A, 128, "SYCL cpu");
                pulldowns[7].addElement(A);
                pulldowns[8].addElement(A);
                memset(&A, 0, sizeof(A));
                openMPposition++;
                snprintf(A, 128, "SYCL gpu");
                pulldowns[7].addElement(A);
                pulldowns[8].addElement(A);
                memset(&A, 0, sizeof(A));
                openMPposition++;
            }
        }

        snprintf(A, 128, "OpenMP");
        pulldowns[7].addElement(A);
        pulldowns[8].addElement(A);
        memset(&A, 0, sizeof(A));

        pulldowns[7].init(window.parentHandle(5), labelWidth, 0, buttonWidth, 1);
        pulldowns[8].init(window.parentHandle(5), 2 * labelWidth + buttonWidth, 0, buttonWidth, 1);
        textBoxes[52].init(window.parentHandle(5), 2 * labelWidth + 2 * buttonWidth, 0, 1, 1);


        sequence.init(parentHandle, buttonCol1, 13, colWidth, 6);
        fitCommand.init(parentHandle, buttonCol1, 21, colWidth, 4);

        buttons[0].init(_T("Run"), parentHandle, buttonCol3, 19, buttonWidth, 1, launchRunThread);
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
        
        textBoxes[44].setLabel(-labelWidth, 0, _T("Batch end"));
        textBoxes[46].setLabel(-labelWidth, 0, _T("Batch steps"));
        pulldowns[4].setLabel(-labelWidth, 0, _T("Propagation"));
        pulldowns[5].setLabel(-labelWidth, 0, _T("Batch mode"));
        pulldowns[6].setLabel(-labelWidth, 0, _T("Batch mode 2"));
        pulldowns[7].setLabel(-labelWidth, 0, _T("Run type:"));
        pulldowns[8].setLabel(-labelWidth, 0, _T(" Offload type:"));
        fitCommand.setLabel(0, -1, _T("Fitting:"));
        sequence.setLabel(0, -1, _T("Sequence:"));
        filePaths[3].overwritePrint("TestFile");
        window.present();

        

        //read the crystal database
        crystalDatabasePtr = (crystalEntry*)calloc(MAX_LOADSTRING, sizeof(crystalEntry));
        wchar_t materialString[128] = { 0 };
        if (crystalDatabasePtr != NULL) {
            //GetCurrentDirectory(MAX_LOADSTRING - 1, programDirectory);
            readCrystalDatabase(crystalDatabasePtr);
            console.wPrint(L"Material database has %i entries:\n", (*crystalDatabasePtr).numberOfEntries);
            for (int i = 0; i < (*crystalDatabasePtr).numberOfEntries; ++i) {
                console.wPrint(L"%2.2i: %s", i, crystalDatabasePtr[i].crystalNameW);
                memset(materialString, 0, 128 * sizeof(wchar_t));
                swprintf(materialString, 128, L"%2.2i: %s", i, crystalDatabasePtr[i].crystalNameW);
                pulldowns[3].addElementW(materialString);
            }
        }
        pulldowns[3].init(parentHandle, textCol2a, 0, 2 * textWidth, 1);
        pulldowns[3].setLabel(-labelWidth, 0, _T("Material"));


        pulldowns[9].addElement(_T("Cobra 1x R5000"));
        pulldowns[9].addElement(_T("Cobra 2x R5000"));
        pulldowns[9].addElement(_T("Cobra 1x V100"));
        pulldowns[9].addElement(_T("Cobra 1x V100"));
        pulldowns[9].addElement(_T("Raven 1x A100"));
        pulldowns[9].addElement(_T("Raven 2x A100"));
        pulldowns[9].addElement(_T("Raven 4x A100"));
        pulldowns[9].init(parentHandle, 0, 24, 2 * buttonWidth, 1);

        char defaultFilename[] = "DefaultValues.ini";
        readInputParametersFile(activeSetPtr, crystalDatabasePtr, defaultFilename);
        setInterfaceValuesToActiveValues();
        g_timeout_add(100, G_SOURCE_FUNC(updateDisplay), NULL);
    }
};
mainGui theGui;

char programDirectory[MAX_LOADSTRING];     // Program working directory (useful if the crystal database has to be reloaded)

///////////////////.
//Definitions over
///////////////////.

void updateDisplay() {
    theGui.console.updateFromBuffer();
    theGui.applyUpdate();
}

void setInterfaceValuesToActiveValues(){
	int i = 0;
    pulse* t = &(*activeSetPtr).pulse1;
    for (int k = 0; k < 2; ++k) {
        theGui.textBoxes[i++].setToDouble((*t).energy);
        theGui.textBoxes[i++].setToDouble(1e-12 * (*t).frequency);
        theGui.textBoxes[i++].setToDouble(1e-12 * (*t).bandwidth);
        theGui.textBoxes[i++].setToDouble((*t).sgOrder);
        theGui.textBoxes[i++].setToDouble(PI * (*t).cep);
        theGui.textBoxes[i++].setToDouble(1e15 * (*t).delay);
        theGui.textBoxes[i++].setToDouble(1e30 * (*t).gdd);
        theGui.textBoxes[i++].setToDouble(1e45 * (*t).tod);
        theGui.textBoxes[i++].setToDouble((*t).phaseMaterial);
        theGui.textBoxes[i++].setToDouble(1e6 * (*t).phaseMaterialThickness);
        theGui.textBoxes[i++].setToDouble(1e6 * (*t).beamwaist);
        theGui.textBoxes[i++].setToDouble(1e6 * (*t).x0);
        theGui.textBoxes[i++].setToDouble(1e6 * (*t).z0);
        theGui.textBoxes[i++].setToDouble(RAD2DEG * (*t).beamAngle);
        theGui.textBoxes[i++].setToDouble(RAD2DEG * (*t).polarizationAngle);
        theGui.textBoxes[i++].setToDouble((*t).circularity);
        t = &(*activeSetPtr).pulse2;
    }

    theGui.pulldowns[0].setValue((*activeSetPtr).pulse1FileType);
    theGui.pulldowns[1].setValue((*activeSetPtr).pulse2FileType);
    theGui.pulldowns[2].setValue((*activeSetPtr).fittingMode);
    theGui.pulldowns[3].setValue((*activeSetPtr).materialIndex);
    
    theGui.textBoxes[i++].setToDouble(RAD2DEG * asin(sin((*activeSetPtr).crystalTheta)));
    theGui.textBoxes[i++].setToDouble(RAD2DEG * asin(sin((*activeSetPtr).crystalPhi)));
    theGui.textBoxes[i++].setToDouble((*activeSetPtr).nonlinearAbsorptionStrength);
    theGui.textBoxes[i++].setToDouble((*activeSetPtr).bandGapElectronVolts);
    theGui.textBoxes[i++].setToDouble(1e-12 * (*activeSetPtr).drudeGamma);
    theGui.textBoxes[i++].setToDouble((*activeSetPtr).effectiveMass);

    if (!(*activeSetPtr).is3D) {
        theGui.textBoxes[i++].setToDouble(1e6 * (*activeSetPtr).spatialWidth);
    }
    else if ((*activeSetPtr).spatialHeight == (*activeSetPtr).spatialWidth
        || (*activeSetPtr).spatialHeight == 1 || (*activeSetPtr).spatialHeight == 0) {
        theGui.textBoxes[i++].setToDouble(1e6 * (*activeSetPtr).spatialWidth);
    }
    else {
        theGui.textBoxes[i++].overwritePrint("%i;%i", (int)(1e6 * (*activeSetPtr).spatialWidth), (int)(1e6 * (*activeSetPtr).spatialHeight));
    }
    theGui.textBoxes[i++].setToDouble(1e6*(*activeSetPtr).rStep);
    theGui.textBoxes[i++].setToDouble(1e15 * (*activeSetPtr).timeSpan);
    theGui.textBoxes[i++].setToDouble(1e15 * (*activeSetPtr).tStep);
    theGui.textBoxes[i++].setToDouble(1e6 * (*activeSetPtr).crystalThickness);
    theGui.textBoxes[i++].setToDouble(1e6 * (*activeSetPtr).propagationStep);
    theGui.textBoxes[i++].setToDouble((*activeSetPtr).batchDestination);
    theGui.textBoxes[i++].setToDouble((*activeSetPtr).batchDestination2);
    theGui.textBoxes[i++].setToDouble((double)(*activeSetPtr).Nsims);
    theGui.textBoxes[i++].setToDouble((double)(*activeSetPtr).Nsims2);
    theGui.pulldowns[4].setValue((*activeSetPtr).symmetryType);
    theGui.pulldowns[5].setValue((*activeSetPtr).batchIndex);
    theGui.pulldowns[6].setValue((*activeSetPtr).batchIndex2);

    if (strlen((*activeSetPtr).sequenceString) > 6) {
        formatSequence((*activeSetPtr).sequenceString, MAX_LOADSTRING);
        theGui.sequence.overwritePrint((*activeSetPtr).sequenceString);
        removeCharacterFromString((*activeSetPtr).sequenceString, MAX_LOADSTRING, '\r');
        removeCharacterFromString((*activeSetPtr).sequenceString, MAX_LOADSTRING, '\n');
    }
    stripLineBreaks((*activeSetPtr).field1FilePath);
    if (!strncmp((*activeSetPtr).field1FilePath, "None",4)) theGui.filePaths[0].overwritePrint((*activeSetPtr).field1FilePath);
    if (!strcmp((*activeSetPtr).field2FilePath, "None")) theGui.filePaths[1].overwritePrint((*activeSetPtr).field2FilePath);
    if (!strcmp((*activeSetPtr).fittingPath, "None")) theGui.filePaths[2].overwritePrint((*activeSetPtr).fittingPath);

    if (!((*activeSetPtr).fittingString[0] == 'N')) {
        insertLineBreaksAfterSemicolons((*activeSetPtr).fittingString, MAX_LOADSTRING);
        theGui.fitCommand.overwritePrint((*activeSetPtr).fittingString);
        removeCharacterFromString((*activeSetPtr).fittingString, MAX_LOADSTRING, '\r');
        removeCharacterFromString((*activeSetPtr).fittingString, MAX_LOADSTRING, '\n');
    }
}


void readParametersFromInterface() {
    int i = 0;
    pulse* t = &(*activeSetPtr).pulse1;
    for (int k = 0; k < 2; ++k) {
        theGui.textBoxes[i++].valueToPointer(&(*t).energy);
        theGui.textBoxes[i++].valueToPointer(1e12, &(*t).frequency);
        theGui.textBoxes[i++].valueToPointer(1e12, &(*t).bandwidth);
        theGui.textBoxes[i++].valueToPointer(&(*t).sgOrder);
        theGui.textBoxes[i++].valueToPointer(1.0/PI, &(*t).cep);
        theGui.textBoxes[i++].valueToPointer(1e-15, &(*t).delay);
        theGui.textBoxes[i++].valueToPointer(1e-30, &(*t).gdd);
        theGui.textBoxes[i++].valueToPointer(1e-45, &(*t).tod);
        theGui.textBoxes[i++].valueToPointer(&(*t).phaseMaterial);
        theGui.textBoxes[i++].valueToPointer(1e-6, &(*t).phaseMaterialThickness);
        theGui.textBoxes[i++].valueToPointer(1e-6, &(*t).beamwaist);
        theGui.textBoxes[i++].valueToPointer(1e-6, &(*t).x0);
        theGui.textBoxes[i++].valueToPointer(1e-6, &(*t).z0);
        theGui.textBoxes[i++].valueToPointer(DEG2RAD, &(*t).beamAngle);
        theGui.textBoxes[i++].valueToPointer(DEG2RAD, &(*t).polarizationAngle);
        theGui.textBoxes[i++].valueToPointer(&(*t).circularity);
        t = &(*activeSetPtr).pulse2;
    }
    
    (*activeSetPtr).pulse1FileType = theGui.pulldowns[0].getValue();
    (*activeSetPtr).pulse2FileType = theGui.pulldowns[1].getValue();
    (*activeSetPtr).fittingMode = theGui.pulldowns[2].getValue();
    (*activeSetPtr).materialIndex = theGui.pulldowns[3].getValue();

    theGui.textBoxes[i++].valueToPointer(DEG2RAD, &(*activeSetPtr).crystalTheta);
    theGui.textBoxes[i++].valueToPointer(DEG2RAD,  &(*activeSetPtr).crystalPhi);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).nonlinearAbsorptionStrength);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).bandGapElectronVolts);
    theGui.textBoxes[i++].valueToPointer(1e12, &(*activeSetPtr).drudeGamma);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).effectiveMass);


    theGui.textBoxes[i++].valueToTwoPointers(1e-6, &(*activeSetPtr).spatialWidth, &(*activeSetPtr).spatialHeight);

    theGui.textBoxes[i++].valueToPointer(1e-6, &(*activeSetPtr).rStep);
    theGui.textBoxes[i++].valueToPointer(1e-15, &(*activeSetPtr).timeSpan);
    theGui.textBoxes[i++].valueToPointer(1e-15, &(*activeSetPtr).tStep);
    theGui.textBoxes[i++].valueToPointer(1e-6, &(*activeSetPtr).crystalThickness);
    theGui.textBoxes[i++].valueToPointer(1e-6, &(*activeSetPtr).propagationStep);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).batchDestination);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).batchDestination2);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).Nsims);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).Nsims2);

    (*activeSetPtr).symmetryType = theGui.pulldowns[4].getValue();
    (*activeSetPtr).batchIndex = theGui.pulldowns[5].getValue();
    (*activeSetPtr).batchIndex2 = theGui.pulldowns[6].getValue();

    theGui.textBoxes[52].valueToPointer(&(*activeSetPtr).NsimsCPU);

    char noneString[] = "None";
    memset((*activeSetPtr).sequenceString, 0, MAX_LOADSTRING);
    theGui.sequence.copyBuffer((*activeSetPtr).sequenceString, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).sequenceString, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).sequenceString, noneString);
    }
    else {
        stripWhiteSpace((*activeSetPtr).sequenceString);
    }
    
    memset((*activeSetPtr).fittingString, 0, MAX_LOADSTRING);
    theGui.fitCommand.copyBuffer((*activeSetPtr).fittingString, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).fittingString, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).fittingString, noneString);
    }
    else {
        stripLineBreaks((*activeSetPtr).fittingString);
    }
    
    memset((*activeSetPtr).field1FilePath, 0, MAX_LOADSTRING);
    theGui.filePaths[0].copyBuffer((*activeSetPtr).field1FilePath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).field1FilePath, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).field1FilePath, noneString);
    }
    stripLineBreaks((*activeSetPtr).field1FilePath);
    memset((*activeSetPtr).field2FilePath, 0, MAX_LOADSTRING);
    theGui.filePaths[1].copyBuffer((*activeSetPtr).field2FilePath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).field2FilePath, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).field2FilePath, noneString);
    }
    stripLineBreaks((*activeSetPtr).field2FilePath);
    memset((*activeSetPtr).outputBasePath, 0, MAX_LOADSTRING);
    theGui.filePaths[3].copyBuffer((*activeSetPtr).outputBasePath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).outputBasePath, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).outputBasePath, noneString);
    }
    stripLineBreaks((*activeSetPtr).outputBasePath);
    
    memset((*activeSetPtr).fittingPath, 0, MAX_LOADSTRING);
    theGui.filePaths[2].copyBuffer((*activeSetPtr).fittingPath, MAX_LOADSTRING);
    if (strnlen_s((*activeSetPtr).fittingPath, MAX_LOADSTRING) == 0) {
        strcpy((*activeSetPtr).fittingPath, noneString);
    }
    stripLineBreaks((*activeSetPtr).fittingPath);

    //derived parameters and cleanup:
    (*activeSetPtr).sellmeierType = 0;
    (*activeSetPtr).axesNumber = 0;
    (*activeSetPtr).Ntime = (size_t)(MIN_GRIDDIM * ceil((*activeSetPtr).timeSpan / (MIN_GRIDDIM * (*activeSetPtr).tStep)));

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
    //printC(L"(x,y) = %lli, %lli\r\n", (*activeSetPtr).Nspace, (*activeSetPtr).Nspace2);
    (*activeSetPtr).Nfreq = (*activeSetPtr).Ntime / 2 + 1;
    (*activeSetPtr).NgridC = (*activeSetPtr).Nfreq * (*activeSetPtr).Nspace * (*activeSetPtr).Nspace2;
    (*activeSetPtr).Ngrid = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * (*activeSetPtr).Nspace2;
    (*activeSetPtr).kStep = TWOPI / ((*activeSetPtr).Nspace * (*activeSetPtr).rStep);
    (*activeSetPtr).fStep = 1.0 / ((*activeSetPtr).Ntime * (*activeSetPtr).tStep);
    (*activeSetPtr).Npropagation = (size_t)round((*activeSetPtr).crystalThickness / (*activeSetPtr).propagationStep);

    (*activeSetPtr).isCylindric = (*activeSetPtr).symmetryType == 1;
    if ((*activeSetPtr).isCylindric) {
        (*activeSetPtr).pulse1.x0 = 0;
        (*activeSetPtr).pulse2.x0 = 0;
        (*activeSetPtr).pulse1.beamAngle = 0;
        (*activeSetPtr).pulse2.beamAngle = 0;
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
}


int insertLineBreaksAfterSemicolons(char* cString, size_t N) {
    size_t i = 0;
    while (i < N - 1) {
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
        }
        i++;
    }
    return 0;
}

int indentForDepth(char* cString, size_t N) {
    size_t i = 0;
    int depth = 0;
    while (i < N - 1) {
        if (cString[i] == '{') depth++;
        if (cString[i + 1] == '}') depth--;

        if (cString[i] == ';' || ((cString[i] == ')' && cString[i + 1] != ';')) || cString[i] == '{') {
            for (int k = 0; k < 2 * depth; k++) {
                if (cString[i + 1] != ' ' && cString[i + 2] != ' ') {
                    memmove(&cString[i + 3], &cString[i + 1], N - i - 3);
                }
                else if (cString[i + 1] != ' ') {
                    memmove(&cString[i + 3], &cString[i + 1], N - i - 3);
                }
                else if (cString[i + 1] == ' ') {
                    memmove(&cString[i + 3], &cString[i + 2], N - i - 2);
                }
                cString[i + 1] = ' ';
                cString[i + 2] = ' ';
                i += 2;
            }
        }
        i++;
    }
    return 0;
}

int insertLineBreaksAfterClosingParenthesis(char* cString, size_t N) {
    size_t i = 0;
    while (i < N - 1) {
        if (cString[i] == ')' && cString[i + 1] != ';' && cString[i + 1] != '{') {
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
        }
        i++;
    }
    return 0;
}

int insertLineBreaksAfterClosingAngle(char* cString, size_t N) {
    size_t i = 0;
    while (i < N - 1) {
        if (cString[i] == '>') {
            memmove(&cString[i + 3], &cString[i + 1], N - i - 3);
            cString[i + 1] = '\r';
            cString[i + 2] = '\n';
            i += 2;
        }
        i++;
    }
    return 0;
}

int insertLineBreaksAfterCurlyBraces(char* cString, size_t N) {
    size_t i = 0;
    while (i < N - 1) {
        if (cString[i] == '{' || cString[i] == '}') {
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
        }
        i++;
    }
    return 0;
}
int formatSequence(char* cString, size_t N) {
    indentForDepth(cString, N);
    insertLineBreaksAfterClosingAngle(cString, N);
    insertLineBreaksAfterClosingParenthesis(cString, N);
    insertLineBreaksAfterSemicolons(cString, N);
    insertLineBreaksAfterCurlyBraces(cString, N);
    return 0;
}

int freeSemipermanentGrids() {
    isGridAllocated = FALSE;
    delete[](*activeSetPtr).ExtOut;
    delete[](*activeSetPtr).EkwOut;
    delete[](*activeSetPtr).totalSpectrum;
    return 0;
}

void checkLibraryAvailability() {
    theGui.console.cPrint("\r\n");
    __try {
        HRESULT hr = __HrLoadAllImportsForDll("cufft64_10.dll");
        if (SUCCEEDED(hr)) {
            CUDAavailable = TRUE;
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
    if (!isGridAllocated) return;
    plotStruct sPlot;
    char* svgBuf = NULL;
    char* svgFilename = NULL;
    FILE* svgFile;

    if (theGui.saveSVG) {
        svgFilename = new char[MAX_LOADSTRING]();
        svgBuf = new char[1024 * 1024]();
    }
    bool logPlot = FALSE;
    if (FALSE) {
        logPlot = FALSE;
    }
    size_t simIndex = 0;// theGui.plotSlider.getIntValue();
    if (simIndex < 0) {
        simIndex = 0;
    }
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }

    bool forceX = FALSE;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = TRUE;
    }
    bool forceY = FALSE;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = TRUE;
    }
    bool overlayTotal = FALSE;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = TRUE;
    }
    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

    size_t cubeMiddle = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * ((*activeSetPtr).Nspace2 / 2);

    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = (*activeSetPtr).tStep / 1e-15;
    sPlot.x0 = -((sPlot.dx * (*activeSetPtr).Ntime) / 2 - sPlot.dx / 2);
    sPlot.data = &(*activeSetPtr).ExtOut[simIndex * (*activeSetPtr).Ngrid * 2 + cubeMiddle + (*activeSetPtr).Ntime * (*activeSetPtr).Nspace / 2];
    sPlot.Npts = (*activeSetPtr).Ntime;
    sPlot.color = LweColor(0, 1, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ey (GV/m)";
    sPlot.unitY = 1e9;
    sPlot.makeSVG = FALSE; // theGui.saveSVG;
    sPlot.svgString = svgBuf;
    sPlot.svgBuffer = 1024 * 1024;

    LwePlot2d(&sPlot);

    if (theGui.saveSVG) {
        memset(svgFilename, 0, MAX_LOADSTRING);
        theGui.filePaths[3].copyBuffer(svgFilename, MAX_LOADSTRING);
        strcat_s(svgFilename, MAX_LOADSTRING, "_Ey.svg");
        if (fopen_s(&svgFile, svgFilename, "w")) return;
        fprintf(svgFile, "%s", svgBuf);
        fclose(svgFile);
        delete[] svgBuf;
        delete[] svgFilename;
    }
}

void drawField2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) return;
    plotStruct sPlot;
    char* svgBuf = NULL;
    char* svgFilename = NULL;
    FILE* svgFile;

    if (theGui.saveSVG) {
        svgFilename = new char[MAX_LOADSTRING]();
        svgBuf = new char[1024 * 1024]();
    }
    bool logPlot = FALSE;
    if (FALSE) {
        logPlot = FALSE;
    }
    size_t simIndex = 0;// theGui.plotSlider.getIntValue();
    if (simIndex < 0) {
        simIndex = 0;
    }
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }

    bool forceX = FALSE;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = TRUE;
    }
    bool forceY = FALSE;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = TRUE;
    }
    bool overlayTotal = FALSE;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = TRUE;
    }
    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

    size_t cubeMiddle = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * ((*activeSetPtr).Nspace2 / 2);

    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = (*activeSetPtr).tStep / 1e-15;
    sPlot.x0 = -((sPlot.dx * (*activeSetPtr).Ntime) / 2 - sPlot.dx / 2);
    sPlot.data = &(*activeSetPtr).ExtOut[(*activeSetPtr).Ngrid + simIndex * (*activeSetPtr).Ngrid * 2 + cubeMiddle + (*activeSetPtr).Ntime * (*activeSetPtr).Nspace / 2];
    sPlot.Npts = (*activeSetPtr).Ntime;
    sPlot.color = LweColor(1, 0, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ex (GV/m)";
    sPlot.unitY = 1e9;
    sPlot.makeSVG = FALSE; // theGui.saveSVG;
    sPlot.svgString = svgBuf;
    sPlot.svgBuffer = 1024 * 1024;

    LwePlot2d(&sPlot);

    if (theGui.saveSVG) {
        memset(svgFilename, 0, MAX_LOADSTRING);
        theGui.filePaths[3].copyBuffer(svgFilename, MAX_LOADSTRING);
        strcat_s(svgFilename, MAX_LOADSTRING, "_Ex.svg");
        if (fopen_s(&svgFile, svgFilename, "w")) return;
        fprintf(svgFile, "%s", svgBuf);
        fclose(svgFile);
        delete[] svgBuf;
        delete[] svgFilename;
    }
}

void drawSpectrum1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) return;
    plotStruct sPlot;
    char* svgBuf = NULL;
    char* svgFilename = NULL;
    FILE* svgFile;

    if (theGui.saveSVG) {
        svgFilename = new char[MAX_LOADSTRING]();
        svgBuf = new char[1024 * 1024]();
    }
    bool logPlot = FALSE;
    if (FALSE) {
        logPlot = FALSE;
    }
    size_t simIndex = 0;// theGui.plotSlider.getIntValue();
    if (simIndex < 0) {
        simIndex = 0;
    }
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }

    bool forceX = FALSE;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = TRUE;
    }
    bool forceY = FALSE;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = TRUE;
    }
    bool overlayTotal = FALSE;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = TRUE;
    }
    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

    size_t cubeMiddle = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * ((*activeSetPtr).Nspace2 / 2);

    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = (*activeSetPtr).fStep / 1e12;
    sPlot.data = &(*activeSetPtr).totalSpectrum[simIndex * 3 * (*activeSetPtr).Nfreq];
    sPlot.Npts = (*activeSetPtr).Nfreq;
    sPlot.logScale = logPlot;
    sPlot.forceYmin = (logPlot || forceY);
    sPlot.forcedYmin = -4.0;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = LweColor(0.5, 0, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Sy (W/Hz)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    if (overlayTotal) {
        sPlot.data2 = &(*activeSetPtr).totalSpectrum[(2 + simIndex * 3) * (*activeSetPtr).Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    sPlot.makeSVG = FALSE; // theGui.saveSVG;
    sPlot.svgString = svgBuf;
    sPlot.svgBuffer = 1024 * 1024;

    LwePlot2d(&sPlot);

    if (theGui.saveSVG) {
        memset(svgFilename, 0, MAX_LOADSTRING);
        theGui.filePaths[3].copyBuffer(svgFilename, MAX_LOADSTRING);
        strcat_s(svgFilename, MAX_LOADSTRING, "_Sx.svg");
        if (fopen_s(&svgFile, svgFilename, "w")) return;
        fprintf(svgFile, "%s", svgBuf);
        fclose(svgFile);
        delete[] svgBuf;
        delete[] svgFilename;
    }
}

void drawSpectrum2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) return;
    plotStruct sPlot;
    char* svgBuf = NULL;
    char* svgFilename = NULL;
    FILE* svgFile;

    if (theGui.saveSVG) {
        svgFilename = new char[MAX_LOADSTRING]();
        svgBuf = new char[1024 * 1024]();
    }
    bool logPlot = FALSE;
    if (FALSE) {
        logPlot = FALSE;
    }
    size_t simIndex = 0;// theGui.plotSlider.getIntValue();
    if (simIndex < 0) {
        simIndex = 0;
    }
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }

    bool forceX = FALSE;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = TRUE;
    }
    bool forceY = FALSE;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = TRUE;
    }
    bool overlayTotal = FALSE;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = TRUE;
    }
    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

    size_t cubeMiddle = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * ((*activeSetPtr).Nspace2 / 2);

    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = (*activeSetPtr).fStep / 1e12;
    sPlot.data = &(*activeSetPtr).totalSpectrum[(1 + simIndex * 3) * (*activeSetPtr).Nfreq];
    sPlot.Npts = (*activeSetPtr).Nfreq;
    sPlot.logScale = logPlot;
    sPlot.forceYmin = (logPlot || forceY);
    sPlot.forcedYmin = -4.0;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = LweColor(1, 0, 0.5, 0.0);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Sx (W/Hz)";
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    if (overlayTotal) {
        sPlot.data2 = &(*activeSetPtr).totalSpectrum[(2 + simIndex * 3) * (*activeSetPtr).Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    sPlot.makeSVG = FALSE; // theGui.saveSVG;
    sPlot.svgString = svgBuf;
    sPlot.svgBuffer = 1024 * 1024;

    LwePlot2d(&sPlot);

	if (theGui.saveSVG) {
		memset(svgFilename, 0, MAX_LOADSTRING);
		theGui.filePaths[3].copyBuffer(svgFilename, MAX_LOADSTRING);
		strcat_s(svgFilename, MAX_LOADSTRING, "_Sx.svg");
		if (fopen_s(&svgFile, svgFilename, "w")) return;
		fprintf(svgFile, "%s", svgBuf);
		fclose(svgFile);
		delete[] svgBuf;
		delete[] svgFilename;
	}
}

int linearRemapZToLogFloat(std::complex<double>* A, int nax, int nay, float* B, int nbx, int nby, double logMin) {
    int i, j;
    float A00;
    float f;

    int nx0, ny0;
    int Ni, Nj;
    for (i = 0; i < nbx; ++i) {
        f = i * (nax / (float)nbx);
        Ni = (int)f;
        nx0 = nay * min(Ni, nax);
        for (j = 0; j < nby; ++j) {
            f = (j * (nay / (float)nby));
            Nj = (int)f;
            ny0 = min(nay, Nj);
            A00 = (float)log10(cModulusSquared(A[ny0 + nx0]) + logMin);
            B[i * nby + j] = A00;
        }
    }
    return 0;
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

void imagePlot(imagePlotStruct* s) {

    int dx = (*s).width;
    int dy = (*s).height;

    size_t plotSize = (size_t)dx * (size_t)dy;
    float* plotarr2 = (float*)malloc(plotSize * sizeof(float));

    switch ((*s).dataType) {
    case 0:
        linearRemapDoubleToFloat((*s).data, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Ntime, plotarr2, (int)dy, (int)dx);
        break;
    case 1:
        std::complex<double> *shiftedFFT = (std::complex<double>*)malloc((*activeSetPtr).Nspace * (*activeSetPtr).Nfreq * sizeof(std::complex<double>));
        fftshiftD2Z((*s).complexData, shiftedFFT, (*activeSetPtr).Nfreq, (*activeSetPtr).Nspace);
        linearRemapZToLogFloat(shiftedFFT, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Nfreq, plotarr2, (int)dy, (int)dx, (*s).logMin);
        free(shiftedFFT);
        break;
    }
    drawArrayAsBitmap((*s).cr, (*s).width, (*s).height, plotarr2, (*s).colorMap);
    free(plotarr2);


}


void drawTimeImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) return;
    imagePlotStruct sPlot;
    char* svgBuf = NULL;
    char* svgFilename = NULL;
    FILE* svgFile;

    bool logPlot = FALSE;
    if (FALSE) {
        logPlot = FALSE;
    }
    size_t simIndex = 0;// theGui.plotSlider.getIntValue();
    if (simIndex < 0) {
        simIndex = 0;
    }
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }

    bool forceX = FALSE;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = TRUE;
    }
    bool forceY = FALSE;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = TRUE;
    }
    bool overlayTotal = FALSE;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = TRUE;
    }
    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

    size_t cubeMiddle = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * ((*activeSetPtr).Nspace2 / 2);

    sPlot.data =
        &(*activeSetPtr).ExtOut[simIndex * (*activeSetPtr).Ngrid * 2 + cubeMiddle];
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 0;
    sPlot.colorMap = 4;
    imagePlot(&sPlot);
}

void drawTimeImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) return;
    imagePlotStruct sPlot;
    char* svgBuf = NULL;
    char* svgFilename = NULL;
    FILE* svgFile;


    bool logPlot = FALSE;
    if (FALSE) {
        logPlot = FALSE;
    }
    size_t simIndex = 0;// theGui.plotSlider.getIntValue();
    if (simIndex < 0) {
        simIndex = 0;
    }
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }

    bool forceX = FALSE;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = TRUE;
    }
    bool forceY = FALSE;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = TRUE;
    }
    bool overlayTotal = FALSE;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = TRUE;
    }
    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

    size_t cubeMiddle = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * ((*activeSetPtr).Nspace2 / 2);

    sPlot.data =
        &(*activeSetPtr).ExtOut[(*activeSetPtr).Ngrid + simIndex * (*activeSetPtr).Ngrid * 2 + cubeMiddle];
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 0;
    sPlot.colorMap = 4;
    imagePlot(&sPlot);
}

void drawFourierImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) return;
    imagePlotStruct sPlot;
    char* svgBuf = NULL;
    char* svgFilename = NULL;
    FILE* svgFile;


    bool logPlot = FALSE;
    if (FALSE) {
        logPlot = FALSE;
    }
    size_t simIndex = 0;// theGui.plotSlider.getIntValue();
    if (simIndex < 0) {
        simIndex = 0;
    }
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }

    bool forceX = FALSE;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = TRUE;
    }
    bool forceY = FALSE;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = TRUE;
    }
    bool overlayTotal = FALSE;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = TRUE;
    }
    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

    size_t cubeMiddle = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * ((*activeSetPtr).Nspace2 / 2);

    sPlot.complexData =
        &(*activeSetPtr).EkwOut[simIndex * (*activeSetPtr).NgridC * 2];
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 1;
    sPlot.logMin = logPlotOffset;
    sPlot.colorMap = 3;
    imagePlot(&sPlot);
}

void drawFourierImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) return;
    imagePlotStruct sPlot;
    char* svgBuf = NULL;
    char* svgFilename = NULL;
    FILE* svgFile;


    bool logPlot = FALSE;
    if (FALSE) {
        logPlot = FALSE;
    }
    size_t simIndex = 0;// theGui.plotSlider.getIntValue();
    if (simIndex < 0) {
        simIndex = 0;
    }
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }

    bool forceX = FALSE;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = TRUE;
    }
    bool forceY = FALSE;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = TRUE;
    }
    bool overlayTotal = FALSE;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = TRUE;
    }
    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

    size_t cubeMiddle = (*activeSetPtr).Ntime * (*activeSetPtr).Nspace * ((*activeSetPtr).Nspace2 / 2);

    sPlot.complexData =
        &(*activeSetPtr).EkwOut[simIndex * (*activeSetPtr).NgridC * 2 + (*activeSetPtr).NgridC];
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 1;
    sPlot.logMin = logPlotOffset;
    sPlot.colorMap = 3;
    imagePlot(&sPlot);
}


void launchRunThread() {
    if(!isRunning) std::thread(simpleSimThread).detach();
}


void mainSimThread() {
    cancellationCalled = FALSE;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();
    HANDLE plotThread = NULL;
    DWORD hplotThread;
    HANDLE cpuThread = NULL;
    DWORD hCpuThread;

    if (isGridAllocated) {
        freeSemipermanentGrids();
    }
    memset(activeSetPtr, 0, sizeof(simulationParameterSet));
    readParametersFromInterface();
    if ((*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 > MAX_SIMULATIONS) {
        theGui.console.cPrint("Too many simulations in batch mode. Must be under %i total.\r\n", MAX_SIMULATIONS);
    }
    (*activeSetPtr).runType = 0;
    allocateGrids(activeSetPtr);
    isGridAllocated = TRUE;
    //setTrackbarLimitsToActiveSet();
    (*activeSetPtr).isFollowerInSequence = FALSE;
    (*activeSetPtr).crystalDatabase = crystalDatabasePtr;
    loadPulseFiles(activeSetPtr);

    if ((*activeSetPtr).sequenceString[0] != 'N') (*activeSetPtr).isInSequence = TRUE;

    configureBatchMode(activeSetPtr);
    int error = 0;
    //run the simulations
    isRunning = TRUE;
    progressCounter = 0;
    int pulldownSelection = theGui.pulldowns[8].getValue();
    auto sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
    auto normalFunction = &solveNonlinearWaveEquationCPU;
    int assignedGPU = 0;
    bool forceCPU = 0;
    int SYCLitems = 0;
    if (syclGPUCount == 0) {
        SYCLitems = (int)SYCLavailable;
    }
    else {
        SYCLitems = 3;
    }
    if (pulldownSelection < cudaGPUCount) {
        sequenceFunction = &solveNonlinearWaveEquationSequence;
        normalFunction = &solveNonlinearWaveEquation;
        assignedGPU = pulldownSelection;
    }
    else if (pulldownSelection == cudaGPUCount && SYCLavailable) {
        sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
        normalFunction = &solveNonlinearWaveEquationSYCL;
    }
    else if (pulldownSelection == cudaGPUCount + 1 && SYCLitems > 1) {
        forceCPU = 1;
        sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
        normalFunction = &solveNonlinearWaveEquationSYCL;
    }
    else if (pulldownSelection == cudaGPUCount + 2 && SYCLitems > 1) {
        assignedGPU = 1;
        sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
        normalFunction = &solveNonlinearWaveEquationSYCL;
    }
    else {
        sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
        normalFunction = &solveNonlinearWaveEquationCPU;
    }

    //cpuThread = CreateThread(NULL, 0, secondaryQueue, &activeSetPtr[(*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 - (*activeSetPtr).NsimsCPU], 0, &hCpuThread);
    for (int j = 0; j < ((*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 - (*activeSetPtr).NsimsCPU); ++j) {

        activeSetPtr[j].runningOnCPU = forceCPU;
        activeSetPtr[j].assignedGPU = assignedGPU;
        if ((*activeSetPtr).isInSequence) {
            error = sequenceFunction(&activeSetPtr[j]);

            if (activeSetPtr[j].memoryError != 0) {
                if (activeSetPtr[j].memoryError == -1) {
                    theGui.console.cPrint(_T("Not enough free GPU memory, sorry.\r\n"), activeSetPtr[j].memoryError);
                }
                else {
                    theGui.console.cPrint(_T("Warning: device memory error (%i).\r\n"), activeSetPtr[j].memoryError);
                }
            }
            if (error) break;
            if (!isPlotting) {
                (*activeSetPtr).plotSim = j;
                //plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);
                //if (plotThread != 0)CloseHandle(plotThread);
            }

        }
        else {
            error = normalFunction(&activeSetPtr[j]);
            if (activeSetPtr[j].memoryError != 0) {
                if (activeSetPtr[j].memoryError == -1) {
                    theGui.console.cPrint(_T("Not enough free GPU memory, sorry.\r\n"), activeSetPtr[j].memoryError);
                }
                else {
                    theGui.console.cPrint(_T("Warning: device memory error (%i).\r\n"), activeSetPtr[j].memoryError);
                }
            }

            if (error) break;
        }

        if (cancellationCalled) {
            theGui.console.cPrint(_T("Warning: series cancelled, stopping after %i simulations.\r\n"), j + 1);
            break;
        }

        if (!isPlotting) {
            (*activeSetPtr).plotSim = j;
            //plotThread = CreateThread(NULL, 0, drawSimPlots, activeSetPtr, 0, &hplotThread);
            //if (plotThread != NULL)CloseHandle(plotThread);
        }
    }

    if ((*activeSetPtr).NsimsCPU != 0 && cpuThread != 0) {
        //WaitForSingleObject(cpuThread, INFINITE);
        //CloseHandle(cpuThread);
    }
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    if (error == 13) {
        theGui.console.cPrint(
            "NaN detected in grid!\r\nTry using a larger spatial/temporal step\r\nor smaller propagation step.\r\nSimulation was cancelled.\r\n");
    }
    else {
        theGui.console.cPrint(_T("Finished after %8.4lf s. \r\n"), 1e-6 *
            (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
    }
    //saveDataSet(activeSetPtr, crystalDatabasePtr, (*activeSetPtr).outputBasePath, FALSE);
    deallocateGrids(activeSetPtr, FALSE);
    isRunning = FALSE;
}


void simpleSimThread() {
    cancellationCalled = FALSE;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();
    if (isGridAllocated) {
        freeSemipermanentGrids();
    }
    memset(activeSetPtr, 0, sizeof(simulationParameterSet));
    readParametersFromInterface();
    if ((*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 > MAX_SIMULATIONS) {
        theGui.console.threadPrint("Too many simulations in batch mode. Must be under %i total.\r\n", MAX_SIMULATIONS);
    }
    (*activeSetPtr).runType = 0;
    allocateGrids(activeSetPtr);
    isGridAllocated = TRUE;
    (*activeSetPtr).isFollowerInSequence = FALSE;
    (*activeSetPtr).crystalDatabase = crystalDatabasePtr;
    //loadPulseFiles(activeSetPtr);

    //if ((*activeSetPtr).sequenceString[0] != 'N') (*activeSetPtr).isInSequence = TRUE;

    configureBatchMode(activeSetPtr);
    int error = 0;
    //run the simulations
    isRunning = TRUE;
    progressCounter = 0;
    
    error = solveNonlinearWaveEquation(activeSetPtr);
    theGui.requestPlotUpdate();
    //Sleep(500);
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    if (error == 13) {
        theGui.console.threadPrint(
            "NaN detected in grid!\r\nTry using a larger spatial/temporal step\r\nor smaller propagation step.\r\nSimulation was cancelled.\r\n");
    }
    else {
        theGui.console.threadPrint(_T("Finished after %8.4lf s.\n"), 1e-6 *
            (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
    }
    
    //saveDataSet(activeSetPtr, crystalDatabasePtr, (*activeSetPtr).outputBasePath, FALSE);
    deallocateGrids(activeSetPtr, FALSE);
    isRunning = FALSE;
}

int main(int argc, char** argv) {
    GtkApplication* app;
    int status;

    app = gtk_application_new("lweTest.wtf", G_APPLICATION_FLAGS_NONE);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    status = g_application_run(G_APPLICATION(app), argc, argv);
    g_object_unref(app);

    return status;
}