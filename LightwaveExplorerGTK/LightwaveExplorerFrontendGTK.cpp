#include "LightwaveExplorerFrontendGTK.h"
#include <thread>
#include <chrono>
#include <locale>
#include <fstream>
#include "../LightwaveExplorerCoreCPU.h"
#include "../LightwaveExplorerCoreCounter.h"
#ifndef CPUONLY
#ifndef NOCUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <nvml.h>
#include "../LightwaveExplorerCore.cuh"
#endif
#ifdef _WIN32
#include "../LightwaveExplorerSYCL/LightwaveExplorerSYCL.h"
#else
#include "LightwaveExplorerDPCPPlib.h"
#endif
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#import<Cocoa/Cocoa.h>
#endif

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
size_t progressCounter = 0;
size_t totalSteps = 0;
#if defined _WIN32 || defined __linux__
const int interfaceThreads = maxN(1, std::thread::hardware_concurrency() / 2);
#else
const int interfaceThreads = std::thread::hardware_concurrency();
#endif
simulationParameterSet* activeSetPtr;       // Main structure containing simulation parameters and pointers
crystalEntry* crystalDatabasePtr;           // Crystal info database
void updateDisplay();

class mainGui {
    bool queueUpdate;
    bool queueSliderUpdate;
    bool queueSliderMove;
    int sliderTarget;
    std::thread threadPoolSim[3];
public:
    LweTextBox textBoxes[54];
    LweButton buttons[16];
    LweButton miniButtons[12];
    LweConsole console;
    LweConsole sequence;
    LweConsole fitCommand;
    LweTextBox filePaths[4];
    LwePulldown pulldowns[10];
    LweDrawBox drawBoxes[8];
    LweDrawBox progressBarBox;
    LweCheckBox checkBoxes[3];
    LweSlider plotSlider;
    LweWindow window;
    LweSpacer spacers[2];
    size_t pathTarget;
    int saveSVG = 0;
    bool isRunning = FALSE;
    bool isPlotting = FALSE;
    bool isGridAllocated = FALSE;
    bool cancellationCalled = FALSE;
    bool loadedDefaults = FALSE;

    mainGui() : queueUpdate(0),
    queueSliderUpdate(0), 
    queueSliderMove(0),
    sliderTarget(0),
    pathTarget(0), 
    saveSVG(0),
    isRunning(0), 
    isPlotting(0), 
    isGridAllocated(0), 
    cancellationCalled(0),
    loadedDefaults(0){}
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

        progressBarBox.queueDraw();
    }
    void requestSliderUpdate() {
        queueSliderUpdate = TRUE;
    }
    void requestSliderMove(int target) {
        queueSliderMove = TRUE;
        sliderTarget = target;
    }
    void updateSlider() {
        if (queueSliderUpdate) {
            plotSlider.setRange(0, (double)(((*activeSetPtr).Nsims * maxN(1,(*activeSetPtr).Nsims2) - 1)));
            queueSliderUpdate = FALSE;
        }
        if (queueSliderMove) {
            plotSlider.setValue(sliderTarget);
            queueSliderMove = FALSE;
        }
    }

    void activate(GtkApplication* app) {
        activeSetPtr = (simulationParameterSet*)calloc(MAX_SIMULATIONS, sizeof(simulationParameterSet));
        int buttonWidth = 4;
        int textWidth = 3;
        int labelWidth = 6;
        int plotWidth = 12;
        int plotHeight = 6;
        int pathChars = 40;
        int colWidth = labelWidth + 2 * textWidth;
        int textCol1a = labelWidth;
        int textCol2a = textCol1a + 2 * textWidth + labelWidth;
        int textCol1b = textCol1a + textWidth;
        int textCol2b = textCol2a + textWidth;
        int buttonCol1 = textCol2a - labelWidth;
        int buttonCol2 = buttonCol1 + buttonWidth;
        int buttonCol3 = buttonCol2 + buttonWidth;
        g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", TRUE, NULL);
        window.init(app, _T("Lightwave Explorer"), 1920, 1080);
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
        filePaths[0].setMaxCharacters(pathChars);
        pulldowns[0].addElement(_T("Synthetic"));
        pulldowns[0].addElement(_T("FROG"));
        pulldowns[0].addElement(_T("EOS"));
        pulldowns[0].init(parentHandle, labelWidth, 16, 2 * textWidth, 1);
        filePaths[0].setLabel(0, -1, _T("Data 1:"));

        filePaths[1].init(parentHandle, 0, 19, colWidth, 1);
        filePaths[1].setMaxCharacters(pathChars);
        filePaths[1].setLabel(0, -1, _T("Data 2:"));
        pulldowns[1].addElement(_T("Synthetic"));
        pulldowns[1].addElement(_T("FROG"));
        pulldowns[1].addElement(_T("EOS"));
        pulldowns[1].init(parentHandle, labelWidth, 18, 2 * textWidth, 1);

        filePaths[2].init(parentHandle, 0, 21, colWidth, 1);
        filePaths[2].setMaxCharacters(pathChars);
        filePaths[2].setLabel(0, -1, _T("Fit data:"));
        pulldowns[2].addElement(_T("Maximize x"));
        pulldowns[2].addElement(_T("Maximize y"));
        pulldowns[2].addElement(_T("Maximize Total"));
        pulldowns[2].addElement(_T("Fit spectrum"));
        pulldowns[2].init(parentHandle, labelWidth, 20, 2 * textWidth, 1);

        filePaths[3].init(parentHandle, buttonCol1, 16, colWidth, 1);
        //filePaths[3].setLabel(0, -1, _T("Output:"));
        filePaths[3].setMaxCharacters(pathChars);

        drawBoxes[0].init(window.parentHandle(2), 0, 0, plotWidth, plotHeight);
        drawBoxes[0].setDrawingFunction(drawTimeImage1);
        drawBoxes[0].setTooltip("Image of the electric field grid: presents a slice of Ey(x,y=0,t), there the horizontal axis is time, and the vertical axis is position");
        drawBoxes[1].init(window.parentHandle(2), 0, plotHeight, plotWidth, plotHeight);
        drawBoxes[1].setDrawingFunction(drawTimeImage2);
        drawBoxes[1].setTooltip("Image of the electric field grid: presents a slice of Ex(x,y=0,t), there the horizontal axis is time, and the vertical axis is position");
        drawBoxes[2].init(window.parentHandle(2), 0, 2 * plotHeight, plotWidth, plotHeight);
        drawBoxes[2].setDrawingFunction(drawField1Plot);
        drawBoxes[2].setTooltip("Plot of the on-axis electric field in the y polarization");
        drawBoxes[3].init(window.parentHandle(2), 0, 3 * plotHeight, plotWidth, plotHeight);
        drawBoxes[3].setDrawingFunction(drawField2Plot);
        drawBoxes[3].setTooltip("Plot of the on-axis electric field in the x polarization");
        drawBoxes[4].init(window.parentHandle(2), plotWidth, 0, plotWidth, plotHeight);
        drawBoxes[4].setDrawingFunction(drawFourierImage1);
        drawBoxes[4].setTooltip("Plot of the electric field grid in momentum-frequency space: Ey(kx,ky=0,f). Is plotted on a logarithmic scale. Vertical axis is transverse momentum kx, and horizontal axis is frequency f.");
        drawBoxes[5].init(window.parentHandle(2), plotWidth, plotHeight, plotWidth, plotHeight);
        drawBoxes[5].setDrawingFunction(drawFourierImage2);
        drawBoxes[5].setTooltip("Plot of the electric field grid in momentum-frequency space: Ex(kx,ky=0,f). Is plotted on a logarithmic scale. Vertical axis is transverse momentum kx, and horizontal axis is frequency f.");
        drawBoxes[6].init(window.parentHandle(2), plotWidth, 2 * plotHeight, plotWidth, plotHeight);
        drawBoxes[6].setDrawingFunction(drawSpectrum1Plot);
        drawBoxes[6].setTooltip("Plot of the energy spectrum of the result, y-polarization.");
        drawBoxes[7].init(window.parentHandle(2), plotWidth, 3 * plotHeight, plotWidth, plotHeight);
        drawBoxes[7].setDrawingFunction(drawSpectrum2Plot);
        drawBoxes[7].setTooltip("Plot of the energy spectrum of the result, x-polarization.");
        progressBarBox.init(window.parentHandle(5), 0, 0, 1, 1);
        progressBarBox.noVerticalExpantion();
        progressBarBox.setDrawingFunction(drawProgress);
        drawBoxes[7].setDrawingFunction(drawSpectrum2Plot);
        textBoxes[48].init(window.parentHandle(4), 2, 0, 2, 1);
        textBoxes[49].init(window.parentHandle(4), 4, 0, 2, 1);
        textBoxes[50].init(window.parentHandle(4), 8, 0, 2, 1);
        textBoxes[51].init(window.parentHandle(4), 10, 0, 2, 1);

        checkBoxes[0].init(_T("Total"), window.parentHandle(4), 12, 0, 1, 1);
        checkBoxes[0].setTooltip("Overlay a plot of the integrated energy spectrum over the two polarization-resolved spectra");
        checkBoxes[1].init(_T("Log"), window.parentHandle(4), 13, 0, 1, 1);
        checkBoxes[1].setTooltip("Plot spectra on a log10 scale");


        pulldowns[4].addElement(_T("2D Cartesian"));
        pulldowns[4].addElement(_T("3D radial symm."));
        pulldowns[4].addElement(_T("3D"));
        pulldowns[4].init(parentHandle, textCol2a, 7, 2 * textWidth, 1);

        char batchModeNames[38][64] = {
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
        "36: Manual",
        "37: i37"
        };

        for (int i = 0; i < 38; ++i) {
            pulldowns[5].addElement(batchModeNames[i]);
            pulldowns[6].addElement(batchModeNames[i]);
        }
        pulldowns[5].init(parentHandle, textCol2a, 8, 2 * textWidth, 1);
        pulldowns[6].init(parentHandle, textCol2a, 9, 2 * textWidth, 1);
        pulldowns[5].setTooltip("Primary batch mode selector: the selected value from the interface will be scanned in a series of simulations, starting from the value entered on the interface, and ending with the batch target set below. The number of simulations is set by the batch steps parameter.");
        pulldowns[6].setTooltip("Secondary batch mode selector - allows a 2D parameter scan. Works in the same way as the primary batch, but uses the values in the right-hand column.");

        int mbRow = 22;
        textBoxes[31].setLabel(-9 ,7,"Sequence:");
#ifdef __APPLE__
        miniButtons[0].init(_T("="), parentHandle, textWidth + 1, mbRow, 2, 1, buttonAddSameCrystal);
        miniButtons[1].init(_T("d"), parentHandle, textWidth + 3, mbRow, 2, 1, buttonAddDefault);
        miniButtons[2].init(_T("r"), parentHandle, textWidth + 5, mbRow, 2, 1, buttonAddRotation);
        miniButtons[3].init(_T("p"), parentHandle, textWidth + 7, mbRow, 2, 1, buttonAddPulse);
        miniButtons[4].init("(", parentHandle, textWidth + 9, mbRow, 2, 1, buttonAddMirror);
        miniButtons[5].init("f", parentHandle, textWidth + 11, mbRow, 2, 1, buttonAddFilter);
        miniButtons[6].init("l", parentHandle, textWidth + 13, mbRow, 2, 1, buttonAddLinear);
        miniButtons[7].init("a", parentHandle, textWidth + 15, mbRow, 2, 1, buttonAddAperture);
        miniButtons[8].init("ff", parentHandle, textWidth + 17, mbRow, 2, 1, buttonAddFarFieldAperture);
        miniButtons[9].init("f", parentHandle, textWidth + 19, mbRow, 2, 1, buttonAddForLoop);
#else
        miniButtons[0].init(_T("\xf0\x9f\x93\xb8"), parentHandle, textWidth + 1, mbRow, 2, 1, buttonAddSameCrystal);
        miniButtons[1].init(_T("\xe2\x99\x8a"), parentHandle, textWidth + 3, mbRow, 2, 1, buttonAddDefault);
        miniButtons[2].init(_T("\xf0\x9f\x92\xab"), parentHandle, textWidth + 5, mbRow, 2, 1, buttonAddRotation);
        miniButtons[3].init(_T("\xf0\x9f\x92\xa1"), parentHandle, textWidth + 7, mbRow, 2, 1, buttonAddPulse);
        miniButtons[4].init("\xf0\x9f\x94\x8e", parentHandle, textWidth + 9, mbRow, 2, 1, buttonAddMirror);
        miniButtons[5].init("\xf0\x9f\x98\x8e", parentHandle, textWidth + 11, mbRow, 2, 1, buttonAddFilter);
        miniButtons[6].init("\xf0\x9f\x93\x8f", parentHandle, textWidth + 13, mbRow, 2, 1, buttonAddLinear);
        miniButtons[7].init("\xf0\x9f\x91\x8c", parentHandle, textWidth + 15, mbRow, 2, 1, buttonAddAperture);
        miniButtons[8].init("\xe2\x9b\xb3", parentHandle, textWidth + 17, mbRow, 2, 1, buttonAddFarFieldAperture);
        miniButtons[9].init("\xf0\x9f\x94\x81", parentHandle, textWidth + 19, mbRow, 2, 1, buttonAddForLoop);

#endif
        miniButtons[0].setTooltip("Make a copy of the crystal currently entered in the interface");
        miniButtons[1].setTooltip("Insert a crystal that will change with the values set on the interface, or modified during a batch calculation");
        miniButtons[2].setTooltip("Rotate the polarization by a specified angle in degrees");
        miniButtons[3].setTooltip("Add a new pulse to the grid; values will be set to duplicate pulse 1 as entered above");
        miniButtons[4].setTooltip("Add a spherical mirror to the beam path, with radius of curvature in meters");
        miniButtons[5].setTooltip("Add a spectral filter to the beam path. Parameters:\n   central frequency (THz)\n   bandwidth (THz)\n   supergaussian order\n   in-band amplitude\n   out-of-band amplitude\n");
        miniButtons[6].setTooltip("Add a linear propagation through the crystal entered on the interface");
        miniButtons[7].setTooltip("Add an aperture to the beam. Parameters:\n   diameter (m)\n   activation parameter\n");
        miniButtons[8].setTooltip("Filter the beam with a far-field aperture. Parameters:\n   opening angle (deg)\n   activation parameter (k)\n   x-angle (deg)\n   y-angle (deg) ");
        miniButtons[9].setTooltip("Add an empty for loop. Parameters:\n   Number of times to execute\n   Variable number in which to put the counter");
        buttons[0].init(_T("Run"), parentHandle, buttonCol3, 15, buttonWidth, 1, launchRunThread);
        buttons[0].setTooltip("Run the simulation as currently entered on the interface. If a sequence is entered in the sequence box below, that will execute, otherwise, a simulation on the input parameters above and to the left in a single medium will be performed.");
        buttons[1].init(_T("Stop"), parentHandle, buttonCol2, 15, buttonWidth, 1, stopButtonCallback);
        buttons[1].setTooltip("Tell a currently-running simulation to stop. It might not stop right away; it will only happen once it reaches a break point");
        buttons[2].init(_T("Script"), parentHandle, 2 * buttonWidth + 1+buttonCol1, 17, textWidth, 1, createRunFile);
        buttons[2].setTooltip("Generate an input file and SLURM script for running the simulation as entered on the selected cluster");
        buttons[3].init(_T("Fit"), parentHandle, buttonCol3, 12, buttonWidth, 1, launchFitThread);
        buttons[3].setTooltip("Run the fitting routine with the above parameters. The mode is set in the pulldown next to the (optional) fitting input data file path.");
        buttons[4].init(_T("Load"), parentHandle, buttonCol1, 15, buttonWidth, 1, loadCallback);
        buttons[4].setTooltip("Load the results of a previous simulation run. You should select the associated .txt file. The parameters will be loaded into the interface, and the data (if it exists) will be plotted.");
        buttons[6].init(_T("Path"), parentHandle, textWidth, 16, textWidth, 1, openFileDialogCallback, 0);
        buttons[7].init(_T("Path"), parentHandle, textWidth, 18, textWidth, 1, openFileDialogCallback, (gpointer)1);
        buttons[8].init(_T("Path"), parentHandle, textWidth, 20, textWidth, 1, openFileDialogCallback, (gpointer)2);
        buttons[9].init(_T("Path"), parentHandle, buttonCol1, 17, textWidth, 1, saveFileDialogCallback, (gpointer)3);
        buttons[9].setTooltip("Sets the base path for the output files");
        buttons[10].init(_T("xlim"), window.parentHandle(4), 0, 0, 1, 1, independentPlotQueue);
        buttons[10].setTooltip("Apply the entered x limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[10].squeeze();
        buttons[11].init(_T("ylim"), window.parentHandle(4), 6, 0, 1, 1, independentPlotQueue);
        buttons[11].setTooltip("Apply the entered y limits to the plot. The two text boxes are for the upper and lower limits applied to the frequency axis. If they are empty, the range will include the whole grid.");
        buttons[11].squeeze();
        checkBoxes[0].setFunction(independentPlotQueue);
        checkBoxes[1].setFunction(independentPlotQueue);
        buttons[12].init(_T("SVG"), window.parentHandle(3), 5, 0, 1, 1, svgCallback);
        buttons[12].setTooltip("Generate SVG files of the four line plots, with filenames based on the base path set above");
        buttons[12].squeeze();
        plotSlider.init(window.parentHandle(3), 0, 0, 4, 1);
        plotSlider.setRange(0.0, 10.0);
        plotSlider.setDigits(0);
        plotSlider.setFunction(independentPlotQueue);

        sequence.init(window.parentHandle(1), 0, 0, 1, 1);
        fitCommand.init(parentHandle, buttonCol1, 13, colWidth, 2);
        console.init(parentHandle, buttonCol1, 18, colWidth, 4);

        checkLibraryAvailability();
        std::string A;
        if (CUDAavailable) {
            pulldowns[7].addElement("CUDA");
            pulldowns[8].addElement("CUDA");
            for (int i = 1; i < cudaGPUCount; ++i) {
                A = Sformat("CUDA {}", i);
                pulldowns[7].addElement(A.c_str());
                pulldowns[8].addElement(A.c_str());
            }
        }
        if (SYCLavailable) {
            A.assign("SYCL");
            pulldowns[7].addElement(A.c_str());
            pulldowns[8].addElement(A.c_str());
            if (syclGPUCount > 0) {

                pulldowns[7].addElement("SYCL cpu");
                pulldowns[8].addElement("SYCL cpu");
                pulldowns[7].addElement("SYCL gpu");
                pulldowns[8].addElement("SYCL gpu");
            }
        }
        pulldowns[7].addElement("OpenMP");
        pulldowns[8].addElement("OpenMP");

        pulldowns[7].init(window.parentHandle(6), 2 + buttonWidth, 0, buttonWidth, 1);
        pulldowns[8].init(window.parentHandle(6), 4 + 2 * buttonWidth, 0, buttonWidth, 1);
        textBoxes[52].init(window.parentHandle(6), 4 + 3 * buttonWidth, 0, 1, 1);
        pulldowns[7].setTooltip("Select the primary method of calculation. The algorithm is the same, but you can run it either on a GPU or CPU depending on your machine");
        pulldowns[8].setTooltip("Select a secondary mode of calculation for offloading jobs from a batch. For example, if the pulldown to the left is set to CUDA and this one is OpenMP, and the number to the right is 2, 2 of the simulations from the batch will be performed on the CPU");


        pulldowns[7].setLabel(-2, 0, _T("Config:"), 8, 2);
        textBoxes[0].setLabel(-labelWidth, 0, _T("Pulse energy (J)"));
        textBoxes[1].setLabel(-labelWidth, 0, _T("Frequency (THz)"));
        textBoxes[2].setLabel(-labelWidth, 0, _T("Bandwidth (THz)"));
        textBoxes[3].setLabel(-labelWidth, 0, _T("SG order"));
        textBoxes[4].setLabel(-labelWidth, 0, _T("CEP/\xcf\x80"));
        textBoxes[5].setLabel(-labelWidth, 0, _T("Delay (fs)"));
        textBoxes[6].setLabel(-labelWidth, 0, _T("GDD (fs\xc2\xb2)"));
        textBoxes[7].setLabel(-labelWidth, 0, _T("TOD (fs\xc2\xb3)"));
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

        fitCommand.setLabel(0, -1, _T("Fitting:"));
        
        //sequence.setLabel(0, -1, _T("Sequence:"), 11, 3);
        filePaths[3].overwritePrint("TestFile");

        GtkCssProvider* textProvider = gtk_css_provider_new();
        gtk_css_provider_load_from_data(textProvider,
            "label, scale { font-family: Arial; font-weight: bold; }\n button, entry, textview { font-family: Arial; font-weight: bold; color: #EEEEEE; background-color: #191919; }", -1);
        gtk_style_context_add_provider_for_display(gdk_display_get_default(), GTK_STYLE_PROVIDER(textProvider), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
        
        GtkCssProvider* buttonShrinker = gtk_css_provider_new();
        gtk_css_provider_load_from_data(buttonShrinker,
            "label, scale, range, button, entry, textview { min-height: 10px; min-width: 10px; }", -1);
        gtk_style_context_add_provider_for_display(gdk_display_get_default(), GTK_STYLE_PROVIDER(buttonShrinker), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

        //read the crystal database
        crystalDatabasePtr = new crystalEntry[MAX_LOADSTRING]();
        std::string materialString;
        if (crystalDatabasePtr != NULL) {
            readCrystalDatabase(crystalDatabasePtr);
            console.cPrint("Material database has {} entries:\n", (*crystalDatabasePtr).numberOfEntries);
            for (int i = 0; i < (*crystalDatabasePtr).numberOfEntries; ++i) {
                console.cPrint("{} : {} \n", i, std::string(crystalDatabasePtr[i].crystalNameW));
                materialString = Sformat("{:2}: {}", i, std::string(crystalDatabasePtr[i].crystalNameW));
                pulldowns[3].addElement(materialString.c_str());
            }
        }
        pulldowns[3].init(parentHandle, textCol2a, 0, 2 * textWidth, 1);
        pulldowns[3].setLabel(-labelWidth, 0, _T("Material"));

        pulldowns[9].addElement(_T("Cobra 1x R5000"));
        pulldowns[9].addElement(_T("Cobra 2x R5000"));
        pulldowns[9].addElement(_T("Cobra 1x V100"));
        pulldowns[9].addElement(_T("Cobra 2x V100"));
        pulldowns[9].addElement(_T("Raven 1x A100"));
        pulldowns[9].addElement(_T("Raven 2x A100"));
        pulldowns[9].addElement(_T("Raven 4x A100"));
        pulldowns[9].init(parentHandle, 1+buttonCol1, 17, 2 * buttonWidth, 1);
        pulldowns[9].squeeze();
        pulldowns[9].setTooltip("Select the cluster and GPU configuration for generating a SLURM script");
        //pulldowns[9].setLabel(-1, 0, "Cluster:", 8, 3);
        

#ifdef __linux__
		if (1 == readInputParametersFile(activeSetPtr, crystalDatabasePtr, "/usr/share/LightwaveExplorer/DefaultValues.ini")) {
			readInputParametersFile(activeSetPtr, crystalDatabasePtr, "DefaultValues.ini");
		}
#elif defined __APPLE__
		uint32_t bufferSize = 1024;
		char sysPathBuffer[1024] = { 0 };
		_NSGetExecutablePath(sysPathBuffer, &bufferSize);
        std::string sysPathFull(sysPathBuffer);
        std::string sysPathIni = sysPathFull.substr(0,sysPathFull.find_last_of("/"));
        sysPathIni.append("/../Resources/DefaultValues.ini");
		if(1 == readInputParametersFile(activeSetPtr, crystalDatabasePtr, sysPathIni.c_str())){
            readInputParametersFile(activeSetPtr, crystalDatabasePtr, "DefaultValues.ini");
        }
#else
		readInputParametersFile(activeSetPtr, crystalDatabasePtr, "DefaultValues.ini");
#endif
        setInterfaceValuesToActiveValues();
        g_timeout_add(100, G_SOURCE_FUNC(updateDisplay), NULL);
        window.present();
    }
};
mainGui theGui;

char programDirectory[MAX_LOADSTRING];     // Program working directory (useful if the crystal database has to be reloaded)

///////////////////.
//Definitions over
///////////////////.

void updateDisplay() {
    theGui.console.updateFromBuffer();
    theGui.updateSlider();
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
        theGui.textBoxes[i++].setToDouble((*t).cep/PI);
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
        theGui.textBoxes[i++].overwritePrint("{};{}", (int)(1e6 * (*activeSetPtr).spatialWidth), (int)(1e6 * (*activeSetPtr).spatialHeight));
    }
    theGui.textBoxes[i++].setToDouble(1e6*(*activeSetPtr).rStep);
    theGui.textBoxes[i++].setToDouble(1e15 * (*activeSetPtr).timeSpan);
    theGui.textBoxes[i++].setToDouble(1e15 * (*activeSetPtr).tStep);
    theGui.textBoxes[i++].setToDouble(1e6 * (*activeSetPtr).crystalThickness);
    theGui.textBoxes[i++].setToDouble(1e9 * (*activeSetPtr).propagationStep);
    theGui.textBoxes[i++].setToDouble((*activeSetPtr).batchDestination);
    theGui.textBoxes[i++].setToDouble((*activeSetPtr).batchDestination2);
    theGui.textBoxes[i++].setToDouble((double)(*activeSetPtr).Nsims);
    theGui.textBoxes[i++].setToDouble((double)(*activeSetPtr).Nsims2);
    theGui.pulldowns[4].setValue((*activeSetPtr).symmetryType);
    theGui.pulldowns[5].setValue((*activeSetPtr).batchIndex);
    theGui.pulldowns[6].setValue((*activeSetPtr).batchIndex2);

    if (std::string((*activeSetPtr).sequenceString).length() > 6) {
        formatSequence((*activeSetPtr).sequenceString, 2*MAX_LOADSTRING);
        theGui.sequence.directOverwritePrint((*activeSetPtr).sequenceString);
        removeCharacterFromString((*activeSetPtr).sequenceString, 2*MAX_LOADSTRING, '\r');
        removeCharacterFromString((*activeSetPtr).sequenceString, 2*MAX_LOADSTRING, '\n');
    }
    stripLineBreaks((*activeSetPtr).field1FilePath);
    if (std::string((*activeSetPtr).field1FilePath).compare("None") != 0) theGui.filePaths[0].overwritePrint((*activeSetPtr).field1FilePath);
    if (std::string((*activeSetPtr).field2FilePath).compare("None") != 0) theGui.filePaths[1].overwritePrint((*activeSetPtr).field2FilePath);
    if (std::string((*activeSetPtr).fittingPath).compare("None") != 0) theGui.filePaths[2].overwritePrint((*activeSetPtr).fittingPath);

    if (!((*activeSetPtr).fittingString[0] == 'N')) {
        insertLineBreaksAfterSemicolons((*activeSetPtr).fittingString, MAX_LOADSTRING);
        theGui.fitCommand.overwritePrint((*activeSetPtr).fittingString);
        removeCharacterFromString((*activeSetPtr).fittingString, MAX_LOADSTRING, '\r');
        removeCharacterFromString((*activeSetPtr).fittingString, MAX_LOADSTRING, '\n');
    }

    if(!theGui.loadedDefaults) theGui.filePaths[3].overwritePrint((*activeSetPtr).outputBasePath);
    theGui.loadedDefaults = TRUE;
}


void readParametersFromInterface() {
    int i = 0;
    pulse* t = &(*activeSetPtr).pulse1;
    for (int k = 0; k < 2; ++k) {
        theGui.textBoxes[i++].valueToPointer(&(*t).energy);
        theGui.textBoxes[i++].valueToPointer(1e12, &(*t).frequency);
        theGui.textBoxes[i++].valueToPointer(1e12, &(*t).bandwidth);
        theGui.textBoxes[i++].valueToPointer(&(*t).sgOrder);
        theGui.textBoxes[i++].valueToPointer(PI, &(*t).cep);
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
    theGui.textBoxes[i++].valueToPointer(1e-9, &(*activeSetPtr).propagationStep);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).batchDestination);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).batchDestination2);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).Nsims);
    theGui.textBoxes[i++].valueToPointer(&(*activeSetPtr).Nsims2);

    (*activeSetPtr).symmetryType = theGui.pulldowns[4].getValue();
    (*activeSetPtr).batchIndex = theGui.pulldowns[5].getValue();
    (*activeSetPtr).batchIndex2 = theGui.pulldowns[6].getValue();
    (*activeSetPtr).runType = theGui.pulldowns[9].getValue();
    theGui.textBoxes[52].valueToPointer(&(*activeSetPtr).NsimsCPU);

    //char noneString[] = "None";
    std::string noneString("None\0");
    std::string s;
    memset((*activeSetPtr).sequenceString, 0, 2*MAX_LOADSTRING);
    theGui.sequence.copyBuffer((*activeSetPtr).sequenceString, 2*MAX_LOADSTRING);
    s.assign((*activeSetPtr).sequenceString);
    if (s.length() == 0) {
        noneString.copy((*activeSetPtr).sequenceString, 2*MAX_LOADSTRING);
    }
    else {
        if((*activeSetPtr).sequenceString[0] != '0')stripWhiteSpace((*activeSetPtr).sequenceString);
    }
    
    memset((*activeSetPtr).fittingString, 0, MAX_LOADSTRING);
    theGui.fitCommand.copyBuffer((*activeSetPtr).fittingString, MAX_LOADSTRING);
    s.assign((*activeSetPtr).fittingString);
    if (s.length() == 0) {
        noneString.copy((*activeSetPtr).fittingString, MAX_LOADSTRING);
    }
    else {
        stripLineBreaks((*activeSetPtr).fittingString);
    }
    
    memset((*activeSetPtr).field1FilePath, 0, MAX_LOADSTRING);
    theGui.filePaths[0].copyBuffer((*activeSetPtr).field1FilePath, MAX_LOADSTRING);
    s.assign((*activeSetPtr).field1FilePath);
    if (s.length() == 0) {
        noneString.copy((*activeSetPtr).field1FilePath, MAX_LOADSTRING);
    }
    stripLineBreaks((*activeSetPtr).field1FilePath);
    memset((*activeSetPtr).field2FilePath, 0, MAX_LOADSTRING);
    theGui.filePaths[1].copyBuffer((*activeSetPtr).field2FilePath, MAX_LOADSTRING);
    s.assign((*activeSetPtr).field2FilePath);
    if (s.length() == 0) {
        noneString.copy((*activeSetPtr).field2FilePath, MAX_LOADSTRING);
    }
    stripLineBreaks((*activeSetPtr).field2FilePath);
    memset((*activeSetPtr).outputBasePath, 0, MAX_LOADSTRING);
    theGui.filePaths[3].copyBuffer((*activeSetPtr).outputBasePath, MAX_LOADSTRING);
    s.assign((*activeSetPtr).outputBasePath);
    if (s.length() == 0) {
        noneString.copy((*activeSetPtr).outputBasePath, MAX_LOADSTRING);
    }
    stripLineBreaks((*activeSetPtr).outputBasePath);
    
    memset((*activeSetPtr).fittingPath, 0, MAX_LOADSTRING);
    theGui.filePaths[2].copyBuffer((*activeSetPtr).fittingPath, MAX_LOADSTRING);
    s.assign((*activeSetPtr).fittingPath);
    if (s.length() == 0) {
        noneString.copy((*activeSetPtr).fittingPath, MAX_LOADSTRING);
    }
    stripLineBreaks((*activeSetPtr).fittingPath);

    //derived parameters and cleanup:
    (*activeSetPtr).sellmeierType = 0;
    (*activeSetPtr).axesNumber = 0;
    (*activeSetPtr).Ntime = (size_t)(MIN_GRIDDIM * round((*activeSetPtr).timeSpan / (MIN_GRIDDIM * (*activeSetPtr).tStep)));
    if ((*activeSetPtr).symmetryType == 2) {
        (*activeSetPtr).is3D = TRUE;
        (*activeSetPtr).spatialWidth = (*activeSetPtr).rStep * (MIN_GRIDDIM * round((*activeSetPtr).spatialWidth / ((*activeSetPtr).rStep * MIN_GRIDDIM)));
        (*activeSetPtr).Nspace = (size_t)round((*activeSetPtr).spatialWidth / (*activeSetPtr).rStep);
        if ((*activeSetPtr).spatialHeight > 0) {
            (*activeSetPtr).spatialHeight = (*activeSetPtr).rStep * (MIN_GRIDDIM * round((*activeSetPtr).spatialHeight / ((*activeSetPtr).rStep * MIN_GRIDDIM)));
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
        (*activeSetPtr).spatialWidth = (*activeSetPtr).rStep * (MIN_GRIDDIM * round((*activeSetPtr).spatialWidth / ((*activeSetPtr).rStep * MIN_GRIDDIM)));
        (*activeSetPtr).Nspace = (size_t)round((*activeSetPtr).spatialWidth / (*activeSetPtr).rStep);
    }

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
    (*activeSetPtr).NsimsCPU = minN((*activeSetPtr).NsimsCPU, (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2);

    (*activeSetPtr).field1IsAllocated = FALSE;
    (*activeSetPtr).field2IsAllocated = FALSE;

    //crystal from database (database must be loaded!)
    (*activeSetPtr).crystalDatabase = crystalDatabasePtr;
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
#if defined CPUONLY
    CUDAavailable = FALSE;
    SYCLavailable = FALSE;
#define solveNonlinearWaveEquationSequence solveNonlinearWaveEquationSequenceCPU
#define solveNonlinearWaveEquation solveNonlinearWaveEquationCPU
#define solveNonlinearWaveEquationSequenceSYCL solveNonlinearWaveEquationSequenceCPU
#define solveNonlinearWaveEquationSYCL solveNonlinearWaveEquationCPU
#define runDlibFitting runDlibFittingCPU
#define runDlibFittingSYCL runDlibFittingCPU
#else

#ifndef NOCUDA
	//Find, count, and name the GPUs
	int CUDAdevice, i;
	cudaGetDeviceCount(&cudaGPUCount);
	cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
	struct cudaDeviceProp activeCUDADeviceProp;
	//if (cuErr == cudaSuccess) {

	if (cudaGPUCount > 0) {
		CUDAavailable = TRUE;
	}
	theGui.console.cPrint("CUDA found {} GPU(s): \n", cudaGPUCount);
	for (i = 0; i < cudaGPUCount; ++i) {
		cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
		theGui.console.cPrint("{}\n", activeCUDADeviceProp.name);
		theGui.console.cPrint("    Memory: {} MB\n    Multiprocessors: {}\n",
			(int)ceil(((float)activeCUDADeviceProp.totalGlobalMem) / 1048576), activeCUDADeviceProp.multiProcessorCount);
	}
#else
#define solveNonlinearWaveEquationSequence solveNonlinearWaveEquationSequenceCPU
#define solveNonlinearWaveEquation solveNonlinearWaveEquationCPU
#define runDlibFitting runDlibFittingCPU
#endif


#ifndef NOSYCL
    SYCLavailable = TRUE;
    char syclDeviceList[MAX_LOADSTRING] = { 0 };
    char syclDefault[MAX_LOADSTRING] = { 0 };
    size_t syclDevices = 0;

    syclDevices = readSYCLDevices(syclDeviceList, syclDefault);

    unsigned char* counts = (unsigned char*)&syclDevices;
    syclGPUCount = (int)counts[1];
    if(syclDevices != 0){
        theGui.console.cPrint("{}",syclDeviceList);
    }
    
#endif

#endif
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

    for (unsigned short int j = 0; j < 256; ++j) {
        switch (cm) {
        case 0:

            colorMap[j][0] = (unsigned char)j;
            colorMap[j][1] = (unsigned char)j;
            colorMap[j][2] = (unsigned char)j;
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
#pragma omp parallel for private(currentValue) num_threads(interfaceThreads)
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
    cairo_surface_finish(cSurface);
    cairo_surface_destroy(cSurface);
    free(pixels);
    return 0;
}

void pathFromDialogBox(GtkDialog* dialog, int response) {
    if (response == GTK_RESPONSE_ACCEPT) {
        GtkFileChooser* chooser = GTK_FILE_CHOOSER(dialog);
        GFile* file = gtk_file_chooser_get_file(chooser);
        std::string s(g_file_get_path(file));
        if (s.substr(s.length() - 4, std::string::npos) == std::string(".txt")) {
            theGui.filePaths[theGui.pathTarget].overwritePrint("{}", s.substr(0,s.length()-4));
        }
        else {
            theGui.filePaths[theGui.pathTarget].overwritePrint("{}", s);
        }
        
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

void changeToBaseNamePath(char* str, size_t maxSize) {
    std::string s(str, maxSize);
    size_t extension = s.find_last_of(".");
    if (extension != std::string::npos) {
        std::string baseName = s.substr(0, extension);
        baseName.append("\0");
        baseName.copy(str, maxSize);
    }
}

void loadFromDialogBox(GtkDialog* dialog, int response) {
    if (response == GTK_RESPONSE_ACCEPT) {
        GtkFileChooser* chooser = GTK_FILE_CHOOSER(dialog);
        GFile* file = gtk_file_chooser_get_file(chooser);
        std::string path(g_file_get_path(file));
        if (isGridAllocated) {
            freeSemipermanentGrids();
            isGridAllocated = FALSE;
        }
        
        int readParameters = readInputParametersFile(activeSetPtr, crystalDatabasePtr, path.c_str());
        allocateGrids(activeSetPtr);
        isGridAllocated = TRUE;
        if (readParameters == 61) {
            size_t extensionLoc = path.find_last_of(".");
            const std::string basePath = path.substr(0, extensionLoc);
            loadSavedFields(activeSetPtr, basePath.c_str());
            setInterfaceValuesToActiveValues();
            theGui.requestSliderUpdate();
            theGui.requestPlotUpdate();

        }
    }
    g_object_unref(dialog);
}

void loadCallback(GtkWidget* widget, gpointer pathTarget) {
    theGui.pathTarget = (size_t)pathTarget;
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new("Open File", theGui.window.windowHandle(), GTK_FILE_CHOOSER_ACTION_OPEN, "Ok", "Cancel");
    g_signal_connect(fileC, "response", G_CALLBACK(loadFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
}

void svgCallback() {
    theGui.saveSVG = 4;
    theGui.requestPlotUpdate();
}

void saveFileDialogCallback(GtkWidget* widget, gpointer pathTarget) {
    theGui.pathTarget = (size_t)pathTarget;
    //get around bug in GTK4 by opening dialog box directly in cocoa on mac
#ifdef __APPLE__
    NSString *filePath;
    NSSavePanel *savePanel = [NSSavePanel savePanel];
    if ([savePanel runModal] == NSModalResponseOK) {
        filePath = [savePanel URL].path;
        theGui.filePaths[theGui.pathTarget].overwritePrint("{}", [filePath UTF8String]);
    }
    return;
#else
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new("Save File", theGui.window.windowHandle(), GTK_FILE_CHOOSER_ACTION_SAVE, "Ok", "Cancel");
    g_signal_connect(fileC, "response", G_CALLBACK(pathFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
#endif

}

void createRunFile() {

    readParametersFromInterface();

    if (isGridAllocated) {
        freeSemipermanentGrids();
    }

    char* fileName = (*activeSetPtr).outputBasePath;
    while (strchr(fileName, '\\') != NULL) {
        fileName = strchr(fileName, '\\');
        fileName++;
    }
    while (strchr(fileName, '/') != NULL) {
        fileName = strchr(fileName, '/');
        fileName++;
    }

    //create SLURM script
    int gpuType = 0;
    int gpuCount = 1;
    switch ((*activeSetPtr).runType) {
    case 0:
        gpuType = 0;
        gpuCount = 1;
        break;
    case 1:
        gpuType = 0;
        gpuCount = 2;
        break;
    case 2:
        gpuType = 1;
        gpuCount = 1;
        break;
    case 3:
        gpuType = 1;
        gpuCount = 2;
        break;
    case 4:
        gpuType = 2;
        gpuCount = 1;
        break;
    case 5:
        gpuType = 2;
        gpuCount = 2;
        break;
    case 6:
        gpuType = 2;
        gpuCount = 4;
        break;
    }
    saveSlurmScript(activeSetPtr, gpuType, gpuCount);


    //create command line settings file
    (*activeSetPtr).runType = 1;
    saveSettingsFile(activeSetPtr);

    theGui.console.tPrint(
        "Run {} on cluster with:\r\nsbatch {}.slurmScript\r\n",
        fileName, fileName);
    isRunning = FALSE;
    isGridAllocated = FALSE;
}

static void buttonAddSameCrystal() {
    if (theGui.textBoxes[34].valueDouble() != 0.0) {
        theGui.sequence.cPrint("plasma({},{},{},{},{},{},{},{},{})\n",
            theGui.pulldowns[3].getValue(), theGui.textBoxes[32].valueDouble(),
            theGui.textBoxes[33].valueDouble(), theGui.textBoxes[34].valueDouble(),
            theGui.textBoxes[35].valueDouble(), theGui.textBoxes[36].valueDouble(),
            theGui.textBoxes[37].valueDouble(), theGui.textBoxes[42].valueDouble(),
            theGui.textBoxes[43].valueDouble());
    }
    else {
        theGui.sequence.cPrint("nonlinear({},{},{},{},{})\n",
            theGui.pulldowns[3].getValue(), theGui.textBoxes[32].valueDouble(),
            theGui.textBoxes[33].valueDouble(), theGui.textBoxes[42].valueDouble(),
            theGui.textBoxes[43].valueDouble());
    }
}

static void buttonAddDefault() {
    theGui.sequence.cPrint("plasma(d,d,d,d,d,d,d,d,d)\n");
}

static void buttonAddMirror() {
    theGui.sequence.cPrint("sphericalMirror(-1.0)\n");
}

static void buttonAddFilter() {
    theGui.sequence.cPrint("filter(130, 20, 4, 1, 0)\n");
}

static void buttonAddLinear() {
    theGui.sequence.cPrint("linear({},{},{},{},{})\n",
        theGui.pulldowns[3].getValue(), theGui.textBoxes[32].valueDouble(),
        theGui.textBoxes[33].valueDouble(), theGui.textBoxes[42].valueDouble(),
        theGui.textBoxes[43].valueDouble());
}

static void buttonAddAperture() {
    theGui.sequence.cPrint("aperture(0.001, 2)\n");
}

static void buttonAddFarFieldAperture() {
    theGui.sequence.cPrint("farFieldAperture(2.0,4000,0,0)\n");
}

static void buttonAddForLoop() {
    theGui.sequence.cPrint("for(10,1){{\n\n}}\n");
}
static void buttonAddPulse() {
    theGui.sequence.cPrint("addPulse({},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})\n",
        theGui.textBoxes[0].valueDouble(),
        theGui.textBoxes[1].valueDouble(),
        theGui.textBoxes[2].valueDouble(),
        theGui.textBoxes[3].valueDouble(),
        theGui.textBoxes[4].valueDouble(),
        theGui.textBoxes[5].valueDouble(),
        theGui.textBoxes[6].valueDouble(),
        theGui.textBoxes[7].valueDouble(),
        theGui.textBoxes[8].valueDouble(),
        theGui.textBoxes[9].valueDouble(),
        theGui.textBoxes[10].valueDouble(),
        theGui.textBoxes[11].valueDouble(),
        theGui.textBoxes[12].valueDouble(),
        theGui.textBoxes[13].valueDouble(),
        0.0,
        theGui.textBoxes[14].valueDouble(),
        theGui.textBoxes[15].valueDouble(),
        theGui.pulldowns[3].getValue(),
        theGui.textBoxes[32].valueDouble(),
        theGui.textBoxes[33].valueDouble());
}

static void buttonAddRotation() {
    theGui.sequence.cPrint("rotate(90)\n");
}

static void activate(GtkApplication* app, gpointer user_data) {
#if defined __linux__ || defined __APPLE__
    setlocale(LC_NUMERIC, "en_US.UTF-8");
#else
    setlocale(LC_NUMERIC, "en_US");
#endif
    theGui.activate(app);
}

int LwePlot2d(plotStruct* inputStruct) {
    plotStruct* s = (plotStruct*)inputStruct;
    if ((*s).Npts == 0) return 1;
    size_t iMin = 0;
    size_t iMax = (*s).Npts;
    std::string svgString;
    cairo_t* cr = (*s).cr;
    cairo_font_extents_t fe;
    memset(&fe, 0, sizeof(cairo_font_extents_t));
    cairo_text_extents_t te;
    memset(&te, 0, sizeof(cairo_text_extents_t));
    double fontSize = 14.0;
    cairo_set_font_size(cr, fontSize);
    cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_font_extents(cr, &fe);
    double x1, y1, x2, y2;
    double width;
    double height;

    double layoutTop = 0.0;
    double layoutBottom = 0.0;
    double layoutLeft = 0.0;
    double layoutRight = 0.0;

    //Fixed parameters affecting aesthetics
    double radius = 2;
    double lineWidth = 1.5;
    double axisSpaceX = 75.0;
    double axisSpaceY = 35.0;
    double axisLabelSpaceX = 21.0;

    //get limits and make the plotting arrays
    double maxY = -1.0e300;
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
        if ((*s).logScale) minY += log10((*s).unitY);
    }
    if ((*s).forceYmax) {
        maxY = (double)(*s).forcedYmax;
        if ((*s).logScale) maxY += log10((*s).unitY);
    }
    if ((*s).forceXmin) {
        minX = (double)(*s).forcedXmin;
    }
    if ((*s).forceXmax) {
        maxX = (double)(*s).forcedXmax;
    }

    //Tickmark labels
    int NyTicks = 3;
    std::string messageBuffer;
    double yTicks1[3] = { maxY, 0.5 * (maxY + minY), minY };
    double xTicks1[3] = { minX + 0.25 * (maxX - minX), minX + 0.5 * (maxX - minX), minX + 0.75 * (maxX - minX) };
    height = (double)(*s).height;
    width = (double)(*s).width;

    //Start SVG file if building one
    auto SVGh = [&](double x) {
        return (int)(15 * x);
    };
    if ((*s).makeSVG) {
        svgString.append("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
        svgString.append(Sformat("<svg width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n",
            width, height, width, height));
        svgString.append(Sformat("<rect fill=\"#{:x}{:x}{:x}\" stroke=\"#000\" x=\"0\" y=\"0\" width=\"{}\" height=\"{}\"/>\n",
            SVGh(0.0f), SVGh(0.0f), SVGh(0.0f), width, height));
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
        if ((*s).makeSVG)svgString.append(Sformat("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#{:x}{:x}{:x}\" stroke-width=\"{}\"/>\n", x1, y1, x2, y2, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), lineWidth));
    };

    auto SVGstdcircle = [&]() {
        if ((*s).makeSVG)svgString.append(Sformat("<circle cx=\"{}\" cy=\"{}\" r=\"{}\" stroke=\"none\" fill=\"#{:x}{:x}{:x}\" />\n", x1, y1, radius, currentColor.rHex(), currentColor.gHex(), currentColor.bHex()));
    };

    auto SVGstartgroup = [&]() {
        if ((*s).makeSVG)svgString.append("<g>\n");
    };

    auto SVGendgroup = [&]() {
        if ((*s).makeSVG)svgString.append("</g>\n");
    };

    auto SVGcentertext = [&]() {
        if ((*s).makeSVG)svgString.append(Sformat("<text font-family=\"Arial\" font-size=\"{}\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" y=\"{}\" text-anchor=\"middle\">\n{}\n</text>\n", fontSize-1, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 0.5 * (layoutLeft + layoutRight), 0.5 * (layoutBottom + layoutTop - te.height), messageBuffer));
    };

    auto SVGlefttext = [&]() {
        if ((*s).makeSVG)svgString.append(Sformat("<text font-family=\"Arial\" font-size=\"{}\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" y=\"{}\">\n{}\n</text>\n", fontSize-1, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), layoutLeft, layoutTop + fontSize, std::string(messageBuffer)));
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
    auto cairoRightText = [&]() {
        currentColor.setCairo(cr);
        cairo_text_extents(cr, messageBuffer.c_str(), &te);
        cairo_move_to(cr, layoutRight-te.width-3, 0.5 * (layoutBottom + layoutTop - te.height));
        cairo_show_text(cr, messageBuffer.c_str());
    };
    auto cairoCenterText = [&]() {
        currentColor.setCairo(cr);
        cairo_text_extents(cr, messageBuffer.c_str(), &te);
        cairo_move_to(cr, 0.5 * (layoutLeft + layoutRight - te.width), 0.5 * (layoutBottom + layoutTop - te.height));
        cairo_show_text(cr, messageBuffer.c_str());
    };

    auto cairoVerticalText = [&]() {
        currentColor.setCairo(cr);
        cairo_text_extents(cr, messageBuffer.c_str(), &te);
        cairo_move_to(cr, 0.0, height);
        cairo_rotate(cr, -3.1415926535897931 / 2);
        cairo_rel_move_to(cr, 0.5 * (layoutLeft + layoutRight - te.width), fontSize);
        cairo_show_text(cr, messageBuffer.c_str());
        cairo_rotate(cr, 3.1415926535897931 / 2);
    };

    currentColor = (*s).textColor;
    //y-tick text labels
	for (int i = 0; i < NyTicks; ++i) {
		double ytVal = yTicks1[i] / (*s).unitY;
        if ((*s).logScale) ytVal = yTicks1[i] - log10((*s).unitY);
		if (abs(ytVal) > 10.0 || abs(ytVal) < 0.01) {
			messageBuffer = Sformat("{:.1e}", ytVal);
		}
		else {
			messageBuffer = Sformat("{:4.4f}", ytVal);
		}

        layoutLeft = axisLabelSpaceX;
        layoutTop = (double)(i * (0.5 * (height)));
        if (i == 2) layoutTop -= 8.0f;
        if (i == 1) layoutTop -= 6.0f;
        layoutBottom = layoutTop + axisSpaceY;
        layoutRight = axisSpaceX;


        cairoRightText();
        SVGlefttext();
    }
    //y-axis name
    if ((*s).yLabel != NULL) {
        messageBuffer = std::string((*s).yLabel);
        layoutLeft = 0;
        layoutTop = height;
        layoutBottom = height + axisSpaceY;
        layoutRight = height;

        cairoVerticalText();
        if ((*s).makeSVG)svgString.append(Sformat("<text font-family=\"Arial\" font-size=\"{}\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" y=\"{}\" text-anchor=\"middle\" transform=\"translate({}, {}) rotate(-90)\">\n{}\n</text>\n", fontSize, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 0.5 * (layoutLeft + layoutRight), layoutTop + fontSize, -(layoutLeft + layoutRight), height, messageBuffer));
    }

    //x-axis name
    if ((*s).xLabel != NULL) {
        layoutLeft = axisSpaceX;
        layoutTop = height + 2.8 * fontSize;
        layoutBottom = height + axisSpaceY;
        layoutRight = axisSpaceX + width;
        messageBuffer.assign((*s).xLabel);
        cairoCenterText();
        SVGcentertext();
    }

    //x-axis tick labels
    for (int i = 0; i < 3; ++i) {
        messageBuffer.assign(Sformat("{}", (int)round(xTicks1[i])));
        layoutLeft = (double)(axisSpaceX + 0.25 * width * ((size_t)(i)+1) - axisSpaceX / 2);
        layoutTop = height+3;
        layoutBottom = height + axisSpaceY;
        layoutRight = layoutLeft + axisSpaceX;

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
        svgString.append("</svg>");
        (*s).SVG.assign(svgString);
    }
    return 0;
}

void drawProgress(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    int x0 = 5;
    width -= x0;
    LweColor black(0.05, 0.05, 0.05, 0.05);
    cairo_rectangle(cr, x0, 0, width, height);
    black.setCairo(cr);
    cairo_fill(cr);

    size_t lengthEstimate = 0;
    if (!(*activeSetPtr).isInFittingMode) {
        lengthEstimate = totalSteps;
    }
    else {
        lengthEstimate = (*activeSetPtr).fittingMaxIterations;
    }

    double newFraction =  0.0;
    if(lengthEstimate) newFraction = minN(1.0, ((double)progressCounter) / (double)lengthEstimate);
    LweColor magenta(1, 0, 1, 0.0);
    LweColor cyan(0, 1, 1, 0);
    cairo_rectangle(cr, x0, height / 3, round(width*newFraction), height / 2);
    magenta.setCairo(cr);
    cairo_fill(cr);
    cairo_rectangle(cr, x0, 2*height/6, round(width * newFraction), height / 4);
    cyan.setCairo(cr);
    cairo_fill(cr);
    return;
}
void drawField1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    plotStruct sPlot;

    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }

    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());

    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
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
    sPlot.makeSVG = saveSVG; // theGui.saveSVG;

    LwePlot2d(&sPlot);

    if (saveSVG) {
        char* svgFilename = new char[MAX_LOADSTRING]();
        theGui.filePaths[3].copyBuffer(svgFilename, MAX_LOADSTRING);
        std::string svgPath(svgFilename);
        svgPath.append("_Ey.svg");
        std::ofstream fs(svgPath);
        fs.write(sPlot.SVG.c_str(),sPlot.SVG.size());
        delete[] svgFilename;
    }
}

void drawField2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    plotStruct sPlot;

    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }

    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());

    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
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
    sPlot.makeSVG = saveSVG;


    LwePlot2d(&sPlot);

    if (saveSVG) {
        char* svgFilename = new char[MAX_LOADSTRING]();
        theGui.filePaths[3].copyBuffer(svgFilename, MAX_LOADSTRING);
        std::string svgPath(svgFilename);
        svgPath.append("_Ex.svg");
        std::ofstream fs(svgPath);
        fs.write(sPlot.SVG.c_str(), sPlot.SVG.size());
        delete[] svgFilename;
    }
}

void drawSpectrum1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    plotStruct sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = FALSE;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = TRUE;
    }
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
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

    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = (*activeSetPtr).fStep / 1e12;
    sPlot.data = &(*activeSetPtr).totalSpectrum[simIndex * 3 * (*activeSetPtr).Nfreq];
    sPlot.Npts = (*activeSetPtr).Nfreq;
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceY;
    //sPlot.forcedYmin = -4.0;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = LweColor(0.5, 0, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Sy (J/THz)";
    sPlot.unitY = 1.0e-12;
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    if (overlayTotal) {
        sPlot.data2 = &(*activeSetPtr).totalSpectrum[(2 + simIndex * 3) * (*activeSetPtr).Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    sPlot.makeSVG = saveSVG;

    LwePlot2d(&sPlot);

    if (saveSVG) {
        char* svgFilename = new char[MAX_LOADSTRING]();
        theGui.filePaths[3].copyBuffer(svgFilename, MAX_LOADSTRING);
        std::string svgPath(svgFilename);
        svgPath.append("Sy.svg");
        std::ofstream fs(svgPath);
        fs.write(sPlot.SVG.c_str(), sPlot.SVG.size());
        delete[] svgFilename;
    }
}

void drawSpectrum2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    plotStruct sPlot;

    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = FALSE;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = TRUE;
    }
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
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
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = (*activeSetPtr).fStep / 1e12;
    sPlot.data = &(*activeSetPtr).totalSpectrum[(1 + simIndex * 3) * (*activeSetPtr).Nfreq];
    sPlot.Npts = (*activeSetPtr).Nfreq;
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceY;
    //sPlot.forcedYmin = -4.0;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = LweColor(1, 0, 0.5, 0.0);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Sx (J/THz)";
    sPlot.unitY = 1.0e-12;
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    if (overlayTotal) {
        sPlot.data2 = &(*activeSetPtr).totalSpectrum[(2 + simIndex * 3) * (*activeSetPtr).Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    sPlot.makeSVG = saveSVG;


    LwePlot2d(&sPlot);

    if (saveSVG) {
        char* svgFilename = new char[MAX_LOADSTRING]();
        theGui.filePaths[3].copyBuffer(svgFilename, MAX_LOADSTRING);
        std::string svgPath(svgFilename);
        svgPath.append("_Sx.svg");
        std::ofstream fs(svgPath);
        fs.write(sPlot.SVG.c_str(), sPlot.SVG.size());
        delete[] svgFilename;
    }
}


int linearRemapZToLogFloatShift(std::complex<double>* A, int nax, int nay, float* B, int nbx, int nby, double logMin) {
    float f;
    int div2 = nax / 2;
    int nx0, ny0;
#pragma omp parallel for private(nx0, ny0, f) num_threads(interfaceThreads)
    for (int i = 0; i < nbx; ++i) {
        f = i * (nax / (float)nbx);
        nx0 = minN((int)f, nax);
        nx0 -= div2 * ((nx0 >= div2) - (nx0 < div2));
        nx0 *= nay;
        for (int j = 0; j < nby; ++j) {
            f = (j * (nay / (float)nby));
            ny0 = minN(nay, (int)f);
            B[i * nby + j] = (float)log10(cModulusSquared(A[ny0 + nx0]) + logMin);
        }
    }
    return 0;
}

int linearRemapZToLogFloat(std::complex<double>* A, int nax, int nay, float* B, int nbx, int nby, double logMin) {
    float A00;
    float f;
    int nx0, ny0;
    int Ni, Nj;
    for (int i = 0; i < nbx; ++i) {
        f = i * (nax / (float)nbx);
        Ni = (int)f;
        nx0 = nay * minN(Ni, nax);
        for (int j = 0; j < nby; ++j) {
            f = (j * (nay / (float)nby));
            Nj = (int)f;
            ny0 = minN(nay, Nj);
            A00 = (float)log10(cModulusSquared(A[ny0 + nx0]) + logMin);
            B[i * nby + j] = A00;
        }
    }
    return 0;
}

int linearRemapDoubleToFloat(double* A, int nax, int nay, float* B, int nbx, int nby) {
    int nx0, ny0;
    int Ni, Nj;
#pragma omp parallel for private(nx0, ny0, Ni, Nj) num_threads(interfaceThreads)
    for (int i = 0; i < nbx; ++i) {
        Ni = (int)(i * (nax / (float)nbx));
        nx0 = nay * minN(Ni, nax);
        for (int j = 0; j < nby; ++j) {
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
        linearRemapZToLogFloatShift((*s).complexData, (int)(*activeSetPtr).Nspace, (int)(*activeSetPtr).Nfreq, plotarr2, (int)dy, (int)dx, (*s).logMin);
        break;
    }
    drawArrayAsBitmap((*s).cr, (*s).width, (*s).height, plotarr2, (*s).colorMap);
    free(plotarr2);
}


void drawTimeImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    imagePlotStruct sPlot;
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
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
    if (!isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    imagePlotStruct sPlot;
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
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
    if (!isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    imagePlotStruct sPlot;
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }
    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

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
    if (!isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    imagePlotStruct sPlot;

    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2) {
        simIndex = 0;
    }

    double logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).timeSpan));
    if ((*activeSetPtr).is3D) {
        logPlotOffset = (double)(1e-4 / ((*activeSetPtr).spatialWidth * (*activeSetPtr).spatialHeight * (*activeSetPtr).timeSpan));
    }

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
    (*activeSetPtr).NsimsCPU = theGui.textBoxes[52].valueInt();
    if(!isRunning) std::thread(mainSimThread, theGui.pulldowns[7].getValue(), theGui.pulldowns[8].getValue()).detach();
}

void launchFitThread() {
    if (!isRunning) std::thread(fittingThread, theGui.pulldowns[7].getValue()).detach();
}

void stopButtonCallback() {
    if (isRunning) {
        cancellationCalled = TRUE;
        for (int i = 0; i < (*activeSetPtr).Nsims; ++i) {
            (*activeSetPtr).statusFlags[i] = 2;
        }
    }
}

void independentPlotQueue(){
    theGui.requestPlotUpdate();
    theGui.applyUpdate();
}

void secondaryQueue(simulationParameterSet* cpuSims, int pulldownSelection, int pulldownSelectionPrimary) {
    if ((*activeSetPtr).NsimsCPU < 1) return;
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
    //launch on CUDA if selected, putting in the correct GPU in multi-gpu systems
    if (pulldownSelection < cudaGPUCount) {
        sequenceFunction = &solveNonlinearWaveEquationSequence;
        normalFunction = &solveNonlinearWaveEquation;
        assignedGPU = pulldownSelection;
    }
    //launch on SYCL, but if the primary queue matches, default to openMP
    else if (pulldownSelection == cudaGPUCount && SYCLitems > 0) {
        if (pulldownSelection == pulldownSelectionPrimary) {
            theGui.console.tPrint("Sorry, can't run two identical SYCL queues\r\n- defaulting to OpenMP for the secondary queue.\r\n");
            sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
            normalFunction = &solveNonlinearWaveEquationCPU;
        }
        else {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
            normalFunction = &solveNonlinearWaveEquationSYCL;
        }
    }
    else if (pulldownSelection == cudaGPUCount + 1 && SYCLitems > 1) {
        if (pulldownSelection == pulldownSelectionPrimary) {
            theGui.console.tPrint("Sorry, can't run two identical SYCL queues\r\n- defaulting to OpenMP for the secondary queue.\r\n");
            sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
            normalFunction = &solveNonlinearWaveEquationCPU;
        }
        else {
            forceCPU = 1;
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
            normalFunction = &solveNonlinearWaveEquationSYCL;
        }
    }
    else if (pulldownSelection == cudaGPUCount + 2 && SYCLitems > 1) {
        if (pulldownSelection == pulldownSelectionPrimary) {
            theGui.console.tPrint("Sorry, can't run two identical SYCL queues\r\n- defaulting to OpenMP for the secondary queue.\r\n");
            sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
            normalFunction = &solveNonlinearWaveEquationCPU;
        }
        else {
            assignedGPU = 1;
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
            normalFunction = &solveNonlinearWaveEquationSYCL;
        }
    }
    else {
        sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
        normalFunction = &solveNonlinearWaveEquationCPU;
    }

    int error = 0;
    if ((*activeSetPtr).isInSequence) {
        for (unsigned int i = 0; i < (*activeSetPtr).NsimsCPU; ++i) {
            cpuSims[i].assignedGPU = assignedGPU;
            cpuSims[i].runningOnCPU = forceCPU;
            error = sequenceFunction(&cpuSims[i]);
            if (error) break;
        }
    }
    else {
        for (unsigned int i = 0; i < (*activeSetPtr).NsimsCPU; ++i) {
            cpuSims[i].assignedGPU = assignedGPU;
            cpuSims[i].runningOnCPU = forceCPU;
            error = normalFunction(&cpuSims[i]);
            if (error) break;
        }
    }

    if (error) {
        theGui.console.tPrint("Encountered error {} in secondary queue.\r\n", error);
        return;
    }

    return;
}

void mainSimThread(int pulldownSelection, int secondPulldownSelection) {
    cancellationCalled = FALSE;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();


    if (isGridAllocated) {
        freeSemipermanentGrids();
        isGridAllocated = FALSE;
    }
    memset(activeSetPtr, 0, sizeof(simulationParameterSet));
    readParametersFromInterface();
    if ((*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 > MAX_SIMULATIONS) {
        theGui.console.tPrint("Too many simulations in batch mode. Must be under {} total.\r\n", MAX_SIMULATIONS);
    }
    (*activeSetPtr).runType = 0;
    allocateGrids(activeSetPtr);
    isGridAllocated = TRUE;
    theGui.requestSliderUpdate();
    (*activeSetPtr).isFollowerInSequence = FALSE;
    (*activeSetPtr).crystalDatabase = crystalDatabasePtr;
    loadPulseFiles(activeSetPtr);

    if ((*activeSetPtr).sequenceString[0] != 'N') (*activeSetPtr).isInSequence = TRUE;
    
    configureBatchMode(activeSetPtr);

    simulationParameterSet* testSet = new simulationParameterSet[(*activeSetPtr).Nsims * (*activeSetPtr).Nsims2]();
    memcpy(testSet, activeSetPtr, (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 * sizeof(simulationParameterSet));
    totalSteps = 0;
    for (int j = 0; j < (*activeSetPtr).Nsims * (*activeSetPtr).Nsims2; j++) {
        if ((*activeSetPtr).isInSequence) {
            testSet[j].progressCounter = &totalSteps;
            testSet[j].runType = -1;
            solveNonlinearWaveEquationSequenceCounter(&testSet[j]);
        }
        else {
            testSet[j].progressCounter = &totalSteps;
            solveNonlinearWaveEquationCounter(&testSet[j]);
        }
    }
    delete[] testSet;
    int error = 0;
    //run the simulations
    isRunning = TRUE;
    progressCounter = 0;
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

    std::thread secondQueueThread(secondaryQueue, 
        &activeSetPtr[(*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 - (*activeSetPtr).NsimsCPU], 
        secondPulldownSelection, pulldownSelection);



    for (int j = 0; j < ((*activeSetPtr).Nsims * (*activeSetPtr).Nsims2 - (*activeSetPtr).NsimsCPU); ++j) {

        activeSetPtr[j].runningOnCPU = forceCPU;
        activeSetPtr[j].assignedGPU = assignedGPU;
        if ((*activeSetPtr).isInSequence) {
            error = sequenceFunction(&activeSetPtr[j]);
            if (activeSetPtr[j].memoryError != 0) {
                if (activeSetPtr[j].memoryError == -1) {
                    theGui.console.tPrint(_T("Not enough free GPU memory, sorry.\n"), activeSetPtr[j].memoryError);
                }
                else {
                    theGui.console.tPrint(_T("Warning: device memory error ({}).\n"), activeSetPtr[j].memoryError);
                }
            }
            if (error) break;
            theGui.requestSliderMove(j);
            independentPlotQueue();
        }
        else {
            error = normalFunction(&activeSetPtr[j]);
            if (activeSetPtr[j].memoryError != 0) {
                if (activeSetPtr[j].memoryError == -1) {
                    theGui.console.tPrint(_T("Not enough free GPU memory, sorry.\n"), activeSetPtr[j].memoryError);
                }
                else {
                    theGui.console.tPrint(_T("Warning: device memory error ({}).\r\n"), activeSetPtr[j].memoryError);
                }
            }
            if (error) break;
        }

        if (cancellationCalled) {
            theGui.console.tPrint(_T("Warning: series cancelled, stopping after {} simulations.\r\n"), j + 1);
            break;
        }
        theGui.requestSliderMove(j);
        independentPlotQueue();
    }

    if (secondQueueThread.joinable()) secondQueueThread.join();
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    if (error == 13) {
        theGui.console.tPrint(
            "NaN detected in grid!\r\nTry using a larger spatial/temporal step\r\nor smaller propagation step.\r\nSimulation was cancelled.\r\n");
    }
    else {
        theGui.console.tPrint(_T("Finished after {:.4} s. \r\n"), 1e-6 *
            (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
    }

    saveDataSet(activeSetPtr);
    deallocateGrids(activeSetPtr, FALSE);
    isRunning = FALSE;
}

void fittingThread(int pulldownSelection) {
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
    if ((*activeSetPtr).Nfitting == 0) {
        theGui.console.tPrint("Couldn't interpret fitting command.\n");
        free((*activeSetPtr).statusFlags);
        free((*activeSetPtr).deffTensor);
        free((*activeSetPtr).loadedField1);
        free((*activeSetPtr).loadedField2);
        return;
    }
    progressCounter = 0;
    (*activeSetPtr).progressCounter = &progressCounter;
    if ((*activeSetPtr).fittingMode == 3) {
        if (loadReferenceSpectrum((*activeSetPtr).fittingPath, activeSetPtr)) {
            theGui.console.tPrint("Could not read reference file!\n");
            free((*activeSetPtr).statusFlags);
            free((*activeSetPtr).deffTensor);
            free((*activeSetPtr).loadedField1);
            free((*activeSetPtr).loadedField2);
            return;
        }
    }

    theGui.console.tPrint("Fitting {} values in mode {} over {} iterations.\nRegion of interest contains {} elements\n",
        (*activeSetPtr).Nfitting, (*activeSetPtr).fittingMode, (*activeSetPtr).fittingMaxIterations, (*activeSetPtr).fittingROIsize);

    int assignedGPU = 0;
    bool forceCPU = 0;
    int SYCLitems = 0;
    if (syclGPUCount == 0) {
        SYCLitems = (int)SYCLavailable;
    }
    else {
        SYCLitems = 3;
    }
    auto fittingFunction = &runDlibFittingCPU;
    if (pulldownSelection < cudaGPUCount) {
        fittingFunction = &runDlibFitting;
        assignedGPU = pulldownSelection;
    }
    else if (pulldownSelection == cudaGPUCount && SYCLavailable) {
        fittingFunction = &runDlibFittingSYCL;
    }
    else if (pulldownSelection == cudaGPUCount + 1 && SYCLitems > 1) {
        forceCPU = 1;
        fittingFunction = &runDlibFittingSYCL;
    }
    else if (pulldownSelection == cudaGPUCount + 2 && SYCLitems > 1) {
        assignedGPU = 1;
        fittingFunction = &runDlibFittingSYCL;
    }
    fittingFunction(activeSetPtr);
    (*activeSetPtr).plotSim = 0;
    (*activeSetPtr).runningOnCPU = forceCPU;
    (*activeSetPtr).assignedGPU = assignedGPU;
    theGui.requestPlotUpdate();
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    theGui.console.tPrint(_T("Finished fitting after {:.4} s.\n"), 1e-6 *
        (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
    saveDataSet(activeSetPtr);
    setInterfaceValuesToActiveValues();
    theGui.console.tPrint("Fitting result:\n (index, value)\n");
    for (int i = 0; i < (*activeSetPtr).Nfitting; ++i) {
        theGui.console.tPrint("{},  {}\n", i, (*activeSetPtr).fittingResult[i]);
    }
    free((*activeSetPtr).statusFlags);
    free((*activeSetPtr).deffTensor);
    free((*activeSetPtr).loadedField1);
    free((*activeSetPtr).loadedField2);
    isRunning = FALSE;
}

int main(int argc, char **argv){
    GtkApplication* app = gtk_application_new("nickkarpowicz.lighwave", (GApplicationFlags)0);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    return g_application_run(G_APPLICATION(app), argc, argv);
}