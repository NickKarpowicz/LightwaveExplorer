#include "LightwaveExplorerFrontendGTK.h"
#include <thread>
#include <chrono>
#include <locale>
#include "../LightwaveExplorerDevices/LightwaveExplorerCoreCPU.h"
#include "../LightwaveExplorerDevices/LightwaveExplorerCoreCounter.h"
#include "../LightwaveExplorerDevices/LightwaveExplorerCoreFP32.cuh"
#include "../LightwaveExplorerDevices/LightwaveExplorerCoreCPUFP32.h"

//conditional includes and definitions
#ifdef __APPLE__
#include <mach-o/dyld.h>
#import<Cocoa/Cocoa.h>
#endif
#ifndef CPUONLY
#ifndef NOCUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <nvml.h>
#include "../LightwaveExplorerCore.cuh"
#endif
#ifdef _WIN32
#include "../LightwaveExplorerDevices/LightwaveExplorerSYCL.h"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
bool isIntelRuntimeInstalled() {
    wchar_t loadBuffer[1024];
    DWORD envcount = GetEnvironmentVariableW(L"INTEL_DEV_REDIST", loadBuffer, 16);
    if (envcount != 0) {
        return true;
    }
    return false;
}
#else
bool isIntelRuntimeInstalled() {
    return true;
}
#include "LightwaveExplorerDPCPPlib.h"
#include "LightwaveExplorerDPCPPlibFP32.h"
#endif
#endif

//Main data structures:
// theSim contains all of the parameters of the current simulation including grid arrays
// theDatabase is the database of crystal properties
std::vector<simulationParameterSet> theSim(1);
crystalDatabase theDatabase;

//Counter atomics
std::atomic_uint32_t progressCounter;
std::atomic_uint32_t totalSteps;

//Limit the number of threads used to draw the interface if the processor supports a lot
const int interfaceThreads = maxN(std::thread::hardware_concurrency() / 2, 2u);

//Main class for controlling the interface
class mainGui {
    bool queueUpdate;
    bool queueSliderUpdate;
    bool queueSliderMove;
    bool queueInterfaceValuesUpdate;
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
    LweCheckBox checkBoxes[4];
    LweSlider plotSlider;
    LweWindow window;
    LweSpacer spacers[2];
    size_t pathTarget;
    int saveSVG = 0;
    bool loadedDefaults = false;
    unsigned int timeoutID;
    mainGui() : queueUpdate(0),
    queueSliderUpdate(0), 
    queueSliderMove(0),
    queueInterfaceValuesUpdate(0),
    sliderTarget(0),
    pathTarget(0), 
    saveSVG(0),
    loadedDefaults(0),
    timeoutID(0){}
    ~mainGui() {}
    void requestPlotUpdate() {
        queueUpdate = true;
    }
    void applyUpdate() {
        if (queueUpdate) {
            queueUpdate = false;
            for (int i = 0; i < 8; ++i) {
                drawBoxes[i].queueDraw();
            }
        }

        progressBarBox.queueDraw();
        if (queueInterfaceValuesUpdate){
            setInterfaceValuesToActiveValues();
            queueInterfaceValuesUpdate = false;
        }
    }
    void requestSliderUpdate() {
        queueSliderUpdate = true;
    }
    void requestSliderMove(int target) {
        queueSliderMove = true;
        sliderTarget = target;
    }
    void requestInterfaceValuesUpdate(){
        queueInterfaceValuesUpdate = true;
    }


    void updateSlider() {
        if (queueSliderUpdate) {
            plotSlider.setRange(0, (double)((theSim[0].Nsims * maxN(theSim[0].Nsims2, 1u) - 1)));
            queueSliderUpdate = false;
        }
        if (queueSliderMove) {
            plotSlider.setValue(sliderTarget);
            queueSliderMove = false;
        }
    }

    void activate(GtkApplication* app) {
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
        g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", true, NULL);
        window.init(app, _T("Lightwave Explorer"), 1600, 980);
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
        pulldowns[2].addElement(_T("Fit spectrum (log)"));
        pulldowns[2].init(parentHandle, labelWidth, 20, 2 * textWidth, 1);

        filePaths[3].init(parentHandle, buttonCol1, 16, colWidth, 1);
        filePaths[3].setMaxCharacters(pathChars);

        drawBoxes[0].init(window.parentHandle(2), 0, 0, plotWidth, plotHeight);
        drawBoxes[0].setDrawingFunction(drawTimeImage1);
        drawBoxes[0].setTooltip("Image of the electric field grid: presents a slice of Ex(x,y=0,t), there the horizontal axis is time, and the vertical axis is position");
        drawBoxes[1].init(window.parentHandle(2), 0, plotHeight, plotWidth, plotHeight);
        drawBoxes[1].setDrawingFunction(drawTimeImage2);
        drawBoxes[1].setTooltip("Image of the electric field grid: presents a slice of Ey(x,y=0,t), there the horizontal axis is time, and the vertical axis is position");
        drawBoxes[2].init(window.parentHandle(2), 0, 2 * plotHeight, plotWidth, plotHeight);
        drawBoxes[2].setDrawingFunction(drawField1Plot);
        drawBoxes[2].setTooltip("Plot of the on-axis electric field in the x polarization");
        drawBoxes[3].init(window.parentHandle(2), 0, 3 * plotHeight, plotWidth, plotHeight);
        drawBoxes[3].setDrawingFunction(drawField2Plot);
        drawBoxes[3].setTooltip("Plot of the on-axis electric field in the y polarization");
        drawBoxes[4].init(window.parentHandle(2), plotWidth, 0, plotWidth, plotHeight);
        drawBoxes[4].setDrawingFunction(drawFourierImage1);
        drawBoxes[4].setTooltip("Plot of the electric field grid in momentum-frequency space: Ex(kx,ky=0,f). Is plotted on a logarithmic scale. Vertical axis is transverse momentum kx, and horizontal axis is frequency f.");
        drawBoxes[5].init(window.parentHandle(2), plotWidth, plotHeight, plotWidth, plotHeight);
        drawBoxes[5].setDrawingFunction(drawFourierImage2);
        drawBoxes[5].setTooltip("Plot of the electric field grid in momentum-frequency space: Ey(kx,ky=0,f). Is plotted on a logarithmic scale. Vertical axis is transverse momentum kx, and horizontal axis is frequency f.");
        drawBoxes[6].init(window.parentHandle(2), plotWidth, 2 * plotHeight, plotWidth, plotHeight);
        drawBoxes[6].setDrawingFunction(drawSpectrum1Plot);
        drawBoxes[6].setTooltip("Plot of the energy spectrum of the result, x-polarization.");
        drawBoxes[7].init(window.parentHandle(2), plotWidth, 3 * plotHeight, plotWidth, plotHeight);
        drawBoxes[7].setDrawingFunction(drawSpectrum2Plot);
        drawBoxes[7].setTooltip("Plot of the energy spectrum of the result, y-polarization.");
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
        checkBoxes[3].init("FP64", window.parentHandle(6), buttonWidth-2, 0, 1, 1);
        checkBoxes[3].setTooltip("Select whether the simulation is performed with 32-bit (unchecked) or 64-bit (checked) floating point numbers");
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
        miniButtons[7].init("\xf0\x9f\x8e\xaf", parentHandle, textWidth + 15, mbRow, 2, 1, buttonAddAperture);
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
        if (theSim[0].CUDAavailable) {
            pulldowns[7].addElement("CUDA");
            pulldowns[8].addElement("CUDA");
            for (int i = 1; i < theSim[0].cudaGPUCount; ++i) {
                A = Sformat("CUDA {}", i);
                pulldowns[7].addElement(A.c_str());
                pulldowns[8].addElement(A.c_str());
            }
        }
        if (theSim[0].SYCLavailable) {
            A.assign("SYCL");
            pulldowns[7].addElement(A.c_str());
            pulldowns[8].addElement(A.c_str());
            if (theSim[0].syclGPUCount > 0) {

                pulldowns[7].addElement("SYCLcpu");
                pulldowns[8].addElement("SYCLcpu");
                pulldowns[7].addElement("SYCLgpu");
                pulldowns[8].addElement("SYCLgpu");
            }
        }
        pulldowns[7].addElement("OpenMP");
        pulldowns[8].addElement("OpenMP");

        pulldowns[7].init(window.parentHandle(6), 2 + buttonWidth, 0, buttonWidth, 1);
        pulldowns[8].init(window.parentHandle(6), 4 + 2 * buttonWidth, 0, buttonWidth, 1);
        textBoxes[52].init(window.parentHandle(6), 4 + 3 * buttonWidth, 0, 1, 1);
        pulldowns[7].setTooltip("Select the primary method of calculation. The algorithm is the same, but you can run it either on a GPU or CPU depending on your machine");
        pulldowns[8].setTooltip("Select a secondary mode of calculation for offloading jobs from a batch. For example, if the pulldown to the left is set to CUDA and this one is OpenMP, and the number to the right is 2, 2 of the simulations from the batch will be performed on the CPU");

        //pulldowns[7].setLabel(-2, 0, _T("Config:"), 8, 2);
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
        std::string materialString;
        for (int i = 0; i < theDatabase.db.size(); ++i) {
            materialString = Sformat("{:2}: {}", i, std::string(theDatabase.db[i].crystalName.c_str()));
            pulldowns[3].addElement(materialString.c_str());
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
        
        //Linux search order:
        // /usr/share/LightwaveExplorer
        // working directory
        //
        //Apple search order:
        // App /Resources folder
        // working directory
#ifdef __linux__
		if (1 == readInputParametersFile(theSim.data(), theDatabase.db.data(), "/usr/share/LightwaveExplorer/DefaultValues.ini")) {
			readInputParametersFile(theSim.data(), theDatabase.db.data(), "DefaultValues.ini");
		}
#elif defined __APPLE__
		uint32_t bufferSize = 1024;
		char sysPathBuffer[1024] = { 0 };
		_NSGetExecutablePath(sysPathBuffer, &bufferSize);
        std::string sysPathFull(sysPathBuffer);
        std::string sysPathIni = sysPathFull.substr(0,sysPathFull.find_last_of("/"));
        sysPathIni.append("/../Resources/DefaultValues.ini");
		if(1 == readInputParametersFile(theSim.data(), theDatabase.db.data(), sysPathIni.c_str())){
            readInputParametersFile(theSim.data(), theDatabase.db.data(), "DefaultValues.ini");
        }
#else
		readInputParametersFile(theSim.data(), theDatabase.db.data(), "DefaultValues.ini");
#endif
        setInterfaceValuesToActiveValues();
        timeoutID = g_timeout_add(50, G_SOURCE_FUNC(updateDisplay), NULL);
        g_signal_connect(window.window, "destroy", G_CALLBACK(destroyMainWindowCallback), NULL);
        window.present();
    }
};
mainGui theGui;


bool updateDisplay() {
    theGui.console.updateFromBuffer();
    theGui.updateSlider();
    theGui.applyUpdate();
    theGui.sequence.paintSequenceText();
    return true;
}

void destroyMainWindowCallback(){
    g_source_remove(theGui.timeoutID);
}

void setInterfaceValuesToActiveValues(){
	int i = 0;
    pulse<double>* t = &theSim[0].pulse1;
    for (int k = 0; k < 2; ++k) {
        theGui.textBoxes[i++].setToDouble((*t).energy);
        theGui.textBoxes[i++].setToDouble(1e-12 * (*t).frequency);
        theGui.textBoxes[i++].setToDouble(1e-12 * (*t).bandwidth);
        theGui.textBoxes[i++].setToDouble((*t).sgOrder);
        theGui.textBoxes[i++].setToDouble((*t).cep/vPi<double>());
        theGui.textBoxes[i++].setToDouble(1e15 * (*t).delay);
        theGui.textBoxes[i++].setToDouble(1e30 * (*t).gdd);
        theGui.textBoxes[i++].setToDouble(1e45 * (*t).tod);
        theGui.textBoxes[i++].setToDouble((*t).phaseMaterial);
        theGui.textBoxes[i++].setToDouble(1e6 * (*t).phaseMaterialThickness);
        theGui.textBoxes[i++].setToDouble(1e6 * (*t).beamwaist);
        theGui.textBoxes[i++].setToDouble(1e6 * (*t).x0);
        theGui.textBoxes[i++].setToDouble(1e6 * (*t).z0);
        theGui.textBoxes[i++].setToDouble(rad2Deg<double>() * (*t).beamAngle);
        theGui.textBoxes[i++].setToDouble(rad2Deg<double>() * (*t).polarizationAngle);
        theGui.textBoxes[i++].setToDouble((*t).circularity);
        t = &theSim[0].pulse2;
    }

    theGui.pulldowns[0].setValue(theSim[0].pulse1FileType);
    theGui.pulldowns[1].setValue(theSim[0].pulse2FileType);
    theGui.pulldowns[2].setValue(theSim[0].fittingMode);
    theGui.pulldowns[3].setValue(theSim[0].materialIndex);
    
    theGui.textBoxes[i++].setToDouble(rad2Deg<double>() * asin(sin(theSim[0].crystalTheta)));
    theGui.textBoxes[i++].setToDouble(rad2Deg<double>() * asin(sin(theSim[0].crystalPhi)));
    theGui.textBoxes[i++].setToDouble(theSim[0].nonlinearAbsorptionStrength);
    theGui.textBoxes[i++].setToDouble(theSim[0].bandGapElectronVolts);
    theGui.textBoxes[i++].setToDouble(1e-12 * theSim[0].drudeGamma);
    theGui.textBoxes[i++].setToDouble(theSim[0].effectiveMass);

    if (!theSim[0].is3D) {
        theGui.textBoxes[i++].setToDouble(1e6 * theSim[0].spatialWidth);
    }
    else if (theSim[0].spatialHeight == theSim[0].spatialWidth
        || theSim[0].spatialHeight == 1 || theSim[0].spatialHeight == 0) {
        theGui.textBoxes[i++].setToDouble(1e6 * theSim[0].spatialWidth);
    }
    else {
        theGui.textBoxes[i++].overwritePrint("{};{}", (int)(1e6 * theSim[0].spatialWidth), (int)(1e6 * theSim[0].spatialHeight));
    }
    theGui.textBoxes[i++].setToDouble(1e6*theSim[0].rStep);
    theGui.textBoxes[i++].setToDouble(1e15 * theSim[0].timeSpan);
    theGui.textBoxes[i++].setToDouble(1e15 * theSim[0].tStep);
    theGui.textBoxes[i++].setToDouble(1e6 * theSim[0].crystalThickness);
    theGui.textBoxes[i++].setToDouble(1e9 * theSim[0].propagationStep);
    theGui.textBoxes[i++].setToDouble(theSim[0].batchDestination);
    theGui.textBoxes[i++].setToDouble(theSim[0].batchDestination2);
    theGui.textBoxes[i++].setToDouble((double)theSim[0].Nsims);
    theGui.textBoxes[i++].setToDouble((double)theSim[0].Nsims2);
    theGui.pulldowns[4].setValue(theSim[0].symmetryType);
    theGui.pulldowns[5].setValue(theSim[0].batchIndex);
    theGui.pulldowns[6].setValue(theSim[0].batchIndex2);
    theGui.sequence.clear();
    if (theSim[0].sequenceString.length() > 6) {
        std::string formattedSequence= theSim[0].sequenceString;
        formatSequence(formattedSequence);
        theGui.sequence.directOverwritePrintSequence(formattedSequence.c_str());
    }
    stripLineBreaks(theSim[0].field1FilePath);
    if (std::string(theSim[0].field1FilePath).compare("None") != 0) theGui.filePaths[0].overwritePrint(theSim[0].field1FilePath);
    if (std::string(theSim[0].field2FilePath).compare("None") != 0) theGui.filePaths[1].overwritePrint(theSim[0].field2FilePath);
    if (std::string(theSim[0].fittingPath).compare("None") != 0) theGui.filePaths[2].overwritePrint(theSim[0].fittingPath);
    theGui.fitCommand.clear();
    if (!(theSim[0].fittingString[0] == 'N')) {
        std::string formattedFit=theSim[0].fittingString;
        insertAfterCharacter(formattedFit,';',std::string("\n"));
        theGui.fitCommand.overwritePrint(formattedFit.c_str());
    }

    if(!theGui.loadedDefaults) theGui.filePaths[3].overwritePrint(theSim[0].outputBasePath);
    theGui.loadedDefaults = true;
}

void readParametersFromInterface() {
    int i = 0;
    pulse<double>* t = &theSim[0].pulse1;
    for (int k = 0; k < 2; ++k) {
        theGui.textBoxes[i++].valueToPointer(&(*t).energy);
        theGui.textBoxes[i++].valueToPointer(1e12, &(*t).frequency);
        theGui.textBoxes[i++].valueToPointer(1e12, &(*t).bandwidth);
        theGui.textBoxes[i++].valueToPointer(&(*t).sgOrder);
        theGui.textBoxes[i++].valueToPointer(vPi<double>(), &(*t).cep);
        theGui.textBoxes[i++].valueToPointer(1e-15, &(*t).delay);
        theGui.textBoxes[i++].valueToPointer(1e-30, &(*t).gdd);
        theGui.textBoxes[i++].valueToPointer(1e-45, &(*t).tod);
        theGui.textBoxes[i++].valueToPointer(&(*t).phaseMaterial);
        theGui.textBoxes[i++].valueToPointer(1e-6, &(*t).phaseMaterialThickness);
        theGui.textBoxes[i++].valueToPointer(1e-6, &(*t).beamwaist);
        theGui.textBoxes[i++].valueToPointer(1e-6, &(*t).x0);
        theGui.textBoxes[i++].valueToPointer(1e-6, &(*t).z0);
        theGui.textBoxes[i++].valueToPointer(deg2Rad<double>(), &(*t).beamAngle);
        theGui.textBoxes[i++].valueToPointer(deg2Rad<double>(), &(*t).polarizationAngle);
        theGui.textBoxes[i++].valueToPointer(&(*t).circularity);
        t = &theSim[0].pulse2;
    }
    
    theSim[0].pulse1FileType = theGui.pulldowns[0].getValue();
    theSim[0].pulse2FileType = theGui.pulldowns[1].getValue();
    theSim[0].fittingMode = theGui.pulldowns[2].getValue();
    theSim[0].materialIndex = theGui.pulldowns[3].getValue();

    theGui.textBoxes[i++].valueToPointer(deg2Rad<double>(), &theSim[0].crystalTheta);
    theGui.textBoxes[i++].valueToPointer(deg2Rad<double>(),  &theSim[0].crystalPhi);
    theGui.textBoxes[i++].valueToPointer(&theSim[0].nonlinearAbsorptionStrength);
    theGui.textBoxes[i++].valueToPointer(&theSim[0].bandGapElectronVolts);
    theGui.textBoxes[i++].valueToPointer(1e12, &theSim[0].drudeGamma);
    theGui.textBoxes[i++].valueToPointer(&theSim[0].effectiveMass);

    theGui.textBoxes[i++].valueToTwoPointers(1e-6, &theSim[0].spatialWidth, &theSim[0].spatialHeight);

    theGui.textBoxes[i++].valueToPointer(1e-6, &theSim[0].rStep);
    theGui.textBoxes[i++].valueToPointer(1e-15, &theSim[0].timeSpan);
    theGui.textBoxes[i++].valueToPointer(1e-15, &theSim[0].tStep);
    theGui.textBoxes[i++].valueToPointer(1e-6, &theSim[0].crystalThickness);
    theGui.textBoxes[i++].valueToPointer(1e-9, &theSim[0].propagationStep);
    theGui.textBoxes[i++].valueToPointer(&theSim[0].batchDestination);
    theGui.textBoxes[i++].valueToPointer(&theSim[0].batchDestination2);
    theGui.textBoxes[i++].valueToPointer(&theSim[0].Nsims);
    theGui.textBoxes[i++].valueToPointer(&theSim[0].Nsims2);

    theSim[0].symmetryType = theGui.pulldowns[4].getValue();
    theSim[0].batchIndex = theGui.pulldowns[5].getValue();
    theSim[0].batchIndex2 = theGui.pulldowns[6].getValue();
    theSim[0].runType = theGui.pulldowns[9].getValue();
    theGui.textBoxes[52].valueToPointer(&theSim[0].NsimsCPU);
    theSim[0].isInSequence = false;
    theGui.sequence.copyBuffer(theSim[0].sequenceString);
    stripWhiteSpace(theSim[0].sequenceString);
    if (theSim[0].sequenceString[0] != 'N' && theSim[0].sequenceString.length() > 5) theSim[0].isInSequence = true;

    theSim[0].isInFittingMode = false;
    theGui.fitCommand.copyBuffer(theSim[0].fittingString);
    stripLineBreaks(theSim[0].fittingString);

    theGui.filePaths[0].copyBuffer(theSim[0].field1FilePath);
    stripLineBreaks(theSim[0].field1FilePath);

    theGui.filePaths[1].copyBuffer(theSim[0].field2FilePath);
    stripLineBreaks(theSim[0].field2FilePath);

    theGui.filePaths[3].copyBuffer(theSim[0].outputBasePath);
    stripLineBreaks(theSim[0].outputBasePath);
    
    theGui.filePaths[2].copyBuffer(theSim[0].fittingPath);
    stripLineBreaks(theSim[0].fittingPath);

    //derived parameters and cleanup:
    theSim[0].sellmeierType = 0;
    theSim[0].axesNumber = 0;
    theSim[0].Ntime = (size_t)(minGridDimension * round(theSim[0].timeSpan / (minGridDimension * theSim[0].tStep)));
    if (theSim[0].symmetryType == 2) {
        theSim[0].is3D = true;
        theSim[0].spatialWidth = theSim[0].rStep * (minGridDimension * round(theSim[0].spatialWidth / (theSim[0].rStep * minGridDimension)));
        theSim[0].Nspace = (size_t)round(theSim[0].spatialWidth / theSim[0].rStep);
        if (theSim[0].spatialHeight > 0) {
            theSim[0].spatialHeight = theSim[0].rStep * (minGridDimension * round(theSim[0].spatialHeight / (theSim[0].rStep * minGridDimension)));
        }
        else {
            theSim[0].spatialHeight = theSim[0].spatialWidth;
        }
        theSim[0].Nspace2 = (size_t)round(theSim[0].spatialHeight / theSim[0].rStep);
    }
    else {
        theSim[0].is3D = false;
        theSim[0].Nspace2 = 1;
        theSim[0].spatialHeight = 0;
        theSim[0].spatialWidth = theSim[0].rStep * (minGridDimension * round(theSim[0].spatialWidth / (theSim[0].rStep * minGridDimension)));
        theSim[0].Nspace = (size_t)round(theSim[0].spatialWidth / theSim[0].rStep);
    }

    theSim[0].Nfreq = theSim[0].Ntime / 2 + 1;
    theSim[0].NgridC = theSim[0].Nfreq * theSim[0].Nspace * theSim[0].Nspace2;
    theSim[0].Ngrid = theSim[0].Ntime * theSim[0].Nspace * theSim[0].Nspace2;
    theSim[0].kStep = twoPi<double>() / (theSim[0].Nspace * theSim[0].rStep);
    theSim[0].fStep = 1.0 / (theSim[0].Ntime * theSim[0].tStep);
    theSim[0].Npropagation = (size_t)round(theSim[0].crystalThickness / theSim[0].propagationStep);

    theSim[0].isCylindric = theSim[0].symmetryType == 1;
    if (theSim[0].isCylindric) {
        theSim[0].pulse1.x0 = 0;
        theSim[0].pulse2.x0 = 0;
        theSim[0].pulse1.beamAngle = 0;
        theSim[0].pulse2.beamAngle = 0;
    }

    if (theSim[0].batchIndex == 0 || theSim[0].Nsims < 1) {
        theSim[0].Nsims = 1;
    }
    if (theSim[0].batchIndex2 == 0 || theSim[0].Nsims2 < 1) {
        theSim[0].Nsims2 = 1;
    }
    theSim[0].NsimsCPU = minN(theSim[0].NsimsCPU, theSim[0].Nsims * theSim[0].Nsims2);

    theSim[0].field1IsAllocated = false;
    theSim[0].field2IsAllocated = false;

    //crystal from database (database must be loaded!)
    theSim[0].crystalDatabase = theDatabase.db.data();
    theSim[0].chi2Tensor = theDatabase.db[theSim[0].materialIndex].d.data();
    theSim[0].chi3Tensor = theDatabase.db[theSim[0].materialIndex].chi3.data();
    theSim[0].nonlinearSwitches = theDatabase.db[theSim[0].materialIndex].nonlinearSwitches.data();
    theSim[0].absorptionParameters = theDatabase.db[theSim[0].materialIndex].absorptionParameters.data();
    theSim[0].sellmeierCoefficients = theDatabase.db[theSim[0].materialIndex].sellmeierCoefficients.data();
    theSim[0].sellmeierType = theDatabase.db[theSim[0].materialIndex].sellmeierType;
    theSim[0].axesNumber = theDatabase.db[theSim[0].materialIndex].axisType;
    theSim[0].progressCounter = &progressCounter;
}

int insertAfterCharacter(std::string& s, char target, std::string appended){
    for(size_t i = 0; i < s.length(); ++i){
        if(s[i] == target){
            s.insert(i+1,appended);
            i += appended.length();
        }
    }
    return 0;
}

int insertAfterCharacterExcept(std::string& s, char target, std::string appended, std::string exclude){
    bool match = false;
    for(size_t i = 0; i < s.length()-1; ++i){
        if(s[i] == target){
            match = false;
            for(int j = 0; j<exclude.length(); j++){
                if(s[i+1] == exclude[j]){
                    match = true;
                }
            }
            if(match){
                ++i;
            }
            else{
                s.insert(i+1,appended);
                i += appended.length();
            }
        }
    }
    return 0;
}

int indentForDepth(std::string& s){
    size_t depth = 0;
    std::string indent("   ");
    for(size_t i = 0; i<s.length()-1; ++i){
        if(s[i]=='{') ++depth;
        if(s[i]=='}' && depth != 0) --depth;
        if(s[i]=='\n' && s[i+1] != '}'){
            for(size_t j = 0; j<depth; j++){
                s.insert(i+1,indent);
                i += indent.length();
            }
        }
    }
    return 0;
}

int formatSequence(std::string& s){
    insertAfterCharacter(s, '>', std::string("\n"));
    insertAfterCharacterExcept(s, ')', std::string("\n"), std::string("{;"));
    insertAfterCharacter(s, ';', std::string("\n"));
    insertAfterCharacter(s, '{', std::string("\n"));
    insertAfterCharacter(s, '}', std::string("\n"));
    indentForDepth(s);
    return 0;
}

int freeSemipermanentGrids() {
    theSim[0].isGridAllocated = false;
    delete[] theSim[0].ExtOut;
    delete[] theSim[0].EkwOut;
    delete[] theSim[0].totalSpectrum;
    return 0;
}

//ifdef guards are in place to only include CUDA/SYCL when they are being used
void checkLibraryAvailability() {   
#if defined CPUONLY
    CUDAavailable = false;
    SYCLavailable = false;
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
	cudaGetDeviceCount(&theSim[0].cudaGPUCount);
	cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
	struct cudaDeviceProp activeCUDADeviceProp;
	//if (cuErr == cudaSuccess) {

	if (theSim[0].cudaGPUCount > 0) {
        theSim[0].CUDAavailable = true;
        if (theSim[0].cudaGPUCount == 1) {
            theGui.console.cPrint("CUDA found a GPU: \n", theSim[0].cudaGPUCount);
        }
        else {
            theGui.console.cPrint("CUDA found {} GPU(s): \n", theSim[0].cudaGPUCount);
        }
        for (i = 0; i < theSim[0].cudaGPUCount; ++i) {
            cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
            theGui.console.cPrint("<span color=\"#66FFFFFF\">   {}\n      Memory: {} MB\n      Multiprocessors: {}</span>\n", 
                activeCUDADeviceProp.name, 
                (int)ceil(((float)activeCUDADeviceProp.totalGlobalMem) / 1048576), 
                activeCUDADeviceProp.multiProcessorCount);
        }
	}

#else
#define solveNonlinearWaveEquationSequence solveNonlinearWaveEquationSequenceCPU
#define solveNonlinearWaveEquation solveNonlinearWaveEquationCPU
#define runDlibFitting runDlibFittingCPU
#endif

#ifndef NOSYCL
    if (isIntelRuntimeInstalled()) {
        theSim[0].SYCLavailable = true;
        char syclDeviceList[1024] = { 0 };
        size_t syclDevices = 0;
        char counts[2] = { 0 };
        readSYCLDevices(counts, syclDeviceList);
        theSim[0].syclGPUCount = (int)counts[1];
        syclDevices = (size_t)counts[0] + (size_t)counts[1];
        if (syclDevices != 0) {
            theGui.console.cPrint("{}", syclDeviceList);
        }
    }
    else {
        theGui.console.cPrint("Not using SYCL because the Intel DPC++\nCompiler Runtime is not installed.\n");
        theSim[0].SYCLavailable = false;
    }
#endif
#endif
}

static constexpr std::array<std::array<unsigned char, 3>, 256> createColormap(const int cm) {
    std::array<std::array<unsigned char, 3>, 256> colorMap{};
    float oneOver255 = 1.0f / 255.0f;
    float nval;
    for (int j = 0; j < 256; ++j) {
        switch (cm) {
        case 0:
            colorMap[j][0] = (unsigned char)j;
            colorMap[j][1] = (unsigned char)j;
            colorMap[j][2] = (unsigned char)j;
            break;
        case 1:
            nval = j * oneOver255;
            colorMap[j][0] = (unsigned char)(255 * cos(vPi<double>() * nval / 2.));
            colorMap[j][1] = (unsigned char)(255 * cos(vPi<double>() * (nval - 0.5)));
            colorMap[j][2] = (unsigned char)(255 * sin(vPi<double>() * nval / 2.));
            break;
        case 2:
            nval = j * oneOver255;
            colorMap[j][0] = (unsigned char)(255 * cos(vPi<double>() * nval / 2.));
            colorMap[j][1] = (unsigned char)(255 * cos(vPi<double>() * (nval - 0.5)));
            colorMap[j][2] = (unsigned char)(255 * sin(vPi<double>() * nval / 2.));
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
    return colorMap;
}

int drawArrayAsBitmap(cairo_t* cr, int Nx, int Ny, float* data, int cm) {
    if (Nx * Ny == 0) return 1;
    
    // creating input
    unsigned char* pixels = new unsigned char[4 * Nx * Ny]();
    if (pixels == nullptr) return 1;

    std::array<std::array<unsigned char, 3>, 256> colorMap = createColormap(cm);
    size_t Ntot = Nx * Ny;
    
    int stride = 4;
    //Find the image maximum and minimum
    float imin = data[0];
    float imax = data[0];
    for (size_t i = 1; i < Ntot; ++i) {
        if (data[i] > imax) imax = data[i];
        if (data[i] < imin) imin = data[i];
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
    delete[] pixels;
    return 0;
}

void pathFromDialogBox(GtkDialog* dialog, int response) {
    if (response == GTK_RESPONSE_ACCEPT) {
        GtkFileChooser* chooser = GTK_FILE_CHOOSER(dialog);
        GFile* file = gtk_file_chooser_get_file(chooser);
        std::string s(g_file_get_path(file));
        if (s.substr(s.length() - 4, std::string::npos) == std::string(".txt") && theGui.pathTarget == 3) {
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
        if (theSim[0].isGridAllocated) {
            freeSemipermanentGrids();
            theSim[0].isGridAllocated = false;
        }
        theGui.sequence.clear();
        theGui.fitCommand.clear();
        theGui.console.cPrint("Reading params\n");
        int readParameters = readInputParametersFile(theSim.data(), theDatabase.db.data(), path.c_str());
        theGui.console.cPrint("Fitting string: {}\n", theSim[0].fittingString);
        allocateGrids(theSim.data());
        theSim[0].isGridAllocated = true;
        if (readParameters == 61) {
            size_t extensionLoc = path.find_last_of(".");
            const std::string basePath = path.substr(0, extensionLoc);
            theGui.console.cPrint("Loading fields\n");
            loadSavedFields(theSim.data(), basePath.c_str());
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

    if (theSim[0].isGridAllocated) {
        freeSemipermanentGrids();
        theSim[0].isGridAllocated = false;
    }

    readParametersFromInterface();
    theSim[0].runType = 0;
    allocateGrids(theSim.data());
    theSim[0].isGridAllocated = true;
    theSim[0].isFollowerInSequence = false;
    theSim[0].crystalDatabase = theDatabase.db.data();

    simulationParameterSet backupSet = theSim[0];
    theSim.assign(theSim[0].Nsims * theSim[0].Nsims2 + 1, backupSet);
    configureBatchMode(theSim.data());

    std::vector<simulationParameterSet> counterData = theSim;
    simulationParameterSet* testSet = counterData.data();
    totalSteps = 0u;
    for (int j = 0; j < theSim[0].Nsims * theSim[0].Nsims2; j++) {
        if (theSim[0].isInSequence) {
            testSet[j].progressCounter = &totalSteps;
            testSet[j].runType = -1;
            solveNonlinearWaveEquationSequenceCounter(&testSet[j]);
        }
        else {
            testSet[j].progressCounter = &totalSteps;
            solveNonlinearWaveEquationCounter(&testSet[j]);
        }
    }

    //create SLURM script
    theSim[0].runType = theGui.pulldowns[9].getValue();
    int gpuType = 0;
    int gpuCount = 1;
    switch (theSim[0].runType) {
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
    double timeEstimate = saveSlurmScript(theSim.data(), gpuType, gpuCount, totalSteps);

    //create command line settings file
    theSim[0].runType = 1;
    saveSettingsFile(theSim.data());

    theGui.console.tPrint(
        "Run {} on cluster with:\nsbatch {}.slurmScript\n",
        getBasename(theSim[0].outputBasePath), getBasename(theSim[0].outputBasePath));
    theGui.console.tPrint("Estimated time to complete: {:.2} hours\n", timeEstimate, (double)totalSteps * theSim[0].Ntime * theSim[0].Nspace * theSim[0].Nspace2, theSim[0].runType);
    theSim[0].isRunning = false;
    theSim[0].isGridAllocated = false;
}

static void buttonAddSameCrystal() {
    if (theGui.textBoxes[34].valueDouble() != 0.0) {
        theGui.sequence.cPrint("plasma({},{},{},{},{},{},{},{},{})\n",
            theGui.pulldowns[3].getValue(), theGui.textBoxes[32].valueDouble(),
            theGui.textBoxes[33].valueDouble(), theGui.textBoxes[34].valueDouble(),
            theGui.textBoxes[35].valueDouble(), theGui.textBoxes[36].valueDouble(),
            theGui.textBoxes[37].valueDouble(), theGui.textBoxes[42].valueDouble(),
            theGui.textBoxes[43].valueDouble());
        theGui.sequence.paintSequenceText();
    }
    else {
        theGui.sequence.cPrint("nonlinear({},{},{},{},{})\n",
            theGui.pulldowns[3].getValue(), theGui.textBoxes[32].valueDouble(),
            theGui.textBoxes[33].valueDouble(), theGui.textBoxes[42].valueDouble(),
            theGui.textBoxes[43].valueDouble());
        theGui.sequence.paintSequenceText();
    }
}

static void buttonAddDefault() {
    theGui.sequence.cPrint("plasma(d,d,d,d,d,d,d,d,d)\n");
    theGui.sequence.paintSequenceText();
}

static void buttonAddMirror() {
    theGui.sequence.cPrint("sphericalMirror(-1.0)\n");
    theGui.sequence.paintSequenceText();
}

static void buttonAddFilter() {
    theGui.sequence.cPrint("filter(130, 20, 4, 1, 0)\n");
    theGui.sequence.paintSequenceText();
}

static void buttonAddLinear() {
    theGui.sequence.cPrint("linear({},{},{},{},{})\n",
        theGui.pulldowns[3].getValue(), theGui.textBoxes[32].valueDouble(),
        theGui.textBoxes[33].valueDouble(), theGui.textBoxes[42].valueDouble(),
        theGui.textBoxes[43].valueDouble());
    theGui.sequence.paintSequenceText();
}

static void buttonAddAperture() {
    theGui.sequence.cPrint("aperture(0.001, 2)\n");
    theGui.sequence.paintSequenceText();
}

static void buttonAddFarFieldAperture() {
    theGui.sequence.cPrint("farFieldAperture(2.0,4000,0,0)\n");
    theGui.sequence.paintSequenceText();
}

static void buttonAddForLoop() {
    theGui.sequence.cPrint("for(10,1){{\n\n}}\n");
    theGui.sequence.paintSequenceText();
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
        0.0,
        theGui.textBoxes[12].valueDouble(),
        theGui.textBoxes[13].valueDouble(),
        0.0,
        theGui.textBoxes[14].valueDouble(),
        theGui.textBoxes[15].valueDouble(),
        theGui.pulldowns[3].getValue(),
        theGui.textBoxes[32].valueDouble(),
        theGui.textBoxes[33].valueDouble());
    theGui.sequence.paintSequenceText();
}

static void buttonAddRotation() {
    theGui.sequence.cPrint("rotate(90)\n");
    theGui.sequence.paintSequenceText();
}

static void activate(GtkApplication* app, gpointer user_data) {
#if defined __linux__ || defined __APPLE__
    setlocale(LC_NUMERIC, "en_US.UTF-8");
#else
    setlocale(LC_NUMERIC, "en_US");
#endif
    theGui.activate(app);
}

void drawProgress(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    int x0 = 5;
    width -= x0;
    LweColor black(0.05, 0.05, 0.05, 0.05);
    cairo_rectangle(cr, x0, 0, width, height);
    black.setCairo(cr);
    cairo_fill(cr);

    size_t lengthEstimate = 0;
    if (!theSim[0].isInFittingMode) {
        lengthEstimate = totalSteps;
    }
    else {
        lengthEstimate = theSim[0].fittingMaxIterations;
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
    if (!theSim[0].isGridAllocated) {
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

    if (simIndex > theSim[0].Nsims * theSim[0].Nsims2) {
        simIndex = 0;
    }

    size_t cubeMiddle = theSim[0].Ntime * theSim[0].Nspace * (theSim[0].Nspace2 / 2);

    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theSim[0].tStep / 1e-15;
    sPlot.x0 = -((sPlot.dx * theSim[0].Ntime) / 2 - sPlot.dx / 2);
    sPlot.data = &theSim[0].ExtOut[simIndex * theSim[0].Ngrid * 2 + cubeMiddle + theSim[0].Ntime * theSim[0].Nspace / 2];
    sPlot.Npts = theSim[0].Ntime;
    sPlot.color = LweColor(0, 1, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ex (GV/m)";
    sPlot.unitY = 1e9;
    sPlot.makeSVG = saveSVG; // theGui.saveSVG;
    LwePlot2d(&sPlot);

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[3].copyBuffer(svgPath);
        svgPath.append("_Ex.svg");
        std::ofstream fs(svgPath);
        fs.write(sPlot.SVG.c_str(),sPlot.SVG.size());
    }
}

void drawField2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim[0].isGridAllocated) {
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

    if (simIndex > theSim[0].Nsims * theSim[0].Nsims2) {
        simIndex = 0;
    }

    size_t cubeMiddle = theSim[0].Ntime * theSim[0].Nspace * (theSim[0].Nspace2 / 2);

    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theSim[0].tStep / 1e-15;
    sPlot.x0 = -((sPlot.dx * theSim[0].Ntime) / 2 - sPlot.dx / 2);
    sPlot.data = &theSim[0].ExtOut[theSim[0].Ngrid + simIndex * theSim[0].Ngrid * 2 + cubeMiddle + theSim[0].Ntime * theSim[0].Nspace / 2];
    sPlot.Npts = theSim[0].Ntime;
    sPlot.color = LweColor(1, 0, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ey (GV/m)";
    sPlot.unitY = 1e9;
    sPlot.makeSVG = saveSVG;

    LwePlot2d(&sPlot);

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[3].copyBuffer(svgPath);
        svgPath.append("_Ey.svg");
        std::ofstream fs(svgPath);
        fs.write(sPlot.SVG.c_str(), sPlot.SVG.size());
    }
}

void drawSpectrum1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim[0].isGridAllocated) {
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
    bool logPlot = false;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = true;
    }
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim[0].Nsims * theSim[0].Nsims2) {
        simIndex = 0;
    }

    bool forceX = false;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }
    bool overlayTotal = false;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = true;
    }

    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theSim[0].fStep / 1e12;
    sPlot.data = &theSim[0].totalSpectrum[simIndex * 3 * theSim[0].Nfreq];
    sPlot.Npts = theSim[0].Nfreq;
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceY;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = LweColor(0.5, 0, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Sx (J/THz)";
    sPlot.unitY = 1.0e-12;
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    if (overlayTotal) {
        sPlot.data2 = &theSim[0].totalSpectrum[(2 + simIndex * 3) * theSim[0].Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    sPlot.makeSVG = saveSVG;

    LwePlot2d(&sPlot);

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[3].copyBuffer(svgPath);
        svgPath.append("_Sx.svg");
        std::ofstream fs(svgPath);
        fs.write(sPlot.SVG.c_str(), sPlot.SVG.size());
    }
}

void drawSpectrum2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim[0].isGridAllocated) {
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
    bool logPlot = false;
    if (theGui.checkBoxes[1].isChecked()) {
        logPlot = true;
    }
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim[0].Nsims * theSim[0].Nsims2) {
        simIndex = 0;
    }

    bool forceX = false;
    double xMin = theGui.textBoxes[48].valueDouble();
    double xMax = theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.textBoxes[50].valueDouble();
    double yMax = theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }
    bool overlayTotal = false;
    if (theGui.checkBoxes[0].isChecked()) {
        overlayTotal = true;
    }
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theSim[0].fStep / 1e12;
    sPlot.data = &theSim[0].totalSpectrum[(1 + simIndex * 3) * theSim[0].Nfreq];
    sPlot.Npts = theSim[0].Nfreq;
    sPlot.logScale = logPlot;
    sPlot.forceYmin = forceY;
    sPlot.forceYmax = forceY;
    sPlot.forcedYmax = yMax;
    if (forceY)sPlot.forcedYmin = yMin;
    sPlot.color = LweColor(1, 0, 0.5, 0.0);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Sy (J/THz)";
    sPlot.unitY = 1.0e-12;
    sPlot.forceXmax = forceX;
    sPlot.forceXmin = forceX;
    sPlot.forcedXmax = xMax;
    sPlot.forcedXmin = xMin;
    if (overlayTotal) {
        sPlot.data2 = &theSim[0].totalSpectrum[(2 + simIndex * 3) * theSim[0].Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    sPlot.makeSVG = saveSVG;

    LwePlot2d(&sPlot);

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[3].copyBuffer(svgPath);
        svgPath.append("_Sy.svg");
        std::ofstream fs(svgPath);
        fs.write(sPlot.SVG.c_str(), sPlot.SVG.size());
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
    float* plotarr2 = new float[plotSize];
    switch ((*s).dataType) {
    case 0:
        linearRemapDoubleToFloat((*s).data, (int)theSim[0].Nspace, (int)theSim[0].Ntime, plotarr2, (int)dy, (int)dx);
        drawArrayAsBitmap((*s).cr, (*s).width, (*s).height, plotarr2, 4);
        break;
    case 1:
        linearRemapZToLogFloatShift((*s).complexData, (int)theSim[0].Nspace, (int)theSim[0].Nfreq, plotarr2, (int)dy, (int)dx, (*s).logMin);
        drawArrayAsBitmap((*s).cr, (*s).width, (*s).height, plotarr2, 3);
        break;
    }
    
    delete[] plotarr2;
}

void drawTimeImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim[0].isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    imagePlotStruct sPlot;
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim[0].Nsims * theSim[0].Nsims2) {
        simIndex = 0;
    }

    size_t cubeMiddle = theSim[0].Ntime * theSim[0].Nspace * (theSim[0].Nspace2 / 2);

    sPlot.data =
        &theSim[0].ExtOut[simIndex * theSim[0].Ngrid * 2 + cubeMiddle];
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 0;
    imagePlot(&sPlot);
}

void drawTimeImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim[0].isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    imagePlotStruct sPlot;
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim[0].Nsims * theSim[0].Nsims2) {
        simIndex = 0;
    }

    size_t cubeMiddle = theSim[0].Ntime * theSim[0].Nspace * (theSim[0].Nspace2 / 2);

    sPlot.data =
        &theSim[0].ExtOut[theSim[0].Ngrid + simIndex * theSim[0].Ngrid * 2 + cubeMiddle];
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 0;
    imagePlot(&sPlot);
}

void drawFourierImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim[0].isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    imagePlotStruct sPlot;
    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim[0].Nsims * theSim[0].Nsims2) {
        simIndex = 0;
    }
    double logPlotOffset = (double)(1e-4 / (theSim[0].spatialWidth * theSim[0].timeSpan));
    if (theSim[0].is3D) {
        logPlotOffset = (double)(1e-4 / (theSim[0].spatialWidth * theSim[0].spatialHeight * theSim[0].timeSpan));
    }

    sPlot.complexData =
        &theSim[0].EkwOut[simIndex * theSim[0].NgridC * 2];
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 1;
    sPlot.logMin = logPlotOffset;
    imagePlot(&sPlot);
}

void drawFourierImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim[0].isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    imagePlotStruct sPlot;

    size_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim[0].Nsims * theSim[0].Nsims2) {
        simIndex = 0;
    }

    double logPlotOffset = (double)(1e-4 / (theSim[0].spatialWidth * theSim[0].timeSpan));
    if (theSim[0].is3D) {
        logPlotOffset = (double)(1e-4 / (theSim[0].spatialWidth * theSim[0].spatialHeight * theSim[0].timeSpan));
    }

    sPlot.complexData =
        &theSim[0].EkwOut[simIndex * theSim[0].NgridC * 2 + theSim[0].NgridC];
    sPlot.area = area;
    sPlot.cr = cr;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 1;
    sPlot.logMin = logPlotOffset;
    imagePlot(&sPlot);
}

void launchRunThread() {
    theSim[0].NsimsCPU = theGui.textBoxes[52].valueInt();
    if(!theSim[0].isRunning) std::thread(mainSimThread, theGui.pulldowns[7].getValue(), theGui.pulldowns[8].getValue(), theGui.checkBoxes[3].isChecked()).detach();
}

void launchFitThread() {
    if (!theSim[0].isRunning) std::thread(fittingThread, theGui.pulldowns[7].getValue(), theGui.checkBoxes[3].isChecked()).detach();
}

void stopButtonCallback() {
    if (theSim[0].isRunning) {
        theSim[0].cancellationCalled = true;
        for (int i = 1; i < theSim[0].Nsims; ++i) {
            theSim[0].cancellationCalled;
        }
    }
}

void independentPlotQueue(){
    theGui.requestPlotUpdate();
    theGui.applyUpdate();
}

void secondaryQueue(simulationParameterSet* cpuSims, int pulldownSelection, int pulldownSelectionPrimary, bool use64bitFloatingPoint) {
    if (theSim[0].NsimsCPU < 1) return;
    auto sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
    auto normalFunction = &solveNonlinearWaveEquationCPU;
    int assignedGPU = 0;
    bool forceCPU = 0;
    [[maybe_unused]]int SYCLitems = 0;
    if (theSim[0].syclGPUCount == 0) {
        SYCLitems = (int)theSim[0].SYCLavailable;
    }
    else {
        SYCLitems = 3;
    }
#ifndef CPUONLY
    //launch on CUDA if selected, putting in the correct GPU in multi-gpu systems
    if (pulldownSelection < theSim[0].cudaGPUCount) {
        sequenceFunction = &solveNonlinearWaveEquationSequence;
        normalFunction = &solveNonlinearWaveEquation;
        assignedGPU = pulldownSelection;
    }
    //launch on SYCL, but if the primary queue matches, default to openMP
    else if (pulldownSelection == theSim[0].cudaGPUCount && SYCLitems > 0) {
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
    else if (pulldownSelection == theSim[0].cudaGPUCount + 1 && SYCLitems > 1) {
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
    else if (pulldownSelection == theSim[0].cudaGPUCount + 2 && SYCLitems > 1) {
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
#endif
    int error = 0;
    if (theSim[0].isInSequence) {
        for (unsigned int i = 0; i < theSim[0].NsimsCPU; ++i) {
            cpuSims[i].assignedGPU = assignedGPU;
            cpuSims[i].runningOnCPU = forceCPU;
            error = sequenceFunction(&cpuSims[i]);
            if (error) break;
        }
    }
    else {
        for (unsigned int i = 0; i < theSim[0].NsimsCPU; ++i) {
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

void mainSimThread(int pulldownSelection, int secondPulldownSelection, bool use64bitFloatingPoint) {
    theSim[0].cancellationCalled = false;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();


    if (theSim[0].isGridAllocated) {
        freeSemipermanentGrids();
        theSim[0].isGridAllocated = false;
    }

    readParametersFromInterface();
    theSim[0].runType = 0;
    allocateGrids(theSim.data());
    theSim[0].isGridAllocated = true;
    theGui.requestSliderUpdate();
    theSim[0].isFollowerInSequence = false;
    theSim[0].crystalDatabase = theDatabase.db.data();
    loadPulseFiles(theSim.data());

    simulationParameterSet backupSet = theSim[0];
    theSim.assign(theSim[0].Nsims * theSim[0].Nsims2+1, backupSet);
    configureBatchMode(theSim.data());

    std::vector<simulationParameterSet> counterData = theSim;
    simulationParameterSet* testSet = counterData.data();
    totalSteps = 0;
    for (int j = 0; j < theSim[0].Nsims * theSim[0].Nsims2; j++) {
        testSet[j] = theSim[j];
        if (theSim[0].isInSequence) {
            testSet[j].progressCounter = &totalSteps;
            testSet[j].runType = -1;
            solveNonlinearWaveEquationSequenceCounter(&testSet[j]);
        }
        else {
            testSet[j].progressCounter = &totalSteps;
            solveNonlinearWaveEquationCounter(&testSet[j]);
        }
    }
    
    theSim[0].runType = 0;
    int error = 0;

    theSim[0].isRunning = true;
    progressCounter = 0;
    auto sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
    auto normalFunction = &solveNonlinearWaveEquationCPU;

    int assignedGPU = 0;
    bool forceCPU = 0;
    [[maybe_unused]]int SYCLitems = 0;
    #ifndef CPUONLY
    if (theSim[0].syclGPUCount == 0) {
        SYCLitems = (int)theSim[0].SYCLavailable;
    }
    else {
        SYCLitems = 3;
    }
    if (pulldownSelection < theSim[0].cudaGPUCount) {
        if (use64bitFloatingPoint) {
            sequenceFunction = &solveNonlinearWaveEquationSequence;
            normalFunction = &solveNonlinearWaveEquation;
        }
        else {
            sequenceFunction = &solveNonlinearWaveEquationSequenceFP32;
            normalFunction = &solveNonlinearWaveEquationFP32;
        }

        assignedGPU = pulldownSelection;
    }
    else if (pulldownSelection == theSim[0].cudaGPUCount && theSim[0].SYCLavailable) {
        if (use64bitFloatingPoint) {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
            normalFunction = &solveNonlinearWaveEquationSYCL;
        }
        else {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
            normalFunction = &solveNonlinearWaveEquationSYCLFP32;
        }

    }
    else if (pulldownSelection == theSim[0].cudaGPUCount + 1 && SYCLitems > 1) {
        forceCPU = 1;
        if (use64bitFloatingPoint) {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
            normalFunction = &solveNonlinearWaveEquationSYCL;
        }
        else {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
            normalFunction = &solveNonlinearWaveEquationSYCLFP32;
        }
    }
    else if (pulldownSelection == theSim[0].cudaGPUCount + 2 && SYCLitems > 1) {
        assignedGPU = 1;
        if (use64bitFloatingPoint) {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
            normalFunction = &solveNonlinearWaveEquationSYCL;
        }
        else {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
            normalFunction = &solveNonlinearWaveEquationSYCLFP32;
        }
    }
    else 
#endif
    {
        if (use64bitFloatingPoint) {
            sequenceFunction = &solveNonlinearWaveEquationSequenceCPU;
            normalFunction = &solveNonlinearWaveEquationCPU;
        }
        else {
            sequenceFunction = &solveNonlinearWaveEquationSequenceCPUFP32;
            normalFunction = &solveNonlinearWaveEquationCPUFP32;
        }
    }

    std::thread secondQueueThread(secondaryQueue, 
        &theSim[theSim[0].Nsims * theSim[0].Nsims2 - theSim[0].NsimsCPU], 
        secondPulldownSelection, pulldownSelection, use64bitFloatingPoint);

    for (int j = 0; j < (theSim[0].Nsims * theSim[0].Nsims2 - theSim[0].NsimsCPU); ++j) {

        theSim[j].runningOnCPU = forceCPU;
        theSim[j].assignedGPU = assignedGPU;
        if (theSim[0].isInSequence) {
            try {
                error = sequenceFunction(&theSim[j]);
            }
            catch (std::exception const& e) {
                std::string errorString=e.what();
                std::erase(errorString,'<');
                std::erase(errorString,'>');
                std::erase(errorString,'&');
                std::erase(errorString,';');
                std::erase(errorString,'{');
                std::erase(errorString,'}');
                theGui.console.tPrint("<span color=\"#FF88FF\">Simulation failed with exception:\n{}</span>\n", errorString);
            }
            if (theSim[j].memoryError != 0) {
                if (theSim[j].memoryError == -1) {
                    theGui.console.tPrint(_T("<span color=\"#FF88FF\">Not enough free GPU memory, sorry.</span>\n"), theSim[j].memoryError);
                }
                else {
                    theGui.console.tPrint(_T("<span color=\"#FF88FF\">Warning: device memory error ({}).</span>\n"), theSim[j].memoryError);
                }
            }
            if (error) break;
            theGui.requestSliderMove(j);
            independentPlotQueue();
        }
        else {
            try {
                error = normalFunction(&theSim[j]);
            } catch (std::exception const& e) {
                std::string errorString=e.what();
                std::erase(errorString,'<');
                std::erase(errorString,'>');
                std::erase(errorString,'&');
                std::erase(errorString,';');
                std::erase(errorString,'{');
                std::erase(errorString,'}');
                theGui.console.tPrint("<span color=\"#FF88FF\">Simulation failed with exception:\n{}</span>\n", errorString);
            }
            
            if (theSim[j].memoryError != 0) {
                if (theSim[j].memoryError == -1) {
                    theGui.console.tPrint(_T("<span color=\"#FF88FF\">Not enough free GPU memory, sorry.</span>\n"), theSim[j].memoryError);
                }
                else {
                    theGui.console.tPrint(_T("<span color=\"#FF88FF\">Warning: device memory error ({}).</span>\r\n"), theSim[j].memoryError);
                }
            }
            if (error) break;
        }
        if (theSim[0].cancellationCalled) {
            theGui.console.tPrint(_T("<span color=\"#FF88FF\">Warning: series cancelled, stopping after {} simulations.</span>\r\n"), j + 1);
            break;
        }
        theGui.requestSliderMove(j);
        independentPlotQueue();
    }

    if (secondQueueThread.joinable()) secondQueueThread.join();
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    if (error == 13) {
        theGui.console.tPrint(
            "<span color=\"#FF88FF\">NaN detected in grid!\nTry using a larger spatial/temporal step\nor smaller propagation step.\nSimulation was cancelled.\n</span>");
    }
    else if (error == 15) {
        theGui.console.tPrint(
            "<span color=\"#FF88FF\">Sorry, that sequence mode has been \nreplaced by the new one. Look in the \ndocumentation for more info. It is a lot \neasier to use now, and hopefully \nit won't take long to set it up. \nSorry about that!\n</span>");
    }
    else if(!error){
        theGui.console.tPrint("<span color=\"#88FFFF\">Finished after {:.4} s. </span>\n", 1e-6 *
            (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
    }
    else {
        theGui.console.tPrint(
            "<span color=\"#FF88FF\">Unhandled error! \nPlease let me know what you \nwere doing when this happened: \nnicholas.karpowicz@mpq.mpg.de</span>");
    }

    saveDataSet(theSim.data());
    deallocateGrids(theSim.data(), false);
    theSim[0].isRunning = false;
}

void fittingThread(int pulldownSelection, bool use64bitFloatingPoint) {
    theSim[0].cancellationCalled = false;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();



    if (theSim[0].isGridAllocated) {
        freeSemipermanentGrids();
        theSim[0].isGridAllocated = false;
    }

    readParametersFromInterface();
    theSim[0].runType = 0;
    allocateGrids(theSim.data());
    theSim[0].isGridAllocated = true;
    theGui.requestSliderUpdate();
    theSim[0].isFollowerInSequence = false;
    theSim[0].crystalDatabase = theDatabase.db.data();
    loadPulseFiles(theSim.data());

    simulationParameterSet backupSet = theSim[0];
    theSim.assign(theSim[0].Nsims * theSim[0].Nsims2 + 1, backupSet);
    configureBatchMode(theSim.data());
    readFittingString(theSim.data());
    if (theSim[0].Nfitting == 0) {
        theGui.console.tPrint("Couldn't interpret fitting command.\n");
        deallocateGrids(theSim.data(), false);
        return;
    }
    progressCounter = 0;
    theSim[0].progressCounter = &progressCounter;
    if (theSim[0].fittingMode == 3) {
        if (loadReferenceSpectrum(theSim[0].fittingPath, theSim.data())) {
            theGui.console.tPrint("Could not read reference file!\n");
            deallocateGrids(theSim.data(), false);
            return;
        }
    }

    theGui.console.tPrint("Fitting {} values in mode {} over {} iterations.\nRegion of interest contains {} elements\n",
        theSim[0].Nfitting, theSim[0].fittingMode, theSim[0].fittingMaxIterations, theSim[0].fittingROIsize);

    int assignedGPU = 0;
    bool forceCPU = 0;
    [[maybe_unused]]int SYCLitems = 0;
    if (theSim[0].syclGPUCount == 0) {
        SYCLitems = (int)theSim[0].SYCLavailable;
    }
    else {
        SYCLitems = 3;
    }

    auto fittingFunction = &runDlibFittingCPU;
    if (!use64bitFloatingPoint) {
        fittingFunction = &runDlibFittingCPUFP32;
    }
#ifndef CPUONLY
    if (pulldownSelection < theSim[0].cudaGPUCount) {
        fittingFunction = &runDlibFitting;
        if (!use64bitFloatingPoint)fittingFunction = &runDlibFittingFP32;
        assignedGPU = pulldownSelection;
    }
    else if (pulldownSelection == theSim[0].cudaGPUCount && theSim[0].SYCLavailable) {
        fittingFunction = &runDlibFittingSYCL;
        if (!use64bitFloatingPoint)fittingFunction = &runDlibFittingSYCLFP32;
    }
    else if (pulldownSelection == theSim[0].cudaGPUCount + 1 && SYCLitems > 1) {
        forceCPU = 1;
        fittingFunction = &runDlibFittingSYCL;
        if (!use64bitFloatingPoint)fittingFunction = &runDlibFittingSYCLFP32;
    }
    else if (pulldownSelection == theSim[0].cudaGPUCount + 2 && SYCLitems > 1) {
        assignedGPU = 1;
        fittingFunction = &runDlibFittingSYCL;
        if (!use64bitFloatingPoint)fittingFunction = &runDlibFittingSYCLFP32;
    }
#endif
    theSim[0].isRunning = true;
    theSim[0].runningOnCPU = forceCPU;
    theSim[0].assignedGPU = assignedGPU;
    fittingFunction(theSim.data());
    theSim[0].plotSim = 0;

    theGui.requestPlotUpdate();
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    theGui.console.tPrint(_T("Finished fitting after {:.4} s.\n"), 1e-6 *
        (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count()));
    theGui.console.tPrint("Fitting result:\n (index, value)\n");
    for (int i = 0; i < theSim[0].Nfitting; ++i) {
        theGui.console.tPrint("{},  {}\n", i, theSim[0].fittingResult[i]);
    }
    saveDataSet(theSim.data());
    deallocateGrids(theSim.data(), false);
    theSim[0].isRunning = false;
}

int main(int argc, char **argv){
    GtkApplication* app = gtk_application_new("nickkarpowicz.lightwave", (GApplicationFlags)0);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    return g_application_run(G_APPLICATION(app), argc, argv);
}