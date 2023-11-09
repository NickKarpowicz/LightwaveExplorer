#include "LightwaveExplorerFrontendGTK.h"

//Main data structures:
// theSim contains all of the parameters of the current simulation including grid arrays
// theDatabase is the database of crystal properties
simulationBatch theSim;
crystalDatabase theDatabase;

//Counter atomics
std::atomic_uint32_t progressCounter{};
std::atomic_uint32_t totalSteps{};

//Main class for controlling the interface
class mainGui {
    bool queueUpdate = false;
    bool queueSliderUpdate = false;
    bool queueSliderMove = false;
    bool queueInterfaceValuesUpdate = false;
    int sliderTarget = 0;
    std::mutex mutex;
public:
    std::array<LweTextBox, 56> textBoxes;
    std::unordered_map<std::string, LweButton> buttons;
    std::unordered_map<std::string, LweButton> miniButtons;
    LweConsole console{};
    LweConsole sequence{};
    LweConsole fitCommand{};
    std::array<LweTextBox,4> filePaths;
    std::unordered_map<std::string, LwePulldown> pulldowns;
    std::array<LweDrawBox,8> drawBoxes;
    LweDrawBox progressBarBox{};
    std::unordered_map<std::string, LweCheckBox> checkBoxes;
    LweSlider plotSlider{};
    LweWindow window{};
    int64_t pathTarget = 0;
    int saveSVG = 0;
    bool loadedDefaults = false;
    bool firstSYCLsimulation = true;
    ~mainGui() { std::lock_guard lastLock(mutex); }
    void requestPlotUpdate() {
        std::unique_lock<std::mutex> lock(mutex);
        queueUpdate = true;
    }
    void applyUpdate() {
        std::unique_lock<std::mutex> lock(mutex);
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
        std::unique_lock<std::mutex> lock(mutex);
        queueSliderUpdate = true;
    }
    void requestSliderMove(int target) {
        std::unique_lock<std::mutex> lock(mutex);
        queueSliderMove = true;
        sliderTarget = target;
    }
    void requestInterfaceValuesUpdate(){
        std::unique_lock<std::mutex> lock(mutex);
        queueInterfaceValuesUpdate = true;
    }

    void updateSlider() {
        std::unique_lock<std::mutex> lock(mutex);
        if (queueSliderUpdate) {
            plotSlider.setRange(0, (double)((theSim.base().Nsims * maxN(theSim.base().Nsims2, 1u) - 1)));
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
        int pulldownWidth = 6;
        int labelWidth = 6;
        int plotWidth = 1;
        int plotHeight = 1;
        int pathChars = 41;
        int colWidth = labelWidth + 2 * textWidth;
        int textCol1a = labelWidth;
        int textCol2a = textCol1a + 2 * textWidth + labelWidth;
        int textCol1b = textCol1a + textWidth;
        int textCol2b = textCol2a + textWidth;
        int buttonCol1 = textCol2a - labelWidth;
        int buttonCol2 = buttonCol1 + buttonWidth;
        int buttonCol3 = buttonCol2 + buttonWidth;
        window.init(app, "Lightwave Explorer", 1400, 800);
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
        pulldowns["pulse1"].addElement(("Synthetic"));
        pulldowns["pulse1"].addElement(("FROG"));
        pulldowns["pulse1"].addElement(("Waveform"));
        pulldowns["pulse1"].addElement("LWE .dat");
        pulldowns["pulse1"].init(parentHandle, labelWidth, 16, pulldownWidth, 1);
        filePaths[0].setLabel(0, -1, ("Data 1:"));

        filePaths[1].init(parentHandle, 0, 19, colWidth, 1);
        filePaths[1].setMaxCharacters(pathChars);
        filePaths[1].setLabel(0, -1, ("Data 2:"));
        pulldowns["pulse2"].addElement(("Synthetic"));
        pulldowns["pulse2"].addElement(("FROG"));
        pulldowns["pulse2"].addElement(("Waveform"));
        pulldowns["pulse2"].addElement("LWE .dat");
        pulldowns["pulse2"].init(parentHandle, labelWidth, 18, pulldownWidth, 1);

        filePaths[2].init(parentHandle, 0, 21, colWidth, 1);
        filePaths[2].setMaxCharacters(pathChars);
        filePaths[2].setLabel(0, -1, ("Fit data:"));
        pulldowns["fit"].addElement(("Maximize x"));
        pulldowns["fit"].addElement(("Maximize y"));
        pulldowns["fit"].addElement(("Maximize Total"));
        pulldowns["fit"].addElement(("Fit spectrum"));
        pulldowns["fit"].addElement(("Fit spectrum (log)"));
        pulldowns["fit"].init(parentHandle, labelWidth, 20, pulldownWidth, 1);

        filePaths[3].init(parentHandle, buttonCol1, 16, colWidth, 1);
        filePaths[3].setMaxCharacters(pathChars);

        drawBoxes[0].init(window.parentHandle(2), 0, 0, plotWidth, plotHeight);
        drawBoxes[0].setDrawingFunction(drawTimeImage1);
        drawBoxes[0].setTooltip(
            "Image of the electric field grid: presents a slice of Ex(x,y=0,t), "
            "there the horizontal axis is time, and the vertical axis is position"
        );
        drawBoxes[1].init(window.parentHandle(2), 0, plotHeight, plotWidth, plotHeight);
        drawBoxes[1].setDrawingFunction(drawTimeImage2);
        drawBoxes[1].setTooltip(
            "Image of the electric field grid: presents a slice of Ey(x,y=0,t), "
            "there the horizontal axis is time, and the vertical axis is position");
        drawBoxes[2].init(window.parentHandle(2), 0, 2 * plotHeight, plotWidth, plotHeight);
        drawBoxes[2].setDrawingFunction(drawField1Plot);
        drawBoxes[2].setTooltip("Plot of the on-axis electric field in the x polarization");
        drawBoxes[3].init(window.parentHandle(2), 0, 3 * plotHeight, plotWidth, plotHeight);
        drawBoxes[3].setDrawingFunction(drawField2Plot);
        drawBoxes[3].setTooltip("Plot of the on-axis electric field in the y polarization");
        drawBoxes[4].init(window.parentHandle(2), plotWidth, 0, plotWidth, plotHeight);
        drawBoxes[4].setDrawingFunction(drawFourierImage1);
        drawBoxes[4].setTooltip(
            "Plot of the electric field grid in momentum-frequency space: Ex(kx,ky=0,f). "
            "Is plotted on a logarithmic scale. Vertical axis is transverse momentum kx, "
            "and horizontal axis is frequency f.");
        drawBoxes[5].init(window.parentHandle(2), plotWidth, plotHeight, plotWidth, plotHeight);
        drawBoxes[5].setDrawingFunction(drawFourierImage2);
        drawBoxes[5].setTooltip("Plot of the electric field grid in momentum-frequency space: "
            "Ey(kx,ky=0,f). Is plotted on a logarithmic scale. Vertical axis is transverse momentum "
            "kx, and horizontal axis is frequency f.");
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

        checkBoxes["Total"].init(("Total"), window.parentHandle(4), 12, 0, 1, 1);
        checkBoxes["Total"].setTooltip("Overlay a plot of the integrated energy spectrum over the two "
            "polarization-resolved spectra");
        checkBoxes["Log"].init(("Log"), window.parentHandle(4), 13, 0, 1, 1);
        checkBoxes["Log"].setTooltip("Plot spectra on a log10 scale");
        checkBoxes["FP64"].init("FP64", window.parentHandle(6), buttonWidth-2, 0, 1, 1);
        checkBoxes["FP64"].setTooltip("Select whether the simulation is performed with 32-bit "
            "(unchecked) or 64-bit (checked) floating point numbers");
        pulldowns["propagator"].addElement(("2D Cartesian"));
        pulldowns["propagator"].addElement(("3D radial symm."));
        pulldowns["propagator"].addElement(("3D"));
        pulldowns["propagator"].addElement(("FDTD 2D"));
        pulldowns["propagator"].addElement(("FDTD 3D"));
        pulldowns["propagator"].init(parentHandle, textCol2a, 7, pulldownWidth, 1);

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
            pulldowns["batch1"].addElement(batchModeNames[i]);
            pulldowns["batch2"].addElement(batchModeNames[i]);
        }
        pulldowns["batch1"].init(parentHandle, textCol2a, 8, pulldownWidth, 1);
        pulldowns["batch2"].init(parentHandle, textCol2a, 9, pulldownWidth, 1);
        pulldowns["batch1"].setTooltip("Primary batch mode selector: the selected value "
            "from the interface will be scanned in a series of simulations, starting "
            "from the value entered on the interface, and ending with the batch target "
            "set below. The number of simulations is set by the batch steps parameter.");
        pulldowns["batch2"].setTooltip("Secondary batch mode selector - allows a 2D parameter "
            "scan. Works in the same way as the primary batch, but uses the values in "
            "the right-hand column.");

        int mbRow = 22;
        textBoxes[31].setLabel(-9 ,7,"Sequence:");
        miniButtons["addSameCrystal"].init(
            "\xf0\x9f\x93\xb8", 
            parentHandle, 
            textWidth + 1, 
            mbRow, 
            2, 
            1, 
            buttonAddSameCrystal);
        miniButtons["addDefault"].init(
            "\xe2\x99\x8a", 
            parentHandle, 
            textWidth + 3, 
            mbRow, 
            2, 
            1, 
            buttonAddDefault);
        miniButtons["addRotation"].init(
            "\xf0\x9f\x92\xab", 
            parentHandle, 
            textWidth + 5, 
            mbRow, 
            2, 
            1, 
            buttonAddRotation);
        miniButtons["addPulse"].init(
            "\xf0\x9f\x92\xa1", 
            parentHandle, 
            textWidth + 7, 
            mbRow, 
            2, 
            1, 
            buttonAddPulse);
        miniButtons["addMirror"].init(
            "\xf0\x9f\x94\x8e", 
            parentHandle, 
            textWidth + 9, 
            mbRow, 
            2, 
            1, 
            buttonAddMirror);
        miniButtons["addFilter"].init(
            "\xf0\x9f\x98\x8e", 
            parentHandle, 
            textWidth + 11, 
            mbRow, 
            2, 
            1, 
            buttonAddFilter);
        miniButtons["addLinear"].init(
            "\xf0\x9f\x93\x8f", 
            parentHandle, 
            textWidth + 13, 
            mbRow, 
            2, 
            1, 
            buttonAddLinear);
        miniButtons["addAperture"].init(
            "\xf0\x9f\x8e\xaf", 
            parentHandle, 
            textWidth + 15, 
            mbRow, 
            2, 
            1, 
            buttonAddAperture);
        miniButtons["addFarFieldAperture"].init(
            "\xe2\x9b\xb3", 
            parentHandle, 
            textWidth + 17, 
            mbRow, 
            2, 
            1, 
            buttonAddFarFieldAperture);
        miniButtons["addForLoop"].init(
            "\xf0\x9f\x94\x81", 
            parentHandle, 
            textWidth + 19, 
            mbRow, 
            2, 
            1, 
            buttonAddForLoop);
        
        miniButtons["addSameCrystal"].setTooltip(
            "Make a copy of the crystal currently entered in the interface");
        miniButtons["addDefault"].setTooltip(
            "Insert a crystal that will change with the values set on the "
            "interface, or modified during a batch calculation");
        miniButtons["addRotation"].setTooltip(
            "Rotate the polarization by a specified angle in degrees");
        miniButtons["addPulse"].setTooltip(
            "Add a new pulse to the grid; values will be set to duplicate "
            "pulse 1 as entered above");
        miniButtons["addMirror"].setTooltip(
            "Add a spherical mirror to the beam path, with radius "
            "of curvature in meters");
        miniButtons["addFilter"].setTooltip(
            "Add a spectral filter to the beam path. "
            "Parameters:\n   central frequency (THz)\n   bandwidth (THz)\n   supergaussian order\n   "
            "in-band amplitude\n   out-of-band amplitude\n");
        miniButtons["addLinear"].setTooltip(
            "Add a linear propagation through the crystal entered on the interface");
        miniButtons["addAperture"].setTooltip(
            "Add an aperture to the beam. Parameters:\n   diameter (m)\n   "
            "activation parameter\n");
        miniButtons["addFarFieldAperture"].setTooltip(
            "Filter the beam with a far-field aperture. Parameters:\n   "
            "opening angle (deg)\n   activation parameter (k)\n   x-angle (deg)\n   y-angle (deg) ");
        miniButtons["addForLoop"].setTooltip(
            "Add an empty for loop. Parameters:\n   "
            "Number of times to execute\n   Variable number in which to put the counter");
        buttons["Run"].init(("Run"), parentHandle, buttonCol3, 15, buttonWidth, 1, launchRunThread);
        buttons["Run"].setTooltip("Run the simulation as currently entered on the "
            "interface. If a sequence is entered in the sequence box below, "
            "that will execute, otherwise, a simulation on the input parameters "
            "above and to the left in a single medium will be performed.");
        buttons["Stop"].init(("Stop"), parentHandle, buttonCol2, 15, buttonWidth, 1, stopButtonCallback);
        buttons["Stop"].setTooltip("Tell a currently-running simulation to stop. "
            "It might not stop right away; it will only happen once it reaches a break point");
        buttons["Script"].init(
            ("Script"), parentHandle, buttonCol3 + 1, 17, textWidth, 1, createRunFile);
        buttons["Script"].setTooltip("Generate an input file and SLURM script for running "
            "the simulation as entered on the selected cluster");
        buttons["Fit"].init(("Fit"), parentHandle, buttonCol3, 12, buttonWidth, 1, launchFitThread);
        buttons["Fit"].setTooltip("Run the fitting routine with the above parameters. "
            "The mode is set in the pulldown next to the (optional) fitting input data file path.");
        buttons["Load"].init(("Load"), parentHandle, buttonCol1, 15, buttonWidth, 1, loadCallback);
        buttons["Load"].setTooltip("Load the results of a previous simulation run. "
            "You should select the associated .txt file. The parameters will be "
            "loaded into the interface, and the data (if it exists) will be plotted.");
        buttons["Field1Path"].init(
            ("Path"), parentHandle, textWidth, 16, textWidth, 1, openFileDialogCallback, 0);
        buttons["Field2Path"].init(
            ("Path"), parentHandle, textWidth, 18, textWidth, 1, openFileDialogCallback, (gpointer)1);
        buttons["FitPath"].init(
            ("Path"), parentHandle, textWidth, 20, textWidth, 1, openFileDialogCallback, (gpointer)2);
        buttons["Path"].init(
            ("Path"), parentHandle, buttonCol1, 17, textWidth, 1, saveFileDialogCallback, (gpointer)3);
        buttons["Path"].setTooltip("Sets the base path for the output files");
        buttons["xlim"].init(("xlim"), window.parentHandle(4), 0, 0, 1, 1, independentPlotQueue);
        buttons["xlim"].setTooltip("Apply the entered x limits to the plot. The two text "
            "boxes are for the upper and lower limits applied to the frequency axis. "
            "If they are empty, the range will include the whole grid.");
        buttons["xlim"].squeeze();
        buttons["ylim"].init(("ylim"), window.parentHandle(4), 6, 0, 1, 1, independentPlotQueue);
        buttons["ylim"].setTooltip("Apply the entered y limits to the plot. The two "
            "text boxes are for the upper and lower limits applied to the frequency "
            "axis. If they are empty, the range will include the whole grid.");
        buttons["ylim"].squeeze();
        checkBoxes["Total"].setFunction(independentPlotQueue);
        checkBoxes["Log"].setFunction(independentPlotQueue);
        buttons["SVG"].init(("SVG"), window.parentHandle(3), 5, 0, 1, 1, svgCallback);
        buttons["SVG"].setTooltip("Generate SVG files of the four line plots, with "
            "filenames based on the base path set above");
        buttons["SVG"].squeeze();
        buttons["expandPlots"].init("\xe2\x86\x94\xef\xb8\x8f",window.parentHandle(3),0,0,1,1,dataPanelCollapseCallback);
        buttons["expandPlots"].setTooltip("Collapse/expand the data entry panel");
        plotSlider.init(window.parentHandle(3), 1, 0, 3, 1);
        plotSlider.setRange(0.0, 10.0);
        plotSlider.setDigits(0);
        plotSlider.setFunction(independentPlotQueue);
        plotSlider.setArrowFunction(sliderResponseToArrows);
        sequence.init(window.parentHandle(1), 0, 0, 1, 1);
        fitCommand.init(parentHandle, buttonCol1, 13, colWidth, 2);
        console.init(parentHandle, buttonCol1, 18, colWidth, 4);
        checkLibraryAvailability();
        
        std::string A;
        if (theSim.base().CUDAavailable) {
            pulldowns["primaryHardware"].addElement("CUDA");
            pulldowns["secondaryHardware"].addElement("CUDA");
            for (int i = 1; i < theSim.base().cudaGPUCount; ++i) {
                A = Sformat("CUDA {}", i);
                pulldowns["primaryHardware"].addElement(A.c_str());
                pulldowns["secondaryHardware"].addElement(A.c_str());
            }
        }
        if (theSim.base().SYCLavailable) {
            A.assign("SYCL");
            pulldowns["primaryHardware"].addElement(A.c_str());
            pulldowns["secondaryHardware"].addElement(A.c_str());
            if (theSim.base().syclGPUCount > 0) {

                pulldowns["primaryHardware"].addElement("SYCLcpu");
                pulldowns["secondaryHardware"].addElement("SYCLcpu");
                pulldowns["primaryHardware"].addElement("SYCLgpu");
                pulldowns["secondaryHardware"].addElement("SYCLgpu");
            }
        }
#if defined _WIN32 || defined __linux__ && not defined CPUONLY
        pulldowns["primaryHardware"].addElement("C++");
        pulldowns["secondaryHardware"].addElement("C++");
#elif defined __APPLE__
        pulldowns["primaryHardware"].addElement("\xef\xa3\xbfGCD");
        pulldowns["secondaryHardware"].addElement("\xef\xa3\xbfGCD");
#endif
#if defined _WIN32 || defined __linux__
        pulldowns["primaryHardware"].addElement("OpenMP");
        pulldowns["secondaryHardware"].addElement("OpenMP");
#endif
        
        pulldowns["primaryHardware"].init(window.parentHandle(6), 2 + buttonWidth, 0, buttonWidth, 1);
        pulldowns["secondaryHardware"].init(window.parentHandle(6), 4 + 2 * buttonWidth, 0, buttonWidth, 1);
        textBoxes[52].init(window.parentHandle(6), 4 + 3 * buttonWidth, 0, 1, 1);
        pulldowns["primaryHardware"].setTooltip("Select the primary method of calculation. "
            "The algorithm is the same, but you can run it either on a "
            "GPU or CPU depending on your machine");
        pulldowns["secondaryHardware"].setTooltip("Select a secondary mode of calculation for "
            "offloading jobs from a batch. For example, if the pulldown to the "
            "left is set to CUDA and this one is OpenMP, and the number to "
            "the right is 2, 2 of the simulations from the batch will be performed on the CPU");

        //pulldowns["primaryHardware"].setLabel(-2, 0, ("Config:"), 8, 2);
        textBoxes[0].setLabel(-labelWidth, 0, ("Pulse energy (J)"));
        textBoxes[1].setLabel(-labelWidth, 0, ("Frequency (THz)"));
        textBoxes[2].setLabel(-labelWidth, 0, ("Bandwidth (THz)"));
        textBoxes[3].setLabel(-labelWidth, 0, ("SG order"));
        textBoxes[4].setLabel(-labelWidth, 0, ("CEP/\xcf\x80"));
        textBoxes[5].setLabel(-labelWidth, 0, ("Delay (fs)"));
        textBoxes[6].setLabel(-labelWidth, 0, ("GDD (fs\xc2\xb2)"));
        textBoxes[7].setLabel(-labelWidth, 0, ("TOD (fs\xc2\xb3)"));
        textBoxes[8].setLabel(-labelWidth, 0, ("Phase material"));
        textBoxes[9].setLabel(-labelWidth, 0, ("Thickness (\xce\xbcm)"));
        textBoxes[10].setLabel(-labelWidth, 0, ("Beamwaist (\xce\xbcm)"));
        textBoxes[11].setLabel(-labelWidth, 0, ("x offset (\xce\xbcm)"));
        textBoxes[12].setLabel(-labelWidth, 0, ("z offset (\xce\xbcm)"));
        textBoxes[13].setLabel(-labelWidth, 0, ("NC angle (deg)"));
        textBoxes[14].setLabel(-labelWidth, 0, ("Polarization (deg)"));
        textBoxes[15].setLabel(-labelWidth, 0, ("Circularity"));
        textBoxes[32].setLabel(-labelWidth, 0, ("Theta, phi (deg)"));
        textBoxes[34].setLabel(-labelWidth, 0, ("NL absorption"));
        textBoxes[36].setLabel(-labelWidth, 0, ("Drude: gamma, m"));
        textBoxes[38].setLabel(-labelWidth, 0, ("Max x, dx (\xce\xbcm)"));
        textBoxes[40].setLabel(-labelWidth, 0, ("Time span, dt (fs)"));
        textBoxes[42].setLabel(-labelWidth, 0, ("Max z, dz (\xce\xbcm,nm)"));
        
        textBoxes[44].setLabel(-labelWidth, 0, ("Batch end"));
        textBoxes[46].setLabel(-labelWidth, 0, ("Batch steps"));
        pulldowns["propagator"].setLabel(-labelWidth, 0, ("Propagation"));
        pulldowns["batch1"].setLabel(-labelWidth, 0, ("Batch mode"));
        pulldowns["batch2"].setLabel(-labelWidth, 0, ("Batch mode 2"));

        fitCommand.setLabel(0, -1, ("Fitting:"));

        filePaths[3].overwritePrint("TestFile");

        //read the crystal database
        std::string materialString;
        for (int i = 0; i < theDatabase.db.size(); ++i) {
            materialString = Sformat(
                "{:2}: {}", i, std::string(theDatabase.db[i].crystalName.c_str()));
            pulldowns["material"].addElement(materialString.c_str());
        }
        
        pulldowns["material"].init(parentHandle, textCol2a, 0, pulldownWidth, 1);
        pulldowns["material"].setLabel(-labelWidth, 0, ("Material"));

        pulldowns["cluster"].addElement("Cobra 1xR5k");
        pulldowns["cluster"].addElement("Cobra 2xR5k");
        pulldowns["cluster"].addElement("Cobra 1xV100");
        pulldowns["cluster"].addElement("Cobra 2xV100");
        pulldowns["cluster"].addElement("Raven 1xA100");
        pulldowns["cluster"].addElement("Raven 2xA100");
        pulldowns["cluster"].addElement("Raven 4xA100");
        pulldowns["cluster"].addElement("Raven NxA100");
        pulldowns["cluster"].init(parentHandle, buttonCol2, 17, buttonWidth + 1, 1);
        pulldowns["cluster"].setTooltip(
            "Select the cluster and GPU configuration for generating a SLURM script");
        
        //Linux search order:
        // ../share/LightwaveExplorer
        // working directory
        //
        //Apple search order:
        // App /Resources folder
        // working directory
#ifdef __linux__
        std::string homePath(std::getenv("HOME"));
        homePath.append("/.LightwaveExplorer/DefaultValues.ini");
        bool firstReadFail = (1 == 
            theSim.sCPU()->readInputParametersFile(theDatabase.db.data(), homePath));
        
        if(firstReadFail){
            char pBuf[256];
            int64_t len = sizeof(pBuf); 
            int bytes = minN(readlink("/proc/self/exe", pBuf, len), len - 1);
            if(bytes >= 0)
                pBuf[bytes] = '\0';
            std::string binPath(pBuf);
            int64_t posPath = binPath.find_last_of("/");
            std::string defaultsPath = 
                binPath.substr(0, posPath).append("/../share/LightwaveExplorer/DefaultValues.ini");

            if (1 == theSim.sCPU()->readInputParametersFile(theDatabase.db.data(), defaultsPath)) {
                theSim.sCPU()->readInputParametersFile(theDatabase.db.data(), "DefaultValues.ini");
            }
        }
        
#elif defined __APPLE__
		uint32_t bufferSize = 1024;
		char sysPathBuffer[1024] = { 0 };
		_NSGetExecutablePath(sysPathBuffer, &bufferSize);
        std::string sysPathFull(sysPathBuffer);
        std::string sysPathIni = sysPathFull.substr(0,sysPathFull.find_last_of("/"));
        sysPathIni.append("/../Resources/DefaultValues.ini");
		if(1 == theSim.sCPU()->readInputParametersFile(theDatabase.db.data(), sysPathIni.c_str())){
            theSim.sCPU()->readInputParametersFile(theDatabase.db.data(), "DefaultValues.ini");
        }
#else
		theSim.sCPU()->readInputParametersFile(theDatabase.db.data(), "DefaultValues.ini");
#endif
        setInterfaceValuesToActiveValues();

        window.connectUpdateFunction(updateDisplay);
        window.connectDestructionFunction(destroyMainWindowCallback);
        window.present();

        //if the only options are SYCL on the CPU and openMP, make the default openMP
        if (!theSim.base().CUDAavailable 
            && theSim.base().SYCLavailable 
            && theSim.base().syclGPUCount == 0) {
            pulldowns["primaryHardware"].setValue(1);
        }
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
    theGui.window.removeUpdateFunction();
   
}

void setInterfaceValuesToActiveValues(){
	int i = 0;
    pulse<double>* t = &theSim.base().pulse1;
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
        t = &theSim.base().pulse2;
    }

    theGui.pulldowns["pulse1"].setValue(theSim.base().pulse1FileType);
    theGui.pulldowns["pulse2"].setValue(theSim.base().pulse2FileType);
    theGui.pulldowns["fit"].setValue(theSim.base().fittingMode);
    theGui.pulldowns["material"].setValue(theSim.base().materialIndex);
    
    theGui.textBoxes[i++].setToDouble(rad2Deg<double>() * asin(sin(theSim.base().crystalTheta)));
    theGui.textBoxes[i++].setToDouble(rad2Deg<double>() * asin(sin(theSim.base().crystalPhi)));
    theGui.textBoxes[i++].setToDouble(theSim.base().nonlinearAbsorptionStrength);
    theGui.textBoxes[i++].setToDouble(theSim.base().bandGapElectronVolts);
    theGui.textBoxes[i++].setToDouble(1e-12 * theSim.base().drudeGamma);
    theGui.textBoxes[i++].setToDouble(theSim.base().effectiveMass);

    if (!theSim.base().is3D) {
        theGui.textBoxes[i++].setToDouble(1e6 * theSim.base().spatialWidth);
    }
    else if (theSim.base().spatialHeight == theSim.base().spatialWidth
        || theSim.base().spatialHeight == 1 || theSim.base().spatialHeight == 0) {
        theGui.textBoxes[i++].setToDouble(1e6 * theSim.base().spatialWidth);
    }
    else {
        theGui.textBoxes[i++].overwritePrint(
            "{};{}", (int)(1e6 * theSim.base().spatialWidth), 
            (int)(1e6 * theSim.base().spatialHeight));
    }
    theGui.textBoxes[i++].setToDouble(1e6*theSim.base().rStep);
    theGui.textBoxes[i++].setToDouble(1e15 * theSim.base().timeSpan);
    theGui.textBoxes[i++].setToDouble(1e15 * theSim.base().tStep);
    theGui.textBoxes[i++].setToDouble(1e6 * theSim.base().crystalThickness);
    theGui.textBoxes[i++].setToDouble(1e9 * theSim.base().propagationStep);
    theGui.textBoxes[i++].setToDouble(theSim.base().batchDestination);
    theGui.textBoxes[i++].setToDouble(theSim.base().batchDestination2);
    theGui.textBoxes[i++].setToDouble((double)theSim.base().Nsims);
    theGui.textBoxes[i++].setToDouble((double)theSim.base().Nsims2);
    theGui.pulldowns["propagator"].setValue(theSim.base().symmetryType);
    theGui.pulldowns["batch1"].setValue(theSim.base().batchIndex);
    theGui.pulldowns["batch2"].setValue(theSim.base().batchIndex2);
    theGui.sequence.clear();
    if (theSim.base().sequenceString.length() > 6) {
        std::string formattedSequence= theSim.base().sequenceString;
        formatSequence(formattedSequence);
        theGui.sequence.directOverwritePrintSequence(formattedSequence.c_str());
    }
    stripLineBreaks(theSim.base().field1FilePath);
    if (std::string(theSim.base().field1FilePath).compare("None") != 0) 
        theGui.filePaths[0].overwritePrint(theSim.base().field1FilePath);
    if (std::string(theSim.base().field2FilePath).compare("None") != 0) 
        theGui.filePaths[1].overwritePrint(theSim.base().field2FilePath);
    if (std::string(theSim.base().fittingPath).compare("None") != 0) 
        theGui.filePaths[2].overwritePrint(theSim.base().fittingPath);
    theGui.fitCommand.clear();
    if (!(theSim.base().fittingString[0] == 'N')) {
        std::string formattedFit=theSim.base().fittingString;
        insertAfterCharacter(formattedFit,';',std::string("\n"));
        theGui.fitCommand.overwritePrint(formattedFit.c_str());
    }

    if(!theGui.loadedDefaults) 
        theGui.filePaths[3].overwritePrint(theSim.base().outputBasePath);
    theGui.loadedDefaults = true;
}

void readParametersFromInterface() {
    int i = 0;
    pulse<double>* t = &theSim.base().pulse1;
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
        t = &theSim.base().pulse2;
    }
    
    theSim.base().pulse1FileType = theGui.pulldowns["pulse1"].getValue();
    theSim.base().pulse2FileType = theGui.pulldowns["pulse2"].getValue();
    theSim.base().fittingMode = theGui.pulldowns["fit"].getValue();
    theSim.base().materialIndex = theGui.pulldowns["material"].getValue();

    theGui.textBoxes[i++].valueToPointer(deg2Rad<double>(), &theSim.base().crystalTheta);
    theGui.textBoxes[i++].valueToPointer(deg2Rad<double>(),  &theSim.base().crystalPhi);
    theGui.textBoxes[i++].valueToPointer(&theSim.base().nonlinearAbsorptionStrength);
    theGui.textBoxes[i++].valueToPointer(&theSim.base().bandGapElectronVolts);
    theGui.textBoxes[i++].valueToPointer(1e12, &theSim.base().drudeGamma);
    theGui.textBoxes[i++].valueToPointer(&theSim.base().effectiveMass);

    if(theSim.base().nonlinearAbsorptionStrength != 0.0 && theSim.base().materialIndex==0){
        theGui.console.tPrint("Note:\nplasma generation doesn't work in vacuum\n");
    }

    theGui.textBoxes[i++].valueToTwoPointers(
        1e-6, &theSim.base().spatialWidth, &theSim.base().spatialHeight);

    theGui.textBoxes[i++].valueToPointer(1e-6, &theSim.base().rStep);
    theGui.textBoxes[i++].valueToPointer(1e-15, &theSim.base().timeSpan);
    theGui.textBoxes[i++].valueToPointer(1e-15, &theSim.base().tStep);
    theGui.textBoxes[i++].valueToPointer(1e-6, &theSim.base().crystalThickness);
    theGui.textBoxes[i++].valueToPointer(1e-9, &theSim.base().propagationStep);
    theGui.textBoxes[i++].valueToPointer(&theSim.base().batchDestination);
    theGui.textBoxes[i++].valueToPointer(&theSim.base().batchDestination2);
    theGui.textBoxes[i++].valueToPointer(&theSim.base().Nsims);
    theGui.textBoxes[i++].valueToPointer(&theSim.base().Nsims2);

    theSim.base().symmetryType = theGui.pulldowns["propagator"].getValue();
    theSim.base().batchIndex = theGui.pulldowns["batch1"].getValue();
    theSim.base().batchIndex2 = theGui.pulldowns["batch2"].getValue();
    //theSim.base().runType = theGui.pulldowns["cluster"].getValue();
    theGui.textBoxes[52].valueToPointer(&theSim.base().NsimsCPU);
    theSim.base().isInSequence = false;
    theGui.sequence.copyBuffer(theSim.base().sequenceString);
    stripWhiteSpace(theSim.base().sequenceString);
    if (theSim.base().sequenceString[0] != 'N' 
        && theSim.base().sequenceString.length() > 5) 
        theSim.base().isInSequence = true;

    theSim.base().isInFittingMode = false;
    theGui.fitCommand.copyBuffer(theSim.base().fittingString);
    stripLineBreaks(theSim.base().fittingString);

    theGui.filePaths[0].copyBuffer(theSim.base().field1FilePath);
    stripLineBreaks(theSim.base().field1FilePath);

    theGui.filePaths[1].copyBuffer(theSim.base().field2FilePath);
    stripLineBreaks(theSim.base().field2FilePath);

    theGui.filePaths[3].copyBuffer(theSim.base().outputBasePath);
    stripLineBreaks(theSim.base().outputBasePath);
    
    theGui.filePaths[2].copyBuffer(theSim.base().fittingPath);
    stripLineBreaks(theSim.base().fittingPath);

    //derived parameters and cleanup:
    theSim.base().sellmeierType = 0;
    theSim.base().axesNumber = 0;
    theSim.base().Ntime = (int64_t)(minGridDimension * round(theSim.base().timeSpan 
        / (minGridDimension * theSim.base().tStep)));
    if (theSim.base().symmetryType == 2) {
        theSim.base().is3D = true;
        theSim.base().isFDTD = false;
        theSim.base().spatialWidth = theSim.base().rStep 
            * (minGridDimension * round(theSim.base().spatialWidth 
                / (theSim.base().rStep * minGridDimension)));
        theSim.base().Nspace = (int64_t)round(theSim.base().spatialWidth 
            / theSim.base().rStep);
        if (theSim.base().spatialHeight > 0) {
            theSim.base().spatialHeight = theSim.base().rStep 
                * (minGridDimension * round(theSim.base().spatialHeight 
                    / (theSim.base().rStep * minGridDimension)));
        }
        else {
            theSim.base().spatialHeight = theSim.base().spatialWidth;
        }
        theSim.base().Nspace2 = (int64_t)round(theSim.base().spatialHeight 
            / theSim.base().rStep);
    }
    else if (theSim.base().symmetryType == 3) {
        theSim.base().is3D = false;
        theSim.base().isFDTD = true;
        theSim.base().Nspace2 = 1;
        theSim.base().spatialHeight = 0;
        theSim.base().spatialWidth = theSim.base().rStep 
            * (minGridDimension * round(theSim.base().spatialWidth 
                / (theSim.base().rStep * minGridDimension)));
        theSim.base().Nspace = (int64_t)round(theSim.base().spatialWidth 
            / theSim.base().rStep);
    }
    else if (theSim.base().symmetryType == 4) {
        theSim.base().is3D = true;
        theSim.base().isFDTD = true;
        theSim.base().spatialWidth = theSim.base().rStep 
            * (minGridDimension * round(theSim.base().spatialWidth 
                / (theSim.base().rStep * minGridDimension)));
        theSim.base().Nspace = (int64_t)round(theSim.base().spatialWidth 
            / theSim.base().rStep);
        if (theSim.base().spatialHeight > 0) {
            theSim.base().spatialHeight = theSim.base().rStep 
                * (minGridDimension * round(theSim.base().spatialHeight 
                    / (theSim.base().rStep * minGridDimension)));
        }
        else {
            theSim.base().spatialHeight = theSim.base().spatialWidth;
        }
        theSim.base().Nspace2 = (int64_t)round(theSim.base().spatialHeight 
            / theSim.base().rStep);
    }
    else {
        theSim.base().is3D = false;
        theSim.base().isFDTD = false;
        theSim.base().Nspace2 = 1;
        theSim.base().spatialHeight = 0;
        theSim.base().spatialWidth = theSim.base().rStep 
            * (minGridDimension * round(theSim.base().spatialWidth 
                / (theSim.base().rStep * minGridDimension)));
        theSim.base().Nspace = (int64_t)round(theSim.base().spatialWidth 
            / theSim.base().rStep);
    }

    theSim.base().Nfreq = theSim.base().Ntime / 2 + 1;
    theSim.base().NgridC = theSim.base().Nfreq * theSim.base().Nspace * theSim.base().Nspace2;
    theSim.base().Ngrid = theSim.base().Ntime * theSim.base().Nspace * theSim.base().Nspace2;
    theSim.base().kStep = twoPi<double>() / (theSim.base().Nspace * theSim.base().rStep);
    theSim.base().fStep = 1.0 / (theSim.base().Ntime * theSim.base().tStep);
    theSim.base().Npropagation = (int64_t)round(theSim.base().crystalThickness 
        / theSim.base().propagationStep);

    theSim.base().isCylindric = theSim.base().symmetryType == 1;
    if (theSim.base().isCylindric) {
        theSim.base().pulse1.x0 = 0;
        theSim.base().pulse2.x0 = 0;
        theSim.base().pulse1.beamAngle = 0;
        theSim.base().pulse2.beamAngle = 0;
    }

    if (theSim.base().batchIndex == 0 || theSim.base().Nsims < 1) {
        theSim.base().Nsims = 1;
    }
    if (theSim.base().batchIndex2 == 0 || theSim.base().Nsims2 < 1) {
        theSim.base().Nsims2 = 1;
    }
    theSim.base().NsimsCPU = minN(theSim.base().NsimsCPU, 
        theSim.base().Nsims * theSim.base().Nsims2);

    theSim.base().field1IsAllocated = false;
    theSim.base().field2IsAllocated = false;

    //crystal from database (database must be loaded!)
    theSim.base().crystalDatabase = theDatabase.db.data();
    theSim.base().chi2Tensor = 
        theDatabase.db[theSim.base().materialIndex].d.data();
    theSim.base().chi3Tensor = 
        theDatabase.db[theSim.base().materialIndex].chi3.data();
    theSim.base().nonlinearSwitches = 
        theDatabase.db[theSim.base().materialIndex].nonlinearSwitches.data();
    theSim.base().absorptionParameters = 
        theDatabase.db[theSim.base().materialIndex].absorptionParameters.data();
    theSim.base().sellmeierCoefficients = 
        theDatabase.db[theSim.base().materialIndex].sellmeierCoefficients.data();
    theSim.base().sellmeierType = 
        theDatabase.db[theSim.base().materialIndex].sellmeierType;
    theSim.base().axesNumber = 
        theDatabase.db[theSim.base().materialIndex].axisType;
    theSim.base().progressCounter = &progressCounter;

    theSim.base().runType = runTypes::normal;
    theSim.base().isFollowerInSequence = false;
    theSim.base().crystalDatabase = theDatabase.db.data();
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

int insertAfterCharacterExcept(
    std::string& s, 
    char target, 
    std::string appended, 
    std::string exclude){
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
    int depth = 0;
    std::string indent("   ");
    for (size_t i = 0; i < s.length() - 1; ++i) {
        if (s[i] == '{') ++depth;
        if (s[i] == '}' && depth != 0) --depth;
        if (s[i] == '\n' && s[i + 1] != '}') {
            for (size_t j = 0; j < depth; j++) {
                s.insert(i + 1, indent);
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


//ifdef guards are in place to only include CUDA/SYCL when they are being used
void checkLibraryAvailability() {   
#if defined CPUONLY
    theSim.base().CUDAavailable = false;
    theSim.base().SYCLavailable = false;
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
	cudaGetDeviceCount(&theSim.base().cudaGPUCount);
	cudaError_t cuErr = cudaGetDevice(&CUDAdevice);
	struct cudaDeviceProp activeCUDADeviceProp;
	//if (cuErr == cudaSuccess) {

	if (theSim.base().cudaGPUCount > 0) {
        theSim.base().CUDAavailable = true;
        if (theSim.base().cudaGPUCount == 1) {
            theGui.console.cPrint("CUDA found a GPU: \n", theSim.base().cudaGPUCount);
        }
        else {
            theGui.console.cPrint("CUDA found {} GPU(s): \n", theSim.base().cudaGPUCount);
        }
        for (i = 0; i < theSim.base().cudaGPUCount; ++i) {
            cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
            theGui.console.cPrint("<span color=\"#66FFFFFF\">   {}\n      "
                "Memory: {} MB\n      Multiprocessors: {}</span>\n", 
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
    bool isIntelRuntimeInstalled = true;
#ifdef _WIN32
    wchar_t loadBuffer[1024];
    DWORD envcount = GetEnvironmentVariableW(L"INTEL_DEV_REDIST", loadBuffer, 16);
    if(envcount==0) isIntelRuntimeInstalled = false;  
#endif
    if (isIntelRuntimeInstalled) {
        theSim.base().SYCLavailable = true;
        char syclDeviceList[1024] = { 0 };
        int64_t syclDevices = 0;
        char counts[2] = { 0 };
        readSYCLDevices(counts, syclDeviceList);
        theSim.base().syclGPUCount = (int)counts[1];
        syclDevices = (int64_t)counts[0] + (int64_t)counts[1];
        if (syclDevices != 0) {
            theGui.console.cPrint("{}", syclDeviceList);
        }
    }
    else {
        theGui.console.cPrint("Not using SYCL because the Intel DPC++\n"
            "Compiler Runtime is not installed.\n");
        theSim.base().SYCLavailable = false;
    }
#endif
#endif
}

void pathFromDialogBox(GtkDialog* dialog, int response) {
    std::unique_lock<std::mutex> GTKlock(GTKmutex);
    if (response == GTK_RESPONSE_ACCEPT) {
        GtkFileChooser* chooser = GTK_FILE_CHOOSER(dialog);
        GFile* file = gtk_file_chooser_get_file(chooser);
        std::string s(g_file_get_path(file));
        GTKlock.unlock();
        if (s.substr(s.length() - 4, std::string::npos) == std::string(".txt") 
            && theGui.pathTarget == 3) {
            theGui.filePaths[theGui.pathTarget].overwritePrint("{}", s.substr(0,s.length()-4));
        }
        else {
            theGui.filePaths[theGui.pathTarget].overwritePrint("{}", s);
        } 
        GTKlock.lock();
    }
    g_object_unref(dialog);
}

void openFileDialog(int64_t pathTarget) {
    std::unique_lock<std::mutex> GTKlock(GTKmutex);
    theGui.pathTarget = pathTarget;
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new(
        "Open File", 
        theGui.window.windowHandle(), 
        GTK_FILE_CHOOSER_ACTION_OPEN, "Ok", "Cancel");
    g_signal_connect(fileC, "response", G_CALLBACK(pathFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
}

void openFileDialogCallback(GtkWidget* widget, gpointer pathTarget) {
    std::unique_lock<std::mutex> GTKlock(GTKmutex);
    theGui.pathTarget = (int64_t)pathTarget;
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new(
        "Open File", 
        theGui.window.windowHandle(), 
        GTK_FILE_CHOOSER_ACTION_OPEN, 
        "Ok", 
        "Cancel");
    g_signal_connect(fileC, "response", G_CALLBACK(pathFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
}

void loadFromDialogBox(GtkDialog* dialog, int response) {
    std::unique_lock<std::mutex> GTKlock(GTKmutex);
    if (response == GTK_RESPONSE_ACCEPT) {
        GtkFileChooser* chooser = GTK_FILE_CHOOSER(dialog);
        GFile* file = gtk_file_chooser_get_file(chooser);
        std::string path(g_file_get_path(file));
        GTKlock.unlock();
        theGui.sequence.clear();
        theGui.fitCommand.clear();
        int readParameters = 
            theSim.base().readInputParametersFile(theDatabase.db.data(), path);
        
        if (readParameters == 61) {
            int64_t extensionLoc = path.find_last_of(".");
            const std::string basePath = path.substr(0, extensionLoc);
            std::string testPath = basePath;
            if(std::filesystem::exists(testPath.append("_Ext.dat"))){
                theSim.configure();
                theSim.base().loadSavedFields(basePath);
            }
            else{
                theSim.configure(false);
                theSim.base().isGridAllocated = false;
            }

            setInterfaceValuesToActiveValues();
            theGui.requestSliderUpdate();
            theGui.requestPlotUpdate();
        }
        GTKlock.lock();
    }
    g_object_unref(dialog);
}

void loadCallback(GtkWidget* widget, gpointer pathTarget) {
    std::unique_lock<std::mutex> GTKlock(GTKmutex);
    theGui.pathTarget = (int64_t)pathTarget;
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new(
        "Open File", 
        theGui.window.windowHandle(), 
        GTK_FILE_CHOOSER_ACTION_OPEN, 
        "Ok", 
        "Cancel");
    g_signal_connect(fileC, "response", G_CALLBACK(loadFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
}

void svgCallback() {
    theGui.saveSVG = 4;
    theGui.requestPlotUpdate();
}

void dataPanelCollapseCallback(){
    theGui.window.toggleSettingsPanel();
}
void saveFileDialogCallback(GtkWidget* widget, gpointer pathTarget) {
    theGui.pathTarget = (int64_t)pathTarget;
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
    std::unique_lock<std::mutex> GTKlock(GTKmutex);
    GtkFileChooserNative* fileC = gtk_file_chooser_native_new(
        "Save File", 
        theGui.window.windowHandle(), 
        GTK_FILE_CHOOSER_ACTION_SAVE, 
        "Ok", 
        "Cancel");
    g_signal_connect(fileC, "response", G_CALLBACK(pathFromDialogBox), NULL);
    gtk_native_dialog_show(GTK_NATIVE_DIALOG(fileC));
#endif

}

void createRunFile() {
    readParametersFromInterface();
    theSim.base().runType = runTypes::normal;
    theSim.base().isGridAllocated = false;
    theSim.base().isFollowerInSequence = false;
    theSim.base().crystalDatabase = theDatabase.db.data();
    theSim.configureCounter();

    std::vector<simulationParameterSet> counterVector = theSim.getParameterVector();
    totalSteps = 0;
    for (int64_t j = 0; j < theSim.base().Nsims * theSim.base().Nsims2; j++) {
        if (theSim.base().isInSequence) {
            counterVector[j].progressCounter = &totalSteps;
            counterVector[j].runType = runTypes::counter;
            solveNonlinearWaveEquationSequenceCounter(&counterVector[j]);
        }
        else {
            counterVector[j].progressCounter = &totalSteps;
            solveNonlinearWaveEquationCounter(&counterVector[j]);
        }
    }

    //create SLURM script
    int cluster = theGui.pulldowns["cluster"].getValue();
    std::string gpuType("ERROR");
    int gpuCount = 1;
    bool arrayMode = false;
    switch (cluster) {
    case 0:
        gpuType.assign("rtx5000");
        gpuCount = 1;
        break;
    case 1:
        gpuType.assign("rtx5000");
        gpuCount = 2;
        break;
    case 2:
        gpuType.assign("v100");
        gpuCount = 1;
        break;
    case 3:
        gpuType.assign("v100");
        gpuCount = 2;
        break;
    case 4:
        gpuType.assign("a100");
        gpuCount = 1;
        break;
    case 5:
        gpuType.assign("a100");
        gpuCount = 2;
        break;
    case 6:
        gpuType.assign("a100");
        gpuCount = 4;
        break;
    case 7:
        gpuType.assign("a100");
        gpuCount = 1;
        arrayMode = true;
        break;
    }
    double timeEstimate = theSim.sCPU()->saveSlurmScript(gpuType, gpuCount, arrayMode, totalSteps);

    //create command line settings file
    
    
    if(arrayMode){
        int simIndex = 0;
        auto& params = theSim.getParameterVector();
        theSim.sCPU()->saveSettingsFile();
        for(int i = 0; i<theSim.sCPU()->Nsims2; ++i){
            for(int j = 0; j<theSim.sCPU()->Nsims; ++j){
                simulationParameterSet arraySim = params[i*theSim.sCPU()->Nsims + j];
                arraySim.Nsims = 1;
                arraySim.Nsims2 = 1;
                arraySim.outputBasePath.append(Sformat("{:04d}",simIndex++));
                arraySim.runType = runTypes::cluster;
                arraySim.batchIndex = 0;
                arraySim.batchIndex2 = 0;
                arraySim.runType = runTypes::cluster;
                arraySim.saveSettingsFile();
            }
        }
        int jobID = 0;
    }
    else{
        theSim.base().runType = runTypes::cluster;
        theSim.sCPU()->saveSettingsFile();
    }
    
    

    theGui.console.tPrint(
        "Run {} on cluster with:\nsbatch {}.slurmScript\n",
        getBasename(theSim.base().outputBasePath), getBasename(theSim.base().outputBasePath));
    theGui.console.tPrint(
        "Upper estimate time to complete: {:.2} hours\n", 
        timeEstimate);
    theSim.base().isRunning = false;
    theSim.base().isGridAllocated = false;
}

static void buttonAddSameCrystal() {
    if (theGui.textBoxes[34].valueDouble() != 0.0) {
        theGui.sequence.cPrint("plasma({},{},{},{},{},{},{},{},{})\n",
            theGui.pulldowns["material"].getValue(), theGui.textBoxes[32].valueDouble(),
            theGui.textBoxes[33].valueDouble(), theGui.textBoxes[34].valueDouble(),
            theGui.textBoxes[35].valueDouble(), theGui.textBoxes[36].valueDouble(),
            theGui.textBoxes[37].valueDouble(), theGui.textBoxes[42].valueDouble(),
            theGui.textBoxes[43].valueDouble());
        theGui.sequence.paintSequenceText();
    }
    else {
        theGui.sequence.cPrint("nonlinear({},{},{},{},{})\n",
            theGui.pulldowns["material"].getValue(), theGui.textBoxes[32].valueDouble(),
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
        theGui.pulldowns["material"].getValue(), theGui.textBoxes[32].valueDouble(),
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
    theGui.sequence.cPrint("addPulse({},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})\n",
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
        theGui.pulldowns["material"].getValue(),
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

    int64_t lengthEstimate = 0;
    if (!theSim.base().isInFittingMode) {
        lengthEstimate = totalSteps;
    }
    else {
        lengthEstimate = theSim.base().fittingMaxIterations;
    }

    double newFraction =  0.0;
    if(lengthEstimate) newFraction = 
        minN(1.0, ((double)progressCounter) / (double)lengthEstimate);
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
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LwePlot sPlot;

    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }

    int64_t simIndex = maxN(0,theGui.plotSlider.getIntValue());

    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = theSim.base().Ntime * theSim.base().Nspace * (theSim.base().Nspace2 / 2);

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[3].copyBuffer(svgPath);
        svgPath.append("_Ex.svg");
        sPlot.SVGPath = svgPath;
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theSim.base().tStep / 1e-15;
    sPlot.x0 = -((sPlot.dx * theSim.base().Ntime) / 2 - sPlot.dx / 2);
    sPlot.data = &theSim.base().ExtOut[
        simIndex * theSim.base().Ngrid * 2 + cubeMiddle + theSim.base().Ntime * theSim.base().Nspace / 2];
    sPlot.Npts = theSim.base().Ntime;
    sPlot.color = LweColor(0, 1, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ex (GV/m)";
    sPlot.unitY = 1e9;
    std::unique_lock dataLock(theSim.mutexes.at(simIndex),std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.plot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        theGui.requestPlotUpdate();
    }
}

void drawField2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LwePlot sPlot;

    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }

    int64_t simIndex = maxN(0,theGui.plotSlider.getIntValue());

    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = 
        theSim.base().Ntime * theSim.base().Nspace * (theSim.base().Nspace2 / 2);

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[3].copyBuffer(svgPath);
        svgPath.append("_Ey.svg");
        sPlot.SVGPath = svgPath;
    }


    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theSim.base().tStep / 1e-15;
    sPlot.x0 = -((sPlot.dx * theSim.base().Ntime) / 2 - sPlot.dx / 2);
    sPlot.data = 
        &theSim.base().ExtOut[
        theSim.base().Ngrid + simIndex * theSim.base().Ngrid * 2 
            + cubeMiddle + theSim.base().Ntime * theSim.base().Nspace / 2];
    sPlot.Npts = theSim.base().Ntime;
    sPlot.color = LweColor(1, 0, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ey (GV/m)";
    sPlot.unitY = 1e9;
    std::unique_lock dataLock(theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.plot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        theGui.requestPlotUpdate();
    }
}

void drawSpectrum1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LwePlot sPlot;
    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkBoxes["Log"].isChecked()) {
        logPlot = true;
    }
    int64_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
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
    if (theGui.checkBoxes["Total"].isChecked()) {
        overlayTotal = true;
    }

    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[3].copyBuffer(svgPath);
        svgPath.append("_Sx.svg");
        sPlot.SVGPath = svgPath;
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theSim.base().fStep / 1e12;
    sPlot.data = &theSim.base().totalSpectrum[simIndex * 3 * theSim.base().Nfreq];
    sPlot.Npts = theSim.base().Nfreq;
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
        sPlot.data2 = &theSim.base().totalSpectrum[(2 + simIndex * 3) * theSim.base().Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    std::unique_lock dataLock(theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.plot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        theGui.requestPlotUpdate();
    }
}

void drawSpectrum2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }

    LwePlot sPlot;

    bool saveSVG = theGui.saveSVG > 0;
    if (saveSVG) {
        theGui.saveSVG--;
    }
    bool logPlot = false;
    if (theGui.checkBoxes["Log"].isChecked()) {
        logPlot = true;
    }
    int64_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
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
    if (theGui.checkBoxes["Total"].isChecked()) {
        overlayTotal = true;
    }
    if (saveSVG) {
        std::string svgPath;
        theGui.filePaths[3].copyBuffer(svgPath);
        svgPath.append("_Sy.svg");
        sPlot.SVGPath = svgPath;
    }

    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theSim.base().fStep / 1e12;
    sPlot.data = &theSim.base().totalSpectrum[(1 + simIndex * 3) * theSim.base().Nfreq];
    sPlot.Npts = theSim.base().Nfreq;
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
        sPlot.data2 = &theSim.base().totalSpectrum[(2 + simIndex * 3) * theSim.base().Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    std::unique_lock dataLock(theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.plot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        theGui.requestPlotUpdate();
    }
}

void drawTimeImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;
    int64_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = theSim.base().Ntime * theSim.base().Nspace * (theSim.base().Nspace2 / 2);

    sPlot.data =
        &theSim.base().ExtOut[simIndex * theSim.base().Ngrid * 2 + cubeMiddle];
    sPlot.dataXdim = theSim.base().Ntime;
    sPlot.dataYdim = theSim.base().Nspace;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.colorMap = 4;
    sPlot.dataType = 0;
    std::unique_lock dataLock(theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.imagePlot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        theGui.requestPlotUpdate();
    }
}

void drawTimeImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;
    int64_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = theSim.base().Ntime * theSim.base().Nspace * (theSim.base().Nspace2 / 2);

    sPlot.data =
    &theSim.base().ExtOut[theSim.base().Ngrid + simIndex * theSim.base().Ngrid * 2 + cubeMiddle];
    sPlot.dataYdim = theSim.base().Nspace;
    sPlot.dataXdim = theSim.base().Ntime;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.colorMap = 4;
    sPlot.dataType = 0;
    std::unique_lock dataLock(theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.imagePlot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        theGui.requestPlotUpdate();
    }
}

void drawFourierImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;
    int64_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }
    double logPlotOffset = (double)(1e-4 / (theSim.base().spatialWidth * theSim.base().timeSpan));
    if (theSim.base().is3D) {
        logPlotOffset = 
            (double)(1e-4 
                / (theSim.base().spatialWidth * theSim.base().spatialHeight * theSim.base().timeSpan));
    }
    sPlot.complexData =
        &theSim.base().EkwOut[simIndex * theSim.base().NgridC * 2];
    sPlot.dataXdim = theSim.base().Nfreq;
    sPlot.dataYdim = theSim.base().Nspace;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 1;
    sPlot.colorMap = 3;
    sPlot.logMin = logPlotOffset;
    std::unique_lock dataLock(theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.imagePlot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        theGui.requestPlotUpdate();
    }
}

void drawFourierImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;

    int64_t simIndex = maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }

    double logPlotOffset = (double)(1e-4 / (theSim.base().spatialWidth * theSim.base().timeSpan));
    if (theSim.base().is3D) {
        logPlotOffset = (double)(1e-4 
            / (theSim.base().spatialWidth * theSim.base().spatialHeight * theSim.base().timeSpan));
    }
    sPlot.complexData =
        &theSim.base().EkwOut[simIndex * theSim.base().NgridC * 2 + theSim.base().NgridC];
    sPlot.dataXdim = theSim.base().Nfreq;
    sPlot.dataYdim = theSim.base().Nspace;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.colorMap = 3;
    sPlot.dataType = 1;
    sPlot.logMin = logPlotOffset;
    std::unique_lock dataLock(theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.imagePlot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        theGui.requestPlotUpdate();
    }
}

void launchRunThread() {
    theSim.base().NsimsCPU = theGui.textBoxes[52].valueInt();
    if(!theSim.base().isRunning) 
        std::thread(
            mainSimThread, 
            theGui.pulldowns["primaryHardware"].getValue(), 
            theGui.pulldowns["secondaryHardware"].getValue(), 
            theGui.checkBoxes["FP64"].isChecked()
        ).detach();
}

void launchFitThread() {
    if (!theSim.base().isRunning) 
        std::thread(
            fittingThread, 
            theGui.pulldowns["primaryHardware"].getValue(), 
            theGui.checkBoxes["FP64"].isChecked()
        ).detach();
}

void stopButtonCallback() {
    if (theSim.base().isRunning) {
        theSim.base().cancellationCalled = true;
        for (int i = 1; i < theSim.base().Nsims; ++i) {
            theSim.base().cancellationCalled = true;
        }
    }
}

void independentPlotQueue(){
    theGui.requestPlotUpdate();
}

bool sliderResponseToArrows(GtkWidget* widget, guint keyValue, guint keyCode, GdkModifierType state, GtkEventControllerKey* eventController) {
    double value = gtk_range_get_value(GTK_RANGE(widget));
    switch (keyValue) {
    case GDK_KEY_Left:
        value -= 1.0;
        theGui.plotSlider.setValue(static_cast<int>(value));
        break;
    case GDK_KEY_Right:
        value += 1.0;
        theGui.plotSlider.setValue(static_cast<int>(value));
        break;
    default: break;
    }
    theGui.requestPlotUpdate();
    return true;
}

void secondaryQueue(
    int64_t cpuSimsIndex, 
    int pulldownSelection, 
    int pulldownSelectionPrimary, 
    bool use64bitFloatingPoint) {
    if (theSim.base().NsimsCPU < 1) return;
    auto sequenceFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSequenceCPU : &solveNonlinearWaveEquationSequenceCPUFP32;
    auto normalFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationCPU : &solveNonlinearWaveEquationCPUFP32;
    int assignedGPU = 0;
    bool forceCPU = false;
    bool useOpenMP = false;
#ifdef CPUONLY
    useOpenMP = true;
#endif
    [[maybe_unused]]int SYCLitems = 0;
    if (theSim.base().syclGPUCount == 0) {
        SYCLitems = (int)theSim.base().SYCLavailable;
    }
    else {
        SYCLitems = 3;
    }
#ifndef CPUONLY
    //launch on CUDA if selected, putting in the correct GPU in multi-gpu systems
    if (pulldownSelection < theSim.base().cudaGPUCount) {
        #ifndef NOCUDA
        sequenceFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSequence : &solveNonlinearWaveEquationSequenceFP32;
        normalFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquation : &solveNonlinearWaveEquationFP32;
        assignedGPU = pulldownSelection;
        #endif
    }
    #ifndef NOSYCL
    //launch on SYCL, but if the primary queue matches, default to openMP
    else if (pulldownSelection == theSim.base().cudaGPUCount && SYCLitems > 0) {
        if (pulldownSelection == pulldownSelectionPrimary) {
            theGui.console.tPrint(
                "Sorry, can't run two identical SYCL queues.\n"
                "Defaulting to OpenMP for the secondary queue.\n");
            sequenceFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSequenceCPU : &solveNonlinearWaveEquationSequenceCPUFP32;
            normalFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationCPU : &solveNonlinearWaveEquationCPUFP32;
        }
        else {
            sequenceFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSequenceSYCL : &solveNonlinearWaveEquationSequenceSYCLFP32;
            normalFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSYCL : &solveNonlinearWaveEquationSYCLFP32;
        }
    }
    #endif
    else if (pulldownSelection == theSim.base().cudaGPUCount + 1 && SYCLitems > 1) {
        if (pulldownSelection == pulldownSelectionPrimary) {
            theGui.console.tPrint(
                "Sorry, can't run two identical SYCL queues.\n"
                "Defaulting to OpenMP for the secondary queue.\n");
            sequenceFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSequenceCPU : &solveNonlinearWaveEquationSequenceCPUFP32;
            normalFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationCPU : &solveNonlinearWaveEquationCPUFP32;
        }
    #ifndef NOSYCL
        else {
            forceCPU = 1;
            sequenceFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSequenceSYCL : &solveNonlinearWaveEquationSequenceSYCLFP32;
            normalFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSYCL : &solveNonlinearWaveEquationSYCLFP32;
        }
    #endif
    }
    else if (pulldownSelection == theSim.base().cudaGPUCount + 2 && SYCLitems > 1) {
        if (pulldownSelection == pulldownSelectionPrimary) {
            theGui.console.tPrint(
                "Sorry, can't run two identical SYCL queues.\n"
                "Defaulting to OpenMP for the secondary queue.\n");
            sequenceFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSequenceCPU : &solveNonlinearWaveEquationSequenceCPUFP32;
            normalFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationCPU : &solveNonlinearWaveEquationCPUFP32;
        }
        #ifndef NOSYCL
        else {
            assignedGPU = 1;
            sequenceFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSequenceSYCL : &solveNonlinearWaveEquationSequenceSYCLFP32;
            normalFunction = use64bitFloatingPoint ? &solveNonlinearWaveEquationSYCL : &solveNonlinearWaveEquationSYCLFP32;
        }
        #endif
    }
    else if (pulldownSelection == (theSim.base().cudaGPUCount + SYCLitems + 1)) {
        useOpenMP = true;
    }
#endif
    int error = 0;
    if (theSim.base().isInSequence) {
        for (int64_t i = cpuSimsIndex; i < theSim.base().NsimsCPU+cpuSimsIndex; ++i) {
            theSim.sCPU()[i].assignedGPU = assignedGPU;
            theSim.sCPU()[i].runningOnCPU = forceCPU;
            theSim.sCPU()[i].useOpenMP = useOpenMP;
            std::lock_guard dataLock(theSim.mutexes.at(i));
            error = sequenceFunction(&theSim.sCPU()[i]);
            if (error) break;
        }
    }
    else {
        for (int64_t i = cpuSimsIndex; i < theSim.base().NsimsCPU+cpuSimsIndex; ++i) {
            theSim.sCPU()[i].assignedGPU = assignedGPU;
            theSim.sCPU()[i].runningOnCPU = forceCPU;
            theSim.sCPU()[i].useOpenMP = useOpenMP;
            std::lock_guard dataLock(theSim.mutexes.at(i));
            error = normalFunction(&theSim.sCPU()[i]);
            if (error) break;
        }
    }

    if (error) {
        theGui.console.tPrint("Encountered error {} in secondary queue.\n", error);
        return;
    }

    return;
}

void mainSimThread(int pulldownSelection, int secondPulldownSelection, bool use64bitFloatingPoint) {
    int error = 0;
    theSim.base().cancellationCalled = false;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

    readParametersFromInterface();
    theSim.configure();
    theGui.requestSliderUpdate();


    std::vector<simulationParameterSet> counterVector = theSim.getParameterVector();
    totalSteps = 0;
    progressCounter = 0;
    try {
        for (int64_t j = 0; j < theSim.base().Nsims * theSim.base().Nsims2; j++) {
            if (theSim.base().isInSequence) {
                counterVector[j].progressCounter = &totalSteps;
                counterVector[j].runType = runTypes::counter;
                solveNonlinearWaveEquationSequenceCounter(&counterVector[j]);
            }
            else {
                counterVector[j].progressCounter = &totalSteps;
                solveNonlinearWaveEquationCounter(&counterVector[j]);
            }
        }
    }
    catch (std::exception const& e) {
        std::string errorString = e.what();
        std::erase(errorString, '<');
        std::erase(errorString, '>');
        std::erase(errorString, '&');
        std::erase(errorString, ';');
        std::erase(errorString, '{');
        std::erase(errorString, '}');
        theGui.console.tPrint(
            "<span color=\"#FF88FF\">Simulation failed with exception:\n{}</span>\n",
            errorString);
        return;
    }

    
    
    theSim.base().isRunning = true;
    auto sequenceFunction = use64bitFloatingPoint ? 
        &solveNonlinearWaveEquationSequenceCPU : &solveNonlinearWaveEquationSequenceCPUFP32;
    auto normalFunction = use64bitFloatingPoint ? 
        &solveNonlinearWaveEquationCPU : &solveNonlinearWaveEquationCPUFP32;
    int assignedGPU = 0;
    bool forceCPU = false;
    bool useOpenMP = false;
#ifdef CPUONLY
    useOpenMP = true;
#endif
    [[maybe_unused]]int SYCLitems = 0;
    #if !defined(CPUONLY)
    if (theSim.base().syclGPUCount == 0) {
        SYCLitems = (int)theSim.base().SYCLavailable;
    }
    else {
        SYCLitems = 3;
    }
    #endif
    #if !defined(CPUONLY) && !defined(NOCUDA)
    if (pulldownSelection < theSim.base().cudaGPUCount) {
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
    #endif
    #ifndef NOSYCL
    if (pulldownSelection == theSim.base().cudaGPUCount && theSim.base().SYCLavailable) {
        if (theGui.firstSYCLsimulation) theGui.console.tPrint(
            "Note: the first time you run SYCL, it will\n"
            "take some time to compile kernels for your\n"
            "device. Subsequent runs will be faster.\n");
        theGui.firstSYCLsimulation = false;
        if (use64bitFloatingPoint) {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
            normalFunction = &solveNonlinearWaveEquationSYCL;
        }
        else {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
            normalFunction = &solveNonlinearWaveEquationSYCLFP32;
        }

    }
    else if (pulldownSelection == theSim.base().cudaGPUCount + 1 && SYCLitems > 1) {
        forceCPU = 1;
        if (theGui.firstSYCLsimulation) theGui.console.tPrint(
            "Note: the first time you run SYCL, it will\n"
            "take some time to compile kernels for your\n"
            "device. Subsequent runs will be faster.\n");
        theGui.firstSYCLsimulation = false;
        if (use64bitFloatingPoint) {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
            normalFunction = &solveNonlinearWaveEquationSYCL;
        }
        else {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
            normalFunction = &solveNonlinearWaveEquationSYCLFP32;
        }
    }
    else if (pulldownSelection == theSim.base().cudaGPUCount + 2 && SYCLitems > 1) {
        assignedGPU = 1;
        if (theGui.firstSYCLsimulation) theGui.console.tPrint(
            "Note: the first time you run SYCL, it will\n"
            "take some time to compile kernels for your\n"
            "device. Subsequent runs will be faster.\n");
        theGui.firstSYCLsimulation = false;
        if (use64bitFloatingPoint) {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
            normalFunction = &solveNonlinearWaveEquationSYCL;
        }
        else {
            sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
            normalFunction = &solveNonlinearWaveEquationSYCLFP32;
        }
    }
    else if (pulldownSelection == (theSim.base().cudaGPUCount + SYCLitems + 1)){
        useOpenMP = true;
    }
    #endif

    std::thread secondQueueThread(secondaryQueue, 
        theSim.base().Nsims * theSim.base().Nsims2 - theSim.base().NsimsCPU, 
        secondPulldownSelection, pulldownSelection, use64bitFloatingPoint);

    for (int j = 0; j < (theSim.base().Nsims * theSim.base().Nsims2 - theSim.base().NsimsCPU); ++j) {
        theSim.sCPU()[j].runningOnCPU = forceCPU;
        theSim.sCPU()[j].assignedGPU = assignedGPU;
        theSim.sCPU()[j].useOpenMP = useOpenMP;
        std::lock_guard dataLock(theSim.mutexes.at(j));
        if (theSim.base().isInSequence) {
            try {
                error = sequenceFunction(&theSim.sCPU()[j]);
            }
            catch (std::exception const& e) {
                std::string errorString=e.what();
                std::erase(errorString,'<');
                std::erase(errorString,'>');
                std::erase(errorString,'&');
                std::erase(errorString,';');
                std::erase(errorString,'{');
                std::erase(errorString,'}');
                theGui.console.tPrint(
                    "<span color=\"#FF88FF\">Simulation failed with exception:\n{}</span>\n", 
                    errorString);
            }
            if (theSim.sCPU()[j].memoryError != 0) {
                if (theSim.sCPU()[j].memoryError == -1) {
                    theGui.console.tPrint((
                        "<span color=\"#FF88FF\">Not enough free GPU memory, sorry.</span>\n"), 
                        theSim.sCPU()[j].memoryError);
                }
                else {
                    theGui.console.tPrint((
                        "<span color=\"#FF88FF\">Warning: device memory error ({}).</span>\n"), 
                        theSim.sCPU()[j].memoryError);
                }
            }
            if (error) break;
            theGui.requestSliderMove(j);
            independentPlotQueue();
        }
        else {
            try {
                error = normalFunction(&theSim.sCPU()[j]);
            } catch (std::exception const& e) {
                std::string errorString=e.what();
                std::erase(errorString,'<');
                std::erase(errorString,'>');
                std::erase(errorString,'&');
                std::erase(errorString,';');
                std::erase(errorString,'{');
                std::erase(errorString,'}');
                theGui.console.tPrint(
                    "<span color=\"#FF88FF\">Simulation failed with exception:\n{}</span>\n", 
                    errorString);
            }
            
            if (theSim.sCPU()[j].memoryError != 0) {
                if (theSim.sCPU()[j].memoryError == -1) {
                    theGui.console.tPrint((
                        "<span color=\"#FF88FF\">Not enough free GPU memory, sorry.</span>\n"), 
                        theSim.sCPU()[j].memoryError);
                }
                else {
                    theGui.console.tPrint((
                        "<span color=\"#FF88FF\">Warning: device memory error ({}).</span>\r\n"), 
                        theSim.sCPU()[j].memoryError);
                }
            }
            if (error) break;
        }
        if (theSim.base().cancellationCalled) {
            theGui.console.tPrint((
                "<span color=\"#FF88FF\">"
                "Warning: series cancelled, stopping\n"
                "after {} simulations.</span>\n"), j + 1);
            break;
        }
        theGui.requestSliderMove(j);
        independentPlotQueue();
    }

    if (secondQueueThread.joinable()) secondQueueThread.join();
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    if (error == 13) {
        theGui.console.tPrint(
            "<span color=\"#FF88FF\">"
            "NaN detected in grid!\n"
            "Try using a larger spatial/temporal step\n"
            "or smaller propagation step.\n"
            "Simulation was cancelled.\n</span>");
    }
    else if (error == 15) {
        theGui.console.tPrint(
            "<span color=\"#FF88FF\">"
            "Sorry, that sequence mode has been \n"
            "replaced by the new one. Look in the \n"
            "documentation for more info. It is a lot\n"
            "easier to use now, and hopefully \n"
            "it won't take long to set it up. \n"
            "Sorry about that!\n</span>");
    }
    else if(!error){
        theGui.console.tPrint(
            "<span color=\"#88FFFF\">Finished after {:.4} s. </span>\n", 1e-6 *
            (double)(std::chrono::duration_cast<std::chrono::microseconds>
                (simulationTimerEnd - simulationTimerBegin).count()));
    }

    theSim.saveDataSet();
    theSim.base().isRunning = false;
}

void fittingThread(int pulldownSelection, bool use64bitFloatingPoint) {
    theSim.base().cancellationCalled = false;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

    readParametersFromInterface();


    theSim.configure();
    theGui.requestSliderUpdate();
    theSim.sCPU()->readFittingString();
    if (theSim.base().Nfitting == 0) {
        theGui.console.tPrint("Couldn't interpret fitting command.\n");
        return;
    }
    progressCounter = 0;
    theSim.base().progressCounter = &progressCounter;
    if (theSim.base().fittingMode == 3) {
        if (theSim.sCPU()->loadReferenceSpectrum()) {
            theGui.console.tPrint("Could not read reference file!\n");
            return;
        }
    }

    theGui.console.tPrint(
        "Fitting {} values in mode {} over {} iterations.\n"
        "Region of interest contains {} elements\n",
        theSim.base().Nfitting, 
        theSim.base().fittingMode, 
        theSim.base().fittingMaxIterations, 
        theSim.base().fittingROIsize);

    int assignedGPU = 0;
    bool forceCPU = false;
    bool useOpenMP = false;
#ifdef CPUONLY
    useOpenMP = true;
#endif
    [[maybe_unused]]int SYCLitems = 0;
    if (theSim.base().syclGPUCount == 0) {
        SYCLitems = (int)theSim.base().SYCLavailable;
    }
    else {
        SYCLitems = 3;
    }

    auto fittingFunction = &runDlibFittingCPU;
    if (!use64bitFloatingPoint) {
        fittingFunction = &runDlibFittingCPUFP32;
    }
#if !defined(CPUONLY) && !defined(NOCUDA)
    if (pulldownSelection < theSim.base().cudaGPUCount) {
        fittingFunction = &runDlibFitting;
        if (!use64bitFloatingPoint)fittingFunction = &runDlibFittingFP32;
        assignedGPU = pulldownSelection;
    }
    #ifndef NOSYCL
    else if (pulldownSelection == theSim.base().cudaGPUCount && theSim.base().SYCLavailable) {
        fittingFunction = &runDlibFittingSYCL;
        if (!use64bitFloatingPoint)fittingFunction = &runDlibFittingSYCLFP32;
    }
    else if (pulldownSelection == theSim.base().cudaGPUCount + 1 && SYCLitems > 1) {
        forceCPU = 1;
        fittingFunction = &runDlibFittingSYCL;
        if (!use64bitFloatingPoint)fittingFunction = &runDlibFittingSYCLFP32;
    }
    else if (pulldownSelection == theSim.base().cudaGPUCount + 2 && SYCLitems > 1) {
        assignedGPU = 1;
        fittingFunction = &runDlibFittingSYCL;
        if (!use64bitFloatingPoint)fittingFunction = &runDlibFittingSYCLFP32;
    }
    #endif
    else if (pulldownSelection == (theSim.base().cudaGPUCount + SYCLitems + 1)) {
        useOpenMP = true;
    }
#endif
    theSim.base().isRunning = true;
    theSim.base().runningOnCPU = forceCPU;
    theSim.base().assignedGPU = assignedGPU;
    theSim.base().useOpenMP = useOpenMP;

    std::unique_lock lock(theSim.mutexes.at(0));
    fittingFunction(theSim.sCPU());
    lock.unlock();
    theSim.base().plotSim = 0;

    theGui.requestPlotUpdate();
    theGui.requestInterfaceValuesUpdate();
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    theGui.console.tPrint(("Finished fitting after {:.4} s.\n"), 1e-6 *
        (double)(std::chrono::duration_cast<std::chrono::microseconds>
            (simulationTimerEnd - simulationTimerBegin).count()));
    theGui.console.tPrint(
        "Fitting result:\n"
        "(index, value)\n");
    for (int i = 0; i < theSim.base().Nfitting; ++i) {
        theGui.console.tPrint("{},  {}\n", i, theSim.base().fittingResult[i]);
    }
    theSim.saveDataSet();
    theSim.base().isRunning = false;
}

int main(int argc, char **argv){
    GtkApplication* app = gtk_application_new("io.github.NickKarpowicz.LightwaveExplorer", static_cast<GApplicationFlags>(0));
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    return g_application_run(G_APPLICATION(app), argc, argv);
}