#include "LightwaveExplorerFrontendQT.h"

class CairoWidget : public QWidget {
    Q_OBJECT
    LWEGui& theGui;
    CairoFunction theFunction;
public:
    CairoWidget(LWEGui& sim, CairoFunction fun, QWidget* parent = nullptr) : 
    QWidget(parent),
    theGui(sim), 
    theFunction(fun) {}
    
protected:
    void paintEvent(QPaintEvent* event) override {
        qreal ratio = devicePixelRatioF();
        QImage image(ratio*size(), QImage::Format_ARGB32);
        image.setDevicePixelRatio(ratio);
        cairo_surface_t* surface = cairo_image_surface_create_for_data(
            image.bits(),
            CAIRO_FORMAT_ARGB32,
            image.width(),
            image.height(),
            image.bytesPerLine());
        cairo_t* cr = cairo_create(surface);
        theFunction(cr,
            image.width(),
            image.height(),
            theGui);
        cairo_destroy(cr);
        cairo_surface_destroy(surface);
        QPainter painter(this);
        painter.drawImage(0, 0, image);
    }

public slots:
    void queueUpdate() {
        update();
    }
};

class GuiMessenger : public QObject {
    Q_OBJECT

public: 
    std::mutex m;
public slots:
    void passString(std::string s){
        std::unique_lock lock(m);
        emit sendText(QString::fromStdString(s));
    }

    void passDrawRequest(){
        emit requestUpdate();
    }

    void passSyncValues(){
        emit requestSyncValues();
    }

    void passSliderUpdate(){
        emit requestSliderUpdate();
    }

    void passSliderPosition(int value){
        emit moveSlider(value);
    }
    void passProgressRange(int min, int max){
        emit sendProgressRange(min, max);
    }
    void passProgressValue(int value){
        emit sendProgressValue(value);
    }

signals:
    void sendText(const QString &text);
    void requestUpdate();
    void requestSyncValues();
    void requestSliderUpdate();
    void moveSlider(int value);
    void sendProgressRange(int min, int max);
    void sendProgressValue(int value);
};

class LWEGui : public QObject {
    Q_OBJECT
    bool isInputRegionHidden=false;
    QThread* messengerThread;
    void populateDatabasePulldown(){
        std::string materialString;
        pulldowns["material"]->clear();
        for (int i = 0; i < theDatabase.db.size(); ++i) {
            materialString = Sformat(
                "{:2}: {}", i, theDatabase.db[i].crystalName);
            pulldowns["material"]->addItem(materialString.c_str());
        }
    }
public:    
//Main data structures:
// theSim contains all of the parameters of the current simulation including grid arrays
// theDatabase is the database of crystal properties
    simulationBatch theSim;
    crystalDatabase theDatabase;

//Extra data structures
    loadedInputData pulse1LoadedData;
    loadedInputData pulse2LoadedData;
    loadedInputData fittingLoadedData;
//Counter atomics
    std::atomic_uint32_t progressCounter{};
    std::atomic_uint32_t totalSteps{};

    //interface elements
    std::unordered_map<std::string, QPushButton*> buttons;
    std::unordered_map<std::string, QLineEdit*> textBoxes;
    std::unordered_map<std::string, QComboBox*> pulldowns;
    std::unordered_map<std::string, QLabel*> labels;
    std::unordered_map<std::string, QCheckBox*> checkboxes;
    std::unordered_map<std::string, CairoWidget*> plots;
    QTextEdit* sequence;
    QTextEdit* console;
    QTextEdit* fitting;
    QProgressBar* progress;
    QWidget *inputRegion;
    QSlider* slider;
    GuiMessenger* messenger;
    bool isMakingSVG = false;
    std::array<std::string,4> SVGStrings;
    template<typename... Args> void cPrint(std::string_view format, Args&&... args) {
        std::string s = Svformat(format, Smake_format_args(args...));
        console->append(s.c_str());
    }
    template<typename... Args> void sPrint(std::string_view format, Args&&... args) {
        std::string s = Svformat(format, Smake_format_args(args...));
        sequence->append(s.c_str());
    }

    void readParametersFromInterface(simulationBatch& sim) {
        auto setToDoubleMultiplier = [](QLineEdit* box, const double multiplier, double& value){;
            value = (box->text()).toDouble()/multiplier;
        };
        auto setToDouble = [](QLineEdit* box, double& value){;
            value = (box->text()).toDouble();
        };
        auto setToInt = [](QLineEdit* box, int& value){
            value = (box->text()).toInt();
        };
        auto setToInt64 = [](QLineEdit* box, int64_t& value){
            value = static_cast<int64_t>((box->text()).toInt());
        };

        setToDouble(textBoxes["Energy1"],sim.base().pulse1.energy);
        setToDoubleMultiplier(textBoxes["Frequency1"],1e-12,sim.base().pulse1.frequency);
        setToDoubleMultiplier(textBoxes["Bandwidth1"],1e-12,sim.base().pulse1.bandwidth);
        setToInt(textBoxes["SGOrder1"],sim.base().pulse1.sgOrder);
        setToDouble(textBoxes["CEP1"],sim.base().pulse1.cep);
        setToDoubleMultiplier(textBoxes["Delay1"],1e15,sim.base().pulse1.delay);
        setToDoubleMultiplier(textBoxes["GDD1"],1e30,sim.base().pulse1.gdd);
        setToDoubleMultiplier(textBoxes["TOD1"],1e45,sim.base().pulse1.tod);
        setToInt(textBoxes["Material1"],sim.base().pulse1.phaseMaterial);
        setToDoubleMultiplier(textBoxes["Thickness1"],1e6,sim.base().pulse1.phaseMaterialThickness);
        setToDoubleMultiplier(textBoxes["Beamwaist1"],1e6,sim.base().pulse1.beamwaist);
        setToDoubleMultiplier(textBoxes["xOffset1"],1e6,sim.base().pulse1.x0);
        setToDoubleMultiplier(textBoxes["zOffset1"],1e6,sim.base().pulse1.z0);
        setToDoubleMultiplier(textBoxes["NCAngle1"],rad2Deg<double>(),sim.base().pulse1.beamAngle);
        setToDoubleMultiplier(textBoxes["Polarization1"],rad2Deg<double>(),sim.base().pulse1.polarizationAngle);
        setToDouble(textBoxes["Circularity1"],sim.base().pulse1.circularity);

        setToDouble(textBoxes["Energy2"],sim.base().pulse2.energy);
        setToDoubleMultiplier(textBoxes["Frequency2"],1e-12,sim.base().pulse2.frequency);
        setToDoubleMultiplier(textBoxes["Bandwidth2"],1e-12,sim.base().pulse2.bandwidth);
        setToInt(textBoxes["SGOrder2"],sim.base().pulse2.sgOrder);
        setToDouble(textBoxes["CEP2"],sim.base().pulse2.cep);
        setToDoubleMultiplier(textBoxes["Delay2"],1e15,sim.base().pulse2.delay);
        setToDoubleMultiplier(textBoxes["GDD2"],1e30,sim.base().pulse2.gdd);
        setToDoubleMultiplier(textBoxes["TOD2"],1e45,sim.base().pulse2.tod);
        setToInt(textBoxes["Material2"],sim.base().pulse2.phaseMaterial);
        setToDoubleMultiplier(textBoxes["Thickness2"],1e6,sim.base().pulse2.phaseMaterialThickness);
        setToDoubleMultiplier(textBoxes["Beamwaist2"],1e6,sim.base().pulse2.beamwaist);
        setToDoubleMultiplier(textBoxes["xOffset2"],1e6,sim.base().pulse2.x0);
        setToDoubleMultiplier(textBoxes["zOffset2"],1e6,sim.base().pulse2.z0);
        setToDoubleMultiplier(textBoxes["NCAngle2"],rad2Deg<double>(),sim.base().pulse2.beamAngle);
        setToDoubleMultiplier(textBoxes["Polarization2"],rad2Deg<double>(),sim.base().pulse2.polarizationAngle);
        setToDouble(textBoxes["Circularity2"],sim.base().pulse2.circularity);
        sim.base().pulse1FileType = pulldowns["pulse1type"]->currentIndex();
        sim.base().pulse2FileType = pulldowns["pulse2type"]->currentIndex();
        sim.base().fittingMode = pulldowns["fit"]->currentIndex();
        sim.base().materialIndex = pulldowns["material"]->currentIndex();
        setToDoubleMultiplier(textBoxes["CrystalTheta"],rad2Deg<double>(),sim.base().crystalTheta);
        setToDoubleMultiplier(textBoxes["CrystalPhi"],rad2Deg<double>(),sim.base().crystalPhi);
        setToDouble(textBoxes["NLAbsorption"],sim.base().nonlinearAbsorptionStrength);
        setToDouble(textBoxes["CrystalBandgap"],sim.base().bandGapElectronVolts);
        setToDoubleMultiplier(textBoxes["DrudeGamma"],1e-12,sim.base().drudeGamma);
        setToDouble(textBoxes["effectiveMass"],sim.base().effectiveMass);
        setToDoubleMultiplier(textBoxes["XSize"],1e6,sim.base().spatialWidth);
        setToDoubleMultiplier(textBoxes["dx"],1e6,sim.base().rStep);
        setToDoubleMultiplier(textBoxes["timeSpan"],1e15,sim.base().timeSpan);
        setToDoubleMultiplier(textBoxes["dt"],1e15,sim.base().tStep);
        setToDoubleMultiplier(textBoxes["ZSize"],1e6,sim.base().crystalThickness);
        setToDoubleMultiplier(textBoxes["dz"],1e9,sim.base().propagationStep);

        sim.base().symmetryType = pulldowns["propagator"]->currentIndex();
        sim.base().batchIndex = pulldowns["batch1"]->currentIndex();
        sim.base().batchIndex2 = pulldowns["batch2"]->currentIndex();
        setToDouble(textBoxes["batch1end"],sim.base().batchDestination);
        setToDouble(textBoxes["batch2end"],sim.base().batchDestination2);
        setToInt64(textBoxes["batch1steps"],sim.base().Nsims);
        setToInt64(textBoxes["batch2steps"],sim.base().Nsims2);
        setToInt64(textBoxes["offload"],sim.base().NsimsCPU);
        
        sim.base().isInSequence = false;
        sim.base().sequenceString = sequence->toPlainText().toStdString();
        stripWhiteSpace(sim.base().sequenceString);
        if (sim.base().sequenceString[0] != 'N' 
            && sim.base().sequenceString.length() > 5) 
            sim.base().isInSequence = true;

        sim.base().isInFittingMode = false;
        sim.base().fittingString = fitting->toPlainText().toStdString();
        stripLineBreaks(sim.base().fittingString);

        sim.base().pulse1LoadedData = pulse1LoadedData;
        sim.base().pulse2LoadedData = pulse2LoadedData;

        //theGui.filePaths[3].copyBuffer(sim.base().outputBasePath);
        // sim.base().outputBasePath = theGui.currentPath;
        // stripLineBreaks(sim.base().outputBasePath);
        // if ((sim.base().outputBasePath.length() > 4 && sim.base().outputBasePath.substr(sim.base().outputBasePath.length() - 4) == ".txt")
        //     || (sim.base().outputBasePath.length() > 4 && sim.base().outputBasePath.substr(sim.base().outputBasePath.length() - 4) == ".zip")) {
        //     sim.base().outputBasePath = sim.base().outputBasePath.substr(0, sim.base().outputBasePath.length() - 4);
        // }

        sim.base().fittingLoadedData = fittingLoadedData;
        
        //derived parameters and cleanup:
        sim.base().sellmeierType = 0;
        sim.base().axesNumber = 0;
        sim.base().Ntime = (int64_t)(minGridDimension * round(sim.base().timeSpan 
            / (minGridDimension * sim.base().tStep)));
        if (sim.base().symmetryType == 2) {
            sim.base().is3D = true;
            sim.base().isFDTD = false;
            sim.base().spatialWidth = sim.base().rStep 
                * (minGridDimension * round(sim.base().spatialWidth 
                    / (sim.base().rStep * minGridDimension)));
            sim.base().Nspace = (int64_t)round(sim.base().spatialWidth 
                / sim.base().rStep);
            if (sim.base().spatialHeight > 0) {
                sim.base().spatialHeight = sim.base().rStep 
                    * (minGridDimension * round(sim.base().spatialHeight 
                        / (sim.base().rStep * minGridDimension)));
            }
            else {
                sim.base().spatialHeight = sim.base().spatialWidth;
            }
            sim.base().Nspace2 = (int64_t)round(sim.base().spatialHeight 
                / sim.base().rStep);
        }
        else if (sim.base().symmetryType == 3) {
            sim.base().is3D = false;
            sim.base().isFDTD = true;
            sim.base().Nspace2 = 1;
            sim.base().spatialHeight = 0;
            sim.base().spatialWidth = sim.base().rStep 
                * (minGridDimension * round(sim.base().spatialWidth 
                    / (sim.base().rStep * minGridDimension)));
            sim.base().Nspace = (int64_t)round(sim.base().spatialWidth 
                / sim.base().rStep);
        }
        else if (sim.base().symmetryType == 4) {
            sim.base().is3D = true;
            sim.base().isFDTD = true;
            sim.base().spatialWidth = sim.base().rStep 
                * (minGridDimension * round(sim.base().spatialWidth 
                    / (sim.base().rStep * minGridDimension)));
            sim.base().Nspace = (int64_t)round(sim.base().spatialWidth 
                / sim.base().rStep);
            if (sim.base().spatialHeight > 0) {
                sim.base().spatialHeight = sim.base().rStep 
                    * (minGridDimension * round(sim.base().spatialHeight 
                        / (sim.base().rStep * minGridDimension)));
            }
            else {
                sim.base().spatialHeight = sim.base().spatialWidth;
            }
            sim.base().Nspace2 = (int64_t)round(sim.base().spatialHeight 
                / sim.base().rStep);
        }
        else {
            sim.base().is3D = false;
            sim.base().isFDTD = false;
            sim.base().Nspace2 = 1;
            sim.base().spatialHeight = 0;
            sim.base().spatialWidth = sim.base().rStep 
                * (minGridDimension * round(sim.base().spatialWidth 
                    / (sim.base().rStep * minGridDimension)));
            sim.base().Nspace = (int64_t)round(sim.base().spatialWidth 
                / sim.base().rStep);
        }

        sim.base().Nfreq = sim.base().Ntime / 2 + 1;
        sim.base().NgridC = sim.base().Nfreq * sim.base().Nspace * sim.base().Nspace2;
        sim.base().Ngrid = sim.base().Ntime * sim.base().Nspace * sim.base().Nspace2;
        sim.base().kStep = twoPi<double>() / (sim.base().Nspace * sim.base().rStep);
        sim.base().fStep = 1.0 / (sim.base().Ntime * sim.base().tStep);
        sim.base().Npropagation = (int64_t)round(sim.base().crystalThickness 
            / sim.base().propagationStep);

        sim.base().isCylindric = sim.base().symmetryType == 1;
        if (sim.base().isCylindric) {
            sim.base().pulse1.x0 = 0;
            sim.base().pulse2.x0 = 0;
            sim.base().pulse1.beamAngle = 0;
            sim.base().pulse2.beamAngle = 0;
        }

        if (sim.base().batchIndex == 0 || sim.base().Nsims < 1) {
            sim.base().Nsims = 1;
        }
        if (sim.base().batchIndex2 == 0 || sim.base().Nsims2 < 1) {
            sim.base().Nsims2 = 1;
        }
        sim.base().NsimsCPU = minN(sim.base().NsimsCPU, 
            sim.base().Nsims * sim.base().Nsims2);

        sim.base().field1IsAllocated = false;
        sim.base().field2IsAllocated = false;

        //crystal from database (database must be loaded!)
        sim.base().crystalDatabase = theDatabase.db.data();
        sim.base().chi2Tensor = 
            theDatabase.db[sim.base().materialIndex].d.data();
        sim.base().chi3Tensor = 
            theDatabase.db[sim.base().materialIndex].chi3.data();
        sim.base().nonlinearSwitches = 
            theDatabase.db[sim.base().materialIndex].nonlinearSwitches.data();
        sim.base().absorptionParameters = 
            theDatabase.db[sim.base().materialIndex].absorptionParameters.data();
        sim.base().sellmeierCoefficients = 
            theDatabase.db[sim.base().materialIndex].sellmeierCoefficients.data();
        sim.base().sellmeierType = 
            theDatabase.db[sim.base().materialIndex].sellmeierType;
        sim.base().axesNumber = 
            theDatabase.db[sim.base().materialIndex].axisType;
        sim.base().progressCounter = &progressCounter;

        sim.base().runType = runTypes::normal;
        sim.base().isFollowerInSequence = false;
        sim.base().crystalDatabase = theDatabase.db.data();
    }

    void setInterfaceValuesToActiveValues(simulationBatch& sim){
        auto setToDouble = [](QLineEdit* box, const double value){
            QString s(Sformat(std::string_view("{:g}"), value).c_str());
            box->setText(s);
        };
        auto setToInt = [](QLineEdit* box, const int value){
            QString s(Sformat(std::string_view("{}"), value).c_str());
            box->setText(s);
        };
        setToDouble(textBoxes["Energy1"],sim.base().pulse1.energy);
        setToDouble(textBoxes["Frequency1"],1e-12*sim.base().pulse1.frequency);
        setToDouble(textBoxes["Bandwidth1"],1e-12*sim.base().pulse1.bandwidth);
        setToDouble(textBoxes["SGOrder1"],sim.base().pulse1.sgOrder);
        setToDouble(textBoxes["CEP1"],sim.base().pulse1.cep);
        setToDouble(textBoxes["Delay1"],1e15*sim.base().pulse1.delay);
        setToDouble(textBoxes["GDD1"],1e30*sim.base().pulse1.gdd);
        setToDouble(textBoxes["TOD1"],1e45*sim.base().pulse1.tod);
        setToInt(textBoxes["Material1"],sim.base().pulse1.phaseMaterial);
        setToDouble(textBoxes["Thickness1"],1e6*sim.base().pulse1.phaseMaterialThickness);
        setToDouble(textBoxes["Beamwaist1"],1e6*sim.base().pulse1.beamwaist);
        setToDouble(textBoxes["xOffset1"],1e6*sim.base().pulse1.x0);
        setToDouble(textBoxes["zOffset1"],1e6*sim.base().pulse1.z0);
        setToDouble(textBoxes["NCAngle1"],rad2Deg<double>()*sim.base().pulse1.beamAngle);
        setToDouble(textBoxes["Polarization1"],rad2Deg<double>()*sim.base().pulse1.polarizationAngle);
        setToDouble(textBoxes["Circularity1"],sim.base().pulse1.circularity);

        setToDouble(textBoxes["Energy2"],sim.base().pulse2.energy);
        setToDouble(textBoxes["Frequency2"],1e-12*sim.base().pulse2.frequency);
        setToDouble(textBoxes["Bandwidth2"],1e-12*sim.base().pulse2.bandwidth);
        setToDouble(textBoxes["SGOrder2"],sim.base().pulse2.sgOrder);
        setToDouble(textBoxes["CEP2"],sim.base().pulse2.cep);
        setToDouble(textBoxes["Delay2"],1e15*sim.base().pulse2.delay);
        setToDouble(textBoxes["GDD2"],1e30*sim.base().pulse2.gdd);
        setToDouble(textBoxes["TOD2"],1e45*sim.base().pulse2.tod);
        setToInt(textBoxes["Material2"],sim.base().pulse2.phaseMaterial);
        setToDouble(textBoxes["Thickness2"],1e6*sim.base().pulse2.phaseMaterialThickness);
        setToDouble(textBoxes["Beamwaist2"],1e6*sim.base().pulse2.beamwaist);
        setToDouble(textBoxes["xOffset2"],1e6*sim.base().pulse2.x0);
        setToDouble(textBoxes["zOffset2"],1e6*sim.base().pulse2.z0);
        setToDouble(textBoxes["NCAngle2"],rad2Deg<double>()*sim.base().pulse2.beamAngle);
        setToDouble(textBoxes["Polarization2"],rad2Deg<double>()*sim.base().pulse2.polarizationAngle);
        setToDouble(textBoxes["Circularity2"],sim.base().pulse2.circularity);

        pulldowns["material"]->setCurrentIndex(sim.base().materialIndex);
        pulldowns["batch1"]->setCurrentIndex(sim.base().batchIndex);
        pulldowns["batch2"]->setCurrentIndex(sim.base().batchIndex2);
        pulldowns["propagator"]->setCurrentIndex(sim.base().symmetryType);

        setToDouble(textBoxes["CrystalTheta"],rad2Deg<double>() * sim.base().crystalTheta);
        setToDouble(textBoxes["CrystalPhi"],rad2Deg<double>() *sim.base().crystalPhi);
        setToDouble(textBoxes["NLAbsorption"],sim.base().nonlinearAbsorptionStrength);
        setToDouble(textBoxes["CrystalBandgap"],sim.base().bandGapElectronVolts);
        setToDouble(textBoxes["DrudeGamma"],1e-12*sim.base().drudeGamma);
        setToDouble(textBoxes["effectiveMass"],sim.base().effectiveMass);
        setToDouble(textBoxes["XSize"],1e6*sim.base().spatialWidth);
        setToDouble(textBoxes["dx"],1e6*sim.base().rStep);
        setToDouble(textBoxes["timeSpan"],1e15*sim.base().timeSpan);
        setToDouble(textBoxes["dt"],1e15*sim.base().tStep);
        setToDouble(textBoxes["ZSize"],1e6*sim.base().crystalThickness);
        setToDouble(textBoxes["dz"],1e9*sim.base().propagationStep);

        setToDouble(textBoxes["batch1end"],sim.base().batchDestination);
        setToDouble(textBoxes["batch2end"],sim.base().batchDestination2);
        setToInt(textBoxes["batch1steps"],sim.base().Nsims);
        setToInt(textBoxes["batch2steps"],sim.base().Nsims2);

        std::string formattedFit=sim.base().fittingString;
        insertAfterCharacter(formattedFit,';',std::string("\n"));
        fitting->setText(QString::fromStdString(formattedFit));

    }
    LWEGui(){
        const int textBoxWidth = 80;
        const int labelWidth = 160;
        const int textBoxHeight = 26;
        const int rowHeight = textBoxHeight+2;
        const int rowWidth = labelWidth + 2*textBoxWidth + 10;
        const int miniButtonWidth = 30;
        const int mainButtonWidth = rowWidth/4;
        const int mainButtonHeight = textBoxHeight;
        const int pulldownContainerWidth = labelWidth+4;

        //Divide the main window into a large expanding upper panel and a control strip at the bottom
        auto squeezeMargins = [&](QBoxLayout* layout){
            layout->setSpacing(0);
            layout->setContentsMargins(1,1,1,1);
        };
        QWidget *windowBody = new QWidget;
        windowBody->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        QVBoxLayout* windowBodyLayOut = new QVBoxLayout(windowBody);
        squeezeMargins(windowBodyLayOut);
        QWidget *mainArea = new QWidget;
        mainArea->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        QWidget *controlStrips = new QWidget;
        controlStrips->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        windowBodyLayOut->addWidget(mainArea);
        windowBodyLayOut->addWidget(controlStrips);

        //Divide the large window area into an input area and a plot area
        inputRegion = new QWidget;
        QWidget *plotRegion = new QWidget;

        QHBoxLayout *mainAreaLayout = new QHBoxLayout(mainArea);
        squeezeMargins(mainAreaLayout);
        inputRegion->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
        plotRegion->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        mainAreaLayout->addWidget(inputRegion);
        mainAreaLayout->addWidget(plotRegion);

        QGridLayout* plotRegionLayout = new QGridLayout(plotRegion);
        plotRegionLayout->setContentsMargins(1,1,1,1);
        plots["TimeImage1"] = new CairoWidget(*this, drawTimeImage1);
        plots["TimeImage1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["TimeImage1"],0,0);
        plots["TimeImage2"] = new CairoWidget(*this, drawTimeImage2);
        plots["TimeImage2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["TimeImage2"],1,0);
        plots["FreqImage1"] = new CairoWidget(*this, drawFourierImage1);
        plots["FreqImage1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["FreqImage1"],0,1);
        plots["FreqImage2"] = new CairoWidget(*this, drawFourierImage2);
        plots["FreqImage2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["FreqImage2"],1,1);

        plots["TimePlot1"] = new CairoWidget(*this, drawField1Plot);
        plots["TimePlot1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["TimePlot1"],2,0);
        plots["TimePlot2"] = new CairoWidget(*this, drawField2Plot);
        plots["TimePlot2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["TimePlot2"],3,0);
        plots["FreqPlot1"] = new CairoWidget(*this, drawSpectrum1Plot);
        plots["FreqPlot1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["FreqPlot1"],2,1);
        plots["FreqPlot2"] = new CairoWidget(*this, drawSpectrum2Plot);
        plots["FreqPlot2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["FreqPlot2"],3,1);

        //Divide the input area into the input grid and the sequence box
        QWidget *inputGrid = new QWidget;
        QWidget *sequenceBox = new QWidget;
        QVBoxLayout *inputRegionLayout = new QVBoxLayout(inputRegion);
        squeezeMargins(inputRegionLayout);
        inputGrid->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sequenceBox->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
        inputRegionLayout->addWidget(inputGrid);
        inputRegionLayout->addWidget(sequenceBox);

        //Divide the input grid into columns
        QWidget *entryColumn1 = new QWidget;
        QWidget *entryColumn2 = new QWidget;
        QHBoxLayout *inputGridLayout = new QHBoxLayout(inputGrid);
        squeezeMargins(inputGridLayout);
        inputGridLayout->addWidget(entryColumn1);
        inputGridLayout->addWidget(entryColumn2);
        QVBoxLayout *entryColumn1Layout = new QVBoxLayout(entryColumn1);
        squeezeMargins(entryColumn1Layout);
        QVBoxLayout *entryColumn2Layout = new QVBoxLayout(entryColumn2);
        squeezeMargins(entryColumn2Layout);

        //Divide the control strip into simulation and plot areas
        QHBoxLayout* controlStripsLayout = new QHBoxLayout(controlStrips);
        squeezeMargins(controlStripsLayout);
        QWidget* simulationControlStrip = new QWidget;
        QWidget* plotControlStrip = new QWidget;
        simulationControlStrip->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        plotControlStrip->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        controlStripsLayout->addWidget(simulationControlStrip);
        controlStripsLayout->addWidget(plotControlStrip);
        simulationControlStrip->setFixedSize(2*rowWidth,rowHeight);

        auto getRowBoxLayout = [&](QBoxLayout* location){
            QWidget *rowBox = new QWidget;
            rowBox->setFixedSize(rowWidth,rowHeight);
            location->addWidget(rowBox);
            QHBoxLayout *rowLayout = new QHBoxLayout(rowBox);
            squeezeMargins(rowLayout);
            return rowLayout;
        };

        auto addTextBoxRow = [&](
            const QString& label,
            const std::string& entry1, 
            const std::string& entry2,  
            QVBoxLayout* location){
                QHBoxLayout* rowLayout = getRowBoxLayout(location);
                labels[entry1] = new QLabel;
                textBoxes[entry1] = new QLineEdit;
                textBoxes[entry2] = new QLineEdit;
                labels[entry1]->setText(label);
                labels[entry1]->setFixedSize(labelWidth,textBoxHeight);
                textBoxes[entry1]->setFixedSize(textBoxWidth,textBoxHeight);
                textBoxes[entry2]->setFixedSize(textBoxWidth,textBoxHeight);
                rowLayout->addWidget(labels[entry1]);
                rowLayout->addWidget(textBoxes[entry1]);
                rowLayout->addWidget(textBoxes[entry2]); 
        };

        auto addPulldownInContainer = [&](int width, int height, QBoxLayout* location, const std::string& entry){
            QWidget* container = new QWidget;
            QHBoxLayout* fitContainerLayout = new QHBoxLayout(container);
            fitContainerLayout->setContentsMargins(0,0,0,0);
            pulldowns[entry] = new QComboBox;
            container->setFixedSize(width,height);
            location->addWidget(container);
            fitContainerLayout->addWidget(pulldowns[entry]);
        };


        //First column
        addTextBoxRow("Pulse energy (J)", "Energy1", "Energy2", entryColumn1Layout);
        addTextBoxRow("Frequency (THz)", "Frequency1", "Frequency2", entryColumn1Layout);
        addTextBoxRow("Bandwidth (THz)", "Bandwidth1", "Bandwidth2", entryColumn1Layout);
        addTextBoxRow("SG order", "SGOrder1", "SGOrder2", entryColumn1Layout);
        addTextBoxRow("CEP/\xcf\x80", "CEP1", "CEP2", entryColumn1Layout);
        addTextBoxRow("Delay (fs)", "Delay1", "Delay2", entryColumn1Layout);
        addTextBoxRow("GDD (fs\xc2\xb2)", "GDD1", "GDD2", entryColumn1Layout);
        addTextBoxRow("TOD (fs\xc2\xb3)", "TOD1", "TOD2", entryColumn1Layout);
        addTextBoxRow("Phase material", "Material1", "Material2", entryColumn1Layout);
        addTextBoxRow("Thickness (\xce\xbcm)", "Thickness1", "Thickness2", entryColumn1Layout);
        addTextBoxRow("Bewamwaist (\xce\xbcm)", "Beamwaist1", "Beamwaist2", entryColumn1Layout);
        addTextBoxRow("x offset (\xce\xbcm)", "xOffset1", "xOffset2", entryColumn1Layout);
        addTextBoxRow("z offset (\xce\xbcm)", "zOffset1", "zOffset2", entryColumn1Layout);
        addTextBoxRow("NC angle (deg)", "NCAngle1", "NCAngle2", entryColumn1Layout);
        addTextBoxRow("Polarization (deg)", "Polarization1", "Polarization2", entryColumn1Layout);
        addTextBoxRow("Circularity", "Circularity1", "Circularity2", entryColumn1Layout);

        QHBoxLayout* pulseTypeRow = getRowBoxLayout(entryColumn1Layout);
        labels["source"] = new QLabel;
        labels["source"]->setText("Source");
        labels["source"]->setFixedWidth(labelWidth);
        pulseTypeRow->addWidget(labels["source"]);
        addPulldownInContainer(textBoxWidth,mainButtonHeight,pulseTypeRow,"pulse1type");
        pulldowns["pulse1type"]->addItem("Synthetic");
        pulldowns["pulse1type"]->addItem("FROG");
        pulldowns["pulse1type"]->addItem("Wave");
        pulldowns["pulse1type"]->addItem("LWE");
        addPulldownInContainer(textBoxWidth,mainButtonHeight,pulseTypeRow,"pulse2type");
        pulldowns["pulse2type"]->addItem("Synthetic");
        pulldowns["pulse2type"]->addItem("FROG");
        pulldowns["pulse2type"]->addItem("Wave");
        pulldowns["pulse2type"]->addItem("LWE");

        QHBoxLayout* loadPulseRow = getRowBoxLayout(entryColumn1Layout);
        labels["sourceFile"] = new QLabel;
        labels["sourceFile"]->setText("Source file");
        labels["sourceFile"]->setFixedWidth(labelWidth);
        loadPulseRow->addWidget(labels["sourceFile"]);
        buttons["LoadPulse1"] = new QPushButton("\xf0\x9f\x93\x82");
        buttons["LoadPulse1"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        loadPulseRow->addWidget(buttons["LoadPulse1"]);
        buttons["LoadPulse2"] = new QPushButton("\xf0\x9f\x93\x82");
        buttons["LoadPulse2"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        loadPulseRow->addWidget(buttons["LoadPulse2"]);

        QHBoxLayout* slurmRow = getRowBoxLayout(entryColumn1Layout);
        labels["slurm"] = new QLabel;
        labels["slurm"]->setText("SLURM script");
        labels["slurm"]->setFixedWidth(labelWidth);
        slurmRow->addWidget(labels["slurm"]);
        addPulldownInContainer(pulldownContainerWidth-miniButtonWidth,mainButtonHeight,slurmRow,"slurm");
        pulldowns["slurm"]->addItem("Cobra 1xR5k");
        pulldowns["slurm"]->addItem("Cobra 2xR5k");
        pulldowns["slurm"]->addItem("Cobra 1xV100");
        pulldowns["slurm"]->addItem("Cobra 2xV100");
        pulldowns["slurm"]->addItem("Raven 1xA100");
        pulldowns["slurm"]->addItem("Raven 2xA100");
        pulldowns["slurm"]->addItem("Raven 4xA100");
        pulldowns["slurm"]->addItem("Raven NxA100");
        
        buttons["saveSLURM"] = new QPushButton("\xf0\x9f\x92\xbe");
        buttons["saveSLURM"]->setFixedSize(miniButtonWidth,mainButtonHeight);
        slurmRow->addWidget(buttons["saveSLURM"]);

        //Second column
        QHBoxLayout* materialRow = getRowBoxLayout(entryColumn2Layout);
        labels["material"] = new QLabel;
        labels["material"]->setText("Material");
        labels["material"]->setFixedWidth(labelWidth-miniButtonWidth);
        materialRow->addWidget(labels["material"]);
        buttons["loadMaterial"] = new QPushButton("\xf0\x9f\x93\x82");
        buttons["loadMaterial"]->setFixedWidth(miniButtonWidth);
        materialRow->addWidget(buttons["loadMaterial"]);
        addPulldownInContainer(pulldownContainerWidth,mainButtonHeight,materialRow,"material");
        populateDatabasePulldown();

        
        addTextBoxRow("Theta, phi (deg)", "CrystalTheta", "CrystalPhi", entryColumn2Layout);
        addTextBoxRow("NL absorption", "NLAbsorption", "CrystalBandgap", entryColumn2Layout);
        addTextBoxRow("Drude: gamma, m", "DrudeGamma", "effectiveMass", entryColumn2Layout);
        addTextBoxRow("Max x, dx (\xce\xbcm)", "XSize", "dx", entryColumn2Layout);
        addTextBoxRow("Time span, dt (fs)", "timeSpan", "dt", entryColumn2Layout);
        addTextBoxRow("Max z, dz (\xce\xbcm, nm)", "ZSize", "dz", entryColumn2Layout);

        QHBoxLayout* propagationRow = getRowBoxLayout(entryColumn2Layout);
        labels["propagator"] = new QLabel;
        labels["propagator"]->setText("Propagation");
        labels["propagator"]->setFixedWidth(labelWidth);
        propagationRow->addWidget(labels["propagator"]);
        addPulldownInContainer(pulldownContainerWidth,mainButtonHeight,propagationRow,"propagator");
        pulldowns["propagator"]->addItem(("2D Cartesian"));
        pulldowns["propagator"]->addItem(("3D radial symmmetry"));
        pulldowns["propagator"]->addItem(("3D"));
        pulldowns["propagator"]->addItem(("FDTD 2D"));
        pulldowns["propagator"]->addItem(("FDTD 3D"));

        QHBoxLayout* batch1Row = getRowBoxLayout(entryColumn2Layout);
        labels["batch1"] = new QLabel;
        labels["batch1"]->setText("Batch mode");
        labels["batch1"]->setFixedWidth(labelWidth);
        batch1Row->addWidget(labels["batch1"]);
        addPulldownInContainer(pulldownContainerWidth,mainButtonHeight,batch1Row,"batch1");

        QHBoxLayout* batch2Row = getRowBoxLayout(entryColumn2Layout);
        labels["batch2"] = new QLabel;
        labels["batch2"]->setText("Batch 2 mode");
        labels["batch2"]->setFixedWidth(labelWidth);
        batch2Row->addWidget(labels["batch2"]);
        addPulldownInContainer(pulldownContainerWidth,mainButtonHeight,batch2Row,"batch2");
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
            pulldowns["batch1"]->addItem(batchModeNames[i]);
            pulldowns["batch2"]->addItem(batchModeNames[i]);
        }

        addTextBoxRow("Batch end", "batch1end", "batch2end", entryColumn2Layout);
        addTextBoxRow("Batch steps", "batch1steps", "batch2steps", entryColumn2Layout);
        
        QHBoxLayout* loadRow = getRowBoxLayout(entryColumn2Layout);
        buttons["Load"] = new QPushButton("Load");
        buttons["Load"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        loadRow->addWidget(buttons["Load"]);
        buttons["Stop"] = new QPushButton("Stop");
        buttons["Stop"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        loadRow->addWidget(buttons["Stop"]);
        buttons["Run"] = new QPushButton("Run");
        buttons["Run"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        loadRow->addWidget(buttons["Run"]);
        buttons["Save"] = new QPushButton("\xf0\x9f\x92\xbe");
        buttons["Save"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        loadRow->addWidget(buttons["Save"]);

        QHBoxLayout* fitRow = getRowBoxLayout(entryColumn2Layout);
        addPulldownInContainer(pulldownContainerWidth,mainButtonHeight,fitRow,"fit");
        pulldowns["fit"]->addItem(("Maximize x"));
        pulldowns["fit"]->addItem(("Maximize y"));
        pulldowns["fit"]->addItem(("Maximize Total"));
        pulldowns["fit"]->addItem(("Fit spectrum"));
        pulldowns["fit"]->addItem(("Fit spectrum (log)"));
        buttons["LoadFitting"] = new QPushButton("\xf0\x9f\x93\x82");
        buttons["LoadFitting"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        fitRow->addWidget(buttons["LoadFitting"]);
        buttons["Fit"] = new QPushButton("Fit");
        buttons["Fit"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        fitRow->addWidget(buttons["Fit"]);

        fitting = new QTextEdit;
        fitting->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        fitting->setFixedSize(rowWidth,2*rowHeight);
        fitting->setTabChangesFocus(true);
        entryColumn2Layout->addWidget(fitting);

        QSpacerItem* consoleSpacer = new QSpacerItem(rowWidth,1);
        entryColumn2Layout->addSpacerItem(consoleSpacer);

        console = new QTextEdit;
        console->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        console->setFixedSize(rowWidth,3*rowHeight-1);
        console->setTabChangesFocus(true);
        entryColumn2Layout->addWidget(console);

        //Divide the sequence box into a row of buttons and the text box
        QVBoxLayout* sequenceBoxLayout = new QVBoxLayout(sequenceBox);
        squeezeMargins(sequenceBoxLayout);
        QWidget* sequenceButtonBox = new QWidget;
        //sequenceButtonBox->setFixedSize(4*textBoxWidth + 2*labelWidth,textBoxHeight);
        sequenceButtonBox->setFixedHeight(textBoxHeight);
        sequenceBoxLayout->addWidget(sequenceButtonBox);
        QHBoxLayout* sequenceButtonBoxLayout = new QHBoxLayout(sequenceButtonBox);
        squeezeMargins(sequenceButtonBoxLayout);
        labels["sequence"] = new QLabel;
        labels["sequence"]->setText("Sequence:");
        labels["sequence"]->setFixedWidth(labelWidth/2);
        sequenceButtonBoxLayout->addWidget(labels["sequence"]);
        auto addMiniButton = [&](const QString& icon, const std::string& entry, const QString& tooltip, std::function<void()> action){
            buttons[entry] = new QPushButton(icon);
            buttons[entry]->setFixedWidth(miniButtonWidth);
            buttons[entry]->setToolTip(tooltip);
            sequenceButtonBoxLayout->addWidget(buttons[entry]);
            QObject::connect(buttons[entry], &QPushButton::clicked, action);
        };

        addMiniButton("\xf0\x9f\x93\xb8", "addSameCrystal", 
        "Add a fixed crystal which matches the values currently entered on the interface",[&](){
            if(textBoxes["NLAbsorption"]->text().toDouble() != 0.0){
                sPrint("plasma({},{},{},{},{},{},{},{},{})",
                pulldowns["material"]->currentIndex(), 
                textBoxes["CrystalTheta"]->text().toDouble(),
                textBoxes["CrystalPhi"]->text().toDouble(), 
                textBoxes["NLAbsorption"]->text().toDouble(),
                textBoxes["CrystalBandgap"]->text().toDouble(), 
                textBoxes["DrudeGamma"]->text().toDouble(),
                textBoxes["effectiveMass"]->text().toDouble(), 
                textBoxes["ZSize"]->text().toDouble(),
                textBoxes["dz"]->text().toDouble());
            }
            else{
                sPrint("nonlinear({},{},{},{},{})",
                pulldowns["material"]->currentIndex(), 
                textBoxes["CrystalTheta"]->text().toDouble(),
                textBoxes["CrystalPhi"]->text().toDouble(), 
                textBoxes["ZSize"]->text().toDouble(),
                textBoxes["dz"]->text().toDouble());
            }

        });
        addMiniButton("\xe2\x99\x8a", "addDefault", "Insert a crystal that will change with the values set on the "
            "interface, or modified during a batch calculation",[&](){
                sPrint("plasma(d,d,d,d,d,d,d,d,d)");
            });
        addMiniButton("\xf0\x9f\x92\xab", "addRotation", "Rotate the polarization by a specified angle in degrees",[&](){
            sPrint("rotate(90)");
        });
        addMiniButton("\xf0\x9f\x92\xa1", "addPulse", "Add a new pulse to the grid; values will be set to duplicate "
            "pulse 1 as entered above",[&](){
                sPrint("addPulse({},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})",
                textBoxes["Energy1"]->text().toDouble(), 
                textBoxes["Frequency1"]->text().toDouble(),
                textBoxes["Bandwidth1"]->text().toDouble(), 
                textBoxes["SGOrder1"]->text().toInt(),
                textBoxes["CEP1"]->text().toDouble(), 
                textBoxes["Delay1"]->text().toDouble(), 
                textBoxes["GDD1"]->text().toDouble(),
                textBoxes["TOD1"]->text().toDouble(),
                textBoxes["Material1"]->text().toInt(), 
                textBoxes["Thickness1"]->text().toDouble(),
                textBoxes["Beamwaist1"]->text().toDouble(),
                textBoxes["xOffset1"]->text().toDouble(),
                0.0,
                textBoxes["zOffset1"]->text().toDouble(),
                textBoxes["NCAngle1"]->text().toDouble(),
                0.0,
                textBoxes["Polarization1"]->text().toDouble(),
                textBoxes["Circularity1"]->text().toDouble(),
                pulldowns["material"]->currentIndex(),
                textBoxes["CrystalTheta"]->text().toDouble(),
                textBoxes["CrystalPhi"]->text().toDouble());
            });
        addMiniButton("\xf0\x9f\x94\x8e", "addMirror", "Add a spherical mirror to the beam path, with radius "
            "of curvature in meters",[&](){
                sPrint("sphericalMirror(-1.0)");
            });
        addMiniButton("\xf0\x9f\x98\x8e", "addFilter", "Add a spectral filter to the beam path. "
            "Parameters:\n   central frequency (THz)\n   bandwidth (THz)\n   supergaussian order\n   "
            "in-band amplitude\n   out-of-band amplitude\n",[&](){
                sPrint("filter(130, 20, 4, 1, 0)");
            });
        addMiniButton("\xf0\x9f\x93\x8f", "addLinear", "Add a linear propagation through the crystal entered on the interface",[&](){
            sPrint("linear({},{},{},{},{})",
                pulldowns["material"]->currentIndex(), 
                textBoxes["CrystalTheta"]->text().toDouble(),
                textBoxes["CrystalPhi"]->text().toDouble(), 
                textBoxes["ZSize"]->text().toDouble(),
                textBoxes["dz"]->text().toDouble());
            });
        addMiniButton("\xf0\x9f\x8e\xaf", "addAperture", "Add an aperture to the beam. Parameters:\n   diameter (m)\n   "
            "activation parameter\n",[&](){
                sPrint("aperture(0.001, 2)");
            });
        addMiniButton("\xe2\x9b\xb3", "addFarFieldAperture", "Filter the beam with a far-field aperture. Parameters:\n   "
            "opening angle (deg)\n   activation parameter (k)\n   x-angle (deg)\n   y-angle (deg) ",[&](){
                sPrint("farFieldAperture(2.0,4000,0,0)");
            });
        addMiniButton("\xf0\x9f\x94\x81", "addForLoop", "Add an empty for loop. Parameters:\n   "
            "Number of times to execute\n   Variable number in which to put the counter",[&](){
                sPrint("for(10,1){{\n\n}}");
            });
        QSpacerItem* miniButtonSpacer = new QSpacerItem(1,1,QSizePolicy::Expanding,QSizePolicy::Fixed);
        sequenceButtonBoxLayout->addSpacerItem(miniButtonSpacer);
        sequence = new QTextEdit;
        sequence->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
        sequence->setFixedWidth(2*rowWidth);
        sequence->setTabChangesFocus(true);
        sequenceBoxLayout->addWidget(sequence);
        
        //Put the control strip below the sequence
        QHBoxLayout* simulationControlStripLayout = new QHBoxLayout(simulationControlStrip);
        squeezeMargins(simulationControlStripLayout);
        progress = new QProgressBar;
        simulationControlStripLayout->addWidget(progress);

        checkboxes["FP64"] = new QCheckBox("FP64");
        simulationControlStripLayout->addWidget(checkboxes["FP64"]);
        addPulldownInContainer(textBoxWidth,mainButtonHeight,simulationControlStripLayout,"primaryHardware");
        addPulldownInContainer(textBoxWidth,mainButtonHeight,simulationControlStripLayout,"secondaryHardware");
        textBoxes["offload"] = new QLineEdit;
        textBoxes["offload"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        simulationControlStripLayout->addWidget(textBoxes["offload"]);

        QHBoxLayout* plotControlStripLayout = new QHBoxLayout(plotControlStrip);
        buttons["collapse"] = new QPushButton("\xe2\x86\x94\xef\xb8\x8f");
        buttons["collapse"]->setFixedSize(miniButtonWidth, textBoxHeight);
        plotControlStripLayout->addWidget(buttons["collapse"]);
        squeezeMargins(plotControlStripLayout);
        slider = new QSlider(Qt::Horizontal);
        plotControlStripLayout->addWidget(slider);
        buttons["SVG"] = new QPushButton("SVG");
        buttons["SVG"]->setFixedSize(miniButtonWidth, textBoxHeight);
        plotControlStripLayout->addWidget(buttons["SVG"]);
        buttons["xlim"] = new QPushButton("xlim");
        buttons["xlim"]->setFixedSize(miniButtonWidth, textBoxHeight);
        plotControlStripLayout->addWidget(buttons["xlim"]);
        textBoxes["xlimStart"] = new QLineEdit;
        textBoxes["xlimStart"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(textBoxes["xlimStart"]);
        textBoxes["xlimStop"] = new QLineEdit;
        textBoxes["xlimStop"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(textBoxes["xlimStop"]);
        textBoxes["ylimStart"] = new QLineEdit;
        buttons["ylim"] = new QPushButton("ylim");
        buttons["ylim"]->setFixedSize(miniButtonWidth, textBoxHeight);
        plotControlStripLayout->addWidget(buttons["ylim"]);
        textBoxes["ylimStart"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(textBoxes["ylimStart"]);
        textBoxes["ylimStop"] = new QLineEdit;
        textBoxes["ylimStop"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(textBoxes["ylimStop"]);
        checkboxes["Total"] = new QCheckBox("Total");
        plotControlStripLayout->addWidget(checkboxes["Total"]);
        checkboxes["Log"] = new QCheckBox("Log");
        plotControlStripLayout->addWidget(checkboxes["Log"]);

        readDefaultValues(theSim, theDatabase);
        setInterfaceValuesToActiveValues(theSim);
        cPrint(checkLibraryAvailability(theSim));

        //populate hardware selectors
        std::string A;
        if (theSim.base().CUDAavailable) {
            pulldowns["primaryHardware"]->addItem("CUDA");
            pulldowns["secondaryHardware"]->addItem("CUDA");
            for (int i = 1; i < theSim.base().cudaGPUCount; ++i) {
                A = Sformat("CUDA {}", i);
                pulldowns["primaryHardware"]->addItem(A.c_str());
                pulldowns["secondaryHardware"]->addItem(A.c_str());
            }
        }
        if (theSim.base().SYCLavailable) {
            A.assign("SYCL");
            pulldowns["primaryHardware"]->addItem(A.c_str());
            pulldowns["secondaryHardware"]->addItem(A.c_str());
            if (theSim.base().syclGPUCount > 0) {

                pulldowns["primaryHardware"]->addItem("SYCLc");
                pulldowns["secondaryHardware"]->addItem("SYCLc");
                pulldowns["primaryHardware"]->addItem("SYCLg");
                pulldowns["secondaryHardware"]->addItem("SYCLg");
            }
        }
#if defined _WIN32 || defined __linux__ && not defined CPUONLY
        pulldowns["primaryHardware"]->addItem("C++");
        pulldowns["secondaryHardware"]->addItem("C++");
#elif defined __APPLE__
        pulldowns["primaryHardware"]->addItem("\xef\xa3\xbfGCD");
        pulldowns["secondaryHardware"]->addItem("\xef\xa3\xbfGCD");
#endif
#if defined _WIN32 || defined __linux__
        pulldowns["primaryHardware"]->addItem("OpenMP");
        pulldowns["secondaryHardware"]->addItem("OpenMP");
#endif

        messengerThread = new QThread;
        messenger = new GuiMessenger;
        messenger->moveToThread(messengerThread);
        QObject::connect(messenger, &GuiMessenger::sendText, console, &QTextEdit::append);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["FreqPlot1"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["FreqPlot2"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["TimePlot1"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["TimePlot2"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["TimeImage1"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["TimeImage2"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["FreqImage1"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["FreqImage2"], &CairoWidget::queueUpdate);

        QObject::connect(slider, &QSlider::valueChanged, plots["FreqPlot1"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["FreqPlot2"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["TimePlot1"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["TimePlot2"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["TimeImage1"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["TimeImage2"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["FreqImage1"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["FreqImage2"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::moveSlider, slider, &QSlider::setValue);
        QObject::connect(messenger, &GuiMessenger::requestSyncValues, this, &LWEGui::queueSyncValues);
        QObject::connect(messenger, &GuiMessenger::requestSliderUpdate, this, &LWEGui::updateSlider);
        QObject::connect(messenger, &GuiMessenger::sendProgressRange, progress, &QProgressBar::setRange);
        QObject::connect(messenger, &GuiMessenger::sendProgressValue, progress, &QProgressBar::setValue);

        windowBody->show();
        connectButtons();
    }
    void connectButtons(){
        QObject::connect(buttons["Run"], &QPushButton::clicked, [&](){
            readParametersFromInterface(theSim);
            theSim.configure();
            updateSlider();
            simulationRun theRun(pulldowns["primaryHardware"]->currentIndex(),checkboxes["FP64"]->isChecked(),theSim);
            simulationRun theOffloadRun(pulldowns["secondaryHardware"]->currentIndex(),checkboxes["FP64"]->isChecked(),theSim);
            std::thread(mainSimThread, std::ref(*this), theRun, theOffloadRun).detach();
        });

        QObject::connect(buttons["Fit"], &QPushButton::clicked, [&](){
            readParametersFromInterface(theSim);
            theSim.configure();
            simulationRun theRun(pulldowns["primaryHardware"]->currentIndex(),checkboxes["FP64"]->isChecked(),theSim);
            std::thread(fittingThread, std::ref(*this), theRun).detach();
        });

        QObject::connect(buttons["Save"], &QPushButton::clicked, [&](){
            theSim.base().outputBasePath = QFileDialog::getSaveFileName(buttons["Save"],"Save LWE result","","LWE Results (*.zip)").toStdString();
            stripLineBreaks(theSim.base().outputBasePath);
            if ((theSim.base().outputBasePath.length() > 4 && theSim.base().outputBasePath.substr(theSim.base().outputBasePath.length() - 4) == ".txt")
                || (theSim.base().outputBasePath.length() > 4 && theSim.base().outputBasePath.substr(theSim.base().outputBasePath.length() - 4) == ".zip")) {
                theSim.base().outputBasePath = theSim.base().outputBasePath.substr(0, theSim.base().outputBasePath.length() - 4);
            }

            messenger->passString("Saving...");
            auto saveLambda = [&](){
                theSim.saveDataSet();
                messenger->passString("done.\n");
            };
            std::thread(saveLambda).detach();
        });
        QObject::connect(buttons["saveSLURM"], &QPushButton::clicked, [&](){
            theSim.base().outputBasePath = QFileDialog::getSaveFileName(buttons["saveSLURM"],"Save cluster script","","LWE Results (*.zip)").toStdString();
            stripLineBreaks(theSim.base().outputBasePath);
            if ((theSim.base().outputBasePath.length() > 4 && theSim.base().outputBasePath.substr(theSim.base().outputBasePath.length() - 4) == ".txt")
                || (theSim.base().outputBasePath.length() > 4 && theSim.base().outputBasePath.substr(theSim.base().outputBasePath.length() - 4) == ".zip")) {
                theSim.base().outputBasePath = theSim.base().outputBasePath.substr(0, theSim.base().outputBasePath.length() - 4);
            }
            createRunFile(*this);
        });

        QObject::connect(buttons["Load"], &QPushButton::clicked, [&](){
            std::string path = QFileDialog::getOpenFileName(buttons["Load"],"Load LWE result","","LWE Results (*.zip);;LWE Inputs (*.txt)").toStdString();
            std::thread([&](std::string path){
                bool isZipFile = (path.length() >= 4 
                && path.substr(path.length()-4)==".zip");
                messenger->passString("Loading...");
                int readParameters =
                theSim.base().readInputParametersFile(theDatabase.db.data(), path);
                theSim.configure();
                std::for_each(theSim.mutexes.begin(), theSim.mutexes.end(), 
                        [](std::mutex& m) {std::lock_guard<std::mutex> lock(m); });
                if (readParameters == 61) {
                    int64_t extensionLoc = path.find_last_of(".");
                    const std::string basePath = path.substr(0, extensionLoc);
                    theSim.base().loadSavedFields(basePath, isZipFile);
                }
                messenger->requestUpdate();
                messenger->passSyncValues();
                messenger->passSliderUpdate();
                messenger->passString("done.");
            },path).detach();
        });

        QObject::connect(buttons["LoadPulse1"], &QPushButton::clicked, [&](){
            std::string path = QFileDialog::getOpenFileName(buttons["LoadPulse1"],"Load field data","","ASCII data (*.*)").toStdString();
            pulse1LoadedData = loadedInputData(path);
            messenger->passString(Sformat("Loaded new file into pulse 1 buffer:\n{}\n", pulse1LoadedData.filePath));
        });

        QObject::connect(buttons["LoadPulse2"], &QPushButton::clicked, [&](){
            std::string path = QFileDialog::getOpenFileName(buttons["LoadPulse2"],"Load field data","","ASCII data (*.*)").toStdString();
            pulse2LoadedData = loadedInputData(path);
            messenger->passString(Sformat("Loaded new file into pulse 2 buffer:\n{}\n", pulse2LoadedData.filePath));
        });

        QObject::connect(buttons["LoadFitting"], &QPushButton::clicked, [&](){
            std::string path = QFileDialog::getOpenFileName(buttons["LoadFitting"],"Load spectral target data","","ASCII data (*.*)").toStdString();
            fittingLoadedData = loadedInputData(path);
            messenger->passString(Sformat("Loaded new fitting spectral target:\n{}\n", fittingLoadedData.filePath));
        });

        QObject::connect(buttons["loadMaterial"], &QPushButton::clicked, [&](){
            std::string path = QFileDialog::getOpenFileName(buttons["LoadMaterial"],"Load crystal database","","ASCII data (*.*)").toStdString();
            theDatabase = crystalDatabase(path);
            populateDatabasePulldown();
            messenger->passString(Sformat("Loaded new crystal database:\n{}\n", path));
        });

        QObject::connect(buttons["Stop"], &QPushButton::clicked, [&](){
            if (theSim.base().isRunning) {
                theSim.base().cancellationCalled = true;
                for (int i = 1; i < theSim.base().Nsims; ++i) {
                    theSim.base().cancellationCalled = true;
                }
            }
        });

        QObject::connect(buttons["collapse"], &QPushButton::clicked, [&](){
            if(isInputRegionHidden){
                inputRegion->show();
                isInputRegionHidden = false;
            }
            else{
                inputRegion->hide();
                isInputRegionHidden = true;
            }
        });

        QObject::connect(buttons["ylim"], &QPushButton::clicked, [&](){
            messenger->passDrawRequest();
        });
        QObject::connect(buttons["xlim"], &QPushButton::clicked, [&](){
            messenger->passDrawRequest();
        });
        QObject::connect(checkboxes["Log"], &QCheckBox::clicked, [&](){
            messenger->passDrawRequest();
        });
        QObject::connect(checkboxes["Total"], &QCheckBox::clicked, [&](){
            messenger->passDrawRequest();
        });
        /*std::ofstream fs(SVGpath);
            fs << SVGstrings[0] << SVGstrings[1] << SVGstrings[2] << SVGstrings[3];*/
        QObject::connect(buttons["SVG"], &QPushButton::clicked, [&](){
            std::string SVGpath = QFileDialog::getSaveFileName(buttons["SVG"],"Save SVG file of plots","","Scalable Vector Graphics (*.svg)").toStdString();
            isMakingSVG = true;
            plots["TimePlot1"]->repaint();
            plots["TimePlot2"]->repaint();
            plots["FreqPlot1"]->repaint();
            plots["FreqPlot2"]->repaint();
            isMakingSVG = false;
            std::ofstream fs(SVGpath);
            fs << SVGStrings[0] << SVGStrings[1] << SVGStrings[2] << SVGStrings[3];
            
        });
        
    }
    public slots:
    void queueSyncValues(){
        setInterfaceValuesToActiveValues(theSim);
    }
    void updateSlider(){
        slider->setMinimum(0);
        slider->setMaximum(theSim.base().Nsims * theSim.base().Nsims2 - 1);
    }
     
};

int main(int argc, char **argv){
    QApplication app(argc, argv);
    LWEGui Gui;
    return app.exec();
}

//ifdef guards are in place to only include CUDA/SYCL when they are being used
std::string checkLibraryAvailability(simulationBatch& theSim) {   
    std::string s;
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
            s.append(Sformat("CUDA found a GPU:\n", theSim.base().cudaGPUCount));
        }
        else {
            s.append(Sformat("CUDA found {} GPU(s):\n", theSim.base().cudaGPUCount));
        }
        for (i = 0; i < theSim.base().cudaGPUCount; ++i) {
            cuErr = cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
            s.append(Sformat("   {}\n", 
                activeCUDADeviceProp.name));
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
    isIntelRuntimeInstalled = LoadLibraryA("pi_win_proxy_loader.dll"); 
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
            s.append(Sformat("{}", syclDeviceList));
        }
    }
    else {
        s.append(Sformat("Not using SYCL because the Intel DPC++\n"
            "Compiler Runtime is not installed.\n"));
        theSim.base().SYCLavailable = false;
    }
#endif
#endif
    return s;
}



void mainSimThread(LWEGui& theGui, simulationRun theRun, simulationRun theOffloadRun) {

    simulationBatch& theSim = theGui.theSim;
    int error = 0;
    theSim.base().cancellationCalled = false;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();
    
    //theGui.requestSliderUpdate();
    std::vector<simulationParameterSet> counterVector = theSim.getParameterVector();
    std::atomic_uint32_t& totalSteps = theGui.totalSteps;
    std::atomic_uint32_t& progressCounter = theGui.progressCounter;
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
        theGui.messenger->passString(Sformat(
            "<span color=\"#FF88FF\">Simulation failed with exception:\n{}</span>\n",
            errorString));
        return;
    }
    theGui.messenger->passProgressRange(0,totalSteps-2);
    
    theSim.base().isRunning = true;
    auto progressThread = [&](){
        while(theSim.base().isRunning){
            theGui.messenger->passProgressValue(theGui.progressCounter);
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
    };
    std::thread advanceProgressBar(progressThread);
    advanceProgressBar.detach();
    auto batchLoop = [&](const int startSim, const int stopSim, const simulationRun& activeRun){
        for (int j = startSim; j < stopSim; ++j) {
            theSim.sCPU()[j].runningOnCPU = activeRun.forceCPU;
            theSim.sCPU()[j].assignedGPU = activeRun.assignedGPU;
            theSim.sCPU()[j].useOpenMP = activeRun.useOpenMP;
                try {
                    if(theSim.base().isInSequence){
                        error = activeRun.sequenceFunction(&theSim.sCPU()[j]);
                    }
                    else{
                        error = activeRun.normalFunction(&theSim.sCPU()[j]);
                    }
                }
                catch (std::exception const& e) {
                    std::string errorString=e.what();
                    std::erase(errorString,'<');
                    std::erase(errorString,'>');
                    std::erase(errorString,'&');
                    std::erase(errorString,';');
                    std::erase(errorString,'{');
                    std::erase(errorString,'}');
                    theGui.messenger->passString(
                        Sformat("Simulation failed with exception:\n{}\n", 
                        errorString));
                }
                if (theSim.sCPU()[j].memoryError != 0) {
                    if (theSim.sCPU()[j].memoryError == -1) {
                        theGui.messenger->passString(
                        Sformat(
                            "Not enough free GPU memory, sorry.\n", 
                            theSim.sCPU()[j].memoryError));
                    }
                    else {
                        theGui.messenger->passString(
                        Sformat(
                            "<span color=\"#FF88FF\">Warning: device memory error ({}).</span>\n", 
                            theSim.sCPU()[j].memoryError));
                    }
                }
                if (error) break;
                theGui.messenger->passSliderPosition(j);
                theGui.messenger->requestUpdate();
            if (theSim.base().cancellationCalled) {
                theGui.messenger->passString(Sformat((
                    "Warning: series cancelled, stopping\n"
                    "after {} simulations.\n"), j + 1));
                break;
            }
        }
    };
    if(theSim.base().NsimsCPU){
        std::thread offloadThread(batchLoop,
        theSim.base().Nsims * theSim.base().Nsims2 - theSim.base().NsimsCPU,
        theSim.base().Nsims * theSim.base().Nsims2,
        theOffloadRun);
        offloadThread.detach();
        batchLoop(0,theSim.base().Nsims * theSim.base().Nsims2 - theSim.base().NsimsCPU,theRun);
        if (offloadThread.joinable()) offloadThread.join();
    }
    else{
        batchLoop(0,theSim.base().Nsims * theSim.base().Nsims2 - theSim.base().NsimsCPU,theRun);
    }
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    if (error == 13) {
        theGui.messenger->passString(
            "NaN detected in grid!\n"
            "Try using a larger spatial/temporal step\n"
            "or smaller propagation step.\n"
            "Simulation was cancelled.");
    }
    else if (error == 15) {
        theGui.messenger->passString(
            "<span color=\"#FF88FF\">"
            "Sorry, that sequence mode has been \n"
            "replaced by the new one. Look in the \n"
            "documentation for more info. It is a lot\n"
            "easier to use now, and hopefully \n"
            "it won't take long to set it up. \n"
            "Sorry about that!\n</span>");
    }
    else if(!error){
        theGui.messenger->passString(Sformat(
            "Finished after {:.4} s.", 1e-6 *
            (double)(std::chrono::duration_cast<std::chrono::microseconds>
                (simulationTimerEnd - simulationTimerBegin).count())));
    }
    theGui.messenger->passProgressValue(theGui.progressCounter);
    theSim.base().isRunning = false;
}
void fittingThread(LWEGui& theGui,  simulationRun theRun) {
    simulationBatch& theSim = theGui.theSim;
    theSim.base().cancellationCalled = false;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

    theSim.sCPU()->readFittingString();
    if (theSim.base().Nfitting == 0) {
        theGui.messenger->passString("Couldn't interpret fitting command.");
        return;
    }
    theGui.progressCounter = 0;
    theSim.base().progressCounter = &theGui.progressCounter;
    if (theSim.base().fittingMode == 3) {
        if (theSim.sCPU()->loadReferenceSpectrum()) {
            theGui.messenger->passString("Could not read reference file!");
            return;
        }
    }

    theGui.messenger->passString(Sformat(
        "Fitting {} values in mode {} over {} iterations.\n"
        "Region of interest contains {} elements\n",
        theSim.base().Nfitting, 
        theSim.base().fittingMode, 
        theSim.base().fittingMaxIterations, 
        theSim.base().fittingROIsize));

    theSim.base().isRunning = true;
    theSim.base().runningOnCPU = theRun.forceCPU;
    theSim.base().assignedGPU = theRun.assignedGPU;
    theSim.base().useOpenMP = theRun.useOpenMP;
    auto progressThread = [&](){
        theGui.messenger->passProgressRange(0,theSim.base().fittingMaxIterations-1);
        while(theSim.base().isRunning){
            theGui.messenger->passProgressValue(theGui.progressCounter);
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
    };
    std::thread advanceProgressBar(progressThread);
    advanceProgressBar.detach();
    std::unique_lock lock(theSim.mutexes.at(0));
    theRun.fittingFunction(theSim.sCPU());
    lock.unlock();
    theSim.base().plotSim = 0;

    theGui.messenger->passDrawRequest();
    //theGui.requestInterfaceValuesUpdate();
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    theGui.messenger->passString(Sformat(("Finished fitting after {:.4} s.\n"), 1e-6 *
        (double)(std::chrono::duration_cast<std::chrono::microseconds>
            (simulationTimerEnd - simulationTimerBegin).count())));
    theGui.messenger->passString(Sformat(
        "Fitting result:\n"
        "(index, value)"));
    for (int i = 0; i < theSim.base().Nfitting; ++i) {
        theGui.messenger->passString(Sformat("{},  {}", i, theSim.base().fittingResult[i]));
    }
    theSim.base().isRunning = false;
    theGui.messenger->requestSyncValues();
}

void readDefaultValues(simulationBatch& sim, crystalDatabase& db){

        //Linux search order:
        // ~/.LightwaveExplorer
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
            sim.sCPU()->readInputParametersFile(db.db.data(), homePath));
        
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

            if (1 == sim.sCPU()->readInputParametersFile(db.db.data(), defaultsPath)) {
                sim.sCPU()->readInputParametersFile(db.db.data(), "DefaultValues.ini");
            }
        }
        
#elif defined __APPLE__
		uint32_t bufferSize = 1024;
		char sysPathBuffer[1024] = { 0 };
		_NSGetExecutablePath(sysPathBuffer, &bufferSize);
        std::string sysPathFull(sysPathBuffer);
        std::string sysPathIni = sysPathFull.substr(0,sysPathFull.find_last_of("/"));
        sysPathIni.append("/../Resources/DefaultValues.ini");
		if(1 == theSim.sCPU()->readInputParametersFile(db.db.data(), sysPathIni.c_str())){
            sim.sCPU()->readInputParametersFile(db.db.data(), "DefaultValues.ini");
        }
#else
		sim.sCPU()->readInputParametersFile(db.db.data(), "DefaultValues.ini");
#endif
}

void drawTimeImage1(cairo_t* cr, int width, int height, LWEGui& theGui) {
    if (!theGui.theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;
    int64_t simIndex = maxN(0,theGui.slider->value());

    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = theGui.theSim.base().Ntime * theGui.theSim.base().Nspace * (theGui.theSim.base().Nspace2 / 2);

    sPlot.data =
        &theGui.theSim.base().ExtOut[simIndex * theGui.theSim.base().Ngrid * 2 + cubeMiddle];
    sPlot.dataXdim = theGui.theSim.base().Ntime;
    sPlot.dataYdim = theGui.theSim.base().Nspace;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.colorMap = 4;
    sPlot.dataType = 0;
    std::unique_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.imagePlot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
    }
}

void drawField1Plot(cairo_t* cr, int width, int height, LWEGui& theGui) {
    if (!theGui.theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LwePlot sPlot;
    int64_t simIndex = maxN(0,theGui.slider->value());

    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = theGui.theSim.base().Ntime * theGui.theSim.base().Nspace * (theGui.theSim.base().Nspace2 / 2);

    sPlot.makeSVG = theGui.isMakingSVG;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theGui.theSim.base().tStep / 1e-15;
    sPlot.x0 = -((sPlot.dx * theGui.theSim.base().Ntime) / 2 - sPlot.dx / 2);
    sPlot.data = &theGui.theSim.base().ExtOut[
        simIndex * theGui.theSim.base().Ngrid * 2 + cubeMiddle + theGui.theSim.base().Ntime * theGui.theSim.base().Nspace / 2];
    sPlot.Npts = theGui.theSim.base().Ntime;
    sPlot.color = LweColor(0, 1, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ex (GV/m)";
    sPlot.unitY = 1e9;
    std::unique_lock dataLock(theGui.theSim.mutexes.at(simIndex),std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG) {
            std::size_t SVGbegin = sPlot.SVGString.find("<svg");
            sPlot.SVGString.insert(SVGbegin,Sformat("<svg width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}"
                "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink="
                "\"http://www.w3.org/1999/xlink\">\n",
                2*width, 2*height, 2*width, 2*height));
            theGui.SVGStrings[0] = sPlot.SVGString;
        }
    }
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
    }
}

void drawField2Plot(cairo_t* cr, int width, int height, LWEGui& theGui) {
    if (!theGui.theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LwePlot sPlot;

    int64_t simIndex = maxN(0,theGui.slider->value());

    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = 
        theGui.theSim.base().Ntime * theGui.theSim.base().Nspace * (theGui.theSim.base().Nspace2 / 2);


    sPlot.makeSVG = theGui.isMakingSVG;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theGui.theSim.base().tStep / 1e-15;
    sPlot.x0 = -((sPlot.dx * theGui.theSim.base().Ntime) / 2 - sPlot.dx / 2);
    sPlot.data = 
        &theGui.theSim.base().ExtOut[
        theGui.theSim.base().Ngrid + simIndex * theGui.theSim.base().Ngrid * 2 
            + cubeMiddle + theGui.theSim.base().Ntime * theGui.theSim.base().Nspace / 2];
    sPlot.Npts = theGui.theSim.base().Ntime;
    sPlot.color = LweColor(1, 0, 1, 1);
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ey (GV/m)";
    sPlot.unitY = 1e9;
    std::unique_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG) {
            std::size_t SVGbegin = sPlot.SVGString.find("width=");
            sPlot.SVGString.insert(SVGbegin,Sformat("x=\"{}\" y=\"{}\" ",
                0, height));
            SVGbegin = sPlot.SVGString.find("<svg");
            theGui.SVGStrings[1] = sPlot.SVGString.substr(SVGbegin);
        }
    }
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
    }
}

void drawSpectrum1Plot(cairo_t* cr, int width, int height, LWEGui& theGui) {
    if (!theGui.theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LwePlot sPlot;

    bool logPlot = false;
    if (theGui.checkboxes["Log"]->isChecked()) {
        logPlot = true;
    }
    int64_t simIndex = maxN(0,theGui.slider->value());
    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    bool forceX = false;
    double xMin = theGui.textBoxes["xlimStart"]->text().toDouble();
    double xMax = theGui.textBoxes["xlimStop"]->text().toDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.textBoxes["ylimStart"]->text().toDouble();
    double yMax = theGui.textBoxes["ylimStop"]->text().toDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }
    bool overlayTotal = false;
    if (theGui.checkboxes["Total"]->isChecked()) {
        overlayTotal = true;
    }

    sPlot.makeSVG = theGui.isMakingSVG;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theGui.theSim.base().fStep / 1e12;
    sPlot.data = &theGui.theSim.base().totalSpectrum[simIndex * 3 * theGui.theSim.base().Nfreq];
    sPlot.Npts = theGui.theSim.base().Nfreq;
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
        sPlot.data2 = &theGui.theSim.base().totalSpectrum[(2 + simIndex * 3) * theGui.theSim.base().Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    std::unique_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG) {
            std::size_t SVGbegin = sPlot.SVGString.find("width=");
            sPlot.SVGString.insert(SVGbegin,Sformat("x=\"{}\" y=\"{}\" ",
                width, 0));
            SVGbegin = sPlot.SVGString.find("<svg");
            theGui.SVGStrings[2] = sPlot.SVGString.substr(SVGbegin);
        }
    }
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
    }
}

void drawSpectrum2Plot(cairo_t* cr, int width, int height, LWEGui& theGui) {
    if (!theGui.theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }

    LwePlot sPlot;
    bool logPlot = false;
    if (theGui.checkboxes["Log"]->isChecked()) {
        logPlot = true;
    }
    int64_t simIndex = maxN(0,theGui.slider->value());
    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    bool forceX = false;
    double xMin = theGui.textBoxes["xlimStart"]->text().toDouble();
    double xMax = theGui.textBoxes["xlimStop"]->text().toDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = theGui.textBoxes["ylimStart"]->text().toDouble();
    double yMax = theGui.textBoxes["ylimStop"]->text().toDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }
    bool overlayTotal = false;
    if (theGui.checkboxes["Total"]->isChecked()) {
        overlayTotal = true;
    }
    sPlot.makeSVG = theGui.isMakingSVG;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theGui.theSim.base().fStep / 1e12;
    sPlot.data = &theGui.theSim.base().totalSpectrum[(1 + simIndex * 3) * theGui.theSim.base().Nfreq];
    sPlot.Npts = theGui.theSim.base().Nfreq;
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
        sPlot.data2 = &theGui.theSim.base().totalSpectrum[(2 + simIndex * 3) * theGui.theSim.base().Nfreq];
        sPlot.ExtraLines = 1;
        sPlot.color2 = LweColor(1.0, 0.5, 0.0, 0);
    }
    std::unique_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG) {
            std::size_t SVGbegin = sPlot.SVGString.find("width=");
            sPlot.SVGString.insert(SVGbegin,Sformat("x=\"{}\" y=\"{}\" ",
                width, height));
            SVGbegin = sPlot.SVGString.find("<svg");
            theGui.SVGStrings[3] = sPlot.SVGString.substr(SVGbegin).append("</svg>");
        }
    }
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
    }
}

void drawTimeImage2(cairo_t* cr, int width, int height, LWEGui& theGui) {
    if (!theGui.theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;
    int64_t simIndex = maxN(0,theGui.slider->value());
    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = theGui.theSim.base().Ntime * theGui.theSim.base().Nspace * (theGui.theSim.base().Nspace2 / 2);

    sPlot.data =
    &theGui.theSim.base().ExtOut[theGui.theSim.base().Ngrid + simIndex * theGui.theSim.base().Ngrid * 2 + cubeMiddle];
    sPlot.dataYdim = theGui.theSim.base().Nspace;
    sPlot.dataXdim = theGui.theSim.base().Ntime;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.colorMap = 4;
    sPlot.dataType = 0;
    std::unique_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.imagePlot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        //theGui.requestPlotUpdate();
    }
}

void drawFourierImage1(cairo_t* cr, int width, int height, LWEGui& theGui) {
    if (!theGui.theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;
    int64_t simIndex = maxN(0,theGui.slider->value());
    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }
    double logPlotOffset = (double)(1e-4 / (theGui.theSim.base().spatialWidth * theGui.theSim.base().timeSpan));
    if (theGui.theSim.base().is3D) {
        logPlotOffset = 
            (double)(1e-4 
                / (theGui.theSim.base().spatialWidth * theGui.theSim.base().spatialHeight * theGui.theSim.base().timeSpan));
    }
    sPlot.complexData =
        &theGui.theSim.base().EkwOut[simIndex * theGui.theSim.base().NgridC * 2];
    sPlot.dataXdim = theGui.theSim.base().Nfreq;
    sPlot.dataYdim = theGui.theSim.base().Nspace;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dataType = 1;
    sPlot.colorMap = 3;
    sPlot.logMin = logPlotOffset;
    std::unique_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.imagePlot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        //theGui.requestPlotUpdate();
    }
}

void drawFourierImage2(cairo_t* cr, int width, int height, LWEGui& theGui) {
    if (!theGui.theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;

    int64_t simIndex = maxN(0,theGui.slider->value());
    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    double logPlotOffset = (double)(1e-4 / (theGui.theSim.base().spatialWidth * theGui.theSim.base().timeSpan));
    if (theGui.theSim.base().is3D) {
        logPlotOffset = (double)(1e-4 
            / (theGui.theSim.base().spatialWidth * theGui.theSim.base().spatialHeight * theGui.theSim.base().timeSpan));
    }
    sPlot.complexData =
        &theGui.theSim.base().EkwOut[simIndex * theGui.theSim.base().NgridC * 2 + theGui.theSim.base().NgridC];
    sPlot.dataXdim = theGui.theSim.base().Nfreq;
    sPlot.dataYdim = theGui.theSim.base().Nspace;
    sPlot.height = height;
    sPlot.width = width;
    sPlot.colorMap = 3;
    sPlot.dataType = 1;
    sPlot.logMin = logPlotOffset;
    std::unique_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) sPlot.imagePlot(cr);
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        //theGui.requestPlotUpdate();
    }
}


void createRunFile(LWEGui& theGui) {
    simulationBatch& theSim = theGui.theSim;
    theGui.readParametersFromInterface(theSim);
    theSim.base().runType = runTypes::normal;
    theSim.base().isGridAllocated = false;
    theSim.base().isFollowerInSequence = false;
    theSim.base().crystalDatabase = theGui.theDatabase.db.data();
    theSim.configureCounter();

    std::vector<simulationParameterSet> counterVector = theSim.getParameterVector();
    theGui.totalSteps = 0;
    for (int64_t j = 0; j < theSim.base().Nsims * theSim.base().Nsims2; j++) {
        counterVector[j].progressCounter = &theGui.totalSteps;
        counterVector[j].runType = runTypes::counter;
        if (theSim.base().isInSequence) {
            solveNonlinearWaveEquationSequenceCounter(&counterVector[j]);
        }
        else {
            solveNonlinearWaveEquationCounter(&counterVector[j]);
        }
    }

    //create SLURM script
    int cluster = theGui.pulldowns["slurm"]->currentIndex();
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
    double timeEstimate = theSim.sCPU()->saveSlurmScript(gpuType, gpuCount, arrayMode, theGui.totalSteps, theSim.parameters, theGui.theDatabase);

    theGui.messenger->passString(Sformat(
        "Run {} on cluster with:\nsbatch {}.slurmScript\n",
        getBasename(theSim.base().outputBasePath), getBasename(theSim.base().outputBasePath)));
    theGui.messenger->passString(Sformat("or\n./lweget.sh {}\n",getBasename(theSim.base().outputBasePath)));
    theGui.messenger->passString(Sformat(
        "Upper estimate time to complete: {:.2} hours\n", 
        timeEstimate));
    theSim.base().isRunning = false;
    theSim.base().isGridAllocated = false;
}

int insertAfterCharacter(std::string& s, char target, std::string appended){
    for(std::size_t i = 0; i < s.length(); ++i){
        if(s[i] == target){
            s.insert(i+1,appended);
            i += appended.length();
        }
    }
    return 0;
}


#include "LightwaveExplorerFrontendQT.moc"