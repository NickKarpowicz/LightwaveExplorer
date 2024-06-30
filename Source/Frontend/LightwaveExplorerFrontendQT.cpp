#include "LightwaveExplorerFrontendQT.h"

class CairoWidget : public QWidget {
    simulationBatch& theSim;
    CairoFunction theFunction;
public:
    CairoWidget(simulationBatch& sim, CairoFunction fun, QWidget* parent = nullptr) : theSim(sim), theFunction(fun), QWidget(parent) {
    }
    
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
            theSim);
        cairo_destroy(cr);
        cairo_surface_destroy(surface);
        QPainter painter(this);
        painter.drawImage(0, 0, image);
    }
};

void drawTimeImage1(cairo_t* cr, int width, int height, simulationBatch& theSim) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;
    int64_t simIndex = 0;//maxN(0,theGui.plotSlider.getIntValue());
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
    }
}

void drawField1Plot(cairo_t* cr, int width, int height, simulationBatch& theSim) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LwePlot sPlot;
    bool saveSVG = false; //theGui.SVGqueue[0];
    int64_t simIndex = 0; //maxN(0,theGui.plotSlider.getIntValue());

    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = theSim.base().Ntime * theSim.base().Nspace * (theSim.base().Nspace2 / 2);

    sPlot.makeSVG = saveSVG;
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
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        // if(saveSVG) {
        //     size_t SVGbegin = sPlot.SVGString.find("<svg");
        //     sPlot.SVGString.insert(SVGbegin,Sformat("<svg width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}"
        //         "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink="
        //         "\"http://www.w3.org/1999/xlink\">\n",
        //         2*width, 2*height, 2*width, 2*height));
        //     theGui.SVGstrings[0] = sPlot.SVGString;
        //     if(saveSVG) theGui.SVGqueue[0] = false;
        // }
    }
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        //theGui.requestPlotUpdate();
    }
}

void drawField2Plot(cairo_t* cr, int width, int height, simulationBatch& theSim) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LwePlot sPlot;

    bool saveSVG = false; //theGui.SVGqueue[1];

    int64_t simIndex = 0;//maxN(0,theGui.plotSlider.getIntValue());

    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = 
        theSim.base().Ntime * theSim.base().Nspace * (theSim.base().Nspace2 / 2);


    sPlot.makeSVG = saveSVG;
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
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        // if(saveSVG) {
        //     size_t SVGbegin = sPlot.SVGString.find("width=");
        //     sPlot.SVGString.insert(SVGbegin,Sformat("x=\"{}\" y=\"{}\" ",
        //         0, height));
        //     SVGbegin = sPlot.SVGString.find("<svg");
        //     theGui.SVGstrings[1] = sPlot.SVGString.substr(SVGbegin);
        //     if(saveSVG) theGui.SVGqueue[1] = false;
        // }
    }
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        //theGui.requestPlotUpdate();
    }
}

void drawSpectrum1Plot(cairo_t* cr, int width, int height, simulationBatch& theSim) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LwePlot sPlot;
    bool saveSVG = false; //theGui.SVGqueue[2];

    bool logPlot = false;
    // if (theGui.checkBoxes["Log"].isChecked()) {
    //     logPlot = true;
    // }
    int64_t simIndex = 0; //maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }

    bool forceX = false;
    double xMin = 0;//theGui.textBoxes[48].valueDouble();
    double xMax = 0;//theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = 0;//theGui.textBoxes[50].valueDouble();
    double yMax = 0;//theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }
    bool overlayTotal = false;
    // if (theGui.checkBoxes["Total"].isChecked()) {
    //     overlayTotal = true;
    // }

    sPlot.makeSVG = saveSVG;
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
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        // if(saveSVG) {
        //     size_t SVGbegin = sPlot.SVGString.find("width=");
        //     sPlot.SVGString.insert(SVGbegin,Sformat("x=\"{}\" y=\"{}\" ",
        //         width, 0));
        //     SVGbegin = sPlot.SVGString.find("<svg");
        //     theGui.SVGstrings[2] = sPlot.SVGString.substr(SVGbegin);
        //     if(saveSVG) theGui.SVGqueue[2] = false;
        // }
    }
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        //theGui.requestPlotUpdate();
    }
}

void drawSpectrum2Plot(cairo_t* cr, int width, int height, simulationBatch& theSim) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }

    LwePlot sPlot;

    bool saveSVG = false;//theGui.SVGqueue[3];

    bool logPlot = false;
    // if (theGui.checkBoxes["Log"].isChecked()) {
    //     logPlot = true;
    // }
    int64_t simIndex = 0;//maxN(0,theGui.plotSlider.getIntValue());
    if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
        simIndex = 0;
    }

    bool forceX = false;
    double xMin = 0;//theGui.textBoxes[48].valueDouble();
    double xMax = 0;//theGui.textBoxes[49].valueDouble();
    if (xMin != xMax && xMax > xMin) {
        forceX = true;
    }
    bool forceY = false;
    double yMin = 0;//theGui.textBoxes[50].valueDouble();
    double yMax = 0;//theGui.textBoxes[51].valueDouble();
    if (yMin != yMax && yMax > yMin) {
        forceY = true;
    }
    bool overlayTotal = false;
    // if (theGui.checkBoxes["Total"].isChecked()) {
    //     overlayTotal = true;
    // }
    sPlot.makeSVG = saveSVG;
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
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        // if(saveSVG) {
        //     size_t SVGbegin = sPlot.SVGString.find("width=");
        //     sPlot.SVGString.insert(SVGbegin,Sformat("x=\"{}\" y=\"{}\" ",
        //         width, height));
        //     SVGbegin = sPlot.SVGString.find("<svg");
        //     theGui.SVGstrings[3] = sPlot.SVGString.substr(SVGbegin).append("</svg>");
        //     if(saveSVG) theGui.SVGqueue[3] = false;
        // }
    }
    else {
        LweColor black(0, 0, 0, 0);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        //theGui.requestPlotUpdate();
    }
}

void drawTimeImage2(cairo_t* cr, int width, int height, simulationBatch& theSim) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;
    int64_t simIndex = 0;//maxN(0,theGui.plotSlider.getIntValue());
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
        //theGui.requestPlotUpdate();
    }
}

void drawFourierImage1(cairo_t* cr, int width, int height, simulationBatch& theSim) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;
    int64_t simIndex = 0;//maxN(0,theGui.plotSlider.getIntValue());
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
        //theGui.requestPlotUpdate();
    }
}

void drawFourierImage2(cairo_t* cr, int width, int height, simulationBatch& theSim) {
    if (!theSim.base().isGridAllocated) {
        LweColor black(0, 0, 0, 0);
        std::unique_lock<std::mutex> GTKlock(GTKmutex);
        cairo_rectangle(cr, 0, 0, width, height);
        black.setCairo(cr);
        cairo_fill(cr);
        return;
    }
    LweImage sPlot;

    int64_t simIndex = 0;// maxN(0,theGui.plotSlider.getIntValue());
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
        //theGui.requestPlotUpdate();
    }
}


class LWEGui {
    //QMainWindow mainWindow;
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
    template<typename... Args> void cPrint(std::string_view format, Args&&... args) {
        std::string s = Svformat(format, Smake_format_args(args...));
        console->append(s.c_str());
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

        sim.base().materialIndex = pulldowns["material"]->currentIndex();
        setToDouble(textBoxes["CrystalTheta"],sim.base().crystalTheta);
        setToDouble(textBoxes["CrystalPhi"],sim.base().crystalPhi);
        setToDouble(textBoxes["NLAbsorption"],sim.base().nonlinearAbsorptionStrength);
        setToDouble(textBoxes["CrystalBandgap"],sim.base().bandGapElectronVolts);
        setToDouble(textBoxes["DrudeGamma"],sim.base().drudeGamma);
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

        setToDouble(textBoxes["CrystalTheta"],sim.base().crystalTheta);
        setToDouble(textBoxes["CrystalPhi"],sim.base().crystalPhi);
        setToDouble(textBoxes["NLAbsorption"],sim.base().nonlinearAbsorptionStrength);
        setToDouble(textBoxes["CrystalBandgap"],sim.base().bandGapElectronVolts);
        setToDouble(textBoxes["DrudeGamma"],sim.base().drudeGamma);
        setToDouble(textBoxes["effectiveMass"],sim.base().effectiveMass);
        setToDouble(textBoxes["XSize"],1e6*sim.base().spatialWidth);
        setToDouble(textBoxes["dx"],1e6*sim.base().rStep);
        setToDouble(textBoxes["timeSpan"],1e15*sim.base().timeSpan);
        setToDouble(textBoxes["dt"],1e15*sim.base().tStep);
        setToDouble(textBoxes["ZSize"],1e6*sim.base().crystalThickness);
        setToDouble(textBoxes["dz"],1e9*sim.base().propagationStep);

    }
    LWEGui(){
        const int textBoxWidth = 80;
        const int labelWidth = 160;
        const int textBoxHeight = 26;
        const int rowHeight = textBoxHeight+2;
        const int rowWidth = labelWidth + 2*textBoxWidth + 10;
        const int miniButtonWidth = 30;
        const int mainButtonWidth = rowWidth/4;
        const int pulldownWidth = 2*textBoxWidth + 30;
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
        QWidget *inputRegion = new QWidget;
        QWidget *plotRegion = new QWidget;

        QHBoxLayout *mainAreaLayout = new QHBoxLayout(mainArea);
        squeezeMargins(mainAreaLayout);
        inputRegion->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
        plotRegion->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        mainAreaLayout->addWidget(inputRegion);
        mainAreaLayout->addWidget(plotRegion);

        QGridLayout* plotRegionLayout = new QGridLayout(plotRegion);
        plots["TimeImage1"] = new CairoWidget(theSim, drawTimeImage1);
        plots["TimeImage1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["TimeImage1"],0,0);
        plots["TimeImage2"] = new CairoWidget(theSim, drawTimeImage2);
        plots["TimeImage2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["TimeImage2"],1,0);
        plots["FreqImage1"] = new CairoWidget(theSim, drawFourierImage1);
        plots["FreqImage1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["FreqImage1"],0,1);
        plots["FreqImage2"] = new CairoWidget(theSim, drawFourierImage2);
        plots["FreqImage2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["FreqImage2"],1,1);

        plots["TimePlot1"] = new CairoWidget(theSim, drawField1Plot);
        plots["TimePlot1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["TimePlot1"],2,0);
        plots["TimePlot2"] = new CairoWidget(theSim, drawField2Plot);
        plots["TimePlot2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["TimePlot2"],3,0);
        plots["FreqPlot1"] = new CairoWidget(theSim, drawSpectrum1Plot);
        plots["FreqPlot1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["FreqPlot1"],2,1);
        plots["FreqPlot2"] = new CairoWidget(theSim, drawSpectrum2Plot);
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

        std::string materialString;
        for (int i = 0; i < theDatabase.db.size(); ++i) {
            materialString = Sformat(
                "{:2}: {}", i, std::string(theDatabase.db[i].crystalName.c_str()));
            pulldowns["material"]->addItem(materialString.c_str());
        }
        
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
        auto addMiniButton = [&](const QString& icon, const std::string& entry, const QString& tooltip){
            buttons[entry] = new QPushButton(icon);
            buttons[entry]->setFixedWidth(miniButtonWidth);
            buttons[entry]->setToolTip(tooltip);
            sequenceButtonBoxLayout->addWidget(buttons[entry]);
        };
        addMiniButton("\xf0\x9f\x93\xb8", "addSameCrystal", "Add a fixed crystal which matches the values currently entered on the interface");
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
        labels["FP64"] = new QLabel("FP64");
        simulationControlStripLayout->addWidget(labels["FP64"]);
        checkboxes["FP64"] = new QCheckBox;
        simulationControlStripLayout->addWidget(checkboxes["FP64"]);
        addPulldownInContainer(textBoxWidth,mainButtonHeight,simulationControlStripLayout,"primaryHardware");
        addPulldownInContainer(textBoxWidth,mainButtonHeight,simulationControlStripLayout,"secondaryHardware");
        textBoxes["offload"] = new QLineEdit;
        textBoxes["offload"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        simulationControlStripLayout->addWidget(textBoxes["offload"]);

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


        windowBody->show();
        //mainSimThread(0,0,false);
        readParametersFromInterface(theSim);
        mainSimThread(theSim,theDatabase,totalSteps,progressCounter,0,0,false);
    } 
};

int main(int argc, char **argv){
    QApplication app(argc, argv);
    LWEGui Gui;
    return app.exec();
}

//ifdef guards are in place to only include CUDA/SYCL when they are being used
std::string checkLibraryAvailability(simulationBatch& theSim) {   
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
    std::string s;
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

void mainSimThread(simulationBatch& theSim, crystalDatabase& theDatabase,std::atomic_uint32_t& totalSteps, std::atomic_uint32_t& progressCounter, int pulldownSelection, int secondPulldownSelection, bool use64bitFloatingPoint) {
    int error = 0;
    theSim.base().cancellationCalled = false;
    auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

    //readParametersFromInterface();
    theSim.configure();
    //theGui.requestSliderUpdate();


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
        // theGui.console.tPrint(
        //     "<span color=\"#FF88FF\">Simulation failed with exception:\n{}</span>\n",
        //     errorString);
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
        // if (theGui.firstSYCLsimulation) theGui.console.tPrint(
        //     "Note: the first time you run SYCL, it will\n"
        //     "take some time to compile kernels for your\n"
        //     "device. Subsequent runs will be faster.\n");
        // theGui.firstSYCLsimulation = false;
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
        // if (theGui.firstSYCLsimulation) theGui.console.tPrint(
        //     "Note: the first time you run SYCL, it will\n"
        //     "take some time to compile kernels for your\n"
        //     "device. Subsequent runs will be faster.\n");
        // theGui.firstSYCLsimulation = false;
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
        // if (theGui.firstSYCLsimulation) theGui.console.tPrint(
        //     "Note: the first time you run SYCL, it will\n"
        //     "take some time to compile kernels for your\n"
        //     "device. Subsequent runs will be faster.\n");
        // theGui.firstSYCLsimulation = false;
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

    // std::thread secondQueueThread(secondaryQueue, 
    //     theSim.base().Nsims * theSim.base().Nsims2 - theSim.base().NsimsCPU, 
    //     secondPulldownSelection, pulldownSelection, use64bitFloatingPoint);
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
                // theGui.console.tPrint(
                //     "<span color=\"#FF88FF\">Simulation failed with exception:\n{}</span>\n", 
                //     errorString);
            }
            if (theSim.sCPU()[j].memoryError != 0) {
                // if (theSim.sCPU()[j].memoryError == -1) {
                //     theGui.console.tPrint((
                //         "<span color=\"#FF88FF\">Not enough free GPU memory, sorry.</span>\n"), 
                //         theSim.sCPU()[j].memoryError);
                // }
                // else {
                //     theGui.console.tPrint((
                //         "<span color=\"#FF88FF\">Warning: device memory error ({}).</span>\n"), 
                //         theSim.sCPU()[j].memoryError);
                // }
            }
            if (error) break;
            //theGui.requestSliderMove(j);
            //independentPlotQueue();
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
                // theGui.console.tPrint(
                //     "<span color=\"#FF88FF\">Simulation failed with exception:\n{}</span>\n", 
                //     errorString);
            }
            
            if (theSim.sCPU()[j].memoryError != 0) {
                if (theSim.sCPU()[j].memoryError == -1) {
                    // theGui.console.tPrint((
                    //     "<span color=\"#FF88FF\">Not enough free GPU memory, sorry.</span>\n"), 
                    //     theSim.sCPU()[j].memoryError);
                }
                else {
                    // theGui.console.tPrint((
                    //     "<span color=\"#FF88FF\">Warning: device memory error ({}).</span>\r\n"), 
                    //     theSim.sCPU()[j].memoryError);
                }
            }
            if (error) break;
        }
        if (theSim.base().cancellationCalled) {
            // theGui.console.tPrint((
            //     "<span color=\"#FF88FF\">"
            //     "Warning: series cancelled, stopping\n"
            //     "after {} simulations.</span>\n"), j + 1);
            break;
        }
        // theGui.requestSliderMove(j);
        // independentPlotQueue();
    }

    //if (secondQueueThread.joinable()) secondQueueThread.join();
    auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
    if (error == 13) {
        // theGui.console.tPrint(
        //     "<span color=\"#FF88FF\">"
        //     "NaN detected in grid!\n"
        //     "Try using a larger spatial/temporal step\n"
        //     "or smaller propagation step.\n"
        //     "Simulation was cancelled.\n</span>");
    }
    else if (error == 15) {
        // theGui.console.tPrint(
        //     "<span color=\"#FF88FF\">"
        //     "Sorry, that sequence mode has been \n"
        //     "replaced by the new one. Look in the \n"
        //     "documentation for more info. It is a lot\n"
        //     "easier to use now, and hopefully \n"
        //     "it won't take long to set it up. \n"
        //     "Sorry about that!\n</span>");
    }
    else if(!error){
        // theGui.console.tPrint(
        //     "<span color=\"#88FFFF\">Finished after {:.4} s. </span>\n", 1e-6 *
        //     (double)(std::chrono::duration_cast<std::chrono::microseconds>
        //         (simulationTimerEnd - simulationTimerBegin).count()));
    }
    theSim.base().isRunning = false;
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