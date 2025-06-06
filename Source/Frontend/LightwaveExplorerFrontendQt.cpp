#include "LightwaveExplorerFrontendQt.h"
#include "cairo.h"
#include <mutex>
#include <regex>
#include <shared_mutex>
void blackoutCairoPlot(cairo_t* cr, const int width, const int height, const bool show_rendering=false){
    LweColor black(0, 0, 0, 0);
    cairo_rectangle(cr, 0, 0, width, height);
    black.setCairo(cr);
    cairo_fill(cr);
    if(show_rendering){
        std::string message("Rendering...");
        cairo_text_extents_t te;
        cairo_text_extents(cr, message.c_str(), &te);
        LweColor white(0.5,0.5,0.5,1.0);
        cairo_move_to(cr, 0.5*width - 0.5*te.x_advance, 0.5*height + 0.5*te.height);
        white.setCairo(cr);
        cairo_show_text(cr, message.c_str());
    }
}
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
        qreal ratio = devicePixelRatio();
        QImage image(ratio * size(), QImage::Format_ARGB32);
        image.setDevicePixelRatio(ratio);
        cairo_surface_t* surface = cairo_image_surface_create_for_data(
            image.bits(),
            CAIRO_FORMAT_ARGB32,
            image.width(),
            image.height(),
            image.bytesPerLine());
        cairo_t* cr = cairo_create(surface);
        cairo_scale(cr, ratio, ratio);
        theFunction(cr,
            ceil(image.width()/ratio),
            ceil(image.height()/ratio),
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

public slots:
    void passString(std::string s){
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
    void passBeamViewRender(){
        emit requestBeamViewRender();
    }
    void passBeamViewUpdate(){
        emit requestBeamViewUpdate();
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
    void requestBeamViewRender();
    void requestBeamViewUpdate();
    void moveSlider(int value);
    void sendProgressRange(int min, int max);
    void sendProgressValue(int value);
};

class SequenceValidator : public QSyntaxHighlighter {
    Q_OBJECT

    struct HighlightingRule{
        QRegularExpression pattern;
        QTextCharFormat format;
    };

    struct ColorSet{
        QColor functionColor;
        QColor numberColor;
        QColor variableColor;
        QColor openParenColor;
        QVector<QColor> matchedParensColors;
        QColor stringColor;
        QColor commentColor;
    };
    const ColorSet darkColors{
        QColor(200,200,255,255),
        QColor(255,128,255,255),
        QColor(128,255,255,255),
        QColor(255,0,0,255),
        {
            QColor(255,0,255,255),
            QColor(0,255,255,255),
            QColor(255,128,0,255)
        },
        QColor(192,128,0,255),
        QColor(64,255,64,255)
    };

    const ColorSet lightColors{
        QColor(64,0,192,255),
        QColor(0,128,0,255),
        QColor(0,64,128,255),
        QColor(255,0,0,255),
        {
            QColor(255,0,255,255),
            QColor(0,200,200,255),
            QColor(255,128,0,255)
        },
        QColor(128,64,0,255),
        QColor(0,192,0,255)
    };

    QList<HighlightingRule> highlightingRulesDark;
    QList<HighlightingRule> highlightingRulesLight;
    QList<HighlightingRule> generateRules(const ColorSet& c){
        QList<HighlightingRule> highlightingRules;
        HighlightingRule rule;
        QTextCharFormat functionFormat;

        //color known function names
        functionFormat.setForeground(c.functionColor);
        functionFormat.setFontWeight(QFont::Bold);
        const QString functionNames[] = {
            "for",
            "plasma",
            "nonlinear",
            "linear",
            "sphericalMirror",
            "parabolicMirror",
            "init",
            "default",
            "rotate",
            "set",
            "plasmaReinject",
            "save",
            "savePlasma",
            "fresnelLoss",
            "aperture",
            "farFieldAperture",
            "farFieldInverseAperture",
            "energy",
            "filter",
            "lorentzian",
            "addPulse",
            "fdtd",
            "fdtdGrid",
            "fdtdGridNearField",
            "applyOptic",
            "applyOpticX",
            "applyOpticY",
            "loadOptic",
            "rotateIntoBiaxial",
            "rotateFromBiaxial",
            "fdtdReflection",
            "polarizer"
        };
        for(const QString &fun : functionNames){
            rule.pattern = QRegularExpression(QString("\\b")+fun+QString("\\b"));
            rule.format = functionFormat;
            highlightingRules.append(rule);
        }

        //color number literals
        QTextCharFormat numberFormat;
        numberFormat.setForeground(c.numberColor);
        rule.pattern = QRegularExpression("[+-]?(\\d+(\\.\\d*)?|\\.\\d+)([eE][+-]?\\d+)?");
        rule.format = numberFormat;
        highlightingRules.append(rule);

        //color the built-in variables (d, iXX, vXX)
        QTextCharFormat variableFormat;
        variableFormat.setForeground(c.variableColor);
        rule.pattern = QRegularExpression("(?<=\\(|\\s|,)d(?=\\s|,|\\))");
        rule.format = variableFormat;
        highlightingRules.append(rule);
        rule.pattern = QRegularExpression("(?<=[\\(,\\s\\+\\*\\-\\/\\^])v\\d{2}(?=[\\s,\\)\\+\\*\\-\\/\\^])");
        highlightingRules.append(rule);
        rule.pattern = QRegularExpression("(?<=[\\(,\\s\\+\\*\\-\\/\\^])i\\d{2}(?=[\\s,\\)\\+\\*\\-\\/\\^])");
        highlightingRules.append(rule);

        //color strings
        QTextCharFormat stringFormat;
        stringFormat.setForeground(c.stringColor);
        rule.pattern = QRegularExpression("\".*\"");
        rule.format = stringFormat;
        highlightingRules.append(rule);

        //color comments
        QTextCharFormat commentFormat;
        commentFormat.setForeground(c.commentColor);
        rule.pattern = QRegularExpression("<.*>");
        rule.format = commentFormat;
        highlightingRules.append(rule);

        return highlightingRules;
    }
public:
    SequenceValidator(QTextDocument* parent = nullptr) : QSyntaxHighlighter(parent){
        highlightingRulesDark = generateRules(
            darkColors
        );
        highlightingRulesLight = generateRules(
            lightColors
        );
    }
protected:
    void highlightBlock(const QString& text) override {
    #ifdef __linux__
        const bool isLightTheme = (QApplication::palette().windowText().color().lightnessF()<0.5);
    #else
        const bool isLightTheme = (QApplication::styleHints()->colorScheme() == Qt::ColorScheme::Light);
    #endif
        QList<HighlightingRule>& highlightingRules = isLightTheme ?
            highlightingRulesLight : highlightingRulesDark;
        const ColorSet& c = isLightTheme ? lightColors : darkColors;
        for (const HighlightingRule &rule : std::as_const(highlightingRules)) {
            QRegularExpressionMatchIterator matchIterator = rule.pattern.globalMatch(text);
            while (matchIterator.hasNext()) {
                QRegularExpressionMatch match = matchIterator.next();
                setFormat(match.capturedStart(), match.capturedLength(), rule.format);
            }
        }

        // Set the colors for matched and unmatched parentheses
        QStack<int> parenthesesStack;
        QTextCharFormat badParenFormat;
        badParenFormat.setForeground(c.openParenColor);
        badParenFormat.setFontWeight(QFont::Bold);
        for (qsizetype i = 0; i < text.length(); ++i) {
            QChar currentChar = text[i];
            if (currentChar == '(') {
                parenthesesStack.push(i);
            } else if (currentChar == ')') {
                if (!parenthesesStack.isEmpty()) {
                    int openingIndex = parenthesesStack.pop();
                    int colorIndex = openingIndex % c.matchedParensColors.size();
                    setFormat(openingIndex, 1, c.matchedParensColors[colorIndex]);
                    setFormat(i, 1, c.matchedParensColors[colorIndex]);
                } else {
                    setFormat(i, 1, badParenFormat);
                }
            }
        }
        while (!parenthesesStack.isEmpty()) {
            int openingIndex = parenthesesStack.pop();
            setFormat(openingIndex, 1, badParenFormat);
        }
    }
private:
};


class LWEGui : public QMainWindow {
    Q_OBJECT
    bool isInputRegionHidden=false;
    QThread* messengerThread;
    void populateDatabasePulldown(){
        std::string materialString;
        pulldowns["material"]->clear();
        for (std::size_t i = 0; i < theDatabase.db.size(); ++i) {
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
    QWidget* inputRegion;
    QWidget* plotControlRegion;
    QSlider* slider;
    LweImage<double> beamView;
    std::shared_mutex beamViewMutex;
    GuiMessenger* messenger;
    std::shared_mutex m;

    bool isMakingSVG = false;
    std::array<std::string,9> SVGStrings;
    template<typename... Args> void cPrint(std::string_view format, Args&&... args) {
        std::string s = Svformat(format, Smake_format_args(args...));
        console->append(s.c_str());
    }
    template<typename... Args> void sPrint(std::string_view format, Args&&... args) {
        std::string s = Svformat(format, Smake_format_args(args...));
        sequence->insertPlainText(s.c_str());
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

        enum class SecondValueDefault {Default_zero = 0, Default_same = 1};
        auto setTwoDoublesIfThereIsASemicolon = [](QLineEdit* box, double& value1, double& value2, SecondValueDefault defaultType, double multiplier = 1.0){
            QString input = box -> text();
            int semicolonIndex = input.indexOf(';');
            if(semicolonIndex == -1){
                value1 = multiplier * input.toDouble();
                switch(defaultType){
                    case SecondValueDefault::Default_zero:
                        value2 = 0.0;
                        break;
                    case SecondValueDefault::Default_same:
                        value2 = value1;
                        break;
                }
            }
            else{
                value1 = multiplier * input.left(semicolonIndex).toDouble();
                value2 = multiplier * input.mid(semicolonIndex+1).toDouble();
            }
        };

        setToDouble(textBoxes["energy1"],sim.base().pulse1.energy);
        setToDoubleMultiplier(textBoxes["frequency1"],1e-12,sim.base().pulse1.frequency);
        setToDoubleMultiplier(textBoxes["bandwidth1"],1e-12,sim.base().pulse1.bandwidth);
        setToInt(textBoxes["SGOrder1"],sim.base().pulse1.sgOrder);
        setToDoubleMultiplier(textBoxes["CEP1"],1.0/vPi<double>(),sim.base().pulse1.cep);
        setToDoubleMultiplier(textBoxes["delay1"],1e15,sim.base().pulse1.delay);
        setToDoubleMultiplier(textBoxes["GDD1"],1e30,sim.base().pulse1.gdd);
        setToDoubleMultiplier(textBoxes["TOD1"],1e45,sim.base().pulse1.tod);
        setToInt(textBoxes["material1"],sim.base().pulse1.phaseMaterial);
        setToDoubleMultiplier(textBoxes["thickness1"],1e6,sim.base().pulse1.phaseMaterialThickness);
        setToDoubleMultiplier(textBoxes["beamwaist1"],1e6,sim.base().pulse1.beamwaist);
        setTwoDoublesIfThereIsASemicolon(textBoxes["xOffset1"], sim.base().pulse1.x0, sim.base().pulse1.y0, SecondValueDefault::Default_zero, 1e-6);
        setToDoubleMultiplier(textBoxes["zOffset1"],1e6,sim.base().pulse1.z0);
        setTwoDoublesIfThereIsASemicolon(textBoxes["NCAngle1"], sim.base().pulse1.beamAngle, sim.base().pulse1.beamAnglePhi, SecondValueDefault::Default_zero, deg2Rad<double>());
        setToDoubleMultiplier(textBoxes["polarization1"],rad2Deg<double>(),sim.base().pulse1.polarizationAngle);
        setToDouble(textBoxes["circularity1"],sim.base().pulse1.circularity);

        setToDouble(textBoxes["energy2"],sim.base().pulse2.energy);
        setToDoubleMultiplier(textBoxes["frequency2"],1e-12,sim.base().pulse2.frequency);
        setToDoubleMultiplier(textBoxes["bandwidth2"],1e-12,sim.base().pulse2.bandwidth);
        setToInt(textBoxes["SGOrder2"],sim.base().pulse2.sgOrder);
        setToDoubleMultiplier(textBoxes["CEP2"],1.0/vPi<double>(),sim.base().pulse2.cep);
        setToDoubleMultiplier(textBoxes["delay2"],1e15,sim.base().pulse2.delay);
        setToDoubleMultiplier(textBoxes["GDD2"],1e30,sim.base().pulse2.gdd);
        setToDoubleMultiplier(textBoxes["TOD2"],1e45,sim.base().pulse2.tod);
        setToInt(textBoxes["material2"],sim.base().pulse2.phaseMaterial);
        setToDoubleMultiplier(textBoxes["thickness2"],1e6,sim.base().pulse2.phaseMaterialThickness);
        setToDoubleMultiplier(textBoxes["beamwaist2"],1e6,sim.base().pulse2.beamwaist);
        setTwoDoublesIfThereIsASemicolon(textBoxes["xOffset2"], sim.base().pulse2.x0, sim.base().pulse2.y0, SecondValueDefault::Default_zero, 1e-6);
        setToDoubleMultiplier(textBoxes["zOffset2"],1e6,sim.base().pulse2.z0);
        setTwoDoublesIfThereIsASemicolon(textBoxes["NCAngle2"], sim.base().pulse2.beamAngle, sim.base().pulse2.beamAnglePhi, SecondValueDefault::Default_zero, deg2Rad<double>());
        setToDoubleMultiplier(textBoxes["polarization2"],rad2Deg<double>(),sim.base().pulse2.polarizationAngle);
        setToDouble(textBoxes["circularity2"],sim.base().pulse2.circularity);

        sim.base().pulse1FileType = pulldowns["pulse1type"]->currentIndex();
        sim.base().pulse2FileType = pulldowns["pulse2type"]->currentIndex();
        sim.base().fittingMode = pulldowns["fit"]->currentIndex();
        sim.base().materialIndex = pulldowns["material"]->currentIndex();
        setToDoubleMultiplier(textBoxes["crystalTheta"],rad2Deg<double>(),sim.base().crystalTheta);
        setToDoubleMultiplier(textBoxes["crystalPhi"],rad2Deg<double>(),sim.base().crystalPhi);
        setTwoDoublesIfThereIsASemicolon(textBoxes["NLAbsorption"],sim.base().nonlinearAbsorptionStrength, sim.base().startingCarrierDensity, SecondValueDefault::Default_zero);
        setToDouble(textBoxes["crystalBandgap"],sim.base().bandGapElectronVolts);
        setToDoubleMultiplier(textBoxes["DrudeGamma"],1e-12,sim.base().drudeGamma);
        setToDouble(textBoxes["effectiveMass"],sim.base().effectiveMass);
        setTwoDoublesIfThereIsASemicolon(textBoxes["xSize"], sim.base().spatialWidth, sim.base().spatialHeight, SecondValueDefault::Default_same, 1e-6);
        setToDoubleMultiplier(textBoxes["dx"],1e6,sim.base().rStep);
        setToDoubleMultiplier(textBoxes["timeSpan"],1e15,sim.base().timeSpan);
        setToDoubleMultiplier(textBoxes["dt"],1e15,sim.base().tStep);
        setToDoubleMultiplier(textBoxes["zSize"],1e6,sim.base().crystalThickness);
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
        sim.base().fittingLoadedData = fittingLoadedData;

        //derived parameters and cleanup:
        sim.base().sellmeierType = 0;
        sim.base().axesNumber = 0;
        sim.base().Ntime = (int64_t)(minGridDimension * round(sim.base().timeSpan
            / (minGridDimension * sim.base().tStep)));
        sim.base().timeSpan = sim.base().tStep * sim.base().Ntime;
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
            theDatabase.db[sim.base().materialIndex].nonlinearSwitches;
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

        enum class NoSemicolonIf{zero = 0, same = 1, same_or_zero = 2};
        auto setToTwoDoubles = [](QLineEdit* box, const double value1, const double value2, NoSemicolonIf oneValueCondition){
            bool setBoth = true;
            switch(oneValueCondition){
                case NoSemicolonIf::zero:
                    if(value2==0.0){
                        setBoth = false;
                    }
                    break;
                case NoSemicolonIf::same:
                    if(value1==value2){
                        setBoth = false;
                    }
                    break;
                case NoSemicolonIf::same_or_zero:
                    if(value1==value2 || value2==0.0){
                        setBoth = false;
                    }
            }
            if(setBoth){
                QString s(Sformat(std::string_view("{:g}; {:g}"), value1, value2).c_str());
                box->setText(s);
            }
            else{
                QString s(Sformat(std::string_view("{:g}"), value1).c_str());
                box->setText(s);
            }
        };

        setToDouble(textBoxes["energy1"],sim.base().pulse1.energy);
        setToDouble(textBoxes["frequency1"],1e-12*sim.base().pulse1.frequency);
        setToDouble(textBoxes["bandwidth1"],1e-12*sim.base().pulse1.bandwidth);
        setToDouble(textBoxes["SGOrder1"],sim.base().pulse1.sgOrder);
        setToDouble(textBoxes["CEP1"],sim.base().pulse1.cep / vPi<double>());
        setToDouble(textBoxes["delay1"],1e15*sim.base().pulse1.delay);
        setToDouble(textBoxes["GDD1"],1e30*sim.base().pulse1.gdd);
        setToDouble(textBoxes["TOD1"],1e45*sim.base().pulse1.tod);
        setToInt(textBoxes["material1"],sim.base().pulse1.phaseMaterial);
        setToDouble(textBoxes["thickness1"],1e6*sim.base().pulse1.phaseMaterialThickness);
        setToDouble(textBoxes["beamwaist1"],1e6*sim.base().pulse1.beamwaist);
        setToTwoDoubles(textBoxes["xOffset1"], 1e6*sim.base().pulse1.x0, 1e6*sim.base().pulse1.y0, NoSemicolonIf::zero);
        setToDouble(textBoxes["zOffset1"],1e6*sim.base().pulse1.z0);
        setToTwoDoubles(textBoxes["NCAngle1"],rad2Deg<double>()*sim.base().pulse1.beamAngle, rad2Deg<double>()*sim.base().pulse1.beamAnglePhi, NoSemicolonIf::zero);
        setToDouble(textBoxes["polarization1"],rad2Deg<double>()*sim.base().pulse1.polarizationAngle);
        setToDouble(textBoxes["circularity1"],sim.base().pulse1.circularity);

        setToDouble(textBoxes["energy2"],sim.base().pulse2.energy);
        setToDouble(textBoxes["frequency2"],1e-12*sim.base().pulse2.frequency);
        setToDouble(textBoxes["bandwidth2"],1e-12*sim.base().pulse2.bandwidth);
        setToDouble(textBoxes["SGOrder2"],sim.base().pulse2.sgOrder);
        setToDouble(textBoxes["CEP2"],sim.base().pulse2.cep / vPi<double>());
        setToDouble(textBoxes["delay2"],1e15*sim.base().pulse2.delay);
        setToDouble(textBoxes["GDD2"],1e30*sim.base().pulse2.gdd);
        setToDouble(textBoxes["TOD2"],1e45*sim.base().pulse2.tod);
        setToInt(textBoxes["material2"],sim.base().pulse2.phaseMaterial);
        setToDouble(textBoxes["thickness2"],1e6*sim.base().pulse2.phaseMaterialThickness);
        setToDouble(textBoxes["beamwaist2"],1e6*sim.base().pulse2.beamwaist);
        setToTwoDoubles(textBoxes["xOffset2"], 1e6*sim.base().pulse2.x0, 1e6*sim.base().pulse2.y0, NoSemicolonIf::zero);
        setToDouble(textBoxes["zOffset2"],1e6*sim.base().pulse2.z0);
        setToTwoDoubles(textBoxes["NCAngle2"],rad2Deg<double>()*sim.base().pulse2.beamAngle, rad2Deg<double>()*sim.base().pulse2.beamAnglePhi, NoSemicolonIf::zero);
        setToDouble(textBoxes["polarization2"],rad2Deg<double>()*sim.base().pulse2.polarizationAngle);
        setToDouble(textBoxes["circularity2"],sim.base().pulse2.circularity);

        pulldowns["material"]->setCurrentIndex(sim.base().materialIndex);
        pulldowns["batch1"]->setCurrentIndex(sim.base().batchIndex);
        pulldowns["batch2"]->setCurrentIndex(sim.base().batchIndex2);
        pulldowns["propagator"]->setCurrentIndex(sim.base().symmetryType);
        pulldowns["pulse1type"]->setCurrentIndex(sim.base().pulse1FileType);
        pulldowns["pulse2type"]->setCurrentIndex(sim.base().pulse2FileType);

        setToDouble(textBoxes["crystalTheta"],rad2Deg<double>() * sim.base().crystalTheta);
        setToDouble(textBoxes["crystalPhi"],rad2Deg<double>() *sim.base().crystalPhi);
        setToTwoDoubles(textBoxes["NLAbsorption"],sim.base().nonlinearAbsorptionStrength, sim.base().startingCarrierDensity, NoSemicolonIf::zero);
        setToDouble(textBoxes["crystalBandgap"],sim.base().bandGapElectronVolts);
        setToDouble(textBoxes["DrudeGamma"],1e-12*sim.base().drudeGamma);
        setToDouble(textBoxes["effectiveMass"],sim.base().effectiveMass);
        setToTwoDoubles(textBoxes["xSize"],1e6*sim.base().spatialWidth, 1e6*sim.base().spatialHeight, NoSemicolonIf::same_or_zero);
        setToDouble(textBoxes["dx"],1e6*sim.base().rStep);
        setToDouble(textBoxes["timeSpan"],1e15*sim.base().timeSpan);
        setToDouble(textBoxes["dt"],1e15*sim.base().tStep);
        setToDouble(textBoxes["zSize"],1e6*sim.base().crystalThickness);
        setToDouble(textBoxes["dz"],1e9*sim.base().propagationStep);

        setToDouble(textBoxes["batch1end"],sim.base().batchDestination);
        setToDouble(textBoxes["batch2end"],sim.base().batchDestination2);
        setToInt(textBoxes["batch1steps"],sim.base().Nsims);
        setToInt(textBoxes["batch2steps"],sim.base().Nsims2);

        std::string formattedFit=sim.base().fittingString;
        insertAfterCharacter(formattedFit,';',std::string("\n"));
        if(formattedFit.size() > 4 && formattedFit.substr(0,4) != "None"){
            fitting->setText(QString::fromStdString(formattedFit));
        }
        else{
            fitting->setText("");
        }

        std::string formattedSequence = sim.base().sequenceString;
        formatSequence(formattedSequence);
        if(formattedSequence.size() > 4 && formattedSequence.substr(4) != "None"){
            sequence->setText(QString::fromStdString(formattedSequence));
        }
        else{
            sequence->setText("");
        }
    }
    LWEGui(){
#if defined(Q_OS_WIN)
        const int textBoxWidth = 78;
        const int textBoxHeight = 26;
        const int miniButtonWidth = 30;
        const int mainButtonHeight = textBoxHeight;
#elif defined(Q_OS_MAC)
        const int textBoxWidth = 80;
        const int textBoxHeight = 26;
        const int miniButtonWidth = 30;
        const int mainButtonHeight = textBoxHeight+6;
#elif defined(Q_OS_LINUX)
        const int textBoxWidth = 78;
        const int textBoxHeight = 26;
        const int miniButtonWidth = 30;
        const int mainButtonHeight = textBoxHeight;
#else
#endif
        const int plotTextBoxWidth = 40;
        const int plotLabelWidth = 80;
        const int labelWidth = 2*textBoxWidth;
        const int rowHeight = textBoxHeight+2;
        const int rowWidth = labelWidth + 2*textBoxWidth + 10;
        const int mainButtonWidth = rowWidth/4;

        const int pulldownContainerWidth = labelWidth+4;
        QFont emojiFont = getEmojiFont();

        //Divide the main window into a large expanding upper panel and a control strip at the bottom
        auto squeezeMargins = [&](QBoxLayout* layout){
            layout->setSpacing(0);
            layout->setContentsMargins(0,0,0,0);
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
        //add a plot control region to show when the input area is collapsed
        plotControlRegion = new QWidget;
        QHBoxLayout *mainAreaLayout = new QHBoxLayout(mainArea);
        squeezeMargins(mainAreaLayout);
        inputRegion->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
        plotRegion->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotControlRegion->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        mainAreaLayout->addWidget(inputRegion);
        mainAreaLayout->addWidget(plotRegion);
        mainAreaLayout->addWidget(plotControlRegion);
        plotControlRegion->setFixedWidth(220);
        plotControlRegion->hide();
        QGridLayout* plotRegionLayout = new QGridLayout(plotRegion);
        plotRegionLayout->setContentsMargins(0,0,0,0);
        plotRegionLayout->setSpacing(2);

        QVBoxLayout* plotControlRegionLayout = new QVBoxLayout(plotControlRegion);

        plots["timeImage1"] = new CairoWidget(*this, drawTimeImage1);
        plots["timeImage1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plots["timeImage1"]->setToolTip(
            "Image of the electric field grid: presents a slice of Ex(x,y=0,t),\n"
            "there the horizontal axis is time, and the vertical axis is position");
        plotRegionLayout->addWidget(plots["timeImage1"],0,0);
        plots["timeImage2"] = new CairoWidget(*this, drawTimeImage2);
        plots["timeImage2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plots["timeImage2"]->setToolTip("Image of the electric field grid: presents a slice of Ey(x,y=0,t),\n"
            "there the horizontal axis is time, and the vertical axis is position");
        plotRegionLayout->addWidget(plots["timeImage2"],1,0);
        plots["freqImage1"] = new CairoWidget(*this, drawFourierImage1);
        plots["freqImage1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plots["freqImage1"]->setToolTip("Plot of the electric field grid in momentum-frequency space: Ex(kx,ky=0,f).\n"
            "Is plotted on a logarithmic scale. Vertical axis is transverse momentum kx,\n"
            "and horizontal axis is frequency f.");
        plotRegionLayout->addWidget(plots["freqImage1"],0,1);
        plots["freqImage2"] = new CairoWidget(*this, drawFourierImage2);
        plots["freqImage2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["freqImage2"],1,1);
        plots["freqImage2"]->setToolTip("Plot of the electric field grid in momentum-frequency space:\n"
            "Ey(kx,ky=0,f). Is plotted on a logarithmic scale. Vertical axis is transverse momentum\n"
            "kx, and horizontal axis is frequency f.");

        plots["beamView"] = new CairoWidget(*this, drawBeamImage);
        plots["beamView"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plotRegionLayout->addWidget(plots["beamView"],2,0,1,2);
        plots["beamView"]->setToolTip("Face-on view of the output light. Depending on the setting in the\n"
        "plot panel (right-hand side when the toggle arrow has been pressed)\n"
        "either near-field at the end of the simulation (position) or far-field (angle)");
        //plots["beamView"]->hide();

        plots["timePlot1"] = new CairoWidget(*this, drawField1Plot);
        plots["timePlot1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plots["timePlot1"]->setToolTip(
            "Plot of the on-axis electric field in the x-polarization. When the polarizations are combined,\n"
        "it will show both x-polarization (cyan) and y-polarization (magenta).");
        plotRegionLayout->addWidget(plots["timePlot1"],4,0);

        plots["timePlot2"] = new CairoWidget(*this, drawField2Plot);
        plots["timePlot2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plots["timePlot2"]->setToolTip(
            "Plot of the on-axis electric field in the y-polarization");
        plotRegionLayout->addWidget(plots["timePlot2"],5,0);
        plots["freqPlot1"] = new CairoWidget(*this, drawSpectrum1Plot);
        plots["freqPlot1"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plots["freqPlot1"]->setToolTip("Plot of the energy spectrum of the result\n x-polarization (purple)\n y-polarization (red)\n total (orange)\n");
        plotRegionLayout->addWidget(plots["freqPlot1"],4,1);
        plots["freqPlot2"] = new CairoWidget(*this, drawSpectrum2Plot);
        plots["freqPlot2"]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        plots["freqPlot2"]->setToolTip(
            "Plot of the energy spectrum of the result, y-polarization.");
        plotRegionLayout->addWidget(plots["freqPlot2"],5,1);

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
        plotControlStrip->setFixedHeight(rowHeight);
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
            QVBoxLayout* location,
            const QString& tooltip){
                QHBoxLayout* rowLayout = getRowBoxLayout(location);
                labels[entry1] = new QLabel;
                labels[entry1]->setToolTip(tooltip);
                textBoxes[entry1] = new QLineEdit;
                textBoxes[entry2] = new QLineEdit;
                textBoxes[entry1]->setToolTip(tooltip);
                textBoxes[entry2]->setToolTip(tooltip);
                labels[entry1]->setToolTip(tooltip);
                labels[entry1]->setText(label);
                labels[entry1]->setFixedSize(labelWidth,textBoxHeight);
                textBoxes[entry1]->setFixedSize(textBoxWidth,textBoxHeight);
                textBoxes[entry2]->setFixedSize(textBoxWidth,textBoxHeight);
                rowLayout->addWidget(labels[entry1]);
                rowLayout->addWidget(textBoxes[entry1]);
                rowLayout->addWidget(textBoxes[entry2]);
        };

        auto addTextBoxRowForPlots = [&](
            const QString& label,
            const std::string& entry1,
            const std::string& entry2,
            const std::string& entry3,
            QVBoxLayout* location,
            const QString& tooltip,
            const float value1,
            const float value2,
            const float value3){
                QHBoxLayout* rowLayout = getRowBoxLayout(location);
                rowLayout->setAlignment(Qt::AlignLeft);
                rowLayout->setSpacing(2);
                labels[entry1] = new QLabel;
                labels[entry1]->setToolTip(tooltip);
                textBoxes[entry1] = new QLineEdit;
                if(entry2 != "none")textBoxes[entry2] = new QLineEdit;
                if(entry3 != "none")textBoxes[entry3] = new QLineEdit;
                textBoxes[entry1]->setToolTip(tooltip);
                if(entry2 != "none")textBoxes[entry2]->setToolTip(tooltip);
                if(entry3 != "none")textBoxes[entry3]->setToolTip(tooltip);
                labels[entry1]->setToolTip(tooltip);
                labels[entry1]->setText(label);
                labels[entry1]->setFixedSize(plotLabelWidth,textBoxHeight);
                textBoxes[entry1]->setFixedSize(plotTextBoxWidth,textBoxHeight);
                if(entry2 != "none")textBoxes[entry2]->setFixedSize(plotTextBoxWidth,textBoxHeight);
                if(entry3 != "none")textBoxes[entry3]->setFixedSize(plotTextBoxWidth,textBoxHeight);
                QString s1(Sformat(std::string_view("{:g}"), value1).c_str());
                textBoxes[entry1]->setText(s1);
                if(entry2 != "none"){
                    QString s2(Sformat(std::string_view("{:g}"), value2).c_str());
                    textBoxes[entry2]->setText(s2);
                }

                if(entry3 != "none"){
                    QString s3(Sformat(std::string_view("{:g}"), value3).c_str());
                    textBoxes[entry3]->setText(s3);
                }
                rowLayout->addWidget(labels[entry1], Qt::AlignLeft);
                rowLayout->addWidget(textBoxes[entry1], Qt::AlignLeft);
                if(entry2 != "none")rowLayout->addWidget(textBoxes[entry2], Qt::AlignLeft);
                if(entry3 != "none")rowLayout->addWidget(textBoxes[entry3], Qt::AlignLeft);
        };

        auto addPulldownInContainer = [&](int width, QBoxLayout* location, const std::string& entry){
            QWidget* container = new QWidget;
            QHBoxLayout* fitContainerLayout = new QHBoxLayout(container);
            container->setFixedSize(width,textBoxHeight);
            fitContainerLayout->setContentsMargins(0,0,0,0);
            pulldowns[entry] = new QComboBox;
            pulldowns[entry]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
            location->addWidget(container);
            fitContainerLayout->addWidget(pulldowns[entry]);
        };


        //First column
        addTextBoxRow("Pulse energy (J)", "energy1", "energy2", entryColumn1Layout,
        "The energy contained in the pulses. \nNote that there are two columns, the first controls the primary pulse.\n"
        " For 2D modes, the intensity will be set as it would be\nfor a round Gaussian beam with that size");
        addTextBoxRow("Frequency (THz)", "frequency1", "frequency2", entryColumn1Layout,
        "The central frequency of the pulse. Note that the simulation in the UPPE\n"
        "modes uses a moving frame centered around the group velocity of the\n"
        "pulse in the first column.");
        addTextBoxRow("Bandwidth (THz)", "bandwidth1", "bandwidth2", entryColumn1Layout,
        "The bandwidth of the pulses.");
        addTextBoxRow("SG order", "SGOrder1", "SGOrder2", entryColumn1Layout,
        "The order of the super-Gaussian spectrum. 2 implies a standard Gaussian, and the spectrum\n"
        " will be increasingly flattop as the order increases. Higher order will lead to\n"
        "additional structure outside of the main pulse in the time domain.");
        addTextBoxRow("CEP/\xcf\x80", "CEP1", "CEP2", entryColumn1Layout,
        "The carrier-envelope phase of the pulse, in units of pi radians.\n"
        "Setting this to 1 will flip the field upside-down.");
        addTextBoxRow("Delay (fs)", "delay1", "delay2", entryColumn1Layout,
        "Temporal delay of the pulses.");
        addTextBoxRow("GDD (fs\xc2\xb2)", "GDD1", "GDD2", entryColumn1Layout,
        "Group-delay dispersion of the pulses, i.e. linear chirp.");
        addTextBoxRow("TOD (fs\xc2\xb3)", "TOD1", "TOD2", entryColumn1Layout,
        "Third-order dispersion of the pulses.");
        addTextBoxRow("Phase material", "material1", "material2", entryColumn1Layout,
        "Select a material number - the numbers come from the Material pulldown menu -\n"
        "to apply to the pulse before it enters the simulation. For example, you can\n"
        "add a given amount of fused silica transmission to the pulse to simulate\n"
        " the effect of scanning a wedge pair in the beam.");
        addTextBoxRow("Thickness (\xce\xbcm)", "thickness1", "thickness2", entryColumn1Layout,
        "This sets the thickness of the linear propagation through the material selected above.");
        addTextBoxRow("Beamwaist (\xce\xbcm)", "beamwaist1", "beamwaist2", entryColumn1Layout,
        "Gaussian beamwaist of the input beam in space. In terms of intensity, this is the 1/e^2 radius");
        addTextBoxRow("x offset (\xce\xbcm)", "xOffset1", "xOffset2", entryColumn1Layout,
        "Offset of the beam in the x-direction (up and down on the screen, transverse to propagation).\n"
        "In 3D radial symmetry mode, this does nothing.");
        addTextBoxRow("z offset (\xce\xbcm)", "zOffset1", "zOffset2", entryColumn1Layout,
        "Offset of the beam in the z-direction (the direction of propagation).");
        addTextBoxRow("NC angle (deg)", "NCAngle1", "NCAngle2", entryColumn1Layout,
        "Angle of the beam with respect to the z-axis (in the x-z plane), i.e. the noncollinear angle.");
        addTextBoxRow("Polarization (deg)", "polarization1", "polarization2", entryColumn1Layout,
        "Polarization angle of the beam. With 0: field points along x-axis, 90: along y-axis.");
        addTextBoxRow("Circularity", "circularity1", "circularity2", entryColumn1Layout,
        "Degree of circularity of the beam. 0: Linear polarization, 1: circular (can be any value in between).");

        QHBoxLayout* pulseTypeRow = getRowBoxLayout(entryColumn1Layout);
        labels["source"] = new QLabel;
        labels["source"]->setText("Source");
        labels["source"]->setFixedWidth(labelWidth);
        pulseTypeRow->addWidget(labels["source"]);
        addPulldownInContainer(textBoxWidth,pulseTypeRow,"pulse1type");
        pulldowns["pulse1type"]->addItem("Synthetic");
        pulldowns["pulse1type"]->addItem("FROG");
        pulldowns["pulse1type"]->addItem("Wave");
        pulldowns["pulse1type"]->addItem("LWE");
        pulldowns["pulse1type"]->setToolTip(
            "Determine whether the pulse is defined by the values set above\n"
            "or using an input file, either:\n"
            "a FROG .speck file\n"
            "a waveform in ASCII format (time|normalized amplitude)\n"
            "A previous LWE result field (a .bin from a result archive)\n");
        addPulldownInContainer(textBoxWidth,pulseTypeRow,"pulse2type");
        pulldowns["pulse2type"]->addItem("Synthetic");
        pulldowns["pulse2type"]->addItem("FROG");
        pulldowns["pulse2type"]->addItem("Wave");
        pulldowns["pulse2type"]->addItem("LWE");
        pulldowns["pulse2type"]->setToolTip(pulldowns["pulse1type"]->toolTip());
        labels["source"]->setToolTip(pulldowns["pulse1type"]->toolTip());
        QHBoxLayout* loadPulseRow = getRowBoxLayout(entryColumn1Layout);
        labels["sourceFile"] = new QLabel;
        labels["sourceFile"]->setText("Source file");
        labels["sourceFile"]->setFixedWidth(labelWidth);
        loadPulseRow->addWidget(labels["sourceFile"]);
        buttons["loadPulse1"] = new QPushButton("\xf0\x9f\x93\x82");
        buttons["loadPulse1"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        buttons["loadPulse1"]->setFont(emojiFont);
        loadPulseRow->addWidget(buttons["loadPulse1"]);
        buttons["loadPulse2"] = new QPushButton("\xf0\x9f\x93\x82");
        buttons["loadPulse2"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        buttons["loadPulse2"]->setFont(emojiFont);
        loadPulseRow->addWidget(buttons["loadPulse2"]);
        buttons["loadPulse1"]->setToolTip(
            "Load the file for describing the electric field, if not using\n"
            "a synthetic pulse.");
        buttons["loadPulse2"]->setToolTip(buttons["loadPulse1"]->toolTip());
        QHBoxLayout* slurmRow = getRowBoxLayout(entryColumn1Layout);
        labels["slurm"] = new QLabel;
        labels["slurm"]->setText("SLURM script");
        labels["slurm"]->setFixedWidth(labelWidth);
        slurmRow->addWidget(labels["slurm"]);
        addPulldownInContainer(pulldownContainerWidth-miniButtonWidth,slurmRow,"slurm");
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
        buttons["saveSLURM"]->setFont(emojiFont);
        slurmRow->addWidget(buttons["saveSLURM"]);
        pulldowns["slurm"]->setToolTip("Generate a SLURM script package\n"
        "which can be run on a GPU cluster, such as the ones\n"
        "managed by the Max Planck Society.");
        labels["slurm"]->setToolTip(pulldowns["slurm"]->toolTip());
        buttons["saveSLURM"]->setToolTip(pulldowns["slurm"]->toolTip());
        //Second column
        QHBoxLayout* materialRow = getRowBoxLayout(entryColumn2Layout);
        labels["material"] = new QLabel;
        labels["material"]->setText("Material");
        labels["material"]->setFixedWidth(labelWidth-miniButtonWidth);
        materialRow->addWidget(labels["material"]);
        buttons["loadMaterial"] = new QPushButton("\xf0\x9f\x93\x82");
        buttons["loadMaterial"]->setFixedSize(miniButtonWidth,mainButtonHeight);
        buttons["loadMaterial"]->setFont(emojiFont);
        materialRow->addWidget(buttons["loadMaterial"]);
        addPulldownInContainer(pulldownContainerWidth,materialRow,"material");
        populateDatabasePulldown();
        pulldowns["material"]->setToolTip(
            "Select the material to used in propagation. If the name is\n"
            "followed by an L, it means the 7-Lorentzian formula was used, so\n"
            "FDTD propagation is possible. If it is followed by G, it uses Gaussian bands\n"
            "for absorption, with Dawson function shapes of the refractive index.");
        labels["material"]->setToolTip(pulldowns["material"]->toolTip());
        buttons["loadMaterial"]->setToolTip("Load a custom crystal database");

        addTextBoxRow("Theta, phi (deg)", "crystalTheta", "crystalPhi", entryColumn2Layout,
        "These angles define the orientation of the nonlinear crystal.");
        addTextBoxRow("NL absorption", "NLAbsorption", "crystalBandgap", entryColumn2Layout,
        "The first parameter is the nonlinear absorption strength inside of the crystal.\n"
        "If it's 0, then there will be no nonlinear absorption or plasma generation. Typical\n"
        "values range from 10^-19 to 10^-18. The second parameter is the bandgap of the crystal\n"
        " in eV. This, plus the frequency of the first pulse, determines the order of multiphoton\n"
        "absorption.");
        addTextBoxRow("Drude: gamma, m", "DrudeGamma", "effectiveMass", entryColumn2Layout,
        "The first parameter is the Drude model momentum relaxation rate, gamma, in units of\n"
        "THz. The second is the reduced effective mass, relative to the electron mass, i.e.\n"
        "a free electron in vacuum would have a value of 1.");
        addTextBoxRow("Max x, dx (\xce\xbcm)", "xSize", "dx", entryColumn2Layout,
        "This gives the parameters of the spatial grid in the direction transverse to beam\n"
        "propagation. The first parameter gives the size of the grid, and the second is the\n"
        "size of each pixel.");
        addTextBoxRow("Time span, dt (fs)", "timeSpan", "dt", entryColumn2Layout,
        "These parameters define the temporal grid.");
        addTextBoxRow("Max z, dz (\xce\xbcm, nm)", "zSize", "dz", entryColumn2Layout,
        "These parameters define the discretization in the z-direction. This has a different\n"
        "meaning in UPPE and FDTD modes. In UPPE, the pulse propagates forward in space along the\n"
        "z-axis - this defines how far to go, and with which propagation step size. In FDTD,\n"
        "the propagation takes place on a spatial grid, and advances in time. Max z there gives\n"
        "the crystal thickness, and dz defines how fine the grid mesh is in the z-direction.");

        QHBoxLayout* propagationRow = getRowBoxLayout(entryColumn2Layout);
        labels["propagator"] = new QLabel;
        labels["propagator"]->setText("Propagation");
        labels["propagator"]->setFixedWidth(labelWidth);
        propagationRow->addWidget(labels["propagator"]);
        addPulldownInContainer(pulldownContainerWidth,propagationRow,"propagator");
        pulldowns["propagator"]->addItem(("2D Cartesian"));
        pulldowns["propagator"]->addItem(("3D radial symmetry"));
        pulldowns["propagator"]->addItem(("3D"));
        pulldowns["propagator"]->addItem(("FDTD 2D"));
        pulldowns["propagator"]->addItem(("FDTD 3D"));
        pulldowns["propagator"]->setToolTip("Choose the propagation mode.\n"
        "The first three use the UPPE, while the last two use FDTD in cartesian\n"
        "coordinates.");
        labels["propagator"]->setToolTip(pulldowns["propagator"]->toolTip());

        QHBoxLayout* batch1Row = getRowBoxLayout(entryColumn2Layout);
        labels["batch1"] = new QLabel;
        labels["batch1"]->setText("Batch mode");
        labels["batch1"]->setFixedWidth(labelWidth);
        batch1Row->addWidget(labels["batch1"]);
        addPulldownInContainer(pulldownContainerWidth,batch1Row,"batch1");

        QHBoxLayout* batch2Row = getRowBoxLayout(entryColumn2Layout);
        labels["batch2"] = new QLabel;
        labels["batch2"]->setText("Batch 2 mode");
        labels["batch2"]->setFixedWidth(labelWidth);
        batch2Row->addWidget(labels["batch2"]);
        addPulldownInContainer(pulldownContainerWidth,batch2Row,"batch2");
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
        pulldowns["batch1"]->setToolTip(
        "The batch mode allows you to perform a series of simulations\n"
        "with different values of a given parameter.\n"
        "You can do two batches at once, creating a 2D dataset of all combinations.");
        pulldowns["batch1"]->setToolTip(pulldowns["batch1"]->toolTip());
        labels["batch1"]->setToolTip(pulldowns["batch1"]->toolTip());
        labels["batch2"]->setToolTip(pulldowns["batch1"]->toolTip());

        addTextBoxRow("Batch end", "batch1end", "batch2end", entryColumn2Layout,
        "These parameters define the targets of a batch calculation. For example, if\n"
        "The pulse has a frequency of 200 THz, and 300 is entered in this box, and\n"
        "'03: frequency 1' is selected from the Batch mode pulldown, the frequency will be\n"
        "scanned from 200 to 300 THz.");
        addTextBoxRow("Batch steps", "batch1steps", "batch2steps", entryColumn2Layout,
        "These give the number of steps to be done in the batch scan. For example, if\n"
        "frequency is scanned from 200 to 300 THz, a value of 11 will perform simulations at:\n"
        "200, 210, 220, ..., 290, 300.");

        QHBoxLayout* loadRow = getRowBoxLayout(entryColumn2Layout);
        buttons["load"] = new QPushButton("Load");
        buttons["load"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        buttons["load"]->setToolTip("Load the results of a previous simulation run.\n"
            "You should select the associated .zip file. The parameters will be\n"
            "loaded into the interface, and the data will be plotted.");
        loadRow->addWidget(buttons["load"]);
        buttons["stop"] = new QPushButton("Stop");
        buttons["stop"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        buttons["stop"]->setToolTip("Tell a currently-running simulation to stop.\n"
            "It might not stop right away; it will only happen once it reaches a break point");
        loadRow->addWidget(buttons["stop"]);
        buttons["run"] = new QPushButton("Run");
        buttons["run"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        buttons["run"]->setToolTip("Run the simulation as currently entered on the\n"
            "interface. If a sequence is entered in the sequence box below,\n"
            "that will execute, otherwise, a simulation on the input parameters\n"
            "above and to the left in a single medium will be performed.");
        loadRow->addWidget(buttons["run"]);
        buttons["save"] = new QPushButton("\xf0\x9f\x92\xbe");
        buttons["save"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        buttons["save"]->setFont(emojiFont);
        buttons["save"]->setToolTip("Save the simulation results.");
        loadRow->addWidget(buttons["save"]);

        QHBoxLayout* fitRow = getRowBoxLayout(entryColumn2Layout);
        addPulldownInContainer(pulldownContainerWidth,fitRow,"fit");
        pulldowns["fit"]->addItem(("Maximize x"));
        pulldowns["fit"]->addItem(("Maximize y"));
        pulldowns["fit"]->addItem(("Maximize Total"));
        pulldowns["fit"]->addItem(("Fit spectrum"));
        pulldowns["fit"]->addItem(("Fit spectrum (log)"));
        pulldowns["fit"]->setToolTip(
            "Select the fitting mode. They will either maximize the power in a certain band\n"
            "and polarization, or match the spectrum to an input file, with weights determined\n"
            "either linearly or logarithmically");
        buttons["loadFitting"] = new QPushButton("\xf0\x9f\x93\x82");
        buttons["loadFitting"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        buttons["loadFitting"]->setFont(emojiFont);
        buttons["loadFitting"]->setToolTip("Load a reference spectrum to fit against");
        fitRow->addWidget(buttons["loadFitting"]);
        buttons["fit"] = new QPushButton("fit");
        buttons["fit"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        buttons["fit"]->setToolTip("Run the fitting routine");
        fitRow->addWidget(buttons["fit"]);

        fitting = new QTextEdit;
        fitting->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        fitting->setFixedSize(rowWidth,2*rowHeight);
        fitting->setTabChangesFocus(true);
        fitting->setToolTip("Here you can enter the commands for the fitting routine.\n"
        "They take the form of a sequence of numbers, in groups of three separated by semicolons\n"
        "The first three numbers are:\n"
        "- start frequency of the range of interest\n"
        "- stop frequency of the ROI\n"
        "- number of iterations to perform\n"
        "The next set gives the parameter to optimize and the range to search:\n"
        "- parameter index (see the batch mode pull downs for the numbers)\n"
        "- minimum value (same units as the GUI uses)\n"
        "- maximum value\n"
        "Example:\n"
        "720e12 780e12 300;\n"
        "29 0 90\n"
        "This will optimize the power in the band from 720-780 THz\n"
        "by changing the phase matching angle theta, in a range from\n"
        "0 to 90 degrees, over 300 iterations.");
        entryColumn2Layout->addWidget(fitting);

        QSpacerItem* consoleSpacer = new QSpacerItem(rowWidth,1);
        entryColumn2Layout->addSpacerItem(consoleSpacer);

        console = new QTextEdit;
        console->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        console->setFixedSize(rowWidth,3*rowHeight-1);
        console->setTabChangesFocus(true);
        console->setFocusPolicy(Qt::NoFocus);
        entryColumn2Layout->addWidget(console);

        //Divide the sequence box into a row of buttons and the text box
        QVBoxLayout* sequenceBoxLayout = new QVBoxLayout(sequenceBox);
        squeezeMargins(sequenceBoxLayout);
        QWidget* sequenceButtonBox = new QWidget;
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
            buttons[entry]->setFixedSize(miniButtonWidth, mainButtonHeight);
            buttons[entry]->setToolTip(tooltip);
            buttons[entry]->setFont(emojiFont);
            sequenceButtonBoxLayout->addWidget(buttons[entry]);
            QObject::connect(buttons[entry], &QPushButton::clicked, action);
        };

        addMiniButton("\xf0\x9f\x93\xb8", "addSameCrystal",
        "Add a fixed crystal which matches the values currently entered on the interface",[&](){
            if(textBoxes["NLAbsorption"]->text().toDouble() != 0.0){
                sPrint("plasma({},{},{},{},{},{},{},{},{})\n",
                pulldowns["material"]->currentIndex(),
                textBoxes["crystalTheta"]->text().toDouble(),
                textBoxes["crystalPhi"]->text().toDouble(),
                textBoxes["NLAbsorption"]->text().toDouble(),
                textBoxes["crystalBandgap"]->text().toDouble(),
                textBoxes["DrudeGamma"]->text().toDouble(),
                textBoxes["effectiveMass"]->text().toDouble(),
                textBoxes["zSize"]->text().toDouble(),
                textBoxes["dz"]->text().toDouble());
            }
            else{
                sPrint("nonlinear({},{},{},{},{})\n",
                pulldowns["material"]->currentIndex(),
                textBoxes["crystalTheta"]->text().toDouble(),
                textBoxes["crystalPhi"]->text().toDouble(),
                textBoxes["zSize"]->text().toDouble(),
                textBoxes["dz"]->text().toDouble());
            }

        });
        addMiniButton("\xe2\x99\x8a", "addDefault", "Insert a crystal that will change with the values set on the "
            "interface, or modified during a batch calculation",[&](){
                sPrint("plasma(d,d,d,d,d,d,d,d,d)\n");
            });
        addMiniButton("\xf0\x9f\x92\xab", "addRotation", "Rotate the polarization by a specified angle in degrees",[&](){
            sPrint("rotate(90)\n");
        });
        addMiniButton("\xf0\x9f\x92\xa1", "addPulse", "Add a new pulse to the grid; values will be set to duplicate "
            "pulse 1 as entered above",[&](){
                sPrint("addPulse({},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})\n",
                textBoxes["energy1"]->text().toDouble(),
                textBoxes["frequency1"]->text().toDouble(),
                textBoxes["bandwidth1"]->text().toDouble(),
                textBoxes["SGOrder1"]->text().toInt(),
                textBoxes["CEP1"]->text().toDouble(),
                textBoxes["delay1"]->text().toDouble(),
                textBoxes["GDD1"]->text().toDouble(),
                textBoxes["TOD1"]->text().toDouble(),
                textBoxes["material1"]->text().toInt(),
                textBoxes["thickness1"]->text().toDouble(),
                textBoxes["beamwaist1"]->text().toDouble(),
                textBoxes["xOffset1"]->text().toDouble(),
                0.0,
                textBoxes["zOffset1"]->text().toDouble(),
                textBoxes["NCAngle1"]->text().toDouble(),
                0.0,
                textBoxes["polarization1"]->text().toDouble(),
                textBoxes["circularity1"]->text().toDouble(),
                pulldowns["material"]->currentIndex(),
                textBoxes["crystalTheta"]->text().toDouble(),
                textBoxes["crystalPhi"]->text().toDouble());
            });
        addMiniButton("\xf0\x9f\x94\x8e", "addMirror", "Add a spherical mirror to the beam path, with radius "
            "of curvature in meters",[&](){
                sPrint("sphericalMirror(-1.0)\n");
            });
        addMiniButton("\xf0\x9f\x98\x8e", "addFilter", "Add a spectral filter to the beam path. "
            "Parameters:\n   central frequency (THz)\n   bandwidth (THz)\n   supergaussian order\n   "
            "in-band amplitude\n   out-of-band amplitude",[&](){
                sPrint("filter(130, 20, 4, 1, 0)\n");
            });
        addMiniButton("\xf0\x9f\x93\x8f", "addLinear", "Add a linear propagation through the crystal entered on the interface",[&](){
            sPrint("linear({},{},{},{},{})\n",
                pulldowns["material"]->currentIndex(),
                textBoxes["crystalTheta"]->text().toDouble(),
                textBoxes["crystalPhi"]->text().toDouble(),
                textBoxes["zSize"]->text().toDouble(),
                textBoxes["dz"]->text().toDouble());
            });
        addMiniButton("\xf0\x9f\x8e\xaf", "addAperture", "Add an aperture to the beam. Parameters:\n   diameter (m)\n   "
            "activation parameter",[&](){
                sPrint("aperture(0.001, 2)\n");
            });
        addMiniButton("\xe2\x9b\xb3", "addFarFieldAperture", "Filter the beam with a far-field aperture. Parameters:\n   "
            "opening angle (deg)\n   activation parameter (k)\n   x-angle (deg)\n   y-angle (deg) ",[&](){
                sPrint("farFieldAperture(2.0,4000,0,0)\n");
            });
        addMiniButton("\xf0\x9f\x94\x81", "addForLoop", "Add an empty for loop. Parameters:\n   "
            "Number of times to execute\n   Variable number in which to put the counter",[&](){
                sPrint("for(10,1){{\n\n}}\n");
            });
        addMiniButton("\xf0\x9f\x97\xba\xef\xb8\x8f", "fdtd", "Run an FDTD simulation, also on a custom grid if you want to.\n"
                "Parameters:\n    "
                "Time divider: factor by which to divide the\n"
                "    time step relative to the time step of the output.\n"
                "    dz (m): custom z-step\n"
                "    front buffer (m): vacuum to insert before the crystal\n"
                "    rear buffer (m): vacuum to insert after the crystal\n"
                "    observation point(m): position at which to\n"
                "        record the time-dependent field.\n"
                "    wait time (s): time to wait before beginning to\n"
                "        record the field in the observation plane\n"
                "    preserve near-field: 0 will remove components which\n"
                "        would not propagate to the far-field. 1 keeps them.\n"
                "    grid index: index in the optics database of a custom grid,\n"
                "        set to d for no custom grid and instead insert a crystal\n"
                "        with the thickness set by the interface.",
                [&](){
                sPrint("fdtd(d,d, 1e-6, 1e-6, 2e-6, 20e-15, 1, d)\n");
            });
        addMiniButton("\xe2\x87\x8c", "fdtdReflection", "Calculate a nonlinear reflection from a surface\n"
            "Parameters:\n"
            "    time step relative to the time step of the output.\n"
            "    dz (m): custom z-step"
            ,[&](){
                    sPrint("fdtdReflection(d,d)\n");
                });
        addMiniButton("\xf0\x9f\xaa\xa9", "applyOptic", "Apply a loaded optic to the field. \nParameters:\n   "
            "Index to the loaded optic",[&](){
                sPrint("applyOptic(0)");
            });
        addMiniButton("\xf0\x9f\x93\x84", "loadOptic", "Load an optic from a file.\n"
            "It should have the format: \nWavelength (nm) | Reflectivity (0.0 to 1.0) | Phase (rad)",[&](){
                std::string path = QFileDialog::getOpenFileName(
                    buttons["loadOptic"],"Load mirror reflectivity and phase","","Plain text file (*.*)").toStdString();
                loadedInputData loadedData(path);
                if(!loadedData.hasData) return;
                theSim.optics.push_back(loadedData);
                cPrint("Loaded {}\ninto index {}",path,theSim.optics.size()-1);
            });
        addMiniButton("\xf0\x9f\x9a\xae", "clearOptics", "Empty out the loaded optics database",[&](){
                theSim.optics.clear();
                cPrint("Optics database now empty.");
            });
        QSpacerItem* miniButtonSpacer = new QSpacerItem(1,1,QSizePolicy::Expanding,QSizePolicy::Fixed);
        sequenceButtonBoxLayout->addSpacerItem(miniButtonSpacer);
        sequence = new QTextEdit;
        sequence->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
        sequence->setFixedWidth(2*rowWidth);
        sequence->setTabChangesFocus(true);
        sequenceBoxLayout->addWidget(sequence);
        sequence->setToolTip("Here you can enter a sequence of events to take place during the simulation\n"
        "The buttons above will enter the commands for a few common things; there are more in the docs.");

        SequenceValidator* contextColors = new SequenceValidator(sequence->document());
    #if QT_VERSION >= QT_VERSION_CHECK(6,5,0)
        QObject::connect(QGuiApplication::styleHints(),
            &QStyleHints::colorSchemeChanged,
            contextColors,
            &QSyntaxHighlighter::rehighlight);
    #endif
        //Put the control strip below the sequence
        QHBoxLayout* simulationControlStripLayout = new QHBoxLayout(simulationControlStrip);
        squeezeMargins(simulationControlStripLayout);
        progress = new QProgressBar;
        simulationControlStripLayout->addSpacerItem(new QSpacerItem(8,1,QSizePolicy::Fixed,QSizePolicy::Fixed));
        simulationControlStripLayout->addWidget(progress);

        checkboxes["FP64"] = new QCheckBox("FP64");
        checkboxes["FP64"]->setToolTip("Determine whether the simulation is run using 32-bit (unchecked)\n"
        "or 64-bit floating point numbers. Checked will be slower but have better precision.");
        simulationControlStripLayout->addSpacerItem(new QSpacerItem(8,1,QSizePolicy::Fixed,QSizePolicy::Fixed));
        simulationControlStripLayout->addWidget(checkboxes["FP64"]);
        addPulldownInContainer(textBoxWidth,simulationControlStripLayout,"primaryHardware");
        addPulldownInContainer(textBoxWidth,simulationControlStripLayout,"secondaryHardware");
        pulldowns["primaryHardware"]->setToolTip("Determine on which hardware/propagation code to run the simulation");
        pulldowns["secondaryHardware"]->setToolTip("Pick a second piece of hardware/propagator to offload a number\n"
        "of simulations in a batch to. For example, if the simulation is running on your GPU, you can give the CPU some work, too");
        textBoxes["offload"] = new QLineEdit;
        textBoxes["offload"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        textBoxes["offload"]->setToolTip("Number of simulations to offload to the second selected hardware.");
        simulationControlStripLayout->addWidget(textBoxes["offload"]);

        QHBoxLayout* plotControlStripLayout = new QHBoxLayout(plotControlStrip);
        squeezeMargins(plotControlStripLayout);
        plotControlStripLayout->addSpacerItem(new QSpacerItem(8,1,QSizePolicy::Fixed,QSizePolicy::Fixed));
        buttons["collapse"] = new QPushButton("\xe2\x86\x94\xef\xb8\x8f");
        buttons["collapse"]->setFixedSize(miniButtonWidth, mainButtonHeight);
        buttons["collapse"]->setToolTip("Collapse the simulation parameters pane and show the plot controls.");
        plotControlStripLayout->addWidget(buttons["collapse"]);

        slider = new QSlider(Qt::Horizontal);
        slider->setFixedHeight(mainButtonHeight);
        plotControlStripLayout->addWidget(slider);
        plotControlStripLayout->addSpacerItem(new QSpacerItem(8,1,QSizePolicy::Fixed,QSizePolicy::Fixed));
        labels["sliderIndex"] = new QLabel("0");
        labels["sliderIndex"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(labels["sliderIndex"]);
        labels["sliderValue"] = new QLabel("0.0");
        labels["sliderValue"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(labels["sliderValue"]);
        buttons["svg"] = new QPushButton("SVG");
        buttons["svg"]->setFixedSize(miniButtonWidth+8, mainButtonHeight);
        buttons["svg"]->setToolTip("Save a .svg file of the plots.");
        plotControlStripLayout->addWidget(buttons["svg"]);
        buttons["xlim"] = new QPushButton("xlim");
        buttons["xlim"]->setFixedSize(miniButtonWidth+8, mainButtonHeight);
        plotControlStripLayout->addWidget(buttons["xlim"]);
        textBoxes["xlimStart"] = new QLineEdit;
        textBoxes["xlimStart"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(textBoxes["xlimStart"]);
        textBoxes["xlimStop"] = new QLineEdit;
        textBoxes["xlimStop"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(textBoxes["xlimStop"]);
        textBoxes["ylimStart"] = new QLineEdit;
        buttons["ylim"] = new QPushButton("ylim");
        buttons["ylim"]->setFixedSize(miniButtonWidth+8, mainButtonHeight);
        plotControlStripLayout->addWidget(buttons["ylim"]);
        textBoxes["ylimStart"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(textBoxes["ylimStart"]);
        textBoxes["ylimStop"] = new QLineEdit;
        textBoxes["ylimStop"]->setFixedSize(textBoxWidth/2,textBoxHeight);
        plotControlStripLayout->addWidget(textBoxes["ylimStop"]);
        checkboxes["Total"] = new QCheckBox("Total");
        checkboxes["Total"]->setToolTip("Overlay the total spectrum.");
        plotControlStripLayout->addSpacerItem(new QSpacerItem(8,1,QSizePolicy::Fixed,QSizePolicy::Fixed));
        checkboxes["Log"] = new QCheckBox("Log");
        checkboxes["Log"]->setToolTip("Plot the spectrum on a log scale. Will look more meaningful if you set y-limits.");
    #ifdef __APPLE__
        plotControlStripLayout->addWidget(checkboxes["Total"], 0, Qt::AlignBottom);
        plotControlStripLayout->addSpacerItem(new QSpacerItem(8,1,QSizePolicy::Fixed,QSizePolicy::Fixed));
        plotControlStripLayout->addWidget(checkboxes["Log"], 0, Qt::AlignBottom);
    #else
        plotControlStripLayout->addWidget(checkboxes["Total"]);
        plotControlStripLayout->addSpacerItem(new QSpacerItem(8,1,QSizePolicy::Fixed,QSizePolicy::Fixed));
        plotControlStripLayout->addWidget(checkboxes["Log"]);
    #endif

    //Layout of plot control area

        checkboxes["plotDark"] = new QCheckBox("Dark mode");
        checkboxes["plotDark"]->setChecked(true);
        checkboxes["plotDark"]->setToolTip("Uncheck to plot with white background");
        plotControlRegionLayout->addWidget(checkboxes["plotDark"]);
        addTextBoxRowForPlots(
            "Font size",
            "fontSize",
            "none",
            "none",
            plotControlRegionLayout,
            "Set font size in the plots\n",
            13, 0, 0);

        checkboxes["showSpaceTimeMomentumFrequency"] = new QCheckBox("Slice images");
        checkboxes["showSpaceTimeMomentumFrequency"]->setChecked(true);
        checkboxes["showSpaceTimeMomentumFrequency"]->setToolTip("Show or hide the space/time and momentum/frequency view.");
        plotControlRegionLayout->addWidget(checkboxes["showSpaceTimeMomentumFrequency"]);
        QObject::connect(checkboxes["showSpaceTimeMomentumFrequency"], &QCheckBox::clicked, [&](){
            if(checkboxes["showSpaceTimeMomentumFrequency"]->isChecked()){
                plots["timeImage1"]->show();
                plots["timeImage2"]->show();
                plots["freqImage1"]->show();
                plots["freqImage2"]->show();
            }
            else{
                plots["timeImage1"]->hide();
                plots["timeImage2"]->hide();
                plots["freqImage1"]->hide();
                plots["freqImage2"]->hide();
            }
        });

        checkboxes["showBeamView"] = new QCheckBox("Beam view");
        checkboxes["showBeamView"]->setChecked(true);
        checkboxes["showBeamView"]->setToolTip("Show or hide the beam view.");
        plotControlRegionLayout->addWidget(checkboxes["showBeamView"]);
        QObject::connect(checkboxes["showBeamView"], &QCheckBox::clicked, [&](){
            if(checkboxes["showBeamView"]->isChecked()){
                plots["beamView"]->show();
            }
            else{
                plots["beamView"]->hide();
            }
        });
        addPulldownInContainer(labelWidth, plotControlRegionLayout, "beamViewMode");
        pulldowns["beamViewMode"]->addItem("Near field");
        pulldowns["beamViewMode"]->addItem("Far field");
        pulldowns["beamViewMode"]->setCurrentIndex(1);
        addTextBoxRowForPlots(
            "Red (f, \xcf\x83, I)",
            "Red_frequency",
            "Red_sigma",
            "Red_strength",
            plotControlRegionLayout,
            "Control the red color channel in the beam view.\n"
            "Central frequency (THz), bandwidth (THz), gain (multiplier)",
            400, 60, 1);
        addTextBoxRowForPlots(
            "Green (f, \xcf\x83, I)",
            "Green_frequency",
            "Green_sigma",
            "Green_strength",
            plotControlRegionLayout,
            "Control the green color channel in the beam view.\n"
            "Central frequency (THz), bandwidth (THz), gain (multiplier)",
            550, 60, 1);
        addTextBoxRowForPlots(
            "Blue (f, \xcf\x83, I)",
            "Blue_frequency",
            "Blue_sigma",
            "Blue_strength",
            plotControlRegionLayout,
            "Control the blue color channel in the beam view.\n"
            "Central frequency (THz), bandwidth (THz), gain (multiplier)",
            680, 90, 1);
        addTextBoxRowForPlots(
                "Farfield \xe2\x88\xa0",
                "farfieldAngle",
                "none",
                "none",
                plotControlRegionLayout,
                "Set the solid angle shown in the farfield view (in degrees)\n",
                5, 0, 0);
        checkboxes["rotate2D"] = new QCheckBox("Rotate 2D around center");
        checkboxes["rotate2D"]->setChecked(true);
        checkboxes["rotate2D"]->setToolTip(
            "Produce images in 2D mode by rotating around the center-of-mass.\n"
            "This will produce ring artifacts for frequencies that don't correspond\n"
            "to a single beam.");
        plotControlRegionLayout->addWidget(checkboxes["rotate2D"]);
        buttons["render"] = new QPushButton("Render");
        buttons["render"]->setFixedSize(mainButtonWidth,mainButtonHeight);
        buttons["render"]->setToolTip("Render the beam view with these parameters");
        plotControlRegionLayout->addWidget(buttons["render"]);
        checkboxes["showWavesAndSpectra"] = new QCheckBox("Show waves and spectra");
        checkboxes["showWavesAndSpectra"]->setChecked(true);
        checkboxes["showWavesAndSpectra"]->setToolTip(
            "show the waveforms and spectra");
        plotControlRegionLayout->addWidget(checkboxes["showWavesAndSpectra"]);
        QObject::connect(checkboxes["showWavesAndSpectra"], &QCheckBox::clicked, [&](){
            if(checkboxes["showWavesAndSpectra"]->isChecked()){
                if(!checkboxes["combinePolarizations"]->isChecked()){
                    plots["timePlot2"]->show();
                    plots["freqPlot2"]->show();
                }
                plots["timePlot1"]->show();
                plots["freqPlot1"]->show();
            }
            else{
                plots["timePlot1"]->hide();
                plots["freqPlot1"]->hide();
                plots["timePlot2"]->hide();
                plots["freqPlot2"]->hide();
            }
        });

        checkboxes["combinePolarizations"] = new QCheckBox("Combine polarizations");
        checkboxes["combinePolarizations"]->setChecked(false);
        checkboxes["combinePolarizations"]->setToolTip(
            "Combine the waveform and spectra plots, showing the two polarization components together");
        plotControlRegionLayout->addWidget(checkboxes["combinePolarizations"]);
        QObject::connect(checkboxes["combinePolarizations"], &QCheckBox::clicked, [&](){
            if(checkboxes["combinePolarizations"]->isChecked()){
                plots["timePlot2"]->hide();
                plots["freqPlot2"]->hide();
            }
            else if(checkboxes["showWavesAndSpectra"]->isChecked()){
                plots["timePlot2"]->show();
                plots["freqPlot2"]->show();
            }
        });


        //squeezeMargins(plotControlRegionLayout);

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
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["freqPlot1"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["beamView"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["freqPlot2"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["timePlot1"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["timePlot2"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["timeImage1"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["timeImage2"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["freqImage1"], &CairoWidget::queueUpdate);
        QObject::connect(messenger, &GuiMessenger::requestUpdate, plots["freqImage2"], &CairoWidget::queueUpdate);

        QObject::connect(slider, &QSlider::valueChanged, plots["freqPlot1"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, this, &LWEGui::renderBeamView);
        QObject::connect(slider, &QSlider::valueChanged, plots["freqPlot2"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["timePlot1"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["timePlot2"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["timeImage1"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["timeImage2"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["freqImage1"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, plots["freqImage2"], &CairoWidget::queueUpdate);
        QObject::connect(slider, &QSlider::valueChanged, this, [&](){
            int i = slider->value();
            labels["sliderIndex"]->setText(QString::number(i));

            if(theSim.base().batchIndex2 == 0){
                if(theSim.base().Nsims < 2) return;
                double batchStart = theSim.base().getByNumberWithMultiplier(theSim.base().batchIndex);
                double batchValue =
                batchStart +
                i * (theSim.base().batchDestination - batchStart)/(theSim.base().Nsims - 1);
                labels["sliderValue"]->setFixedWidth(60);
                labels["sliderValue"]->setText(QString::fromStdString(Sformat("{:.3g}",batchValue)));
            }
            else{
                if(theSim.base().Nsims < 2 && theSim.base().Nsims2 < 2) return;
                double batchStart = theSim.base().getByNumberWithMultiplier(theSim.base().batchIndex);
                double batchValue =
                batchStart +
                (i % theSim.base().Nsims) * (theSim.base().batchDestination - batchStart)/(theSim.base().Nsims - 1);

                double batchStart2 = theSim.base().getByNumberWithMultiplier(theSim.base().batchIndex2);
                double batchValue2 =
                batchStart2 +
                (i / theSim.base().Nsims) * (theSim.base().batchDestination2 - batchStart2)/(theSim.base().Nsims2 - 1);
                labels["sliderValue"]->setFixedWidth(120);
                labels["sliderValue"]->setText(QString::fromStdString(Sformat("{:.3g}, {:.3g}",batchValue,batchValue2)));
            }

        });

        QObject::connect(messenger, &GuiMessenger::moveSlider, slider, &QSlider::setValue);
        QObject::connect(messenger, &GuiMessenger::requestSyncValues, this, &LWEGui::queueSyncValues);
        QObject::connect(messenger, &GuiMessenger::requestSliderUpdate, this, &LWEGui::updateSlider);
        QObject::connect(messenger, &GuiMessenger::requestBeamViewUpdate, this, &LWEGui::updateBeamView);
        QObject::connect(messenger, &GuiMessenger::requestBeamViewRender, this, &LWEGui::renderBeamView);
        QObject::connect(messenger, &GuiMessenger::sendProgressRange, progress, &QProgressBar::setRange);
        QObject::connect(messenger, &GuiMessenger::sendProgressValue, progress, &QProgressBar::setValue);

        windowBody->show();
        connectButtons();
    }
    void connectButtons(){
        QObject::connect(buttons["run"], &QPushButton::clicked, [&](){
            if(theSim.base().isRunning) return;
            std::unique_lock guiLock(m);
            readParametersFromInterface(theSim);
            theSim.configure();
            updateSlider();
            simulationRun theRun(pulldowns["primaryHardware"]->currentIndex(),checkboxes["FP64"]->isChecked(),theSim);
            simulationRun theOffloadRun(pulldowns["secondaryHardware"]->currentIndex(),checkboxes["FP64"]->isChecked(),theSim);
            std::thread(mainSimThread, std::ref(*this), theRun, theOffloadRun).detach();
        });

        QObject::connect(buttons["fit"], &QPushButton::clicked, [&](){
            if(theSim.base().isRunning) return;
            std::unique_lock guiLock(m);
            readParametersFromInterface(theSim);
            theSim.configure();
            simulationRun theRun(pulldowns["primaryHardware"]->currentIndex(),checkboxes["FP64"]->isChecked(),theSim);
            std::thread(fittingThread, std::ref(*this), theRun).detach();
        });

        QObject::connect(buttons["save"], &QPushButton::clicked, [&](){
            if(!theSim.base().isGridAllocated) return;
            std::string path = QFileDialog::getSaveFileName(buttons["save"],"Save LWE result","","LWE Results (*.zip)").toStdString();
            if(path.empty()) return;
            theSim.base().outputBasePath = path;
            stripLineBreaks(theSim.base().outputBasePath);
            if ((theSim.base().outputBasePath.length() > 4 && theSim.base().outputBasePath.substr(theSim.base().outputBasePath.length() - 4) == ".txt")
                || (theSim.base().outputBasePath.length() > 4 && theSim.base().outputBasePath.substr(theSim.base().outputBasePath.length() - 4) == ".zip")) {
                theSim.base().outputBasePath = theSim.base().outputBasePath.substr(0, theSim.base().outputBasePath.length() - 4);
            }

            messenger->passString("Saving...");
            auto saveLambda = [&](){
                std::unique_lock lock(m);
                theSim.saveDataSet();
                messenger->passString("done.\n");
                messenger->passDrawRequest();
            };
            std::thread(saveLambda).detach();
        });
        QObject::connect(buttons["saveSLURM"], &QPushButton::clicked, [&](){
            if(theSim.base().isRunning) return;
            std::string path = QFileDialog::getSaveFileName(buttons["saveSLURM"],"Save cluster script","","LWE Results (*.zip)").toStdString();
            if(path.empty()) return;
            theSim.base().outputBasePath = path;
            stripLineBreaks(theSim.base().outputBasePath);
            if ((theSim.base().outputBasePath.length() > 4 && theSim.base().outputBasePath.substr(theSim.base().outputBasePath.length() - 4) == ".txt")
                || (theSim.base().outputBasePath.length() > 4 && theSim.base().outputBasePath.substr(theSim.base().outputBasePath.length() - 4) == ".zip")) {
                theSim.base().outputBasePath = theSim.base().outputBasePath.substr(0, theSim.base().outputBasePath.length() - 4);
            }
            theSim.base().optics = theSim.optics;
            createRunFile(*this);
        });

        QObject::connect(buttons["load"], &QPushButton::clicked, [&](){
            if(theSim.base().isRunning) return;
            std::string path = QFileDialog::getOpenFileName(buttons["load"],"Load LWE result","","LWE Results (*.zip);;LWE Inputs (*.txt)").toStdString();
            if(path.empty()) return;
            std::thread([&](std::string path){
                std::unique_lock lock(m);
                bool isZipFile = (path.length() >= 4
                && path.substr(path.length()-4)==".zip");
                messenger->passString("Loading...");
                int readParameters =
                theSim.base().readInputParametersFile(theDatabase.db.data(), path);
                if(isZipFile) theSim.loadOptics(path);
                theSim.configure(true);
                std::for_each(theSim.mutexes.begin(), theSim.mutexes.end(),
                        [](std::shared_mutex& m) {std::lock_guard<std::shared_mutex> lock(m); });
                if (readParameters == 61 && theSim.base().isGridAllocated) {
                    int64_t extensionLoc = path.find_last_of(".");
                    const std::string basePath = path.substr(0, extensionLoc);
                    theSim.base().loadSavedFields(basePath, isZipFile);
                    if(theSim.base().pulse1LoadedData.hasData){
                        pulse1LoadedData = theSim.base().pulse1LoadedData;
                        messenger->passString("Loaded stored pulse 1 file into buffer.");
                    }
                    if(theSim.base().pulse2LoadedData.hasData){
                        pulse2LoadedData = theSim.base().pulse2LoadedData;
                        messenger->passString("Loaded stored pulse 2 file into buffer.");
                    }
                }
                messenger->passDrawRequest();
                messenger->passSyncValues();
                messenger->passSliderUpdate();
                messenger->passBeamViewRender();
                messenger->passString("done.");
            },path).detach();
        });

        QObject::connect(buttons["loadPulse1"], &QPushButton::clicked, [&](){
            std::string path = QFileDialog::getOpenFileName(buttons["loadPulse1"],"Load field data","","ASCII data (*.*)").toStdString();
            if(path.empty()) return;
            pulse1LoadedData = loadedInputData(path);
            messenger->passString(Sformat("Loaded new file into pulse 1 buffer:\n{}\n", pulse1LoadedData.filePath));
        });

        QObject::connect(buttons["loadPulse2"], &QPushButton::clicked, [&](){
            std::string path = QFileDialog::getOpenFileName(buttons["loadPulse2"],"Load field data","","ASCII data (*.*)").toStdString();
            if(path.empty()) return;
            pulse2LoadedData = loadedInputData(path);
            messenger->passString(Sformat("Loaded new file into pulse 2 buffer:\n{}\n", pulse2LoadedData.filePath));
        });

        QObject::connect(buttons["loadFitting"], &QPushButton::clicked, [&](){
            std::string path = QFileDialog::getOpenFileName(buttons["loadFitting"],"Load spectral target data","","ASCII data (*.*)").toStdString();
            if(path.empty()) return;
            fittingLoadedData = loadedInputData(path);
            messenger->passString(Sformat("Loaded new fitting spectral target:\n{}\n", fittingLoadedData.filePath));
        });

        QObject::connect(buttons["loadMaterial"], &QPushButton::clicked, [&](){
            std::string path = QFileDialog::getOpenFileName(buttons["LoadMaterial"],"Load crystal database","","ASCII data (*.*)").toStdString();
            if(path.empty()) return;
            theDatabase = crystalDatabase(path);
            populateDatabasePulldown();
            messenger->passString(Sformat("Loaded new crystal database:\n{}\n", path));
        });

        QObject::connect(buttons["stop"], &QPushButton::clicked, [&](){
            if (theSim.base().isRunning) {
                for(simulationParameterSet& sim : theSim.parameters){
                    sim.cancellationCalled = true;
                }
            }
        });

        QObject::connect(buttons["collapse"], &QPushButton::clicked, [&](){
            if(isInputRegionHidden){
                inputRegion->show();
                isInputRegionHidden = false;

                plotControlRegion->hide();
            }
            else{
                inputRegion->hide();
                isInputRegionHidden = true;
                plotControlRegion->show();
            }
        });

        QObject::connect(buttons["render"], &QPushButton::clicked, [&](){
            messenger->passBeamViewRender();
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

        QObject::connect(buttons["svg"], &QPushButton::clicked, [&](){
            std::string SVGpath = QFileDialog::getSaveFileName(buttons["svg"],"Save SVG file of plots","","Scalable Vector Graphics (*.svg)").toStdString();
            if(SVGpath.empty()) return;
            isMakingSVG = true;


            //get dimensions lambda
            auto getDimensions = [](const std::string& SVGstr) {
                std::regex width_regex("width=\"(\\d+)\"");
                std::regex height_regex("height=\"(\\d+)\"");
                std::smatch width_match, height_match;

                std::regex_search(SVGstr, width_match, width_regex);
                std::regex_search(SVGstr, height_match, height_regex);
                int width = std::stoi(width_match[1].str());
                int height = std::stoi(height_match[1].str());
                return std::pair<int,int>(width, height);
            };
            for(auto& str : SVGStrings){
                str.clear();
            }
            for(auto& pair : plots){
                if(pair.second->isVisible()){
                    pair.second->repaint();
                }
            }

            const bool showSlices = checkboxes["showSpaceTimeMomentumFrequency"]->isChecked();
            const bool showBeamView = checkboxes["showBeamView"]->isChecked();
            const bool showWaves = checkboxes["showWavesAndSpectra"] ->isChecked();
            const bool combineWaves = checkboxes["combinePolarizations"]->isChecked();
            bool startedSVG = false;
            int totalHeight = 0;
            int totalWidth = 0;
            int height_slice = 0;
            int width_slice = 0;
            int height_wave = 0;
            int width_wave = 0;
            int firstPlot = -1;
            int verticalElements = 0;
            int horizontalElements = 2;

            if(showSlices){
                startedSVG = true;
                firstPlot = 0;
                auto width_height = getDimensions(SVGStrings[0]);
                width_slice = width_height.first;
                height_slice = width_height.second;
                totalWidth = 2 * width_slice;
                totalHeight = 2 * height_slice;
                verticalElements += 2;
            }
            if(showBeamView){
                if(!startedSVG){
                    firstPlot = 8;
                    horizontalElements = 1;
                }
                startedSVG = true;
                auto width_height = getDimensions(SVGStrings[8]);
                totalWidth = width_height.first;
                totalHeight += width_height.second;
                verticalElements++;
            }
            if(showWaves){
                if(!startedSVG) firstPlot = 4;
                startedSVG = true;
                auto width_height = getDimensions(SVGStrings[4]);
                horizontalElements = 2;
                width_wave = width_height.first;
                height_wave = width_height.second;
                totalWidth = 2*width_wave;
                totalHeight += height_wave;
                if(!combineWaves){
                    totalHeight += height_wave;
                    verticalElements++;
                }
                verticalElements++;
            }

            std::size_t SVGbegin = SVGStrings[firstPlot].find("<svg");
            SVGStrings[firstPlot].insert(SVGbegin,Sformat("<svg cool=\"1\" width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}"
                "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink="
                "\"http://www.w3.org/1999/xlink\">\n",
                totalWidth, totalHeight, totalWidth, totalHeight));


            auto appendPlot = [&](int index, int x, int y){
                SVGbegin = SVGStrings[index].find("width=");
                SVGStrings[index].insert(SVGbegin,Sformat("x=\"{}\" y=\"{}\" ",
                x * (totalWidth/horizontalElements), y * (totalHeight/verticalElements)));
                SVGbegin = SVGStrings[index].find("<svg");
                SVGStrings[index] = SVGStrings[index].substr(SVGbegin);

            };

            int y_ind = 0;
            if(showSlices){
                appendPlot(1,0,1);
                appendPlot(2,1,0);
                appendPlot(3,1,1);

                y_ind = 2;
            }
            if(showSlices && showBeamView){
                appendPlot(8,0,y_ind);
            }
            if(showBeamView) y_ind++;
            if(showWaves && (showSlices || showBeamView)){
                appendPlot(4,0,y_ind);
            }
            if(showWaves){
                appendPlot(6,1,y_ind);
                y_ind++;
                if(!combineWaves){
                    appendPlot(5,0,y_ind);
                    appendPlot(7,1,y_ind);
                }
            }


            isMakingSVG = false;
            std::ofstream fs(SVGpath);

            fs << SVGStrings[0];
            fs << SVGStrings[1];
            fs << SVGStrings[2];
            fs << SVGStrings[3];
            fs << SVGStrings[8];
            fs << SVGStrings[4];
            fs << SVGStrings[5];
            fs << SVGStrings[6];
            fs << SVGStrings[7];


            fs << "</svg>";
        });

    }
    public slots:
    void queueSyncValues(){
        setInterfaceValuesToActiveValues(theSim);
    }
    void updateSlider(){
        slider->setMinimum(0);
        slider->setMaximum(theSim.base().Nsims * theSim.base().Nsims2 - 1);
        labels["sliderIndex"]->setText(QString::number(slider->value()));
    }
    void renderBeamView(){
        if(plots["beamView"]->isVisible()){
            int64_t simIndex = maxN(0,slider->value());
            VisualizationType type = static_cast<VisualizationType>(pulldowns["beamViewMode"]->currentIndex());
            if (simIndex > theSim.base().Nsims * theSim.base().Nsims2) {
                simIndex = 0;
            }
            beamView.height = theSim.base().Nspace;
            beamView.width = theSim.base().Nspace;
            beamView.hasFullSizeRenderedImage = true;
            VisualizationConfig config(&(theSim.base()), type);
            config.red_f0 = 1e12 * textBoxes["Red_frequency"]->text().toDouble();
            config.green_f0 = 1e12 * textBoxes["Green_frequency"]->text().toDouble();
            config.blue_f0 = 1e12 * textBoxes["Blue_frequency"]->text().toDouble();
            config.red_sigma = 1e12 * textBoxes["Red_sigma"]->text().toDouble();
            config.green_sigma = 1e12 * textBoxes["Green_sigma"]->text().toDouble();
            config.blue_sigma = 1e12 * textBoxes["Blue_sigma"]->text().toDouble();
            config.red_amplitude = 1e-11 * textBoxes["Red_strength"]->text().toDouble();
            config.green_amplitude = 1e-11 * textBoxes["Green_strength"]->text().toDouble();
            config.blue_amplitude = 1e-11 * textBoxes["Blue_strength"]->text().toDouble();
            config.maxAngle = deg2Rad<float>() * textBoxes["farfieldAngle"]->text().toDouble();
            config.dTheta = 2*config.maxAngle/config.Nx;
            config.simIndex = simIndex;
            config.result_pixels = &beamView.pixels;
            config.rotate2D = checkboxes["rotate2D"]->isChecked();
            std::thread(drawBeamThread,std::ref(*this),config).detach();
        }
    }
    void updateBeamView(){
        plots["beamView"]->repaint();
    }

};

int main(int argc, char **argv){
    QApplication app(argc, argv);
    app.setApplicationName("Lightwave Explorer");
    #ifdef __linux__
    app.setDesktopFileName("io.github.NickKarpowicz.LightwaveExplorer");
    #endif
    //Win10 dark theme via https://forum.qt.io/topic/101391/windows-10-dark-theme/4
#ifdef Q_OS_WIN
    QSettings settings("HKEY_CURRENT_USER\\Software\\Microsoft\\Windows\\CurrentVersion\\Themes\\Personalize",QSettings::NativeFormat);
    if(QSysInfo::productVersion().toInt() < 11){
        app.setStyle(QStyleFactory::create("Fusion"));
        QPalette darkPalette;
        QColor darkColor = QColor(45,45,45);
        QColor disabledColor = QColor(127,127,127);
        darkPalette.setColor(QPalette::Window, darkColor);
        darkPalette.setColor(QPalette::WindowText, Qt::white);
        darkPalette.setColor(QPalette::Base, QColor(18,18,18));
        darkPalette.setColor(QPalette::AlternateBase, darkColor);
        darkPalette.setColor(QPalette::ToolTipBase, Qt::white);
        darkPalette.setColor(QPalette::ToolTipText, Qt::white);
        darkPalette.setColor(QPalette::Text, Qt::white);
        darkPalette.setColor(QPalette::Disabled, QPalette::Text, disabledColor);
        darkPalette.setColor(QPalette::Button, darkColor);
        darkPalette.setColor(QPalette::ButtonText, Qt::white);
        darkPalette.setColor(QPalette::Disabled, QPalette::ButtonText, disabledColor);
        darkPalette.setColor(QPalette::BrightText, Qt::red);
        darkPalette.setColor(QPalette::Link, QColor(42, 130, 218));

        darkPalette.setColor(QPalette::Highlight, QColor(42, 130, 218));
        darkPalette.setColor(QPalette::HighlightedText, Qt::black);
        darkPalette.setColor(QPalette::Disabled, QPalette::HighlightedText, disabledColor);

        app.setPalette(darkPalette);

        app.setStyleSheet("QToolTip { color: #ffffff; background-color: #2a82da; border: 1px solid white; }");
    }
#endif

    LWEGui Gui;
    return app.exec();
}


void setSYCLvars(std::string& s) {
#ifdef _WIN32
    wchar_t loadBuffer[1024];
    DWORD envcount = GetEnvironmentVariableW(L"SYCL_CACHE_PERSISTENT", loadBuffer, 16);
    if (envcount == 0) {
        STARTUPINFO si;
        PROCESS_INFORMATION pi;

        ZeroMemory(&si, sizeof(si));
        si.cb = sizeof(si);
        ZeroMemory(&pi, sizeof(pi));
        wchar_t setSYCLpersistent[] = L"setx SYCL_CACHE_PERSISTENT 1";
        // Start the child process.
        CreateProcess(NULL,   // No module name (use command line)
            setSYCLpersistent,        // Command line
            NULL,           // Process handle not inheritable
            NULL,           // Thread handle not inheritable
            false,          // Set handle inheritance to false
            CREATE_NO_WINDOW,              // No creation flags
            NULL,           // Use parent's environment block
            NULL,           // Use parent's starting directory
            &si,            // Pointer to STARTUPINFO structure
            &pi);           // Pointer to PROCESS_INFORMATION structure

        // Wait until child process exits.
        WaitForSingleObject(pi.hProcess, INFINITE);

        // Close process and thread handles.
        CloseHandle(pi.hProcess);
        CloseHandle(pi.hThread);
        s.append("Note: SYCL cache was disabled, now it's enabled.\n"
        "I recommend you restart LWE to avoid long delays when you start a simulation.\n"
        "Now only the first run will be slow.\n");
    }
#endif
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

#ifdef _WIN32
    auto checkDLL = [](const char* name){
        HMODULE hModule = LoadLibraryA(name);
        if(hModule){
            FreeLibrary(hModule);
            return true;
        }
        return false;
    };
#else
    [[maybe_unused]] auto checkDLL = []([[maybe_unused]] const char* name){
        return true;
    };
#endif

#ifndef NOCUDA
    if(checkDLL("nvml.dll")){
        int CUDAdevice;
        cudaGetDeviceCount(&theSim.base().cudaGPUCount);
        cudaGetDevice(&CUDAdevice);
        struct cudaDeviceProp activeCUDADeviceProp;

        if (theSim.base().cudaGPUCount > 0) {
            theSim.base().CUDAavailable = true;
            if (theSim.base().cudaGPUCount == 1) {
                s.append(Sformat("CUDA found a GPU:\n", theSim.base().cudaGPUCount));
            }
            else {
                s.append(Sformat("CUDA found {} GPU(s):\n", theSim.base().cudaGPUCount));
            }
            for (int i = 0; i < theSim.base().cudaGPUCount; ++i) {
                cudaGetDeviceProperties(&activeCUDADeviceProp, CUDAdevice);
                s.append(Sformat("   {}\n",
                    activeCUDADeviceProp.name));
            }
        }
    }


#else
#define solveNonlinearWaveEquationSequence solveNonlinearWaveEquationSequenceCPU
#define solveNonlinearWaveEquation solveNonlinearWaveEquationCPU
#define runDlibFitting runDlibFittingCPU
#endif

#ifndef NOSYCL
    bool isIntelRuntimeInstalled = checkDLL("sycl8.dll");
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
        setSYCLvars(s);
    }
    else {
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

    std::vector<simulationParameterSet> counterVector = theSim.getParameterVector();
    std::atomic_uint32_t& totalSteps = theGui.totalSteps;
    std::atomic_uint32_t& progressCounter = theGui.progressCounter;
    totalSteps = 0;
    progressCounter = 0;

    try {
        for (int64_t j = 0; j < theSim.base().Nsims * theSim.base().Nsims2; j++) {
            totalSteps += 1;
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
    theGui.messenger->passProgressRange(0,static_cast<int>(totalSteps));
    theSim.base().isRunning = true;
    auto progressThread = [&](){
        while(theSim.base().isRunning){
            theGui.messenger->passProgressValue(static_cast<int>(theGui.progressCounter));
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
            std::unique_lock lock(theSim.mutexes.at(j));
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
            theGui.progressCounter += 1;
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
    theSim.base().isRunning = false;
    theGui.messenger->passProgressValue(static_cast<int>(theGui.progressCounter));
    theGui.messenger->passBeamViewRender();
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
		if(1 == sim.sCPU()->readInputParametersFile(db.db.data(), sysPathIni.c_str())){
            sim.sCPU()->readInputParametersFile(db.db.data(), "DefaultValues.ini");
        }
#else
		sim.sCPU()->readInputParametersFile(db.db.data(), "DefaultValues.ini");
#endif
}

void drawBeamThread(LWEGui& theGui, VisualizationConfig config){
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    std::shared_lock dataLock(theGui.theSim.mutexes.at(config.simIndex), std::try_to_lock);
    if(!(dataLock.owns_lock()) || !(guiLock.owns_lock())) return;
    std::unique_lock imageLock(theGui.beamViewMutex);
    if(&theGui.theSim.base() != config.sCPU) return;
    renderVisualizationCPU(config);
    theGui.messenger->requestUpdate();
}
void drawBeamImage(cairo_t* cr, int width, int height, LWEGui& theGui) {
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    if (!theGui.theSim.base().isGridAllocated) {
        blackoutCairoPlot(cr,width,height);
        return;
    }
    if (!(guiLock.owns_lock())) {
        blackoutCairoPlot(cr,width,height,true);
        return;
    }
    int64_t simIndex = maxN(0,theGui.slider->value());
    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    LwePlot<double> sPlot;
    if(!(theGui.checkboxes["plotDark"]->isChecked())){
        sPlot.backgroundColor =  LweColor(1,1,1,0);
        sPlot.axisColor = LweColor(0.0, 0.0, 0.0, 0);
        sPlot.textColor = sPlot.axisColor;
    }
    sPlot.fontSize = theGui.textBoxes["fontSize"]->text().toDouble();
    sPlot.height = height;
    sPlot.width = width;
    if(theGui.pulldowns["beamViewMode"]->currentIndex() == 1){
        sPlot.xLabel = "\xce\xb8 (deg)";
        sPlot.yLabel = "\xcf\x86 (deg)";
        sPlot.forcedXmax = theGui.textBoxes["farfieldAngle"]->text().toDouble();
        sPlot.forcedXmin = -sPlot.forcedXmax.value();
        sPlot.forcedYmin = theGui.textBoxes["farfieldAngle"]->text().toDouble();
        sPlot.forcedYmax = -sPlot.forcedYmin.value();
    }
    else{
        sPlot.xLabel = "x (\xce\xbcm)";
        sPlot.yLabel = "y (\xce\xbcm)";
        sPlot.forcedXmax = 0.5e6 * theGui.theSim.base().Nspace * theGui.theSim.base().rStep;
        sPlot.forcedXmin = -sPlot.forcedXmax.value();
        sPlot.forcedYmin = 0.5e6 * theGui.theSim.base().Nspace * theGui.theSim.base().rStep;
        sPlot.forcedYmax = -sPlot.forcedYmin.value();
    }

    sPlot.makeSVG = theGui.isMakingSVG;
    std::shared_lock imageLock(theGui.beamViewMutex, std::try_to_lock);
    if (imageLock.owns_lock()) {
        sPlot.image = &theGui.beamView;
        sPlot.plot(cr);
        if(theGui.isMakingSVG){
            theGui.SVGStrings[8] = sPlot.SVGString;
        }
    }
    else blackoutCairoPlot(cr, width, height, true);
}

void drawTimeImage1(cairo_t* cr, int width, int height, LWEGui& theGui) {
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    if (!theGui.theSim.base().isGridAllocated || !(guiLock.owns_lock())) {
        blackoutCairoPlot(cr,width,height);
        return;
    }

    int64_t simIndex = maxN(0,theGui.slider->value());
    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }
    int64_t cubeMiddle = theGui.theSim.base().Ntime * theGui.theSim.base().Nspace * (theGui.theSim.base().Nspace2 / 2);

    LweImage<double> image;
    image.data = DataSlice<double,2>{&theGui.theSim.base().ExtOut[simIndex * theGui.theSim.base().Ngrid * 2 + cubeMiddle],
        {static_cast<size_t>(theGui.theSim.base().Ntime), static_cast<size_t>(theGui.theSim.base().Nspace)}};

    image.height = height;
    image.width = width;
    image.colorMap = ColorMap::cyan_magenta;


    LwePlot<double> sPlot;
    if(!(theGui.checkboxes["plotDark"]->isChecked())){
        sPlot.backgroundColor =  LweColor(1,1,1,0);
        sPlot.axisColor = LweColor(0.0, 0.0, 0.0, 0);
        sPlot.textColor = sPlot.axisColor;
    }
    sPlot.fontSize = theGui.textBoxes["fontSize"]->text().toDouble();
    sPlot.height = height;
    sPlot.width = width;
    sPlot.image = &image;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "x (mm)";
    sPlot.unitY = 1e3;
    sPlot.forcedXmax = 0.5e15 * theGui.theSim.base().Ntime * theGui.theSim.base().tStep;
    sPlot.forcedXmin = -sPlot.forcedXmax.value();
    sPlot.forcedYmin = 0.5e3 * theGui.theSim.base().spatialWidth;
    sPlot.forcedYmax = -sPlot.forcedYmin.value();
    sPlot.makeSVG = theGui.isMakingSVG;
    std::shared_lock dataLock(theGui.theSim.mutexes.at(simIndex),std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG){
            theGui.SVGStrings[0] = sPlot.SVGString;
        }
    }
    else blackoutCairoPlot(cr, width, height);
}

void drawField1Plot(cairo_t* cr, int width, int height, LWEGui& theGui) {
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    if (!theGui.theSim.base().isGridAllocated || !(guiLock.owns_lock())) {
        blackoutCairoPlot(cr,width,height);
        return;
    }
    LwePlot<double> sPlot;
    int64_t simIndex = maxN(0,theGui.slider->value());

    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = theGui.theSim.base().Ntime * theGui.theSim.base().Nspace * (theGui.theSim.base().Nspace2 / 2);
    if(!(theGui.checkboxes["plotDark"]->isChecked())){
        sPlot.backgroundColor =  LweColor(1,1,1,0);
        sPlot.axisColor = LweColor(0.0, 0.0, 0.0, 0);
        sPlot.textColor = sPlot.axisColor;
    }
    sPlot.makeSVG = theGui.isMakingSVG;
    sPlot.fontSize = theGui.textBoxes["fontSize"]->text().toDouble();
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theGui.theSim.base().tStep * 1e15;
    sPlot.x0 = -0.5e15 * (theGui.theSim.base().Ntime) * theGui.theSim.base().tStep;
    sPlot.forcedXmin = sPlot.x0;
    sPlot.forcedXmax = -sPlot.x0;
    sPlot.data.push_back(&theGui.theSim.base().ExtOut[
        simIndex * theGui.theSim.base().Ngrid * 2 + cubeMiddle + theGui.theSim.base().Ntime * theGui.theSim.base().Nspace / 2]);

    sPlot.Npts = theGui.theSim.base().Ntime;
    sPlot.lineColors = {LweColor(0, 1, 1, 1)};
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ex (GV/m)";
    if(theGui.checkboxes["combinePolarizations"]->isChecked()){
        sPlot.lineColors.push_back(LweColor(1, 0, 1, 1));
        sPlot.data.push_back(&theGui.theSim.base().ExtOut[
        theGui.theSim.base().Ngrid + simIndex * theGui.theSim.base().Ngrid * 2
            + cubeMiddle + theGui.theSim.base().Ntime * theGui.theSim.base().Nspace / 2]);
        sPlot.yLabel = "E (GV/m)";
    }
    sPlot.unitY = 1e9;
    std::shared_lock dataLock(theGui.theSim.mutexes.at(simIndex),std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG) {
            theGui.SVGStrings[4] = sPlot.SVGString;
        }
    }
    else blackoutCairoPlot(cr, width, height);
}

void drawField2Plot(cairo_t* cr, int width, int height, LWEGui& theGui) {
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    if (!theGui.theSim.base().isGridAllocated || !(guiLock.owns_lock())) {
        blackoutCairoPlot(cr,width,height);
        return;
    }
    LwePlot<double> sPlot;

    int64_t simIndex = maxN(0,theGui.slider->value());

    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle =
        theGui.theSim.base().Ntime * theGui.theSim.base().Nspace * (theGui.theSim.base().Nspace2 / 2);
    if(!(theGui.checkboxes["plotDark"]->isChecked())){
        sPlot.backgroundColor =  LweColor(1,1,1,0);
        sPlot.axisColor = LweColor(0.0, 0.0, 0.0, 0);
        sPlot.textColor = sPlot.axisColor;
    }
    sPlot.makeSVG = theGui.isMakingSVG;
    sPlot.fontSize = theGui.textBoxes["fontSize"]->text().toDouble();
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theGui.theSim.base().tStep / 1e-15;
    sPlot.x0 = -0.5e15 * (theGui.theSim.base().Ntime) * theGui.theSim.base().tStep;
    sPlot.forcedXmin = sPlot.x0;
    sPlot.forcedXmax = -sPlot.x0;
    sPlot.data.push_back(
        &theGui.theSim.base().ExtOut[
        theGui.theSim.base().Ngrid + simIndex * theGui.theSim.base().Ngrid * 2
            + cubeMiddle + theGui.theSim.base().Ntime * theGui.theSim.base().Nspace / 2]);
    sPlot.Npts = theGui.theSim.base().Ntime;
    sPlot.lineColors = {LweColor(1, 0, 1, 1)};
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "Ey (GV/m)";
    sPlot.unitY = 1e9;
    std::shared_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG) {
            theGui.SVGStrings[5] = sPlot.SVGString;
        }
    }
    else blackoutCairoPlot(cr, width, height);

}

void drawSpectrum1Plot(cairo_t* cr, int width, int height, LWEGui& theGui) {
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    if (!theGui.theSim.base().isGridAllocated || !(guiLock.owns_lock())) {
        blackoutCairoPlot(cr,width,height);
        return;
    }
    LwePlot<double> sPlot;

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
    if(!(theGui.checkboxes["plotDark"]->isChecked())){
        sPlot.backgroundColor =  LweColor(1,1,1,0);
        sPlot.axisColor = LweColor(0.0, 0.0, 0.0, 0);
        sPlot.textColor = sPlot.axisColor;
    }
    sPlot.makeSVG = theGui.isMakingSVG;
    sPlot.fontSize = theGui.textBoxes["fontSize"]->text().toDouble();
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theGui.theSim.base().fStep / 1e12;
    sPlot.yLabel = "Sx (J/THz)";
    if (theGui.checkboxes["Total"]->isChecked()) {
        sPlot.data.push_back(&theGui.theSim.base().totalSpectrum[(2 + simIndex * 3) * theGui.theSim.base().Nfreq]);
        sPlot.lineColors.push_back(LweColor(1.0, 0.5, 0.0, 0));
    }
    if (theGui.checkboxes["combinePolarizations"]->isChecked()) {
        sPlot.data.push_back(&theGui.theSim.base().totalSpectrum[(1 + simIndex * 3) * theGui.theSim.base().Nfreq]);
        sPlot.lineColors.push_back(LweColor(1, 0, 0.5, 0.0));
        sPlot.yLabel = "S (J/THz)";
    }
    sPlot.data.push_back(
        &theGui.theSim.base().totalSpectrum[simIndex * 3 * theGui.theSim.base().Nfreq]);
    sPlot.Npts = theGui.theSim.base().Nfreq;
    sPlot.logScale = theGui.checkboxes["Log"]->isChecked();;
    sPlot.lineColors.push_back(LweColor(0.5, 0, 1, 1));
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";

    sPlot.unitY = 1.0e-12;
    if(forceX){
        sPlot.forcedXmax = xMax;
        sPlot.forcedXmin = xMin;
    }
    if (forceY){
        sPlot.forcedYmin = yMin;
        sPlot.forcedYmax = yMax;
    }
    std::shared_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG) {
            theGui.SVGStrings[6] = sPlot.SVGString;
        }
    }
    else blackoutCairoPlot(cr, width, height);
}

void drawSpectrum2Plot(cairo_t* cr, int width, int height, LWEGui& theGui) {
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    if (!theGui.theSim.base().isGridAllocated || !(guiLock.owns_lock())) {
        blackoutCairoPlot(cr,width,height);
        return;
    }
    LwePlot<double> sPlot;

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
    if(!(theGui.checkboxes["plotDark"]->isChecked())){
        sPlot.backgroundColor =  LweColor(1,1,1,0);
        sPlot.axisColor = LweColor(0.0, 0.0, 0.0, 0);
        sPlot.textColor = sPlot.axisColor;
    }
    sPlot.makeSVG = theGui.isMakingSVG;
    sPlot.fontSize = theGui.textBoxes["fontSize"]->text().toDouble();
    sPlot.height = height;
    sPlot.width = width;
    sPlot.dx = theGui.theSim.base().fStep / 1e12;
    if (theGui.checkboxes["Total"]->isChecked()) {
        sPlot.data.push_back(&theGui.theSim.base().totalSpectrum[(2 + simIndex * 3) * theGui.theSim.base().Nfreq]);
        sPlot.lineColors = {LweColor(1.0, 0.5, 0.0, 0), LweColor(1, 0, 0.5, 0.0)};
    }
    else{
        sPlot.lineColors = {LweColor(1, 0, 0.5, 0.0)};
    }
    sPlot.data.push_back(&theGui.theSim.base().totalSpectrum[(1 + simIndex * 3) * theGui.theSim.base().Nfreq]);
    sPlot.Npts = theGui.theSim.base().Nfreq;
    sPlot.logScale = theGui.checkboxes["Log"]->isChecked();
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "Sy (J/THz)";
    sPlot.unitY = 1.0e-12;
    if(forceX){
        sPlot.forcedXmax = xMax;
        sPlot.forcedXmin = xMin;
    }
    if (forceY){
        sPlot.forcedYmin = yMin;
        sPlot.forcedYmax = yMax;
    }

    std::shared_lock dataLock(theGui.theSim.mutexes.at(simIndex), std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG) {
            theGui.SVGStrings[7] = sPlot.SVGString;
        }
    }
    else blackoutCairoPlot(cr, width, height);
}

void drawTimeImage2(cairo_t* cr, int width, int height, LWEGui& theGui) {
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    if (!theGui.theSim.base().isGridAllocated || !(guiLock.owns_lock())) {
        blackoutCairoPlot(cr,width,height);
        return;
    }

    int64_t simIndex = maxN(0,theGui.slider->value());
    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    int64_t cubeMiddle = theGui.theSim.base().Ntime * theGui.theSim.base().Nspace * (theGui.theSim.base().Nspace2 / 2);

    LweImage<double> image;
    image.data = DataSlice<double,2>{
        &theGui.theSim.base().ExtOut[simIndex * theGui.theSim.base().Ngrid * 2 + theGui.theSim.base().Ngrid + cubeMiddle],
        {static_cast<size_t>(theGui.theSim.base().Ntime),
        static_cast<size_t>(theGui.theSim.base().Nspace)}};
    image.height = height;
    image.width = width;
    image.colorMap = ColorMap::cyan_magenta;


    LwePlot<double> sPlot;
    if(!(theGui.checkboxes["plotDark"]->isChecked())){
        sPlot.backgroundColor =  LweColor(1,1,1,0);
        sPlot.axisColor = LweColor(0.0, 0.0, 0.0, 0);
        sPlot.textColor = sPlot.axisColor;
    }
    sPlot.fontSize = theGui.textBoxes["fontSize"]->text().toDouble();
    sPlot.height = height;
    sPlot.width = width;
    sPlot.image = &image;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Time (fs)";
    sPlot.yLabel = "x (mm)";
    sPlot.unitY = 1e3;
    sPlot.forcedXmax = 0.5e15 * theGui.theSim.base().Ntime * theGui.theSim.base().tStep;
    sPlot.forcedXmin = -sPlot.forcedXmax.value();
    sPlot.forcedYmin = 0.5e3 * theGui.theSim.base().spatialWidth;
    sPlot.forcedYmax = -sPlot.forcedYmin.value();
    sPlot.makeSVG = theGui.isMakingSVG;
    std::shared_lock dataLock(theGui.theSim.mutexes.at(simIndex),std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG){
            theGui.SVGStrings[1] = sPlot.SVGString;
        }
    }
    else blackoutCairoPlot(cr, width, height);
}

void drawFourierImage1(cairo_t* cr, int width, int height, LWEGui& theGui) {
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    if (!theGui.theSim.base().isGridAllocated || !(guiLock.owns_lock())) {
        blackoutCairoPlot(cr,width,height);
        return;
    }
    LweImage<double> image;
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
    image.complexData =
        DataSlice<std::complex<double>,2>{
            &theGui.theSim.base().EkwOut[simIndex * theGui.theSim.base().NgridC * 2],
            {static_cast<size_t>(theGui.theSim.base().Nfreq),
            static_cast<size_t>(theGui.theSim.base().Nspace)}};
    image.height = height;
    image.width = width;
    image.colorMap = ColorMap::purple;
    image.logMin = logPlotOffset;
    LwePlot<double> sPlot;
    if(!(theGui.checkboxes["plotDark"]->isChecked())){
        sPlot.backgroundColor =  LweColor(1,1,1,0);
        sPlot.axisColor = LweColor(0.0, 0.0, 0.0, 0);
        sPlot.textColor = sPlot.axisColor;
    }
    sPlot.fontSize = theGui.textBoxes["fontSize"]->text().toDouble();
    sPlot.height = height;
    sPlot.width = width;
    sPlot.image = &image;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "k\xe2\x82\x93 (1/\xce\xbcm)";
    sPlot.unitY = 1;
    sPlot.forcedXmax = 1e-12*(theGui.theSim.base().Nfreq-1) * theGui.theSim.base().fStep;
    sPlot.forcedXmin = 0;
    sPlot.forcedYmin = -0.5e-6 * theGui.theSim.base().Nspace * theGui.theSim.base().kStep;
    sPlot.forcedYmax = -sPlot.forcedYmin.value();
    sPlot.makeSVG = theGui.isMakingSVG;
    std::shared_lock dataLock(theGui.theSim.mutexes.at(simIndex),std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG){
            theGui.SVGStrings[2] = sPlot.SVGString;
        }
    }
    else blackoutCairoPlot(cr, width, height);
}

void drawFourierImage2(cairo_t* cr, int width, int height, LWEGui& theGui) {
    std::shared_lock guiLock(theGui.m, std::try_to_lock);
    if (!theGui.theSim.base().isGridAllocated || !(guiLock.owns_lock())) {
        blackoutCairoPlot(cr,width,height);
        return;
    }
    LweImage<double> image;

    int64_t simIndex = maxN(0,theGui.slider->value());
    if (simIndex > theGui.theSim.base().Nsims * theGui.theSim.base().Nsims2) {
        simIndex = 0;
    }

    double logPlotOffset = (double)(1e-4 / (theGui.theSim.base().spatialWidth * theGui.theSim.base().timeSpan));
    if (theGui.theSim.base().is3D) {
        logPlotOffset = (double)(1e-4
            / (theGui.theSim.base().spatialWidth * theGui.theSim.base().spatialHeight * theGui.theSim.base().timeSpan));
    }
    image.complexData = DataSlice<std::complex<double>,2>{
        &theGui.theSim.base().EkwOut[simIndex * theGui.theSim.base().NgridC * 2 + theGui.theSim.base().NgridC],
        {static_cast<size_t>(theGui.theSim.base().Nfreq),
        static_cast<size_t>(theGui.theSim.base().Nspace)}};
    image.height = height;
    image.width = width;
    image.colorMap = ColorMap::purple;
    image.logMin = logPlotOffset;
    LwePlot<double> sPlot;
    if(!(theGui.checkboxes["plotDark"]->isChecked())){
        sPlot.backgroundColor =  LweColor(1,1,1,0);
        sPlot.axisColor = LweColor(0.0, 0.0, 0.0, 0);
        sPlot.textColor = sPlot.axisColor;
    }
    sPlot.fontSize = theGui.textBoxes["fontSize"]->text().toDouble();
    sPlot.height = height;
    sPlot.width = width;
    sPlot.image = &image;
    sPlot.axisColor = LweColor(0.8, 0.8, 0.8, 0);
    sPlot.xLabel = "Frequency (THz)";
    sPlot.yLabel = "k\xe2\x82\x93 (1/\xce\xbcm)";
    sPlot.unitY = 1;
    sPlot.forcedXmax = 1e-12*(theGui.theSim.base().Nfreq-1) * theGui.theSim.base().fStep;
    sPlot.forcedXmin = 0;
    sPlot.forcedYmin = -0.5e-6 * theGui.theSim.base().Nspace * theGui.theSim.base().kStep;
    sPlot.forcedYmax = -sPlot.forcedYmin.value();
    sPlot.makeSVG = theGui.isMakingSVG;
    std::shared_lock dataLock(theGui.theSim.mutexes.at(simIndex),std::try_to_lock);
    if (dataLock.owns_lock()) {
        sPlot.plot(cr);
        if(theGui.isMakingSVG){
            theGui.SVGStrings[3] = sPlot.SVGString;
        }
    }
    else blackoutCairoPlot(cr, width, height);
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

int insertAfterCharacterExcept(
    std::string& s,
    char target,
    std::string appended,
    std::string exclude){
    bool match = false;
    for(std::size_t i = 0; i < s.length()-1; ++i){
        if(s[i] == target){
            match = false;
            for(std::size_t j = 0; j<exclude.length(); j++){
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
    std::size_t depth = 0;
    std::string indent("   ");
    for (std::size_t i = 0; i < s.length() - 1; ++i) {
        if (s[i] == '{') ++depth;
        if (s[i] == '}' && depth != 0) --depth;
        if (s[i] == '\n' && s[i + 1] != '}') {
            for (std::size_t j = 0; j < depth; j++) {
                s.insert(i + 1, indent);
                i += indent.length();
            }
        }
    }
    return 0;
}

int formatSequence(std::string& s){
    if(s.size()==0) return 0;
    insertAfterCharacter(s, '>', std::string("\n"));
    insertAfterCharacterExcept(s, ')', std::string("\n"), std::string("^{;"));
    insertAfterCharacter(s, ';', std::string("\n"));
    insertAfterCharacter(s, '{', std::string("\n"));
    insertAfterCharacter(s, '}', std::string("\n"));
    indentForDepth(s);
    return 0;
}

#include "LightwaveExplorerFrontendQt.moc"
