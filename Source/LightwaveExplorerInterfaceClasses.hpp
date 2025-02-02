#pragma once
#include "DataStructures.hpp"
#include "LightwaveExplorerUtilities.h"
#include <atomic>

class loadedInputData {
    public:
    std::string fileContents;
    std::string filePath;
    bool hasData = false;
    loadedInputData(){}
    loadedInputData(const std::string& path){
        filePath = path;
        std::ifstream file(filePath);
        if(file.fail()) return;
        fileContents = std::string((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        hasData = fileContents.size() > 1;
    }


    //method to output the saved data to a complex spectrum (e.g. reflectivity or transmission)
    //assuming the data format is: wavelength (nm) | R (0.0 to 1.0) | phase(rad)
    template<typename deviceFP>
    std::vector<std::complex<deviceFP>> toComplexSpectrum(const int64_t Nfreq, const double fStep) {
        std::vector<std::complex<deviceFP>> complexReflectivity(Nfreq, std::complex<deviceFP>(0.0,0.0));
        if (!hasData) {
            return complexReflectivity;
        }

        std::stringstream fs(fileContents);
        int64_t maxFileSize = 16384;
        int64_t currentRow = 0;
        constexpr double c = 1e9 * lightC<double>();
        struct Reflectivity{
            double f;
            std::complex<double> R;
        };

        std::vector<Reflectivity> loadedReflectivities;
        loadedReflectivities.reserve(8192);
        loadedReflectivities.push_back({0.0,{0.0,0.0}});
        double minFrequency{};
        double maxFrequency{};
        while (fs.good() && currentRow < maxFileSize) {
            double wavelength;
            double R;
            double phase;
            fs >> wavelength >> R >> phase;
            double f = c/wavelength;
            //store in complex as {R, phase}
            loadedReflectivities.push_back({
                f,
                std::complex<double>(R,phase)});
            if (currentRow == 0) {
                minFrequency = f;
                maxFrequency = f;
            }
            else {
                maxFrequency = maxN(maxFrequency, f);
                minFrequency = minN(minFrequency, f);
            }
            currentRow++;
        }

        std::sort(
            loadedReflectivities.begin(),
            loadedReflectivities.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.f < rhs.f;
            }
        );

        double currentFrequency = 0;
        double df;
        for (int64_t i = 1; i < Nfreq; i++) {
            currentFrequency = i * fStep;
            if ((currentFrequency > minFrequency) && (currentFrequency < maxFrequency)) {
                //find the first frequency greater than the current value
                int64_t j = 0;
                while (loadedReflectivities[j].f < currentFrequency) {
                    j++;
                }
                //linear interpolation
                df = loadedReflectivities[j].f - loadedReflectivities[j-1].f;
                double t = (currentFrequency - loadedReflectivities[j-1].f)/df;
                std::complex<double> currentValue = 
                    (1.0 - t)* loadedReflectivities[j-1].R 
                    + t * loadedReflectivities[j].R;

                //put in complex representation
                currentValue = 
                    std::sqrt(std::abs(currentValue.real())) * 
                    std::exp(std::complex<double>(0.0,currentValue.imag()));
                complexReflectivity[i] = std::complex<deviceFP>(
                    static_cast<deviceFP>(currentValue.real()),
                    static_cast<deviceFP>(currentValue.imag()));
            }
            else{
                complexReflectivity[i] = {};
            }
        }
        return complexReflectivity;
    }

    loadedInputData(const std::string& zipPath, const std::string& filename){
        if(zipContainsFile(zipPath, filename)){
            std::vector<char> data;
            zipIntoMemory(zipPath,filename,data);
            if(data.size()<1) return;
            fileContents = std::string(data.begin(),data.end());
            filePath = filename;
            hasData = true;
        }
    }
};


//Simulation parameter class containing the complete description of the running simulation
//intended only to be present on the CPU, as certain contents (std::array, std::string) can not
//be assumed to be implemented on the device.
class simulationParameterSet {
public:
    double rStep = 0;
    double tStep = 0;
    double fStep = 0;
    double kStep = 0;
    double propagationStep = 0;
    int64_t Npropagation = 0;
    int64_t Ntime = 0;
    int64_t Nfreq = 0;
    int64_t Nspace = 0;
    int64_t Nspace2 = 0;
    int64_t Ngrid = 0;
    int64_t NgridC = 0;
    int64_t Nsims = 0;
    int64_t Nsims2 = 0;
    std::atomic_uint32_t* progressCounter = 0;
    int64_t NsimsCPU = 0;
    pulse<double> pulse1;
    pulse<double> pulse2;
    double spatialWidth = 0;
    double spatialHeight = 0;
    double timeSpan = 0;
    int materialIndex = 0;
    int materialIndexAlternate = 0;
    double bandGapElectronVolts = 0;
    double effectiveMass = 0;
    double nonlinearAbsorptionStrength = 0;
    double startingCarrierDensity = 0;
    double drudeGamma = 0;
    double crystalTheta = 0;
    double crystalPhi = 0;
    double crystalThickness = 0;
    double* chi2Tensor = 0;
    double* chi3Tensor = 0;
    double* sellmeierCoefficients = 0;
    int sellmeierType = 0;
    int axesNumber = 0;
    NonlinearPropertyFlags nonlinearSwitches = {};
    bool isCylindric = 0;
    bool is3D = 0;
    bool isFDTD = 0;
    int symmetryType = 0;
    bool useOpenMP = true;
    //loaded FROG/EOS fields
    std::complex<double>* loadedField1 = 0;
    std::complex<double>* loadedField2 = 0;
    double* loadedFullGrid1 = 0;
    double* loadedFullGrid2 = 0;
    bool field1IsAllocated = 0;
    bool field2IsAllocated = 0;
    int pulse1FileType = 0;
    int pulse2FileType = 0;
    loadedInputData pulse1LoadedData;
    loadedInputData pulse2LoadedData;
    loadedInputData fittingLoadedData;
    int pulsetype = 0;
    double* ExtOut = 0;
    std::complex<double>* EkwOut = 0;
    double* totalSpectrum = 0;
    int memoryError = 0;
    int assignedGPU = 0;
    int plotSim = 0;
    crystalEntry* crystalDatabase = 0;
    int batchIndex = 0;
    int batchIndex2 = 0;
    double batchDestination = 0;
    double batchDestination2 = 0;
    std::string outputBasePath;
    runTypes runType = runTypes::normal;
    bool runningOnCPU = 0;

    //sequence
    bool isInSequence = 0;
    bool isFollowerInSequence = 0;
    bool isReinjecting = 0;
    bool forceLinear = 0;
    std::string sequenceString;
    double i37 = 0.0;
    int64_t batchLoc1 = 0;
    int64_t batchLoc2 = 0;
    std::vector<loadedInputData> optics{};

    //fitting
    bool isInFittingMode = false;
    std::string fittingString;
    std::array<double, 256> fittingArray = {};
    double* fittingReference = 0;
    int Nfitting = 0;
    int fittingMode = 0;
    int fittingMaxIterations = 0;
    int64_t fittingROIstart = 0;
    int64_t fittingROIstop = 0;
    int64_t fittingROIsize = 0;
    std::array<double, 64> fittingResult = {};
    std::array<double, 64> fittingError = {};

    //Status
    bool isRunning = false;
    bool isGridAllocated = false;
    bool cancellationCalled = false;
    bool CUDAavailable = false;
    bool SYCLavailable = false;
    int cudaGPUCount = 0;
    int syclGPUCount = 0;

	std::array<double, 38> multipliers = { 0,
        1, 1, 1e12, 1e12,
        1e12, 1e12, vPi<double>(), vPi<double>(),
        1e-15, 1e-15, 1e-30, 1e-30,
        1e-45, 1e-45, 1e-6, 1e-6,
        1e-6, 1e-6,
        1e-6, 1e-6, 1e-6, 1e-6,
        deg2Rad<double>(), deg2Rad<double>(), deg2Rad<double>(), deg2Rad<double>(),
        1, 1, deg2Rad<double>(), deg2Rad<double>(),
        1, 1e12, 1, 1e-6,
        1e-9, 1, 1 };

    [[nodiscard]] constexpr double getByNumberWithMultiplier(const std::size_t index) {
        if (index == 0 || index == 36 || index >= multipliers.size()) return 0.0;
        return  getByNumber(index) / multipliers[index];
    }
    
    constexpr double getByNumber(const std::size_t index) {
        switch (index) {
        case 0:
            return 0.0;
        case 1:
            return pulse1.energy;
        case 2:
            return pulse2.energy;
        case 3:
            return pulse1.frequency;
        case 4:
            return pulse2.frequency;
        case 5:
            return pulse1.bandwidth;
        case 6:
            return pulse2.bandwidth;
        case 7:
            return pulse1.cep;
        case 8:
            return pulse2.cep;
        case 9:
            return pulse1.delay;
        case 10:
            return pulse2.delay;
        case 11:
            return pulse1.gdd;
        case 12:
            return pulse2.gdd;
        case 13:
            return pulse1.tod;
        case 14:
            return pulse2.tod;
        case 15:
            return pulse1.phaseMaterialThickness;
        case 16:
            return pulse2.phaseMaterialThickness;
        case 17:
            return pulse1.beamwaist;
        case 18:
            return pulse2.beamwaist;
        case 19:
            return pulse1.x0;
        case 20:
            return pulse2.x0;
        case 21:
            return pulse1.z0;
        case 22:
            return pulse2.z0;
        case 23:
            return pulse1.beamAngle;
        case 24:
            return pulse2.beamAngle;
        case 25:
            return pulse1.polarizationAngle;
        case 26:
            return pulse2.polarizationAngle;
        case 27:
            return pulse1.circularity;
        case 28:
            return pulse2.circularity;
        case 29:
            return crystalTheta;
        case 30:
            return crystalPhi;
        case 31:
            return nonlinearAbsorptionStrength;
        case 32:
            return drudeGamma;
        case 33:
            return effectiveMass;
        case 34:
            return crystalThickness;
        case 35:
            return propagationStep;
        case 36:
            return 0.0;
        case 37:
            return i37;
        default:
            return 0.0;
        };
    }
    void setByNumber(const int64_t index, const double value);
    void setByNumberWithMultiplier(const std::size_t index, const double value);
    int loadSavedFields(const std::string& outputBase, bool isZipFile);
    int loadReferenceSpectrum();
    int readInputParametersFile(crystalEntry* crystalDatabasePtr, const std::string filePath);
    std::string settingsString();
    int saveSettingsFile();
    double saveSlurmScript(const std::string& gpuType, int gpuCount, bool useJobArray, int64_t totalSteps, std::vector<simulationParameterSet>& params, const class crystalDatabase& db);
    int readFittingString();


    template <typename deviceFP, typename C>
    void initializeDeviceParameters(deviceParameterSet<deviceFP, C>* s) {
        (*s).Ntime = Ntime;
        (*s).Nspace = Nspace;
        (*s).Nspace2 = Nspace2;
        (*s).is3D = is3D;
        (*s).Nfreq = ((*s).Ntime / 2 + 1);
        (*s).Ngrid = (*s).Ntime * (*s).Nspace * (*s).Nspace2;
        (*s).NgridC = (*s).Nfreq * (*s).Nspace * (*s).Nspace2; //size of the positive frequency side of the grid
        (*s).fftNorm = static_cast<deviceFP>(1.0 / (*s).Ngrid);
        (*s).dt = static_cast<deviceFP>(tStep);
        (*s).dx = static_cast<deviceFP>(rStep);
        (*s).dk1 = static_cast<deviceFP>(twoPi<double>() / (Nspace * rStep));
        (*s).dk2 = static_cast<deviceFP>(twoPi<double>() / (Nspace2 * rStep));
        (*s).fStep = static_cast<deviceFP>(fStep);
        (*s).Nsteps = static_cast<int64_t>(round(crystalThickness / propagationStep));
        (*s).h = static_cast<deviceFP>(crystalThickness / ((*s).Nsteps)); //adjust step size so that thickness can be varied continuously by fitting
        (*s).axesNumber = axesNumber;
        (*s).sellmeierType = sellmeierType;
        (*s).crystalPhi = static_cast<deviceFP>(crystalPhi);
        (*s).crystalTheta = static_cast<deviceFP>(crystalTheta);
        (*s).f0 = static_cast<deviceFP>(pulse1.frequency);
        (*s).Nthread = threadsPerBlock;
        (*s).Nblock = static_cast<int>((*s).Ngrid / threadsPerBlock);
        (*s).NblockC = static_cast<int>((*s).NgridC / threadsPerBlock);
        (*s).isCylindric = isCylindric;
        (*s).forceLinear = forceLinear;
        (*s).isNonLinear = (nonlinearSwitches.hasChi2 || nonlinearSwitches.hasChi3);
        (*s).isUsingMillersRule = (crystalDatabase[materialIndex].nonlinearReferenceFrequencies[0]) != 0.0;

        if (nonlinearAbsorptionStrength > 0. || startingCarrierDensity > 0) {
            (*s).hasPlasma = true;
            (*s).isNonLinear = true;
        }
        else {
            (*s).hasPlasma = false;
        }

        if ((*s).forceLinear) {
            (*s).hasPlasma = false;
            (*s).isNonLinear = false;
        }
    }

    template <typename deviceFP, typename C>
    void fillRotationMatricies(deviceParameterSet<deviceFP, C>* s) {
        double cosT = cos(crystalTheta);
        double sinT = sin(crystalTheta);
        double cosP = cos(crystalPhi);
        double sinP = sin(crystalPhi);
        double forward[9] =
        { cosT * cosP, sinP, -sinT * cosP, -sinP * cosT, cosP, sinP * sinT, sinT, 0, cosT };

        //reverse direction (different order of operations)
        double backward[9] =
        { cosT * cosP, -sinP * cosT, sinT, sinP, cosP, 0, -sinT * cosP, sinP * sinT, cosT };

        for (int64_t i = 0; i < 9; i++) {
            (*s).rotationForward[i] = (deviceFP)forward[i];
            (*s).rotationBackward[i] = (deviceFP)backward[i];
        }
    }

    template<typename deviceFP, typename C>
    void finishConfiguration(deviceParameterSet<deviceFP,C>* s) {
        int64_t beamExpansionFactor = 1;
        if ((*s).isCylindric) {
            beamExpansionFactor = 2;
        }
        //second polarization grids are to pointers within the first polarization
        //to have contiguous memory
        (*s).gridETime2 = (*s).gridETime1 + (*s).Ngrid;
        (*s).workspace2 = (*s).workspace1 + (*s).NgridC;
        (*s).gridPolarizationTime2 = (*s).gridPolarizationTime1 + (*s).Ngrid;
        (*s).workspace2P = (*s).workspace1 + beamExpansionFactor * (*s).NgridC;
        (*s).k2 = (*s).k1 + (*s).NgridC;
        (*s).chiLinear2 = (*s).chiLinear1 + (*s).Nfreq;
        (*s).fieldFactor2 = (*s).fieldFactor1 + (*s).Nfreq;
        (*s).inverseChiLinear2 = (*s).inverseChiLinear1 + (*s).Nfreq;
        (*s).gridRadialLaplacian2 = (*s).gridRadialLaplacian1 + (*s).Ngrid;
        (*s).gridPropagationFactor1Rho2 = (*s).gridPropagationFactor1Rho1 + (*s).NgridC;
        (*s).gridPolarizationFactor2 = (*s).gridPolarizationFactor1 + (*s).NgridC;
        (*s).gridEFrequency1Next2 = (*s).gridEFrequency1Next1 + (*s).NgridC;
        (*s).gridPropagationFactor2 = (*s).gridPropagationFactor1 + (*s).NgridC;
        (*s).gridEFrequency2 = (*s).gridEFrequency1 + (*s).NgridC;

        PlasmaParameters<double> plasmaParameters = {};
        plasmaParameters.fieldExponent = static_cast<int>(ceil(bandGapElectronVolts * 241.79893e12 / pulse1.frequency) - 2);
        plasmaParameters.nonlinearAbsorption = nonlinearAbsorptionStrength; //nonlinear absorption strength parameter
        plasmaParameters.drudeGamma = drudeGamma; //gamma
        plasmaParameters.initialDensity = startingCarrierDensity;
        if (nonlinearAbsorptionStrength > 0.) {
            plasmaParameters.integrationFactor = tStep * tStep
                * 2.817832e-08 / (1.6022e-19 * bandGapElectronVolts * effectiveMass); // (dt^2)*e* e / (m * band gap));
        }
        else {
            plasmaParameters.integrationFactor = 0.0;
        }

        for (int j = 0; j < 18; ++j) {
            (*s).chi2Tensor[j] = static_cast<deviceFP>(2e-12 * chi2Tensor[j]); //go from d in pm/V to chi2 in m/V
            if (j > 8) (*s).chi2Tensor[j] *= 2.0; //multiply cross-terms by 2 for consistency with convention
        }

        (*s).nonlinearSwitches = nonlinearSwitches;

        for (int64_t i = 0; i < 81; i++) {
            (*s).chi3Tensor[i] = static_cast<deviceFP>(chi3Tensor[i]);
        }


        (*s).plasmaParameters = plasmaParameters;
    }
};

class simulationBatch {
    std::vector<double> Ext;
    std::vector<std::complex<double>> Ekw;
    std::vector<std::complex<double>> loadedField1;
    std::vector<std::complex<double>> loadedField2;
    std::vector<double> loadedFullGrid1;
    std::vector<double> loadedFullGrid2;
    std::vector<double> fitReference;
    std::vector<double> totalSpectrum;
    int64_t Nsimstotal = 0;
    int64_t Nsims = 0;
    int64_t Nsims2 = 0;
    int64_t Nfreq = 0;
    int64_t Ngrid = 0;
    int64_t NgridC = 0;
public:
    std::vector<loadedInputData> optics;
    std::vector<simulationParameterSet> parameters;
    std::vector<std::mutex> mutexes = std::vector<std::mutex>(1);

    simulationBatch() {
        parameters = std::vector<simulationParameterSet>(1);
    }
    ~simulationBatch() {
        std::for_each(mutexes.begin(), mutexes.end(), 
            [](std::mutex& m) {std::lock_guard<std::mutex> lock(m); });
    }
    void configure(bool allocateFields=true);
    void configureCounter();
    void loadPulseFiles();
    void loadOptics(const std::string& zipPath);
    int saveDataSet();
    [[nodiscard]] double* getExt(int64_t i) {
        return &Ext.data()[i * Ngrid * 2];
    }
    [[nodiscard]] std::complex<double>* getEkw(int64_t i) {
        return &Ekw.data()[i * NgridC * 2];
    }
    [[nodiscard]] double* getTotalSpectrum(int64_t i) {
        return &totalSpectrum.data()[i * 3 * Nfreq];
    }
    [[nodiscard]] std::vector<simulationParameterSet>& getParameterVector() {
        return parameters;
    }
    [[nodiscard]] simulationParameterSet* sCPU() {
        return parameters.data();
    }

    [[nodiscard]] simulationParameterSet& base() {
        return parameters[0];
    }
};

int				loadFrogSpeck(const loadedInputData& frogFilePath, std::complex<double>* Egrid, const int64_t Ntime, const double fStep, const double gateLevel);
int             loadWaveformFile(const loadedInputData& waveformFile, std::complex<double>* outputGrid, const int64_t Ntime, const double fStep);
int             loadSavedGridFile(const loadedInputData& file, std::vector<double>& outputGrid, int64_t Ngrid);
int             loadSavedGridFileMultiple(const loadedInputData& file, std::vector<double>& outputGrid, int64_t Ngrid, int64_t Nsims);