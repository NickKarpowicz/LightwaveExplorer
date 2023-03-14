#pragma once
#define maxN(a,b)            (((a) > (b)) ? (a) : (b))
#define minN(a,b)            (((a) < (b)) ? (a) : (b))
#include <complex>
#include <cstring>
#include <vector>
#include <fstream>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#define MAX_LOADSTRING 1024

#ifndef LWEFLOATINGPOINT
#define LWEFLOATINGPOINT 32
#endif
#if LWEFLOATINGPOINT==32
#define LWEFLOATINGPOINTTYPE float
#define ANGLETOLERANCE 1e-12f
#define TWOPI 6.2831853071795862f
#define PI 3.1415926535897931f
#define INVSQRTPI 0.5641895835477563f
#define DEG2RAD 1.7453292519943295e-02f
#define RAD2DEG 57.2957795130823229f
#define LIGHTC 2.99792458e8f
#define EPS0 8.8541878128e-12f
#define SIXTH 0.1666666666666667f
#define THIRD 0.3333333333333333f
#define KLORENTZIAN 3182.607353999257f //(e * e / (epsilon_o * m_e)
#else
#define LWEFLOATINGPOINTTYPE double
#define ANGLETOLERANCE 1e-12
#define TWOPI 6.2831853071795862
#define PI 3.1415926535897931
#define INVSQRTPI 0.5641895835477563
#define DEG2RAD 1.7453292519943295e-02
#define RAD2DEG 57.2957795130823229
#define LIGHTC 2.99792458e8
#define EPS0 8.8541878128e-12
#define SIXTH 0.1666666666666667
#define THIRD 0.3333333333333333
#define KLORENTZIAN 3182.607353999257 //(e * e / (epsilon_o * m_e)
#endif
#include <complex>
#include <cstring>
#ifdef __CUDACC__
#include <cufft.h>
#include <thrust/complex.h>
#include <fftw3_mkl.h>
typedef thrust::complex<LWEFLOATINGPOINTTYPE> deviceComplex;
typedef LWEFLOATINGPOINTTYPE deviceFP;
#elif defined RUNONSYCL
#include <oneapi/dpl/complex>
#include <oneapi/dpl/cmath>
typedef oneapi::dpl::complex<LWEFLOATINGPOINTTYPE> deviceComplex;
typedef LWEFLOATINGPOINTTYPE deviceFP;
#include <sycl/sycl.hpp>
#elif defined CPUONLY
#include <fftw3.h>
typedef std::complex<LWEFLOATINGPOINTTYPE> deviceComplex;
typedef LWEFLOATINGPOINTTYPE deviceFP;
#else
#include <fftw3_mkl.h>
typedef std::complex<LWEFLOATINGPOINTTYPE> deviceComplex;
typedef LWEFLOATINGPOINTTYPE deviceFP;
#endif

#define THREADS_PER_BLOCK 32
#define MIN_GRIDDIM 8
#define FALSE 0
#define TRUE 1
#define MAX_LOADSTRING 1024


class deviceParameterSet {
public:
    deviceComplex* workspace1 = 0;
    deviceComplex* workspace2 = 0;
    deviceComplex* workspace2P = 0;
    deviceComplex* gridETemp1 = 0;
    deviceComplex* gridETemp2 = 0;
    deviceComplex* gridEFrequency1 = 0;
    deviceComplex* gridEFrequency2 = 0;
    deviceComplex* gridPropagationFactor1 = 0;
    deviceComplex* gridPropagationFactor1Rho1 = 0;
    deviceComplex* gridPropagationFactor1Rho2 = 0;
    deviceComplex* gridPolarizationFactor1 = 0;
    deviceComplex* gridPolarizationFrequency1 = 0;
    deviceComplex* gridPropagationFactor2 = 0;
    deviceComplex* gridPolarizationFactor2 = 0;
    deviceComplex* gridPolarizationFrequency2 = 0;
    deviceComplex* gridEFrequency1Next1 = 0;
    deviceComplex* gridEFrequency1Next2 = 0;
    deviceComplex* gridPlasmaCurrentFrequency1 = 0;
    deviceComplex* gridPlasmaCurrentFrequency2 = 0;
    deviceComplex* chiLinear1 = 0;
    deviceComplex* chiLinear2 = 0;
    deviceFP* inverseChiLinear1 = 0;
    deviceFP* inverseChiLinear2 = 0;
    deviceFP* fieldFactor1 = 0;
    deviceFP* fieldFactor2 = 0;
    deviceComplex* k1 = 0;
    deviceComplex* k2 = 0;
    deviceComplex n0 = 0.0;
    deviceFP* gridRadialLaplacian1 = 0;
    deviceFP* gridRadialLaplacian2 = 0;
    deviceFP* gridETime1 = 0;
    deviceFP* gridETime2 = 0;
    deviceFP* gridPolarizationTime1 = 0;
    deviceFP* gridPolarizationTime2 = 0;
    deviceFP* expGammaT = 0;
    deviceFP* gridPlasmaCurrent1 = 0;
    deviceFP* gridPlasmaCurrent2 = 0;

    //fixed length arrays
    deviceFP firstDerivativeOperation[6] = { 0 };
    deviceFP plasmaParameters[6] = { 0 }; //[dt^2 * e^2/m * nonlinearAbsorptionStrength, gamma] 
    deviceFP chi2Tensor[18] = { 0 };
    deviceFP chi3Tensor[81] = { 0 };
    deviceFP absorptionParameters[6] = { 0 };
    deviceFP rotationForward[9] = { 0 };
    deviceFP rotationBackward[9] = { 0 };
    int nonlinearSwitches[4] = { 0 };

    bool isCylindric = 0;
    bool is3D = 0;
    bool hasPlasma = 0;
    bool isNonLinear = 0;
    bool isUsingMillersRule = 0;
    bool forceLinear = 0;
    size_t Ntime = 0;
    size_t Nfreq = 0;
    size_t Nspace = 0;
    size_t Nspace2 = 0;
    size_t Ngrid = 0;
    size_t NgridC = 0;
    deviceFP fftNorm = 0;
    int axesNumber = 0;
    int sellmeierType = 0;
    deviceFP crystalTheta;
    deviceFP crystalPhi;
    deviceFP f0 = 0;
    deviceFP fStep = 0;
    deviceFP dt = 0;
    deviceFP dx = 0;
    deviceFP dk1 = 0;
    deviceFP dk2 = 0;
    deviceFP h = 0;
    size_t Nsteps = 0;
    int Nthread = 0;
    int NblockC = 0;
    int Nblock = 0;
};


class crystalEntry {
public:
    std::string crystalName;
    int axisType = 0;
    int sellmeierType = 0;
    int nonlinearSwitches[4] = { 0 };
    double sellmeierCoefficients[66] = { 0 };
    std::string sellmeierReference;
    double d[18] = { 0 };
    std::string dReference;
    double chi3[81] = { 0 };
    std::string chi3Reference;
    double absorptionParameters[6] = { 0 };
    char spectralFile[512] = { 0 };
    double nonlinearReferenceFrequencies[7] = { 0 };
};

class crystalDatabase {
public:
    std::vector<crystalEntry> db;

    crystalDatabase() {
#ifdef __APPLE__
        uint32_t bufferSize = 1024;
        char sysPath[1024] = { 0 };
        _NSGetExecutablePath(sysPath, &bufferSize);
        std::string macPath(sysPath);
        size_t posPath = macPath.find_last_of("/");
        std::string databasePath = macPath.substr(0, posPath).append("/../Resources/CrystalDatabase.txt");
        std::ifstream fs(databasePath);
        if (!fs.is_open()) {
            fs.open("CrystalDatabase.txt");
        }
#elif defined __linux__
        std::ifstream fs("/usr/share/LightwaveExplorer/CrystalDatabase.txt");
        if (!fs.is_open()) {
            fs.open("CrystalDatabase.txt");
        }
#else
        std::ifstream fs("CrystalDatabase.txt");
#endif
        
        std::string line;
        if (!fs.is_open())return;
        while (!fs.eof() && fs.good()) {
            crystalEntry newEntry;
            std::getline(fs, line);//Name:

            std::getline(fs, line);
            newEntry.crystalName = line;

            std::getline(fs, line); //Type:
            fs >> newEntry.axisType;
            std::getline(fs, line);

            std::getline(fs, line); //Sellmeier equation:
            fs >> newEntry.sellmeierType;
            std::getline(fs, line);

            std::getline(fs, line); //1st axis coefficients:
            for (int k = 0; k < 22; ++k) {
                fs >> newEntry.sellmeierCoefficients[k];
            }
            std::getline(fs, line);

            std::getline(fs, line); //2nd axis coefficients:
            for (int k = 0; k < 22; ++k) {
                fs >> newEntry.sellmeierCoefficients[k + 22];
            }
            std::getline(fs, line);

            std::getline(fs, line); //3rd axis coefficients:
            for (int k = 0; k < 22; ++k) {
                fs >> newEntry.sellmeierCoefficients[k + 44];
            }
            std::getline(fs, line);

            std::getline(fs, line); //Sellmeier reference:
            std::getline(fs, line);
            newEntry.sellmeierReference = line;

            std::getline(fs, line); // chi2 type:
            fs >> newEntry.nonlinearSwitches[0];
            std::getline(fs, line);

            std::getline(fs, line); //d:
            for (int k = 0; k < 18; ++k) {
                fs >> newEntry.d[3 * (k % 6) + (k / 6)]; //column-order!
            }
            std::getline(fs, line);

            std::getline(fs, line); //d reference:
            std::getline(fs, line);
            newEntry.dReference = line;

            std::getline(fs, line); //chi3 type:
            fs >> newEntry.nonlinearSwitches[1];
            std::getline(fs, line);

            std::getline(fs, line); //chi3:
            memset(newEntry.chi3, 0, 81 * sizeof(double));
            switch (newEntry.nonlinearSwitches[1]) {
            case 0: //no chi3, skip all three lines
                std::getline(fs, line);
                std::getline(fs, line);
                std::getline(fs, line);
                break;
            case 1: //read full chi3
                for (int k = 0; k < 81; ++k) {
                    fs >> newEntry.chi3[k]; //row-order!
                }
                std::getline(fs, line);
                break;

            case 2: //assume centrosymmetric, just read chi3_1111, then skip
                fs >> newEntry.chi3[0];
                std::getline(fs, line);
                std::getline(fs, line);
                std::getline(fs, line);
                break;
            }

            std::getline(fs, line); //chi3 reference:
            std::getline(fs, line);
            newEntry.chi3Reference = line;

            std::getline(fs, line); //Spectral file:
            std::getline(fs, line);
            line.copy(newEntry.spectralFile, 512);

            std::getline(fs, line); //Nonlinear reference frequencies:
            for (int k = 0; k < 7; ++k) {
                fs >> newEntry.nonlinearReferenceFrequencies[k];
            }
            std::getline(fs, line);
            std::getline(fs, line); //~~~crystal end~~~
            if(fs.good())db.push_back(newEntry);
        }
    }
};

template <typename T>
class pulse {
public:
    T energy;
    T frequency;
    T bandwidth;
    int sgOrder;
    T cep;
    T delay;
    T gdd;
    T tod;
    int phaseMaterial;
    T phaseMaterialThickness;
    T beamwaist;
    T x0;
    T y0;
    T z0;
    T beamAngle;
    T polarizationAngle;
    T beamAnglePhi;
    T circularity;
    T pulseSum;

    pulse() : energy(), 
        frequency(),
        bandwidth(),
        sgOrder(),
        cep(),
        delay(),
        gdd(),
        tod(),
        phaseMaterial(),
        phaseMaterialThickness(),
        beamwaist(),
        x0(),
        y0(),
        z0(),
        beamAngle(),
        polarizationAngle(),
        beamAnglePhi(),
        circularity(),
        pulseSum(){}

    template<typename U>
    pulse(pulse<U>& other) : energy((T)other.energy),
        frequency((T)other.frequency),
        bandwidth((T)other.bandwidth),
        sgOrder(other.sgOrder),
        cep((T)other.cep),
        delay((T)other.delay),
        gdd((T)other.gdd),
        tod((T)other.tod),
        phaseMaterial(other.phaseMaterial),
        phaseMaterialThickness((T)other.phaseMaterialThickness),
        beamwaist((T)other.beamwaist),
        x0((T)other.x0),
        y0((T)other.y0),
        z0((T)other.z0),
        beamAngle((T)other.beamAngle),
        polarizationAngle((T)other.polarizationAngle),
        beamAnglePhi((T)other.beamAnglePhi),
        circularity((T)other.circularity),
        pulseSum((T)other.pulseSum) {}

    template <typename U>
    pulse& operator=(const pulse<U>& other) {
        energy = (T)other.energy;
        frequency = (T)other.frequency;
        bandwidth = (T)other.bandwidth;
        sgOrder = other.sgOrder;
        cep = (T)other.cep;
        delay = (T)other.delay;
        gdd = (T)other.gdd;
        tod = (T)other.tod;
        phaseMaterial = other.phaseMaterial;
        phaseMaterialThickness = (T)other.phaseMaterialThickness;
        beamwaist = (T)other.beamwaist;
        x0 = (T)other.x0;
        y0 = (T)other.y0;
        z0 = (T)other.z0;
        beamAngle = (T)other.beamAngle;
        polarizationAngle = (T)other.polarizationAngle;
        beamAnglePhi = (T)other.beamAnglePhi;
        circularity = (T)other.circularity;
        pulseSum = (T)other.pulseSum;
        return *this;
    }
};


//Simulation parameter struct to pass to the simulations running in threads
class simulationParameterSet {
public:
    double rStep = 0;
    double tStep = 0;
    double fStep = 0;
    double kStep = 0;
    double propagationStep = 0;
    size_t Npropagation = 0;
    size_t Ntime = 0;
    size_t Nfreq = 0;
    size_t Nspace = 0;
    size_t Nspace2 = 0;
    size_t Ngrid = 0;
    size_t NgridC = 0;
    size_t Nsims = 0;
    size_t Nsims2 = 0;
    size_t* progressCounter = 0;
    size_t NsimsCPU = 0;
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
    double drudeGamma = 0;
    double crystalTheta = 0;
    double crystalPhi = 0;
    double crystalThickness = 0;
    double* chi2Tensor = 0;
    double* deffTensor = 0;
    double* chi3Tensor = 0;
    double* sellmeierCoefficients = 0;
    double* absorptionParameters = 0;
    int sellmeierType = 0;
    int axesNumber = 0;
    int* nonlinearSwitches = 0;
    bool isCylindric = 0;
    bool is3D = 0;
    int symmetryType = 0;

    //loaded FROG/EOS fields
    std::complex<double>* loadedField1 = 0;
    std::complex<double>* loadedField2 = 0;
    bool field1IsAllocated = 0;
    bool field2IsAllocated = 0;
    int pulse1FileType = 0;
    int pulse2FileType = 0;
    char field1FilePath[MAX_LOADSTRING] = { 0 };
    char field2FilePath[MAX_LOADSTRING] = { 0 };

    int pulsetype = 0;
    double* ExtOut = 0;
    std::complex<double>* EkwOut = 0;
    double* totalSpectrum = 0;
    int* statusFlags = 0;
    int memoryError = 0;
    int assignedGPU = 0;
    int plotSim = 0;
    crystalEntry* crystalDatabase = 0;
    int batchIndex = 0;
    int batchIndex2 = 0;
    double batchDestination = 0;
    double batchDestination2 = 0;
    char outputBasePath[MAX_LOADSTRING] = { 0 };
    int runType = 0;
    bool runningOnCPU = 0;

    //sequence
    bool isInSequence = 0;
    bool isFollowerInSequence = 0;
    bool isReinjecting = 0;
    bool forceLinear = 0;
    char sequenceString[2*MAX_LOADSTRING] = { 0 };
    double i37 = 0.0;
    size_t batchLoc1 = 0;
    size_t batchLoc2 = 0;

    //fitting
    bool isInFittingMode;
    char fittingString[1024] = { 0 };
    char fittingPath[1024] = { 0 };
    double fittingArray[256] = { 0 };
    double* fittingReference = 0;
    int Nfitting = 0;
    int fittingMode = 0;
    int fittingMaxIterations = 0;
    size_t fittingROIstart = 0;
    size_t fittingROIstop = 0;
    size_t fittingROIsize = 0;
    double fittingResult[64];
    double fittingError[64];
};

int             loadSavedFields(simulationParameterSet* sCPU, const char* outputBase);
int             removeCharacterFromString(char* cString, size_t N, char removedChar);
int				fftshiftZ(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int             fftshiftD2Z(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int				fftshiftAndFilp(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int             loadReferenceSpectrum(char* spectrumPath, simulationParameterSet* sCPU);
int             readFittingString(simulationParameterSet* sCPU);
int             saveSettingsFile(simulationParameterSet* sCPU);
void            unixNewLine(FILE* iostream);
double          saveSlurmScript(simulationParameterSet* sCPU, int gpuType, int gpuCount, size_t totalSteps);
int				loadFrogSpeck(char* frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel);
double          cModulusSquared(std::complex<double>complexNumber);
int             allocateGrids(simulationParameterSet* sCPU);
int             deallocateGrids(simulationParameterSet* sCPU, bool alsoDeleteDisplayItems);
int             configureBatchMode(simulationParameterSet* sCPU);
int             saveDataSet(simulationParameterSet* sCPU);
int             readInputParametersFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, const char* filePath);
int             loadPulseFiles(simulationParameterSet* sCPU);
int             skipFileUntilCharacter(FILE* fstream, char target);
int             copyParamsIntoStrings(char parameterBlock[22][256], const char* cc, int n);
void            applyOp(char op, double* result, double* readout);
double          parameterStringToDouble(const char* pString, double* iBlock, double* vBlock);
std::string     getBasename(char* fullPath);
void            stripWhiteSpace(char* sequenceString, size_t bufferSize);
void            stripLineBreaks(char* sequenceString, size_t bufferSize);
int             interpretParameters(std::string cc, int n, double *iBlock, double *vBlock, double *parameters, bool* defaultMask);
