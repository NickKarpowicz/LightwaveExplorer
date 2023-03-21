#pragma once
#include <complex>
#include <vector>
#include <array>
#include <fstream>
#include <atomic>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

static const unsigned int threadsPerBlock = 32;
static const unsigned int minGridDimension = 8;
static const size_t pathArrayLength = 1024;

#ifndef LWEFLOATINGPOINT
#define LWEFLOATINGPOINT 64
#endif
#if LWEFLOATINGPOINT==32
#define LWEFLOATINGPOINTTYPE float
#else
#define LWEFLOATINGPOINTTYPE double
#endif

//ugly macro prefix because CUDA doesn't like my nice constexprs
#ifdef __CUDACC__
#define hostOrDevice __host__ __device__
#else
#define hostOrDevice
#endif

//variadic template to constexpr the product of a bunch of values
//in a way that keeps Xe Graphics happy (no doubles)
//convention: if there are multiple types as inputs, 
//return type is the type of the first argument
template<typename T>
hostOrDevice static constexpr T constProd(T x) {
    return x;
}
template<typename T, typename... Args>
hostOrDevice static constexpr T constProd(T x, Args... args) {
    return (T)(x * constProd(args...));
}

template <typename T, typename U>
hostOrDevice static constexpr inline T maxN(T a, U b) {
    return (((a) > (b)) ? (a) : (b));
}

template <typename T, typename U>
hostOrDevice static constexpr inline T minN(T a, U b) {
    return (((a) < (b)) ? (a) : (b));
}
//pi to stupid digits
template <typename T>
hostOrDevice static constexpr T vPi() {
    return (T)3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384L;
}
template <typename T>
hostOrDevice static constexpr T twoPi() {
    return (T)(2.0 * vPi<T>());
}
template <typename T>
hostOrDevice static constexpr T angleTolerance() {
    return (T)1e-12;
}
template <typename T>
hostOrDevice static constexpr T invSqrtPi() {
    return (T)(1.0/sqrt(vPi<T>()));
}
template <typename T>
hostOrDevice static constexpr T deg2Rad() {
    return (T)(vPi<T>()/180.0);
}
template <typename T>
hostOrDevice static constexpr T rad2Deg() {
    return (T)(180.0/vPi<T>());
}
template <typename T>
hostOrDevice static constexpr T lightC() {
    return (T)2.99792458e8;
}
template <typename T>
hostOrDevice static constexpr T eps0() {
    return (T)8.8541878128e-12;
}
template <typename T>
hostOrDevice static constexpr T sixth() {
    return (T)((T)1.0/6.0);
}
template <typename T>
hostOrDevice static constexpr T third() {
    return (T)((T)1.0/3.0);
}
template <typename T>
hostOrDevice static constexpr T kLorentzian() {
    return (T)3182.607353999257;
}

#include <complex>
#include <cstring>
#ifdef __CUDACC__
#include <cufft.h>
#include <thrust/complex.h>
#include <fftw3_mkl.h>
#elif defined RUNONSYCL
#include <oneapi/dpl/complex>
#include <oneapi/dpl/cmath>
#include <sycl/sycl.hpp>
#elif defined CPUONLY
#include <fftw3.h>
#else
#include <fftw3_mkl.h>
#endif

//Enum for determining the FFT type:
// D2Z: real to complex (time to frequency)
// Z2D: complex to real (f to t)
// D2Z_1D: real-to-complex, on the time/frequency axis only (to a grid in space vs. frequency)
// Z2D_1D: complex to real, time/frequency axis only (to space vs. time)
// D2Z_Polarization: in the case of cylindrical symmetry, a double-sized fft in the de-interlaced
//                   representation. If plasma calculations are also being done, their FFT will
//                   be included here as a batch. Thus, this is the most expensive operation
//                   in such a propagation.
enum class deviceFFT : int {
    D2Z = 0,
    Z2D = 1,
    D2Z_1D = 2,
    Z2D_1D = 3,
    D2Z_Polarization = 4
};

//Determine the type of data transfer - not necessary on all devices, but should be specified
//consistently in the abstracted functions.
// ToDevice: source is the host, destination is device
// ToHost: source is the device, destination is host
// OnDevice: device is both source and destination
// there is no on-host memory transfer to discourage memcpy-like operations where not required...
enum class copyType : int {
    ToDevice = 1,
    ToHost =  2,
    OnDevice = 3
};

//class holding the device data structures
//note that it uses c-style arrays-- this is for compatibility
//with all of the platforms involved, and because it is transferred
//to the device with a memcpy-like operation, so constructors
//would not be called.
template <typename deviceFP, typename deviceComplex>
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

//Class which holds a single entry in the crystal database
class crystalEntry {
public:
    std::string crystalName;
    int axisType = 0;
    int sellmeierType = 0;
    std::array<int,4> nonlinearSwitches = {};
    std::array<double,66> sellmeierCoefficients = {};
    std::string sellmeierReference;
    std::array<double,18> d = {};
    std::string dReference;
    std::array<double,81> chi3 = {};
    std::string chi3Reference;
    std::array<double,6> absorptionParameters = {};
    std::string spectralFile;
    std::array<double,7> nonlinearReferenceFrequencies = {};
};

//Crystal database class; primarily holds a std::vector of crystalEntry elements
//comprising the database, plus method for loading the database from the file
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
            newEntry.spectralFile = line;
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

//templated class for describing a pulse in various floating point representations
//copyable between representations (required for strict FP32 mode)
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
    size_t Npropagation = 0;
    size_t Ntime = 0;
    size_t Nfreq = 0;
    size_t Nspace = 0;
    size_t Nspace2 = 0;
    size_t Ngrid = 0;
    size_t NgridC = 0;
    size_t Nsims = 0;
    size_t Nsims2 = 0;
    std::atomic_uint32_t* progressCounter = 0;
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
    std::string field1FilePath;
    std::string field2FilePath;

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
    std::string outputBasePath;
    int runType = 0;
    bool runningOnCPU = 0;

    //sequence
    bool isInSequence = 0;
    bool isFollowerInSequence = 0;
    bool isReinjecting = 0;
    bool forceLinear = 0;
    std::string sequenceString;
    double i37 = 0.0;
    size_t batchLoc1 = 0;
    size_t batchLoc2 = 0;

    //fitting
    bool isInFittingMode = false;
    std::string fittingString;
    std::string fittingPath;
    std::array<double, 256> fittingArray = {};
    double* fittingReference = 0;
    int Nfitting = 0;
    int fittingMode = 0;
    int fittingMaxIterations = 0;
    size_t fittingROIstart = 0;
    size_t fittingROIstop = 0;
    size_t fittingROIsize = 0;
    std::array<double, 64> fittingResult = {};
    std::array<double, 64> fittingError = {};
};

int             loadSavedFields(simulationParameterSet* sCPU, const char* outputBase);
int             removeCharacterFromString(char* cString, size_t N, char removedChar);
void            removeCharacterFromString(std::string& s, char removedChar);
int				fftshiftZ(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int             fftshiftD2Z(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int				fftshiftAndFilp(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int             loadReferenceSpectrum(std::string spectrumPath, simulationParameterSet* sCPU);
int             readFittingString(simulationParameterSet* sCPU);
int             saveSettingsFile(const simulationParameterSet* sCPU);
double          saveSlurmScript(simulationParameterSet* sCPU, int gpuType, int gpuCount, size_t totalSteps);
int				loadFrogSpeck(std::string frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel);
double          cModulusSquared(const std::complex<double>& x);
int             allocateGrids(simulationParameterSet* sCPU);
int             deallocateGrids(simulationParameterSet* sCPU, bool alsoDeleteDisplayItems);
int             configureBatchMode(simulationParameterSet* sCPU);
int             saveDataSet(simulationParameterSet* sCPU);
int             readInputParametersFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, const char* filePath);
int             loadPulseFiles(simulationParameterSet* sCPU);
int             skipFileUntilCharacter(FILE* fstream, char target);
void            applyOp(char op, double* result, double* readout);
double          parameterStringToDouble(std::string& ss, double* iBlock, double* vBlock);
std::string     getBasename(char* fullPath);
std::string     getBasename(const std::string& fullPath);
void            stripWhiteSpace(char* sequenceString, size_t bufferSize);
void            stripWhiteSpace(std::string& s);
void            stripLineBreaks(std::string& s);
int             interpretParameters(std::string cc, int n, double *iBlock, double *vBlock, double *parameters, bool* defaultMask);
