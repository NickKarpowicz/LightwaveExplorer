#pragma once
#include <complex>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <atomic>
#include <mutex>
#include <algorithm>
#include "../LightwaveExplorerDevices/LightwaveExplorerHelpers.h"
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#ifdef __linux__
#include <unistd.h>
#endif

static const unsigned int threadsPerBlock = 64;
static const unsigned int minGridDimension = 8;


std::string     getBasename(const std::string& fullPath);
int				loadFrogSpeck(const std::string& frogFilePath, std::complex<double>* Egrid, const int64_t Ntime, const double fStep, const double gateLevel);
int             loadSavedGridFile(const std::string& filePath, std::vector<double>& outputGrid, int64_t Ngrid);
int             loadSavedGridFileMultiple(const std::string& filePath, std::vector<double>& outputGrid, int64_t Ngrid, int64_t Nsims);
int             loadWaveformFile(const std::string& filePath, std::complex<double>* outputGrid, const int64_t Ntime, const double fStep);
double          cModulusSquared(const std::complex<double>& x);
void            applyOp(const char op, double* result, const double* readout);
double          parameterStringToDouble(const std::string& ss, const double* iBlock, const double* vBlock);
void            stripWhiteSpace(std::string& s);
void            stripLineBreaks(std::string& s);
int             interpretParameters(const std::string& cc, const int n, const double *iBlock, const double *vBlock, double *parameters, bool* defaultMask);


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

template <typename deviceFP>
class maxwellPoint {
public:
    deviceFP x{};
    deviceFP y{};
    deviceFP z{};
    hostOrDevice inline deviceFP& operator()(const int index) {
        switch (index) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: return x;
        }
    }
    hostOrDevice inline void operator+=(
        const maxwellPoint<deviceFP>& other) {
        x += other.x;
        y += other.y;
        z += other.z;
    }
    hostOrDevice inline void operator-=(
        const maxwellPoint<deviceFP>& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
    }
    hostOrDevice inline void operator*=(
        const maxwellPoint<deviceFP>& other) {
        x *= other.x;
        y *= other.y;
        z *= other.z;
    }
    hostOrDevice inline void operator/=(
        const maxwellPoint<deviceFP>& other) {
        x /= other.x;
        y /= other.y;
        z /= other.z;
    }

    hostOrDevice inline void operator+=(
        const deviceFP other) {
        x += other;
        y += other;
        z += other;
    }
    hostOrDevice inline void operator-=(
        const deviceFP other) {
        x -= other;
        y -= other;
        z -= other;
    }
    hostOrDevice inline void operator*=(
        const deviceFP other) {
        x *= other;
        y *= other;
        z *= other;
    }
    hostOrDevice inline void operator/=(
        const deviceFP other) {
        x /= other;
        y /= other;
        z /= other;
    }

    hostOrDevice inline maxwellPoint<deviceFP> operator*(
        const maxwellPoint<deviceFP>& other) const {
        return maxwellPoint<deviceFP>{
                x * other.x,
                y * other.y,
                z * other.z};
    }
    hostOrDevice inline maxwellPoint<deviceFP> operator*(
        const deviceFP other) const {
        return maxwellPoint<deviceFP>{
                x * other,
                y * other,
                z * other};
    }
    hostOrDevice inline friend maxwellPoint<deviceFP> operator*(const deviceFP a, const maxwellPoint<deviceFP>& b) {
        return maxwellPoint<deviceFP>{
                a * b.x,
                a * b.y,
                a * b.z};
    }

    hostOrDevice inline maxwellPoint<deviceFP> operator/(
        const maxwellPoint<deviceFP>& other) const {
        return maxwellPoint<deviceFP>{
                x / other.x,
                y / other.y,
                z / other.z};
    }
    hostOrDevice inline maxwellPoint<deviceFP> operator/(
        const deviceFP other) const {
        return maxwellPoint<deviceFP>{
                x / other,
                y / other,
                z / other};
    }
    hostOrDevice inline friend maxwellPoint<deviceFP> operator/(const deviceFP a, const maxwellPoint<deviceFP>& b) {
        return maxwellPoint<deviceFP>{
                a / b.x,
                a / b.y,
                a / b.z};
    }

    hostOrDevice inline maxwellPoint<deviceFP> operator+(
        const maxwellPoint<deviceFP>& other) const {
        return maxwellPoint<deviceFP>{
                x + other.x,
                y + other.y,
                z + other.z};
    }
    hostOrDevice inline maxwellPoint<deviceFP> operator+(
        const deviceFP other) const {
        return maxwellPoint<deviceFP>{
                x + other,
                y + other,
                z + other};
    }
    hostOrDevice inline friend maxwellPoint<deviceFP> operator+(const deviceFP a, const maxwellPoint<deviceFP>& b) {
        return maxwellPoint<deviceFP>{
                a + b.x,
                a + b.y,
                a + b.z};
    }

    hostOrDevice inline maxwellPoint<deviceFP> operator-(
        const maxwellPoint<deviceFP>& other) const {
        return maxwellPoint<deviceFP>{
                x - other.x,
                y - other.y,
                z - other.z};
    }
    hostOrDevice inline maxwellPoint<deviceFP> operator-(
        const deviceFP other) const {
        return maxwellPoint<deviceFP>{
                x - other,
                y - other,
                z - other};
    }
    hostOrDevice inline friend maxwellPoint<deviceFP> operator-(const deviceFP a, const maxwellPoint<deviceFP>& b) {
        return maxwellPoint<deviceFP>{
                a - b.x,
                a - b.y,
                a - b.z};
    }
};

template <typename deviceFP>
class maxwellKPoint {
public:
    maxwellPoint<deviceFP> kE;
    maxwellPoint<deviceFP> kH;
};

template <typename deviceFP>
class oscillator {
public:
    maxwellPoint<deviceFP> J;
    maxwellPoint<deviceFP> P;

    //note that I only defined the three operations I need rather than the
    //full set. Might be worth filling in everything later.
    hostOrDevice void operator+=(
        const oscillator<deviceFP>& other) {
        J += other.J;
        P += other.P;
    }
    hostOrDevice inline oscillator<deviceFP> operator*(
        const deviceFP other) const {
        return oscillator<deviceFP>{
                J * other,
                P * other,
        };
    }
    hostOrDevice inline oscillator<deviceFP> operator+(
        const oscillator<deviceFP>& other) const {
        return oscillator<deviceFP>{
                J + other.J,
                P + other.P
        };
    }
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
    int64_t Ntime = 0;
    int64_t Nfreq = 0;
    int64_t Nspace = 0;
    int64_t Nspace2 = 0;
    int64_t Ngrid = 0;
    int64_t NgridC = 0;
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
    int64_t Nsteps = 0;
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
        #include <mach-o/dyld.h>
        uint32_t bufferSize = 1024;
        char sysPath[1024] = { 0 };
        _NSGetExecutablePath(sysPath, &bufferSize);
        std::string macPath(sysPath);
        int64_t posPath = macPath.find_last_of("/");
        std::string databasePath = macPath.substr(0, posPath).append("/../Resources/CrystalDatabase.txt");
        std::ifstream fs(databasePath);
        if (!fs.is_open()) {
            fs.open("CrystalDatabase.txt");
        }
#elif defined __linux__
        char pBuf[256];
        int64_t len = sizeof(pBuf); 
        int bytes = minN(readlink("/proc/self/exe", pBuf, len), len - 1);
        if(bytes >= 0)
            pBuf[bytes] = '\0';
        std::string binPath(pBuf);
        int64_t posPath = binPath.find_last_of("/");
        std::string databasePath = binPath.substr(0, posPath).append("/../share/LightwaveExplorer/CrystalDatabase.txt");
        std::ifstream fs(databasePath);
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
    double drudeGamma = 0;
    double crystalTheta = 0;
    double crystalPhi = 0;
    double crystalThickness = 0;
    double* chi2Tensor = 0;
    double* chi3Tensor = 0;
    double* sellmeierCoefficients = 0;
    double* absorptionParameters = 0;
    int sellmeierType = 0;
    int axesNumber = 0;
    int* nonlinearSwitches = 0;
    bool isCylindric = 0;
    bool is3D = 0;
    bool isFDTD = 0;
    int symmetryType = 0;

    //loaded FROG/EOS fields
    std::complex<double>* loadedField1 = 0;
    std::complex<double>* loadedField2 = 0;
    double* loadedFullGrid1 = 0;
    double* loadedFullGrid2 = 0;
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
    int64_t batchLoc1 = 0;
    int64_t batchLoc2 = 0;

    //fitting
    bool isInFittingMode = false;
    std::string fittingString;
    std::string fittingPath;
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

    [[nodiscard]] constexpr double getByNumberWithMultiplier(const size_t index) {
        if (index == 0 || index == 36 || index >= multipliers.size()) return 0.0;
        return  getByNumber(index) / multipliers[index];
    }
    
    constexpr double getByNumber(const size_t index) {
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
    void setByNumberWithMultiplier(const size_t index, const double value);
    int loadSavedFields(const std::string& outputBase);
    int loadReferenceSpectrum();
    int readInputParametersFile(crystalEntry* crystalDatabasePtr, const std::string filePath);
    int saveSettingsFile();
    double saveSlurmScript(int gpuType, int gpuCount, int64_t totalSteps);
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
        (*s).fftNorm = (deviceFP)1.0 / (*s).Ngrid;
        (*s).dt = (deviceFP)tStep;
        (*s).dx = (deviceFP)rStep;
        (*s).dk1 = (deviceFP)(twoPi<double>() / (Nspace * rStep));
        (*s).dk2 = (deviceFP)(twoPi<double>() / (Nspace2 * rStep));
        (*s).fStep = (deviceFP)fStep;
        (*s).Nsteps = (int64_t)round(crystalThickness / propagationStep);
        (*s).h = (deviceFP)crystalThickness / ((*s).Nsteps); //adjust step size so that thickness can be varied continuously by fitting
        (*s).axesNumber = axesNumber;
        (*s).sellmeierType = sellmeierType;
        (*s).crystalPhi = (deviceFP)crystalPhi;
        (*s).crystalTheta = (deviceFP)crystalTheta;
        (*s).f0 = (deviceFP)pulse1.frequency;
        (*s).Nthread = threadsPerBlock;
        (*s).Nblock = (int)((*s).Ngrid / threadsPerBlock);
        (*s).NblockC = (int)((*s).NgridC / threadsPerBlock);
        (*s).isCylindric = isCylindric;
        (*s).forceLinear = forceLinear;
        (*s).isNonLinear = (nonlinearSwitches[0] + nonlinearSwitches[1]) > 0;
        (*s).isUsingMillersRule = (crystalDatabase[materialIndex].nonlinearReferenceFrequencies[0]) != 0.0;

        if (nonlinearAbsorptionStrength > 0.) {
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

        double firstDerivativeOperation[6] = { -1. / 60.,  3. / 20., -3. / 4.,  3. / 4.,  -3. / 20., 1. / 60. };
        for (int64_t i = 0; i < 6; ++i) {
            firstDerivativeOperation[i] *= (-2.0 / ((*s).dx));
        }

        //set nonlinearSwitches[3] to the number of photons needed to overcome bandgap
        nonlinearSwitches[3] = (int)ceil(bandGapElectronVolts * 241.79893e12 / pulse1.frequency) - 2;
        double plasmaParametersCPU[6] = { 0 };


        plasmaParametersCPU[0] = nonlinearAbsorptionStrength; //nonlinear absorption strength parameter
        plasmaParametersCPU[1] = drudeGamma; //gamma
        if (nonlinearAbsorptionStrength > 0.) {
            plasmaParametersCPU[2] = tStep * tStep
                * 2.817832e-08 / (1.6022e-19 * bandGapElectronVolts * effectiveMass); // (dt^2)*e* e / (m * band gap));
        }
        else {
            plasmaParametersCPU[2] = 0;
        }

        for (int j = 0; j < 18; ++j) {
            (*s).chi2Tensor[j] = (deviceFP)(2e-12 * chi2Tensor[j]); //go from d in pm/V to chi2 in m/V
            if (j > 8) (*s).chi2Tensor[j] *= 2.0; //multiply cross-terms by 2 for consistency with convention
        }

        for (int i = 0; i < 4; ++i) {
            (*s).nonlinearSwitches[i] = nonlinearSwitches[i];
        }

        for (int64_t i = 0; i < 81; i++) {
            (*s).chi3Tensor[i] = (deviceFP)chi3Tensor[i];
        }

        for (int64_t i = 0; i < 6; i++) {
            (*s).absorptionParameters[i] = (deviceFP)absorptionParameters[i];
        }

        for (int64_t i = 0; i < 6; i++) {
            (*s).plasmaParameters[i] = (deviceFP)plasmaParametersCPU[i];
        }

        for (int64_t i = 0; i < 6; i++) {
            (*s).firstDerivativeOperation[i] = (deviceFP)firstDerivativeOperation[i];
        }
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
    std::vector<simulationParameterSet> parameters;
public:
    std::vector<std::mutex> mutexes = std::vector<std::mutex>(1);

    simulationBatch() {
        parameters = std::vector<simulationParameterSet>(1);
    }
    ~simulationBatch() {
        std::for_each(mutexes.begin(), mutexes.end(), 
            [](std::mutex& m) {std::lock_guard<std::mutex> lock(m); });
    }
    void configure();
    void loadPulseFiles();
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

template <typename deviceFP, typename E, typename H, typename O>
class maxwellCalculation {
public:
    E* Egrid{};
    E* EgridNext{};
    E* EgridEstimate{};
    E* EgridEstimate2{};
    H* Hgrid{};
    H* HgridNext{};
    H* HgridEstimate{};
    H* HgridEstimate2{};
    O* materialGrid{};
    O* materialGridNext{};
    O* materialGridEstimate{};
    O* materialGridEstimate2{};
    int64_t* materialIndexMap{};
    maxwellPoint<deviceFP> sellmeierEquations[22][8]{};
    maxwellPoint<deviceFP> chi3[27][8]{};
    maxwellPoint<deviceFP> chi2[6][8]{};
    int nonlinearAbsorptionOrder[8]{};
    deviceFP kNonlinearAbsorption[8]{};
    deviceFP kDrude[8]{};
    deviceFP gammaDrude[8]{};
    deviceFP kCarrierGeneration[8]{};
    deviceFP rotateForward[9][8]{};
    deviceFP rotateBackward[9][8]{};
    bool hasChi2[8]{};
    bool hasFullChi3[8]{};
    bool hasSingleChi3[8]{};
    bool hasPlasma[8]{};
    deviceFP* inOutEy{};
    deviceFP* inOutEx{};
    deviceFP* inputExFFT{};
    deviceFP* inputEyFFT{};
    deviceFP omegaStep{};
    deviceFP xyStep{};
    deviceFP zStep{};
    deviceFP tStep{};
    deviceFP frontBuffer{};
    deviceFP backBuffer{};
    deviceFP crystalThickness{};
    deviceFP inverseXyStep{};
    deviceFP inverseZStep{};
    deviceFP omegaMax{};
    int64_t observationPoint{};
    int64_t waitFrames{};
    int64_t Nx{};
    int64_t Ny{};
    int64_t Nz{};
    int64_t Nt{};
    int64_t Ngrid{};
    int64_t NMaterialGrid{};
    int Noscillators{};
    int64_t NtIO{};
    int64_t frequencyLimit{};
    int64_t tGridFactor=1;
    int64_t materialStart{};
    int64_t materialStop{};
    maxwellCalculation<deviceFP, E, H, O>* deviceCopy = nullptr;
    
    void fillRotationMatricies(double crystalTheta, double crystalPhi, int64_t crystalNumber) {
        double cosT = cos(crystalTheta);
        double sinT = sin(crystalTheta);
        double cosP = cos(crystalPhi);
        double sinP = sin(crystalPhi);
        double forward[9] =
        { cosT * cosP, sinP, -sinT * cosP, 
            -sinP * cosT, cosP, sinP * sinT, 
            sinT, 0.0, cosT };

        //reverse direction (different order of operations)
        double backward[9] =
        { cosT * cosP, -sinP * cosT, sinT, 
            sinP, cosP, 0.0, 
            -sinT * cosP, sinP * sinT, cosT };

        for (int64_t i = 0; i < 9; i++) {
            rotateForward[i][crystalNumber] = static_cast<deviceFP>(forward[i]);
            rotateBackward[i][crystalNumber] = static_cast<deviceFP>(backward[i]);
        }
    }
    maxwellCalculation(simulationParameterSet* s, int64_t timeFactor, deviceFP zStep_in, deviceFP frontBuffer_in, deviceFP backBuffer_in) {
        frontBuffer = frontBuffer_in;
        backBuffer = backBuffer_in;
        crystalThickness = (*s).crystalThickness;
        zStep = zStep_in;
        Nx = (*s).Nspace;
        Ny = ((*s).is3D) ? (*s).Nspace2 : 1;
        Nz = (frontBuffer + backBuffer + crystalThickness) / zStep;
        Nz = minGridDimension * (Nz / minGridDimension + (Nz % minGridDimension > 0));
        NtIO = (*s).Ntime;
        xyStep = (*s).rStep;
        tStep = (*s).tStep / timeFactor;
        omegaMax = 0.1 * twoPi<deviceFP>()*lightC<deviceFP>() / zStep;
        omegaStep = twoPi<deviceFP>() * s->fStep;
        frequencyLimit = minN(static_cast<int64_t>(omegaMax / omegaStep),s->Nfreq);
        inverseXyStep = 1.0 / xyStep;
        inverseZStep = 1.0 / zStep;
        materialStart = frontBuffer / zStep;
        materialStop = materialStart + (crystalThickness / zStep);
        observationPoint = materialStop + 10;
        tGridFactor = timeFactor;
        Ngrid = Nz * Ny * Nx;
        fillRotationMatricies((*s).crystalTheta, (*s).crystalPhi, 0);
    }

    
};
