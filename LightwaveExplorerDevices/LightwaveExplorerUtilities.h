#pragma once
#include <complex>
#include <cstring>
#ifdef __CUDACC__
#include <cufft.h>
#include <thrust/complex.h>
#include <fftw3_mkl.h>
#define deviceLib thrust
typedef thrust::complex<double> deviceComplex;
typedef double deviceFP;
#elif defined RUNONSYCL
#include <oneapi/dpl/complex>
#include <oneapi/dpl/cmath>
typedef oneapi::dpl::complex<double> deviceComplex;
typedef double deviceFP;
#define deviceLib oneapi::dpl
#include <sycl/sycl.hpp>
#elif defined CPUONLY
#include <fftw3.h>
typedef std::complex<double> deviceComplex;
typedef double deviceFP;
#define deviceLib std
#else
#include <fftw3_mkl.h>
typedef std::complex<double> deviceComplex;
typedef double deviceFP;
#define deviceLib std
#endif

#define THREADS_PER_BLOCK 32
#define MIN_GRIDDIM 8
#define ANGLETOLERANCE 1e-12
#define FALSE 0
#define TRUE 1
#define MAX_LOADSTRING 1024
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
#define maxN(a,b)            (((a) > (b)) ? (a) : (b))
#define minN(a,b)            (((a) < (b)) ? (a) : (b))

typedef struct crystalEntry {
    char crystalNameW[256] = { 0 };
    int axisType = 0;
    int sellmeierType = 0;
    int nonlinearSwitches[4] = { 0 };
    double sellmeierCoefficients[66] = { 0 };
    char sellmeierReference[512] = { 0 };
    double d[18] = { 0 };
    char dReference[512] = { 0 };
    double chi3[81] = { 0 };
    char chi3Reference[512] = { 0 };
    double absorptionParameters[6] = { 0 };
    char spectralFile[512] = { 0 };
    double spectralData[2048] = { 0 };
    double nonlinearReferenceFrequencies[7] = { 0 };
    int numberOfEntries = 0;
} crystalEntry;

typedef struct pulse {
    double energy;
    double frequency;
    double bandwidth;
    int sgOrder;
    double cep;
    double delay;
    double gdd;
    double tod;
    int phaseMaterial;
    double phaseMaterialThickness;
    double beamwaist;
    double x0;
    double y0;
    double z0;
    double beamAngle;
    double polarizationAngle;
    double beamAnglePhi;
    double circularity;
    double pulseSum;
} pulse;

//Simulation parameter struct to pass to the simulations running in threads
typedef struct simulationParameterSet {
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
    pulse pulse1;
    pulse pulse2;
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
    double neref = 0;
    double noref = 0;
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
    double* InfoVec = 0;
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
    double fittingPrecision = 0;
    double* fittingReference = 0;
    int Nfitting = 0;
    int fittingMode = 0;
    int fittingMaxIterations = 0;
    size_t fittingROIstart = 0;
    size_t fittingROIstop = 0;
    size_t fittingROIsize = 0;
    double fittingResult[64];
    double fittingError[64];

} simulationParameterSet;

typedef struct deviceParameterSet {
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
    deviceFP* J0 = 0;
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
} deviceParameterSet;

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
int             readCrystalDatabase(crystalEntry* db);
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
