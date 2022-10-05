#pragma once
#include <complex>
#ifdef __CUDACC__
#include <cufft.h>
#include <thrust/complex.h>
#define deviceLib thrust
#define deviceComplex thrust::complex<double>
#elif defined RUNONSYCL
#include <fftw3_mkl.h>
#define deviceComplex
#define deviceLib std
#define deviceComplex std::complex<double>
#else
#include <fftw3_mkl.h>
#define deviceComplex std::complex<double>
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
#define DEG2RAD 1.7453292519943295e-02
#define LIGHTC 2.99792458e8
#define EPS0 8.8541878128e-12
#define SIXTH 0.1666666666666667
#define THIRD 0.3333333333333333
#define KLORENTZIAN 3182.607353999257 //(e * e / (epsilon_o * m_e)
#define maxN(a,b)            (((a) > (b)) ? (a) : (b))
#define minN(a,b)            (((a) < (b)) ? (a) : (b))

typedef struct crystalEntry {
    wchar_t crystalNameW[256] = { 0 };
    int axisType = 0;
    int sellmeierType = 0;
    int nonlinearSwitches[4] = { 0 };
    double sellmeierCoefficients[66] = { 0 };
    wchar_t sellmeierReference[512] = { 0 };
    double d[18] = { 0 };
    wchar_t dReference[512] = { 0 };
    double chi3[81] = { 0 };
    wchar_t chi3Reference[512] = { 0 };
    double absorptionParameters[6] = { 0 };
    wchar_t spectralFile[512] = { 0 };
    double spectralData[2048] = { 0 };
    double nonlinearReferenceFrequencies[7] = { 0 };
    int numberOfEntries = 0;
} crystalEntry;

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
    double pulse1measEnergy = 0;

    //spatial beam properties = 0;
    double beamwaist1 = 0;
    double beamwaist2 = 0;
    double z01 = 0;
    double z02 = 0;
    double x01 = 0;
    double x02 = 0;
    double y01 = 0;
    double y02 = 0;
    double propagationAngle1 = 0;
    double propagationAngle2 = 0;
    double propagationAnglePhi1 = 0;
    double propagationAnglePhi2 = 0;
    bool isCylindric = 0;
    bool is3D = 0;
    int symmetryType = 0;

    //spectral/temporal field properties
    double pulseEnergy1 = 0;
    double pulseEnergy2 = 0;
    double frequency1 = 0;
    double frequency2 = 0;
    int sgOrder1 = 0;
    int sgOrder2 = 0;
    double bandwidth1 = 0;
    double bandwidth2 = 0;
    double cephase1 = 0;
    double cephase2 = 0;
    double delay1 = 0;
    double delay2 = 0;
    double gdd1 = 0;
    double gdd2 = 0;
    double tod1 = 0;
    double tod2 = 0;
    double phaseMaterialThickness1 = 0;
    double phaseMaterialThickness2 = 0;
    int phaseMaterialIndex1 = 0;
    int phaseMaterialIndex2 = 0;

    //loaded FROG/EOS fields
    std::complex<double>* loadedField1 = 0;
    std::complex<double>* loadedField2 = 0;
    bool field1IsAllocated = 0;
    bool field2IsAllocated = 0;
    int pulse1FileType = 0;
    int pulse2FileType = 0;
    char field1FilePath[MAX_LOADSTRING] = { 0 };
    char field2FilePath[MAX_LOADSTRING] = { 0 };


    //polarization properties
    double polarizationAngle1 = 0;
    double polarizationAngle2 = 0;
    double circularity1 = 0;
    double circularity2 = 0;

    int pulsetype = 0;
    double* InfoVec = 0;
    double* ExtOut = 0;
    std::complex<double>* EkwOut = 0;
    double* totalSpectrum = 0;
    int* imdone = 0;
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
    char sequenceString[MAX_LOADSTRING] = { 0 };
    double sequenceArray[MAX_LOADSTRING] = { 0 };
    int Nsequence = 0;

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

typedef struct cudaParameterSet {
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
    double* inverseChiLinear1 = 0;
    double* inverseChiLinear2 = 0;
    double* fieldFactor1 = 0;
    double* fieldFactor2 = 0;
    deviceComplex* k1 = 0;
    deviceComplex* k2 = 0;

    double* gridRadialLaplacian1 = 0;
    double* gridRadialLaplacian2 = 0;
    double* gridETime1 = 0;
    double* gridETime2 = 0;
    double* gridPolarizationTime1 = 0;
    double* gridPolarizationTime2 = 0;
    double* expGammaT = 0;
    double* gridPlasmaCurrent1 = 0;
    double* gridPlasmaCurrent2 = 0;

    //fixed length arrays
    double firstDerivativeOperation[6] = { 0 };
    double plasmaParameters[6] = { 0 }; //[dt^2 * e^2/m * nonlinearAbsorptionStrength, gamma] 
    double chi2Tensor[18] = { 0 };
    double chi3Tensor[81] = { 0 };
    double absorptionParameters[6] = { 0 };
    double rotationForward[9] = { 0 };
    double rotationBackward[9] = { 0 };
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
    double fftNorm = 0;
    int axesNumber = 0;
    int sellmeierType = 0;
    double f0 = 0;
    double fStep = 0;
    double dt = 0;
    double dx = 0;
    double dk1 = 0;
    double dk2 = 0;
    double h = 0;
    size_t Nsteps = 0;
    int Nthread = 0;
    int NblockC = 0;
    int Nblock = 0;
} cudaParameterSet;

int             loadSavedFields(simulationParameterSet* sCPU, char* outputBase);
int             removeCharacterFromString(char* cString, size_t N, char removedChar);
int				fftshiftZ(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int             fftshiftD2Z(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int				fftshiftAndFilp(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int             loadReferenceSpectrum(char* spectrumPath, simulationParameterSet* sCPU);
int             readFittingString(simulationParameterSet* sCPU);
int             saveSettingsFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr);
void            unixNewLine(FILE* iostream);
int             saveSlurmScript(simulationParameterSet* sCPU, int gpuType, int gpuCount);
int				loadFrogSpeck(char* frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel);
double          cModulusSquared(std::complex<double>complexNumber);
int             allocateGrids(simulationParameterSet* sCPU);
int             deallocateGrids(simulationParameterSet* sCPU, bool alsoDeleteDisplayItems);
int             readCrystalDatabase(crystalEntry* db);
int             readSequenceString(simulationParameterSet* sCPU);
int             configureBatchMode(simulationParameterSet* sCPU);
int             saveDataSet(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, char* outputbase, bool saveInputs);
int             readInputParametersFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, char* filePath);
int             loadPulseFiles(simulationParameterSet* sCPU);
int             skipFileUntilCharacter(FILE* fstream, char target);

