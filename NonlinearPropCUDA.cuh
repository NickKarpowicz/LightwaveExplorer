#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <complex>
#include <cuComplex.h>
#include "cufft.h"

#define MAX_LOADSTRING 1024

typedef struct crystalEntry {
    wchar_t crystalNameW[256];
    int axisType;
    int sellmeierType;
    int nonlinearSwitches[4];
    double sellmeierCoefficients[66];
    wchar_t sellmeierReference[512];
    double d[18];
    wchar_t dReference[512];
    double chi3[81];
    wchar_t chi3Reference[512];
    double absorptionParameters[6];
    wchar_t spectralFile[512];
    double spectralData[2048];
    double nonlinearReferenceFrequencies[7];
    int numberOfEntries;
} crystalEntry;

//Simulation parameter struct to pass to the simulations running in threads
typedef struct simulationParameterSet {
    double rStep;
    double tStep;
    double fStep;
    double kStep;
    double propagationStep;
    size_t Npropagation;
    size_t Ntime;
    size_t Nfreq;
    size_t Nspace;
    size_t Nspace2;
    size_t Ngrid;
    size_t NgridC;
    size_t Nsims;
    size_t Nsims2;
    size_t *progressCounter;
    double spatialWidth;
    double spatialHeight;
    double timeSpan;
    int materialIndex;
    int materialIndexAlternate;
    double bandGapElectronVolts;
    double effectiveMass;
    double nonlinearAbsorptionStrength;
    double drudeGamma;
    double crystalTheta;
    double crystalPhi;
    double crystalThickness;
    double* chi2Tensor;
    double* deffTensor;
    double* chi3Tensor;
    double* sellmeierCoefficients;
    double* absorptionParameters;
    int sellmeierType;
    int axesNumber;
    double neref;
    double noref;
    int* nonlinearSwitches;
    double pulse1measEnergy;

    //spatial beam properties;
    double beamwaist1;
    double beamwaist2;
    double z01;
    double z02;
    double x01;
    double x02;
    double y01;
    double y02;
    double propagationAngle1;
    double propagationAngle2;
    double propagationAnglePhi1;
    double propagationAnglePhi2;
    bool isCylindric;
    bool is3D;
    int symmetryType;

    //spectral/temporal field properties
    double pulseEnergy1;
    double pulseEnergy2;
    double frequency1;
    double frequency2;
    int sgOrder1;
    int sgOrder2;
    double bandwidth1;
    double bandwidth2;
    double cephase1;
    double cephase2;
    double delay1;
    double delay2;
    double gdd1;
    double gdd2;
    double tod1;
    double tod2;
    double phaseMaterialThickness1;
    double phaseMaterialThickness2;
    int phaseMaterialIndex1;
    int phaseMaterialIndex2;

    //loaded FROG/EOS fields
    std::complex<double>* loadedField1;
    std::complex<double>* loadedField2;
    bool field1IsAllocated = 0;
    bool field2IsAllocated = 0;
    int pulse1FileType;
    int pulse2FileType;
    char field1FilePath[MAX_LOADSTRING];
    char field2FilePath[MAX_LOADSTRING];


    //polarization properties
    double polarizationAngle1;
    double polarizationAngle2;
    double circularity1;
    double circularity2;


    int pulsetype;
    double* InfoVec;
    double* ExtOut;
    std::complex<double>* EkwOut;
    double* totalSpectrum;
    int* imdone;
    int memoryError;
    int assignedGPU;
    int plotSim;
    crystalEntry* crystalDatabase;
    int batchIndex;
    int batchIndex2;
    double batchDestination;
    double batchDestination2;
    char outputBasePath[MAX_LOADSTRING];
    int runType;

    //sequence
    bool isInSequence;
    bool isFollowerInSequence;
    char sequenceString[MAX_LOADSTRING];
    double sequenceArray[MAX_LOADSTRING];
    int Nsequence;

    //fitting
    bool isInFittingMode;
    char fittingString[1024];
    char fittingPath[1024];
    double fittingArray[1024];
    double fittingPrecision;
    double* fittingReference;
    int Nfitting;
    int fittingMode;
    int fittingMaxIterations;
    size_t fittingROIstart;
    size_t fittingROIstop;
    size_t fittingROIsize;

} simulationParameterSet;

typedef struct cudaParameterSet {
    cudaStream_t CUDAStream;

    cuDoubleComplex* workspace1;
    cuDoubleComplex* workspace2;
    cuDoubleComplex* workspace2P;
	cuDoubleComplex* gridETemp1;
	cuDoubleComplex* gridETemp2;
	cuDoubleComplex* gridEFrequency1;
	cuDoubleComplex* gridEFrequency2;
	cuDoubleComplex* gridPropagationFactor1;
	cuDoubleComplex* gridPropagationFactor1Rho1;
	cuDoubleComplex* gridPropagationFactor1Rho2;
	cuDoubleComplex* gridPolarizationFactor1;
	cuDoubleComplex* gridPolarizationFrequency1;
	cuDoubleComplex* gridPropagationFactor2;
	cuDoubleComplex* gridPolarizationFactor2;
	cuDoubleComplex* gridPolarizationFrequency2;
	cuDoubleComplex* gridEFrequency1Next1;
	cuDoubleComplex* gridEFrequency1Next2;
	cuDoubleComplex* gridPlasmaCurrentFrequency1;
	cuDoubleComplex* gridPlasmaCurrentFrequency2;
    cuDoubleComplex* chiLinear1;
    cuDoubleComplex* chiLinear2;
	cuDoubleComplex* k1;
	cuDoubleComplex* k2;

    double* gridRadialLaplacian1;
    double* gridRadialLaplacian2;
    double* gridETime1;
    double* gridETime2;
    double* gridPolarizationTime1;
    double* gridPolarizationTime2;
    double* expGammaT;
    double* gridPlasmaCurrent1;
    double* gridPlasmaCurrent2;

    //fixed length arrays
    double firstDerivativeOperation[6];
    double plasmaParameters[6]; //[dt^2 * e^2/m * nonlinearAbsorptionStrength, gamma] 
	double chi2Tensor[9];
	double* chi3Tensor;
	double absorptionParameters[6];
    double rotationForward[9];
    double rotationBackward[9];
	int nonlinearSwitches[4];
	
    cufftHandle fftPlanZ2D;
    cufftHandle fftPlanD2Z;
    cufftHandle doublePolfftPlan;
    bool isCylindric;
    bool is3D;
    bool hasPlasma;
    bool isNonLinear;
    bool isUsingMillersRule;
	size_t Ntime;
    size_t Nfreq;
	size_t Nspace;
    size_t Nspace2;
	size_t Ngrid;
    size_t NgridC;
	int axesNumber;
	int sellmeierType;
	double f0;
	double fStep;
	double dt;
	double dx;
	double h;
	size_t Nsteps;
	int Nthread;
    int NblockC;
	int Nblock;
} cudaParameterSet;

unsigned long	solveNonlinearWaveEquation(void* lpParam);
int				runRK4Step(cudaParameterSet* sH, cudaParameterSet* sD, uint8_t stepNumber);
int				prepareElectricFieldArrays(simulationParameterSet* s, cudaParameterSet* sc);
int				calcEffectiveChi2Tensor(double* defftensor, double* dtensor, double theta, double phi);
int				fftshiftZ(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int             fftshiftD2Z(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int				fftshiftAndFilp(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
std::complex<double> sellmeier(std::complex<double>* ne, std::complex<double>* no, double* a, double f, double theta, double phi, int type, int eqn);
int				preparePropagation2DCartesian(simulationParameterSet* s, cudaParameterSet sc);
int				preparePropagation3DCylindric(simulationParameterSet* s, cudaParameterSet sc);
int             preparePropagation3D(simulationParameterSet* s, cudaParameterSet sc);
int				loadFrogSpeck(char* frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel, int fieldIndex);
int				rotateField(simulationParameterSet* s, double rotationAngle);
double          cModulusSquared(std::complex<double>complexNumber);
int             allocateGrids(simulationParameterSet* sCPU);
int             readCrystalDatabase(crystalEntry* db);
int             readSequenceString(simulationParameterSet* sCPU);
int             configureBatchMode(simulationParameterSet* sCPU);
int             saveDataSet(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, char* outputbase, bool saveInputs);
int             resolveSequence(int currentIndex, simulationParameterSet* s, crystalEntry* db);
int             readInputParametersFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, char* filePath);
int             loadPulseFiles(simulationParameterSet* sCPU);
unsigned long   solveNonlinearWaveEquationSequence(void* lpParam);
int             saveSettingsFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr);
void            unixNewLine(FILE* iostream);
int             saveSlurmScript(simulationParameterSet* sCPU, int gpuType, int gpuCount);
int             loadSavedFields(simulationParameterSet* sCPU, char* outputBase, bool GPUisPresent);
int             getTotalSpectrum(simulationParameterSet* sCPU, cudaParameterSet* sc);
unsigned long   runFitting(simulationParameterSet* sCPU);
void            runFittingIteration(int* m, int* n, double* fittingValues, double* fittingFunction);
int             readFittingString(simulationParameterSet* sCPU);
int             loadReferenceSpectrum(char* spectrumPath, simulationParameterSet* sCPU);
int             removeCharacterFromString(char* cString, size_t N, char removedChar);
int             skipFileUntilCharacter(FILE* fstream, char target);
int             applyLinearPropagation(simulationParameterSet* s, int materialIndex, double thickness);
int             fillRotationMatricies(simulationParameterSet* sCPU, cudaParameterSet* s);
int             deallocateCudaParameterSet(cudaParameterSet* s);
int             initializeCudaParameterSet(simulationParameterSet* sCPU, cudaParameterSet* s);