#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <complex>
#include <cuComplex.h>
#include "cufft.h"

//Simulation parameter struct to pass to the simulations running in threads
struct simulationParameterSet {
    double rStep;
    double tStep;
    double fStep;
    double kStep;
    double propagationStep;
    size_t Npropagation;
    size_t Ntime;
    size_t Nspace;
    size_t Ngrid;
    size_t Nsims;
    double spatialWidth;
    double timeSpan;
    int materialIndex;
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
    std::complex<double>* refractiveIndex1;
    std::complex<double>* refractiveIndex2;
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
    double propagationAngle1;
    double propagationAngle2;
    bool isCylindric;
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

    //loaded FROG/EOS fields
    std::complex<double>* loadedField1;
    std::complex<double>* loadedField2;
    bool field1IsAllocated = 0;
    bool field2IsAllocated = 0;



    //polarization properties
    double polarizationAngle1;
    double polarizationAngle2;
    double circularity1;
    double circularity2;


    int pulsetype;
    double* InfoVec;
    std::complex<double>* Ext;
    std::complex<double>* Ekw;
    std::complex<double>* ExtOut;
    std::complex<double>* EkwOut;
    int* imdone;
    int memoryError;
    int plotSim;
    struct crystalEntry* crystalDatabase;
    int batchIndex;
    double batchDestination;
    char* outputBasePath;

    //sequence
    bool isInSequence;
    bool isFollowerInSequence;
    char* sequenceString;
    double* sequenceArray;
    int Nsequence;
};

struct crystalEntry {
    wchar_t crystalNameW[256];
    int axisType;
    int sellmeierType;
    int nonlinearSwitches[4];
    double sellmeierCoefficients[66];
    wchar_t sellmeierReference[512];
    double d[18];
    wchar_t dReference[512];
    double chi3[48];
    wchar_t chi3Reference[512];
    double absorptionParameters[6];
    wchar_t spectralFile[512];
    double spectralData[2048];
    int numberOfEntries;
};

struct cudaParameterSet {
	cuDoubleComplex* gridETime1;
	cuDoubleComplex* gridETime2;
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
	cuDoubleComplex* gridRadialLaplacian1;
	cuDoubleComplex* gridRadialLaplacian2;
	cuDoubleComplex* gridPlasmaCurrentFrequency1;
	cuDoubleComplex* gridPlasmaCurrentFrequency2;
	cuDoubleComplex* k1;
	cuDoubleComplex* k2;

	cuDoubleComplex* ne;
	cuDoubleComplex* no;
	bool isCylindric;
	bool hasPlasma;
	double* gridPolarizationTime1; //future optimization: this could be double rather than complex, if I can figure out the data layout of D2Z cufft
	double* gridPolarizationTime2;
	double* firstDerivativeOperation;
	double* plasmaParameters; //[dt^2 * e^2/m * nonlinearAbsorptionStrength, gamma] 
	double* chi2Tensor;
	double* chi3Tensor;
	double* gridPlasmaCurrent1;
	double* gridPlasmaCurrent2;
	double* absorptionParameters;
	double* expGammaT;
	int* nonlinearSwitches;
	long long* propagationInts;
	cufftHandle fftPlan;
	cufftHandle polfftPlan;
	bool isNonLinear;
	size_t Ntime;
	size_t Nspace;
	size_t Ngrid;
	int axesNumber;
	int sellmeierType;
	double f0;
	double fStep;
	double dt;
	double dx;
	double h;
	long long Nsteps;
	int Nthread;
	int Nblock;
};

unsigned long	solveNonlinearWaveEquation(void* lpParam);
int				runRK4Step(struct cudaParameterSet s, int stepNumber);
int				prepareElectricFieldArrays(struct simulationParameterSet* s, struct cudaParameterSet* sc);
int				calcEffectiveChi2Tensor(double* defftensor, double* dtensor, double theta, double phi);
int				fftshiftZ(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
int				fftshiftAndFilp(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2);
std::complex<double> sellmeier(std::complex<double>* ne, std::complex<double>* no, double* a, double f, double theta, double phi, int type, int eqn);
int				preparePropagation2DCartesian(struct simulationParameterSet* s, struct cudaParameterSet sc);
int				preparePropagation3DCylindric(struct simulationParameterSet* s, struct cudaParameterSet sc);
double			findWalkoffAngles(struct simulationParameterSet* s, double dk, double f, double tol);
int				loadFrogSpeck(char* frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel, int fieldIndex);
int				plotDataXY(double* X, double* Y, double minX, double maxX, double minY, double maxY, int N, int plotSizeX, int plotSizeY, double lineWidth, double markerWidth, double* plotGrid, double* xTicks, int NxTicks, double* yTicks, int NyTicks);
int				rotateField(struct simulationParameterSet* s, double rotationAngle);
double          cModulusSquared(std::complex<double>complexNumber);
int             allocateGrids(struct simulationParameterSet* sCPU);
int             readCrystalDatabase(struct crystalEntry* db);
int             readSequenceString(struct simulationParameterSet* sCPU);
int             configureBatchMode(struct simulationParameterSet* sCPU);
int             saveDataSet(struct simulationParameterSet* sCPU, struct crystalEntry* crystalDatabasePtr, char* outputbase);
int             resolveSequence(int currentIndex, struct simulationParameterSet* s, struct crystalEntry* db);
int             readInputParametersFile(struct simulationParameterSet* sCPU, struct crystalEntry* crystalDatabasePtr, char* filePath);