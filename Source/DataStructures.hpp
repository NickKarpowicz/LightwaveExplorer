#pragma once
#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <complex>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#ifdef __linux__
#include <unistd.h>
#endif
#include "LightwaveExplorerHelpers.h"

//Forward declarations of interface classes
class simulationParameterSet;
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

enum class runTypes : int {
    counter = -1,
    normal = 0,
    cluster = 1
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

template<typename deviceFP>
class PlasmaParameters{
    public:
    deviceFP nonlinearAbsorption = {};
    deviceFP bandgap = {};
    deviceFP drudeGamma = {};
    deviceFP effectiveMass = {};
    deviceFP initialDensity = {};
    deviceFP integrationFactor = {};
    deviceFP energyFactor = {};
    int fieldExponent = {};

    PlasmaParameters() = default;

    template<typename otherFP>
    PlasmaParameters(const PlasmaParameters<otherFP>& other){
        nonlinearAbsorption = static_cast<deviceFP>(other.nonlinearAbsorption);
        bandgap = static_cast<deviceFP>(other.bandgap);
        drudeGamma = static_cast<deviceFP>(other.drudeGamma);
        effectiveMass = static_cast<deviceFP>(other.effectiveMass);
        initialDensity = static_cast<deviceFP>(other.initialDensity);
        integrationFactor = static_cast<deviceFP>(other.integrationFactor);
        energyFactor = static_cast<deviceFP>(other.energyFactor);
        fieldExponent = other.fieldExponent;
    }
};

class NonlinearPropertyFlags{
    public:
    bool hasChi2 = false;
    bool hasChi3 = false;
    bool assumeCentrosymmetric = false;
};

class LWEDevice{
    public:
    simulationParameterSet* cParams;
	int memoryStatus = -1;
	bool configuredFFT = false;
    virtual int deviceCalloc(void** ptr, size_t N, size_t elementSize) = 0;
    virtual void deviceMemset(void* ptr, int value, size_t count) = 0;
    virtual void deviceMemcpyImplementation(
		void* dst, 
		const void* src, 
		size_t count, 
		copyType kind) = 0;
    virtual void deviceFree(void* block) = 0;
    virtual void fft(const void* input, void* output, deviceFFT type) = 0;
    virtual void reset(simulationParameterSet* sCPU) = 0;
    
    void deviceMemcpy(
		void* dst, 
		const void* src, 
		size_t count, 
		copyType kind) {
			deviceMemcpyImplementation(dst, src, count, kind);
	}

    void deviceMemcpy(
		void* dst, 
		void* src, 
		size_t count, 
		copyType kind) {
            if(src==dst) return;
			deviceMemcpyImplementation(dst, src, count, kind);
	}

	void deviceMemcpy(
		double* dst, 
		const float* src, 
		size_t count, 
		copyType kind) {
            float* copyBuffer = new float[count / sizeof(double)];
            deviceMemcpyImplementation(
                static_cast<void*>(copyBuffer), 
                static_cast<const void*>(src), 
                count/2, 
                kind);
            for (size_t i = 0; i < count / sizeof(double); i++) {
                dst[i] = copyBuffer[i];
            }
            delete[] copyBuffer;
	}
	
	void deviceMemcpy(
		float* dst, 
		const double* src, 
		size_t count, 
		copyType kind) {
            float* copyBuffer = new float[count / sizeof(double)];
            for (size_t i = 0; i < count / sizeof(double); i++) {
                copyBuffer[i] = (float)src[i];
            }
            deviceMemcpyImplementation(
                static_cast<void*>(dst), 
                static_cast<void*>(copyBuffer), 
                count / 2, 
                kind);
            delete[] copyBuffer;
	}

    void deviceMemcpy(
		std::complex<double>* dst, 
		const std::complex<float>* src, 
		size_t count, 
		copyType kind) {
		std::complex<float>* copyBuffer = 
			new std::complex<float>[count / sizeof(std::complex<double>)];
		deviceMemcpyImplementation(
			static_cast<void*>(copyBuffer), 
			static_cast<const void*>(src), 
			count/2, 
			kind);
		for (size_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			dst[i] = std::complex<double>(copyBuffer[i].real(), copyBuffer[i].imag());
		}
		delete[] copyBuffer;
	}

	void deviceMemcpy(
		std::complex<float>* dst, 
		const std::complex<double>* src, 
		size_t count, 
		copyType kind) {
		std::complex<float>* copyBuffer = 
			new std::complex<float>[count / sizeof(std::complex<double>)];
		
		for (size_t i = 0; i < count / sizeof(std::complex<double>); i++) {
			copyBuffer[i] = 
				std::complex<float>((float)src[i].real(), (float)src[i].imag());
		}
		deviceMemcpyImplementation(
			static_cast<void*>(dst), 
			static_cast<void*>(copyBuffer), 
			count / 2, 
			kind);
		delete[] copyBuffer;
	}
};

template<typename T>
class LWEBuffer{
    LWEDevice* d = nullptr;

public:
    T* buffer = nullptr;
    size_t count = 0;
    size_t bytes = 0;
    LWEBuffer(){}
    LWEBuffer(const LWEBuffer&) = delete;
    LWEBuffer& operator=(const LWEBuffer&) = delete;
    LWEBuffer(LWEDevice* d, const size_t N, const size_t elementSize = sizeof(T)) : 
    d(d),
    count(maxN(N,static_cast<size_t>(1))),
    bytes(count*elementSize)
    {
        int error = d->deviceCalloc((void**)&buffer, count, elementSize);
        if(error){
            throw std::runtime_error("Couldn't allocate enough memory.");
        }
    }
    LWEBuffer(LWEDevice* d, const std::vector<T> v) : 
    d(d),
    count(v.size()),
    bytes(count * sizeof(T)){
        int error = d->deviceCalloc((void**)&buffer, count, sizeof(T));
        if(error){
            throw std::runtime_error("Couldn't allocate enough memory.");
        }
        d->deviceMemcpy(buffer, v.data(), count, copyType::ToDevice);
    }
    
    ~LWEBuffer(){
        if(bytes) d->deviceFree(buffer);
    }
    T* device_ptr() const {
        if(bytes == 0 || buffer == nullptr) throw std::runtime_error("Attempted to access empty LWEBuffer");
        return buffer;
    }

    void initialize_to_zero(){
        if(bytes) d->deviceMemset(buffer, 0, bytes);
    }

    void resize(int64_t newCount){
        if(newCount == count && buffer != nullptr){
            initialize_to_zero();
            return;
        }
        if(buffer != nullptr) d->deviceFree(buffer);
        int error = d->deviceCalloc((void**)&buffer, newCount, sizeof(T));
        if(error){
            throw std::runtime_error("Resize of LWE buffer failed with error");
        }
        count = newCount;
        bytes = count * sizeof(T);
    }

    void deallocate(){
        if(buffer != nullptr) d->deviceFree(buffer);
        buffer = nullptr;
        count = 0;
        bytes = 0;
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
    deviceFP* gridBiaxialDelta = 0;
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
    PlasmaParameters<deviceFP> plasmaParameters = {}; //[dt^2 * e^2/m * nonlinearAbsorptionStrength, gamma] 
    deviceFP chi2Tensor[18] = { 0 };
    deviceFP chi3Tensor[81] = { 0 };
    deviceFP rotationForward[9] = { 0 };
    deviceFP rotationBackward[9] = { 0 };
    NonlinearPropertyFlags nonlinearSwitches{};

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
    NonlinearPropertyFlags nonlinearSwitches = {};
    std::array<double,66> sellmeierCoefficients = {};
    std::string sellmeierReference;
    std::array<double,18> d = {};
    std::string dReference;
    std::array<double,81> chi3 = {};
    std::string chi3Reference;
    std::string spectralFile;
    std::array<double,7> nonlinearReferenceFrequencies = {};
    std::array<double,132> offDiagonalCoefficients = {};
};

//Crystal database class; primarily holds a std::vector of crystalEntry elements
//comprising the database, plus method for loading the database from the file
class crystalDatabase {
public:
    std::vector<crystalEntry> db;
    std::string path = "CrystalDatabase.txt";
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
        else{
            path = databasePath;
        }
#elif defined __linux__

        std::string databasePath(std::getenv("HOME"));
        databasePath.append("/.LightwaveExplorer/CrystalDatabase.txt");
        std::ifstream fs(databasePath);

        if(!fs.is_open()){
            char pBuf[256];
            int64_t len = sizeof(pBuf); 
            int bytes = minN(readlink("/proc/self/exe", pBuf, len), len - 1);
            if(bytes >= 0)
                pBuf[bytes] = '\0';
            std::string binPath(pBuf);
            int64_t posPath = binPath.find_last_of("/");
            databasePath = binPath.substr(0, posPath).append("/../share/LightwaveExplorer/CrystalDatabase.txt");
            fs.open(databasePath);
            if (!fs.is_open()) {
                fs.open("CrystalDatabase.txt");
            }
            else{
                path = databasePath;
            }
        }
        
#else
        std::ifstream fs("CrystalDatabase.txt");
#endif
        loadFromFilestream(fs);
    }
    crystalDatabase(std::string customPath){
        std::ifstream fs(customPath);
        loadFromFilestream(fs);
        path = customPath;
    }
    void loadFromFilestream(std::ifstream& fs){
        std::string line;
        if (!fs.is_open())return;
        int tempInt = 0;
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

            std::getline(fs, line); //Sellmeier reference:, or if it has Monoclinic values, read them
            if(line=="Off-diagonal susceptibility coefficients:"){
                for(int k = 0; k<132; ++k){
                    fs >> newEntry.offDiagonalCoefficients[k];
                }
                std::getline(fs,line);
                std::getline(fs,line);
            }
            std::getline(fs, line);
            newEntry.sellmeierReference = line;

            std::getline(fs, line); // chi2 type:
            fs >> tempInt;
            newEntry.nonlinearSwitches.hasChi2 = tempInt == 1;
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
            fs >> tempInt;
            std::getline(fs, line);
            std::getline(fs, line); //chi3:
            switch (tempInt) {
            case 0: //no chi3, skip all three lines
                std::getline(fs, line);
                std::getline(fs, line);
                std::getline(fs, line);
                newEntry.nonlinearSwitches.hasChi3 = false;
                newEntry.nonlinearSwitches.assumeCentrosymmetric = false;
                break;
            case 1: //read full chi3
                for (int k = 0; k < 81; ++k) {
                    fs >> newEntry.chi3[k]; //row-order!
                }
                std::getline(fs, line);
                newEntry.nonlinearSwitches.hasChi3 = true;
                newEntry.nonlinearSwitches.assumeCentrosymmetric = false;
                break;

            case 2: //assume centrosymmetric, just read chi3_1111, then skip
                fs >> newEntry.chi3[0];
                std::getline(fs, line);
                std::getline(fs, line);
                std::getline(fs, line);

                newEntry.nonlinearSwitches.hasChi3 = true;
                newEntry.nonlinearSwitches.assumeCentrosymmetric = true;
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


