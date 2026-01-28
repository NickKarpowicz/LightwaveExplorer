#pragma once
#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <complex>
#include <cstdint>
#include <concepts>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#ifdef __linux__
#include <unistd.h>
#endif
#include "LightwaveExplorerHelpers.h"

template<typename T>
concept FPType = std::same_as<T, float> || std::same_as<T, double>;

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

template <FPType T>
class maxwellPoint {
public:
    T x{};
    T y{};
    T z{};
    hostOrDevice inline T& operator()(const int index) {
        switch (index) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: return x;
        }
    }
    hostOrDevice inline void operator+=(
        const maxwellPoint<T>& other) {
        x += other.x;
        y += other.y;
        z += other.z;
    }
    hostOrDevice inline void operator-=(
        const maxwellPoint<T>& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
    }
    hostOrDevice inline void operator*=(
        const maxwellPoint<T>& other) {
        x *= other.x;
        y *= other.y;
        z *= other.z;
    }
    hostOrDevice inline void operator/=(
        const maxwellPoint<T>& other) {
        x /= other.x;
        y /= other.y;
        z /= other.z;
    }

    hostOrDevice inline void operator+=(
        const T other) {
        x += other;
        y += other;
        z += other;
    }
    hostOrDevice inline void operator-=(
        const T other) {
        x -= other;
        y -= other;
        z -= other;
    }
    hostOrDevice inline void operator*=(
        const T other) {
        x *= other;
        y *= other;
        z *= other;
    }
    hostOrDevice inline void operator/=(
        const T other) {
        x /= other;
        y /= other;
        z /= other;
    }

    hostOrDevice inline maxwellPoint<T> operator*(
        const maxwellPoint<T>& other) const {
        return maxwellPoint<T>{
                x * other.x,
                y * other.y,
                z * other.z};
    }
    hostOrDevice inline maxwellPoint<T> operator*(
        const T other) const {
        return maxwellPoint<T>{
                x * other,
                y * other,
                z * other};
    }
    hostOrDevice inline friend maxwellPoint<T> operator*(const T a, const maxwellPoint<T>& b) {
        return maxwellPoint<T>{
                a * b.x,
                a * b.y,
                a * b.z};
    }

    hostOrDevice inline maxwellPoint<T> operator/(
        const maxwellPoint<T>& other) const {
        return maxwellPoint<T>{
                x / other.x,
                y / other.y,
                z / other.z};
    }
    hostOrDevice inline maxwellPoint<T> operator/(
        const T other) const {
        return maxwellPoint<T>{
                x / other,
                y / other,
                z / other};
    }
    hostOrDevice inline friend maxwellPoint<T> operator/(const T a, const maxwellPoint<T>& b) {
        return maxwellPoint<T>{
                a / b.x,
                a / b.y,
                a / b.z};
    }

    hostOrDevice inline maxwellPoint<T> operator+(
        const maxwellPoint<T>& other) const {
        return maxwellPoint<T>{
                x + other.x,
                y + other.y,
                z + other.z};
    }
    hostOrDevice inline maxwellPoint<T> operator+(
        const T other) const {
        return maxwellPoint<T>{
                x + other,
                y + other,
                z + other};
    }
    hostOrDevice inline friend maxwellPoint<T> operator+(const T a, const maxwellPoint<T>& b) {
        return maxwellPoint<T>{
                a + b.x,
                a + b.y,
                a + b.z};
    }

    hostOrDevice inline maxwellPoint<T> operator-(
        const maxwellPoint<T>& other) const {
        return maxwellPoint<T>{
                x - other.x,
                y - other.y,
                z - other.z};
    }
    hostOrDevice inline maxwellPoint<T> operator-(
        const T other) const {
        return maxwellPoint<T>{
                x - other,
                y - other,
                z - other};
    }
    hostOrDevice inline friend maxwellPoint<T> operator-(const T a, const maxwellPoint<T>& b) {
        return maxwellPoint<T>{
                a - b.x,
                a - b.y,
                a - b.z};
    }
};

template <typename T>
class maxwellKPoint {
public:
    maxwellPoint<T> kE;
    maxwellPoint<T> kH;
};

template <typename T>
class oscillator {
public:
    maxwellPoint<T> J;
    maxwellPoint<T> P;

    //note that I only defined the three operations I need rather than the
    //full set. Might be worth filling in everything later.
    hostOrDevice void operator+=(
        const oscillator<T>& other) {
        J += other.J;
        P += other.P;
    }
    hostOrDevice inline oscillator<T> operator*(
        const T other) const {
        return oscillator<T>{
                J * other,
                P * other,
        };
    }
    hostOrDevice inline oscillator<T> operator+(
        const oscillator<T>& other) const {
        return oscillator<T>{
                J + other.J,
                P + other.P
        };
    }
};

template<FPType T>
class PlasmaParameters{
    public:
    T nonlinearAbsorption = {};
    T bandgap = {};
    T drudeGamma = {};
    T effectiveMass = {};
    T initialDensity = {};
    T integrationFactor = {};
    T energyFactor = {};
    int fieldExponent = {};

    PlasmaParameters() = default;

    template<FPType otherFP>
    PlasmaParameters(const PlasmaParameters<otherFP>& other){
        nonlinearAbsorption = static_cast<T>(other.nonlinearAbsorption);
        bandgap = static_cast<T>(other.bandgap);
        drudeGamma = static_cast<T>(other.drudeGamma);
        effectiveMass = static_cast<T>(other.effectiveMass);
        initialDensity = static_cast<T>(other.initialDensity);
        integrationFactor = static_cast<T>(other.integrationFactor);
        energyFactor = static_cast<T>(other.energyFactor);
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
    bool visualizationOnly = false;
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

    void resize(size_t newCount){
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
template <FPType T, typename deviceComplex>
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
    T* gridBiaxialDelta = 0;
    deviceComplex* chiLinear1 = 0;
    deviceComplex* chiLinear2 = 0;
    T* inverseChiLinear1 = 0;
    T* inverseChiLinear2 = 0;
    T* fieldFactor1 = 0;
    T* fieldFactor2 = 0;
    deviceComplex* k1 = 0;
    deviceComplex* k2 = 0;
    deviceComplex n0 = 0.0;
    T* gridRadialLaplacian1 = 0;
    T* gridRadialLaplacian2 = 0;
    T* gridETime1 = 0;
    T* gridETime2 = 0;
    T* gridPolarizationTime1 = 0;
    T* gridPolarizationTime2 = 0;
    T* expGammaT = 0;
    T* gridPlasmaCurrent1 = 0;
    T* gridPlasmaCurrent2 = 0;

    //fixed length arrays
    PlasmaParameters<T> plasmaParameters = {}; //[dt^2 * e^2/m * nonlinearAbsorptionStrength, gamma]
    T chi2Tensor[18] = { 0 };
    T chi3Tensor[81] = { 0 };
    T rotationForward[9] = { 0 };
    T rotationBackward[9] = { 0 };
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
    T fftNorm = 0;
    int axesNumber = 0;
    int sellmeierType = 0;
    T crystalTheta;
    T crystalPhi;
    T f0 = 0;
    T fStep = 0;
    T dt = 0;
    T dx = 0;
    T dk1 = 0;
    T dk2 = 0;
    T h = 0;
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
            if(fs.good())db.push_back(newEntry);
            std::getline(fs, line); //~~~crystal end~~~
        }
    }
};

enum class BeamBasis {
    laguerre,
    hermite
};

template <FPType T, int number_of_modes, int max_expansion_order>
struct BeamSpecification {
    BeamBasis basis = BeamBasis::hermite;
    int relevant_modes = 1;
    int relevant_expansion = 1;
    std::uint8_t m[number_of_modes] = {};
    std::uint8_t l[number_of_modes] = {};
    T weight[number_of_modes] = { 1 };
    T phase[number_of_modes] = {};
    T waist[number_of_modes][max_expansion_order] = {};
    T rotation[number_of_modes][max_expansion_order] = {};
    T x_offset[number_of_modes][max_expansion_order] = {};
    T y_offset[number_of_modes][max_expansion_order] = {};
    T z_offset[number_of_modes][max_expansion_order] = {};
    T angle_x[number_of_modes][max_expansion_order] = {};
    T angle_y[number_of_modes][max_expansion_order] = {};

    BeamSpecification() = default;

    template<FPType otherFP>
    BeamSpecification(const BeamSpecification<otherFP, number_of_modes, max_expansion_order>& other){
        basis = other.basis;
        relevant_modes = other.relevant_modes;
        relevant_expansion = other.relevant_expansion;
        for(int i = 0; i < number_of_modes; i++){
            m[i] = other.m[i];
            l[i] = other.l[i];
            weight[i] = static_cast<T>(other.weight[i]);
            phase[i] = static_cast<T>(other.phase[i]);
            for(int j = 0; j< max_expansion_order; j++){
                waist[i][j] = static_cast<T>(other.waist[i][j]);
                rotation[i][j] = static_cast<T>(other.rotation[i][j]);
                x_offset[i][j] = static_cast<T>(other.x_offset[i][j]);
                y_offset[i][j] = static_cast<T>(other.y_offset[i][j]);
                z_offset[i][j] = static_cast<T>(other.z_offset[i][j]);
                angle_x[i][j] = static_cast<T>(other.angle_x[i][j]);
                angle_y[i][j] = static_cast<T>(other.angle_x[i][j]);
            }
        }
    }

    template<FPType otherFP>
    BeamSpecification& operator=(const BeamSpecification<otherFP, number_of_modes, max_expansion_order>& other){
        basis = other.basis;
        relevant_modes = other.relevant_modes;
        relevant_expansion = other.relevant_expansion;
        for(int i = 0; i < number_of_modes; i++){
            m[i] = other.m[i];
            l[i] = other.l[i];
            weight[i] = static_cast<T>(other.weight[i]);
            phase[i] = static_cast<T>(other.phase[i]);
            for(int j = 0; j< max_expansion_order; j++){
                waist[i][j] = static_cast<T>(other.waist[i][j]);
                rotation[i][j] = static_cast<T>(other.rotation[i][j]);
                x_offset[i][j] = static_cast<T>(other.x_offset[i][j]);
                y_offset[i][j] = static_cast<T>(other.y_offset[i][j]);
                z_offset[i][j] = static_cast<T>(other.z_offset[i][j]);
                angle_x[i][j] = static_cast<T>(other.angle_x[i][j]);
                angle_y[i][j] = static_cast<T>(other.angle_y[i][j]);
            }
        }
    }

    BeamSpecification(const std::string& descriptor, const BeamBasis b){
        basis = b;
        std::vector<std::vector<T>> data = parse_string_to_vecs<T>(descriptor);
        std::cout << "modes: " << data.size() << '\n';
        std::cout << "elems: " << data[0].size() << std::endl;
        relevant_modes = std::min(number_of_modes, static_cast<int>(data.size()));
        //TODO: validation: should have minimum element number for valid beam
        for(int mode_idx=0; mode_idx<relevant_modes; mode_idx++){
            int current_expansion = std::min((static_cast<int>(data[mode_idx].size()) - 2) / 6, max_expansion_order);
            relevant_expansion = std::max(relevant_expansion, current_expansion);
            l[mode_idx] = static_cast<std::uint8_t>(data[mode_idx][0]);
            m[mode_idx] = static_cast<std::uint8_t>(data[mode_idx][1]);
            weight[mode_idx] = static_cast<T>(data[mode_idx][2]);
            phase[mode_idx] = static_cast<T>(data[mode_idx][3]);
            for(int expansion_idx=0; expansion_idx<current_expansion; expansion_idx++){
                waist[mode_idx][expansion_idx] = 1e-6f * data[mode_idx][4 + 6*expansion_idx];
                rotation[mode_idx][expansion_idx] = data[mode_idx][4 + 6*expansion_idx + 1];
                std::cout << "waist: " << waist[mode_idx][expansion_idx] << '\n';
                x_offset[mode_idx][expansion_idx] = 1e-6f * data[mode_idx][4 + 6*expansion_idx+2];
                y_offset[mode_idx][expansion_idx] = 1e-6f * data[mode_idx][4 + 6*expansion_idx+3];
                z_offset[mode_idx][expansion_idx] = 1e-6f * data[mode_idx][4 + 6*expansion_idx+4];
                angle_x[mode_idx][expansion_idx] = data[mode_idx][4 + 6*expansion_idx+5];
                angle_y[mode_idx][expansion_idx] = data[mode_idx][4 + 6*expansion_idx+6];
            }
        }


    }
};

//templated class for describing a pulse in various floating point representations
//copyable between representations (required for strict FP32 mode)
template <FPType T>
class Pulse {
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
    T polarizationAngle;
    T circularity;
    T pulseSum;
    BeamSpecification<T, 16, 4> beam_spec;
    T beamwaist;
    T x_offset;
    T y_offset;
    T z_offset;
    T angle_x_offset;
    T angle_y_offset;
    Pulse() : energy(),
        frequency(),
        bandwidth(),
        sgOrder(),
        cep(),
        delay(),
        gdd(),
        tod(),
        phaseMaterial(),
        phaseMaterialThickness(),
        polarizationAngle(),
        circularity(),
        pulseSum(),
        beam_spec(),
        beamwaist(),
        x_offset(),
        y_offset(),
        z_offset(),
        angle_x_offset(),
        angle_y_offset() {}

    template<FPType U>
    Pulse(Pulse<U>& other) : energy((T)other.energy),
        frequency((T)other.frequency),
        bandwidth((T)other.bandwidth),
        sgOrder(other.sgOrder),
        cep((T)other.cep),
        delay((T)other.delay),
        gdd((T)other.gdd),
        tod((T)other.tod),
        phaseMaterial(other.phaseMaterial),
        phaseMaterialThickness((T)other.phaseMaterialThickness),
        polarizationAngle((T)other.polarizationAngle),
        circularity((T)other.circularity),
        pulseSum((T)other.pulseSum),
        beam_spec(other.beam_spec),
        beamwaist((T)other.beamwaist),
        x_offset((T)other.x_offset),
        y_offset((T)other.y_offset),
        z_offset((T)other.z_offset),
        angle_x_offset((T)other.angle_x_offset),
        angle_y_offset((T)other.angle_y_offset) {}

    template <FPType U>
    Pulse& operator=(const Pulse<U>& other) {
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
        polarizationAngle = (T)other.polarizationAngle;
        circularity = (T)other.circularity;
        pulseSum = (T)other.pulseSum;
        beam_spec = other.beam_spec;
        beamwaist = (T)other.beamwaist;
        x_offset = (T)other.x_offset;
        y_offset = (T)other.y_offset;
        z_offset = (T)other.z_offset;
        angle_x_offset = (T)other.angle_x_offset;
        angle_y_offset = (T)other.angle_y_offset;
        return *this;
    }
};
