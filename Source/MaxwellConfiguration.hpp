#pragma once
#include "DataStructures.hpp"
#include "LightwaveExplorerInterfaceClasses.hpp"
constexpr int64_t NmaterialMax = 8;
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
    char* materialMap{};
    int64_t* oscillatorIndexMap{};
    int64_t materialKeys[NmaterialMax]{};
    deviceFP materialTheta[NmaterialMax]{};
    deviceFP materialPhi[NmaterialMax]{};

    maxwellPoint<deviceFP> sellmeierEquations[22][NmaterialMax]{};
    maxwellPoint<deviceFP> chi3[27][NmaterialMax]{};
    maxwellPoint<deviceFP> chi2[6][NmaterialMax]{};
    int nonlinearAbsorptionOrder[NmaterialMax]{};
    deviceFP kNonlinearAbsorption[NmaterialMax]{};
    deviceFP kDrude[NmaterialMax]{};
    deviceFP gammaDrude[NmaterialMax]{};
    deviceFP kCarrierGeneration[NmaterialMax]{};
    deviceFP startingCarriers[NmaterialMax]{};
    deviceFP rotateForward[9][NmaterialMax]{};
    deviceFP rotateBackward[9][NmaterialMax]{};
    bool hasChi2[NmaterialMax]{};
    bool hasFullChi3[NmaterialMax]{};
    bool hasSingleChi3[NmaterialMax]{};
    bool hasPlasma[NmaterialMax]{};
    bool hasMaterialMap{};
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
    int64_t fftSize{};
    int64_t Ninjection{};
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
        Ninjection = NtIO * timeFactor;
        fftSize = NtIO / 2 + 1;
        tGridFactor = timeFactor;
        Ngrid = Nz * Ny * Nx;
        fillRotationMatricies((*s).crystalTheta, (*s).crystalPhi, 0);
    }
};