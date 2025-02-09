#include "LightwaveExplorerTrilingual.h"
#include "LightwaveExplorerInterfaceClasses.hpp"
#include <dlib/global_optimization.h>
#include "Kernels.hpp"

using namespace kernelNamespace;

namespace hostFunctions{
	static simulationParameterSet* fittingSet;
	static ActiveDevice* dFit;

	//forward declarations when useful
	static unsigned long int solveFDTD(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		int64_t tFactor, 
		deviceFP dz, 
		deviceFP frontBuffer, 
		deviceFP backBuffer,
		deviceFP observationPoint = 0.0,
		int gridIndex = -1,
		bool preserveNearField = false,
		int64_t manualWaitFrames = -1);

	static std::complex<double> hostSellmeierFunc(
		double ls, 
		const double omega, 
		const double* a, 
		const int eqn) {
		const double omega2 = omega * omega;
		double realPart;
		std::complex<double> compPart;
		switch (eqn) {
		case 0:
			if (ls == -a[3] || ls == -a[6] || ls == -a[9] || ls == -a[12]) return std::complex<double>{};
			realPart = a[0]
				+ (a[1] + a[2] * ls) / (ls + a[3])
				+ (a[4] + a[5] * ls) / (ls + a[6])
				+ (a[7] + a[8] * ls) / (ls + a[9])
				+ (a[10] + a[11] * ls) / (ls + a[12])
				+ a[13] * ls
				+ a[14] * ls * ls
				+ a[15] * ls * ls * ls;
			compPart = kLorentzian<double>() * a[16] 
				/ std::complex<double>(a[17] - omega2, a[18] * omega)
				+ kLorentzian<double>() * a[19] 
				/ std::complex<double>(a[20] - omega2, a[21] * omega);
			return std::sqrt(maxN(realPart, 0.9f) + compPart);
		case 1:
			//up to 7 Lorentzian lines
			compPart = a[1] / std::complex<double>(a[2] - omega2, a[3] * omega)
				+ a[4] / std::complex<double>(a[5] - omega2, a[6] * omega)
				+ a[7] / std::complex<double>(a[8] - omega2, a[9] * omega)
				+ a[10] / std::complex<double>(a[11] - omega2, a[12] * omega)
				+ a[13] / std::complex<double>(a[14] - omega2, a[15] * omega)
				+ a[16] / std::complex<double>(a[17] - omega2, a[18] * omega)
				+ a[19] / std::complex<double>(a[20] - omega2, a[21] * omega);
			compPart *= kLorentzian<double>();
			compPart += a[0];
			return std::complex<double>((
				std::sqrt(compPart)).real(), 
				-std::abs((std::sqrt(compPart)).imag()));

		case 100:
		{
			if (ls == -a[3] 
				|| ls == -a[6] 
				|| ls == -a[9] 
				|| ls == -a[12]) return std::complex<double>{};
		}
		realPart = a[0]
			+ (a[1] + a[2] * ls) / (ls + a[3])
			+ (a[4] + a[5] * ls) / (ls + a[6])
			+ (a[7] + a[8] * ls) / (ls + a[9])
			+ (a[10] + a[11] * ls) / (ls + a[12])
			+ a[13] * ls
			+ a[14] * ls * ls
			+ a[15] * ls * ls * ls;
		//"real-valued equation has no business being < 1" - causality
		return std::complex<double>(std::sqrt(maxN(realPart, 0.9f)), 0.0f);
		}
		return std::complex<double>(1.0,0.0);
	};
	static int getTotalSpectrum(ActiveDevice& d) {
		simulationParameterSet* sCPU = d.cParams;
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;

		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.fft((*sc).gridETime1, (*sc).workspace1, deviceFFT::D2Z_1D);
		if ((*sc).is3D) {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum3DKernel{ d.dParamsDevice });
		}
		else if ((*sc).isCylindric) {
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel{ d.dParamsDevice });
		}
		else {
			//uncomment and change logic if I want to use the square spectra
			//d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrum2DSquareKernel, d.dParamsDevice);
			d.deviceLaunch((unsigned int)(*sCPU).Nfreq, 1u, totalSpectrumKernel{ d.dParamsDevice });
		}
		
		d.deviceMemcpy((double*)(*sCPU).totalSpectrum, 
			(deviceFP*)(*sc).gridPolarizationTime1, 
			3 * (*sCPU).Nfreq * sizeof(double), copyType::ToHost);

		//apply normalization to result of 3D calculation for numerical precision (value may not be
		//represtentable as a float)
		if ((*sCPU).runType != runTypes::counter) {
			double volumeElement = constProd(lightC<double>(), 2 * eps0<double>())
					* (*sCPU).rStep * (*sCPU).tStep * (*sCPU).tStep;
			if((*sCPU).is3D){
				volumeElement *= (*sCPU).rStep;
			} 
			else{
				volumeElement *= vPi<double>();
			}
			for (int64_t i = 0; i < 3 * (*sCPU).Nfreq; i++) {
				(*sCPU).totalSpectrum[i] *= volumeElement;
			}
		}

		return 0;
	}

	static int forwardHankel(ActiveDevice& d, deviceFP* in, deviceComplex* out) {
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		d.deviceLaunch(
			(*sc).Nblock, 
			(*sc).Nthread, 
			hankelKernel{ 
				d.dParamsDevice, 
				in, 
				(deviceFP*)(*sc).workspace1 });
		d.fft((*sc).workspace1, out, deviceFFT::D2Z_1D);
		return 0;
	}
	static int backwardHankel(ActiveDevice& d, deviceComplex* in, deviceFP* out) {
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		d.fft(in, (*sc).workspace1, deviceFFT::Z2D_1D);
		d.deviceLaunch(
			(*sc).Nblock, 
			(*sc).Nthread, 
			inverseHankelKernel{ 
				d.dParamsDevice, 
				(deviceFP*)(*sc).workspace1, 
				out });
		return 0;
	}

	static int addPulseToFieldArrays(
		ActiveDevice& d, 
		pulse<double>& pCPU, 
		const bool useLoadedField, 
		const std::complex<double>* loadedFieldIn) {
		if(pCPU.energy == 0.0) return 0;
		simulationParameterSet* s = d.cParams;
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		deviceParameterSet<deviceFP, deviceComplex>* scDevice = d.dParamsDevice;
		pulse<deviceFP>* p;
		d.deviceCalloc((void**)&p, 1, sizeof(pulse<deviceFP>));
		pulse<deviceFP> devpCPU = pCPU;

		d.deviceMemcpy(
			d.dParamsDevice, 
			sc, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);

		deviceFP* materialPhase;
		deviceComplex* loadedField;

		d.deviceCalloc((void**)&loadedField, (*sc).Ntime, sizeof(deviceComplex));

		//get the material phase
		deviceFP* materialCoefficients, * sellmeierPropagationMedium;

		d.deviceCalloc((void**)&materialCoefficients, 66, sizeof(deviceFP));
		d.deviceCalloc((void**)&sellmeierPropagationMedium, 66, sizeof(deviceFP));
		d.deviceCalloc((void**)&materialPhase, (*s).Nfreq, sizeof(deviceFP));
		d.deviceMemcpy(
			materialCoefficients, 
			(*s).crystalDatabase[pCPU.phaseMaterial].sellmeierCoefficients.data(), 
			66 * sizeof(double), 
			copyType::ToDevice);
		d.deviceMemcpy(
			sellmeierPropagationMedium, 
			(*s).crystalDatabase[(*s).materialIndex].sellmeierCoefficients.data(), 
			66 * sizeof(double), 
			copyType::ToDevice);
		d.deviceLaunch(
			(unsigned int)(*s).Nfreq, 
			1, 
			materialPhaseKernel { 
				(deviceFP)(*s).fStep,
				(*s).Ntime, 
				materialCoefficients, 
				(deviceFP)pCPU.frequency, 
				(deviceFP)pCPU.phaseMaterialThickness,
				(*s).crystalDatabase[pCPU.phaseMaterial].sellmeierType,
				materialPhase });

		deviceFP* pulseSum = &materialCoefficients[0];

		if (useLoadedField) {
			d.deviceMemcpy(
				loadedField, 
				loadedFieldIn, 
				(*s).Nfreq * sizeof(std::complex<double>), 
				copyType::ToDevice);
		}
		d.deviceMemset(pulseSum, 0, sizeof(deviceFP));
		d.deviceMemset((*sc).workspace1, 0, 2 * (*sc).NgridC * sizeof(deviceComplex));
		d.deviceMemcpy(
			p, 
			&devpCPU, 
			sizeof(pulse<deviceFP>), 
			copyType::ToDevice);
		if ((*sc).is3D) {
			d.deviceLaunch(
				(*sc).Nblock / 2, 
				(*sc).Nthread, 
				beamGenerationKernel3D{
					(*sc).workspace1, 
					p, 
					pulseSum, 
					scDevice, 
					useLoadedField, 
					loadedField, 
					materialPhase,
					sellmeierPropagationMedium });
		}
		else {
			d.deviceLaunch(
				(*sc).Nblock / 2, 
				(*sc).Nthread, 
				beamGenerationKernel2D{
					(*sc).workspace1, 
					p, 
					pulseSum, 
					scDevice, 
					useLoadedField, 
					loadedField, 
					materialPhase,
					sellmeierPropagationMedium });
		}

		d.fft((*sc).workspace1, (*sc).gridPolarizationTime1, deviceFFT::Z2D_1D);

		d.deviceLaunch(
			2 * (*sc).Nblock, 
			(*sc).Nthread, 
			beamNormalizeKernel{ 
				scDevice, 
				pulseSum, 
				(*sc).gridPolarizationTime1, 
				(deviceFP)pCPU.energy });

		//add the pulses
		d.deviceLaunch(
			2 * (*sc).Nblock, 
			(*sc).Nthread, 
			addDoubleArraysKernel{ 
				(*sc).gridETime1, 
				(deviceFP*)(*sc).gridPolarizationTime1 });

		//fft onto frequency grid
		d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);

		d.deviceFree(materialPhase);
		d.deviceFree(materialCoefficients);
		d.deviceFree(sellmeierPropagationMedium);
		d.deviceFree(loadedField);
		d.deviceFree(p);
		return 0;
	}

	static int addPreviousGridToFieldArrays(
		ActiveDevice& d, 
		double* loadedField) {
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		
		d.deviceMemcpy(
			(deviceFP*)(*sc).gridPolarizationTime1, 
			loadedField, 
			(*sc).Ngrid * 2 * sizeof(double), 
			copyType::ToDevice);

		//add the pulses
		d.deviceLaunch(
			2 * (*sc).Nblock, 
			(*sc).Nthread, 
			addDoubleArraysKernel{ 
				(*sc).gridETime1, 
				(deviceFP*)(*sc).gridPolarizationTime1 });

		//fft onto frequency grid
		d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);

		return 0;
	}
	
	static int prepareElectricFieldArrays(ActiveDevice& d) {

		simulationParameterSet* s = d.cParams;
		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		
		d.deviceMemcpy(
			d.dParamsDevice, 
			sc, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);
		deviceParameterSet<deviceFP, deviceComplex>* scDevice = d.dParamsDevice;
		
		if (!(*s).isFollowerInSequence || (*s).isReinjecting) {
			if (!(*s).isReinjecting) {
				d.deviceMemset((*sc).gridETime1, 0, 2 * (*sc).Ngrid * sizeof(deviceFP));
			}
			else {
				d.deviceMemcpy(
					(*sc).gridETime1, 
					(*s).ExtOut, 
					2 * (*s).Ngrid * sizeof(double), 
					copyType::ToDevice);
				d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);
			}
			
			if (d.cParams->pulse1FileType == 3) {
				addPreviousGridToFieldArrays(d, d.cParams->loadedFullGrid1);
			}
			else {
				addPulseToFieldArrays(
					d, 
					d.cParams->pulse1, 
					d.cParams->field1IsAllocated, 
					d.cParams->loadedField1);
			}
			
			if (d.cParams->pulse2FileType == 3) {
				addPreviousGridToFieldArrays(d, d.cParams->loadedFullGrid2);
			}
			else {
				addPulseToFieldArrays(
					d, 
					d.cParams->pulse2, 
					d.cParams->field2IsAllocated, 
					d.cParams->loadedField2);
			}
			
		}
		else {
			d.deviceMemcpy(
				(*sc).gridETime1, 
				(*s).ExtOut, 
				2 * (*s).Ngrid * sizeof(double), 
				copyType::ToDevice);
			d.fft((*sc).gridETime1, (*sc).gridEFrequency1, deviceFFT::D2Z);
		}
		
		if((*sc).axesNumber==2){
			d.deviceLaunch(
				(*sc).Nblock / 2,
				(*sc).Nthread,
				biaxialRotationKernel {d.dParamsDevice,(*sc).gridEFrequency1,true}
			);
		}

		//Copy the field into the temporary array
		d.deviceMemcpy(
			(void*)(*sc).gridEFrequency1Next1, 
			(void*)(*sc).gridEFrequency1, 
			2 * (*sc).NgridC * sizeof(deviceComplex), 
			copyType::OnDevice);


		

		//set the propagation grids how they should be at the beginning of the next step
		d.deviceLaunch(
			(unsigned int)((*sc).NgridC / minGridDimension), 
			2 * minGridDimension, 
			multiplicationKernelCompactDoubleVector{
				(*sc).fieldFactor1, 
				(*sc).gridEFrequency1Next1, 
				(*sc).workspace1, 
				scDevice });

		return 0;
	}

	static int applyFresnelLoss(
		ActiveDevice& d, 
		simulationParameterSet* s, 
		deviceParameterSet<deviceFP, 
		deviceComplex>& sc, 
		int materialIndex1, 
		int materialIndex2) {
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(
			sellmeierCoefficientsAugmentedCPU, 
			(*s).crystalDatabase[materialIndex1].sellmeierCoefficients.data(), 
			66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		deviceFP* sellmeierCoefficients1;
		deviceFP* sellmeierCoefficients2;
		d.deviceCalloc((void**)&sellmeierCoefficients1, 74, sizeof(deviceFP));
		d.deviceCalloc((void**)&sellmeierCoefficients2, 74, sizeof(deviceFP));
		d.deviceMemcpy(
			sellmeierCoefficients1, 
			sellmeierCoefficientsAugmentedCPU, 
			(66 + 8) * sizeof(double), 
			copyType::ToDevice);

		memcpy(
			sellmeierCoefficientsAugmentedCPU, 
			(*s).crystalDatabase[materialIndex2].sellmeierCoefficients.data(), 
			66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		d.deviceMemcpy(
			sellmeierCoefficients2, 
			sellmeierCoefficientsAugmentedCPU, 
			(66 + 8) * sizeof(double), 
			copyType::ToDevice);
		d.deviceMemcpy(
			sc.gridEFrequency1, 
			(*s).EkwOut, 
			2 * (*s).NgridC * sizeof(std::complex<double>), 
			copyType::ToDevice);

		//transform final result
		d.fft(sc.gridEFrequency1, sc.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(2 * sc.Nblock, sc.Nthread, multiplyByConstantKernelD{
			sc.gridETime1, static_cast<deviceFP>(1.0 / sc.Ngrid) });
		//copy the field arrays from the GPU to CPU memory
		d.deviceMemcpy(
			(*s).ExtOut, 
			sc.gridETime1, 
			2 * (*s).Ngrid * sizeof(double), 
			copyType::ToHost);
		d.deviceMemcpy(
			(*s).EkwOut, 
			sc.gridEFrequency1, 
			2 * (*s).Ngrid * sizeof(std::complex<double>), 
			copyType::ToHost);

		d.deviceFree(sellmeierCoefficients1);
		d.deviceFree(sellmeierCoefficients2);

		return 0;
	}

	static int applyFilter(ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		const double f0, 
		const double bandwidth, 
		const double order, 
		const double inBandAmplitude, 
		const double outOfBandAmplitude) {

		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			filterKernel{
				sDevice,
				static_cast<deviceFP>(1.0e12 * f0),
				static_cast<deviceFP>(1.0e12 * bandwidth),
				static_cast<int>(round(order)),
				static_cast<deviceFP>(inBandAmplitude),
				static_cast<deviceFP>(outOfBandAmplitude) });

		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);

		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);

		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		getTotalSpectrum(d);

		return 0;
	}

	static int applyOptic(
		ActiveDevice& d,
		simulationParameterSet& sCPU,
		const int index,
		const bool applyX,
		const bool applyY){
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			sCPU.ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		std::vector<std::complex<deviceFP>> complexReflectivity = 
			sCPU.optics[index].toComplexSpectrum<deviceFP>(sCPU.Nfreq,sCPU.fStep);
		d.deviceMemcpy(
			d.deviceStruct.gridPropagationFactor1, 
			complexReflectivity.data(), 
			sCPU.Nfreq*sizeof(deviceComplex), 
			copyType::ToDevice
		);
		d.deviceLaunch(
			d.deviceStruct.Nblock/2,
			d.deviceStruct.Nthread,
			applyOpticKernel{d.dParamsDevice,d.deviceStruct.gridPropagationFactor1, applyX, applyY}
		);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			sCPU.EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		d.deviceMemcpy(
			sCPU.ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * sCPU.Ngrid * sizeof(double), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyLorenzian(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		const double amplitude, 
		const double f0, 
		const double gamma, 
		const double radius, 
		const double order) {

		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, lorentzianSpotKernel{
			sDevice,
			(deviceFP)amplitude,
			(deviceFP)(1.0e12 * f0),
			(deviceFP)(1.0e12 * gamma),
			(deviceFP)radius,
			(deviceFP)order });
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ntime) });
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		getTotalSpectrum(d);

		return 0;
	}

	static int applyAperatureFarFieldHankel(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		double diameter, 
		double activationParameter, 
		double xOffset, 
		double yOffset) {
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		forwardHankel(d, d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			apertureFarFieldKernelHankel{
				sDevice,
				(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter),
				(deviceFP)activationParameter,
				(deviceFP)(deg2Rad<deviceFP>() * xOffset),
				(deviceFP)(deg2Rad<deviceFP>() * yOffset) });
		backwardHankel(d, d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyInverseAperatureFarFieldHankel(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		double diameter, 
		double activationParameter, 
		double xOffset, 
		double yOffset) {
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		forwardHankel(d, d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			inverseApertureFarFieldKernelHankel{
				sDevice,
				(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter),
				(deviceFP)activationParameter,
				(deviceFP)(deg2Rad<deviceFP>() * xOffset),
				(deviceFP)(deg2Rad<deviceFP>() * yOffset) });
		backwardHankel(d, d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyAperatureFarField(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		double diameter, 
		double activationParameter, 
		double xOffset, 
		double yOffset) {
		if ((*sCPU).isCylindric) {
			return applyAperatureFarFieldHankel(d, sCPU, diameter, activationParameter, xOffset, yOffset);
		}
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, apertureFarFieldKernel{
			sDevice,
			(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter),
			(deviceFP)activationParameter,
			(deviceFP)(deg2Rad<deviceFP>() * xOffset),
			(deviceFP)(deg2Rad<deviceFP>() * yOffset) });

		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);

		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);

		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		getTotalSpectrum(d);
		return 0;
	}

	static int applyInverseAperatureFarField(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		double diameter, 
		double activationParameter, 
		double xOffset, 
		double yOffset) {
		if ((*sCPU).isCylindric) {
			return applyInverseAperatureFarFieldHankel(d, sCPU, diameter, activationParameter, xOffset, yOffset);
		}
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(d.deviceStruct.Nblock / 2, d.deviceStruct.Nthread, inverseApertureFarFieldKernel{
			sDevice,
			(deviceFP)(0.5 * deg2Rad<deviceFP>() * diameter),
			(deviceFP)activationParameter,
			(deviceFP)(deg2Rad<deviceFP>() * xOffset),
			(deviceFP)(deg2Rad<deviceFP>() * yOffset) });

		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);

		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);

		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		getTotalSpectrum(d);
		return 0;
	}

	static int applyAperature(
		ActiveDevice& d, 
		const simulationParameterSet* sCPU, 
		const double diameter, 
		const double activationParameter) {
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);

		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceLaunch(
			d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			apertureKernel{
				sDevice,
				(deviceFP)(0.5 * diameter),
				(deviceFP)(activationParameter) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applySphericalMirror(
		ActiveDevice& d, 
		const simulationParameterSet* sCPU, 
		deviceParameterSet<deviceFP, deviceComplex>& s, 
		const double ROC) {

		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceMemcpy(
			sDevice, 
			&s, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);

		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			sphericalMirrorKernel{ sDevice, (deviceFP)ROC });
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.deviceLaunch(
			2 * d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			multiplyByConstantKernelD{ 
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ntime )});
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyParabolicMirror(ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		deviceParameterSet<deviceFP, deviceComplex>& s, 
		const double focus) {

		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;
		d.deviceMemcpy(
			d.deviceStruct.gridETime1, 
			(*sCPU).ExtOut, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToDevice);
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z_1D);
		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			parabolicMirrorKernel{ sDevice, (deviceFP)focus });
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D_1D);
		d.deviceLaunch(
			2 * d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			multiplyByConstantKernelD { 
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ntime) });
		d.fft(d.deviceStruct.gridETime1, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		getTotalSpectrum(d);
		return 0;
	}

	static int applyLinearPropagation(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		const int materialIndex, 
		const double thickness) {

		if (d.hasPlasma) {
			simulationParameterSet sCopy = *sCPU;
			sCopy.nonlinearAbsorptionStrength = 0.0;
			d.reset(&sCopy);
		}

		d.deviceMemcpy(
			d.deviceStruct.gridEFrequency1, 
			(*sCPU).EkwOut, 
			d.deviceStruct.NgridC * 2 * sizeof(std::complex<double>), 
			copyType::ToDevice);

		deviceFP* sellmeierCoefficients = (deviceFP*)d.deviceStruct.gridEFrequency1Next1;
		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[74] = { 0 };
		memcpy(
			sellmeierCoefficientsAugmentedCPU, 
			(*sCPU).crystalDatabase[materialIndex].sellmeierCoefficients.data(), 
			66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*sCPU).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*sCPU).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*sCPU).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*sCPU).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*sCPU).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*sCPU).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		d.deviceMemcpy(
			sellmeierCoefficients, 
			sellmeierCoefficientsAugmentedCPU, 
			(66 + 8) * sizeof(double), 
			copyType::ToDevice);
		d.deviceStruct.axesNumber = (*sCPU).crystalDatabase[materialIndex].axisType;
		d.deviceStruct.sellmeierType = (*sCPU).crystalDatabase[materialIndex].sellmeierType;
		deviceParameterSet<deviceFP, deviceComplex>* sDevice = d.dParamsDevice;

		d.deviceLaunch(
			d.deviceStruct.Nblock / 2, 
			d.deviceStruct.Nthread, 
			applyLinearPropagationKernel{ 
				sellmeierCoefficients, 
				(deviceFP)thickness, 
				sDevice });
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			d.deviceStruct.NgridC * 2 * sizeof(std::complex<double>), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(
			2 * d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			multiplyByConstantKernelD{ 
				d.deviceStruct.gridETime1, 
				(deviceFP)(1.0 / d.deviceStruct.Ngrid) });

		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);
		getTotalSpectrum(d);

		return 0;
	}

	static int preparePropagationGrids(
		ActiveDevice& d) {

		deviceParameterSet<deviceFP, deviceComplex>* sc = d.s;
		simulationParameterSet* s = d.cParams;
		deviceFP* sellmeierCoefficients = (deviceFP*)(*sc).gridEFrequency1Next1;
		//construct augmented sellmeier coefficients used in the kernel to find the walkoff angles
		double sellmeierCoefficientsAugmentedCPU[79];
		memcpy(
			sellmeierCoefficientsAugmentedCPU, 
			(*s).sellmeierCoefficients, 
			66 * (sizeof(double)));
		sellmeierCoefficientsAugmentedCPU[66] = (*s).crystalTheta;
		sellmeierCoefficientsAugmentedCPU[67] = (*s).crystalPhi;
		sellmeierCoefficientsAugmentedCPU[68] = (*s).axesNumber;
		sellmeierCoefficientsAugmentedCPU[69] = (*s).sellmeierType;
		sellmeierCoefficientsAugmentedCPU[70] = (*s).kStep;
		sellmeierCoefficientsAugmentedCPU[71] = (*s).fStep;
		sellmeierCoefficientsAugmentedCPU[72] = 1.0e-12;
		memcpy(
			sellmeierCoefficientsAugmentedCPU + 72, 
			(*s).crystalDatabase[(*s).materialIndex].nonlinearReferenceFrequencies.data(), 
			7 * sizeof(double));
		d.deviceMemcpy(
			sellmeierCoefficients, 
			sellmeierCoefficientsAugmentedCPU, 
			79 * sizeof(double), 
			copyType::ToDevice);

		//prepare the propagation grids
		deviceParameterSet<deviceFP, deviceComplex>* sD = d.dParamsDevice;
		d.deviceMemcpy(
			sD, 
			sc, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToDevice);
		d.deviceLaunch(
			(unsigned int)maxN((*sc).Ntime / 2, 82), 
			1, 
			getChiLinearKernel{ sD, sellmeierCoefficients });
		if ((*s).is3D) {
			d.deviceLaunch(
				static_cast<unsigned int>((*sc).Nblock / 2u), 
				static_cast<unsigned int>((*sc).Nthread), 
				prepare3DGridsKernel{ sellmeierCoefficients, sD });
		}
		else if ((*s).isCylindric) {
			d.deviceLaunch(
				static_cast<unsigned int>((*sc).Nblock / 2u), 
				static_cast<unsigned int>((*sc).Nthread), 
				prepareCylindricGridsKernel{ sellmeierCoefficients, sD });
		}
		else {
			d.deviceLaunch(
				static_cast<unsigned int>((*sc).Nblock / 2u), 
				static_cast<unsigned int>((*sc).Nthread), 
				prepareCartesianGridsKernel{ sellmeierCoefficients, sD });
		}
		d.deviceMemcpy(
			sc, 
			sD, 
			sizeof(deviceParameterSet<deviceFP, deviceComplex>), 
			copyType::ToHost);
		return 0;
	}

	//Rotate the field on the GPU
	//Allocates memory and copies from CPU, then copies back to CPU and deallocates
	// - inefficient but the general principle is that only the CPU memory is preserved
	// after simulations finish... and this only runs at the end of the simulation
	static int rotateField(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		const double rotationAngle) {

		deviceComplex* Ein1 = d.deviceStruct.gridEFrequency1;
		deviceComplex* Ein2 = d.deviceStruct.gridEFrequency2;
		deviceComplex* Eout1 = d.deviceStruct.gridEFrequency1Next1;
		deviceComplex* Eout2 = d.deviceStruct.gridEFrequency1Next2;

		//retrieve/rotate the field from the CPU memory
		d.deviceMemcpy(
			Ein1, 
			(*sCPU).EkwOut, 
			2 * (*sCPU).NgridC * sizeof(std::complex<double>), 
			copyType::ToDevice);

		if(rotationAngle == deg2Rad<deviceFP>()*90.0){
			d.deviceLaunch(
				static_cast<unsigned int>(d.deviceStruct.NgridC / minGridDimension), 
				minGridDimension, 
				rotateField90Kernel{ Ein1, Ein2, Eout1, Eout2});
		}
		else if(rotationAngle == deg2Rad<deviceFP>()*180.0){
			d.deviceLaunch(
				static_cast<unsigned int>(d.deviceStruct.NgridC / minGridDimension), 
				minGridDimension, 
				rotateField180Kernel{ Ein1, Ein2, Eout1, Eout2});
		}
		else {
			d.deviceLaunch(
				static_cast<unsigned int>(d.deviceStruct.NgridC / minGridDimension), 
				minGridDimension, 
				rotateFieldKernel{ Ein1, Ein2, Eout1, Eout2, (deviceFP)rotationAngle });
		}
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			Eout1, 
			2 * (*sCPU).NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);

		//transform back to time
		d.fft(Eout1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(
			2 * d.deviceStruct.Nblock, 
			d.deviceStruct.Nthread, 
			multiplyByConstantKernelD{ 
				d.deviceStruct.gridETime1, 
				static_cast<deviceFP>(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);

		//update spectrum
		getTotalSpectrum(d);
		return 0;
	}

	static void rotateBiaxial(
		ActiveDevice& d,
		simulationParameterSet* sCPU,
		int materialIndex,
		double theta,
		double phi,
		double frequency,
		bool forward
		){
		double ls = 1e6*lightC<double>()/frequency;
		ls *= ls;
		double omega = twoPi<double>()*frequency;
		if((*sCPU).crystalDatabase[materialIndex].axisType != 2) return;
		std::complex<double> na = hostSellmeierFunc(ls, omega, &((*sCPU).crystalDatabase[materialIndex].sellmeierCoefficients[0]), (*sCPU).crystalDatabase[materialIndex].sellmeierType);
		std::complex<double> nb = hostSellmeierFunc(ls, omega, &((*sCPU).crystalDatabase[materialIndex].sellmeierCoefficients[22]), (*sCPU).crystalDatabase[materialIndex].sellmeierType);
		std::complex<double> nc = hostSellmeierFunc(ls, omega, &((*sCPU).crystalDatabase[materialIndex].sellmeierCoefficients[44]), (*sCPU).crystalDatabase[materialIndex].sellmeierType);
		na *= na;
		nb *= nb;
		nc *= nc;
		double cp = std::cos(phi);
		double sp = std::sin(phi);
		double ct = std::cos(theta);
		double st = std::sin(theta);

		double delta = (na == nb) ? 
		double{} : 
		0.5 * std::atan(std::sin(2*phi) * ct / 
		(((1.0/nb.real() - 1.0/nc.real())/(1.0/na.real() - 1.0/nb.real())) * st*st - cp*cp * ct*ct + sp*sp));
		if(!forward) delta *= -1;
		rotateField(d, sCPU, delta);
	}

	static void savePlasma(
		ActiveDevice& d, 
		int64_t saveloc, 
		int densityLoc){

		const deviceParameterSet<deviceFP, deviceComplex>* sH = d.s; 
		const deviceParameterSet<deviceFP, deviceComplex>* sD = d.dParamsDevice;
		//prepare the grids
		preparePropagationGrids(d);
		prepareElectricFieldArrays(d);

		//Clear out the memory to which the data will be saved
		d.deviceMemset(d.deviceStruct.gridEFrequency1,0,4*sizeof(deviceFP)*d.cParams->NgridC);

		//run the plasma kernels
		d.deviceLaunch(
			(*sH).Nblock, 
			(*sH).Nthread, 
			plasmaCurrentKernel_twoStage_A{ sD });
		d.deviceLaunch(
			(unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 
			minGridDimension, 
			plasmaCurrentKernel_SaveOutput{ sD });

		//save to memory
		d.deviceMemcpy(
			d.cParams->ExtOut + 2*d.cParams->Ngrid*saveloc, 
			d.deviceStruct.gridPolarizationTime1, 
			2 * d.deviceStruct.Ngrid * sizeof(double), 
			copyType::ToHost);

		if (densityLoc == 1) {
			d.deviceMemcpy(
				d.cParams->ExtOut + 2 * d.cParams->Ngrid * saveloc, 
				reinterpret_cast<deviceFP*>(d.deviceStruct.gridEFrequency1), 
				d.deviceStruct.Ngrid * sizeof(double), 
				copyType::ToHost);
		}
		else if (densityLoc == 2) {
			d.deviceMemcpy(
				d.cParams->ExtOut + 2 * d.cParams->Ngrid * saveloc + d.cParams->Ngrid, 
				reinterpret_cast<deviceFP*>(d.deviceStruct.gridEFrequency1), 
				d.deviceStruct.Ngrid * sizeof(double), 
				copyType::ToHost);
		}
		
	}
//function to run a RK4 time step
//stepNumber is the sub-step index, from 0 to 3
	static int runRK4Step(
		ActiveDevice& d, 
		const uint8_t stepNumber){

		const deviceParameterSet<deviceFP, deviceComplex>* sH = d.s; 
		const deviceParameterSet<deviceFP, deviceComplex>* sD = d.dParamsDevice;

		// Beam with symmetry around z axis:
		// Nonlinear polarization and plasma use expanded grid
		// Radial laplacian uses standard grid
		// two possible FFT shapes (radial Laplacian always performed)
		if((*sH).isCylindric){
			//ifft to time domain 
			d.fft((*sH).workspace1, (*sH).gridETime1, deviceFFT::Z2D);

			//Nonlinear polarization and plasma current are fft-ed in a batch
			//from separate (de-interlaced) time-domain grids.
			//assumption: no plasma without other nonlinearities
			if((*sH).isNonLinear){
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					nonlinearPolarizationKernel{ sD });
				if((*sH).hasPlasma){
					d.deviceLaunch(
						(*sH).Nblock, 
						(*sH).Nthread, 
						plasmaCurrentKernel_twoStage_A{ sD });
					
					//CUDA and other platforms perform very 
					//differently for different versions of this kernel, use optimum per platform
					#ifdef __CUDACC__
					d.deviceLaunch(
						(unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 
						2 * minGridDimension, 
						plasmaCurrentKernel_twoStage_B{ sD });
					#else
					d.deviceLaunch(
						(unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 
						minGridDimension, 
						plasmaCurrentKernel_twoStage_B_simultaneous{ sD });
					#endif
					
					d.deviceLaunch(
						(*sH).Nblock, 
						(*sH).Nthread, 
						expandCylindricalBeam{ sD });
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFT::D2Z_Polarization);
					d.deviceLaunch(
						(*sH).Nblock / 2, 
						(*sH).Nthread, 
						updateKwithPlasmaKernelCylindric{ sD });
				}
				else{
					d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFT::D2Z_Polarization);
					d.deviceLaunch(
						(*sH).Nblock / 2, 
						(*sH).Nthread, 
						updateKwithPolarizationKernelCylindric{ sD });
				}
			}
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				radialLaplacianKernel{ sD });
			d.fft((*sH).gridRadialLaplacian1, (*sH).workspace1, deviceFFT::D2Z);

			switch (stepNumber) {
			case 0:
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					rkKernel0Cylindric{ sD });
				break;
			case 1:
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					rkKernel1Cylindric{ sD });
				break;
			case 2:
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					rkKernel2Cylindric{ sD });
				break;
			case 3:
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					rkKernel3Cylindric{ sD });
				break;
			}
			return 0;
		}
		// 2D and 3D cartesian
		// Only one type of FFT
		// currently nonlinear polarization and plasma ffts are not batched, could give
		// minor speed boost by combining them, but requires additional memory, so
		// nonlinear polarization and plasma are fft-ed separately to accommodate larger
		// systems.
		else if ((*sH).isNonLinear) {
			//perform inverse FFT to get time-space electric field
			if((*sH).axesNumber==2){
				//undo delta rotation if biaxial
				d.deviceLaunch(
					(*sH).Nblock / 2,
					(*sH).Nthread,
					biaxialRotationKernel {sD,(*sH).workspace1,false}
				);
			}
			d.fft((*sH).workspace1, (*sH).gridETime1, deviceFFT::Z2D);
			//Plasma/multiphoton absorption
			if ((*sH).hasPlasma) {
				d.deviceLaunch(
					(*sH).Nblock, 
					(*sH).Nthread, 
					plasmaCurrentKernel_twoStage_A{ sD });
				d.deviceLaunch(
					(unsigned int)(((*sH).Nspace2 * (*sH).Nspace) / minGridDimension), 
					minGridDimension, 
					plasmaCurrentKernel_twoStage_B_simultaneous{ sD });
				d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFT::D2Z);
				d.deviceLaunch(
					(*sH).Nblock / 2, 
					(*sH).Nthread, 
					updateKwithPlasmaKernel{ sD });
			}
			//Nonlinear polarization
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				nonlinearPolarizationKernel{ sD });
			d.fft((*sH).gridPolarizationTime1, (*sH).workspace1, deviceFFT::D2Z);
			//rotate nonlinear polarization into the delta frame if biaxial
			if((*sH).axesNumber==2){
				d.deviceLaunch(
					(*sH).Nblock / 2,
					(*sH).Nthread,
					biaxialRotationKernel {sD,(*sH).workspace1,true}
				);
			}
		}

		//advance an RK4 step
		switch (stepNumber) {
		case 0:
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				rkKernel0{ sD });
			break;
		case 1:
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				rkKernel1{ sD });
			break;
		case 2:
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				rkKernel2{ sD });
			break;
		case 3:
			d.deviceLaunch(
				(*sH).Nblock, 
				(*sH).Nthread, 
				rkKernel3{ sD });
			break;
		}
		return 0;
	}

	static unsigned long int solveNonlinearWaveEquationWithDevice(
		ActiveDevice& d, 
		simulationParameterSet* sCPU) {

		if (sCPU->isFDTD) {
			return solveFDTD(d, sCPU, 5, 0, 1e-6, 1e-6);
		}
		//prepare the propagation arrays
		preparePropagationGrids(d);
		prepareElectricFieldArrays(d);

		deviceFP* canaryPointer = 
			&d.deviceStruct.gridETime1[
				d.deviceStruct.Ntime / 2 
					+ d.deviceStruct.Ntime 
					* (d.deviceStruct.Nspace / 2 
						+ d.deviceStruct.Nspace 
						* (d.deviceStruct.Nspace2 / 2))];
		//Core propagation loop
		for (int64_t i = 0; i < d.deviceStruct.Nsteps; ++i) {

			//RK4
			runRK4Step(d, 0);
			runRK4Step(d, 1);
			runRK4Step(d, 2);
			runRK4Step(d, 3);

			//periodically check if the simulation diverged or was cancelled
			if ((*sCPU).cancellationCalled) break;
			if (i % 10 == 0 && d.isTheCanaryPixelNaN(canaryPointer)) break;
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
		}
		if ((*sCPU).isInFittingMode && !(*sCPU).isInSequence)(*(*sCPU).progressCounter)++;

		if(d.deviceStruct.axesNumber==2){
			d.deviceLaunch(
				d.s->Nblock / 2,
				d.s->Nthread,
				biaxialRotationKernel {d.dParamsDevice,d.deviceStruct.gridEFrequency1,false}
			);
		}
		//take final spectra and transfer the results to the CPU
		d.deviceMemcpy(
			(*sCPU).EkwOut, 
			d.deviceStruct.gridEFrequency1, 
			2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
			copyType::ToHost);
		d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
		d.deviceLaunch(
			(int)(d.deviceStruct.Ngrid / minGridDimension), 
			2 * minGridDimension, 
			multiplyByConstantKernelD{ 
				d.deviceStruct.gridETime1, 
				static_cast<deviceFP>(1.0 / d.deviceStruct.Ngrid) });
		d.deviceMemcpy(
			(*sCPU).ExtOut, 
			d.deviceStruct.gridETime1, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToHost);
		getTotalSpectrum(d);

		return 13 * d.isTheCanaryPixelNaN(canaryPointer);
	}

	template<typename maxwellType>
	static unsigned int freeFDTD(
		ActiveDevice& d, 
		maxwellType& maxCalc) {
		unsigned int errorValue = 13 * d.isTheCanaryPixelNaN(&(maxCalc.Egrid[0].y));
		d.deviceFree(maxCalc.Egrid);
		d.deviceFree(maxCalc.EgridEstimate);
		d.deviceFree(maxCalc.EgridEstimate2);
		d.deviceFree(maxCalc.EgridNext);
		d.deviceFree(maxCalc.Hgrid);
		d.deviceFree(maxCalc.HgridEstimate);
		d.deviceFree(maxCalc.HgridEstimate2);
		d.deviceFree(maxCalc.HgridNext);
		d.deviceFree(maxCalc.materialGrid);
		d.deviceFree(maxCalc.materialGridNext);
		d.deviceFree(maxCalc.materialGridEstimate);
		d.deviceFree(maxCalc.materialGridEstimate2);
		d.deviceFree(maxCalc.deviceCopy);
		if(maxCalc.hasMaterialMap){
			d.deviceFree(maxCalc.materialMap);
			d.deviceFree(maxCalc.oscillatorIndexMap);
		}
		return errorValue;
	}


	template<typename maxwellType>
	static void calculateFDTDParameters(
		const simulationParameterSet* sCPU, 
		maxwellType& maxCalc,
		int mapValue = 0){

		if(mapValue == 0){
			double n0 = hostSellmeierFunc(
			0, twoPi<double>() * sCPU->pulse1.frequency, 
			(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			double nm1 = hostSellmeierFunc(
				0, twoPi<double>() * (-2e11 + sCPU->pulse1.frequency), 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			double np1 = hostSellmeierFunc(
				0, twoPi<double>() * (2e11 + sCPU->pulse1.frequency), 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			double nGroup = n0 + sCPU->pulse1.frequency * (np1 - nm1) / 4.0e11;

			if(maxCalc.waitFrames == -1){
				maxCalc.waitFrames = 
					(maxCalc.frontBuffer + nGroup * maxCalc.crystalThickness + 10 * maxCalc.zStep) 
					/ (lightC<double>() * maxCalc.tStep);
				maxCalc.waitFrames = maxCalc.tGridFactor * (maxCalc.waitFrames / maxCalc.tGridFactor);
			}

			maxCalc.Nt = maxCalc.waitFrames + (*sCPU).Ntime * maxCalc.tGridFactor;
		}
		

		//copy the crystal info
		if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].axisType == 0) {
			//isotropic
			for (int i = 0; i < 22; i++) {
				maxCalc.sellmeierEquations[i][mapValue].x = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
				maxCalc.sellmeierEquations[i][mapValue].y = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
				maxCalc.sellmeierEquations[i][mapValue].z = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
			}
		}
		else if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].axisType == 1) {
			//uniaxial
			for (int i = 0; i < 22; i++) {
				maxCalc.sellmeierEquations[i][mapValue].x = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
				maxCalc.sellmeierEquations[i][mapValue].y = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]);
				maxCalc.sellmeierEquations[i][mapValue].z = 
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i+22]);
			}
		}
		else {
			//biaxial
			for (int i = 0; i < 22; i++) {
				maxCalc.sellmeierEquations[i][mapValue] = maxwellPoint<deviceFP>{
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i]),
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i+22]),
					static_cast<deviceFP>(
						(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients[i+44])
				};
			}
		}
		
		if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearSwitches.hasChi2) 
			maxCalc.hasChi2[mapValue] = true;
		if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearSwitches.hasChi3) {
			if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearSwitches.assumeCentrosymmetric) {
				maxCalc.hasSingleChi3[mapValue] = true;
			}
			else{
				maxCalc.hasFullChi3[mapValue] = true;
			}
		}
			

		//perform millers rule normalization on the nonlinear coefficients while copying the values
		double millersRuleFactorChi2 = 2e-12;
		double millersRuleFactorChi3 = 1.0;
		if ((*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[0]) {
			double wRef = 
				twoPi<double>() * 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[0];
			double nRef = hostSellmeierFunc(
				0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi2 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[1];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi2 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[2];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi2 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[3];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi3 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[4];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi3 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[5];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi3 /= nRef * nRef - 1.0;
			wRef = twoPi<double>() 
				* (*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].nonlinearReferenceFrequencies[6];
			nRef = hostSellmeierFunc(0, wRef, 
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].sellmeierCoefficients.data(), 1).real();
			millersRuleFactorChi3 /= nRef * nRef - 1.0;
		}
		for (int i = 0; i < 6; i++) {
			maxCalc.chi2[i][mapValue] = maxwellPoint<deviceFP>{
				static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].d[3 * i] * millersRuleFactorChi2),
				static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].d[3 * i + 1] * millersRuleFactorChi2),
				static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].d[3 * i + 2] * millersRuleFactorChi2) };
			if (i > 2) maxCalc.chi2[i][mapValue] *= 2.0;
		}
		for (int i = 0; i < 27; i++) {
			maxCalc.chi3[i][mapValue].x = static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].chi3[i] * millersRuleFactorChi3);
			maxCalc.chi3[i][mapValue].y = static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].chi3[i+27] * millersRuleFactorChi3);
			maxCalc.chi3[i][mapValue].z = static_cast<deviceFP>(
				(*sCPU).crystalDatabase[maxCalc.materialKeys[mapValue]].chi3[i+54] * millersRuleFactorChi3);
		}

		//count the nonzero oscillators
		if(!maxCalc.hasMaterialMap){
			for (int i = 0; i < 7; i++) {
				if (maxCalc.sellmeierEquations[1 + i * 3][mapValue].x > 0.0) maxCalc.Noscillators++;
				else break;
			}

			//collect plasma properties
			if ((*sCPU).nonlinearAbsorptionStrength > 0.0 || (*sCPU).startingCarrierDensity > 0.0) {
				maxCalc.Noscillators++;
				maxCalc.hasPlasma[mapValue] = true;
				maxCalc.kCarrierGeneration[mapValue] = 2.0 / ((*sCPU).bandGapElectronVolts);
				maxCalc.kDrude[mapValue] = - elCharge<double>() 
					/ ((*sCPU).effectiveMass * elMass<double>());
				maxCalc.gammaDrude[mapValue] = (*sCPU).drudeGamma;
				maxCalc.kNonlinearAbsorption[mapValue] =  0.5 * (*sCPU).nonlinearAbsorptionStrength;
				maxCalc.startingCarriers[mapValue] = (*sCPU).startingCarrierDensity;
				maxCalc.nonlinearAbsorptionOrder[mapValue] = static_cast<int>(
					std::ceil(eVtoHz<double>() * (*sCPU).bandGapElectronVolts 
						/ (*sCPU).pulse1.frequency)) - 1;
			}
			maxCalc.NMaterialGrid =
				(maxCalc.materialStop - maxCalc.materialStart)
				* maxCalc.Nx * maxCalc.Ny * maxCalc.Noscillators;
		}
	}

	static void prepareFDTD(
		ActiveDevice& d, 
		const simulationParameterSet* sCPU, 
		maxwell3D& maxCalc,
		const int gridIndex = -1) {

		//Check if there is a material map and allocate/load it if necessary
		if (gridIndex != -1 && (*sCPU).runType != runTypes::counter && (*sCPU).runType != runTypes::cluster) {
			//throw std::runtime_error(std::string("path string:\n").append(materialMapPath));
			maxCalc.hasMaterialMap = false;
			int64_t NmaterialPoints = 0;
			std::vector<int8_t> materialMapCPU(maxCalc.Ngrid, 0);
			std::vector<int64_t> oscillatorIndexMapCPU(maxCalc.Ngrid, 0);
			std::stringstream fs(sCPU->optics[gridIndex].fileContents);

			if (fs.good()) {
				auto moveToColon = [&]() {
					char x = 0;
					while (x != ':' && fs.good()) {
						fs >> x;
					}
					return 0;
					};
				//First line: materials
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.materialKeys[i];
					//optimization: count oscillators
				}
				maxCalc.Noscillators = 7;
				//Second line: theta
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.materialTheta[i];
				}
				//Third line: phi
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.materialPhi[i];
				}
				//Fourth: Nonlinear absorption (NOTE: FIX PLASMA PARAMS)
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.kNonlinearAbsorption[i];
					maxCalc.hasPlasma[i] = (maxCalc.kNonlinearAbsorption[i] > 0.0) || (maxCalc.startingCarriers[i] > 0);
				}
				//Fifth: bandgap
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.kDrude[i];
				}
				//Sixth: gamma
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.gammaDrude[i];
				}
				//Seventh: effective mass
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.kCarrierGeneration[i];
				}
				//Eighth: Initial carriers
				moveToColon();
				for (int i = 0; i < NmaterialMax; i++) {
					fs >> maxCalc.startingCarriers[i];
				}
				//9th: Dimensions
				int64_t fileNx{};
				int64_t fileNy{};
				int64_t fileNz{};
				moveToColon();
				fs >> fileNx >> fileNy >> fileNz;
				if(fileNx != maxCalc.Nx) throw std::runtime_error("Nx of grid doesn't match file: "+std::to_string(fileNx)+" vs. "+std::to_string(maxCalc.Nx));
				if(fileNy != maxCalc.Ny) throw std::runtime_error("Ny of grid doesn't match file: "+std::to_string(fileNy)+" vs. "+std::to_string(maxCalc.Ny));
				if(fileNz != maxCalc.Nz) throw std::runtime_error("Nz of grid doesn't match file: "+std::to_string(fileNz)+" vs. "+std::to_string(maxCalc.Nz));
				//Eighth: start data
				int64_t gridCount{};
				while (fs.good() && gridCount < maxCalc.Ngrid) {
					int currentInt;
					fs >> currentInt;
					materialMapCPU[gridCount] = static_cast<int8_t>(currentInt);
					if (materialMapCPU[gridCount] > 0) {
						oscillatorIndexMapCPU[gridCount] = NmaterialPoints;
						NmaterialPoints++;
					}
					gridCount++;
				}

				d.deviceCalloc((void**)&(maxCalc.materialMap),
					maxCalc.Ngrid, sizeof(char));
				d.deviceMemcpy(
					(void*)maxCalc.materialMap,
					(void*)materialMapCPU.data(),
					maxCalc.Ngrid * sizeof(char),
					copyType::ToDevice);
				d.deviceCalloc((void**)&(maxCalc.oscillatorIndexMap),
					maxCalc.Ngrid, sizeof(int64_t));
				d.deviceMemcpy(
					(void*)maxCalc.oscillatorIndexMap,
					(void*)oscillatorIndexMapCPU.data(),
					maxCalc.Ngrid * sizeof(int64_t),
					copyType::ToDevice);
				maxCalc.hasMaterialMap = true;
				maxCalc.NMaterialGrid = NmaterialPoints * maxCalc.Noscillators;

				for(int i = 0; i<NmaterialMax; i++){
					calculateFDTDParameters(sCPU, maxCalc, i);
				}
				if(maxCalc.observationPoint==0){
					maxCalc.observationPoint = maxCalc.Nz - 10;
				}
				
			}
			else{
				std::runtime_error("Failed to load material map.\n");
			}

		}
		else{
			maxCalc.materialKeys[0] = (*sCPU).materialIndex;
			calculateFDTDParameters(sCPU, maxCalc);
		}

		
		//Make sure that the time grid is populated and do a 1D (time) FFT onto the frequency grid
		//make both grids available through the maxCalc class
		d.deviceMemcpy(
			d.s->gridETime1, 
			(*sCPU).ExtOut, 
			2 * (*sCPU).Ngrid * sizeof(double), 
			copyType::ToDevice);
		maxCalc.inOutEx = d.s->gridETime1;
		maxCalc.inOutEy = d.s->gridETime2;
		d.fft(d.s->gridETime1, d.s->gridEFrequency1, deviceFFT::D2Z_1D);
		maxCalc.inputExFFT = reinterpret_cast<deviceFP*>(d.s->gridEFrequency1);
		maxCalc.inputEyFFT = reinterpret_cast<deviceFP*>(d.s->gridEFrequency2);
		d.deviceMemset(maxCalc.inOutEx, 0, 2*(*sCPU).Ngrid * sizeof(deviceFP));

		//allocate the new memory needed for the maxwell calculation
		d.deviceCalloc((void**)&(maxCalc.Egrid), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.EgridEstimate), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.EgridNext), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.EgridEstimate2), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.Hgrid), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.HgridEstimate), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.HgridNext), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.HgridEstimate2), 
			maxCalc.Ngrid, sizeof(maxwellPoint<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.materialGrid), 
			maxCalc.NMaterialGrid, sizeof(oscillator<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.materialGridEstimate), 
			maxCalc.NMaterialGrid, sizeof(oscillator<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.materialGridEstimate2), 
			maxCalc.NMaterialGrid, sizeof(oscillator<deviceFP>));
		d.deviceCalloc((void**)&(maxCalc.materialGridNext), 
			maxCalc.NMaterialGrid, sizeof(oscillator<deviceFP>));

		//make a device copy of the maxCalc class
		maxwell3D* maxCalcDevice{};
		d.deviceCalloc((void**)&maxCalcDevice, 1, sizeof(maxwell3D));
		d.deviceMemcpy(
			(void*)maxCalcDevice, 
			(void*)&maxCalc, 
			sizeof(maxwell3D), 
			copyType::ToDevice);
		maxCalc.deviceCopy = maxCalcDevice;

		//if there is an initial plasma population, set the values in the material grid
		d.deviceLaunch(
				maxCalc.Ngrid / 64, 
				64, 
				maxwellSetInitialCarrierDensity{ maxCalc.deviceCopy});
	}

	static unsigned long int solveFDTD(
		ActiveDevice& d, 
		simulationParameterSet* sCPU, 
		int64_t tFactor, 
		deviceFP dz, 
		deviceFP frontBuffer, 
		deviceFP backBuffer,
		deviceFP observationPoint,
		int gridIndex,
		bool preserveNearField,
		int64_t manualWaitFrames) {
		
		//initialize the grid if necessary
		if (!sCPU->isFollowerInSequence) {
			simulationParameterSet sCPUcopy = *sCPU;

			sCPUcopy.materialIndex = 0;
			sCPUcopy.crystalTheta = 0.0;
			sCPUcopy.crystalPhi = 0.0;
			sCPUcopy.crystalThickness = 0.0;
			sCPUcopy.propagationStep = 1e-9;

			sCPUcopy.nonlinearAbsorptionStrength = 0.0;
			sCPUcopy.chi2Tensor = sCPUcopy.crystalDatabase[0].d.data();
			sCPUcopy.chi3Tensor = sCPUcopy.crystalDatabase[0].chi3.data();
			sCPUcopy.nonlinearSwitches = 
				sCPUcopy.crystalDatabase[0].nonlinearSwitches;
			sCPUcopy.sellmeierCoefficients = 
				sCPUcopy.crystalDatabase[0].sellmeierCoefficients.data();

			sCPUcopy.sellmeierType = sCPUcopy.crystalDatabase[0].sellmeierType;
			sCPUcopy.axesNumber = 0;
			sCPUcopy.isFDTD = false;
			d.reset(&sCPUcopy);
			solveNonlinearWaveEquationWithDevice(d, &sCPUcopy);
			d.reset(sCPU);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			(*sCPU).isFollowerInSequence = true;
		}
		if (dz == 0.0) dz = (*sCPU).propagationStep;
		//generate the FDTD data structure and prepare the device
		maxwell3D maxCalc = maxwell3D(sCPU, tFactor, dz, frontBuffer, backBuffer);
		if(observationPoint != 0.0){
			maxCalc.observationPoint = static_cast<int>(round(observationPoint/dz));
			if(maxCalc.observationPoint < 1 || maxCalc.observationPoint >= maxCalc.Nz){
				throw std::runtime_error("Invalid observation point in FDTD\n");
			}
		}
		maxCalc.waitFrames = manualWaitFrames;

		prepareFDTD(d, sCPU, maxCalc, gridIndex);
		
		//RK loop
		for (int64_t i = 0; i < maxCalc.Nt; i++) {
			d.deviceLaunch(
				maxCalc.Ngrid / 64, 
				64, 
				maxwellRKkernel0{ maxCalc.deviceCopy, i });
			d.deviceLaunch(
				maxCalc.Ngrid / 64, 
				64, 
				maxwellRKkernel1{ maxCalc.deviceCopy, i });
			d.deviceLaunch(
				maxCalc.Ngrid / 64,
				64, 
				maxwellRKkernel2{ maxCalc.deviceCopy, i });
			d.deviceLaunch(
				maxCalc.Ngrid / 64, 
				64, 
				maxwellRKkernel3{ maxCalc.deviceCopy, i });
			if (i % maxCalc.tGridFactor == 0 
				&& (i >= maxCalc.waitFrames)) {
				d.deviceLaunch(
					(maxCalc.Nx * maxCalc.Ny) / minGridDimension, 
					minGridDimension, 
					maxwellSampleGrid{ 
						maxCalc.deviceCopy, 
						(i - maxCalc.waitFrames) / maxCalc.tGridFactor });
			}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			if (i % 20 == 0 && d.isTheCanaryPixelNaN(&(maxCalc.Egrid[0].y))) break;
			if ((*sCPU).cancellationCalled) break;
		}
		

		
		
		//correct far-field amplitudes for vectorial effects
		if(preserveNearField){
			d.deviceMemcpy(
				(*sCPU).ExtOut, 
				maxCalc.inOutEx, 
				2*(*sCPU).Ngrid * sizeof(double), 
				copyType::ToHost);
			d.fft(maxCalc.inOutEx, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
			d.deviceMemcpy(
				(*sCPU).EkwOut,
				d.deviceStruct.gridEFrequency1,
				2 * d.deviceStruct.NgridC * sizeof(std::complex<double>),
				copyType::ToHost);
		}
		else {
			d.fft(maxCalc.inOutEx, d.deviceStruct.gridEFrequency1, deviceFFT::D2Z);
			if((*sCPU).is3D){
				d.deviceLaunch(
					d.deviceStruct.Nblock / 2, 
					d.deviceStruct.Nthread,
					correctFDTDAmplitudesKernel{ d.dParamsDevice });
			}
			else{
				d.deviceLaunch(
					d.deviceStruct.Nblock / 2, 
					d.deviceStruct.Nthread,
					correctFDTDAmplitudesKernel2D{ d.dParamsDevice });
			}
			d.deviceMemcpy(
				(*sCPU).EkwOut,
				d.deviceStruct.gridEFrequency1,
				2 * d.deviceStruct.NgridC * sizeof(std::complex<double>),
				copyType::ToHost);
			d.fft(d.deviceStruct.gridEFrequency1, d.deviceStruct.gridETime1, deviceFFT::Z2D);
			//transfer result to CPU memory and take spectrum
			d.deviceMemcpy(
				(*sCPU).ExtOut, 
				d.deviceStruct.gridETime1, 
				2*(*sCPU).Ngrid * sizeof(double), 
				copyType::ToHost);
		}

		getTotalSpectrum(d);

		//free device memory		
		return freeFDTD(d, maxCalc);
	}

	static constexpr unsigned int functionID(const char* s) {
		unsigned int result = 0;
		for (int i = 0; s[i]; ++i) {
			result += s[i];
		}
		return result;
	}

	//Dispatcher of the sequence mode. 
	// New functions go here, and should have a unique ID (chances of a conflict are small, and 
	// will produce a compile-time error.
	// Functions cannot start with a number or the string "None".
	static int interpretCommand(
		const std::string& cc, 
		const double* iBlock, 
		double* vBlock, 
		ActiveDevice& d, 
		simulationParameterSet *sCPU) {
		if (cc.size() == 0) return 0;
		crystalEntry* db = (*sCPU).crystalDatabase;
		int error = 0;
		double parameters[32] = {};
		bool defaultMask[32] = {};
		auto functionNameEnd = cc.find_first_of('(');
		if (functionNameEnd ==std::string::npos) throw std::runtime_error(std::string("Didn't find opening parenthesis in\n").append(cc));
		switch (functionID(cc.substr(0, functionNameEnd).c_str())) {
		case functionID("rotate"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			rotateField(d, sCPU, deg2Rad<deviceFP>() * parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("rotateIntoBiaxial"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			rotateBiaxial(
				d, 
				sCPU, 
				defaultMask[0] ? (*sCPU).materialIndex : static_cast<int>(parameters[0]), 
				defaultMask[1] ? (*sCPU).crystalTheta : deg2Rad<deviceFP>() * parameters[1],
				defaultMask[2] ? (*sCPU).crystalPhi : deg2Rad<deviceFP>() * parameters[2],
				defaultMask[3] ? (*sCPU).pulse1.frequency : 1e12 * parameters[3],
				true);
			break;
		case functionID("rotateFromBiaxial"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			rotateBiaxial(
				d, 
				sCPU, 
				defaultMask[0] ? (*sCPU).materialIndex : static_cast<int>(parameters[0]), 
				defaultMask[1] ? (*sCPU).crystalTheta : deg2Rad<deviceFP>() * parameters[1],
				defaultMask[2] ? (*sCPU).crystalPhi : deg2Rad<deviceFP>() * parameters[2],
				defaultMask[3] ? (*sCPU).pulse1.frequency : 1e12 * parameters[3],
				false);
			break;
		case functionID("set"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			if(parameters[0]<100 && parameters[0] >= 0){ 
				vBlock[static_cast<int>(parameters[0])] = parameters[1];
			}
			else throw std::runtime_error("set() index must be less than 100\n");
			break;
		case functionID("plasmaReinject"):
			(*sCPU).isReinjecting = true;
			[[fallthrough]];
		case functionID("plasma"):
		{
			interpretParameters(cc, 9, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
			if (!defaultMask[3])(*sCPU).nonlinearAbsorptionStrength = parameters[3];
			if (!defaultMask[4])(*sCPU).bandGapElectronVolts = parameters[4];
			if (!defaultMask[5])(*sCPU).drudeGamma = parameters[5];
			if (!defaultMask[6])(*sCPU).effectiveMass = parameters[6];
			if (!defaultMask[7])(*sCPU).crystalThickness = 1e-6 * parameters[7];
			if (!defaultMask[8])(*sCPU).propagationStep = 1e-9 * parameters[8];
			(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d.data();
			(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3.data();
			(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();

			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU);
			(*sCPU).isFollowerInSequence = true;
		}
			break;
		case functionID("nonlinear"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
			if (!defaultMask[3])(*sCPU).crystalThickness = 1e-6 * parameters[3];
			if (!defaultMask[4])(*sCPU).propagationStep = 1e-9 * parameters[4];

			(*sCPU).nonlinearAbsorptionStrength = 0.0;
			(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d.data();
			(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3.data();
			(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();

			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU);
			(*sCPU).isFollowerInSequence = true;
			break;
		case functionID("fdtd"):
			interpretParameters(cc, 8, iBlock, vBlock, parameters, defaultMask);
			{
				int64_t timeFactor = 5;
				if(!defaultMask[0]) timeFactor = static_cast<int64_t>(parameters[0]);

				double dz = sCPU->propagationStep;
				if(!defaultMask[1]) dz = parameters[1];

				double frontBuffer = 1e-6;
				if(!defaultMask[2]) frontBuffer = parameters[2];

				double backBuffer = 1e-6;
				if(!defaultMask[3]) backBuffer = parameters[3];

				double observationPoint = 0.0;
				if(!defaultMask[4]) observationPoint = parameters[4];

				double dt = (sCPU->tStep)/timeFactor;
				int64_t waitFrames = -1;
				if(!defaultMask[5]) waitFrames = static_cast<int64_t>(round(parameters[5]/dt)); 

				bool preserveNearField = false;
				if(!defaultMask[6]) preserveNearField = (parameters[6]==1.0);

				int gridIndex = -1;
				if(!defaultMask[7]) gridIndex = static_cast<int>(parameters[7]);

				error = solveFDTD(
					d, 
					sCPU, 
					timeFactor, 
					dz, 
					frontBuffer, 
					backBuffer,
					observationPoint,
					gridIndex,
					preserveNearField,
					waitFrames);
			}
			break;
		case functionID("fdtdReflection"):
    		interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
    		{
				int64_t timeFactor = 5;
				if(!defaultMask[0]) timeFactor = static_cast<int64_t>(parameters[0]);

				double dz = sCPU->propagationStep;
				if(!defaultMask[1]) dz = parameters[1];

				//Set the front buffer size to fit the input field
                double offset = 1e-6;
                double frontBuffer = 0.5 * (sCPU->timeSpan)*lightC<double>();

				//Set the wait time so that the first moment of the input field
				//reaches the observation plane after reflecting
				double dt = (sCPU->tStep)/timeFactor;
                double waitTime = (offset + 2 * frontBuffer)/lightC<double>();
				int64_t waitFrames = static_cast<int64_t>(round(waitTime/dt));

                error = solveFDTD(
    			d,
    			sCPU,
    			timeFactor,
    			dz,
    			frontBuffer + offset,
    			offset,
                offset,
                -1,
                true,
                waitFrames);
            }
            break;
		case functionID("default"):
			d.reset(sCPU);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU);
			break;
		case functionID("save"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			{
				int64_t saveLoc = (int64_t)parameters[0];
				if (saveLoc < (*sCPU).Nsims && saveLoc != 0 && (*sCPU).runType != runTypes::counter) {
					memcpy(
						&(*sCPU).ExtOut[saveLoc * (*sCPU).Ngrid * 2], 
						(*sCPU).ExtOut, 
						2 * (*sCPU).Ngrid * sizeof(double));
					memcpy(
						&(*sCPU).EkwOut[saveLoc * (*sCPU).NgridC * 2], 
						(*sCPU).EkwOut, 
						2 * (*sCPU).NgridC * sizeof(std::complex<double>));
					memcpy(
						&(*sCPU).totalSpectrum[saveLoc * 3 * (*sCPU).Nfreq], 
						(*sCPU).totalSpectrum, 
						3 * (*sCPU).Nfreq * sizeof(double));
				}
				else if ((*sCPU).runType != runTypes::counter) {
					throw std::runtime_error(
						std::string("Attempted out-of-bounds save() to index ")
						.append(std::to_string(saveLoc)).append("\n"));
				}
			}
			break;
		case functionID("savePlasma"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			{
				int64_t saveLoc = (int64_t)parameters[0];
				int64_t plasmaLoc = (int)parameters[1];
				if (saveLoc < (*sCPU).Nsims 
					&& saveLoc != 0 
					&& (*sCPU).runType != runTypes::counter) 
					savePlasma(d, saveLoc, plasmaLoc);
				else if ((*sCPU).runType == runTypes::counter){
					return error;
				}
				else {
					throw std::runtime_error(
						std::string("Attempted out-of-bounds savePlasma() to index ")
						.append(std::to_string(saveLoc)).append("\n"));
				}
			}
			break;
		case functionID("init"):
			(*sCPU).materialIndex = 0;
			(*sCPU).crystalTheta = 0.0;
			(*sCPU).crystalPhi = 0.0;
			(*sCPU).crystalThickness = 0;
			(*sCPU).propagationStep = 1e-9;

			(*sCPU).nonlinearAbsorptionStrength = 0.0;
			(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d.data();
			(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3.data();
			(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();

			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			error = solveNonlinearWaveEquationWithDevice(d, sCPU);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			(*sCPU).isFollowerInSequence = true;
			break;

		case functionID("linear"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if ((*sCPU).isCylindric) {
				if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
				if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
				if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
				if (!defaultMask[3])(*sCPU).crystalThickness = 1e-6 * parameters[3];
				if (!defaultMask[4])(*sCPU).propagationStep = 1e-9 * parameters[4];
				(*sCPU).nonlinearAbsorptionStrength = 0.0;
				(*sCPU).chi2Tensor = db[(*sCPU).materialIndex].d.data();
				(*sCPU).chi3Tensor = db[(*sCPU).materialIndex].chi3.data();
				(*sCPU).nonlinearSwitches = db[(*sCPU).materialIndex].nonlinearSwitches;
				(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();
				(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
				(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
				(*sCPU).forceLinear = true;
				d.reset(sCPU);
				error = solveNonlinearWaveEquationWithDevice(d, sCPU);
				(*sCPU).isFollowerInSequence = true;
			}
			else {
				if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
				if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
				if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
				if (!defaultMask[3])(*sCPU).crystalThickness = 1e-6 * parameters[3];
				if (d.hasPlasma) {
					(*sCPU).nonlinearAbsorptionStrength = 0.0;
					(*sCPU).startingCarrierDensity = 0.0;
					(*sCPU).forceLinear = true;
				}
				(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();
				(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
				(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
				d.reset(sCPU);
				applyLinearPropagation(d, sCPU, (*sCPU).materialIndex, (*sCPU).crystalThickness);
				if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			}

			break;
		case functionID("fresnelLoss"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			if (!defaultMask[0])(*sCPU).materialIndex = (int)parameters[0];
			if (!defaultMask[1])(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[1];
			if (!defaultMask[2])(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[2];
			d.reset(sCPU);
			applyFresnelLoss(d, sCPU, d.deviceStruct,
				(int)parameters[4],
				(int)parameters[5]);
			break;
		case functionID("sphericalMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applySphericalMirror(d, sCPU, d.deviceStruct, parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("parabolicMirror"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyParabolicMirror(d, sCPU, d.deviceStruct, parameters[0]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("aperture"):
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyAperature(d, sCPU,
				parameters[0],
				parameters[1]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("farFieldAperture"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			(*sCPU).materialIndex = 0;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();
			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			applyAperatureFarField(d, sCPU,
				parameters[0],
				parameters[1],
				parameters[2],
				parameters[3]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("farFieldInverseAperture"):
			interpretParameters(cc, 4, iBlock, vBlock, parameters, defaultMask);
			(*sCPU).materialIndex = 0;
			(*sCPU).sellmeierCoefficients = db[(*sCPU).materialIndex].sellmeierCoefficients.data();
			(*sCPU).sellmeierType = db[(*sCPU).materialIndex].sellmeierType;
			(*sCPU).axesNumber = db[(*sCPU).materialIndex].axisType;
			d.reset(sCPU);
			applyInverseAperatureFarField(d, sCPU,
				parameters[0],
				parameters[1],
				parameters[2],
				parameters[3]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("energy"):
			{
			if ((*sCPU).runType == runTypes::counter) break;
			interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
			int targetVar = (int)parameters[0];
			int spectrumType = (int)parameters[1];
			double energy = 0.0;
			for(int i = 0; i < (*sCPU).Nfreq; i++){
				energy += (*sCPU).totalSpectrum[i + (*sCPU).Nfreq * spectrumType];
			}
			energy *= (*sCPU).fStep;
			if(targetVar<100) vBlock[targetVar] = energy;
			}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("filter"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyFilter(d, sCPU,
				parameters[0],
				parameters[1],
				parameters[2],
				parameters[3],
				parameters[4]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("lorentzian"):
			interpretParameters(cc, 5, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			applyLorenzian(
				d, 
				sCPU, 
				parameters[0], 
				parameters[1], 
				parameters[2], 
				parameters[3], 
				parameters[4]);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("applyOptic"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			if((*sCPU).optics.size() < static_cast<std::size_t>(parameters[0])+1){
				throw std::runtime_error("applyOptic: optic is not loaded");
			}
			applyOptic(
				d, 
				*sCPU, 
				static_cast<int>(parameters[0]),
				true,
				true);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("applyOpticX"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			if((*sCPU).optics.size() < static_cast<std::size_t>(parameters[0])+1){
				throw std::runtime_error("applyOpticX: optic is not loaded");
			}
			applyOptic(
				d, 
				*sCPU, 
				static_cast<int>(parameters[0]),
				true,
				false);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("applyOpticY"):
			interpretParameters(cc, 1, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			if((*sCPU).optics.size() < static_cast<std::size_t>(parameters[0])+1){
				throw std::runtime_error("applyOpticY: optic is not loaded");
			}
			applyOptic(
				d, 
				*sCPU, 
				static_cast<int>(parameters[0]),
				false,
				true);
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		case functionID("loadOptic"):{
				loadedInputData file(cc.substr(cc.find('"')+1, cc.find('"', cc.find('"') + 1) - cc.find('"') - 1));
				if(file.hasData) (*sCPU).optics.push_back(file);
				else throw std::runtime_error("loadOptic failed, check the file path. If you are using the Flatpak version, load the optic using the load optic button.");
			}
			break;
		case functionID("addPulse"):
			if ((*sCPU).runType == runTypes::counter) break;
		{
			interpretParameters(cc, 21, iBlock, vBlock, parameters, defaultMask);
			d.reset(sCPU);
			d.deviceMemcpy(
				d.deviceStruct.gridETime1, 
				(*sCPU).ExtOut, 
				2 * d.deviceStruct.Ngrid * sizeof(double), 
				copyType::ToDevice);
			d.deviceMemcpy(
				d.deviceStruct.gridEFrequency1, 
				(*sCPU).EkwOut, 
				2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
				copyType::ToDevice);

			pulse<double> p;
			p = sCPU->pulse1;
			p.energy = parameters[0];
			p.frequency = 1e12 * parameters[1];
			p.bandwidth = 1e12 * parameters[2];
			p.sgOrder = (int)parameters[3];
			p.cep = parameters[4] * vPi<deviceFP>();
			p.delay = 1e-15 * parameters[5];
			p.gdd = 1e-30 * parameters[6];
			p.tod = 1e-45 * parameters[7];
			p.phaseMaterial = (int)parameters[8];
			p.phaseMaterialThickness = 1e-6 * parameters[9];
			p.beamwaist = 1e-6 * parameters[10];
			p.x0 = 1e-6 * parameters[11];
			p.y0 = 1e-6 * parameters[12];
			p.z0 = 1e-6 *parameters[13];
			p.beamAngle = deg2Rad<deviceFP>() * parameters[14];
			p.beamAnglePhi = deg2Rad<deviceFP>() * parameters[15];
			p.polarizationAngle = deg2Rad<deviceFP>() * parameters[16];
			p.circularity = parameters[17];
			(*sCPU).materialIndex = (int)parameters[18];
			(*sCPU).crystalTheta = deg2Rad<deviceFP>() * parameters[19];
			(*sCPU).crystalPhi = deg2Rad<deviceFP>() * parameters[20];

			addPulseToFieldArrays(d, p, false, NULL);
			d.deviceMemcpy(
				(*sCPU).EkwOut, 
				d.deviceStruct.gridEFrequency1, 
				2 * d.deviceStruct.NgridC * sizeof(std::complex<double>), 
				copyType::ToHost);
			d.deviceMemcpy(
				(*sCPU).ExtOut, 
				d.deviceStruct.gridETime1, 
				2 * (*sCPU).Ngrid * sizeof(double), 
				copyType::ToHost);

			getTotalSpectrum(d);
		}
			if (!(*sCPU).isInFittingMode)(*(*sCPU).progressCounter)++;
			break;
		
		case functionID("for"): {
				interpretParameters(cc, 2, iBlock, vBlock, parameters, defaultMask);
				int counter = (int)parameters[0];
				int targetVar = (int)parameters[1];
				std::string currentString = cc.substr(cc.find_first_of('{') + 1, std::string::npos);
				std::string forStartString = currentString;
				vBlock[targetVar] = 0.0;
				for (int i = 0; i < counter; i++) {

					while (currentString.length() > 0 && currentString.at(0) != '}') {
						if (currentString.at(0) == '<') {
							currentString =
								currentString.substr(currentString.find_first_of('>'), std::string::npos);
							if (currentString.length() > 0) currentString =
								currentString.substr(1, std::string::npos);
						}
						if (currentString.at(0) == '{') {
							currentString = currentString.substr(1, std::string::npos);
							while (currentString.find_first_of('{') != std::string::npos
								&& currentString.find_first_of('{') < currentString.find_first_of('}')) {
								currentString =
									currentString.substr(currentString.find_first_of('}'), std::string::npos);
								currentString = currentString.substr(1, std::string::npos);
							}
							currentString =
								currentString.substr(currentString.find_first_of('}'), std::string::npos);
							if (currentString.length() < 5) break;
							currentString = currentString.substr(1, std::string::npos);
						}
						error = interpretCommand(currentString, iBlock, vBlock, d, sCPU);
						currentString = 
							currentString.substr(findParenthesesClosure(currentString)+1, std::string::npos);
						if (error || (*sCPU).cancellationCalled) break;
					}
					++vBlock[targetVar];
					currentString = forStartString;
					if (error || (*sCPU).cancellationCalled) break;
				}
			}
			break;
		default: throw std::runtime_error(std::string("Could not interpret:\n").append(cc));
		}
		return error;
	}
	
	static int solveSequenceWithDevice(
		ActiveDevice& d, 
		simulationParameterSet* sCPU) {

		int error = 0;
		//if it starts with 0, it's an old sequence; quit
		if ((*sCPU).sequenceString[0] == '0') {
			return 15;
		}

		//main text interpreter
		double iBlock[100] = { 0.0 };

		for (int k = 1; k < 38; k++) {
			iBlock[k] = (*sCPU).getByNumberWithMultiplier(k);
		}

		double vBlock[100] = { 0.0 };
		std::string currentString((*sCPU).sequenceString);
		simulationParameterSet backupSet = *sCPU;
		//shortest command is for(), if there's only 4 characters left, it can only
		//be whitespace or other trailing symbols, and even then something is wrong,
		//since for should be followed by a loop... next shortest is init(), and that
		//certainly doesn't belong at the end. So anything under six characters=bail
		size_t minLength = 5;
		while (currentString.length() > minLength) {
			//skip curly braces (for loops should have been handled by interpretCommand() already)
			if (currentString.at(0) == '{') {
				currentString = currentString.substr(1, std::string::npos);
				while (currentString.find_first_of('{') != std::string::npos
					&& currentString.find_first_of('{') < currentString.find_first_of('}')) {
					currentString = 
						currentString.substr(currentString.find_first_of('}'), std::string::npos);
					currentString = currentString.substr(1, std::string::npos);
				}
				currentString = 
					currentString.substr(currentString.find_first_of('}'), std::string::npos);
				if (currentString.length() < minLength) break;
				currentString = currentString.substr(1, std::string::npos);
			}
			//skip angle brackets (comments)
			while (currentString.at(0) == '<') {
				currentString = 
					currentString.substr(currentString.find_first_of('>'), std::string::npos);
				if (currentString.length() < minLength) break;
				currentString = currentString.substr(1, std::string::npos);
			}
			if (error || (*sCPU).cancellationCalled) break;
			error = interpretCommand(currentString, iBlock, vBlock, d, sCPU);
			if (error || (*sCPU).cancellationCalled) break;
			currentString = 
				currentString.substr(findParenthesesClosure(currentString), std::string::npos);
			if (currentString.length() < minLength) break;

			currentString = currentString.substr(1, std::string::npos);
			backupSet.isFollowerInSequence = (*sCPU).isFollowerInSequence;
			backupSet.optics = (*sCPU).optics;
			*sCPU = backupSet;
		}
		*sCPU = backupSet;
		return error;
	}

	// helper function for fitting mode, runs the simulation and returns difference from the desired outcome.
	static double getResidual(
		const dlib::matrix<double, 0, 1>& x){

		double result = 0.0;
		for (int i = 0; i < (*fittingSet).Nfitting; ++i) {
			(*fittingSet).setByNumberWithMultiplier((int64_t)(*fittingSet).fittingArray[3 * i], x(i));
		}

		ActiveDevice& d = *dFit;
		d.cParams = fittingSet;
		d.reset(fittingSet);

		if ((*fittingSet).isInSequence) {
			(*fittingSet).isFollowerInSequence = false;
			solveSequenceWithDevice(d, fittingSet);
			(*(*fittingSet).progressCounter)++;
		}
		else {
			solveNonlinearWaveEquationWithDevice(d, fittingSet);
		}

		//maximize total spectrum in ROI
		if ((*fittingSet).fittingMode < 3) {
			for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
				result += 
					(*fittingSet).totalSpectrum[(*fittingSet).fittingMode 
					* (*fittingSet).Nfreq + (*fittingSet).fittingROIstart + i];
			}
			return result;
		}

		//mode 3 & 4: match total spectrum to reference given in ascii file
		double a;
		double maxSim = 0.0;
		double maxRef = 0.0;
		double* simSpec = &(*fittingSet).totalSpectrum[
			2 * (*fittingSet).Nfreq + (*fittingSet).fittingROIstart];
		double* refSpec = &(*fittingSet).fittingReference[(*fittingSet).fittingROIstart];
		for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
			maxSim = maxN(maxSim, simSpec[i]);
			maxRef = maxN(maxRef, refSpec[i]);
		}

		if (maxSim == 0) {
			maxSim = 1;
		}
		if (maxRef == 0) {
			maxRef = 1;
		}
		result = 0.0;

		if ((*fittingSet).fittingMode == 4) {
			for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
				a = log10(refSpec[i] / maxRef) - log10(simSpec[i] / maxSim);
				if (!isnan(a)) result += a * a;
			}
		}
		else {
			for (int i = 0; i < (*fittingSet).fittingROIsize; ++i) {
				a = (refSpec[i] / maxRef) - (simSpec[i] / maxSim);
				result += a * a;
			}
		}
		return sqrt(result);
	}
}
using namespace hostFunctions;

//Main (non sequence) solver. New device typeps should have a unique definition of the name
// e.g. solveNonlinearWaveEquationSYCL so that the correct one may be called. That's why it's
// a preprocessor definition here.
unsigned long solveNonlinearWaveEquationX(simulationParameterSet* sCPU) {
	ActiveDevice d(sCPU);
	if (d.memoryStatus) return 1;
	return solveNonlinearWaveEquationWithDevice(d, sCPU);
}

// Main function for running a sequence
unsigned long solveNonlinearWaveEquationSequenceX(simulationParameterSet* sCPU) {
	if ((*sCPU).batchIndex == 36 && (*sCPU).batchLoc1 != 0) return 0;
	ActiveDevice d(sCPU);
	if (d.memoryStatus) return 1;
	return solveSequenceWithDevice(d, sCPU);
}

//run in fitting mode
unsigned long runDlibFittingX(simulationParameterSet* sCPU) {
	simulationParameterSet sCPUcurrent = *sCPU;
	fittingSet = &sCPUcurrent;
	ActiveDevice d(fittingSet);
	dFit = &d;
	dlib::matrix<double, 0, 1> parameters;
	parameters.set_size((*sCPU).Nfitting);
	dlib::matrix<double, 0, 1> lowerBounds;
	lowerBounds.set_size((*sCPU).Nfitting);
	dlib::matrix<double, 0, 1> upperBounds;
	upperBounds.set_size((*sCPU).Nfitting);
	for (int i = 0; i < (*sCPU).Nfitting; ++i) {
		parameters(i) = (*sCPU).getByNumber((int64_t)round((*sCPU).fittingArray[3 * i]));
		lowerBounds(i) = (*sCPU).fittingArray[3 * i + 1];
		upperBounds(i) = (*sCPU).fittingArray[3 * i + 2];
	}

	dlib::function_evaluation result;
	if ((*sCPU).fittingMode != 3) {
		result = dlib::find_max_global(
			getResidual, 
			lowerBounds, 
			upperBounds, 
			dlib::max_function_calls((*sCPU).fittingMaxIterations));
	}
	else {
		result = dlib::find_min_global(
			getResidual, 
			lowerBounds, 
			upperBounds, 
			dlib::max_function_calls((*sCPU).fittingMaxIterations));
	}

	for (int i = 0; i < (*sCPU).Nfitting; ++i) {
		(*sCPU).setByNumberWithMultiplier((int64_t)round((*sCPU).fittingArray[3 * i]), result.x(i));
		(*sCPU).fittingResult[i] = result.x(i);
	}

	std::atomic_uint32_t fitCounter{ 0 };
	std::atomic_uint32_t* originalCounter = (*sCPU).progressCounter;
	(*sCPU).progressCounter = &fitCounter;
	d.cParams = sCPU;
	d.reset(sCPU);
	if ((*sCPU).isInSequence) {
		solveSequenceWithDevice(d, sCPU);
	}
	else {
		solveNonlinearWaveEquationWithDevice(d, sCPU);
	}

	(*sCPU).progressCounter = originalCounter;
	dFit = nullptr;
	return 0;
}



