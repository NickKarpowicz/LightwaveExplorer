#include "LightwaveExplorerTrilingual.h"	
#include "LightwaveExplorerUtilities.h"
#include <iostream>
#include <chrono>
#include <thread>

#ifdef RUNONCUDA
	#include "LightwaveExplorerCore.cuh"
	#include "Devices/LightwaveExplorerCoreFP32.cuh"
#elif defined CPUONLY
	#include "Devices/LightwaveExplorerCoreCPU.h"
	#include "Devices/LightwaveExplorerCoreCPUFP32.h"
#elif defined RUNONSYCL
	#include "Devices/LightwaveExplorerSYCL.h"
	#include "Devices/LightwaveExplorerSYCLFP32.h"
#endif

int main(int argc, char* argv[]){
	char* filepath = argv[1];
	std::atomic_uint32_t progressCounter = 0;
	int CUDAdeviceCount = 1;
#ifdef __CUDACC__
	if (hardwareCheck(&CUDAdeviceCount)) return 1;
#endif
	if (argc < 2) {
		std::cout << "No input file specified." << std::endl;
		return 2;
	}

	// allocate databases, main structs
	simulationBatch theSim;
	crystalDatabase db;
	theSim.base().crystalDatabase = db.db.data();
	theSim.base().progressCounter = &progressCounter;
	// read crystal database
	if (db.db.size() == 0) {
		std::cout << "Could not read crystal database." << std::endl;
		return 12;
	}
	std::cout << "Read " << db.db.size() << " crystal database entries:" << std::endl;
	for (int j = 0; j < db.db.size(); ++j) {
		std::cout << "Material " << j << " name: " << db.db[j].crystalName << std::endl;
	}

	// read from settings file
	if (theSim.sCPU()->readInputParametersFile(db.db.data(), filepath) == 1) {
		std::cout << "Could not read input file." << std::endl;
		return 13;
	}

	// read optics files if there are any
	int opticCount = 0;
	while(true){
		loadedInputData newOptic(theSim.base().outputBasePath + "_optic" + std::to_string(opticCount) + ".txt");
		if(newOptic.hasData) theSim.optics.push_back(newOptic);
		else break;
		opticCount++;
	}
	
	if ((theSim.sCPU()->sequenceString[0] != 'N') && theSim.sCPU()->sequenceString[0] != 0) theSim.sCPU()->isInSequence = true;
	theSim.configure();
	simulationParameterSet* sCPU = theSim.sCPU();
	auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

	if ((*sCPU).Nfitting != 0) {
		std::cout << "Running optimization for" << (*sCPU).fittingMaxIterations << "iterations..." << std::endl;

		runDlibFittingX(sCPU);

		auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
		std::cout << "Finished after" <<
			1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count())
			<< "s" << std::endl;
		theSim.saveDataSet();

		std::cout << "Optimization result:" << std::endl << "(index, value)" << std::endl;
		for (int i = 0; i < (*sCPU).Nfitting; ++i) {
			std::cout << i << (*sCPU).fittingResult[i] << std::endl;
		}
		return 0;
	}
	// run simulations
	std::thread* threadBlock = new std::thread[(*sCPU).Nsims * (*sCPU).Nsims2];
	size_t maxThreads = minN(CUDAdeviceCount, (*sCPU).Nsims * (*sCPU).Nsims2);
	for (int j = 0; j < (*sCPU).Nsims * (*sCPU).Nsims2; ++j) {

		sCPU[j].assignedGPU = j % CUDAdeviceCount;
		if (j >= maxThreads) {
			if (threadBlock[j - maxThreads].joinable()) {
				threadBlock[j - maxThreads].join();
			}
		}

		if ((*sCPU).isInSequence) {
			threadBlock[j] = std::thread(solveNonlinearWaveEquationSequenceX, &sCPU[j]);
		}
		else {
			threadBlock[j] = std::thread(solveNonlinearWaveEquationX, &sCPU[j]);
		}

	}
	for (int i = 0; i < (*sCPU).Nsims * (*sCPU).Nsims2; ++i) {
		if (sCPU[i].memoryError > 0) {
			std::cout << "Warning: device memory error " << sCPU[i].memoryError << std::endl;
		}
		if (threadBlock[i].joinable()) {
			threadBlock[i].join();
		}
	}

	auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
	std::cout << "Finished after " <<
		1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count())
		<< "s" << std::endl;
	theSim.saveDataSet();
	delete[] threadBlock;
	return 0;
}
