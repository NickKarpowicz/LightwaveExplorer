#include "LightwaveExplorerTrilingual.h"	
#include "LightwaveExplorerUtilities.h"
#include <iostream>

	#ifndef PLATFORMTYPE
	#include "../LightwaveExplorerCore.cuh"
	#include <chrono>
	#include <thread>
	#endif
//main function - if included in the GUI, this should have a different name
// than main() - this one applies when running in command line mode (e.g. on
// the clusters.)
int main(int argc, char* argv[]){
	char* filepath = argv[1];
	int i, j;

	size_t progressCounter = 0;
	int CUDAdeviceCount = 1;

	if (hardwareCheck(&CUDAdeviceCount)) return 1;

	if (argc < 2) {
		std::cout << "No input file specified." << std::endl;
		return 2;
	}

	// allocate databases, main structs
	simulationParameterSet initializationStruct;
	memset(&initializationStruct, 0, sizeof(simulationParameterSet));
	crystalDatabase db;
	initializationStruct.crystalDatabase = db.db.data();
	initializationStruct.progressCounter = &progressCounter;
	// read crystal database
	if (db.db.size() == 0) {
		std::cout << "Could not read crystal database." << std::endl;
		return 12;
	}
	std::cout << "Read " << db.db.size() << "crystal database entries :" << std::endl;
	for (j = 0; j < db.db.size(); ++j) {
		std::cout << "Material " << j << "name: " << db.db[j].crystalName.c_str() << std::endl;
	}

	// read from settings file
	if (readInputParametersFile(&initializationStruct, db.db.data(), filepath) == 1) {
		std::cout << "Could not read input file." << std::endl;
		return 13;
	}
	size_t Ntotal = minN(1, initializationStruct.Nsims * minN(1, initializationStruct.Nsims2));
	simulationParameterSet* sCPU = new simulationParameterSet[Ntotal];
	memcpy(sCPU, &initializationStruct, sizeof(simulationParameterSet));
	allocateGrids(sCPU);
	if (loadPulseFiles(sCPU) == 1) {
		std::cout << "Could not read pulse file." << std::endl;
		deallocateGrids(sCPU, true);
		delete[] sCPU;
		return 14;
	}

	if (((*sCPU).sequenceString[0] != 'N') && (*sCPU).sequenceString[0] != 0) (*sCPU).isInSequence = true;
	configureBatchMode(sCPU);
	readFittingString(sCPU);
	auto simulationTimerBegin = std::chrono::high_resolution_clock::now();

	if ((*sCPU).Nfitting != 0) {
		std::cout << "Running optimization for" << (*sCPU).fittingMaxIterations << "iterations..." << std::endl;

		runDlibFitting(sCPU);

		auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
		std::cout << "Finished after" <<
			1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count())
			<< "s" << std::endl;
		saveDataSet(sCPU);

		std::cout << "Optimization result:" << std::endl << "(index, value)" << std::endl;
		for (int i = 0; i < (*sCPU).Nfitting; ++i) {
			std::cout << i << (*sCPU).fittingResult[i] << std::endl;
		}

		deallocateGrids(sCPU, true);
		delete[] sCPU;
		return 0;
	}
	// run simulations
	std::thread* threadBlock = new std::thread[Ntotal];
	size_t maxThreads = minN(CUDAdeviceCount, Ntotal);
	for (j = 0; j < Ntotal; ++j) {

		sCPU[j].assignedGPU = j % CUDAdeviceCount;
		if (j >= maxThreads) {
			if (threadBlock[j - maxThreads].joinable()) {
				threadBlock[j - maxThreads].join();
			}
		}

		if ((*sCPU).isInSequence) {
			threadBlock[j] = std::thread(solveNonlinearWaveEquationSequence, &sCPU[j]);
		}
		else {
			threadBlock[j] = std::thread(solveNonlinearWaveEquation, &sCPU[j]);
		}

	}
	for (i = 0; i < (*sCPU).Nsims * (*sCPU).Nsims2; ++i) {
		if (sCPU[i].memoryError > 0) {
			std::cout << "Warning: device memory error " << sCPU[i].memoryError << std::endl;
		}
		if (threadBlock[i].joinable()) {
			threadBlock[i].join();
		}
	}

	auto simulationTimerEnd = std::chrono::high_resolution_clock::now();
	std::cout << "Finished after" <<
		1e-6 * (double)(std::chrono::duration_cast<std::chrono::microseconds>(simulationTimerEnd - simulationTimerBegin).count())
		<< "s" << std::endl;
	saveDataSet(sCPU);
	delete[] threadBlock;
	deallocateGrids(sCPU, true);
	delete[] sCPU;
	return 0;
}
