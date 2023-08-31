#include "LightwaveExplorerUtilities.h"
#include <iostream>
#ifdef CPUONLY
#include <fftw3.h>
#else
#include <fftw3_mkl.h>
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

int simulationParameterSet::loadSavedFields(const std::string& outputBase) {
	std::string Epath = outputBase;
	Epath.append("_Ext.dat");
	std::ifstream Efile(Epath, std::ios::binary);
	if (Efile.is_open()) {
		Efile.read(reinterpret_cast<char*>(ExtOut), 2 * (Ngrid * Nsims * Nsims2) * sizeof(double));
	}
	else return 1;

	std::string Spath = outputBase;
	Spath.append("_spectrum.dat");
	std::ifstream Sfile(Spath, std::ios::binary);
	if (Sfile.is_open()) Sfile.read(
		reinterpret_cast<char*>(totalSpectrum), 
		Nsims * Nsims2 * 3 * Nfreq * sizeof(double));

	fftw_plan fftwPlanD2Z;
	if (is3D) {
		const int fftwSizes[] = { (int)Nspace2, (int)Nspace, (int)Ntime };
		fftwPlanD2Z = fftw_plan_many_dft_r2c(
			3, 
			fftwSizes, 
			2, 
			ExtOut, 
			NULL, 
			1, 
			(int)Ngrid, 
			(fftw_complex*)EkwOut, 
			NULL, 
			1, 
			(int)NgridC, 
			FFTW_MEASURE);
	}
	else {
		const int fftwSizes[] = { (int)Nspace, (int)Ntime };
		fftwPlanD2Z = fftw_plan_many_dft_r2c(
			2, 
			fftwSizes, 
			2, 
			ExtOut, 
			NULL, 
			1, 
			(int)Ngrid, 
			(fftw_complex*)EkwOut, 
			NULL, 
			1, 
			(int)NgridC, 
			FFTW_MEASURE);
	}

	for (int64_t i = 0; i < (Nsims * Nsims2); i++) {
		fftw_execute_dft_r2c(
			fftwPlanD2Z, 
			&ExtOut[2 * i * Ngrid], 
			(fftw_complex*) &EkwOut[2 * i * NgridC]);
	}
	fftw_destroy_plan(fftwPlanD2Z);

	return 0;
}


//note to self: replace raw pointers with std::vector
int simulationParameterSet::loadReferenceSpectrum() {
	std::ifstream fs(fittingPath);
	if (fs.fail()) {
		return 1;
	}
	int64_t maxFileSize = 16384;
	int64_t currentRow = 0;
	constexpr double c = 1e9 * lightC<double>();
	std::vector<double> loadedWavelengths(1);
	loadedWavelengths.reserve(8192);
	std::vector<double> loadedFrequencies(1);
	loadedFrequencies.reserve(8192);
	std::vector<double> loadedIntensities(1);
	loadedIntensities.reserve(8192);
	double maxWavelength{};
	double minWavelength{};
	double lastRead;
	while (fs.good() && currentRow < maxFileSize) {
		//fs >> loadedWavelengths[currentRow] >> loadedIntensities[currentRow];
		fs >> lastRead;
		loadedWavelengths.push_back(lastRead);
		fs >> lastRead;
		loadedIntensities.push_back(lastRead);
		if (currentRow == 0) {
			maxWavelength = loadedWavelengths[currentRow];
			minWavelength = loadedWavelengths[currentRow];
		}
		else {
			maxWavelength = maxN(maxWavelength, loadedWavelengths[currentRow]);
			minWavelength = minN(minWavelength, loadedWavelengths[currentRow]);
		}
		//rescale to frequency spacing
		loadedIntensities[currentRow] *= 
			loadedWavelengths[currentRow] * loadedWavelengths[currentRow];
		loadedFrequencies.push_back(c / loadedWavelengths[currentRow]);
		currentRow++;
	}
	int64_t sizeData = currentRow - 1;
	int64_t i, j;

	double maxFrequency = c / minWavelength;
	double minFrequency = c / maxWavelength;
	double currentFrequency = 0;
	double df;

	for (i = 1; i < Nfreq; i++) {
		currentFrequency = i * fStep;
		if ((currentFrequency > minFrequency) && (currentFrequency < maxFrequency)) {
			//find the first frequency greater than the current value
			j = sizeData - 1;
			while ((loadedFrequencies[j] <= currentFrequency) && (j > 2)) {
				j--;
			}
			df = loadedFrequencies[j] - loadedFrequencies[j - 1];
			fittingReference[i] =
				(loadedIntensities[j - 1] * (loadedFrequencies[j] - currentFrequency)
					+ loadedIntensities[j] * (currentFrequency - loadedFrequencies[j - 1]))
				/ df; 
			//linear interpolation
		}
	}
	return 0;
}

double simulationParameterSet::saveSlurmScript(int gpuType, int gpuCount, int64_t totalSteps) {
	std::string outputFile=outputBasePath;
	outputFile.append(".slurmScript");
	std::ofstream fs(outputFile, std::ios::binary);
	if (fs.fail()) return 1;

	//Estimate the time required on the cluster, proportional to number of grid points x steps
	double timeEstimate = static_cast<double>((totalSteps / (Nsims * Nsims2)) * Nspace * Ntime) * ceil(static_cast<double>(Nsims * Nsims2) / gpuCount);
	//3D
	if (symmetryType == 2) {
		timeEstimate *= Nspace2;
		timeEstimate *= 1.0e-11;
	}
	//Radial symmetry
	else if (symmetryType == 1) {
		timeEstimate *= 2.0e-11;
	}
	//2D
	else {
		timeEstimate *= 1.0e-11;
	}
	//Plasma doubles estimate
	if (nonlinearAbsorptionStrength != 0.0) timeEstimate *= 2.1;
	if (gpuType != 2) timeEstimate *= 8; //if it's not an A100, assume its slow.
	timeEstimate += 20.0; //fixed offset for loading .etc
	timeEstimate /= 3600.0; //convert to hours
	
	if (fittingString[0] != 0 && fittingString[0] != 'N') {
		readFittingString();
		if (fittingMaxIterations > 0) timeEstimate *= fittingMaxIterations;
	}
	std::string baseName = getBasename(outputBasePath);

	fs << "#!/bin/bash -l" << '\x0A';
	fs << "#SBATCH -o ./tjob.out.%j" << '\x0A';
	fs << "#SBATCH -e ./tjob.err.%j" << '\x0A';
	fs << "#SBATCH -D ./" << '\x0A';
	fs << "#SBATCH -J lightwave" << '\x0A';
	fs << "#SBATCH --constraint=\"gpu\"" << '\x0A';
	if (gpuType == 0) {
		fs << "#SBATCH --gres=gpu:rtx5000:" << minN(gpuCount, 2) << '\x0A';
	}
	if (gpuType == 1) {
		fs << "#SBATCH --gres=gpu:v100:" << minN(gpuCount, 2) << '\x0A';
	}
	if (gpuType == 2) {
		fs << "#SBATCH --gres=gpu:a100:" << minN(gpuCount, 4) << '\x0A';
		fs << "#SBATCH --cpus-per-task=" << 1 + minN(gpuCount, 4) << '\x0A';
	}
	fs << "#SBATCH --mem=" << 
		8192 + (18 * sizeof(double) * Ngrid * maxN(Nsims, 1u)) / 1048576 << "M\x0A";
	fs << "#SBATCH --nodes=1" << '\x0A';
	fs << "#SBATCH --ntasks-per-node=1" << '\x0A';
	fs << "#SBATCH --time=" << (int)ceil(1.5 * timeEstimate) << ":00:00" << '\x0A';
	fs << "module purge" << '\x0A';
	fs << "module load cuda/11.6" << '\x0A';
	fs << "module load mkl/2022.1" << '\x0A';
	fs << "export LD_LIBRARY_PATH=$MKL_HOME/lib/intel64:$LD_LIBRARY_PATH" << '\x0A';
	if (gpuType == 0 || gpuType == 1) {
		fs << "srun ./lwe " << baseName << ".input > " << baseName << ".out\x0A";
	}
	if (gpuType == 2) {
		fs << "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}" << '\x0A';
		fs << "srun ./lwe " << baseName << ".input > " << baseName << ".out\x0A";
	}

	return timeEstimate;
}

int simulationParameterSet::saveSettingsFile() {
	crystalEntry *crystalDatabasePtr = crystalDatabase;
	std::string outputPath(outputBasePath);
	if (runType > 0) {
		outputPath.append(".input");
	}
	else {
		outputPath.append(".txt");
	}
	std::ofstream fs(outputPath, std::ios::binary);
	if (fs.fail()) return -1;

	std::string baseName = getBasename(outputBasePath);
	std::string referenceBaseName = getBasename(fittingPath);
	std::string pulse1BaseName = getBasename(field1FilePath);
	std::string pulse2BaseName = getBasename(field2FilePath);
	fs.precision(15);

	fs << "Pulse energy 1 (J): " << pulse1.energy << '\x0A';
	fs << "Pulse energy 2 (J): " << pulse2.energy << '\x0A';
	fs << "Frequency 1 (Hz): " << pulse1.frequency << '\x0A';
	fs << "Frequency 2 (Hz): " << pulse2.frequency << '\x0A';
	fs << "Bandwidth 1 (Hz): " << pulse1.bandwidth << '\x0A';
	fs << "Bandwidth 2 (Hz): " << pulse2.bandwidth << '\x0A';
	fs << "SG order 1: " << pulse1.sgOrder << '\x0A';
	fs << "SG order 2: " << pulse2.sgOrder << '\x0A';
	fs << "CEP 1 (rad): " << pulse1.cep << '\x0A';
	fs << "CEP 2 (rad): " << pulse2.cep << '\x0A';
	fs << "Delay 1 (s): " << pulse1.delay << '\x0A';
	fs << "Delay 2 (s): " << pulse2.delay << '\x0A';
	fs << "GDD 1 (s^-2): " << pulse1.gdd << '\x0A';
	fs << "GDD 2 (s^-2): " << pulse2.gdd << '\x0A';
	fs << "TOD 1 (s^-3): " << pulse1.tod << '\x0A';
	fs << "TOD 2 (s^-3): " << pulse2.tod << '\x0A';
	fs << "Phase material 1 index: " << pulse1.phaseMaterial << '\x0A';
	fs << "Phase material 2 index: " << pulse2.phaseMaterial << '\x0A';
	fs << "Phase material thickness 1 (mcr.): " << pulse1.phaseMaterialThickness << '\x0A';
	fs << "Phase material thickness 2 (mcr.): " << pulse2.phaseMaterialThickness << '\x0A';
	fs << "Beam mode placeholder: " << 0 << '\x0A';
	fs << "Beamwaist 1 (m): " << pulse1.beamwaist << '\x0A';
	fs << "Beamwaist 2 (m): " << pulse2.beamwaist << '\x0A';
	fs << "x offset 1 (m): " << pulse1.x0 << '\x0A';
	fs << "x offset 2 (m): " << pulse2.x0 << '\x0A';
	fs << "y offset 1 (m): " << pulse1.y0 << '\x0A';
	fs << "y offset 2 (m): " << pulse2.y0 << '\x0A';
	fs << "z offset 1 (m): " << pulse1.z0 << '\x0A';
	fs << "z offset 2 (m): " << pulse2.z0 << '\x0A';
	fs << "NC angle 1 (rad): " << pulse1.beamAngle << '\x0A';
	fs << "NC angle 2 (rad): " << pulse2.beamAngle << '\x0A';
	fs << "NC angle phi 1 (rad): " << pulse1.beamAnglePhi << '\x0A';
	fs << "NC angle phi 2 (rad): " << pulse2.beamAnglePhi << '\x0A';
	fs << "Polarization 1 (rad): " << pulse1.polarizationAngle << '\x0A';
	fs << "Polarization 2 (rad): " << pulse2.polarizationAngle << '\x0A';
	fs << "Circularity 1: " << pulse1.circularity << '\x0A';
	fs << "Circularity 2: " << pulse2.circularity << '\x0A';
	fs << "Material index: " << materialIndex << '\x0A';
	fs << "Alternate material index: " << materialIndexAlternate << '\x0A';
	fs << "Crystal theta (rad): " << crystalTheta << '\x0A';
	fs << "Crystal phi (rad): " << crystalPhi << '\x0A';
	fs << "Grid width (m): " << spatialWidth << '\x0A';
	fs << "Grid height (m): " << spatialHeight << '\x0A';
	fs << "dx (m): " << rStep << '\x0A';
	fs << "Time span (s): " << timeSpan << '\x0A';
	fs << "dt (s): " << tStep << '\x0A';
	fs << "Thickness (m): " << crystalThickness << '\x0A';
	fs << "dz (m): " << propagationStep << '\x0A';
	fs << "Nonlinear absorption parameter: " << nonlinearAbsorptionStrength << '\x0A';
	fs << "Band gap (eV): " << bandGapElectronVolts << '\x0A';
	fs << "Effective mass (relative): " << effectiveMass << '\x0A';
	fs << "Drude gamma (Hz): " << drudeGamma << '\x0A';
	fs << "Propagation mode: " << symmetryType << '\x0A';
	fs << "Batch mode: " << batchIndex << '\x0A';
	fs << "Batch destination: " << batchDestination << '\x0A';
	fs << "Batch steps: " << Nsims << '\x0A';
	fs << "Batch mode 2: " << batchIndex2 << '\x0A';
	fs << "Batch destination 2: " << batchDestination2 << '\x0A';
	fs << "Batch steps 2: " << Nsims2 << '\x0A';
	fs << "Sequence: " << sequenceString << '\x0A';
	fs << "Fitting: " << fittingString << '\x0A';
	fs << "Fitting mode: " << fittingMode << '\x0A';
	if (runType > 0) { //don't include full path if making a cluster script
		fs << "Output base path: " << baseName << '\x0A';
		fs << "Field 1 from file type: " << pulse1FileType << '\x0A';
		fs << "Field 2 from file type: " << pulse2FileType << '\x0A';
		fs << "Field 1 file path: " << pulse1BaseName << '\x0A';
		fs << "Field 2 file path: " << pulse2BaseName << '\x0A';
		fs << "Fitting reference file path: " << referenceBaseName << '\x0A';
	}
	else {
		fs << "Output base path: " << outputBasePath << '\x0A';
		fs << "Field 1 from file type: " << pulse1FileType << '\x0A';
		fs << "Field 2 from file type: " << pulse2FileType << '\x0A';
		fs << "Field 1 file path: " << field1FilePath << '\x0A';
		fs << "Field 2 file path: " << field2FilePath << '\x0A';
		fs << "Fitting reference file path: " << fittingPath << '\x0A';
	}

	fs << "Material name: " << crystalDatabasePtr[materialIndex].crystalName << '\x0A';
	fs << "Sellmeier reference: " << crystalDatabasePtr[materialIndex].sellmeierReference << '\x0A';
	fs << "Chi2 reference: " << crystalDatabasePtr[materialIndex].dReference << '\x0A';
	fs << "Chi3 reference: " << crystalDatabasePtr[materialIndex].chi3Reference << '\x0A';
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 22; ++j) {
			fs << crystalDatabasePtr[materialIndex].sellmeierCoefficients[i * 22 + j];
			if (j < 21) fs << ',';
		}
		fs << '\x0A';
	}
	fs << "Code version: 2023.08";
	fs << '\x0A';

	return 0;
}

void simulationParameterSet::setByNumber(const int64_t index, const double value) {
	switch (index) {
	case 0:
		return;
	case 1:
		pulse1.energy = value; return;
	case 2:
		pulse2.energy = value; return;
	case 3:
		pulse1.frequency = value; return;
	case 4:
		pulse2.frequency = value; return;
	case 5:
		pulse1.bandwidth = value; return;
	case 6:
		pulse2.bandwidth = value; return;
	case 7:
		pulse1.cep = value; return;
	case 8:
		pulse2.cep = value; return;
	case 9:
		pulse1.delay = value; return;
	case 10:
		pulse2.delay = value; return;
	case 11:
		pulse1.gdd = value; return;
	case 12:
		pulse2.gdd = value; return;
	case 13:
		pulse1.tod = value; return;
	case 14:
		pulse2.tod = value; return;
	case 15:
		pulse1.phaseMaterialThickness = value; return;
	case 16:
		pulse2.phaseMaterialThickness = value; return;
	case 17:
		pulse1.beamwaist = value; return;
	case 18:
		pulse2.beamwaist = value; return;
	case 19:
		pulse1.x0 = value; return;
	case 20:
		pulse2.x0 = value; return;
	case 21:
		pulse1.z0 = value; return;
	case 22:
		pulse2.z0 = value; return;
	case 23:
		pulse1.beamAngle = value; return;
	case 24:
		pulse2.beamAngle = value; return;
	case 25:
		pulse1.polarizationAngle = value; return;
	case 26:
		pulse2.polarizationAngle = value; return;
	case 27:
		pulse1.circularity = value; return;
	case 28:
		pulse2.circularity = value; return;
	case 29:
		crystalTheta = value; return;
	case 30:
		crystalPhi = value; return;
	case 31:
		nonlinearAbsorptionStrength = value; return;
	case 32:
		drudeGamma = value; return;
	case 33:
		effectiveMass = value; return;
	case 34:
		crystalThickness = value; return;
	case 35:
		propagationStep = value; return;
	case 36:
		return;
	case 37:
		i37 = value; return;
	default:
		return;
	}
}

void simulationParameterSet::setByNumberWithMultiplier(
	const size_t index, 
	const double value) {
	if (index > multipliers.size()) return;
	setByNumber(index, value * multipliers[index]);
}


int simulationParameterSet::readInputParametersFile(
	crystalEntry* crystalDatabasePtr, 
	const std::string filePath) {
	std::ifstream fs(filePath);
	std::string line;

	if (fs.fail()) return 1;

	auto moveToColon = [&]() {
		char x = 0;
		while (x != ':' && fs.good()) {
			fs >> x;
		}
		return 0;
	};

	moveToColon();
	fs >> pulse1.energy;
	moveToColon();
	fs >> pulse2.energy;
	moveToColon();
	fs >> pulse1.frequency;
	moveToColon();
	fs >> pulse2.frequency;
	moveToColon();
	fs >> pulse1.bandwidth;
	moveToColon();
	fs >> pulse2.bandwidth;
	moveToColon();
	fs >> pulse1.sgOrder;
	moveToColon();
	fs >> pulse2.sgOrder;
	moveToColon();
	fs >> pulse1.cep;
	moveToColon();
	fs >> pulse2.cep;
	moveToColon();
	fs >> pulse1.delay;
	moveToColon();
	fs >> pulse2.delay;
	moveToColon();
	fs >> pulse1.gdd;
	moveToColon();
	fs >> pulse2.gdd;
	moveToColon();
	fs >> pulse1.tod;
	moveToColon();
	fs >> pulse2.tod;
	moveToColon();
	fs >> pulse1.phaseMaterial;
	moveToColon();
	fs >> pulse2.phaseMaterial;
	moveToColon();
	fs >> pulse1.phaseMaterialThickness;
	moveToColon();
	fs >> pulse2.phaseMaterialThickness;
	moveToColon();
	moveToColon();
	fs >> pulse1.beamwaist;
	moveToColon();
	fs >> pulse2.beamwaist;
	moveToColon();
	fs >> pulse1.x0;
	moveToColon();
	fs >> pulse2.x0;
	moveToColon();
	fs >> pulse1.y0;
	moveToColon();
	fs >> pulse2.y0;
	moveToColon();
	fs >> pulse1.z0;
	moveToColon();
	fs >> pulse2.z0;
	moveToColon();
	fs >> pulse1.beamAngle;
	moveToColon();
	fs >> pulse2.beamAngle;
	moveToColon();
	fs >> pulse1.beamAnglePhi;
	moveToColon();
	fs >> pulse2.beamAnglePhi;
	moveToColon();
	fs >> pulse1.polarizationAngle;
	moveToColon();
	fs >> pulse2.polarizationAngle;
	moveToColon();
	fs >> pulse1.circularity;
	moveToColon();
	fs >> pulse2.circularity;
	moveToColon();
	fs >> materialIndex;
	moveToColon();
	fs >> materialIndexAlternate;
	moveToColon();
	fs >> crystalTheta;
	moveToColon();
	fs >> crystalPhi;
	moveToColon();
	fs >> spatialWidth;
	moveToColon();
	fs >> spatialHeight;
	moveToColon();
	fs >> rStep;
	moveToColon();
	fs >> timeSpan;
	moveToColon();
	fs >> tStep;
	moveToColon();
	fs >> crystalThickness;
	moveToColon();
	fs >> propagationStep;
	moveToColon();
	fs >> nonlinearAbsorptionStrength;
	moveToColon();
	fs >> bandGapElectronVolts;
	moveToColon();
	fs >> effectiveMass;
	moveToColon();
	fs >> drudeGamma;
	moveToColon();
	fs >> symmetryType;
	moveToColon();
	fs >> batchIndex;
	moveToColon();
	fs >> batchDestination;
	moveToColon();
	fs >> Nsims;
	moveToColon();
	fs >> batchIndex2;
	moveToColon();
	fs >> batchDestination2;
	moveToColon();
	fs >> Nsims2;

	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());

	sequenceString = line;
	
	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());

	fittingString = line;
	moveToColon();
	fs >> fittingMode;

	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());

	outputBasePath = line;
	moveToColon();
	fs >> pulse1FileType;
	moveToColon();
	fs >> pulse2FileType;

	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());

	field1FilePath = line;
	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());

	field2FilePath = line;
	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());

	fittingPath = line;
	removeCharacterFromString(field1FilePath, '\r');
	removeCharacterFromString(field1FilePath, '\n');
	removeCharacterFromString(field2FilePath, '\r');
	removeCharacterFromString(field2FilePath, '\n');
	removeCharacterFromString(fittingPath, '\r');
	removeCharacterFromString(fittingPath, '\n');
	removeCharacterFromString(fittingString, '\r');
	removeCharacterFromString(fittingString, '\n');
	removeCharacterFromString(sequenceString, '\r');
	removeCharacterFromString(sequenceString, '\n');
	removeCharacterFromString(outputBasePath, '\r');
	removeCharacterFromString(outputBasePath, '\n');

	//derived parameters and cleanup:
	sellmeierType = 0;
	axesNumber = 0;
	Ntime = (int64_t)(minGridDimension * round(timeSpan / (minGridDimension * tStep)));
	Nfreq = Ntime / 2 + 1;
	Nspace = (int64_t)(minGridDimension * round(spatialWidth / (minGridDimension * rStep)));
	Nspace2 = (int64_t)(minGridDimension * round(spatialHeight / (minGridDimension * rStep)));
	Ngrid = Ntime * Nspace;
	NgridC = Nfreq * Nspace;
	kStep = twoPi<double>() / (Nspace * rStep);
	fStep = 1.0 / (Ntime * tStep);
	Npropagation = (int64_t)round(crystalThickness / propagationStep);

	isCylindric = symmetryType == 1;
	is3D = symmetryType == 2 || symmetryType == 4;
	isFDTD = symmetryType == 3 || symmetryType == 4;
	if (isCylindric) {
		pulse1.x0 = 0;
		pulse2.x0 = 0;
		pulse1.beamAngle = 0;
		pulse2.beamAngle = 0;
	}
	if (is3D) {
		Ngrid = Ntime * Nspace * Nspace2;
		NgridC = Nfreq * Nspace * Nspace2;
	}
	else {
		Nspace2 = 1;
	}
	if (batchIndex == 0 || Nsims < 1) {
		Nsims = 1;
	}
	if (batchIndex2 == 0 || Nsims2 < 1) {
		Nsims2 = 1;
	}

	field1IsAllocated = false;
	field2IsAllocated = false;

	//crystal from database (database must be loaded!)
	crystalDatabase = crystalDatabasePtr;
	chi2Tensor = crystalDatabasePtr[materialIndex].d.data();
	chi3Tensor = crystalDatabasePtr[materialIndex].chi3.data();
	nonlinearSwitches = crystalDatabasePtr[materialIndex].nonlinearSwitches.data();
	absorptionParameters = crystalDatabasePtr[materialIndex].absorptionParameters.data();
	sellmeierCoefficients = crystalDatabasePtr[materialIndex].sellmeierCoefficients.data();
	sellmeierType = crystalDatabasePtr[materialIndex].sellmeierType;
	axesNumber = crystalDatabasePtr[materialIndex].axisType;

	if (fs.good())return 61;
	else return -1;
}

void simulationBatch::configure() {
	Nfreq = parameters[0].Nfreq;
	Nsims = parameters[0].Nsims;
	Nsims2 = parameters[0].Nsims2;
	Nsimstotal = Nsims * Nsims2;
	Ngrid = parameters[0].Ngrid;
	NgridC = parameters[0].NgridC;
	simulationParameterSet base = parameters[0];
	parameters.resize(Nsimstotal, base);
	std::for_each(mutexes.begin(), mutexes.end(),
		[](std::mutex& m) {std::lock_guard<std::mutex> lock(m); });
	Ext = std::vector<double>(Nsimstotal * Ngrid * 2, 0.0);
	Ekw = std::vector<std::complex<double>>(
		Nsimstotal * NgridC * 2, 
		std::complex<double>(0.0, 0.0));
	mutexes = std::vector<std::mutex>(Nsimstotal);
	totalSpectrum = std::vector<double>(Nfreq * Nsimstotal * 3);

	if (parameters[0].pulse1FileType == 1 || parameters[0].pulse1FileType == 2) {
		loadedField1 = std::vector<std::complex<double>>(Nfreq, std::complex<double>(0.0, 0.0));
	}
	if (parameters[0].pulse2FileType == 1 || parameters[0].pulse2FileType == 2) {
		loadedField2 = std::vector<std::complex<double>>(Nfreq, std::complex<double>(0.0, 0.0));
	}
	if (parameters[0].pulse1FileType == 3) {
		loadedFullGrid1 = std::vector<double>(2 * Ngrid, 0.0);
	}
	if (parameters[0].pulse2FileType == 3) {
		loadedFullGrid2 = std::vector<double>(2 * Ngrid, 0.0);
	}
	if (parameters[0].fittingMode > 2) {
		fitReference = std::vector<double>(Nfreq, 0.0);
	}

	//configure
	double step1 = (parameters[0].batchDestination 
		- parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex)) 
		/ (Nsims - 1);
	double step2 = 0.0;
	if (Nsims2 > 0) {
		step2 = (parameters[0].batchDestination2 
			- parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex2)) 
			/ (Nsims2 - 1);
	}
	
	parameters[0].ExtOut = Ext.data();
	parameters[0].EkwOut = Ekw.data();
	parameters[0].totalSpectrum = totalSpectrum.data();
	parameters[0].loadedField1 = loadedField1.data();
	parameters[0].loadedField2 = loadedField2.data();
	parameters[0].loadedFullGrid1 = loadedFullGrid1.data();
	parameters[0].loadedFullGrid2 = loadedFullGrid2.data();
	parameters[0].fittingReference = fitReference.data();
	parameters[0].isGridAllocated = true;
	loadPulseFiles();

	for (int64_t i = 0; i < Nsims2; i++) {
		int64_t currentRow = i * Nsims;

		if (currentRow > 0) {
			parameters[currentRow] = parameters[0];
		}
		if (Nsims2 > 1) {
			parameters[currentRow].setByNumberWithMultiplier(
				parameters[0].batchIndex2, 
				parameters[0].getByNumberWithMultiplier(
					parameters[0].batchIndex2) + i * step2);
		}

		for (int64_t j = 0; j < Nsims; j++) {

			if (j > 0) {
				parameters[j + currentRow] = parameters[currentRow];
			}

			parameters[j + currentRow].batchLoc1 = j;
			parameters[j + currentRow].batchLoc2 = i;
			parameters[j + currentRow].ExtOut = getExt((j + currentRow));
			parameters[j + currentRow].EkwOut = getEkw((j + currentRow));
			parameters[j + currentRow].totalSpectrum = getTotalSpectrum((j + currentRow));
			parameters[j + currentRow].isFollowerInSequence = false;
			parameters[j + currentRow].setByNumberWithMultiplier(
				parameters[0].batchIndex, 
				parameters[0].getByNumberWithMultiplier(
					parameters[0].batchIndex) + j * step1);
		}
	}
}


void simulationBatch::configureCounter() {
	Nfreq = parameters[0].Nfreq;
	Nsims = parameters[0].Nsims;
	Nsims2 = parameters[0].Nsims2;
	Nsimstotal = Nsims * Nsims2;
	Ngrid = parameters[0].Ngrid;
	NgridC = parameters[0].NgridC;
	simulationParameterSet base = parameters[0];
	parameters.resize(Nsimstotal, base);
	std::for_each(mutexes.begin(), mutexes.end(),
		[](std::mutex& m) {std::lock_guard<std::mutex> lock(m); });
	mutexes = std::vector<std::mutex>(Nsimstotal);

	//configure
	double step1 = (parameters[0].batchDestination
		- parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex))
		/ (Nsims - 1);
	double step2 = 0.0;
	if (Nsims2 > 0) {
		step2 = (parameters[0].batchDestination2
			- parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex2))
			/ (Nsims2 - 1);
	}

	parameters[0].isGridAllocated = false;


	for (int64_t i = 0; i < Nsims2; i++) {
		int64_t currentRow = i * Nsims;

		if (currentRow > 0) {
			parameters[currentRow] = parameters[0];
		}
		if (Nsims2 > 1) {
			parameters[currentRow].setByNumberWithMultiplier(
				parameters[0].batchIndex2,
				parameters[0].getByNumberWithMultiplier(
					parameters[0].batchIndex2) + i * step2);
		}

		for (int64_t j = 0; j < Nsims; j++) {

			if (j > 0) {
				parameters[j + currentRow] = parameters[currentRow];
			}

			parameters[j + currentRow].batchLoc1 = j;
			parameters[j + currentRow].batchLoc2 = i;
			parameters[j + currentRow].isFollowerInSequence = false;
			parameters[j + currentRow].setByNumberWithMultiplier(
				parameters[0].batchIndex,
				parameters[0].getByNumberWithMultiplier(
					parameters[0].batchIndex) + j * step1);
		}
	}
}

void simulationBatch::loadPulseFiles() {
	//pulse type specifies if something has to be loaded to describe the pulses, or if they should be
	//synthesized later. 1: FROG .speck format; 2: Time-domain waveform; 3: Previous LWE result
	int frogLines = 0;
	if (parameters[0].pulse1FileType == 1) {
		frogLines = loadFrogSpeck(
			parameters[0].field1FilePath, 
			loadedField1.data(), 
			parameters[0].Ntime, 
			parameters[0].fStep, 
			0.0);
		parameters[0].field1IsAllocated = (frogLines > 1);
	}
	if (parameters[0].pulse2FileType == 1) {
		frogLines = loadFrogSpeck(
			parameters[0].field2FilePath, 
			loadedField2.data(), 
			parameters[0].Ntime, 
			parameters[0].fStep, 
			0.0);
		parameters[0].field1IsAllocated = (frogLines > 1);
	}

	if (parameters[0].pulse1FileType == 2) {
		frogLines = loadWaveformFile(
			parameters[0].field1FilePath, 
			loadedField1.data(), 
			parameters[0].Ntime, 
			parameters[0].fStep);
		parameters[0].field1IsAllocated = (frogLines > 1);
	}
	if (parameters[0].pulse2FileType == 2) {
		frogLines = loadWaveformFile(
			parameters[0].field2FilePath, 
			loadedField2.data(), 
			parameters[0].Ntime, 
			parameters[0].fStep);
		parameters[0].field2IsAllocated = (frogLines > 1);
	}

	if (parameters[0].pulse1FileType == 3) {
		parameters[0].field1IsAllocated = loadSavedGridFile(
			parameters[0].field1FilePath, 
			loadedFullGrid1, 
			parameters[0].Ngrid);
	}
	if (parameters[0].pulse2FileType == 3) {
		parameters[0].field2IsAllocated = loadSavedGridFile(
			parameters[0].field2FilePath, 
			loadedFullGrid2, 
			parameters[0].Ngrid);
	}
}

int simulationBatch::saveDataSet() {
	parameters[0].saveSettingsFile();
	std::for_each(mutexes.begin(), mutexes.end(),
		[&](std::mutex& m) { m.lock(); });
	std::string Epath=parameters[0].outputBasePath;
	Epath.append("_Ext.dat");
	std::ofstream Efile(Epath, std::ios::binary);
	if (Efile.is_open()) Efile.write(
		reinterpret_cast<char*>(Ext.data()), 
		Ext.size() * sizeof(double));

	std::string Spath=parameters[0].outputBasePath;
	Spath.append("_spectrum.dat");
	std::ofstream Sfile(Spath, std::ios::binary);
	if (Sfile.is_open()) Sfile.write(
		reinterpret_cast<char*>(totalSpectrum.data()), 
		totalSpectrum.size() * sizeof(double));
	std::for_each(mutexes.begin(), mutexes.end(),
		[&](std::mutex& m) { m.unlock(); });
	return 0;
}


int simulationParameterSet::readFittingString() {
	removeCharacterFromString(fittingString, '\r');
	removeCharacterFromString(fittingString, '\n');
	removeCharacterFromString(fittingString, '\t');
	std::stringstream ss(fittingString);
	double ROIbegin, ROIend;
	int maxIterations = 0;
	int fittingCount = 0;
	ss >> ROIbegin >> ROIend >> maxIterations;
	ss.ignore(fittingString.length(), ';');

	fittingROIstart = (int64_t)(ROIbegin / fStep);
	fittingROIstop = (int64_t)minN(ROIend / fStep, Ntime / 2);
	fittingROIsize = minN(maxN(fittingROIstop - fittingROIstart, 1u), Ntime / 2);
	fittingMaxIterations = maxIterations;

	while (ss.good()) {
		ss >> 
			fittingArray[fittingCount] >> 
			fittingArray[fittingCount + 1] >> 
			fittingArray[fittingCount + 2];
		if (ss.good()) fittingCount += 3;
		ss.ignore(fittingString.length(), ';');
	}

	Nfitting = fittingCount / 3;
	isInFittingMode = ((Nfitting > 0) && (maxIterations > 0));

	if (!isInFittingMode) {
		std::string noneString("None.");
		fittingString = noneString;
	}

	return 0;
}

int removeCharacterFromStringSkippingChars(
	std::string& s, 
	const char removedChar, 
	const char startChar, 
	const char endChar) {
	bool removing = true;
	for (size_t i = 0; i < s.length(); ++i) {
		if (s[i] == removedChar && removing) {
			s.erase(i,1);
			--i;
		}
		if (s[i] == startChar) removing = false;
		if (s[i] == endChar) removing = true;
	}
	return 0;
}

void stripWhiteSpace(std::string& s) {
	removeCharacterFromStringSkippingChars(s, ' ', '<', '>');
	removeCharacterFromString(s, '\r');
	removeCharacterFromString(s, '\n');
	removeCharacterFromString(s, '\t');
}

void stripLineBreaks(std::string& s) {
	removeCharacterFromString(s, '\r');
	removeCharacterFromString(s, '\n');
}

int interpretParameters(
	const std::string& cc, 
	const int numberParams, 
	const double *iBlock, 
	const double *vBlock, 
	double *parameters, 
	bool* defaultMask){
		std::string arguments = cc.substr(cc.find_first_of('(')+1, std::string::npos);
		// pattern: search for a , or ) depending on how many have been found
		// and the value of numberParams. 
		// If an ) is encountered while searching for , throw "too few"
		// If an , is encountered while searching for ) throw "too many"
		int numberFound = 0;
		size_t startArgument = 0;
		std::vector<std::string> argTable;
		argTable.reserve(numberParams);
		char expectedDelimiter = ',';
		char wrongDelimiter = ')';
		int openParens = 0;

		for(int i = 0; i<arguments.size(); ++i){
			if(numberFound == numberParams-1) {
				expectedDelimiter = ')';
				wrongDelimiter = ',';
			}
			if(arguments[i] == '(') openParens++;
			if(arguments[i] == expectedDelimiter && openParens == 0){
				if(i != startArgument){
					argTable.push_back(arguments.substr(startArgument,i-startArgument));
					numberFound++;
					startArgument = i+1;
				}
				else{
					throw std::runtime_error("Malformed argument\n");
				}
				if(expectedDelimiter==')') break;
			}
			else if(openParens > 0 && arguments[i]==')'){
				openParens--;
			}
			else if(arguments[i] == wrongDelimiter){
				throw std::runtime_error(std::string("Wrong number of arguments\n").append(cc).append("\n"));
			}
		}

		for(int i = 0; i<argTable.size(); ++i){
			if(argTable[i].at(0)=='d'){
				defaultMask[i] = true;
			}
			else{
				parameters[i] = parameterStringToDouble(argTable[i],iBlock,vBlock);
			}
		}
		return 0;
}

size_t findParenthesesClosure(std::string& a){
	int nParen = 0;
	bool foundFirstParen = false;
	for(int i = 0; i<a.size(); ++i){
		if(a[i]=='('){
			nParen++;
			foundFirstParen = true;
		}
		if(a[i]==')'){
			nParen--;
			if(nParen==0 && foundFirstParen){
				return i;
			}
			else if(nParen < 0 || !foundFirstParen){
				throw std::runtime_error(std::string("Weird parenthesis in:\n").append(a).append("\n"));
			}
		}
	}
	throw std::runtime_error(std::string("Weird parenthesis in:\n").append(a).append("\n"));
}

double parameterStringToDouble(
const std::string& ss, 
	const double* iBlock, 
	const double* vBlock) {
	std::vector<size_t> operatorTable;
	std::vector<double> numberTable;
	int lastOperator = -1;
	bool lastNumberWasNotParenthesized = true;
	for(size_t i = 0; i<ss.size(); ++i){
		if(
		  (ss[i] == '+'
		|| ss[i] == '-'
		|| ss[i] == '*'
		|| ss[i] == '/'
		|| ss[i] == '^')
		&& i != 0){
			operatorTable.push_back(i);
			if(lastNumberWasNotParenthesized)numberTable.push_back(std::stod(ss.substr(lastOperator+1,i-lastOperator-1)));
			lastOperator = i;
			lastNumberWasNotParenthesized = true;
		}
		else if(ss[i] == '('){
			std::string parenString = ss.substr(i, std::string::npos);
			parenString = parenString.substr(1,
				findParenthesesClosure(
					parenString)-1);
			numberTable.push_back(parameterStringToDouble(parenString, iBlock, vBlock));
			lastOperator = i + parenString.size();
			i += parenString.size()+1;
			lastNumberWasNotParenthesized = false;
		}
		else if(ss[i] == 'e'){
			if(ss[i] == '+' || ss[i] == '-') i += 2;
			else i++;
		}
		else if(ss[i] == 'v'){
			int ind = std::stoi(ss.substr(i+1,2));
			numberTable.push_back(vBlock[ind]);
			lastOperator = i + 2;
			i += 2;
			lastNumberWasNotParenthesized = false;
		}
		else if(ss[i] == 'i'){
			int ind = std::stoi(ss.substr(i+1,2));
			numberTable.push_back(iBlock[ind]);
			lastOperator = i + 2;
			i += 2;
			lastNumberWasNotParenthesized = false;
		}
	}
	if(lastOperator == -1) {
		return std::stod(ss);
	}
	if(lastOperator+1 < ss.size()){
			numberTable.push_back(std::stod(ss.substr(lastOperator+1,std::string::npos)));
		}

	auto applyOperators = [&](char firstOp, char secondOp){
		for(int i = 0; i < operatorTable.size(); ++i){
			if(ss[operatorTable[i]] == firstOp || ss[operatorTable[i]] == secondOp){
				switch(ss[operatorTable[i]]){
				case '^':
					numberTable[i] = std::pow(numberTable[i],numberTable[i+1]);
					break;
				case '*':
					numberTable[i] = numberTable[i] * numberTable[i+1];
					break;
				case '/':
					numberTable[i] = numberTable[i] / numberTable[i+1];
					break;
				case '+':
					numberTable[i] = numberTable[i] + numberTable[i+1];
					break;
				case '-':
					numberTable[i] = numberTable[i] - numberTable[i+1];
					break;
				}
				numberTable.erase(numberTable.begin() + i + 1);
				operatorTable.erase(operatorTable.begin()+i);
				i--;
			}
		}
	};
	applyOperators('^',0);
	applyOperators('*','/');
	applyOperators('+','-');
	if(std::isnan(numberTable[0])) throw std::runtime_error(
		std::string("NaN encountered when evaluating:\n").append(ss).append("\n"));
	return numberTable[0];
}

std::string getBasename(const std::string& fullPath) {
	std::string pathString = fullPath;
	int64_t positionOfName = pathString.find_last_of("/\\");
	if (positionOfName == std::string::npos) return pathString;
	return pathString.substr(positionOfName + 1);
}

//calculates the squared modulus of a complex number
inline double cModulusSquared(const std::complex<double>& x) {
	return x.real()*x.real() + x.imag()*x.imag();
}

int loadFrogSpeck(
	const std::string& frogFilePath, 
	std::complex<double>* Egrid, 
	const int64_t Ntime, 
	const double fStep, 
	const double gateLevel) {
	std::string line;
	std::ifstream fs(frogFilePath);
	if (fs.fail()) return -1;
	const int maxFileSize = 16384;
	double wavelength, R, phi, complexX, complexY, f, f0, f1;
	double fmax{};

	constexpr double cNanometers = 1e9 * lightC<double>(); 
	double df{};
	double fmin{};
	int currentRow{};
	std::vector<std::complex<double>> E;
	E.reserve(8192);

	while (fs.good() && currentRow < maxFileSize) {
		fs >> wavelength;
		fs >> R;
		fs >> phi;
		fs >> complexX;
		fs >> complexY;
		std::getline(fs, line);

		//get the complex field from the data
		E.push_back(std::complex<double>(complexX, complexY));
		//keep track of the frequency step of the grid 
		//(running sum, divide by number of rows at end to get average)
		if (currentRow > 0) df += cNanometers / wavelength - fmax;

		//keep track of the highest frequency in the data
		fmax = cNanometers / wavelength;

		//store the lowest frequency in the data
		if (currentRow == 0) fmin = fmax;

		currentRow++;
	}

	//return an error if nothing was loaded
	if (currentRow == 0) {
		return -1;
	}

	df /= currentRow; //average frequency step

	//interpolate the FROG data onto the simulation grid

	//fill the simulation grid based on the data
	for (int i = 0; i < Ntime / 2 + 1; i++) {

		//frequency grid used in the simulation
		f = i * fStep;

		int k0 = (int)floor((f - fmin) / df);
		int k1 = (int)ceil((f - fmin) / df);
		if (k0 < 0 || k1 >= currentRow) {
			Egrid[i] = 0; //field is zero outside of data range
		}
		else {
			f0 = fmin + k0 * df;
			f1 = fmin + k1 * df;
			Egrid[i] = (E[k0] * (f1 - f) + E[k1] * (f - f0)) / df; //linear interpolation
			Egrid[i] *= (abs(Egrid[i]) > gateLevel);
		}
	}
	return currentRow;
}

int loadWaveformFile(
	const std::string& filePath, 
	std::complex<double>* outputGrid,
	const int64_t Ntime, 
	const double fStep) {
	std::vector<double> Ein;
	Ein.reserve(8192);
	std::ifstream fs(filePath);
	std::string line;

	//read the waveform file: assumption is that the first column 
	//is time and second is the waveform.
	double dataT;
	double dataE;
	double dataDeltaT;
	double maxE2 = 0.0;
	int maxLoc = 0;
	int lineCount = 0;
	while (fs.good()) {
		fs >> dataT;
		fs >> dataE;
		if (dataE * dataE > maxE2) {
			maxE2 = dataE * dataE;
			maxLoc = lineCount;
		}
		std::getline(fs, line);
		if (lineCount == 0) dataDeltaT = dataT;
		if (lineCount == 1) dataDeltaT = dataT - dataDeltaT;
		Ein.push_back(dataE);
		lineCount++;
	}
	if (lineCount == 0) return 0;

	//frequency grid of the data
	const double df = 1.0 / (dataDeltaT * lineCount);
	const int64_t NfreqData = lineCount / 2 + 1;

	//FFT the waveform onto a frequency grid
	std::vector<std::complex<double>> fftOfEin(NfreqData + 1, 0.0);
	fftw_plan fftwPlanD2Z = fftw_plan_dft_r2c_1d(
		lineCount, Ein.data(), 
		reinterpret_cast<fftw_complex*>(fftOfEin.data()), 
		FFTW_ESTIMATE);
	fftw_execute_dft_r2c(
		fftwPlanD2Z, 
		Ein.data(), 
		reinterpret_cast<fftw_complex*>(fftOfEin.data()));
	fftw_destroy_plan(fftwPlanD2Z);

	//apply a time shift so that the frequency-domain solution 
	//oscillates slowly (will be undone after interpolation)
	const std::complex<double> timeShift = 
		std::complex<double>(0.0, 1.0) * twoPi<double>() 
		* df * static_cast<double>(maxLoc) * dataDeltaT;
	const std::complex<double> timeShiftResult = 
		std::complex<double>(0.0, -1.0) * twoPi<double>() 
		* static_cast<double>(maxLoc - lineCount/2) * dataDeltaT;
	for (int i = 0; i < NfreqData; i++) {
		fftOfEin[i] *= std::exp(timeShift * static_cast<double>(i));
	}

	//Interpolate in the frequency domain
	for (int i = 0; i < Ntime / 2 + 1; i++) {
		//frequency grid used in the simulation
		double f = i * fStep;

		int k0 = (int)floor(f / df);
		int k1 = (int)ceil(f / df);
		if (k0 < 0 || k1 >= NfreqData) {
			outputGrid[i] = {}; //field is zero outside of data range
		}
		else {
			double f0 = k0 * df;
			double f1 = k1 * df;
			//linear interpolation
			outputGrid[i] = (fftOfEin[k0] * (f1 - f) + fftOfEin[k1] * (f - f0)) / df; 
			outputGrid[i] *= std::exp(timeShiftResult * f);
		}
	}

	return lineCount;
}

int loadSavedGridFile(
	const std::string& filePath, 
	std::vector<double>& outputGrid, 
	const int64_t Ngrid) {
	std::ifstream Efile(filePath, std::ios::binary);
	outputGrid.resize(Ngrid);
	if (Efile.is_open()) {
		Efile.read(
			reinterpret_cast<char*>(outputGrid.data()), 
			2 * Ngrid * sizeof(double));
		return 0;
	}
	else return 1;
}

int loadSavedGridFileMultiple(
	const std::string& filePath, 
	std::vector<double>& outputGrid, 
	const int64_t Ngrid, 
	const int64_t Nsims) {
	outputGrid.resize(Ngrid * Nsims);
	std::ifstream Efile(filePath, std::ios::binary);
	if (Efile.is_open()) {
		Efile.read(reinterpret_cast<char*>(outputGrid.data()), 2 * Ngrid * Nsims * sizeof(double));
		return 0;
	}
	else return 1;
}