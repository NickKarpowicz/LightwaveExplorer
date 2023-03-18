#include <cstdlib>
#include <string.h>
#include <complex>
#include <fstream>
#include "LightwaveExplorerUtilities.h"


int readFittingString(simulationParameterSet* sCPU) {
	std::string sIn((*sCPU).fittingString);
#ifndef __CUDACC__
	std::erase(sIn, '\r');
	std::erase(sIn, '\n');
	std::erase(sIn, '\t');
	memset((*sCPU).fittingString, 0, MAX_LOADSTRING);
	sIn.copy((*sCPU).fittingString, MAX_LOADSTRING - 1);
#endif


	std::stringstream ss(sIn);
	double ROIbegin, ROIend;
	int maxIterations = 0;
	int fittingCount = 0;
	ss >> ROIbegin >> ROIend >> maxIterations;
	ss.ignore(MAX_LOADSTRING, ';');

	(*sCPU).fittingROIstart = (size_t)(ROIbegin / (*sCPU).fStep);
	(*sCPU).fittingROIstop = (size_t)minN(ROIend / (*sCPU).fStep, (*sCPU).Ntime / 2);
	(*sCPU).fittingROIsize = minN(maxN((*sCPU).fittingROIstop - (*sCPU).fittingROIstart, 1u), (*sCPU).Ntime / 2);
	(*sCPU).fittingMaxIterations = maxIterations;

	while (ss.good()) {
		ss >> (*sCPU).fittingArray[fittingCount] >> (*sCPU).fittingArray[fittingCount + 1] >> (*sCPU).fittingArray[fittingCount + 2];
		if (ss.good()) fittingCount += 3;
		ss.ignore(MAX_LOADSTRING, ';');
	}

	(*sCPU).Nfitting = fittingCount / 3;
	(*sCPU).isInFittingMode = (((*sCPU).Nfitting > 0) && (maxIterations > 0));

	if (!(*sCPU).isInFittingMode) {
		std::string noneString("None.\0");
		noneString.copy((*sCPU).fittingString, MAX_LOADSTRING);
	}

	return 0;
}

int skipFileUntilCharacter(FILE* fstream, char target) {
	int c = 0;
	while (c != EOF && c != target) {
		c = fgetc(fstream);
	}
	return 0;
}

int removeCharacterFromString(char* cString, size_t N, char removedChar) {
	size_t i = 0;
	size_t r = 0;
	if(cString[0] == 0) return 0;
	while (i < N - 1) {
		if (cString[i] == removedChar) {
			memmove(&cString[i], &cString[i + 1], N - i - r - 1);
			cString[N - r - 1] = 0;
			r++;
		}
		else {
			i++;
		}
	}
	if (cString[N - 1] == removedChar) {
		cString[N - 1] = '\0';
	}
	return 0;
}




int removeCharacterFromStringSkippingChars(std::string& s, char removedChar, char startChar, char endChar) {
	bool removing = TRUE;
	for (size_t i = 0; i < s.length(); ++i) {
		if (s[i] == removedChar && removing) {
			s.erase(i,1);
			--i;
		}
		if (s[i] == startChar) removing = FALSE;
		if (s[i] == endChar) removing = TRUE;
	}
	return 0;
}


void stripWhiteSpace(char* sequenceString, size_t bufferSize) {
	std::string s(sequenceString);
	removeCharacterFromStringSkippingChars(s, ' ', '<', '>');
#ifndef __CUDACC__
	std::erase(s, '\r');
	std::erase(s, '\n');
	std::erase(s, '\t');
#endif


	memset(sequenceString, 0, bufferSize);
	s.copy(sequenceString, bufferSize - 1);
}

void stripLineBreaks(char* sequenceString, size_t bufferSize) {
	std::string s(sequenceString);
#ifndef __CUDACC__
	std::erase(s, '\r');
	std::erase(s, '\n');
#endif
	memset(sequenceString, 0, bufferSize);
	s.copy(sequenceString, bufferSize-1);
}


int interpretParameters(std::string cc, int n, double *iBlock, double *vBlock, double *parameters, bool* defaultMask){
	std::string ccSegment = cc.substr(cc.find_first_of('(')+1, std::string::npos);

	std::stringstream ss(ccSegment);
	std::string line;
	int i = 0;
	while(i<n && ss.good()){
		std::getline(ss, line, ',');
		if(line.at(0)=='d'){
			defaultMask[i] = TRUE;
		}
		else{
			parameters[i] = parameterStringToDouble(line.c_str(),iBlock,vBlock);
		}
		i++;
	}

	return 0;
}

int copyParamsIntoStrings(char parameterBlock[22][256], const char* cc, int n) {
	int loc = 0;
	//scan to opening parenthesis
	for (;;) {
		if (cc[loc] == '(') break;
		if (cc[loc] == 0) return 1;
		loc++;
	}

	loc++;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < 64; j++) {
			if (cc[loc] == ',' || cc[loc] == ')') {
				j = 64;
			}
			else {
				parameterBlock[i][j] = cc[loc];
			}
			loc++;
		}
	}
	return loc;
}

void applyOp(char op, double* result, double* readout) {
	switch (op) {
	case '*':
		*result *= *readout;
		return;
	case '/':
		*result /= *readout;
		return;
	case '-':
		*result -= *readout;
		return;
	case '+':
		*result += *readout;
		return;
	case '^':
		if (*readout < 1.0 && *result < 0.0) {
			*result = 0.0;
		}
		else {
			*result = pow(*result, *readout);
		}
		return;
	}
}

double parameterStringToDouble(const char* pString, double* iBlock, double* vBlock) {
	std::string ss(pString);

	auto nextInt = [&](std::string iStr, int location) {
		std::stringstream s(iStr.substr(location));
		int a;
		s >> a;
		return a;
	};

	auto nextDouble = [&](std::string iStr, int location) {
		std::stringstream s(iStr.substr(location));
		double a;
		s >> a;
		return a;
	};

	int loc = 0;
	double result = 0.0;
	double readout = 0.0;
	int ind = 0;
	bool previousCharWasOp = 0;
	char lastOp = 0;
	while (loc < (ss.length()) && ss[loc]!=')') {
		if (!previousCharWasOp) {
			if (ss[loc] == 'v') {
				++loc;
				ind = nextInt(ss,loc);
				loc += 2;
				if (ind < 100) result = vBlock[ind];
			}
			else if (ss[loc] == 'i') {
				++loc;
				ind = nextInt(ss, loc);
				loc += 2;
				if (ind < 100) result = iBlock[ind];
			}

			else if (ss[loc] == '*'
				|| ss[loc] == '-'
				|| ss[loc] == '+'
				|| ss[loc] == '/'
				|| ss[loc] == '^') {
				previousCharWasOp = 1;
				lastOp = ss[loc];
				loc++;
			}
			else {
				result = nextDouble(ss, loc);
				while (!(ss[loc] == '*'
					|| ss[loc] == '-'
					|| ss[loc] == '+'
					|| ss[loc] == '/'
					|| ss[loc] == '^')) {
					if (loc >= (ss.length()-1) || ss[loc]==')') {
						return result;
					}
					if (ss[loc] == 'e') {
						if (ss.at(loc+1) == '-' || ss.at(loc+1) == '+') loc += 2;
					}
					else {
						loc++;
					}
				}
			}
		}
		else {
			if (ss[loc] == 'v') {
				++loc;
				ind = nextInt(ss, loc);
				loc += 2;
				if (ind < 100)readout = vBlock[ind];
				applyOp(lastOp, &result, &readout);
				previousCharWasOp = 0;
			}
			else if (ss[loc] == 'i') {
				++loc;
				ind = nextInt(ss, loc);
				loc += 2;
				if (ind < 100) readout = iBlock[ind];
				applyOp(lastOp, &result, &readout);
				previousCharWasOp = 0;
			}
			else {
				readout = nextDouble(ss, loc);
				applyOp(lastOp, &result, &readout);
				previousCharWasOp = 0;
				while (!(ss[loc] == '*'
					|| ss[loc] == '-'
					|| ss[loc] == '+'
					|| ss[loc] == '/'
					|| ss[loc] == '^')) {
					if (loc >= (ss.length()-1) || ss[loc]==')') {
						return result;
					}
					if (ss[loc] == 'e') {
						if (ss.at(loc+1) == '-' || ss.at(loc+1) == '+') loc += 2;
					}
					else {
						loc++;
					}
				}
			}
		}

	}
	return result;
}

//c implementation of fftshift, working on complex double precision
//A is the input array, B is the output
//dim1: column length
//dim2: row length
int fftshiftZ(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2) {
	long long i, j;
	long long div1 = dim1 / 2;
	long long div2 = dim2 / 2;
	//Quadrant 1
	for (i = 0; i < div1; i++) {
		for (j = 0; j < div2; j++) {
			B[i + dim1 * j] = A[i + div1 + dim1 * (j + div2)];
		}
	}
	//Quadrant 2
	for (i = 0; i < div1; i++) {
		for (j = div2; j < dim2; j++) {
			B[i + dim1 * j] = A[i + div1 + dim1 * (j - div2)];
		}
	}
	//Quadrant 3
	for (i = div1; i < dim1; i++) {
		for (j = 0; j < div2; j++) {
			B[i + dim1 * j] = A[i - div1 + dim1 * (j + div2)];
		}
	}
	//Quadrant 4
	for (i = div1; i < dim1; i++) {
		for (j = div2; j < dim2; j++) {
			B[i + dim1 * j] = A[i - div1 + dim1 * (j - div2)];
		}
	}
	return 0;
}

//c implementation of fftshift, working on complex double precision
//A is the input array, B is the output
//dim1: column length
//dim2: row length
int fftshiftD2Z(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2) {
	long long j;
	long long div2 = dim2 / 2;
	//Quadrant 1

	for (j = 0; j < div2; j++) {
		memcpy(&B[dim1 * j], &A[dim1 * (j + div2)], dim1 * 2 * sizeof(double));
	}
	//Quadrant 2

	for (j = div2; j < dim2; j++) {
		memcpy(&B[dim1 * j], &A[dim1 * (j - div2)], dim1 * 2 * sizeof(double));
	}


	return 0;
}

//same as fftshiftZ, but flips the output array columns
int fftshiftAndFilp(std::complex<double>* A, std::complex<double>* B, long long dim1, long long dim2) {
	long long i, j;
	long long div1 = dim1 / 2;
	long long div2 = dim2 / 2;
	//Quadrant 1
	for (i = 0; i < div1; i++) {
		for (j = 0; j < div2; j++) {
			B[(dim1 - i - 1) + dim1 * j] = A[i + div1 + dim1 * (j + div2)];
		}
	}
	//Quadrant 2
	for (i = 0; i < div1; i++) {
		for (j = div2; j < dim2; j++) {
			B[(dim1 - i - 1) + dim1 * j] = A[i + div1 + dim1 * (j - div2)];
		}
	}
	//Quadrant 3
	for (i = div1; i < dim1; i++) {
		for (j = 0; j < div2; j++) {
			B[(dim1 - i - 1) + dim1 * j] = A[i - div1 + dim1 * (j + div2)];
		}
	}
	//Quadrant 4
	for (i = div1; i < dim1; i++) {
		for (j = div2; j < dim2; j++) {
			B[(dim1 - i - 1) + dim1 * j] = A[i - div1 + dim1 * (j - div2)];
		}
	}
	return 0;
}

int loadReferenceSpectrum(char* spectrumPath, simulationParameterSet* sCPU) {
	std::ifstream fs(spectrumPath);
	if (fs.fail()) {
		return 1;
	}
	size_t maxFileSize = 16384;
	size_t currentRow = 0;
	double c = 1e9 * lightC<double>();
	double* loadedWavelengths = new double[maxFileSize]();
	double* loadedFrequencies = new double[maxFileSize]();
	double* loadedIntensities = new double[maxFileSize]();
	if (loadedWavelengths == NULL) {
		return 1;
	}
	double maxWavelength = 0;
	double minWavelength = 0;

	while (fs.good() && currentRow < maxFileSize) {
		fs >> loadedWavelengths[currentRow] >> loadedIntensities[currentRow];
		if (currentRow == 0) {
			maxWavelength = loadedWavelengths[currentRow];
			minWavelength = loadedWavelengths[currentRow];
		}
		else {
			maxWavelength = maxN(maxWavelength, loadedWavelengths[currentRow]);
			minWavelength = minN(minWavelength, loadedWavelengths[currentRow]);
		}
		//rescale to frequency spacing
		loadedIntensities[currentRow] *= loadedWavelengths[currentRow] * loadedWavelengths[currentRow];
		loadedFrequencies[currentRow] = c / loadedWavelengths[currentRow];
		currentRow++;
	}
	size_t sizeData = currentRow - 1;
	size_t i, j;

	double maxFrequency = c / minWavelength;
	double minFrequency = c / maxWavelength;
	double currentFrequency = 0;
	double df;

	for (i = 1; i < (*sCPU).Nfreq; i++) {
		currentFrequency = i * (*sCPU).fStep;
		if ((currentFrequency > minFrequency) && (currentFrequency < maxFrequency)) {
			//find the first frequency greater than the current value
			j = sizeData - 1;
			while ((loadedFrequencies[j] <= currentFrequency) && (j > 2)) {
				j--;
			}
			df = loadedFrequencies[j] - loadedFrequencies[j - 1];
			(*sCPU).fittingReference[i] =
				(loadedIntensities[j - 1] * (loadedFrequencies[j] - currentFrequency)
					+ loadedIntensities[j] * (currentFrequency - loadedFrequencies[j - 1])) / df; //linear interpolation
		}
	}

	delete[] loadedWavelengths;
	delete[] loadedIntensities;
	delete[] loadedFrequencies;

	return 0;
}

int loadSavedFields(simulationParameterSet* sCPU, const char* outputBase) {
	std::string Epath(outputBase);
	Epath.append("_Ext.dat");
	std::ifstream Efile(Epath, std::ios::binary);
	if (Efile.is_open()) {
		Efile.read(reinterpret_cast<char*>((*sCPU).ExtOut), 2 * ((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2) * sizeof(double));
	}
	else return 1;

	std::string Spath(outputBase);
	Spath.append("_spectrum.dat");
	std::ifstream Sfile(Spath, std::ios::binary);
	if (Sfile.is_open()) Sfile.read(reinterpret_cast<char*>((*sCPU).totalSpectrum), (*sCPU).Nsims * (*sCPU).Nsims2 * 3 * (*sCPU).Nfreq * sizeof(double));

	fftw_plan fftwPlanD2Z;
	if ((*sCPU).is3D) {
		const int fftwSizes[] = { (int)(*sCPU).Nspace2, (int)(*sCPU).Nspace, (int)(*sCPU).Ntime };
		fftwPlanD2Z = fftw_plan_many_dft_r2c(3, fftwSizes, 2, (*sCPU).ExtOut, NULL, 1, (int)(*sCPU).Ngrid, (fftw_complex*)(*sCPU).EkwOut, NULL, 1, (int)(*sCPU).NgridC, FFTW_MEASURE);
	}
	else {
		const int fftwSizes[] = { (int)(*sCPU).Nspace, (int)(*sCPU).Ntime };
		fftwPlanD2Z = fftw_plan_many_dft_r2c(2, fftwSizes, 2, (*sCPU).ExtOut, NULL, 1, (int)(*sCPU).Ngrid, (fftw_complex*)(*sCPU).EkwOut, NULL, 1, (int)(*sCPU).NgridC, FFTW_MEASURE);
	}

	for (size_t i = 0; i < ((*sCPU).Nsims * (*sCPU).Nsims2); i++) {
		fftw_execute_dft_r2c(fftwPlanD2Z, &(*sCPU).ExtOut[2 * i * (*sCPU).Ngrid], (fftw_complex*) & (*sCPU).EkwOut[2 * i * (*sCPU).NgridC]);
	}
	fftw_destroy_plan(fftwPlanD2Z);

	return 0;
}

std::string getBasename(char* fullPath) {
	std::string pathString(fullPath);
	std::size_t positionOfName = pathString.find_last_of("/\\");
	if (positionOfName == std::string::npos) return pathString;
	return pathString.substr(positionOfName + 1);
}

double saveSlurmScript(simulationParameterSet* sCPU, int gpuType, int gpuCount, size_t totalSteps) {
	std::string outputFile((*sCPU).outputBasePath);
	outputFile.append(".slurmScript");
	std::ofstream fs(outputFile, std::ios::binary);
	if (fs.fail()) return 1;

	//Estimate the time required on the cluster, proportional to number of grid points x steps
	double timeEstimate = (double)(totalSteps * (*sCPU).Nspace * (*sCPU).Ntime) * ceil(((float)(*sCPU).Nsims)/gpuCount);
	if ((*sCPU).symmetryType == 2) {
		timeEstimate *= (*sCPU).Nspace2;
		timeEstimate *= 1.2e-09;
	}
	else if ((*sCPU).symmetryType == 1) {
		timeEstimate *= 4.7e-9;
	}
	else {
		timeEstimate *= 4.6e-9;
	}
	if ((*sCPU).nonlinearAbsorptionStrength != 0.0) timeEstimate *= 2.1;
	timeEstimate /= 3600.0;
	if (gpuType != 2) timeEstimate *= 8;
	if ((*sCPU).fittingString[0] != 0 && (*sCPU).fittingString[0] != 'N') {
		readFittingString(sCPU);
		if ((*sCPU).fittingMaxIterations > 0) timeEstimate *= (*sCPU).fittingMaxIterations;
	}
	std::string baseName = getBasename((*sCPU).outputBasePath);

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
		fs << "#SBATCH --cpus-per-task=" << 2 * minN(gpuCount, 4) << '\x0A';
	}
	fs << "#SBATCH --mem=" << 8192 + (18 * sizeof(double) * (*sCPU).Ngrid * maxN((*sCPU).Nsims, 1u)) / 1048576 << "M\x0A";
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

//print a linefeed without a carriage return so that linux systems don't complain
//about impure scripts from DOS machines
//fopen() should be called with "wb"
void unixNewLine(FILE* iostream) {
	char LF = '\x0A';
	fwrite(&LF, sizeof(char), 1, iostream);
}

int saveSettingsFile(simulationParameterSet* sCPU) {
	crystalEntry *crystalDatabasePtr = (*sCPU).crystalDatabase;
	std::string outputPath((*sCPU).outputBasePath);
	if ((*sCPU).runType > 0) {
		outputPath.append(".input");
	}
	else {
		outputPath.append(".txt");
	}
	std::ofstream fs(outputPath, std::ios::binary);
	if (fs.fail()) return -1;

	std::string baseName = getBasename((*sCPU).outputBasePath);
	std::string referenceBaseName = getBasename((*sCPU).fittingPath);
	std::string pulse1BaseName = getBasename((*sCPU).field1FilePath);
	std::string pulse2BaseName = getBasename((*sCPU).field2FilePath);
	fs.precision(15);

	fs << "Pulse energy 1 (J): " << (*sCPU).pulse1.energy << '\x0A';
	fs << "Pulse energy 2 (J): " << (*sCPU).pulse2.energy << '\x0A';
	fs << "Frequency 1 (Hz): " << (*sCPU).pulse1.frequency << '\x0A';
	fs << "Frequency 2 (Hz): " << (*sCPU).pulse2.frequency << '\x0A';
	fs << "Bandwidth 1 (Hz): " << (*sCPU).pulse1.bandwidth << '\x0A';
	fs << "Bandwidth 2 (Hz): " << (*sCPU).pulse2.bandwidth << '\x0A';
	fs << "SG order 1: " << (*sCPU).pulse1.sgOrder << '\x0A';
	fs << "SG order 2: " << (*sCPU).pulse2.sgOrder << '\x0A';
	fs << "CEP 1 (rad): " << (*sCPU).pulse1.cep << '\x0A';
	fs << "CEP 2 (rad): " << (*sCPU).pulse2.cep << '\x0A';
	fs << "Delay 1 (s): " << (*sCPU).pulse1.delay << '\x0A';
	fs << "Delay 2 (s): " << (*sCPU).pulse2.delay << '\x0A';
	fs << "GDD 1 (s^-2): " << (*sCPU).pulse1.gdd << '\x0A';
	fs << "GDD 2 (s^-2): " << (*sCPU).pulse2.gdd << '\x0A';
	fs << "TOD 1 (s^-3): " << (*sCPU).pulse1.tod << '\x0A';
	fs << "TOD 2 (s^-3): " << (*sCPU).pulse2.tod << '\x0A';
	fs << "Phase material 1 index: " << (*sCPU).pulse1.phaseMaterial << '\x0A';
	fs << "Phase material 2 index: " << (*sCPU).pulse2.phaseMaterial << '\x0A';
	fs << "Phase material thickness 1 (mcr.): " << (*sCPU).pulse1.phaseMaterialThickness << '\x0A';
	fs << "Phase material thickness 2 (mcr.): " << (*sCPU).pulse2.phaseMaterialThickness << '\x0A';
	fs << "Beam mode placeholder: " << 0 << '\x0A';
	fs << "Beamwaist 1 (m): " << (*sCPU).pulse1.beamwaist << '\x0A';
	fs << "Beamwaist 2 (m): " << (*sCPU).pulse2.beamwaist << '\x0A';
	fs << "x offset 1 (m): " << (*sCPU).pulse1.x0 << '\x0A';
	fs << "x offset 2 (m): " << (*sCPU).pulse2.x0 << '\x0A';
	fs << "y offset 1 (m): " << (*sCPU).pulse1.y0 << '\x0A';
	fs << "y offset 2 (m): " << (*sCPU).pulse2.y0 << '\x0A';
	fs << "z offset 1 (m): " << (*sCPU).pulse1.z0 << '\x0A';
	fs << "z offset 2 (m): " << (*sCPU).pulse2.z0 << '\x0A';
	fs << "NC angle 1 (rad): " << (*sCPU).pulse1.beamAngle << '\x0A';
	fs << "NC angle 2 (rad): " << (*sCPU).pulse2.beamAngle << '\x0A';
	fs << "NC angle phi 1 (rad): " << (*sCPU).pulse1.beamAnglePhi << '\x0A';
	fs << "NC angle phi 2 (rad): " << (*sCPU).pulse2.beamAnglePhi << '\x0A';
	fs << "Polarization 1 (rad): " << (*sCPU).pulse1.polarizationAngle << '\x0A';
	fs << "Polarization 2 (rad): " << (*sCPU).pulse2.polarizationAngle << '\x0A';
	fs << "Circularity 1: " << (*sCPU).pulse1.circularity << '\x0A';
	fs << "Circularity 2: " << (*sCPU).pulse2.circularity << '\x0A';
	fs << "Material index: " << (*sCPU).materialIndex << '\x0A';
	fs << "Alternate material index: " << (*sCPU).materialIndexAlternate << '\x0A';
	fs << "Crystal theta (rad): " << (*sCPU).crystalTheta << '\x0A';
	fs << "Crystal phi (rad): " << (*sCPU).crystalPhi << '\x0A';
	fs << "Grid width (m): " << (*sCPU).spatialWidth << '\x0A';
	fs << "Grid height (m): " << (*sCPU).spatialHeight << '\x0A';
	fs << "dx (m): " << (*sCPU).rStep << '\x0A';
	fs << "Time span (s): " << (*sCPU).timeSpan << '\x0A';
	fs << "dt (s): " << (*sCPU).tStep << '\x0A';
	fs << "Thickness (m): " << (*sCPU).crystalThickness << '\x0A';
	fs << "dz (m): " << (*sCPU).propagationStep << '\x0A';
	fs << "Nonlinear absorption parameter: " << (*sCPU).nonlinearAbsorptionStrength << '\x0A';
	fs << "Band gap (eV): " << (*sCPU).bandGapElectronVolts << '\x0A';
	fs << "Effective mass (relative): " << (*sCPU).effectiveMass << '\x0A';
	fs << "Drude gamma (Hz): " << (*sCPU).drudeGamma << '\x0A';
	fs << "Propagation mode: " << (*sCPU).symmetryType << '\x0A';
	fs << "Batch mode: " << (*sCPU).batchIndex << '\x0A';
	fs << "Batch destination: " << (*sCPU).batchDestination << '\x0A';
	fs << "Batch steps: " << (*sCPU).Nsims << '\x0A';
	fs << "Batch mode 2: " << (*sCPU).batchIndex2 << '\x0A';
	fs << "Batch destination 2: " << (*sCPU).batchDestination2 << '\x0A';
	fs << "Batch steps 2: " << (*sCPU).Nsims2 << '\x0A';
	fs << "Sequence: " << (*sCPU).sequenceString << '\x0A';
	fs << "Fitting: " << (*sCPU).fittingString << '\x0A';
	fs << "Fitting mode: " << (*sCPU).fittingMode << '\x0A';
	if ((*sCPU).runType > 0) { //don't include full path if making a cluster script
		fs << "Output base path: " << baseName << '\x0A';
		fs << "Field 1 from file type: " << (*sCPU).pulse1FileType << '\x0A';
		fs << "Field 2 from file type: " << (*sCPU).pulse2FileType << '\x0A';
		fs << "Field 1 file path: " << pulse1BaseName << '\x0A';
		fs << "Field 2 file path: " << pulse2BaseName << '\x0A';
		fs << "Fitting reference file path: " << referenceBaseName << '\x0A';
	}
	else {
		fs << "Output base path: " << (*sCPU).outputBasePath << '\x0A';
		fs << "Field 1 from file type: " << (*sCPU).pulse1FileType << '\x0A';
		fs << "Field 2 from file type: " << (*sCPU).pulse2FileType << '\x0A';
		fs << "Field 1 file path: " << (*sCPU).field1FilePath << '\x0A';
		fs << "Field 2 file path: " << (*sCPU).field2FilePath << '\x0A';
		fs << "Fitting reference file path: " << (*sCPU).fittingPath << '\x0A';
	}

	fs << "Material name: " << crystalDatabasePtr[(*sCPU).materialIndex].crystalName << '\x0A';
	fs << "Sellmeier reference: " << crystalDatabasePtr[(*sCPU).materialIndex].sellmeierReference << '\x0A';
	fs << "Chi2 reference: " << crystalDatabasePtr[(*sCPU).materialIndex].dReference << '\x0A';
	fs << "Chi3 reference: " << crystalDatabasePtr[(*sCPU).materialIndex].chi3Reference << '\x0A';
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 22; ++j) {
			fs << crystalDatabasePtr[(*sCPU).materialIndex].sellmeierCoefficients[i * 22 + j];
			if (j < 21) fs << ',';
		}
		fs << '\x0A';
	}
	fs << "Code version: 0.71,  March 18, 2023";
	fs << '\x0A';

	return 0;
}


int loadPulseFiles(simulationParameterSet* sCPU) {

	//pulse type specifies if something has to be loaded to describe the pulses, or if they should be
	//synthesized later. 1: FROG .speck format; 2: EOS (not implemented yet)
	int frogLines = 0;
	int errCount = 0;
	if ((*sCPU).pulse1FileType == 1) {
		frogLines = loadFrogSpeck((*sCPU).field1FilePath, (*sCPU).loadedField1, (*sCPU).Ntime, (*sCPU).fStep, 0.0);
		if (frogLines > 1) {
			(*sCPU).field1IsAllocated = TRUE;
		}
		else {
			(*sCPU).field1IsAllocated = FALSE;
			errCount++;
		}
	}

	if ((*sCPU).pulse2FileType == 1) {
		frogLines = loadFrogSpeck((*sCPU).field2FilePath, (*sCPU).loadedField2, (*sCPU).Ntime, (*sCPU).fStep, 0.0);
		if (frogLines > 1) {
			(*sCPU).field2IsAllocated = TRUE;
		}
		else {
			(*sCPU).field2IsAllocated = FALSE;
			errCount++;
		}
	}
	return errCount;
}


int readInputParametersFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, const char* filePath) {
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
	fs >> (*sCPU).pulse1.energy;
	moveToColon();
	fs >> (*sCPU).pulse2.energy;
	moveToColon();
	fs >> (*sCPU).pulse1.frequency;
	moveToColon();
	fs >> (*sCPU).pulse2.frequency;
	moveToColon();
	fs >> (*sCPU).pulse1.bandwidth;
	moveToColon();
	fs >> (*sCPU).pulse2.bandwidth;
	moveToColon();
	fs >> (*sCPU).pulse1.sgOrder;
	moveToColon();
	fs >> (*sCPU).pulse2.sgOrder;
	moveToColon();
	fs >> (*sCPU).pulse1.cep;
	moveToColon();
	fs >> (*sCPU).pulse2.cep;
	moveToColon();
	fs >> (*sCPU).pulse1.delay;
	moveToColon();
	fs >> (*sCPU).pulse2.delay;
	moveToColon();
	fs >> (*sCPU).pulse1.gdd;
	moveToColon();
	fs >> (*sCPU).pulse2.gdd;
	moveToColon();
	fs >> (*sCPU).pulse1.tod;
	moveToColon();
	fs >> (*sCPU).pulse2.tod;
	moveToColon();
	fs >> (*sCPU).pulse1.phaseMaterial;
	moveToColon();
	fs >> (*sCPU).pulse2.phaseMaterial;
	moveToColon();
	fs >> (*sCPU).pulse1.phaseMaterialThickness;
	moveToColon();
	fs >> (*sCPU).pulse2.phaseMaterialThickness;
	moveToColon();
	moveToColon();
	fs >> (*sCPU).pulse1.beamwaist;
	moveToColon();
	fs >> (*sCPU).pulse2.beamwaist;
	moveToColon();
	fs >> (*sCPU).pulse1.x0;
	moveToColon();
	fs >> (*sCPU).pulse2.x0;
	moveToColon();
	fs >> (*sCPU).pulse1.y0;
	moveToColon();
	fs >> (*sCPU).pulse2.y0;
	moveToColon();
	fs >> (*sCPU).pulse1.z0;
	moveToColon();
	fs >> (*sCPU).pulse2.z0;
	moveToColon();
	fs >> (*sCPU).pulse1.beamAngle;
	moveToColon();
	fs >> (*sCPU).pulse2.beamAngle;
	moveToColon();
	fs >> (*sCPU).pulse1.beamAnglePhi;
	moveToColon();
	fs >> (*sCPU).pulse2.beamAnglePhi;
	moveToColon();
	fs >> (*sCPU).pulse1.polarizationAngle;
	moveToColon();
	fs >> (*sCPU).pulse2.polarizationAngle;
	moveToColon();
	fs >> (*sCPU).pulse1.circularity;
	moveToColon();
	fs >> (*sCPU).pulse2.circularity;
	moveToColon();
	fs >> (*sCPU).materialIndex;
	moveToColon();
	fs >> (*sCPU).materialIndexAlternate;
	moveToColon();
	fs >> (*sCPU).crystalTheta;
	moveToColon();
	fs >> (*sCPU).crystalPhi;
	moveToColon();
	fs >> (*sCPU).spatialWidth;
	moveToColon();
	fs >> (*sCPU).spatialHeight;
	moveToColon();
	fs >> (*sCPU).rStep;
	moveToColon();
	fs >> (*sCPU).timeSpan;
	moveToColon();
	fs >> (*sCPU).tStep;
	moveToColon();
	fs >> (*sCPU).crystalThickness;
	moveToColon();
	fs >> (*sCPU).propagationStep;
	moveToColon();
	fs >> (*sCPU).nonlinearAbsorptionStrength;
	moveToColon();
	fs >> (*sCPU).bandGapElectronVolts;
	moveToColon();
	fs >> (*sCPU).effectiveMass;
	moveToColon();
	fs >> (*sCPU).drudeGamma;
	moveToColon();
	fs >> (*sCPU).symmetryType;
	moveToColon();
	fs >> (*sCPU).batchIndex;
	moveToColon();
	fs >> (*sCPU).batchDestination;
	moveToColon();
	fs >> (*sCPU).Nsims;
	moveToColon();
	fs >> (*sCPU).batchIndex2;
	moveToColon();
	fs >> (*sCPU).batchDestination2;
	moveToColon();
	fs >> (*sCPU).Nsims2;

	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());
	memset((*sCPU).sequenceString, 0, 2 * MAX_LOADSTRING);
	line.copy((*sCPU).sequenceString, 2*MAX_LOADSTRING);
	
	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());
	memset((*sCPU).fittingString, 0, MAX_LOADSTRING);
	line.copy((*sCPU).fittingString, MAX_LOADSTRING);
	
	moveToColon();
	fs >> (*sCPU).fittingMode;

	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());
	line.copy((*sCPU).outputBasePath, MAX_LOADSTRING);

	moveToColon();
	fs >> (*sCPU).pulse1FileType;
	moveToColon();
	fs >> (*sCPU).pulse2FileType;

	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());
	line.copy((*sCPU).field1FilePath, MAX_LOADSTRING);

	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());
	line.copy((*sCPU).field2FilePath, MAX_LOADSTRING);

	moveToColon();
	std::getline(fs, line);
	line.erase(line.begin());
	line.copy((*sCPU).fittingPath, MAX_LOADSTRING);

	removeCharacterFromString((*sCPU).field1FilePath, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).field1FilePath, MAX_LOADSTRING, '\n');
	removeCharacterFromString((*sCPU).field2FilePath, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).field2FilePath, MAX_LOADSTRING, '\n');
	removeCharacterFromString((*sCPU).fittingPath, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).fittingPath, MAX_LOADSTRING, '\n');
	removeCharacterFromString((*sCPU).fittingString, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).fittingString, MAX_LOADSTRING, '\n');
	removeCharacterFromString((*sCPU).sequenceString, 2*MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).sequenceString, 2*MAX_LOADSTRING, '\n');
	removeCharacterFromString((*sCPU).outputBasePath, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).outputBasePath, MAX_LOADSTRING, '\n');

	//derived parameters and cleanup:
	(*sCPU).sellmeierType = 0;
	(*sCPU).axesNumber = 0;
	(*sCPU).Ntime = (size_t)(MIN_GRIDDIM * round((*sCPU).timeSpan / (MIN_GRIDDIM * (*sCPU).tStep)));
	(*sCPU).Nfreq = (*sCPU).Ntime / 2 + 1;
	(*sCPU).Nspace = (size_t)(MIN_GRIDDIM * round((*sCPU).spatialWidth / (MIN_GRIDDIM * (*sCPU).rStep)));
	(*sCPU).Nspace2 = (size_t)(MIN_GRIDDIM * round((*sCPU).spatialHeight / (MIN_GRIDDIM * (*sCPU).rStep)));
	(*sCPU).Ngrid = (*sCPU).Ntime * (*sCPU).Nspace;
	(*sCPU).NgridC = (*sCPU).Nfreq * (*sCPU).Nspace;
	(*sCPU).kStep = twoPi<double>() / ((*sCPU).Nspace * (*sCPU).rStep);
	(*sCPU).fStep = 1.0 / ((*sCPU).Ntime * (*sCPU).tStep);
	(*sCPU).Npropagation = (size_t)round((*sCPU).crystalThickness / (*sCPU).propagationStep);

	(*sCPU).isCylindric = (*sCPU).symmetryType == 1;
	(*sCPU).is3D = (*sCPU).symmetryType == 2;
	if ((*sCPU).isCylindric) {
		(*sCPU).pulse1.x0 = 0;
		(*sCPU).pulse2.x0 = 0;
		(*sCPU).pulse1.beamAngle = 0;
		(*sCPU).pulse2.beamAngle = 0;
	}
	if ((*sCPU).is3D) {
		(*sCPU).Ngrid = (*sCPU).Ntime * (*sCPU).Nspace * (*sCPU).Nspace2;
		(*sCPU).NgridC = (*sCPU).Nfreq * (*sCPU).Nspace * (*sCPU).Nspace2;
	}
	else {
		(*sCPU).Nspace2 = 1;
	}
	if ((*sCPU).batchIndex == 0 || (*sCPU).Nsims < 1) {
		(*sCPU).Nsims = 1;
	}
	if ((*sCPU).batchIndex2 == 0 || (*sCPU).Nsims2 < 1) {
		(*sCPU).Nsims2 = 1;
	}

	(*sCPU).field1IsAllocated = FALSE;
	(*sCPU).field2IsAllocated = FALSE;

	//crystal from database (database must be loaded!)
	(*sCPU).crystalDatabase = crystalDatabasePtr;
	(*sCPU).chi2Tensor = crystalDatabasePtr[(*sCPU).materialIndex].d;
	(*sCPU).chi3Tensor = crystalDatabasePtr[(*sCPU).materialIndex].chi3;
	(*sCPU).nonlinearSwitches = crystalDatabasePtr[(*sCPU).materialIndex].nonlinearSwitches;
	(*sCPU).absorptionParameters = crystalDatabasePtr[(*sCPU).materialIndex].absorptionParameters;
	(*sCPU).sellmeierCoefficients = crystalDatabasePtr[(*sCPU).materialIndex].sellmeierCoefficients;
	(*sCPU).sellmeierType = crystalDatabasePtr[(*sCPU).materialIndex].sellmeierType;
	(*sCPU).axesNumber = crystalDatabasePtr[(*sCPU).materialIndex].axisType;

	if (fs.good())return 61;
	else return -1;
}

int saveDataSet(simulationParameterSet* sCPU) {
	saveSettingsFile(sCPU);

	std::string Epath((*sCPU).outputBasePath);
	Epath.append("_Ext.dat");
	std::ofstream Efile(Epath, std::ios::binary);
	if (Efile.is_open()) Efile.write(reinterpret_cast<char*>((*sCPU).ExtOut), 2 * ((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2) * sizeof(double));

	std::string Spath((*sCPU).outputBasePath);
	Spath.append("_spectrum.dat");
	std::ofstream Sfile(Spath, std::ios::binary);
	if (Sfile.is_open()) Sfile.write(reinterpret_cast<char*>((*sCPU).totalSpectrum), (*sCPU).Nsims * (*sCPU).Nsims2 * 3 * (*sCPU).Nfreq * sizeof(double));
	
	return 0;
}


int configureBatchMode(simulationParameterSet* sCPU) {
	size_t i, j, currentRow;
	if ((*sCPU).batchIndex == 0 || (*sCPU).Nsims == 1 || (*sCPU).batchIndex > 37 || (*sCPU).batchIndex2 > 37) {
		return 0;
	}

	//pointers to values that can be scanned in batch mode
	double* targets[38] = { 0,
		&(*sCPU).pulse1.energy, &(*sCPU).pulse2.energy, &(*sCPU).pulse1.frequency, &(*sCPU).pulse2.frequency,
		&(*sCPU).pulse1.bandwidth, &(*sCPU).pulse2.bandwidth, &(*sCPU).pulse1.cep, &(*sCPU).pulse2.cep,
		&(*sCPU).pulse1.delay, &(*sCPU).pulse2.delay, &(*sCPU).pulse1.gdd, &(*sCPU).pulse2.gdd,
		&(*sCPU).pulse1.tod, &(*sCPU).pulse2.tod, &(*sCPU).pulse1.phaseMaterialThickness, &(*sCPU).pulse2.phaseMaterialThickness,
		&(*sCPU).pulse1.beamwaist, &(*sCPU).pulse2.beamwaist,
		&(*sCPU).pulse1.x0, &(*sCPU).pulse2.x0, &(*sCPU).pulse1.z0, &(*sCPU).pulse2.z0,
		&(*sCPU).pulse1.beamAngle, &(*sCPU).pulse2.beamAngle, &(*sCPU).pulse1.polarizationAngle, &(*sCPU).pulse2.polarizationAngle,
		&(*sCPU).pulse1.circularity, &(*sCPU).pulse2.circularity, &(*sCPU).crystalTheta, &(*sCPU).crystalPhi,
		&(*sCPU).nonlinearAbsorptionStrength, &(*sCPU).drudeGamma, &(*sCPU).effectiveMass, &(*sCPU).crystalThickness,
		&(*sCPU).propagationStep, &(*sCPU).i37, &(*sCPU).i37};

	//multipliers to the Batch end value from the interface
	// (e.g. frequency in THz requires 1e12 multiplier)
	double multipliers[38] = { 0,
		1, 1, 1e12, 1e12,
		1e12, 1e12, vPi<double>(), vPi<double>(),
		1e-15, 1e-15, 1e-30, 1e-30,
		1e-45, 1e-45, 1e-6, 1e-6,
		1e-6, 1e-6,
		1e-6, 1e-6, 1e-6, 1e-6,
		deg2Rad<double>(), deg2Rad<double>(), deg2Rad<double>(), deg2Rad<double>(),
		1, 1, deg2Rad<double>(), deg2Rad<double>(),
		1, 1e12, 1, 1e-6,
		1e-9, 1.0, 1.0 };

	//Configure the struct array if in a batch
	for (i = 0; i < (*sCPU).Nsims2; i++) {
		currentRow = i * (*sCPU).Nsims;
		
		if (currentRow > 0) {
			memcpy(&sCPU[currentRow], sCPU, sizeof(simulationParameterSet));
		}
		if ((*sCPU).Nsims2 > 1) {
			*((double*)((simulationParameterSet*)targets[(*sCPU).batchIndex2] + currentRow)) +=
				i * (multipliers[(*sCPU).batchIndex2] * (*sCPU).batchDestination2 - *targets[(*sCPU).batchIndex2])
				/ ((*sCPU).Nsims2 - 1);
		}

		for (j = 0; j < (*sCPU).Nsims; j++) {

			if (j > 0) {
				memcpy(&sCPU[j + currentRow], &sCPU[currentRow], sizeof(simulationParameterSet));
			}

			if ((*sCPU).deffTensor != NULL) {
				sCPU[j + currentRow].deffTensor = &(*sCPU).deffTensor[9 * (j + currentRow)];;
			}
			sCPU[j + currentRow].batchLoc1 = j;
			sCPU[j + currentRow].batchLoc2 = i;
			sCPU[j + currentRow].ExtOut = &(*sCPU).ExtOut[(j + currentRow) * (*sCPU).Ngrid * 2];
			sCPU[j + currentRow].EkwOut = &(*sCPU).EkwOut[(j + currentRow) * (*sCPU).NgridC * 2];
			sCPU[j + currentRow].totalSpectrum = &(*sCPU).totalSpectrum[(j + currentRow) * (*sCPU).Nfreq * 3];
			sCPU[j + currentRow].isFollowerInSequence = FALSE;

			// To add new modes, append values to the two arrays above, and to the combobox in the UI.
			// Cast the pointer to the original value to a pointer to a struct, 
			// increment, recast to a pointer to double and resolve then add j times the scan step size.
			*((double*)((simulationParameterSet*)targets[(*sCPU).batchIndex] + (j + currentRow))) +=
				j * (multipliers[(*sCPU).batchIndex] * (*sCPU).batchDestination - *targets[(*sCPU).batchIndex])
				/ ((*sCPU).Nsims - 1);
		}
	}

	return 0;
}


int allocateGrids(simulationParameterSet* sCPU) {
	(*sCPU).loadedField1 = new std::complex<double>[(*sCPU).Ntime]();
	(*sCPU).loadedField2 = new std::complex<double>[(*sCPU).Ntime]();
	(*sCPU).ExtOut = new double[(*sCPU).Ngrid * 2 * (*sCPU).Nsims * (*sCPU).Nsims2]();
	(*sCPU).EkwOut = new std::complex<double>[(*sCPU).NgridC * 2 * (*sCPU).Nsims * (*sCPU).Nsims2]();
	(*sCPU).deffTensor = new double[9 * (*sCPU).Nsims * (*sCPU).Nsims2]();
	(*sCPU).totalSpectrum = new double[(*sCPU).Nsims * (*sCPU).Nsims2 * (*sCPU).Nfreq * 3]();
	(*sCPU).statusFlags = new int[(*sCPU).Nsims * (*sCPU).Nsims2]();
	(*sCPU).fittingReference = new double[(*sCPU).Nfreq]();
	return 0;
}

int deallocateGrids(simulationParameterSet* sCPU, bool alsoDeleteDisplayItems) {
	delete[] (*sCPU).loadedField1;
	delete[] (*sCPU).loadedField2;
	delete[] (*sCPU).deffTensor;
	delete[] (*sCPU).statusFlags;
	delete[] (*sCPU).fittingReference;
	if (alsoDeleteDisplayItems) {
		delete[](*sCPU).ExtOut;
		delete[](*sCPU).EkwOut;
		delete[](*sCPU).totalSpectrum;
	}
	return 0;
}

//calculates the squard modulus of a complex number, under the assumption that the
//machine's complex number format is interleaved doubles.
//c forced to run in c++ for nostalgia reasons
double cModulusSquared(std::complex<double>complexNumber) {
	double* xy = (double*)&complexNumber;
	return xy[0] * xy[0] + xy[1] * xy[1];
}


int loadFrogSpeck(char* frogFilePath, std::complex<double>* Egrid, long long Ntime, double fStep, double gateLevel) {
	std::string line;
	std::ifstream fs(frogFilePath);
	if (fs.fail()) return -1;
	int maxFileSize = 16384;
	double wavelength, R, phi, complexX, complexY, f, f0, f1;
	double fmax = 0.0;
	int i, k0, k1;
	double c = 1e9 * lightC<double>(); //for conversion of wavelength in nm to frequency
	double df = 0;
	double fmin = 0;
	int currentRow = 0;
	std::complex<double>* E = new std::complex<double>[maxFileSize]();
	if (E == NULL) {
		return -2;
	}

	while (fs.good() && currentRow < maxFileSize) {
		fs >> wavelength;
		fs >> R;
		fs >> phi;
		fs >> complexX;
		fs >> complexY;
		std::getline(fs, line);

		//get the complex field from the data
		E[currentRow].real(complexX);
		E[currentRow].imag(complexY);

		//keep track of the frequency step of the grid (running sum, divide by number of rows at end to get average)
		if (currentRow > 0) df += c / wavelength - fmax;

		//keep track of the highest frequency in the data
		fmax = c / wavelength;

		//store the lowest frequency in the data
		if (currentRow == 0) fmin = fmax;

		currentRow++;
	}

	//return an error if nothing was loaded
	if (currentRow == 0) {
		delete[] E;
		return -1;
	}

	df /= currentRow; //average frequency step

	//interpolate the FROG data onto the simulation grid

	//fill the simulation grid based on the data
	for (i = 0; i < Ntime / 2 + 1; i++) {

		//frequency grid used in the simulation
		f = i * fStep;

		k0 = (int)floor((f - fmin) / df);
		k1 = (int)ceil((f - fmin) / df);
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

	delete[] E;
	return currentRow;
}
