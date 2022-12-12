#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex>
#include "LightwaveExplorerUtilities.h"



int readFittingString(simulationParameterSet* sCPU) {
	//read the fitting string (if there is one), convert it into an array if it exists
	char fittingString[MAX_LOADSTRING];
	double ROIbegin;
	double ROIend;
	strcpy_s(fittingString, MAX_LOADSTRING, (*sCPU).fittingString);
	char* nextToken = NULL;
	char* tokToken = strtok_s(fittingString, ";", &nextToken);
	bool paramsRead = (3 == sscanf_s(fittingString, "%lf %lf %d",
		&ROIbegin, &ROIend, &(*sCPU).fittingMaxIterations));
	(*sCPU).fittingROIstart = (size_t)(ROIbegin / (*sCPU).fStep);
	(*sCPU).fittingROIstop = (size_t)minN(ROIend / (*sCPU).fStep, (*sCPU).Ntime / 2);
	(*sCPU).fittingROIsize = minN(maxN(1, (*sCPU).fittingROIstop - (*sCPU).fittingROIstart), (*sCPU).Ntime / 2);
	int fittingCount = 0;
	tokToken = strtok_s(NULL, ";", &nextToken);
	int lastread = 3;
	while (tokToken != NULL && lastread == 3) {
		lastread = sscanf_s(tokToken, "%lf %lf %lf",
			&(*sCPU).fittingArray[fittingCount], &(*sCPU).fittingArray[fittingCount + 1],
			&(*sCPU).fittingArray[fittingCount + 2]);
		if (lastread > 0) {
			fittingCount += lastread;
		}
		tokToken = strtok_s(NULL, ";", &nextToken);
	}
	(*sCPU).Nfitting = fittingCount / 3;
	(*sCPU).isInFittingMode = (((*sCPU).Nfitting) > 0 && paramsRead);

	if (!(*sCPU).isInFittingMode) {
		char nopeString[] = "None.";
		strcpy_s((*sCPU).fittingString, MAX_LOADSTRING, nopeString);
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

int removeCharacterFromStringSkippingAngleBrackets(char* cString, size_t N, char removedChar) {
	size_t i = 0;
	size_t r = 0;
	while (i < N - 1) {
		if (cString[i] == removedChar) {
			memmove(&cString[i], &cString[i + 1], N - i - r - 1);
			cString[N - r - 1] = 0;
			r++;
		}
		else if (cString[i] == '<') {
			i += findClosingAngleBracket(&cString[i]) - &cString[i];
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

void stripWhiteSpace(char* sequenceString) {
	removeCharacterFromStringSkippingAngleBrackets(sequenceString, strlen(sequenceString), ' ');
	removeCharacterFromString(sequenceString, strlen(sequenceString), '\n');
	removeCharacterFromString(sequenceString, strlen(sequenceString), '\t');
	removeCharacterFromString(sequenceString, strlen(sequenceString), '\r');
}

void stripLineBreaks(char* sequenceString) {
	removeCharacterFromString(sequenceString, strlen(sequenceString), '\n');
	removeCharacterFromString(sequenceString, strlen(sequenceString), '\r');
}

char* findClosingParenthesis(const char* s) {
	char* b = (char*)s;
	for (;;) {
		if (b[0] == 0) return NULL;
		if (b[0] == ')') {
			if (b[1] == ';') return b + 1;
			return b;
		}
		else {
			b++;
		}
	}
}

char* findClosingCurlyBracket(const char* s) {
	char* b = (char*)s;
	for (;;) {
		if (b[0] == 0) return NULL;
		if (b[0] == '}') {
			if (b[1] == ';') return b + 1;
			return b;
		}
		else {
			b++;
		}
	}
}

char* findClosingAngleBracket(const char* s) {
	char* b = (char*)s;
	for (;;) {
		if (b[0] == 0) return NULL;
		if (b[0] == '>') {
			if (b[1] == ';') return b + 1;
			return b;
		}
		else {
			b++;
		}
	}
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
		*result = pow(*result, *readout);
		return;
	}
}

double parameterStringToDouble(const char* pString, double* iBlock, double* vBlock) {
	int loc = 0;
	double result = 0.0;
	double readout = 0.0;
	int ind = 0;
	bool previousCharWasOp = 0;
	char lastOp = 0;
	while (loc < 255) {
		if (pString[loc] == 0) return result;

		if (!previousCharWasOp) {
			if (pString[loc] == 'v') {
				++loc;
				sscanf_s(&pString[loc], "%d", &ind);
				loc += 2;
				if (ind < 100) result = vBlock[ind];
			}
			else if (pString[loc] == 'i') {
				++loc;
				sscanf_s(&pString[loc], "%d", &ind);
				loc += 2;
				if (ind < 100) result = iBlock[ind];
			}

			else if (pString[loc] == '*'
				|| pString[loc] == '-'
				|| pString[loc] == '+'
				|| pString[loc] == '/'
				|| pString[loc] == '^') {
				previousCharWasOp = 1;
				lastOp = pString[loc];
				loc++;
			}
			else {
				sscanf_s(&pString[loc], "%lf", &result);
				while (!(pString[loc] == '*'
					|| pString[loc] == '-'
					|| pString[loc] == '+'
					|| pString[loc] == '/'
					|| pString[loc] == '^')) {
					if (pString[loc] == 0) {
						return result;
					}
					if (pString[loc] == 'e') {
						if (pString[loc + 1] == '-' || pString[loc + 1] == '+') loc += 3;
					}
					else {
						loc++;
					}
				}
			}
		}
		else {
			if (pString[loc] == 'v') {
				++loc;
				sscanf_s(&pString[loc], "%d", &ind);
				loc += 2;
				if (ind < 100)readout = vBlock[ind];
				applyOp(lastOp, &result, &readout);
				previousCharWasOp = 0;
			}
			else if (pString[loc] == 'i') {
				++loc;
				sscanf_s(&pString[loc], "%d", &ind);
				loc += 2;
				if (ind < 100) readout = iBlock[ind];
				applyOp(lastOp, &result, &readout);
				previousCharWasOp = 0;
			}
			else {
				sscanf_s(&pString[loc], "%lf", &readout);
				applyOp(lastOp, &result, &readout);
				previousCharWasOp = 0;
				while (!(pString[loc] == '*'
					|| pString[loc] == '-'
					|| pString[loc] == '+'
					|| pString[loc] == '/'
					|| pString[loc] == '^')) {
					if (pString[loc] == 0) {
						return result;
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
	
	FILE* fp;// = fopen(spectrumPath, "r");
	if (fopen_s(&fp, spectrumPath, "r")) {
		printf("Could not read reference file\r\n");
		return 1;
	}
	size_t maxFileSize = 16384;
	size_t currentRow = 0;
	double c = 1e9 * LIGHTC;
	double* loadedWavelengths = (double*)calloc(8192*3, sizeof(double));
	double* loadedFrequencies = loadedWavelengths + 8192;
	double* loadedIntensities = loadedWavelengths + 16384;;
	if (loadedWavelengths == NULL) {
		return 1;
	}
	double maxWavelength = 0;
	double minWavelength = 0;

	while (fscanf_s(fp, "%lf %lf", &loadedWavelengths[currentRow], &loadedIntensities[currentRow]) == 2 && currentRow < maxFileSize) {
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
	//memset((*sCPU).fittingArray, 0, (*sCPU).Nfreq * sizeof(double));
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
	fclose(fp);
	free(loadedWavelengths);

	return 0;
}



int loadSavedFields(simulationParameterSet* sCPU, char* outputBase) {
	char outputpath[MAX_LOADSTRING] = { 0 };
	size_t writeSize = 2 * ((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2);



	//read fields as binary
	FILE* ExtOutFile;
	strcpy_s(outputpath, MAX_LOADSTRING, outputBase);
	strcat_s(outputpath, MAX_LOADSTRING, "_Ext.dat");
	//ExtOutFile = fopen(outputpath, "rb");
	if (fopen_s(&ExtOutFile,outputpath,"rb")) return 1;
	
	fread_s((*sCPU).ExtOut, writeSize*sizeof(double), sizeof(double), writeSize, ExtOutFile);
	fclose(ExtOutFile);

	FILE* spectrumFile;
	strcpy_s(outputpath, MAX_LOADSTRING, outputBase);
	strcat_s(outputpath, MAX_LOADSTRING, "_spectrum.dat");
	if (fopen_s(&spectrumFile, outputpath, "rb")) return 1;
	fread_s((*sCPU).totalSpectrum, (*sCPU).Nsims * (*sCPU).Nsims2 * 3 * (*sCPU).Nfreq * sizeof(double), sizeof(double), (*sCPU).Nsims * (*sCPU).Nsims2 * 3 * (*sCPU).Nfreq, spectrumFile);
	fclose(spectrumFile);

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


int saveSlurmScript(simulationParameterSet* sCPU, int gpuType, int gpuCount) {
	FILE* textfile;
	char outputpath[MAX_LOADSTRING] = { 0 };

	char* fileName = (*sCPU).outputBasePath;
	while (strchr(fileName, '\\') != NULL) {
		fileName = strchr(fileName, '\\');
		fileName++;
	}

	strcpy_s(outputpath, MAX_LOADSTRING, (*sCPU).outputBasePath);
	strcat_s(outputpath, MAX_LOADSTRING, ".slurmScript");
	if (fopen_s(&textfile, outputpath, "wb")) return 1;
	fprintf(textfile, "#!/bin/bash -l"); unixNewLine(textfile);
	fprintf(textfile, "#SBATCH -o ./tjob.out.%%j"); unixNewLine(textfile);
	fprintf(textfile, "#SBATCH -e ./tjob.err.%%j"); unixNewLine(textfile);
	fprintf(textfile, "#SBATCH -D ./"); unixNewLine(textfile);
	fprintf(textfile, "#SBATCH -J lightwave");  unixNewLine(textfile);
	fprintf(textfile, "#SBATCH --constraint=\"gpu\""); unixNewLine(textfile);
	if (gpuType == 0) {
		fprintf(textfile, "#SBATCH --gres=gpu:rtx5000:%i", minN(gpuCount, 2)); unixNewLine(textfile);
	}
	if (gpuType == 1) {
		fprintf(textfile, "#SBATCH --gres=gpu:v100:%i", minN(gpuCount, 2)); unixNewLine(textfile);
	}
	if (gpuType == 2) {
		fprintf(textfile, "#SBATCH --gres=gpu:a100:%i", minN(gpuCount, 4)); unixNewLine(textfile);
		fprintf(textfile, "#SBATCH --cpus-per-task=%i", 2 * minN(gpuCount, 4)); unixNewLine(textfile);
	}
	fprintf(textfile, "#SBATCH --mem=%zuM", 8192 + (18 * sizeof(double) * (*sCPU).Ngrid * maxN(1, (*sCPU).Nsims)) / 1048576);
	unixNewLine(textfile);
	fprintf(textfile, "#SBATCH --nodes=1"); unixNewLine(textfile);
	fprintf(textfile, "#SBATCH --ntasks-per-node=1"); unixNewLine(textfile);
	fprintf(textfile, "#SBATCH --time=03:00:00"); unixNewLine(textfile);
	fprintf(textfile, "module purge"); unixNewLine(textfile);
	fprintf(textfile, "module load cuda/11.6"); unixNewLine(textfile);
	fprintf(textfile, "module load mkl/2022.1"); unixNewLine(textfile);
	fprintf(textfile, "export LD_LIBRARY_PATH=$MKL_HOME/lib/intel64:$LD_LIBRARY_PATH"); unixNewLine(textfile);
	if (gpuType == 0 || gpuType == 1) {
		fprintf(textfile, "srun ./lwe %s.input > %s.out", fileName, fileName); unixNewLine(textfile);
	}
	if (gpuType == 2) {
		fprintf(textfile, "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}"); unixNewLine(textfile);
		fprintf(textfile, "srun ./lwe %s.input > %s.out", fileName, fileName); unixNewLine(textfile);
	}
	fclose(textfile);
	return 0;
}

//print a linefeed without a carriage return so that linux systems don't complain
//about impure scripts from DOS machines
//fopen() should be called with "wb"
void unixNewLine(FILE* iostream) {
	char LF = '\x0A';
	fwrite(&LF, sizeof(char), 1, iostream);
}


int saveSettingsFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr) {
	int j, k;
	FILE* textfile;
	wchar_t wideStringConversionBuffer[MAX_LOADSTRING] = { 0 };
	char outputpath[MAX_LOADSTRING] = { 0 };
	strcpy_s(outputpath, MAX_LOADSTRING, (*sCPU).outputBasePath);
	if ((*sCPU).runType > 0) {
		strcat_s(outputpath, MAX_LOADSTRING, ".input");
	}
	else {
		strcat_s(outputpath, MAX_LOADSTRING, ".txt");
	}

	if(fopen_s(&textfile, outputpath, "w")) return 1;
	fwprintf(textfile, L"Pulse energy 1 (J): %14.14e\nPulse energy 2 (J): %14.14e\nFrequency 1 (Hz): %14.14e\n",
		(*sCPU).pulse1.energy, (*sCPU).pulse2.energy, (*sCPU).pulse1.frequency);
	fwprintf(textfile, L"Frequency 2 (Hz): %14.14e\nBandwidth 1 (Hz): %14.14e\nBandwidth 2 (Hz): %14.14e\n",
		(*sCPU).pulse2.frequency, (*sCPU).pulse1.bandwidth, (*sCPU).pulse2.bandwidth);
	fwprintf(textfile, L"SG order 1: %i\nSG order 2: %i\nCEP 1 (rad): %14.14e\nCEP 2 (rad): %14.14e\nDelay 1 (s): %14.14e\nDelay 2 (s): %14.14e\nGDD 1 (s^-2): %14.14e\nGDD 2 (s^-2): %14.14e\nTOD 1 (s^-3): %14.14e\nTOD 2 (s^-3): %14.14e\n",
		(*sCPU).pulse1.sgOrder, (*sCPU).pulse2.sgOrder, (*sCPU).pulse1.cep, (*sCPU).pulse2.cep, (*sCPU).pulse1.delay, (*sCPU).pulse2.delay, (*sCPU).pulse1.gdd, (*sCPU).pulse2.gdd, (*sCPU).pulse1.tod, (*sCPU).pulse2.tod);
	fwprintf(textfile, L"Phase material 1 index: %i\nPhase material 2 index: %i\nPhase material thickness 1 (mcr.): %14.14e\nPhase material thickness 2 (mcr.): %14.14e\n",
		(*sCPU).pulse1.phaseMaterial, (*sCPU).pulse2.phaseMaterial, (*sCPU).pulse1.phaseMaterialThickness, (*sCPU).pulse2.phaseMaterialThickness);
	fwprintf(textfile, L"Beam mode placeholder: 0\n");
	fwprintf(textfile, L"Beamwaist 1 (m): %14.14e\nBeamwaist 2 (m): %14.14e\nx offset 1 (m): %14.14e\nx offset 2 (m): %14.14e\n",
		(*sCPU).pulse1.beamwaist, (*sCPU).pulse2.beamwaist, (*sCPU).pulse1.x0, (*sCPU).pulse2.x0);
	fwprintf(textfile, L"y offset 1 (m): %14.14e\ny offset 2 (m): %14.14e\n",
		(*sCPU).pulse1.y0, (*sCPU).pulse2.y0);
	fwprintf(textfile, L"z offset 1 (m): %14.14e\nz offset 2 (m): %14.14e\n",
		(*sCPU).pulse1.z0, (*sCPU).pulse2.z0);
	fwprintf(textfile, L"NC angle 1 (rad): %14.14e\nNC angle 2 (rad): %14.14e\n",
		(*sCPU).pulse1.beamAngle, (*sCPU).pulse2.beamAngle);
	fwprintf(textfile, L"NC angle phi 1 (rad): %14.14e\nNC angle phi 2 (rad): %14.14e\n",
		(*sCPU).pulse1.beamAnglePhi, (*sCPU).pulse2.beamAnglePhi);
	fwprintf(textfile, L"Polarization 1 (rad): %14.14e\nPolarization 2 (rad): %14.14e\nCircularity 1: %14.14e\nCircularity 2: %14.14e\n",
		(*sCPU).pulse1.polarizationAngle, (*sCPU).pulse2.polarizationAngle, (*sCPU).pulse1.circularity, (*sCPU).pulse2.circularity);
	fwprintf(textfile, L"Material index: %i\nAlternate material index: %i\n",
		(*sCPU).materialIndex, (*sCPU).materialIndexAlternate);
	fwprintf(textfile, L"Crystal theta (rad): %14.14e\nCrystal phi (rad): %14.14e\nGrid width (m): %14.14e\nGrid height (m): %14.14e\ndx (m): %14.14e\nTime span (s): %14.14e\ndt (s): %14.14e\nThickness (m): %14.14e\ndz (m): %14.14e\n",
		(*sCPU).crystalTheta, (*sCPU).crystalPhi, (*sCPU).spatialWidth, (*sCPU).spatialHeight, (*sCPU).rStep, (*sCPU).timeSpan, (*sCPU).tStep, (*sCPU).crystalThickness, (*sCPU).propagationStep);
	fwprintf(textfile, L"Nonlinear absorption parameter: %14.14e\nBand gap (eV): %14.14e\nEffective mass (relative): %14.14e\nDrude gamma (Hz): %14.14e\n",
		(*sCPU).nonlinearAbsorptionStrength, (*sCPU).bandGapElectronVolts, (*sCPU).effectiveMass, (*sCPU).drudeGamma);
	fwprintf(textfile, L"Propagation mode: %i\n",
		(*sCPU).symmetryType);
	fwprintf(textfile, L"Batch mode: %i\nBatch destination: %14.14e\nBatch steps: %lli\n",
		(*sCPU).batchIndex, (*sCPU).batchDestination, (*sCPU).Nsims);
	fwprintf(textfile, L"Batch mode 2: %i\nBatch destination 2: %14.14e\nBatch steps 2: %lli\n",
		(*sCPU).batchIndex2, (*sCPU).batchDestination2, (*sCPU).Nsims2);
	size_t convertedCount;
	mbstowcs_s(&convertedCount, wideStringConversionBuffer, (*sCPU).sequenceString, MAX_LOADSTRING);
	fwprintf(textfile, L"Sequence: %ls\n", wideStringConversionBuffer);
	mbstowcs_s(&convertedCount, wideStringConversionBuffer, (*sCPU).fittingString, MAX_LOADSTRING);
	fwprintf(textfile, L"Fitting: %ls\n", wideStringConversionBuffer);
	fwprintf(textfile, L"Fitting mode: %i\n", (*sCPU).fittingMode);
	

	if ((*sCPU).runType > 0) {
		char* fileName = (*sCPU).outputBasePath;
		while (strchr(fileName, '\\') != NULL) {
			fileName = strchr(fileName, '\\');
			fileName++;
		}
		mbstowcs_s(&convertedCount, wideStringConversionBuffer, fileName, MAX_LOADSTRING);
		wideStringConversionBuffer[strlen(fileName)] = L'\0';
		fwprintf(textfile, L"Output base path: %ls\n", wideStringConversionBuffer);
	}
	else {
		mbstowcs_s(&convertedCount, wideStringConversionBuffer, (*sCPU).outputBasePath, MAX_LOADSTRING);
		fwprintf(textfile, L"Output base path: %ls\n", wideStringConversionBuffer);
	}

	fwprintf(textfile, L"Field 1 from file type: %i\nField 2 from file type: %i\n", (*sCPU).pulse1FileType, (*sCPU).pulse2FileType);
	mbstowcs_s(&convertedCount, wideStringConversionBuffer, (*sCPU).field1FilePath, MAX_LOADSTRING);
	fwprintf(textfile, L"Field 1 file path: %ls\n", wideStringConversionBuffer);
	mbstowcs_s(&convertedCount, wideStringConversionBuffer, (*sCPU).field2FilePath, MAX_LOADSTRING);
	fwprintf(textfile, L"Field 2 file path: %ls\n", wideStringConversionBuffer);
	mbstowcs_s(&convertedCount, wideStringConversionBuffer, (*sCPU).fittingPath, MAX_LOADSTRING);
	fwprintf(textfile, L"Fitting reference file path: %ls\n", wideStringConversionBuffer);

	fwprintf(textfile, L"Material name: %ls\nSellmeier reference: %ls\nChi2 reference: %ls\nChi3 reference: %ls\n", crystalDatabasePtr[(*sCPU).materialIndex].crystalNameW, crystalDatabasePtr[(*sCPU).materialIndex].sellmeierReference, crystalDatabasePtr[(*sCPU).materialIndex].dReference, crystalDatabasePtr[(*sCPU).materialIndex].chi3Reference);
	fwprintf(textfile, L"Sellmeier coefficients: \n");
	for (j = 0; j < 3; j++) {
		for (k = 0; k < 22; k++) {
			fwprintf(textfile, L"%14.14e ", crystalDatabasePtr[(*sCPU).materialIndex].sellmeierCoefficients[j * 22 + k]);
		}
		fwprintf(textfile, L"\n");
	}
	fwprintf(textfile, L"Code version: 0.5 September 16, 2022\n");

	fclose(textfile);
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


int readInputParametersFile(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, char* filePath) {
	FILE* textfile;
	if (fopen_s(&textfile, filePath, "r")) {
		return 1;
	}
	//read parameters using fscanf:
	//recipe for programming: copy/paste the block of fprintf statements in the saveDataSet() function,
	//then find/replace:
	// fprintf->fscanf
	// (*CPU). -> &(*CPU).
	// %e -> %lf
	// &(*sCPU).sequenceString -> (*sCPU).sequenceString
	// &(*sCPU).outputBasePath -> (*sCPU).outputBasePath
	int readValueCount = 0;
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.energy);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.energy);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.frequency);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.frequency);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.bandwidth);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.bandwidth);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%i", &(*sCPU).pulse1.sgOrder);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%i", &(*sCPU).pulse2.sgOrder);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.cep);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.cep);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.delay);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.delay);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.gdd);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.gdd);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.tod);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.tod);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%d", &(*sCPU).pulse1.phaseMaterial);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%d", &(*sCPU).pulse2.phaseMaterial);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.phaseMaterialThickness);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.phaseMaterialThickness);
	skipFileUntilCharacter(textfile, ':');
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.beamwaist);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.beamwaist);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.x0);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.x0);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.y0);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.y0);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.z0);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.z0);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.beamAngle);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.beamAngle);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.beamAnglePhi);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.beamAnglePhi);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.polarizationAngle);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.polarizationAngle);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse1.circularity);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).pulse2.circularity);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%i", &(*sCPU).materialIndex);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%i", &(*sCPU).materialIndexAlternate);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).crystalTheta);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).crystalPhi);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).spatialWidth);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).spatialHeight);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).rStep);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).timeSpan);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).tStep);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).crystalThickness);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).propagationStep);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).nonlinearAbsorptionStrength);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).bandGapElectronVolts);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).effectiveMass);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).drudeGamma);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%i", &(*sCPU).symmetryType);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%i", &(*sCPU).batchIndex);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).batchDestination);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%zu", &(*sCPU).Nsims);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%i", &(*sCPU).batchIndex2);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%lf", &(*sCPU).batchDestination2);
	skipFileUntilCharacter(textfile, ':');
	readValueCount += fscanf_s(textfile, "%zu", &(*sCPU).Nsims2);
	readValueCount += fscanf_s(textfile, "\nSequence: ");
	fgets((*sCPU).sequenceString, MAX_LOADSTRING, textfile);
	readValueCount += fscanf_s(textfile, "Fitting: ");
	fgets((*sCPU).fittingString, MAX_LOADSTRING, textfile);
	readValueCount += fscanf_s(textfile, "Fitting mode : %i\n", &(*sCPU).fittingMode);
	readValueCount += fscanf_s(textfile, "Output base path: ");
	fgets((*sCPU).outputBasePath, MAX_LOADSTRING, textfile);
	readValueCount += fscanf_s(textfile, "Field 1 from file type: %i\nField 2 from file type: %i\n",
		&(*sCPU).pulse1FileType, &(*sCPU).pulse2FileType);
	readValueCount += fscanf_s(textfile, "Field 1 file path: ");
	fgets((*sCPU).field1FilePath, MAX_LOADSTRING, textfile);
	readValueCount += fscanf_s(textfile, "Field 2 file path: ");
	fgets((*sCPU).field2FilePath, MAX_LOADSTRING, textfile);
	readValueCount += fscanf_s(textfile, "Fitting reference file path: ");
	fgets((*sCPU).fittingPath, MAX_LOADSTRING, textfile);

	removeCharacterFromString((*sCPU).field1FilePath, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).field1FilePath, MAX_LOADSTRING, '\n');
	removeCharacterFromString((*sCPU).field2FilePath, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).field2FilePath, MAX_LOADSTRING, '\n');
	removeCharacterFromString((*sCPU).fittingPath, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).fittingPath, MAX_LOADSTRING, '\n');
	removeCharacterFromString((*sCPU).fittingString, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).fittingString, MAX_LOADSTRING, '\n');
	removeCharacterFromString((*sCPU).sequenceString, MAX_LOADSTRING, '\r');
	removeCharacterFromString((*sCPU).sequenceString, MAX_LOADSTRING, '\n');
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
	(*sCPU).kStep = TWOPI / ((*sCPU).Nspace * (*sCPU).rStep);
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
	(*sCPU).chi2Tensor = crystalDatabasePtr[(*sCPU).materialIndex].d;
	(*sCPU).chi3Tensor = crystalDatabasePtr[(*sCPU).materialIndex].chi3;
	(*sCPU).nonlinearSwitches = crystalDatabasePtr[(*sCPU).materialIndex].nonlinearSwitches;
	(*sCPU).absorptionParameters = crystalDatabasePtr[(*sCPU).materialIndex].absorptionParameters;
	(*sCPU).sellmeierCoefficients = crystalDatabasePtr[(*sCPU).materialIndex].sellmeierCoefficients;
	(*sCPU).sellmeierType = crystalDatabasePtr[(*sCPU).materialIndex].sellmeierType;
	(*sCPU).axesNumber = crystalDatabasePtr[(*sCPU).materialIndex].axisType;

	fclose(textfile);
	return readValueCount;
}

int saveDataSet(simulationParameterSet* sCPU, crystalEntry* crystalDatabasePtr, char* outputbase, bool saveInputs) {

	saveSettingsFile(sCPU, crystalDatabasePtr);
	char outputpath[MAX_LOADSTRING] = { 0 };
	double matlabpadding[1024] = { 0 };

	char* outputbaseVar = strrchr(outputbase, '\\');
	if (!outputbaseVar) {
		outputbaseVar = outputbase;
	}
	else {
		outputbaseVar++;
	}
	
	//write field as binary
	FILE* ExtOutFile;
	size_t writeSize = 2 * ((*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2);
	strcpy_s(outputpath, MAX_LOADSTRING, outputbase);
	strcat_s(outputpath, MAX_LOADSTRING, "_Ext.dat");
	if (fopen_s(&ExtOutFile, outputpath, "wb")) return 1;
	fwrite((*sCPU).ExtOut, sizeof(double), writeSize, ExtOutFile);
	fwrite(matlabpadding, sizeof(double), 1024, ExtOutFile);
	fclose(ExtOutFile);

	//Save the spectrum
	FILE* totalSpectrumFile;
	strcpy_s(outputpath, MAX_LOADSTRING, outputbase);
	strcat_s(outputpath, MAX_LOADSTRING, "_spectrum.dat");
	if(fopen_s(&totalSpectrumFile, outputpath, "wb")) return 1;
	fwrite((*sCPU).totalSpectrum, sizeof(double), 3 * (*sCPU).Nfreq * (*sCPU).Nsims * (*sCPU).Nsims2, totalSpectrumFile);
	fwrite(matlabpadding, sizeof(double), 1024, totalSpectrumFile);
	fclose(totalSpectrumFile);

	FILE* matlabfile;
	strcpy_s(outputpath, MAX_LOADSTRING, outputbase);
	strcat_s(outputpath, MAX_LOADSTRING, ".m");
	if(fopen_s(&matlabfile, outputpath, "w")) return 1;

	if (saveInputs) {
		fprintf(matlabfile, "fid = fopen('%s_ExtIn.dat','rb'); \n", outputbaseVar);
		fprintf(matlabfile, "%s_ExtIn = fread(fid, %zu, 'double'); \n", outputbaseVar, 2 * (*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2);
		fprintf(matlabfile, "%s_ExtIn = reshape(%s_ExtIn,[%zu %zu %zu]); \n", outputbaseVar, outputbaseVar, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims * (*sCPU).Nsims2);
		fprintf(matlabfile, "fclose(fid); \n");
	}

	fprintf(matlabfile, "fid = fopen('%s_Ext.dat','rb'); \n", outputbaseVar);
	fprintf(matlabfile, "%s_Ext = fread(fid, %zu, 'double'); \n", outputbaseVar, 2 * (*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2);
	fprintf(matlabfile, "%s_Ext = reshape(%s_Ext,[%zu %zu %zu]); \n", outputbaseVar, outputbaseVar, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims * (*sCPU).Nsims2);
	fprintf(matlabfile, "fclose(fid); \n");
	fprintf(matlabfile, "fid = fopen('%s_spectrum.dat','rb'); \n", outputbaseVar);
	fprintf(matlabfile, "%s_spectrum = fread(fid, %zu, 'double'); \n", outputbaseVar, 3 * (*sCPU).Nfreq * (*sCPU).Nsims * (*sCPU).Nsims2);
	fprintf(matlabfile, "%s_spectrum = reshape(%s_spectrum,[%zu %d %zu]); \n", outputbaseVar, outputbaseVar, (*sCPU).Nfreq, 3, (*sCPU).Nsims * (*sCPU).Nsims2);
	fprintf(matlabfile, "fclose(fid); \n");
	fprintf(matlabfile, "dt = %e;\ndz = %e;\ndx = %e;\ndf = %e;\n", (*sCPU).tStep, (*sCPU).propagationStep, (*sCPU).rStep, (*sCPU).fStep);
	fclose(matlabfile);

	//write a python script for loading the output fields in a proper shape
	char scriptfilename[MAX_LOADSTRING];
	strcpy_s(scriptfilename, MAX_LOADSTRING, outputbase);
	strcat_s(scriptfilename, MAX_LOADSTRING, ".py");
	FILE* scriptfile;
	if(fopen_s(&scriptfile, scriptfilename, "w")) return 1;
	fprintf(scriptfile, "#!/usr/bin/python\nimport numpy as np\n");
	fprintf(scriptfile, "dt = %e\ndz = %e\ndx = %e\ndf = %e\n", (*sCPU).tStep, (*sCPU).propagationStep, (*sCPU).rStep, (*sCPU).fStep);
	if (saveInputs) {
		fprintf(scriptfile, "%s_ExtIn = np.reshape(np.fromfile(\"", outputbaseVar);
		fprintf(scriptfile, "%s_ExtIn.dat", outputbaseVar);
		fprintf(scriptfile, "\",dtype=np.double)[0:%zu],(%zu,%zu,%zu),order='F')\n", 2 * (*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims * (*sCPU).Nsims2);
	}
	fprintf(scriptfile, "%s_Ext = np.reshape(np.fromfile(\"", outputbaseVar);
	fprintf(scriptfile, "%s_Ext.dat", outputbaseVar);
	fprintf(scriptfile, "\",dtype=np.double)[0:%zu],(%zu,%zu,%zu),order='F')\n", 2 * (*sCPU).Ngrid * (*sCPU).Nsims * (*sCPU).Nsims2, (*sCPU).Ntime, (*sCPU).Nspace, 2 * (*sCPU).Nsims * (*sCPU).Nsims2);
	fprintf(scriptfile, "%s_spectrum = np.reshape(np.fromfile(\"", outputbaseVar);
	fprintf(scriptfile, "%s_spectrum.dat", outputbaseVar);
	fprintf(scriptfile, "\",dtype=np.double)[0:%zu],(%zu,%i,%zu),order='F')\n", 3 * (*sCPU).Nfreq * (*sCPU).Nsims * (*sCPU).Nsims2, (*sCPU).Nfreq, 3, (*sCPU).Nsims * (*sCPU).Nsims2);
	fclose(scriptfile);
	return 0;
}


int configureBatchMode(simulationParameterSet* sCPU) {
	size_t i, j, currentRow;
	if ((*sCPU).batchIndex == 0 || (*sCPU).Nsims == 1 || (*sCPU).batchIndex > 35 || (*sCPU).batchIndex2 > 35) {
		return 0;
	}

	//pointers to values that can be scanned in batch mode
	double* targets[36] = { 0,
		&(*sCPU).pulse1.energy, &(*sCPU).pulse2.energy, &(*sCPU).pulse1.frequency, &(*sCPU).pulse2.frequency,
		&(*sCPU).pulse1.bandwidth, &(*sCPU).pulse2.bandwidth, &(*sCPU).pulse1.cep, &(*sCPU).pulse2.cep,
		&(*sCPU).pulse1.delay, &(*sCPU).pulse2.delay, &(*sCPU).pulse1.gdd, &(*sCPU).pulse2.gdd,
		&(*sCPU).pulse1.tod, &(*sCPU).pulse2.tod, &(*sCPU).pulse1.phaseMaterialThickness, &(*sCPU).pulse2.phaseMaterialThickness,
		&(*sCPU).pulse1.beamwaist, &(*sCPU).pulse2.beamwaist,
		&(*sCPU).pulse1.x0, &(*sCPU).pulse2.x0, &(*sCPU).pulse1.z0, &(*sCPU).pulse2.z0,
		&(*sCPU).pulse1.beamAngle, &(*sCPU).pulse2.beamAngle, &(*sCPU).pulse1.polarizationAngle, &(*sCPU).pulse2.polarizationAngle,
		&(*sCPU).pulse1.circularity, &(*sCPU).pulse2.circularity, &(*sCPU).crystalTheta, &(*sCPU).crystalPhi,
		&(*sCPU).nonlinearAbsorptionStrength, &(*sCPU).drudeGamma, &(*sCPU).effectiveMass, &(*sCPU).crystalThickness,
		&(*sCPU).propagationStep };

	//multipliers to the Batch end value from the interface
	// (e.g. frequency in THz requires 1e12 multiplier)
	double multipliers[36] = { 0,
		1, 1, 1e12, 1e12,
		1e12, 1e12, PI, PI,
		1e-15, 1e-15, 1e-30, 1e-30,
		1e-45, 1e-45, 1e-6, 1e-6,
		1e-6, 1e-6,
		1e-6, 1e-6, 1e-6, 1e-6,
		DEG2RAD, DEG2RAD, DEG2RAD, DEG2RAD,
		1, 1, DEG2RAD, DEG2RAD,
		1, 1e12, 1, 1e-6,
		1e-9 };

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


int readSequenceString(simulationParameterSet* sCPU) {
	//read the sequence string (if there is one), convert it into an array if it exists
	char sequenceString[MAX_LOADSTRING];
	strcpy_s(sequenceString, MAX_LOADSTRING, (*sCPU).sequenceString);
	char* nextToken = NULL;
	char* tokToken = strtok_s(sequenceString, ";", &nextToken);
	int sequenceCount = sscanf_s(sequenceString, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&(*sCPU).sequenceArray[0], &(*sCPU).sequenceArray[1], &(*sCPU).sequenceArray[2],
		&(*sCPU).sequenceArray[3], &(*sCPU).sequenceArray[4], &(*sCPU).sequenceArray[5],
		&(*sCPU).sequenceArray[6], &(*sCPU).sequenceArray[7], &(*sCPU).sequenceArray[8],
		&(*sCPU).sequenceArray[9], &(*sCPU).sequenceArray[10]);

	tokToken = strtok_s(NULL, ";", &nextToken);
	int lastread = sequenceCount;
	while (tokToken != NULL && lastread == 11) {
		lastread = sscanf_s(tokToken, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&(*sCPU).sequenceArray[sequenceCount], &(*sCPU).sequenceArray[sequenceCount + 1],
			&(*sCPU).sequenceArray[sequenceCount + 2], &(*sCPU).sequenceArray[sequenceCount + 3],
			&(*sCPU).sequenceArray[sequenceCount + 4], &(*sCPU).sequenceArray[sequenceCount + 5],
			&(*sCPU).sequenceArray[sequenceCount + 6], &(*sCPU).sequenceArray[sequenceCount + 7],
			&(*sCPU).sequenceArray[sequenceCount + 8], &(*sCPU).sequenceArray[sequenceCount + 9],
			&(*sCPU).sequenceArray[sequenceCount + 10]);
		if (lastread > 0) {
			sequenceCount += lastread;
		}
		tokToken = strtok_s(NULL, ";", &nextToken);
	}
	(*sCPU).Nsequence = sequenceCount / 11;
	(*sCPU).isInSequence = ((*sCPU).Nsequence > 0);

	if (!(*sCPU).isInSequence) {
		char nopeString[] = "None.";
		strcpy_s((*sCPU).sequenceString, MAX_LOADSTRING, nopeString);
	}
	return 0;
}


int readCrystalDatabase(crystalEntry* db) {
	int i = 0;
	double* fd;
	FILE* fp;
	if (fopen_s(&fp, "CrystalDatabase.txt", "r")) {
		return -2;
	}
	wchar_t lineBuffer[MAX_LOADSTRING] = { 0 };

	//read the entries line
	int readErrors = 0;

	while (readErrors == 0 && !feof(fp) && i < MAX_LOADSTRING) {
		readErrors += 0 != fwscanf_s(fp, L"Name:\n");
		fgetws(db[i].crystalNameW, 256, fp);
		readErrors += 1 != fwscanf_s(fp, L"Type:\n%d\n", &db[i].axisType);
		readErrors += 1 != fwscanf_s(fp, L"Sellmeier equation:\n%d\n", &db[i].sellmeierType);
		fd = &db[i].sellmeierCoefficients[0];
		readErrors += 22 != fwscanf_s(fp, L"1st axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			&fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
		fd = &db[i].sellmeierCoefficients[22];
		if (db[i].sellmeierType == 0 && db[i].sellmeierCoefficients[16] == 0.0 && db[i].sellmeierCoefficients[19] == 0.0) {
			db[i].sellmeierType = 100; //use a different equation, 2, for eqn 0 with no imaginary components.
		}
		readErrors += 22 != fwscanf_s(fp, L"2nd axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			&fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
		fd = &db[i].sellmeierCoefficients[44];
		readErrors += 22 != fwscanf_s(fp, L"3rd axis coefficients:\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			&fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6], &fd[7], &fd[8], &fd[9], &fd[10], &fd[11], &fd[12], &fd[13], &fd[14], &fd[15], &fd[16], &fd[17], &fd[18], &fd[19], &fd[20], &fd[21]);
		readErrors += 0 != fwscanf_s(fp, L"Sellmeier reference:\n");
		fgetws(db[i].sellmeierReference, 512, fp);
		readErrors += 1 != fwscanf_s(fp, L"chi2 type:\n%d\n", &db[i].nonlinearSwitches[0]);
		fd = &db[i].d[0];
		readErrors += 18 != fwscanf_s(fp, L"d:\n%lf %lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf %lf\n",
			&fd[0], &fd[3], &fd[6], &fd[9], &fd[12], &fd[15],
			&fd[1], &fd[4], &fd[7], &fd[10], &fd[13], &fd[16],
			&fd[2], &fd[5], &fd[8], &fd[11], &fd[14], &fd[17]);
		readErrors += 0 != fwscanf_s(fp, L"d reference:\n");
		fgetws(db[i].dReference, 512, fp);
		readErrors += 1 != fwscanf_s(fp, L"chi3 type:\n%d\nchi3:\n", &db[i].nonlinearSwitches[1]);
		//handle chi3 in a flexible way to avoid making the user write 81 zeroes when not needed
		fd = &db[i].chi3[0];
		memset(fd, 0, 81 * sizeof(double));
		if (db[i].nonlinearSwitches[1] != 1 && db[i].nonlinearSwitches[1] != 2) {
			fgetws(lineBuffer, MAX_LOADSTRING, fp);
			fgetws(lineBuffer, MAX_LOADSTRING, fp);
			fgetws(lineBuffer, MAX_LOADSTRING, fp);
		}
		else if (db[i].nonlinearSwitches[1] == 1) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 27; k++) {
					readErrors += 1 != fwscanf_s(fp, L"%lf", &fd[j + 3 * k]);
				}
				readErrors += fwscanf_s(fp, L"\n");
			}
		}
		else if (db[i].nonlinearSwitches[1] == 2) {
			readErrors += 1 != fwscanf_s(fp, L"%lf", &fd[0]);
			fgetws(lineBuffer, MAX_LOADSTRING, fp);
			fgetws(lineBuffer, MAX_LOADSTRING, fp);
			fgetws(lineBuffer, MAX_LOADSTRING, fp);
		}
		readErrors += 0 != fwscanf_s(fp, L"chi3 reference:\n");
		fgetws(db[i].chi3Reference, 512, fp);
		readErrors += 0 != fwscanf_s(fp, L"Spectral file:\n");
		fgetws(db[i].spectralFile, 512, fp);
		fd = db[i].nonlinearReferenceFrequencies;
		readErrors += 0 != fwscanf_s(fp, L"Nonlinear reference frequencies:\n");
		readErrors += 7 != fwscanf_s(fp, L"%lf %lf %lf %lf %lf %lf %lf\n",
			&fd[0], &fd[1], &fd[2], &fd[3], &fd[4], &fd[5], &fd[6]);
		readErrors += 0 != fwscanf_s(fp, L"~~~crystal end~~~\n");
		if (readErrors == 0) i++;
	}
	db[0].numberOfEntries = i;
	fclose(fp);
	return i;
}

int allocateGrids(simulationParameterSet* sCPU) {
	(*sCPU).loadedField1 = new std::complex<double>[(*sCPU).Ntime];
	(*sCPU).loadedField2 = new std::complex<double>[(*sCPU).Ntime];
	(*sCPU).ExtOut = new double[(*sCPU).Ngrid * 2 * (*sCPU).Nsims * (*sCPU).Nsims2];
	(*sCPU).EkwOut = new std::complex<double>[(*sCPU).NgridC * 2 * (*sCPU).Nsims * (*sCPU).Nsims2];
	(*sCPU).deffTensor = new double[9 * (*sCPU).Nsims * (*sCPU).Nsims2];
	(*sCPU).totalSpectrum = new double[(*sCPU).Nsims * (*sCPU).Nsims2 * (*sCPU).Nfreq * 3];
	(*sCPU).statusFlags = new int[(*sCPU).Nsims * (*sCPU).Nsims2];
	(*sCPU).fittingReference = new double[(*sCPU).Nfreq];
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
	FILE* fp;
	int maxFileSize = 16384;
	double wavelength, R, phi, complexX, complexY, f, f0, f1;
	double fmax = 0.0;
	int i, k0, k1;
	double c = 1e9 * LIGHTC; //for conversion of wavelength in nm to frequency
	double df = 0;
	double fmin = 0;
	int currentRow = 0;
	std::complex<double>* E = (std::complex<double>*)calloc(maxFileSize, sizeof(std::complex<double>));
	if (E == NULL) {
		return -2;
	}
	//read the data
	if (fopen_s(&fp, frogFilePath, "r")) {
		free(E);
		return -1;
	}
	while (fscanf_s(fp, "%lf %lf %lf %lf %lf", &wavelength, &R, &phi, &complexX, &complexY) == 5 && currentRow < maxFileSize) {
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
	fclose(fp);

	//return an error if nothing was loaded
	if (currentRow == 0) {
		free(E);
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

	free(E);
	return currentRow;
}
