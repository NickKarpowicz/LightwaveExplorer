#include "LightwaveExplorerUtilities.h"
#include "LightwaveExplorerInterfaceClasses.hpp"
#include "LightwaveExplorerInterfaceClasses.cpp"
#include <iostream>
#include <sstream>


//remove whitespace and line breaks, except inside of protected blocks demarcated with characters in the
//startChars and endChars strings
//note that ANY character in the start or stop string will do.
int removeCharacterFromStringSkippingChars(
	std::string& s,
	const char removedChar,
	const std::string& startChars,
	const std::string& endChars) {
	bool removing = true;
	for (std::size_t i = 0; i < s.length(); ++i) {
		if (s[i] == removedChar && removing) {
			s.erase(i,1);
			--i;
		}
		if (removing && startChars.find(s[i]) != std::string::npos && endChars.find(s[i]) != std::string::npos){
			removing = !removing;
		}
		else if (startChars.find(s[i]) != std::string::npos) removing = false;
		else if (endChars.find(s[i]) != std::string::npos) removing = true;
	}
	return 0;
}

void stripWhiteSpace(std::string& s) {
	removeCharacterFromStringSkippingChars(s, ' ', "<\"", ">\"");
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
		std::size_t startArgument = 0;
		std::vector<std::string> argTable;
		argTable.reserve(numberParams);
		char expectedDelimiter = ',';
		char wrongDelimiter = ')';
		int openParens = 0;

		for(std::size_t i = 0; i<arguments.size(); ++i){
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

		for(int i = 0; i<static_cast<int>(argTable.size()); ++i){
			if(argTable[i].at(0)=='d'){
				defaultMask[i] = true;
			}
			else{
				parameters[i] = parameterStringToDouble(argTable[i],iBlock,vBlock);
			}
		}
		return 0;
}

std::size_t findParenthesesClosure(std::string& a){
	int nParen = 0;
	bool foundFirstParen = false;
	for(std::size_t i = 0; i<a.size(); ++i){
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
	std::vector<std::size_t> operatorTable;
	std::vector<double> numberTable;
	int lastOperator = -1;
	bool lastNumberWasNotParenthesized = true;
	for(std::size_t i = 0; i<ss.size(); ++i){
		if(
		  (ss[i] == '+'
		|| ss[i] == '-'
		|| ss[i] == '*'
		|| ss[i] == '/'
		|| ss[i] == '^')
		&& i != 0){
			operatorTable.push_back(i);
			if(lastNumberWasNotParenthesized)numberTable.push_back(std::stod(ss.substr(lastOperator+1,i-lastOperator-1)));
			lastOperator = static_cast<int>(i);
			lastNumberWasNotParenthesized = true;
		}
		else if(ss[i] == '('){
			std::string parenString = ss.substr(i, std::string::npos);
			parenString = parenString.substr(1,
				findParenthesesClosure(
					parenString)-1);
			numberTable.push_back(parameterStringToDouble(parenString, iBlock, vBlock));
			lastOperator = static_cast<int>(i + parenString.size());
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
			lastOperator = static_cast<int>(i + 2);
			i += 2;
			lastNumberWasNotParenthesized = false;
		}
		else if(ss[i] == 'i'){
			int ind = std::stoi(ss.substr(i+1,2));
			numberTable.push_back(iBlock[ind]);
			lastOperator = static_cast<int>(i + 2);
			i += 2;
			lastNumberWasNotParenthesized = false;
		}
	}
	if(lastOperator == -1) {
		return std::stod(ss);
	}
	if(lastOperator+1 < static_cast<int>(ss.size())){
			numberTable.push_back(std::stod(ss.substr(lastOperator+1,std::string::npos)));
		}

	auto applyOperators = [&](char firstOp, char secondOp){
		for(std::size_t i = 0; i < operatorTable.size(); ++i){
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
	std::size_t positionOfName = pathString.find_last_of("/\\");
	if (positionOfName == std::string::npos) return pathString;
	return pathString.substr(positionOfName + 1);
}

//calculates the squared modulus of a complex number
inline double cModulusSquared(const std::complex<double>& x) {
	return x.real()*x.real() + x.imag()*x.imag();
}

int loadFrogSpeck(
	const loadedInputData& frogSpeckData,
	std::complex<double>* Egrid,
	const int64_t Ntime,
	const double fStep,
	const double gateLevel) {
	std::string line;
	std::stringstream fs(frogSpeckData.fileContents);
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
	const loadedInputData& waveformFile,
	std::complex<double>* outputGrid,
	const int64_t Ntime,
	const double fStep) {
	std::vector<double> Ein;
	Ein.reserve(8192);
	std::stringstream fs(waveformFile.fileContents);
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
	// fftw_plan fftwPlanD2Z = fftw_plan_dft_r2c_1d(
	// 	lineCount, Ein.data(),
	// 	reinterpret_cast<fftw_complex*>(fftOfEin.data()),
	// 	FFTW_ESTIMATE);
	// fftw_execute_dft_r2c(
	// 	fftwPlanD2Z,
	// 	Ein.data(),
	// 	reinterpret_cast<fftw_complex*>(fftOfEin.data()));
	// fftw_destroy_plan(fftwPlanD2Z);
	//
	pocketfft::r2c(
                        {static_cast<size_t>(lineCount)},
                        {sizeof(float)},
                        {sizeof(std::complex<double>)},
                        {0}, pocketfft::FORWARD,
                        Ein.data(), fftOfEin.data(),
                        1.0);

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
	const loadedInputData& file,
	std::vector<double>& outputGrid,
	const int64_t Ngrid) {
	outputGrid.resize(Ngrid);
	if (file.hasData) {
		std::stringstream Efile(file.fileContents);
		Efile.read(
			reinterpret_cast<char*>(outputGrid.data()),
			2 * Ngrid * sizeof(double));
		return 0;
	}
	else return 1;
}

int loadSavedGridFileMultiple(
	const loadedInputData& file,
	std::vector<double>& outputGrid,
	const int64_t Ngrid,
	const int64_t Nsims) {
	outputGrid.resize(Ngrid * Nsims);
	if (file.hasData) {
		std::stringstream Efile(file.fileContents);
		Efile.read(reinterpret_cast<char*>(outputGrid.data()), 2 * Ngrid * Nsims * sizeof(double));
		return 0;
	}
	else return 1;
}
