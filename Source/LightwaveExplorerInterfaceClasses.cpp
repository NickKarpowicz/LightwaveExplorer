#include "LightwaveExplorerInterfaceClasses.hpp"
#include "LightwaveExplorerUtilities.h"
#ifdef USEFFTW
#include <fftw3.h>
#else
#include <fftw3_mkl.h>
#endif
#include <algorithm>

int simulationParameterSet::loadSavedFields(const std::string &outputBase,
                                            bool isZipFile) {
  if (isZipFile) {
    std::string zipPath = outputBase + ".zip";
    std::string trueBase = zipGetBasename(zipPath);
    zipIntoMemory(zipPath, trueBase + "_Ext.dat", ExtOut,
                  2 * (Ngrid * Nsims * Nsims2) * sizeof(double));
    zipIntoMemory(zipPath, trueBase + "_spectrum.dat", totalSpectrum,
                  Nsims * Nsims2 * 3 * Nfreq * sizeof(double));
    pulse1LoadedData = loadedInputData(zipPath, trueBase + "_pulse1.dat");
    pulse2LoadedData = loadedInputData(zipPath, trueBase + "_pulse2.dat");
    fittingLoadedData =
        loadedInputData(zipPath, trueBase + "_fittingTarget.dat");
  } else {
    std::string Epath = outputBase;
    Epath.append("_Ext.dat");
    std::ifstream Efile(Epath, std::ios::binary);
    if (Efile.is_open()) {
      Efile.read(reinterpret_cast<char *>(ExtOut),
                 2 * (Ngrid * Nsims * Nsims2) * sizeof(double));
    } else
      return 1;

    std::string Spath = outputBase;
    Spath.append("_spectrum.dat");
    std::ifstream Sfile(Spath, std::ios::binary);
    if (Sfile.is_open())
      Sfile.read(reinterpret_cast<char *>(totalSpectrum),
                 Nsims * Nsims2 * 3 * Nfreq * sizeof(double));
  }

  fftw_plan fftwPlanD2Z;
  if (is3D) {
    const int fftwSizes[] = {(int)Nspace2, (int)Nspace, (int)Ntime};
    fftwPlanD2Z = fftw_plan_many_dft_r2c(3, fftwSizes, 2, ExtOut, 0, 1,
                                         (int)Ngrid, (fftw_complex *)EkwOut, 0,
                                         1, (int)NgridC, FFTW_ESTIMATE);
  } else {
    const int fftwSizes[] = {(int)Nspace, (int)Ntime};
    fftwPlanD2Z = fftw_plan_many_dft_r2c(2, fftwSizes, 2, ExtOut, 0, 1,
                                         (int)Ngrid, (fftw_complex *)EkwOut, 0,
                                         1, (int)NgridC, FFTW_ESTIMATE);
  }

  for (int64_t i = 0; i < (Nsims * Nsims2); i++) {
    fftw_execute_dft_r2c(fftwPlanD2Z, &ExtOut[2 * i * Ngrid],
                         (fftw_complex *)&EkwOut[2 * i * NgridC]);
  }
  fftw_destroy_plan(fftwPlanD2Z);

  return 0;
}

int simulationParameterSet::loadReferenceSpectrum() {
  if (!fittingLoadedData.hasData) {
    return 1;
  }
  std::stringstream fs(fittingLoadedData.fileContents);

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
    // fs >> loadedWavelengths[currentRow] >> loadedIntensities[currentRow];
    fs >> lastRead;
    loadedWavelengths.push_back(lastRead);
    fs >> lastRead;
    loadedIntensities.push_back(lastRead);
    if (currentRow == 0) {
      maxWavelength = loadedWavelengths[currentRow];
      minWavelength = loadedWavelengths[currentRow];
    } else {
      maxWavelength = maxN(maxWavelength, loadedWavelengths[currentRow]);
      minWavelength = minN(minWavelength, loadedWavelengths[currentRow]);
    }
    // rescale to frequency spacing
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
    if ((currentFrequency > minFrequency) &&
        (currentFrequency < maxFrequency)) {
      // find the first frequency greater than the current value
      j = sizeData - 1;
      while ((loadedFrequencies[j] <= currentFrequency) && (j > 2)) {
        j--;
      }
      df = loadedFrequencies[j] - loadedFrequencies[j - 1];
      fittingReference[i] =
          (loadedIntensities[j - 1] *
               (loadedFrequencies[j] - currentFrequency) +
           loadedIntensities[j] *
               (currentFrequency - loadedFrequencies[j - 1])) /
          df;
      // linear interpolation
    }
  }
  return 0;
}

double simulationParameterSet::saveSlurmScript(
    const std::string &gpuType, int gpuCount, bool useJobArray,
    int64_t totalSteps, std::vector<simulationParameterSet> &params,
    const class crystalDatabase &db) {
  std::string scriptPath = getBasename(outputBasePath) + ".slurmScript";
  std::string Zpath = outputBasePath + ".zip";
  std::stringstream fs;
  int64_t voxel_count = (totalSteps / (Nsims * Nsims2)) * Nspace * Ntime;
  // Estimate the time required on the cluster, proportional to number of grid
  // points x steps
  double timeEstimate =
      static_cast<double>(voxel_count) *
      ceil(static_cast<double>(Nsims * Nsims2) / gpuCount);
  // 3D
  if (symmetryType == 2) {
    timeEstimate *= Nspace2;
    timeEstimate *= 1.0e-11;
  }
  // Radial symmetry
  else if (symmetryType == 1) {
    timeEstimate *= 2.0e-11;
  } else if (symmetryType == 4) {
    timeEstimate =
        1e-11 *
        static_cast<double>(voxel_count * Nspace2) * (crystalThickness / propagationStep) *
        ceil(static_cast<double>(Nsims * Nsims2) / gpuCount);
  }
  // 2D
  else {
    timeEstimate *= 1.0e-11;
  }
  // Plasma doubles estimate
  if (nonlinearAbsorptionStrength != 0.0)
    timeEstimate *= 2.1;
  if (gpuType != "a100")
    timeEstimate *= 8;  // if it's not an A100, assume its slow.
  timeEstimate *= 1.25; // safety margin
  timeEstimate += 60.0; // fixed offset for loading .etc
  int timeEstimateHours =
      static_cast<int>(timeEstimate / 3600.0); // convert to hours
  int timeEstimateMinutes =
      static_cast<int>(timeEstimate / 60.0) - 60 * timeEstimateHours;
  if (timeEstimateHours >= 24) {
    timeEstimateHours = 24;
    timeEstimateMinutes = 0;
  }
  if (timeEstimateHours == 0 && timeEstimateMinutes < 15) {
    timeEstimateMinutes = 15;
  }

  if (fittingString[0] != 0 && fittingString[0] != 'N') {
    readFittingString();
    if (fittingMaxIterations > 0)
      timeEstimate *= fittingMaxIterations;
  }
  std::string baseName = getBasename(outputBasePath);
  int memoryMB = static_cast<int>(
      (18 * sizeof(double) * Ngrid * maxN(Nsims, 1u)) / 1048576);
  // TODO: needs a better estimate for FDTD mode

  if (useJobArray)
    memoryMB /= static_cast<int>(Nsims);
  memoryMB += 8192; // base level
  if (useJobArray)
    gpuCount = 1;

  fs << "#!/bin/bash -l" << '\x0A';
  fs << "#SBATCH -o ./tjob.out.%j" << '\x0A';
  fs << "#SBATCH -e ./tjob.err.%j" << '\x0A';
  fs << "#SBATCH -D ./" << '\x0A';
  fs << "#SBATCH -J lightwave" << '\x0A';
  fs << "#SBATCH --constraint=\"gpu\"" << '\x0A';
  fs << "#SBATCH --gres=gpu:" << gpuType << ":" << gpuCount << '\x0A';
  fs << "#SBATCH --cpus-per-task=" << 1 + gpuCount << '\x0A';

  fs << "#SBATCH --mem=" << memoryMB << "M\x0A";
  fs << "#SBATCH --nodes=1" << '\x0A';
  fs << "#SBATCH --ntasks-per-node=1" << '\x0A';
  fs << "#SBATCH --time=" << timeEstimateHours << ':' << timeEstimateMinutes
     << ":00" << '\x0A';
  if (useJobArray) {
    fs << "#SBATCH --array=0-" << Nsims * Nsims2 - 1 << '\x0A';
  }

  fs << "module purge" << '\x0A';
  fs << "module load cuda/12.2" << '\x0A';
  fs << "module load mkl/2024.0" << '\x0A';
  fs << "export LD_LIBRARY_PATH=$MKL_HOME/lib/intel64:$LD_LIBRARY_PATH"
     << '\x0A';
  fs << "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}" << '\x0A';

  if (useJobArray) {
    fs << "file_id=$(printf \"%05d\" $SLURM_ARRAY_TASK_ID)\x0a";
    fs << "base_name=\"" << baseName << "$file_id\"" << '\x0A';
  } else {
    fs << "base_name=\"" << baseName << "\"" << '\x0A';
  }

  fs << "srun ../lwe $base_name.input > $base_name.out\x0A";
  // optionally upload to a webdav location, if a token is provided
  fs << "if [ -f ../webdav_token.txt ]; then" << '\x0A';
  fs << "    webdav_token=$(<../webdav_token.txt)" << '\x0A';
  fs << "    webdav_url=$(<../webdav_url.txt)" << '\x0A';
  fs << "    curl --user $webdav_token:nopass $webdav_url --upload-file "
        "\"$base_name\".out"
     << '\x0A';
  fs << "    curl --user $webdav_token:nopass $webdav_url --upload-file "
        "\"$base_name\".zip"
     << '\x0A';
  fs << "    rm \"$base_name\".out" << '\x0A';
  fs << "    rm \"$base_name\".zip" << '\x0A';
  fs << "fi" << '\x0A';
  std::string script = fs.str();
  mz_zip_archive zip = {};
  mz_zip_writer_init_file(&zip, Zpath.c_str(), 0);
  mz_zip_writer_add_mem(&zip, getBasename(scriptPath).c_str(), script.c_str(),
                        script.size(), MZ_DEFAULT_COMPRESSION);

  // scan the sequence for quotes, include the data if something is found
  std::size_t startPosition = sequenceString.find_first_of('\"');
  std::string workingString = sequenceString;
  std::string modifiedString = sequenceString;
  while (startPosition != std::string::npos) {
    workingString = workingString.substr(startPosition + 1);
    std::size_t endPosition = workingString.find_first_of('\"');
    loadedInputData newData(workingString.substr(0, endPosition));
    std::size_t pos = modifiedString.find(newData.filePath);
    modifiedString.replace(pos, newData.filePath.length(),
                           getBasename(newData.filePath));
    mz_zip_writer_add_mem(&zip, getBasename(newData.filePath).c_str(),
                          newData.fileContents.c_str(),
                          newData.fileContents.size(), MZ_DEFAULT_COMPRESSION);
    workingString = workingString.substr(endPosition + 1);
    startPosition = workingString.find_first_of('\"');
  }

  if (useJobArray) {
    int simIndex = 0;
    std::string settings = settingsString();
    std::string base = getBasename(outputBasePath);
    std::string mainSettings = base + ".txt";
    mz_zip_writer_add_mem(&zip, getBasename(mainSettings).c_str(),
                          settings.c_str(), settings.size(),
                          MZ_DEFAULT_COMPRESSION);
    int64_t loopNsims2 = Nsims2;
    int64_t loopNsims = Nsims;
    std::vector<double> i37(Nsims);
    std::vector<double> i37_2(Nsims2);
    if (batchIndex == 37) {
      for (int64_t i = 0; i < Nsims; ++i) {
        i37[i] = i * batchDestination / (Nsims - 1);
      }
    }
    if (batchIndex2 == 37) {
      for (int64_t i = 0; i < Nsims2; ++i) {
        i37_2[i] = i * batchDestination2 / (Nsims2 - 1);
      }
    }
    for (int64_t i = 0; i < loopNsims2; ++i) {
      for (int64_t j = 0; j < loopNsims; ++j) {
        simulationParameterSet arraySim = params[i * loopNsims + j];
        arraySim.Nsims = 1;
        arraySim.Nsims2 = 1;
        std::ostringstream oss;
        oss << std::setw(5) << std::setfill('0') << simIndex++;
        arraySim.outputBasePath.append(oss.str());
        arraySim.runType = runTypes::cluster;
        arraySim.batchIndex = 0;
        arraySim.batchIndex2 = 0;
        arraySim.runType = runTypes::cluster;
        arraySim.sequenceString = modifiedString;
        if (batchIndex == 37) {
          std::ostringstream oss;
          oss << std::setprecision(15) << i37[j];
          std::string replacement = oss.str();
          std::size_t pos = arraySim.sequenceString.find("i37");
          while (pos != std::string::npos) {
            arraySim.sequenceString.replace(pos, 3, replacement);
            pos = arraySim.sequenceString.find("i37", pos + 3);
          }
        }
        if (batchIndex2 == 37) {
          std::ostringstream oss;
          oss << std::setprecision(15) << i37_2[i];
          std::string replacement = oss.str();
          std::size_t pos = arraySim.sequenceString.find("i37");
          while (pos != std::string::npos) {
            arraySim.sequenceString.replace(pos, 3, replacement);
            pos = arraySim.sequenceString.find("i37", pos + 3);
          }
        }
        settings = arraySim.settingsString();
        std::string currentPath =
            getBasename(arraySim.outputBasePath) + ".input";
        mz_zip_writer_add_mem(&zip, getBasename(currentPath).c_str(),
                              settings.c_str(), settings.size(),
                              MZ_DEFAULT_COMPRESSION);
      }
    }
  } else {
    runType = runTypes::cluster;
    std::string settings = settingsString();
    std::string mainSettings = getBasename(outputBasePath) + ".input";
    mz_zip_writer_add_mem(&zip, getBasename(mainSettings).c_str(),
                          settings.c_str(), settings.size(),
                          MZ_DEFAULT_COMPRESSION);
  }
  loadedInputData dbData(db.path);
  mz_zip_writer_add_mem(&zip, "CrystalDatabase.txt",
                        dbData.fileContents.c_str(), dbData.fileContents.size(),
                        MZ_DEFAULT_COMPRESSION);
  std::string FittingTargetPath = outputBasePath + "_fittingTarget.dat";
  std::string Pulse1Path = outputBasePath + "_pulse1.dat";
  std::string Pulse2Path = outputBasePath + "_pulse2.dat";
  if (pulse1LoadedData.hasData) {
    mz_zip_writer_add_mem(&zip, getBasename(Pulse1Path).c_str(),
                          pulse1LoadedData.fileContents.c_str(),
                          pulse1LoadedData.fileContents.size(),
                          MZ_DEFAULT_COMPRESSION);
  }
  if (pulse2LoadedData.hasData) {
    mz_zip_writer_add_mem(&zip, getBasename(Pulse2Path).c_str(),
                          pulse2LoadedData.fileContents.c_str(),
                          pulse2LoadedData.fileContents.size(),
                          MZ_DEFAULT_COMPRESSION);
  }
  if (fittingLoadedData.hasData) {
    mz_zip_writer_add_mem(&zip, getBasename(FittingTargetPath).c_str(),
                          fittingLoadedData.fileContents.c_str(),
                          fittingLoadedData.fileContents.size(),
                          MZ_DEFAULT_COMPRESSION);
  }
  if (!optics.empty()) {
    for (int i = 0; i < static_cast<int>(optics.size()); i++) {
      std::string opticPath =
          outputBasePath + "_optic" + std::to_string(i) + ".txt";
      mz_zip_writer_add_mem(
          &zip, getBasename(opticPath).c_str(), optics[i].fileContents.c_str(),
          optics[i].fileContents.size(), MZ_DEFAULT_COMPRESSION);
    }
  }

  mz_zip_writer_finalize_archive(&zip);
  mz_zip_writer_end(&zip);

  return timeEstimate / 3600.0;
}

std::string simulationParameterSet::settingsString() {
  std::string baseName = getBasename(outputBasePath);
  std::string referenceBaseName = getBasename(fittingLoadedData.filePath);
  std::string pulse1BaseName = getBasename(pulse1LoadedData.filePath);
  std::string pulse2BaseName = getBasename(pulse2LoadedData.filePath);
  std::stringstream fs;
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
  fs << "Phase material thickness 1 (mcr.): " << pulse1.phaseMaterialThickness
     << '\x0A';
  fs << "Phase material thickness 2 (mcr.): " << pulse2.phaseMaterialThickness
     << '\x0A';
  fs << "Beam mode placeholder: " << 0 << '\x0A';
  fs << "Beamwaist 1 (m): " << pulse1.beam_spec.waist[0][0] << '\x0A';
  fs << "Beamwaist 2 (m): " << pulse2.beam_spec.waist[0][0] << '\x0A';
  fs << "x offset 1 (m): " << pulse1.beam_spec.x_offset[0][0] << '\x0A';
  fs << "x offset 2 (m): " << pulse2.beam_spec.x_offset[0][0] << '\x0A';
  fs << "y offset 1 (m): " << pulse1.beam_spec.y_offset[0][0] << '\x0A';
  fs << "y offset 2 (m): " << pulse2.beam_spec.y_offset[0][0] << '\x0A';
  fs << "z offset 1 (m): " << pulse1.beam_spec.z_offset[0][0] << '\x0A';
  fs << "z offset 2 (m): " << pulse2.beam_spec.z_offset[0][0] << '\x0A';
  fs << "NC angle 1 (rad): " << pulse1.beam_spec.angle_x[0][0] << '\x0A';
  fs << "NC angle 2 (rad): " << pulse2.beam_spec.angle_x[0][0] << '\x0A';
  fs << "NC angle phi 1 (rad): " << pulse1.beam_spec.angle_y[0][0] << '\x0A';
  fs << "NC angle phi 2 (rad): " << pulse2.beam_spec.angle_y[0][0] << '\x0A';
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
  fs << "Nonlinear absorption parameter: " << nonlinearAbsorptionStrength
     << '\x0A';
  fs << "Initial carrier density (m^-3): " << startingCarrierDensity << '\x0A';
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
  if (runType ==
      runTypes::cluster) { // don't include full path if making a cluster script
    fs << "Output base path: " << baseName << '\x0A';
    fs << "Field 1 from file type: " << pulse1FileType << '\x0A';
    fs << "Field 2 from file type: " << pulse2FileType << '\x0A';
    fs << "Field 1 file path: " << pulse1BaseName << '\x0A';
    fs << "Field 2 file path: " << pulse2BaseName << '\x0A';
    fs << "Fitting reference file path: " << referenceBaseName << '\x0A';
  } else {
    fs << "Output base path: " << outputBasePath << '\x0A';
    fs << "Field 1 from file type: " << pulse1FileType << '\x0A';
    fs << "Field 2 from file type: " << pulse2FileType << '\x0A';
    fs << "Field 1 file path: " << pulse1LoadedData.filePath << '\x0A';
    fs << "Field 2 file path: " << pulse2LoadedData.filePath << '\x0A';
    fs << "Fitting reference file path: " << fittingLoadedData.filePath
       << '\x0A';
  }

  fs << "Material name: " << crystalDatabase[materialIndex].crystalName
     << '\x0A';
  fs << "Sellmeier reference: "
     << crystalDatabase[materialIndex].sellmeierReference << '\x0A';
  fs << "Chi2 reference: " << crystalDatabase[materialIndex].dReference
     << '\x0A';
  fs << "Chi3 reference: " << crystalDatabase[materialIndex].chi3Reference
     << '\x0A';
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 22; ++j) {
      fs << crystalDatabase[materialIndex].sellmeierCoefficients[i * 22 + j];
      if (j < 21)
        fs << ',';
    }
    fs << '\x0A';
  }
  fs << "Code version: 2025.6";
  fs << '\x0A';
  return fs.str();
}

int simulationParameterSet::saveSettingsFile() {
  std::string outputPath(outputBasePath);
  if (runType == runTypes::cluster) {
    outputPath.append(".input");
  } else {
    outputPath.append(".txt");
  }
  std::ofstream fs(outputPath, std::ios::binary);
  if (fs.fail())
    return -1;
  fs << settingsString();
  return 0;
}

void simulationParameterSet::setByNumber(const int64_t index,
                                         const double value) {
  switch (index) {
  case 0:
    return;
  case 1:
    pulse1.energy = value;
    return;
  case 2:
    pulse2.energy = value;
    return;
  case 3:
    pulse1.frequency = value;
    return;
  case 4:
    pulse2.frequency = value;
    return;
  case 5:
    pulse1.bandwidth = value;
    return;
  case 6:
    pulse2.bandwidth = value;
    return;
  case 7:
    pulse1.cep = value;
    return;
  case 8:
    pulse2.cep = value;
    return;
  case 9:
    pulse1.delay = value;
    return;
  case 10:
    pulse2.delay = value;
    return;
  case 11:
    pulse1.gdd = value;
    return;
  case 12:
    pulse2.gdd = value;
    return;
  case 13:
    pulse1.tod = value;
    return;
  case 14:
    pulse2.tod = value;
    return;
  case 15:
    pulse1.phaseMaterialThickness = value;
    return;
  case 16:
    pulse2.phaseMaterialThickness = value;
    return;
  case 17:
    pulse1.beam_spec.waist[0][0] = value;
    return;
  case 18:
    pulse2.beam_spec.waist[0][0] = value;
    return;
  case 19:
    pulse1.beam_spec.x_offset[0][0] = value;
    return;
  case 20:
    pulse2.beam_spec.x_offset[0][0] = value;
    return;
  case 21:
    pulse1.beam_spec.z_offset[0][0] = value;
    return;
  case 22:
    pulse2.beam_spec.z_offset[0][0] = value;
    return;
  case 23:
    pulse1.beam_spec.angle_x[0][0] = value;
    return;
  case 24:
    pulse2.beam_spec.angle_x[0][0] = value;
    return;
  case 25:
    pulse1.polarizationAngle = value;
    return;
  case 26:
    pulse2.polarizationAngle = value;
    return;
  case 27:
    pulse1.circularity = std::clamp(value, -1.0, 1.0);
    return;
  case 28:
    pulse2.circularity = std::clamp(value, -1.0, 1.0);
    return;
  case 29:
    crystalTheta = value;
    return;
  case 30:
    crystalPhi = value;
    return;
  case 31:
    nonlinearAbsorptionStrength = value;
    return;
  case 32:
    drudeGamma = value;
    return;
  case 33:
    effectiveMass = value;
    return;
  case 34:
    crystalThickness = value;
    return;
  case 35:
    propagationStep = value;
    return;
  case 36:
    return;
  case 37:
    i37 = value;
    return;
  default:
    return;
  }
}

void simulationParameterSet::setByNumberWithMultiplier(const std::size_t index,
                                                       const double value) {
  if (index > multipliers.size())
    return;
  setByNumber(index, value * multipliers[index]);
}

int simulationParameterSet::readInputParametersFile(
    crystalEntry *crystalDatabasePtr, const std::string filePath) {
  std::string contents;
  if (filePath.length() >= 4 &&
      filePath.substr(filePath.length() - 4) == ".zip") {
    std::vector<char> dataVector;
    std::string textPath = zipGetBasename(filePath);
    textPath = textPath + ".txt";
    zipIntoMemory(filePath, textPath, dataVector);
    dataVector.push_back(0);
    contents = std::string(dataVector.data(), dataVector.size());
  } else {
    std::ifstream file(filePath);
    if (file.fail())
      return 1;
    contents = std::string((std::istreambuf_iterator<char>(file)),
                           std::istreambuf_iterator<char>());
  }

  std::string line;
  std::stringstream fs(contents);
  if (fs.fail())
    return 1;

  auto moveToColon = [&fs]() {
    fs.ignore(std::numeric_limits<std::streamsize>::max(), ':');
  };

  // Check if the next line in the stream contains the label expectedLabel
  // followed by a colon if it does, return the number that follows as a number
  // if there's no colon, or the label isn't what was expected,
  // return the stream to the starting position and return the default value
  auto checkLineName = [&fs](const std::string &expectedLabel,
                             double defaultValue) {
    std::string line;
    std::streampos startingPosition = fs.tellg();
    std::getline(fs, line);
    std::size_t colon = line.find(':');
    if (colon == std::string::npos) {
      fs.seekg(startingPosition);
      return defaultValue;
    }
    std::string label = line.substr(0, colon);
    std::string value = line.substr(colon + 1);
    std::size_t firstNonWhitespace = value.find_first_not_of(" \t");
    if (firstNonWhitespace != std::string::npos && firstNonWhitespace != 0) {
      value.erase(0, firstNonWhitespace);
    }
    value.erase(value.find_last_not_of(" \t") + 1);
    if (label == expectedLabel) {
      return std::stod(value);
    }
    fs.seekg(startingPosition);
    return defaultValue;
  };

  pulse1.energy = checkLineName("Pulse energy 1 (J)", 0.0);
  pulse2.energy = checkLineName("Pulse energy 2 (J)", 0.0);
  pulse1.frequency = checkLineName("Frequency 1 (Hz)", 0.0);
  pulse2.frequency = checkLineName("Frequency 2 (Hz)", 0.0);
  pulse1.bandwidth = checkLineName("Bandwidth 1 (Hz)", 0.0);
  pulse2.bandwidth = checkLineName("Bandwidth 2 (Hz)", 0.0);
  pulse1.sgOrder = checkLineName("SG order 1", 6);
  pulse2.sgOrder = checkLineName("SG order 2", 6);
  pulse1.cep = checkLineName("CEP 1 (rad)", 0.0);
  pulse2.cep = checkLineName("CEP 2 (rad)", 0.0);
  pulse1.delay = checkLineName("Delay 1 (s)", 0.0);
  pulse2.delay = checkLineName("Delay 2 (s)", 0.0);
  pulse1.gdd = checkLineName("GDD 1 (s^-2)", 0.0);
  pulse2.gdd = checkLineName("GDD 2 (s^-2)", 0.0);
  pulse1.tod = checkLineName("TOD 1 (s^-3)", 0.0);
  pulse2.tod = checkLineName("TOD 2 (s^-3)", 0.0);
  pulse1.phaseMaterial = checkLineName("Phase material 1 index", 0.0);
  pulse2.phaseMaterial = checkLineName("Phase material 2 index", 0.0);
  pulse1.phaseMaterialThickness =
      checkLineName("Phase material thickness 1 (mcr.)", 0.0);
  pulse2.phaseMaterialThickness =
      checkLineName("Phase material thickness 2 (mcr.)", 0.0);
  checkLineName("Beam mode placeholder", 0.0);
  pulse1.beam_spec.waist[0][0] = checkLineName("Beamwaist 1 (m)", 0.0);
  pulse2.beam_spec.waist[0][0] = checkLineName("Beamwaist 2 (m)", 0.0);
  pulse1.beam_spec.x_offset[0][0] = checkLineName("x offset 1 (m)", 0.0);
  pulse2.beam_spec.x_offset[0][0] = checkLineName("x offset 2 (m)", 0.0);
  pulse1.beam_spec.y_offset[0][0] = checkLineName("y offset 1 (m)", 0.0);
  pulse2.beam_spec.y_offset[0][0] = checkLineName("y offset 2 (m)", 0.0);
  pulse1.beam_spec.z_offset[0][0] = checkLineName("z offset 1 (m)", 0.0);
  pulse2.beam_spec.z_offset[0][0] = checkLineName("z offset 2 (m)", 0.0);
  pulse1.beam_spec.angle_x[0][0] = checkLineName("NC angle 1 (rad)", 0.0);
  pulse2.beam_spec.angle_x[0][0] = checkLineName("NC angle 2 (rad)", 0.0);
  pulse1.beam_spec.angle_y[0][0] = checkLineName("NC angle phi 1 (rad)", 0.0);
  pulse2.beam_spec.angle_y[0][0] = checkLineName("NC angle phi 2 (rad)", 0.0);
  pulse1.polarizationAngle = checkLineName("Polarization 1 (rad)", 0.0);
  pulse2.polarizationAngle = checkLineName("Polarization 2 (rad)", 0.0);
  pulse1.circularity = checkLineName("Circularity 1", 0.0);
  pulse1.circularity = std::clamp(pulse1.circularity, -1.0, 1.0);
  pulse2.circularity = checkLineName("Circularity 2", 0.0);
  pulse2.circularity = std::clamp(pulse2.circularity, -1.0, 1.0);
  materialIndex = checkLineName("Material index", 0.0);
  materialIndexAlternate = checkLineName("Alternate material index", 0.0);
  crystalTheta = checkLineName("Crystal theta (rad)", 0.0);
  crystalPhi = checkLineName("Crystal phi (rad)", 0.0);
  spatialWidth = checkLineName("Grid width (m)", 0.0);
  spatialHeight = checkLineName("Grid height (m)", 0.0);
  rStep = checkLineName("dx (m)", 0.0);
  timeSpan = checkLineName("Time span (s)", 0.0);
  tStep = checkLineName("dt (s)", 0.0);
  crystalThickness = checkLineName("Thickness (m)", 0.0);
  propagationStep = checkLineName("dz (m)", 0.0);
  nonlinearAbsorptionStrength =
      checkLineName("Nonlinear absorption parameter", 0.0);
  startingCarrierDensity = checkLineName("Initial carrier density (m^-3)", 0.0);
  bandGapElectronVolts = checkLineName("Band gap (eV)", 0.0);
  effectiveMass = checkLineName("Effective mass (relative)", 0.0);
  drudeGamma = checkLineName("Drude gamma (Hz)", 0.0);
  symmetryType = checkLineName("Propagation mode", 0.0);
  batchIndex = checkLineName("Batch mode", 0.0);
  batchDestination = checkLineName("Batch destination", 0.0);
  Nsims = checkLineName("Batch steps", 0.0);
  batchIndex2 = checkLineName("Batch mode 2", 0.0);
  batchDestination2 = checkLineName("Batch destination 2", 0.0);
  Nsims2 = checkLineName("Batch steps 2", 0.0);

  moveToColon();

  std::getline(fs, line);
  if(line.length()>0) line.erase(line.begin());

  sequenceString = line;

  moveToColon();
  std::getline(fs, line);

  if(line.length()>0) line.erase(line.begin());

  fittingString = line;
  moveToColon();
  fs >> fittingMode;

  moveToColon();
  std::getline(fs, line);
  if(line.length()>0) line.erase(line.begin());

  outputBasePath = line;
  moveToColon();
  fs >> pulse1FileType;
  moveToColon();
  fs >> pulse2FileType;
  moveToColon();
  std::getline(fs, line);
  if(line.length()>0) line.erase(line.begin());
  removeCharacterFromString(line, '\r');
  removeCharacterFromString(line, '\n');
  pulse1LoadedData = loadedInputData(line);
  moveToColon();
  std::getline(fs, line);
  if(line.length()>0) line.erase(line.begin());
  removeCharacterFromString(line, '\r');
  removeCharacterFromString(line, '\n');
  pulse2LoadedData = loadedInputData(line);
  moveToColon();
  std::getline(fs, line);
  if(line.length()>0) line.erase(line.begin());
  removeCharacterFromString(line, '\r');
  removeCharacterFromString(line, '\n');
  fittingLoadedData = loadedInputData(line);

  removeCharacterFromString(sequenceString, '\r');
  removeCharacterFromString(sequenceString, '\n');
  removeCharacterFromString(outputBasePath, '\r');
  removeCharacterFromString(outputBasePath, '\n');

  // derived parameters and cleanup:
  sellmeierType = 0;
  axesNumber = 0;
  Ntime = (int64_t)(minGridDimension *
                    round(timeSpan / (minGridDimension * tStep)));
  timeSpan = Ntime * tStep; // force timeSpan to be consistent with Ntime
  Nfreq = Ntime / 2 + 1;
  Nspace = (int64_t)(minGridDimension *
                     round(spatialWidth / (minGridDimension * rStep)));
  spatialWidth = Nspace * rStep;
  Nspace2 = (int64_t)(minGridDimension *
                      round(spatialHeight / (minGridDimension * rStep)));
  spatialHeight = Nspace2 * rStep;
  Ngrid = Ntime * Nspace;
  NgridC = Nfreq * Nspace;
  kStep = twoPi<double>() / (Nspace * rStep);
  fStep = 1.0 / (Ntime * tStep);
  Npropagation = (int64_t)round(crystalThickness / propagationStep);

  isCylindric = symmetryType == 1;
  is3D = symmetryType == 2 || symmetryType == 4;
  isFDTD = symmetryType == 3 || symmetryType == 4;
  if (isCylindric) {
    pulse1.beam_spec.x_offset[0][0] = 0.0;
    pulse2.beam_spec.x_offset[0][0] = 0.0;
    pulse1.beam_spec.angle_x[0][0] = 0.0;
    pulse2.beam_spec.angle_x[0][0] = 0.0;
  }
  if (is3D) {
    Ngrid = Ntime * Nspace * Nspace2;
    NgridC = Nfreq * Nspace * Nspace2;
  } else {
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

  // crystal from database (database must be loaded!)
  crystalDatabase = crystalDatabasePtr;
  chi2Tensor = crystalDatabasePtr[materialIndex].d.data();
  chi3Tensor = crystalDatabasePtr[materialIndex].chi3.data();
  nonlinearSwitches = crystalDatabasePtr[materialIndex].nonlinearSwitches;
  sellmeierCoefficients =
      crystalDatabasePtr[materialIndex].sellmeierCoefficients.data();
  sellmeierType = crystalDatabasePtr[materialIndex].sellmeierType;
  axesNumber = crystalDatabasePtr[materialIndex].axisType;

  if (fs.good())
    return 61;
  else
    return -1;
}

void simulationBatch::configure(bool allocateFields) {
  Nfreq = parameters[0].Nfreq;
  Nsims = parameters[0].Nsims;
  Nsims2 = parameters[0].Nsims2;
  Nsimstotal = Nsims * Nsims2;
  Ngrid = parameters[0].Ngrid;
  NgridC = parameters[0].NgridC;
  simulationParameterSet base = parameters[0];
  parameters.resize(Nsimstotal, base);
  std::for_each(mutexes.begin(), mutexes.end(), [](std::shared_mutex &m) {
    std::lock_guard<std::shared_mutex> lock(m);
  });
  if (allocateFields) {
    Ext = std::vector<double>(Nsimstotal * Ngrid * 2, 0.0);
    Ekw = std::vector<std::complex<double>>(Nsimstotal * NgridC * 2,
                                            std::complex<double>(0.0, 0.0));
  }

  mutexes = std::vector<std::shared_mutex>(Nsimstotal);
  totalSpectrum = std::vector<double>(Nfreq * Nsimstotal * 3);

  if (parameters[0].pulse1FileType == 1 || parameters[0].pulse1FileType == 2) {
    loadedField1 = std::vector<std::complex<double>>(
        Nfreq, std::complex<double>(0.0, 0.0));
  }
  if (parameters[0].pulse2FileType == 1 || parameters[0].pulse2FileType == 2) {
    loadedField2 = std::vector<std::complex<double>>(
        Nfreq, std::complex<double>(0.0, 0.0));
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

  // configure
  double step1 = (Nsims > 1) ? (parameters[0].batchDestination -
                                parameters[0].getByNumberWithMultiplier(
                                    parameters[0].batchIndex)) /
                                   (Nsims - 1)
                             : 0.0;
  double step2 = (Nsims2 > 1) ? (parameters[0].batchDestination2 -
                                 parameters[0].getByNumberWithMultiplier(
                                     parameters[0].batchIndex2)) /
                                    (Nsims2 - 1)
                              : 0.0;

  parameters[0].ExtOut = Ext.data();
  parameters[0].EkwOut = Ekw.data();
  parameters[0].totalSpectrum = totalSpectrum.data();
  parameters[0].loadedField1 = loadedField1.data();
  parameters[0].loadedField2 = loadedField2.data();
  parameters[0].loadedFullGrid1 = loadedFullGrid1.data();
  parameters[0].loadedFullGrid2 = loadedFullGrid2.data();
  parameters[0].fittingReference = fitReference.data();
  parameters[0].isGridAllocated = true;
  parameters[0].optics = optics;
  loadPulseFiles();

  for (int64_t i = 0; i < Nsims2; i++) {
    int64_t currentRow = i * Nsims;

    if (currentRow > 0) {
      parameters[currentRow] = parameters[0];
    }
    if (Nsims2 > 1) {
      parameters[currentRow].setByNumberWithMultiplier(
          parameters[0].batchIndex2,
          parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex2) +
              i * step2);
    }

    for (int64_t j = 0; j < Nsims; j++) {

      if (j > 0) {
        parameters[j + currentRow] = parameters[currentRow];
      }

      parameters[j + currentRow].batchLoc1 = j;
      parameters[j + currentRow].batchLoc2 = i;
      parameters[j + currentRow].ExtOut = getExt((j + currentRow));
      parameters[j + currentRow].EkwOut = getEkw((j + currentRow));
      parameters[j + currentRow].totalSpectrum =
          getTotalSpectrum((j + currentRow));
      parameters[j + currentRow].isFollowerInSequence = false;
      parameters[j + currentRow].cancellationCalled = false;
      parameters[j + currentRow].setByNumberWithMultiplier(
          parameters[0].batchIndex,
          parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex) +
              j * step1);
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
  std::for_each(mutexes.begin(), mutexes.end(), [](std::shared_mutex &m) {
    std::lock_guard<std::shared_mutex> lock(m);
  });
  mutexes = std::vector<std::shared_mutex>(Nsimstotal);

  // configure
  double step1 =
      (parameters[0].batchDestination -
       parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex)) /
      (Nsims - 1);
  double step2 = 0.0;
  if (Nsims2 > 0) {
    step2 =
        (parameters[0].batchDestination2 -
         parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex2)) /
        (Nsims2 - 1);
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
          parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex2) +
              i * step2);
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
          parameters[0].getByNumberWithMultiplier(parameters[0].batchIndex) +
              j * step1);
    }
  }
}

void simulationBatch::loadPulseFiles() {
  // pulse type specifies if something has to be loaded to describe the pulses,
  // or if they should be synthesized later. 1: FROG .speck format; 2:
  // Time-domain waveform; 3: Previous LWE result
  int frogLines = 0;
  if (parameters[0].pulse1FileType == 1) {
    frogLines =
        loadFrogSpeck(parameters[0].pulse1LoadedData, loadedField1.data(),
                      parameters[0].Ntime, parameters[0].fStep, 0.0);
    parameters[0].field1IsAllocated = (frogLines > 1);
  }
  if (parameters[0].pulse2FileType == 1) {
    frogLines =
        loadFrogSpeck(parameters[0].pulse2LoadedData, loadedField2.data(),
                      parameters[0].Ntime, parameters[0].fStep, 0.0);
    parameters[0].field2IsAllocated = (frogLines > 1);
  }

  if (parameters[0].pulse1FileType == 2) {
    frogLines =
        loadWaveformFile(parameters[0].pulse1LoadedData, loadedField1.data(),
                         parameters[0].Ntime, parameters[0].fStep);
    parameters[0].field1IsAllocated = (frogLines > 1);
  }
  if (parameters[0].pulse2FileType == 2) {
    frogLines =
        loadWaveformFile(parameters[0].pulse2LoadedData, loadedField2.data(),
                         parameters[0].Ntime, parameters[0].fStep);
    parameters[0].field2IsAllocated = (frogLines > 1);
  }

  if (parameters[0].pulse1FileType == 3) {
    parameters[0].field1IsAllocated = loadSavedGridFile(
        parameters[0].pulse1LoadedData, loadedFullGrid1, parameters[0].Ngrid);
  }
  if (parameters[0].pulse2FileType == 3) {
    parameters[0].field2IsAllocated = loadSavedGridFile(
        parameters[0].pulse2LoadedData, loadedFullGrid2, parameters[0].Ngrid);
  }
}
void simulationBatch::loadOptics(const std::string &zipPath) {
  bool keepLoading = true;
  std::string base = zipGetBasename(zipPath) + "_optic";
  int ind = 0;
  while (keepLoading) {
    std::string currentFile = base + std::to_string(ind) + ".txt";
    keepLoading = zipContainsFile(zipPath, currentFile);
    if (!keepLoading)
      break;
    if (ind == 0)
      optics.clear();
    loadedInputData file(zipPath, currentFile);
    optics.push_back(file);
    ind++;
  }
}

int simulationBatch::saveDataSet() {
  std::for_each(mutexes.begin(), mutexes.end(),
                [&](std::shared_mutex &m) { m.lock(); });

  std::string Epath = parameters[0].outputBasePath + "_Ext.dat";
  std::string Spath = parameters[0].outputBasePath + "_spectrum.dat";
  std::string Zpath = parameters[0].outputBasePath + ".zip";
  std::string Tpath = parameters[0].outputBasePath + ".txt";
  std::string FittingTargetPath =
      parameters[0].outputBasePath + "_fittingTarget.dat";
  std::string Pulse1Path = parameters[0].outputBasePath + "_pulse1.dat";
  std::string Pulse2Path = parameters[0].outputBasePath + "_pulse2.dat";
  std::string outputText = parameters[0].settingsString();
  mz_zip_archive zip = {};
  mz_zip_writer_init_file(&zip, Zpath.c_str(), 0);
  mz_zip_writer_add_mem(&zip, getBasename(Tpath).c_str(), outputText.c_str(),
                        outputText.size(), MZ_DEFAULT_COMPRESSION);
  mz_zip_writer_add_mem(&zip, getBasename(Epath).c_str(), Ext.data(),
                        sizeof(double) * Ext.size(), MZ_DEFAULT_COMPRESSION);
  mz_zip_writer_add_mem(&zip, getBasename(Spath).c_str(), totalSpectrum.data(),
                        sizeof(double) * totalSpectrum.size(),
                        MZ_DEFAULT_COMPRESSION);
  if (parameters[0].pulse1LoadedData.hasData) {
    mz_zip_writer_add_mem(&zip, getBasename(Pulse1Path).c_str(),
                          parameters[0].pulse1LoadedData.fileContents.c_str(),
                          parameters[0].pulse1LoadedData.fileContents.size(),
                          MZ_DEFAULT_COMPRESSION);
  }
  if (parameters[0].pulse2LoadedData.hasData) {
    mz_zip_writer_add_mem(&zip, getBasename(Pulse2Path).c_str(),
                          parameters[0].pulse2LoadedData.fileContents.c_str(),
                          parameters[0].pulse2LoadedData.fileContents.size(),
                          MZ_DEFAULT_COMPRESSION);
  }
  if (parameters[0].fittingLoadedData.hasData) {
    mz_zip_writer_add_mem(&zip, getBasename(FittingTargetPath).c_str(),
                          parameters[0].fittingLoadedData.fileContents.c_str(),
                          parameters[0].fittingLoadedData.fileContents.size(),
                          MZ_DEFAULT_COMPRESSION);
  }
  if (parameters[0].optics.size() > 0) {
    for (int i = 0; i < static_cast<int>(optics.size()); i++) {
      std::string opticPath =
          parameters[0].outputBasePath + "_optic" + std::to_string(i) + ".txt";
      mz_zip_writer_add_mem(
          &zip, getBasename(opticPath).c_str(), optics[i].fileContents.c_str(),
          optics[i].fileContents.size(), MZ_DEFAULT_COMPRESSION);
    }
  }
  mz_zip_writer_finalize_archive(&zip);
  mz_zip_writer_end(&zip);

  std::for_each(mutexes.begin(), mutexes.end(),
                [&](std::shared_mutex &m) { m.unlock(); });
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
    ss >> fittingArray[fittingCount] >> fittingArray[fittingCount + 1] >>
        fittingArray[fittingCount + 2];
    if (ss.good())
      fittingCount += 3;
    ss.ignore(fittingString.length(), ';');
  }

  Nfitting = fittingCount / 3;
  isInFittingMode = ((Nfitting > 0) && (maxIterations > 0));

  return 0;
}
