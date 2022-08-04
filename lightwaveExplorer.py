import numpy as np
import re
import os

class lightwaveExplorerResult:
    def __init__(self):
        self.rStep = 0
        self.tStep = 0
        self.fStep = 0
        self.kStep = 0
        self.propagationStep = 0
        self.Npropagation = 0
        self.Ntime = 0
        self.Nspace = 0
        self.Ngrid = 0
        self.Nsims = 0
        self.spatialWidth = 0
        self.timeSpan = 0
        self.materialIndex = 0
        self.bandGapElectronVolts = 0
        self.effectiveMass = 0
        self.nonlinearAbsorptionStrength = 0
        self.drudeGamma = 0
        self.crystalTheta = 0
        self.crystalPhi = 0
        self.crystalThickness = 0
        self.beamwaist1 = 0
        self.beamwaist2 = 0
        self.z01 = 0
        self.z02 = 0
        self.x01 = 0
        self.x02 = 0
        self.propagationAngle1 = 0
        self.propagationAngle2 = 0
        self.symmetryType = 0
        self.pulseEnergy1 = 0
        self.pulseEnergy2 = 0
        self.frequency1 = 0
        self.frequency2 = 0
        self.sgOrder1 = 0
        self.sgOrder2 = 0
        self.bandwidth1 = 0
        self.bandwidth2 = 0
        self.cephase1 = 0
        self.cephase2 = 0
        self.delay1 = 0
        self.delay2 = 0
        self.gdd1 = 0
        self.gdd2 = 0
        self.tod1 = 0
        self.tod2 = 0
        self.phaseMaterialThickness1 = 0
        self.phaseMaterialThickness2 = 0
        self.phaseMaterialIndex = 0
        self.polarizationAngle1 = 0
        self.polarizationAngle2 = 0
        self.circularity1 = 0
        self.circularity2 = 0
        self.pulsetype = 0
        self.Ext = 0
        self.batchIndex = 0
        self.batchDestination = 0
        self.batchStart = 0
        self.batchVector = 0
        self.Ext_x = 0
        self.Ext_y = 0
        self.spectrumTotal = 0
        self.spectrum_x = 0
        self.spectrum_y = 0
        self.frequencyVector = 0
        self.timeVector = 0
        self.spaceVector = 0

def readLine(line: str):
    rr = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)
    return float(rr[-1])
    

def load(filePath: str):
    s = lightwaveExplorerResult()

    #read in the values from the settings file
    settingsFile = open(filePath, "r")
    lines = settingsFile.readlines()
    s.pulseEnergy1 = readLine(lines[0])
    s.pulseEnergy2 = readLine(lines[1])
    s.frequency1 = readLine(lines[2])
    s.frequency2 = readLine(lines[3])
    s.bandwidth1 = readLine(lines[4])
    s.bandwidth2 = readLine(lines[5])
    s.sgOrder1 = readLine(lines[6])
    s.sgOrder2 = readLine(lines[6])
    s.cephase1 = readLine(lines[7])
    s.cephase2 = readLine(lines[8])
    s.delay1 = readLine(lines[9])
    s.delay2 = readLine(lines[10])
    s.gdd1 = readLine(lines[11])
    s.gdd2 = readLine(lines[12])
    s.tod1 = readLine(lines[13])
    s.tod2 = readLine(lines[14])
    s.phaseMaterialIndex = readLine(lines[15])
    s.phaseMaterialThickness1 = readLine(lines[16])
    s.phaseMaterialThickness2 = readLine(lines[17])
    s.beamwaist1 = readLine(lines[18])
    s.beamwaist2 = readLine(lines[19])
    s.x01 = readLine(lines[20])
    s.x02 = readLine(lines[21])
    s.z01 = readLine(lines[22])
    s.z02 = readLine(lines[23])
    s.propagationAngle1 = readLine(lines[24])
    s.propagationAngle2 = readLine(lines[25])
    s.polarizationAngle1 = readLine(lines[26])
    s.polarizationAngle2 = readLine(lines[27])
    s.circularity1 = readLine(lines[28])
    s.circularity2 = readLine(lines[29])
    s.materialIndex = readLine(lines[30])
    s.crystalTheta = readLine(lines[31])
    s.crystalPhi = readLine(lines[32])
    s.spatialWidth = readLine(lines[33])
    s.rStep = readLine(lines[34])
    s.timeSpan = readLine(lines[35])
    s.tStep = readLine(lines[36])
    s.crystalThickness = readLine(lines[37])
    s.propagationStep = readLine(lines[38])
    s.nonlinearAbsorptionStrength = readLine(lines[39])
    s.bandGapElectronVolts = readLine(lines[40])
    s.effectiveMass = readLine(lines[41])
    s.drudeGamma = readLine(lines[42])
    s.symmetryType = readLine(lines[43])
    s.batchIndex = readLine(lines[44])
    s.batchDestination = readLine(lines[45])
    s.Nsims = int(readLine(lines[46]))
    s.Ntime = int(np.round(s.timeSpan/s.tStep))
    s.Nspace = int(np.round(s.spatialWidth/s.rStep))
    s.Ngrid = int(s.Ntime*s.Nspace)
    
    #now load the output data from binary format. Note that this will fail if you're using wrong-endian CPUs
    fileBase = os.path.splitext(filePath)
    s.Ext = np.reshape(np.fromfile(fileBase[0]+"_Ext.dat",dtype=np.double)[0:(2*s.Ngrid*s.Nsims)],(s.Ntime,s.Nspace,2*s.Nsims),order='F')
    s.Ext_x = np.squeeze(s.Ext[:,:,0:(2*s.Nsims):2])
    s.Ext_y = np.squeeze(s.Ext[:,:,1:(2*s.Nsims + 1):2])
    
    s.spectrum = np.reshape(np.fromfile(fileBase[0]+"_spectrum.dat",dtype=np.double)[0:3*s.Ntime*s.Nsims],(s.Ntime,3,s.Nsims),order='F')
    s.spectrumTotal = np.squeeze(s.spectrum[0:int(s.Ntime/2),2,:]).T
    s.spectrum_x = np.squeeze(s.spectrum[0:int(s.Ntime/2),0,:]).T
    s.spectrum_y = np.squeeze(s.spectrum[0:int(s.Ntime/2),1,:]).T
    
    #make scale vector corresponding to the batch scan and correct units of the scan
    if s.batchIndex == 0:
        s.batchStart = 0
    elif s.batchIndex == 1:
        s.batchStart = s.pulseEnergy1
    elif s.batchIndex == 2:
        s.batchStart = s.pulseEnergy2
    elif s.batchIndex == 3:
        s.batchStart = s.frequency1
        s.batchDestination *= 1e12
    elif s.batchIndex == 4:
        s.batchStart = s.frequency2
        s.batchDestination *= 1e12
    elif s.batchIndex == 5:
        s.batchStart = s.frequency1
        s.batchDestination *= 1e12
    elif s.batchIndex == 6:
        s.batchStart = s.frequency2
        s.batchDestination *= 1e12
    elif s.batchIndex == 7:
        s.batchStart = s.cephase1
        s.batchDestination *= np.pi
    elif s.batchIndex == 8:
        s.batchStart = s.cephase2
        s.batchDestination *= np.pi
    elif s.batchIndex == 9:
        s.batchStart = s.delay1
        s.batchDestination *= 1e-15
    elif s.batchIndex == 10:
        s.batchStart = s.delay2
        s.batchDestination *= 1e-15
    elif s.batchIndex == 11:
        s.batchStart = s.gdd1
        s.batchDestination *= 1e-30
    elif s.batchIndex == 12:
        s.batchStart = s.gdd2
        s.batchDestination *= 1e-30
    elif s.batchIndex == 13:
        s.batchStart = s.tod1
        s.batchDestination *= 1e-45
    elif s.batchIndex == 14:
        s.batchStart = s.tod2
        s.batchDestination *= 1e-45   
    elif s.batchIndex == 15:
        s.batchStart = s.phaseMaterialThickness1
        s.batchDestination *= 1e-6
    elif s.batchIndex == 16:
        s.batchStart = s.phaseMaterialThickness2
        s.batchDestination *= 1e-6
    elif s.batchIndex == 17:
        s.batchStart = s.beamwaist1
        s.batchDestination *= 1e-6
    elif s.batchIndex == 18:
        s.batchStart = s.beamwaist2
        s.batchDestination *= 1e-6
    elif s.batchIndex == 19:
        s.batchStart = s.x01
        s.batchDestination *= 1e-6
    elif s.batchIndex == 20:
        s.batchStart = s.x02
        s.batchDestination *= 1e-6
    elif s.batchIndex == 21:
        s.batchStart = s.z01
        s.batchDestination *= 1e-6
    elif s.batchIndex == 22:
        s.batchStart = s.z02
        s.batchDestination *= 1e-6
    elif s.batchIndex == 23:
        s.batchStart = s.propagationAngle1
        s.batchDestination *= (np.pi/180)
    elif s.batchIndex == 24:
        s.batchStart = s.propagationAngle2
        s.batchDestination *= (np.pi/180)
    elif s.batchIndex == 25:
        s.batchStart = s.polarizationAngle1
        s.batchDestination *=(np.pi/180)
    elif s.batchIndex == 26:
        s.batchStart = s.polarizationAngle2
        s.batchDestination *= (np.pi/180)
    elif s.batchIndex == 27:
        s.batchStart = s.circularity1
    elif s.batchIndex == 28:
        s.batchStart = s.circularity2
    elif s.batchIndex == 29:
        s.batchStart = s.crystalTheta
        s.batchDestination *= (np.pi/180)
    elif s.batchIndex == 30:
        s.batchStart = s.crystalPhi
        s.batchDestination *= (np.pi/180)
    elif s.batchIndex == 31:
        s.batchStart = s.nonlinearAbsorptionStrength
    elif s.batchIndex == 32:
        s.batchStart = s.drudeGamma
        s.batchDestination *= 1e12
    elif s.batchIndex == 33:
        s.batchStart = s.effectiveMass
    elif s.batchIndex == 34:
        s.batchStart = s.crystalThickness
        s.batchDestination *= 1e-6
    elif s.batchIndex == 35:
        s.batchStart = s.propagationStep
        s.batchDestination *= 1e-9
        
    #that was fun, wasn't it? I think we all just had a good time. In c++ I just added an offset to a pointer. Now make the scale vector.
    s.batchVector = np.linspace(s.batchStart,s.batchDestination,s.Nsims)

    s.timeVector = s.tStep*np.arange(0,s.Ntime)
    s.frequencyVector = np.fft.fftfreq(s.Ntime, d=s.tStep)
    s.frequencyVectorSpectrum = s.frequencyVector[0:int(s.Ntime/2)]
    s.spaceVector = s.rStep * (np.arange(0,s.Nspace) - s.Nspace/2) + 0.25 * s.rStep
    return s

#Takes a result of a calculation and extracts the electro-optic sampling signal
#The simulation must be done in a way that includes mixing of the signal and local oscillator
#either by fully simulating the waveplate, or by rotating the system of coordinates in
#a sequence.
# Example sequence:
#   0 7 -1 -1 0 0 0 0 -1 -1 45;
#   0 0 -1 -1 0 0 0 0 0 -1 0;
# This propagates through crystal with index 7 with thickness, angles and step defined by the interface,
# then rotates the system by 45 degrees. The last 0-length propagation through vacuum re-evaluates the spectra.
# The spectral response of the detection system can be adjusted in three ways:
# 1: application of a supergaussian bandpass filter with a three element list:
#   bandpass=[centralFrequency, sigma, superGaussianOrder]
# 2: application of a measured transmission of the filters and beamline optics:
#   filterTransmissionNanometers = np.array([wavelengthsInNanometers, transmissionIntensity])
# 3: application of the detector's spectral response (assumed identical for the two photodiodes)
#   same format as 2
def EOS(s: lightwaveExplorerResult, bandpass=None, filterTransmissionNanometers=None, detectorResponseNanometers=None):
    c = 2.99792458e8
    totalResponse = 1.

    bandpass = np.array(bandpass)
    filterTransmissionNanometers = np.array(filterTransmissionNanometers)
    detectorResponseNanometers = np.array(detectorResponseNanometers)

    if bandpass.any() != None:
        bandpassFilter = np.exp(-(s.frequencyVectorSpectrum-bandpass[0])**bandpass[2]/(2*bandpass[1]**bandpass[2]))
        totalResponse *= bandpassFilter

    if filterTransmissionNanometers.any() != None:
        sortIndicies = np.argsort(1e9*c/filterTransmissionNanometers[0,:])
        dataFrequencyAxis = np.array([1e9*c/filterTransmissionNanometers[sortIndicies,0], filterTransmissionNanometers[sortIndicies,1]])
        totalResponse  *= np.interp(s.frequencyVectorSpectrum, dataFrequencyAxis[0,:], dataFrequencyAxis[1,:])
    
    if detectorResponseNanometers.any() != None:
        sortIndicies = np.argsort(1e9*c/detectorResponseNanometers[0,:])
        dataFrequencyAxis = np.array([1e9*c/detectorResponseNanometers[0,sortIndicies], detectorResponseNanometers[1,sortIndicies]])
        totalResponse  *= np.interp(s.frequencyVectorSpectrum, dataFrequencyAxis[0,:], dataFrequencyAxis[1,:])

    EOSsignal = np.sum((totalResponse*(s.spectrum_x-s.spectrum_y)), axis=1)
    return EOSsignal;