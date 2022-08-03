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

def readLine(line: str):
    rr = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)
    return float(rr[-1])
    

def load(filePath: str):
    loadedData = lightwaveExplorerResult()
    settingsFile = open(filePath, "r")
    lines = settingsFile.readlines()
    loadedData.pulseEnergy1 = readLine(lines[0])
    loadedData.pulseEnergy2 = readLine(lines[1])
    loadedData.frequency1 = readLine(lines[2])
    loadedData.frequency2 = readLine(lines[3])
    loadedData.bandwidth1 = readLine(lines[4])
    loadedData.bandwidth2 = readLine(lines[5])
    loadedData.sgOrder1 = readLine(lines[6])
    loadedData.sgOrder2 = readLine(lines[6])
    loadedData.cephase1 = readLine(lines[7])
    loadedData.cephase2 = readLine(lines[8])
    loadedData.delay1 = readLine(lines[9])
    loadedData.delay2 = readLine(lines[10])
    loadedData.gdd1 = readLine(lines[11])
    loadedData.gdd2 = readLine(lines[12])
    loadedData.tod1 = readLine(lines[13])
    loadedData.tod2 = readLine(lines[14])
    loadedData.phaseMaterialIndex = readLine(lines[15])
    loadedData.phaseMaterialThickness1 = readLine(lines[16])
    loadedData.phaseMaterialThickness2 = readLine(lines[17])
    loadedData.beamwaist1 = readLine(lines[18])
    loadedData.beamwaist2 = readLine(lines[19])
    loadedData.x01 = readLine(lines[20])
    loadedData.x02 = readLine(lines[21])
    loadedData.z01 = readLine(lines[22])
    loadedData.z02 = readLine(lines[23])
    loadedData.propagationAngle1 = readLine(lines[24])
    loadedData.propagationAngle2 = readLine(lines[25])
    loadedData.polarizationAngle1 = readLine(lines[26])
    loadedData.polarizationAngle2 = readLine(lines[27])
    loadedData.circularity1 = readLine(lines[28])
    loadedData.circularity2 = readLine(lines[29])
    loadedData.materialIndex = readLine(lines[30])
    loadedData.crystalTheta = readLine(lines[31])
    loadedData.crystalPhi = readLine(lines[32])
    loadedData.spatialWidth = readLine(lines[33])
    loadedData.rStep = readLine(lines[34])
    loadedData.timeSpan = readLine(lines[35])
    loadedData.tStep = readLine(lines[36])
    loadedData.crystalThickness = readLine(lines[37])
    loadedData.propagationStep = readLine(lines[38])
    loadedData.nonlinearAbsorptionStrength = readLine(lines[39])
    loadedData.bandGapElectronVolts = readLine(lines[40])
    loadedData.effectiveMass = readLine(lines[41])
    loadedData.drudeGamma = readLine(lines[42])
    loadedData.symmetryType = readLine(lines[43])
    loadedData.batchIndex = readLine(lines[44])
    loadedData.batchDestination = readLine(lines[45])
    loadedData.Nsims = int(readLine(lines[46]))
    loadedData.Ntime = int(np.round(loadedData.timeSpan/loadedData.tStep))
    loadedData.Nspace = int(np.round(loadedData.spatialWidth/loadedData.rStep))
    loadedData.Ngrid = int(loadedData.Ntime*loadedData.Nspace)
    
    #now load the output data from binary format. Note that this will fail if you're using wrong-endian CPUs
    fileBase = os.path.splitext(filePath)
    loadedData.Ext = np.reshape(np.fromfile(fileBase[0]+"_Ext.dat",dtype=np.double)[0:(2*loadedData.Ngrid*loadedData.Nsims)],(loadedData.Ntime,loadedData.Nspace,2*loadedData.Nsims),order='F')
    loadedData.Ext_x = np.squeeze(loadedData.Ext[:,:,0:2:-1])
    loadedData.Ext_y = np.squeeze(loadedData.Ext[:,:,1:2:-1])
    
    loadedData.spectrum = np.reshape(np.fromfile(fileBase[0]+"_spectrum.dat",dtype=np.double)[0:3*loadedData.Ntime*loadedData.Nsims],(loadedData.Ntime,3,loadedData.Nsims),order='F')
    loadedData.spectrumTotal = np.squeeze(loadedData.spectrum[:,2,:])
    loadedData.spectrum_x = np.squeeze(loadedData.spectrum[:,0,:])
    loadedData.spectrum_y = np.squeeze(loadedData.spectrum[:,1,:])
    
    #make scale vector corresponding to the batch scan and correct units of the scan
    if loadedData.batchIndex == 0:
        loadedData.batchStart = 0
    elif loadedData.batchIndex == 1:
        loadedData.batchStart = loadedData.pulseEnergy1
    elif loadedData.batchIndex == 2:
        loadedData.batchStart = loadedData.pulseEnergy2
    elif loadedData.batchIndex == 3:
        loadedData.batchStart = loadedData.frequency1
        loadedData.batchDestination *= 1e12
    elif loadedData.batchIndex == 4:
        loadedData.batchStart = loadedData.frequency2
        loadedData.batchDestination *= 1e12
    elif loadedData.batchIndex == 5:
        loadedData.batchStart = loadedData.frequency1
        loadedData.batchDestination *= 1e12
    elif loadedData.batchIndex == 6:
        loadedData.batchStart = loadedData.frequency2
        loadedData.batchDestination *= 1e12
    elif loadedData.batchIndex == 7:
        loadedData.batchStart = loadedData.cephase1
        loadedData.batchDestination *= np.pi
    elif loadedData.batchIndex == 8:
        loadedData.batchStart = loadedData.cephase2
        loadedData.batchDestination *= np.pi
    elif loadedData.batchIndex == 9:
        loadedData.batchStart = loadedData.delay1
        loadedData.batchDestination *= 1e-15
    elif loadedData.batchIndex == 10:
        loadedData.batchStart = loadedData.delay2
        loadedData.batchDestination *= 1e-15
    elif loadedData.batchIndex == 11:
        loadedData.batchStart = loadedData.gdd1
        loadedData.batchDestination *= 1e-30
    elif loadedData.batchIndex == 12:
        loadedData.batchStart = loadedData.gdd2
        loadedData.batchDestination *= 1e-30
    elif loadedData.batchIndex == 13:
        loadedData.batchStart = loadedData.tod1
        loadedData.batchDestination *= 1e-45
    elif loadedData.batchIndex == 14:
        loadedData.batchStart = loadedData.tod2
        loadedData.batchDestination *= 1e-45   
    elif loadedData.batchIndex == 15:
        loadedData.batchStart = loadedData.phaseMaterialThickness1
        loadedData.batchDestination *= 1e-6
    elif loadedData.batchIndex == 16:
        loadedData.batchStart = loadedData.phaseMaterialThickness2
        loadedData.batchDestination *= 1e-6
    elif loadedData.batchIndex == 17:
        loadedData.batchStart = loadedData.beamwaist1
        loadedData.batchDestination *= 1e-6
    elif loadedData.batchIndex == 18:
        loadedData.batchStart = loadedData.beamwaist2
        loadedData.batchDestination *= 1e-6
    elif loadedData.batchIndex == 19:
        loadedData.batchStart = loadedData.x01
        loadedData.batchDestination *= 1e-6
    elif loadedData.batchIndex == 20:
        loadedData.batchStart = loadedData.x02
        loadedData.batchDestination *= 1e-6
    elif loadedData.batchIndex == 21:
        loadedData.batchStart = loadedData.z01
        loadedData.batchDestination *= 1e-6
    elif loadedData.batchIndex == 22:
        loadedData.batchStart = loadedData.z02
        loadedData.batchDestination *= 1e-6
    elif loadedData.batchIndex == 23:
        loadedData.batchStart = loadedData.propagationAngle1
        loadedData.batchDestination *= (np.pi/180)
    elif loadedData.batchIndex == 24:
        loadedData.batchStart = loadedData.propagationAngle2
        loadedData.batchDestination *= (np.pi/180)
    elif loadedData.batchIndex == 25:
        loadedData.batchStart = loadedData.polarizationAngle1
        loadedData.batchDestination *=(np.pi/180)
    elif loadedData.batchIndex == 26:
        loadedData.batchStart = loadedData.polarizationAngle2
        loadedData.batchDestination *= (np.pi/180)
    elif loadedData.batchIndex == 27:
        loadedData.batchStart = loadedData.circularity1
    elif loadedData.batchIndex == 28:
        loadedData.batchStart = loadedData.circularity2
    elif loadedData.batchIndex == 29:
        loadedData.batchStart = loadedData.crystalTheta
        loadedData.batchDestination *= (np.pi/180)
    elif loadedData.batchIndex == 30:
        loadedData.batchStart = loadedData.crystalPhi
        loadedData.batchDestination *= (np.pi/180)
    elif loadedData.batchIndex == 31:
        loadedData.batchStart = loadedData.nonlinearAbsorptionStrength
    elif loadedData.batchIndex == 32:
        loadedData.batchStart = loadedData.drudeGamma
        loadedData.batchDestination *= 1e12
    elif loadedData.batchIndex == 33:
        loadedData.batchStart = loadedData.effectiveMass
    elif loadedData.batchIndex == 34:
        loadedData.batchStart = loadedData.crystalThickness
        loadedData.batchDestination *= 1e-6
    elif loadedData.batchIndex == 35:
        loadedData.batchStart = loadedData.propagationStep
        loadedData.batchDestination *= 1e-9
        
        
        
    loadedData.batchVector = np.linspace(loadedData.batchStart,loadedData.batchDestination,loadedData.Nsims)
    return loadedData
    