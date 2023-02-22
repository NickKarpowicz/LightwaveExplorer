import numpy as np
import re
import os
from scipy.optimize import least_squares
from scipy.special import dawsn
import yaml
import urllib.request

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
        self.Nspace2 = 0
        self.Nfreq = 0
        self.Ngrid = 0
        self.Nsims = 1
        self.Nsims2 = 1
        self.spatialWidth = 0
        self.timeSpan = 0
        self.materialIndex = 0
        self.materialIndexAlternate = 0
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
        self.y01 = 0
        self.y02 = 0
        self.x01 = 0
        self.x02 = 0
        self.propagationAngle1 = 0
        self.propagationAngle2 = 0
        self.propagationAnglePhi1 = 0
        self.propagationAnglePhi2 = 0
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
        self.phaseMaterialIndex1 = 0
        self.phaseMaterialIndex2 = 0
        self.polarizationAngle1 = 0
        self.polarizationAngle2 = 0
        self.circularity1 = 0
        self.circularity2 = 0
        self.pulsetype = 0
        self.Ext = 0
        self.batchIndex = 0
        self.batchIndex2 = 0
        self.batchDestination = 0
        self.batchDestination2 = 0
        self.batchStart = 0
        self.batchStart2 = 0
        self.batchVector = 0
        self.batchVector2 = 0
        self.Ext_x = 0
        self.Ext_y = 0
        self.spectrumTotal = 0
        self.spectrum_x = 0
        self.spectrum_y = 0
        self.frequencyVector = 0
        self.frequencyVectorSpectrum = 0
        self.timeVector = 0
        self.spaceVector = 0

def readLine(line: str):
    rr = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)
    return float(rr[-1])
    
#fwhm
# gives the full-width at (height) maximum (default 0.5)
# copied from Stackoverflow...
def fwhm(x: np.ndarray, y: np.ndarray, height: float = 0.5) -> float:
    heightLevel = np.max(y) * height
    indexMax = np.argmax(y)
    xLower = np.interp(heightLevel, y[:indexMax], x[:indexMax])
    xUpper = np.interp(heightLevel, np.flip(y[indexMax:]), np.flip(x[indexMax:]))
    return xUpper - xLower

def norma(v: np.ndarray):
    return v/v.max()

def printSellmeier(sc):
    s = np.array2string(sc, formatter={'float_kind': '{0:.6g}'.format}).replace('\n','').replace('[','').replace(']','')
    print(s)

def sellmeier(wavelengthMicrons, a, equationType: int):
    w = 2 * np.pi * 2.99792458e8 / (1e-6*wavelengthMicrons)
    ls = wavelengthMicrons**2
    k = 3182.607353999257
    def rectangularBand(w,w0,width,height):
        if width != 0:
            scaledAxis = (w-w0)/width
        else:
            scaledAxis = (w-w0)
        imagPart = np.zeros(w.size)
        imagPart[np.abs(scaledAxis)<0.5] = -height
        realPart = np.zeros(w.size) + 1.0
        realPart[scaledAxis!=0.5] = -height * np.log(np.abs((scaledAxis+0.5)/(scaledAxis-0.5)))/np.pi
        return realPart + 1j*imagPart
    def gaussianBand(w,w0,width,height):
        if width == 0.0:
            return 0 + 1j*0
        scaledF = (w-w0)/(np.sqrt(2.0) * width)
        realPart = -dawsn(scaledF)/np.sqrt(np.pi)
        imagPart = np.exp(-scaledF*scaledF)
        return height * (realPart - 1j * imagPart)

    if equationType == 0:
        n = (a[0] + (a[1] + a[2] * ls) / (ls + a[3]) + (a[4] + a[5] * ls) / (ls + a[6]) + (a[7] + a[8] * ls) / (ls + a[9]) + (a[10] + a[11] * ls) / (ls + a[12]) + a[13] * ls + a[14] * ls * ls + a[15] * ls * ls * ls) 
        n[n<0] = 1
        n = n + k * a[16] / ((a[17] - w ** 2) +  (a[18] * w) * 1j)
        n += k * a[19] / ((a[20] - w ** 2) +  (a[21] * w) * 1j)
    elif equationType == 1:
        n = 1.0
        n = k * a[0] / ((a[1] - w ** 2) +  (a[2] * w) * 1j)
        n += k * a[3] / ((a[4] - w ** 2) +  (a[5] * w) * 1j)
        n += k * a[6] / ((a[7] - w ** 2) +  (a[8] * w) * 1j)
        n += k * a[9] / ((a[10] - w ** 2) +  (a[11] * w) * 1j)
        n += k * a[12] / ((a[13] - w ** 2) +  (a[14] * w) * 1j)
        n += k * a[15] / ((a[16] - w ** 2) +  (a[17] * w) * 1j)
        n += k * a[18] / ((a[19] - w ** 2) +  (a[20] * w) * 1j)
    elif equationType == 2:
        n = a[0]
        for i in range(0,6):
            n += gaussianBand(w,a[i*3+1],a[2 + i*3],a[3+i*3])

    return np.sqrt(n)

def load(filePath: str, loadFieldArray=True):
    s = lightwaveExplorerResult()

    #read in the values from the settings file
    settingsFile = open(filePath, "r")
    lines = settingsFile.readlines()
    i = 0
    s.pulseEnergy1 = readLine(lines[i])
    i+=1
    s.pulseEnergy2 = readLine(lines[i])
    i+=1
    s.frequency1 = readLine(lines[i])
    i+=1
    s.frequency2 = readLine(lines[i])
    i+=1
    s.bandwidth1 = readLine(lines[i])
    i+=1
    s.bandwidth2 = readLine(lines[i])
    i+=1
    s.sgOrder1 = readLine(lines[i])
    i+=1
    s.sgOrder2 = readLine(lines[i])
    i+=1
    s.cephase1 = readLine(lines[i])
    i+=1
    s.cephase2 = readLine(lines[i])
    i+=1
    s.delay1 = readLine(lines[i])
    i+=1
    s.delay2 = readLine(lines[i])
    i+=1
    s.gdd1 = readLine(lines[i])
    i+=1
    s.gdd2 = readLine(lines[i])
    i+=1
    s.tod1 = readLine(lines[i])
    i+=1
    s.tod2 = readLine(lines[i])
    i+=1
    s.phaseMaterialIndex1 = readLine(lines[i])
    i+=1
    s.phaseMaterialIndex2 = readLine(lines[i])
    i+=1
    s.phaseMaterialThickness1 = readLine(lines[i])
    i+=1
    s.phaseMaterialThickness2 = readLine(lines[i])
    i+=2
    #NOTE THAT THERE IS A PLACEHOLDER HERE FOR EXTRA BEAM MODE
    s.beamwaist1 = readLine(lines[i])
    i+=1
    s.beamwaist2 = readLine(lines[i])
    i+=1
    s.x01 = readLine(lines[i])
    i+=1
    s.x02 = readLine(lines[i])
    i+=1
    s.y01 = readLine(lines[i])
    i+=1
    s.y02 = readLine(lines[i])
    i+=1
    s.z01 = readLine(lines[i])
    i+=1
    s.z02 = readLine(lines[i])
    i+=1
    s.propagationAngle1 = readLine(lines[i])
    i+=1
    s.propagationAngle2 = readLine(lines[i])
    i+=1
    s.propagationAnglePhi1 = readLine(lines[i])
    i+=1
    s.propagationAnglePhi2 = readLine(lines[i])
    i+=1
    s.polarizationAngle1 = readLine(lines[i])
    i+=1
    s.polarizationAngle2 = readLine(lines[i])
    i+=1
    s.circularity1 = readLine(lines[i])
    i+=1
    s.circularity2 = readLine(lines[i])
    i+=1
    s.materialIndex = readLine(lines[i])
    i+=1
    s.materialIndexAlternate = readLine(lines[i])
    i+=1
    s.crystalTheta = readLine(lines[i])
    i+=1
    s.crystalPhi = readLine(lines[i])
    i+=1
    s.spatialWidth = readLine(lines[i])
    i+=1
    s.spatialHeight = readLine(lines[i])
    i+=1
    s.rStep = readLine(lines[i])
    i+=1
    s.timeSpan = readLine(lines[i])
    i+=1
    s.tStep = readLine(lines[i])
    i+=1
    s.crystalThickness = readLine(lines[i])
    i+=1
    s.propagationStep = readLine(lines[i])
    i+=1
    s.nonlinearAbsorptionStrength = readLine(lines[i])
    i+=1
    s.bandGapElectronVolts = readLine(lines[i])
    i+=1
    s.effectiveMass = readLine(lines[i])
    i+=1
    s.drudeGamma = readLine(lines[i])
    i+=1
    s.symmetryType = readLine(lines[i])
    i+=1
    s.batchIndex = readLine(lines[i])
    i+=1
    s.batchDestination = readLine(lines[i])
    i+=1
    s.Nsims = int(readLine(lines[i]))
    i+=1
    testString = lines[i]
    if testString[0] == 'B':
        s.batchIndex2 = readLine(lines[i])
        i+=1
        s.batchDestination2 = readLine(lines[i])
        i+=1
        s.Nsims2 = int(readLine(lines[i]))
        i+=1
    MIN_GRIDDIM=8
    s.Ntime = int(MIN_GRIDDIM * np.round(s.timeSpan / (MIN_GRIDDIM * s.tStep)))
    s.Nfreq = int(s.Ntime/2 + 1)
    s.fStep = 1.0/(s.Ntime * s.tStep)
    s.Nspace = int(MIN_GRIDDIM * np.round(s.spatialWidth / (MIN_GRIDDIM * s.rStep)))
    if s.symmetryType == 2:
        s.Nspace2 = int(MIN_GRIDDIM * np.round(s.spatialHeight / (MIN_GRIDDIM * s.rStep)))
    else:
        s.Nspace2 = 1
    s.Ngrid = int(s.Ntime*s.Nspace*s.Nspace2)

    
    #now load the output data from binary format. Note that this will fail if you're using wrong-endian CPUs
    fileBase = os.path.splitext(filePath)
    if loadFieldArray: 
        if s.symmetryType == 2:
            s.Ext = np.reshape(np.fromfile(fileBase[0]+"_Ext.dat",dtype=np.double)[0:(2*s.Ngrid*s.Nsims*s.Nsims2)],(s.Ntime,s.Nspace,s.Nspace2, 2*s.Nsims,s.Nsims2),order='F')
            s.Ext_x = np.squeeze(s.Ext[:,:,:,0:(2*s.Nsims):2,:])
            s.Ext_y = np.squeeze(s.Ext[:,:,:,1:(2*s.Nsims + 1):2,:])
        else:
            s.Ext = np.reshape(np.fromfile(fileBase[0]+"_Ext.dat",dtype=np.double)[0:(2*s.Ngrid*s.Nsims*s.Nsims2)],(s.Ntime,s.Nspace, 2*s.Nsims, s.Nsims2),order='F')
            s.Ext_x = np.squeeze(s.Ext[:,:,0:(2*s.Nsims):2,:])
            s.Ext_y = np.squeeze(s.Ext[:,:,1:(2*s.Nsims + 1):2,:])
    
    s.spectrum = np.reshape(np.fromfile(fileBase[0]+"_spectrum.dat",dtype=np.double)[0:3*s.Nfreq*s.Nsims*s.Nsims2],(s.Nfreq,3,s.Nsims,s.Nsims2),order='F')
    s.spectrumTotal = np.squeeze(s.spectrum[:,2,:,:]).T
    s.spectrum_x = np.squeeze(s.spectrum[:,0,:,:]).T
    s.spectrum_y = np.squeeze(s.spectrum[:,1,:,:]).T
    
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
    elif s.batchIndex == 36:
        s.batchStart = 0
    elif s.batchIndex == 37:
        s.batchStart = 0
        
    #that was fun, wasn't it? I think we all just had a good time. In c++ I just added an offset to a pointer. Now make the scale vector.
    s.batchVector = np.linspace(s.batchStart,s.batchDestination,s.Nsims)

    #make scale vector corresponding to the batch scan and correct units of the scan
    if s.batchIndex2 == 0:
        s.batchStart2 = 0
    elif s.batchIndex2 == 1:
        s.batchStart2 = s.pulseEnergy1
    elif s.batchIndex2 == 2:
        s.batchStart2 = s.pulseEnergy2
    elif s.batchIndex2 == 3:
        s.batchStart2 = s.frequency1
        s.batchDestination2 *= 1e12
    elif s.batchIndex2 == 4:
        s.batchStart2 = s.frequency2
        s.batchDestination2 *= 1e12
    elif s.batchIndex2 == 5:
        s.batchStart2 = s.frequency1
        s.batchDestination2 *= 1e12
    elif s.batchIndex2 == 6:
        s.batchStart2 = s.frequency2
        s.batchDestination2 *= 1e12
    elif s.batchIndex2 == 7:
        s.batchStart2 = s.cephase1
        s.batchDestination2 *= np.pi
    elif s.batchIndex2 == 8:
        s.batchStart2 = s.cephase2
        s.batchDestination2 *= np.pi
    elif s.batchIndex2 == 9:
        s.batchStart2 = s.delay1
        s.batchDestination2 *= 1e-15
    elif s.batchIndex2 == 10:
        s.batchStart2 = s.delay2
        s.batchDestination2 *= 1e-15
    elif s.batchIndex2 == 11:
        s.batchStart2 = s.gdd1
        s.batchDestination2 *= 1e-30
    elif s.batchIndex2 == 12:
        s.batchStart2 = s.gdd2
        s.batchDestination2 *= 1e-30
    elif s.batchIndex2 == 13:
        s.batchStart2 = s.tod1
        s.batchDestination2 *= 1e-45
    elif s.batchIndex2 == 14:
        s.batchStart2 = s.tod2
        s.batchDestination2 *= 1e-45   
    elif s.batchIndex2 == 15:
        s.batchStart2 = s.phaseMaterialThickness1
        s.batchDestination2 *= 1e-6
    elif s.batchIndex2 == 16:
        s.batchStart2 = s.phaseMaterialThickness2
        s.batchDestination2 *= 1e-6
    elif s.batchIndex2 == 17:
        s.batchStart2 = s.beamwaist1
        s.batchDestination2 *= 1e-6
    elif s.batchIndex2 == 18:
        s.batchStart2 = s.beamwaist2
        s.batchDestination2 *= 1e-6
    elif s.batchIndex2 == 19:
        s.batchStart2 = s.x01
        s.batchDestination2 *= 1e-6
    elif s.batchIndex2 == 20:
        s.batchStart2 = s.x02
        s.batchDestination2 *= 1e-6
    elif s.batchIndex2 == 21:
        s.batchStart2 = s.z01
        s.batchDestination2 *= 1e-6
    elif s.batchIndex2 == 22:
        s.batchStart2 = s.z02
        s.batchDestination2 *= 1e-6
    elif s.batchIndex2 == 23:
        s.batchStart2 = s.propagationAngle1
        s.batchDestination2 *= (np.pi/180)
    elif s.batchIndex2 == 24:
        s.batchStart2 = s.propagationAngle2
        s.batchDestination2 *= (np.pi/180)
    elif s.batchIndex2 == 25:
        s.batchStart2 = s.polarizationAngle1
        s.batchDestination2 *=(np.pi/180)
    elif s.batchIndex2 == 26:
        s.batchStart2 = s.polarizationAngle2
        s.batchDestination2 *= (np.pi/180)
    elif s.batchIndex2 == 27:
        s.batchStart2 = s.circularity1
    elif s.batchIndex2 == 28:
        s.batchStart2 = s.circularity2
    elif s.batchIndex2 == 29:
        s.batchStart2 = s.crystalTheta
        s.batchDestination2 *= (np.pi/180)
    elif s.batchIndex2 == 30:
        s.batchStart2 = s.crystalPhi
        s.batchDestination2 *= (np.pi/180)
    elif s.batchIndex2 == 31:
        s.batchStart2 = s.nonlinearAbsorptionStrength
    elif s.batchIndex2 == 32:
        s.batchStart2 = s.drudeGamma
        s.batchDestination2 *= 1e12
    elif s.batchIndex2 == 33:
        s.batchStart2 = s.effectiveMass
    elif s.batchIndex2 == 34:
        s.batchStart2 = s.crystalThickness
        s.batchDestination2 *= 1e-6
    elif s.batchIndex2 == 35:
        s.batchStart2 = s.propagationStep
        s.batchDestination2 *= 1e-9
    elif s.batchIndex2 == 36:
        s.batchStart = 0
    elif s.batchIndex2 == 37:
        s.batchStart = 0
    #that was fun, wasn't it? I think we all just had a good time. In c++ I just added an offset to a pointer. Now make the scale vector.
    s.batchVector2 = np.linspace(s.batchStart2,s.batchDestination2,s.Nsims2)

    s.timeVector = s.tStep*np.arange(0,s.Ntime)
    s.frequencyVector = np.fft.fftfreq(s.Ntime, d=s.tStep)
    s.frequencyVectorSpectrum = s.frequencyVector[0:s.Nfreq]
    s.frequencyVectorSpectrum[-1] *= -1
    s.spaceVector = s.rStep * (np.arange(0,s.Nspace) - s.Nspace/2) + 0.25 * s.rStep
    return s

def getRII_object(url):
    riiurl = urllib.request.urlopen(url)
    return  yaml.safe_load(riiurl)

def getSellmeierFromRII(url):
    RIIobj = getRII_object(url)
    types = np.array([x["type"] for x in RIIobj["DATA"]])
    retrievedCoeffs = np.zeros(22)

    if 'formula 1' in types:
        firstTrue = ((types=='formula 2').cumsum().cumsum()==1).argmax()
        RIIcoeffs = np.fromstring(RIIobj["DATA"][firstTrue]["coefficients"], dtype=float, sep=' ')
        Nresonances = int(np.floor((RIIcoeffs.size-1)/2))
        retrievedCoeffs[0] = RIIcoeffs[0] + 1
        for i in range(0,min(7,Nresonances)):
            retrievedCoeffs[i*3 + 2] = RIIcoeffs[1 + 2*i]
            retrievedCoeffs[i*3 + 3] = -(RIIcoeffs[2 + 2*i]**2)
    elif 'formula 2' in types:
        firstTrue = ((types=='formula 2').cumsum().cumsum()==1).argmax()
        RIIcoeffs = np.fromstring(RIIobj["DATA"][firstTrue]["coefficients"], dtype=float, sep=' ')
        Nresonances = int(np.floor((RIIcoeffs.size-1)/2))
        retrievedCoeffs[0] = RIIcoeffs[0] + 1
        for i in range(0,min(7,Nresonances)):
            retrievedCoeffs[i*3 + 2] = RIIcoeffs[1 + 2*i]
            retrievedCoeffs[i*3 + 3] = -RIIcoeffs[2 + 2*i]
    elif 'formula 4' in types:
        firstTrue = ((types=='formula 4').cumsum().cumsum()==1).argmax()
        RIIcoeffs = np.fromstring(RIIobj["DATA"][firstTrue]["coefficients"], dtype=float, sep=' ')
        Nresonances = int(np.floor((RIIcoeffs.size-1)/4))
        retrievedCoeffs[0] = RIIcoeffs[0]
        for i in range(0,min(Nresonances,2)):
            if RIIcoeffs[2 + 4*i] == 0:
                retrievedCoeffs[i*3 + 1] = RIIcoeffs[1 + 4*i]
            else:
                retrievedCoeffs[i*3 + 2] = RIIcoeffs[1 + 4*i]
            retrievedCoeffs[i*3 + 3] = -(RIIcoeffs[3 + 4*i]**RIIcoeffs[4+4*i])
        if RIIcoeffs.size > 5:
            Npoly = int(np.floor(int(RIIcoeffs.size - 5)/2))
            for i in range(0,Npoly):
                if RIIcoeffs[5 + 2*i + 1] == 2.0:
                    retrievedCoeffs[13] = RIIcoeffs[5 + 2*i]
                if RIIcoeffs[5 + 2*i + 1] == 4.0:
                    retrievedCoeffs[14] = RIIcoeffs[5 + 2*i]
                if RIIcoeffs[5 + 2*i + 1] == 6.0:
                    retrievedCoeffs[15] = RIIcoeffs[5 + 2*i]
    else:
        print("Sorry, this formula hasn't been implemented yet.")
    return retrievedCoeffs

def getTabulatedDataFromRII(url):
    RIIobj = getRII_object(url)
    types = np.array([x["type"] for x in RIIobj["DATA"]])
    if 'tabulated nk' in types:
        firstTrue = ((types=='tabulated nk').cumsum().cumsum()==1).argmax()
        nkData = np.array([np.array(i.split(' ')).astype(float) for i in str(RIIobj["DATA"][firstTrue]["data"]).strip().split('\n')])
    elif 'tabulated k' in types:
        firstTrue = ((types=='tabulated k').cumsum().cumsum()==1).argmax()
        kData = np.array([np.array(i.split(' ')).astype(float) for i in str(RIIobj["DATA"][firstTrue]["data"]).strip().split('\n')])
        nkData = np.zeros((kData[:,0].size, 3))
        nkData[:,0] = kData[:,0]
        nkData[:,2] = kData[:,1]
        sellCoeffs = getSellmeierFromRII(url)
        nkData[:,1] = np.real(sellmeier(nkData[:,0],sellCoeffs,0))
    elif 'tabulated n' in types:
        firstTrue = ((types=='tabulated n').cumsum().cumsum()==1).argmax()
        nData = np.array([np.array(i.split(' ')).astype(float) for i in str(RIIobj["DATA"][firstTrue]["data"]).strip().split('\n')])
        nkData = np.zeros((kData[:,0].size, 3))
        nkData[:,0] = nData[:,0]
        nkData[:,1] = nData[:,1]
        
    return nkData
        
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

    #make everything numpy arrays
    bandpass = np.array(bandpass)
    filterTransmissionNanometers = np.array(filterTransmissionNanometers)
    detectorResponseNanometers = np.array(detectorResponseNanometers)

    #resolve the various filters
    if bandpass.any() != None:
        bandpassFilter = np.exp(-(s.frequencyVectorSpectrum-bandpass[0])**bandpass[2]/(2*bandpass[1]**bandpass[2]))
        totalResponse *= bandpassFilter

    if filterTransmissionNanometers.any() != None:
        sortIndicies = np.argsort(1e9*c/filterTransmissionNanometers[0,:])
        dataFrequencyAxis = np.array([1e9*c/filterTransmissionNanometers[0,sortIndicies], filterTransmissionNanometers[1,sortIndicies]])
        totalResponse  *= np.interp(s.frequencyVectorSpectrum, dataFrequencyAxis[0,:], dataFrequencyAxis[1,:])
    
    if detectorResponseNanometers.any() != None:
        sortIndicies = np.argsort(1e9*c/detectorResponseNanometers[0,:])
        dataFrequencyAxis = np.array([1e9*c/detectorResponseNanometers[0,sortIndicies], detectorResponseNanometers[1,sortIndicies]])
        totalResponse  *= np.interp(s.frequencyVectorSpectrum, dataFrequencyAxis[0,:], dataFrequencyAxis[1,:])

    #EOS signal is the integral of the difference between the two spectra, multiplied by the total spectral response
    if s.Nsims2>1:
        EOSsignal = np.array([np.sum((totalResponse*(np.squeeze(s.spectrum_x[i,:,:]-s.spectrum_y[i,:,:]))), axis=1) for i in range(s.Nsims2)])
    else:
        EOSsignal = np.sum((totalResponse*(s.spectrum_x-s.spectrum_y)), axis=1)
    return EOSsignal

def sellmeierFit(wavelengthMicrons, startingCoefficients, activeElements, eqnType: int, nTarget, fitImaginary: bool):
    def expandCoeffs(x):
        x1 = np.zeros(22)
        for i in range(0,np.size(startingCoefficients)):
            x1[activeElements[i]] = x[activeElements[i]]
        return x1
    def fun_nforx(x):
        return sellmeier(wavelengthMicrons, expandCoeffs(x), eqnType)
    def fun_residual(x):
        nx = fun_nforx(x)
        if fitImaginary:
            returnVals = np.append(np.real(nTarget - nx), np.imag(nTarget - nx))
        else:
            returnVals = np.real(nTarget - nx)
        return returnVals

    res = least_squares(fun_residual, startingCoefficients, gtol=None, xtol=None, ftol = 1e-12, max_nfev=16384)
    print("Resulting coefficients:")
    printSellmeier(expandCoeffs(res.x))
    return expandCoeffs(res.x), fun_nforx(res.x)
