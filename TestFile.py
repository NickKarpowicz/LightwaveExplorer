#!/usr/bin/python
import numpy as np
dt = 5.000000e-16
dz = 2.500000e-08
dx = 1.000000e-06
TestFile_ExtIn = np.reshape(np.fromfile("TestFile_ExtIn.dat",dtype=np.double)[0:1048576],(1024,512,2),order='F')
TestFile_ExtOut = np.reshape(np.fromfile("TestFile_ExtOut.dat",dtype=np.double)[0:1048576],(1024,512,2),order='F')
