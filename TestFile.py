#!/usr/bin/python
import numpy as np
dt = 7.000000e-16
dz = 2.500000e-08
dx = 6.200000e-07
TestFile_ExtIn = np.reshape(np.fromfile("TestFile_ExtIn.dat",dtype=np.double)[0:327680],(512,320,2),order='F')
TestFile_ExtOut = np.reshape(np.fromfile("TestFile_ExtOut.dat",dtype=np.double)[0:327680],(512,320,2),order='F')
