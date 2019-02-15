#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 12:39:07 2018
This module is used for reading gObs and gBar data from a file and plotting it.
The outputs of the MDAR method in my Tester.java file should be inputs here.
@author: Virinchi Rallabhandi
"""
import matplotlib.pyplot as plt
import numpy as np
gBarFile = open("D:/University/ICRAR Studentship/Data/gBar23.txt", "r")
gBar = gBarFile.readlines()
gBarFile.close()
gObsFile = open("D:/University/ICRAR Studentship/Data/gObs23.txt", "r")
gObs = gObsFile.readlines()
gObsFile.close()
for i in range(len(gBar)):
    gBar[i] = float(gBar[i])
    gObs[i] = float(gObs[i])
mcGaughRange = []
mcGaughFit = []
unity= []
for i in range(1000):
    g = 10**-11.7 + 10**-12*float(i)
    mcGaughRange.append(np.log10(g))
    unity.append(np.log10(g))
    mcGaughFit.append(np.log10(g/(1 - np.e**(-1*(g/1.2/10**(-10))**0.5))))
plt.figure()
plt.plot(gBar, gObs, "co")
plt.plot(mcGaughRange, mcGaughFit, "k-")
plt.plot(mcGaughRange, unity, "b--")
plt.xlabel("log(gBar (m/s^2))")
plt.ylabel("log(gObs (m/s^2))")
plt.legend(["Points from simulated galaxies", "McGaugh et al. fit", "gObs = gBar"])
plt.title("MDAR compiled from many simulated galaxies")
plt.savefig("D:/University/ICRAR Studentship/Data/MDAR compiled 23")