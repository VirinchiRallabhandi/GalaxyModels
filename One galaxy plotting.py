#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:03:08 2018
This module is used for plotting various features, e.g. gas density profile,
rotation curve, halo density profile etc., from a single model galaxy.
The output of the oneGalaxy method of Tester.java should be an input here. 
@author: Virinchi Rallabhandi
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
dataFile = open("D:/University/ICRAR Studentship/Data/215a.txt", "r")
G = 6.67408*10**-11
v0 = 310000.0
z0 = (0.45*v0/100000.0 - 0.14)*3.086*10**19
#This will unfortunately have to be supplied manually everytime because I
#can't be bothered recording extra information to file
data = dataFile.readlines()
dataFile.close()
radius = []
logRadius = []
logGasDensity = []
gasDensity = []
stellarDensity = []
baryonDensity = []
logBaryonDensity = []
haloDensity = []
logHaloDensity = []
rotationCurve = []
accelerationS = []
accelerationG = []
velocityS = []
velocityG = []
velocityH = []
velocityB = []
gBar = []
gObs = []
sigmaSR = []
for line in data:
    numbers = line.split()
    radius.append(float(numbers[0]))
    gasDensity.append(float(numbers[2]))
    stellarDensity.append(float(numbers[1]))
    baryonDensity.append(float(numbers[1]) + float(numbers[2]))
    sigmaSR.append(1.67*np.sqrt(2*np.pi*G*float(numbers[1])*z0))
    haloDensity.append(float(numbers[3]))
    if float(numbers[0]) > 0 and float(numbers[2]) > 0 and float(numbers[3]) > 0:
        logRadius.append(np.log10(float(numbers[0])))
        logGasDensity.append(np.log10(float(numbers[2])))
        logHaloDensity.append(np.log10(float(numbers[3])))
        logBaryonDensity.append(np.log10(float(numbers[1]) + float(numbers[2])))
    accelerationS.append(float(numbers[4]))
    accelerationG.append(float(numbers[5]))
    gb = float(numbers[6])
    go = float(numbers[8])
    if gb > 0 and go > 0:
        gBar.append(np.log10(gb))
        gObs.append(np.log10(go))
    velocityS.append(np.sqrt(float(numbers[4])*float(numbers[0])))
    vg = float(numbers[5])*float(numbers[0])
    if vg > 0:
        velocityG.append(np.sqrt(vg))
    else:
        velocityG.append(0.0)
    velocityB.append(np.sqrt(float(numbers[6])*float(numbers[0])))
    rotationCurve.append(np.sqrt(float(numbers[8])*float(numbers[0])))
    vh = float(numbers[7])*float(numbers[0])
    if vh > 0:
        velocityH.append(np.sqrt(vh))
    else:
        velocityH.append(0.0)
mcGaughRange = []
mcGaughFit = []
unity= []
for i in range(1000):
    g = 10**-11.6 + 10**-12*float(i)
    mcGaughRange.append(np.log10(g))
    unity.append(np.log10(g))
    mcGaughFit.append(np.log10(g/(1 - np.e**(-1*(g/1.2/10**(-10))**0.5))))
"""The dark matter density at a radius, r, is rho = rho0/(1 + (r/a)^2)
   This method is used for fitting the deduced dark matter halo"""
def coredIsothermalSphere(r, rho0, a):
    return rho0/(1 + (r/a)**2)
haloForFit = np.zeros(len(haloDensity) - 100)
radiiForFit = np.zeros(len(radius) - 100)
logHaloForFit = np.zeros(len(logHaloDensity) - 100)
logRadiiForFit = np.zeros(len(logRadius) - 100)
for i in range(len(haloForFit)):
    haloForFit[i] = haloDensity[i + 50]*10**20
    radiiForFit[i] = radius[i + 50]/10**20
for i in range(len(logHaloForFit)):
    logHaloForFit[i] = logHaloDensity[i + 100]
    logRadiiForFit[i] = logRadius[i + 100]
fit = np.polyfit(logRadiiForFit, logHaloForFit, 1)
logFitHalo = np.zeros(len(logRadius))
for i in range(len(logFitHalo)):
    logFitHalo[i] = fit[0]*logRadius[i] + fit[1]
coreFit, covariance = opt.curve_fit(coredIsothermalSphere, radiiForFit, haloForFit)
rho0Core = coreFit[0]/10**20
aCore = coreFit[1]*10**20
corePoints = np.zeros(len(radius))
for i in range(len(radius)):
    corePoints[i] = coredIsothermalSphere(radius[i], rho0Core, aCore)
rho0Core = float(int(rho0Core*10**26))/10**26
aCore = float(int(aCore/10**15))*10**15
radiiForPlot = np.zeros(len(radius) - 50)
haloForPlot = np.zeros(len(radiiForPlot))
coreForPlot = np.zeros(len(radiiForPlot))
for i in range(len(radiiForPlot)):
    radiiForPlot[i] = radius[i + 50]
    haloForPlot[i] = haloDensity[i + 50]
    coreForPlot[i] = corePoints[i + 50]
for i in range(int(len(gBar)*0.25)):
    del gBar[0]
    del gObs[0]
    del gObs[len(gObs) - 1]
    del gBar[len(gBar) - 1]
plt.figure()
plt.plot(radius, sigmaSR)
plt.xlabel("Radius (m)")
plt.ylabel("Velocity (m/s)")
plt.title("Deduced stellar radial velocity dispersion")
#The stellar radial velocity dispersion plot should be ignored. It contains and
#old version of the stellar radial velocity equation.
plt.figure()
plt.plot(radius, gasDensity)
plt.xlabel("Radius (m)")
plt.ylabel("Gas density (kg/m^2)")
plt.title("Deduced gas surface mass density")
plt.figure()
plt.plot(logRadius, logGasDensity)
plt.xlabel("log(radius (m))")
plt.ylabel("log(gas density (kg/m^2))")
plt.title("Deduced gas surface mass density")
plt.figure()
plt.plot(radius, stellarDensity)
plt.xlabel("Radius (m)")
plt.ylabel("Stellar density (kg/m^2)")
plt.title("Stellar surface mass density")
plt.figure()
plt.plot(radius, baryonDensity)
plt.xlabel("Radius (m)")
plt.ylabel("Baryon density (kg/m^2)")
plt.title("Total baryon density")
plt.figure()
plt.plot(logRadius, logBaryonDensity)
plt.xlabel("log(radius (m))")
plt.ylabel("log(baryon density (kg/m^2))")
plt.title("Total baryon density")
plt.figure()
plt.plot(radius, accelerationS)
plt.xlabel("Radius (m)")
plt.ylabel("Acceleration (m/s^2)")
plt.title("Accleration due to the stellar disk")
plt.figure()
plt.plot(radius, accelerationG)
plt.xlabel("Radius (m)")
plt.ylabel("Acceleration (m/s^2)")
plt.title("Acceleration due to the gas disk")
plt.figure()
plt.plot(radius, haloDensity, "bo")
plt.plot(radius, corePoints, "r-")
plt.xlabel("Radius (m)")
plt.ylabel("Density (kg/m^3)")
plt.legend(["Deduced from model"
            , str(rho0Core)+"/"+"(1 + (r/"+str(aCore)+")^2)"])
plt.title("Dark matter density from rotation curve")
#plt.savefig("/Users/virinchirallabhandi/Desktop/Dark matter density 4")
plt.figure()
plt.plot(logRadius, logHaloDensity, "b-")
plt.plot(logRadius, logFitHalo, "r-")
plt.xlabel("log(radius (m))")
plt.ylabel("log(halo density(kg/m^3))")
plt.legend(["Deduced halo density", str(fit[0])+"log(R) + "+str(fit[1])])
plt.title("Dark matter density from rotation curve")
plt.figure()
plt.plot(radiiForPlot, haloForPlot, "bo")
plt.plot(radiiForPlot, coreForPlot, "r-")
plt.xlabel("Radius (m)")
plt.ylabel("Density (kg/m^3)")
plt.legend(["Deduced from model"
            , str(rho0Core)+"/"+"(1 + (r/"+str(aCore)+")^2)"])
plt.title("Dark matter density from rotation curve")
plt.figure()
plt.plot(radius, velocityS, "r-")
plt.plot(radius, velocityG, "g-")
plt.plot(radius, velocityH, "k-")
plt.plot(radius, velocityB, "m-")
plt.plot(radius, rotationCurve, "c-")
plt.xlabel("Radius (m)")
plt.ylabel("Velocity (m/s)")
plt.legend(["Stellar potential only", "Gas potential only",
    "Halo potential only", "Baryon potential only", "All matter"])
plt.title("Rotation curve")
#plt.savefig("/Users/virinchirallabhandi/Desktop/Rotation curves 4")
plt.figure()
plt.plot(gBar, gObs, "bo")
plt.plot(mcGaughRange, mcGaughFit, "r-")
plt.plot(mcGaughRange, unity, "g-")
plt.xlabel("log(gBar (m/s^2))")
plt.ylabel("log(gObs (m/s^2))")
plt.legend(["Points from simulated galaxies", "McGaugh et al. fit", "gObs = gBar"])
plt.title("MDAR from one galaxy")