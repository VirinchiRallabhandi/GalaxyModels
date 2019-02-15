#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 14:46:35 2018
Good luck if you're trying to understand this module 
@author: Virinchi Rallabhandi
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from mpl_toolkits import mplot3d
"""Gives the circular velocity at a radius r as given by the URC with
   parameters v0, a and beta. rOpt is taken to be 3.2r0"""
def URC(r, lam):
    rOpt = 13*np.sqrt(lam)
    #rOpt is now in units of kpc
    Kpc = 3.086*10**19
    r = r/Kpc
    #r is now also in units of kpc
    v1 = np.sqrt(0.8 + 0.49*np.log10(lam) + 0.75*pow(np.e, -0.4*lam)/(
        0.47 + 2.25*pow(lam, 0.4)))
    v1 = 200*pow(lam, 0.41)/v1
    v = (0.72 + 0.44*np.log10(lam))*1.97*pow(r/rOpt, 1.22)/pow(
        pow(r/rOpt, 2) + 0.61, 1.43)
    v = v + 1.6*pow(np.e, -0.4*lam)*pow(r/rOpt, 2)/(pow(r/rOpt, 2)
        + 2.25*pow(lam, 0.4))
    return v1*np.sqrt(v)*1000
"""The derivative of the URC with respect to radius. The additional parameter
   is v = URC(r, v0, a, beta. It is provided to improve time efficiency."""
def dURCdr(r, lam):
    rOpt = 13*pow(lam, 0.5)
    #rOpt is now in units of kpc
    Kpc = 3.086*10**19
    r = r/Kpc
    #r is now also in units of kpc
    v1 = np.sqrt(0.8 + 0.49*np.log10(lam) + 0.75*pow(np.e, -0.4*lam)/(
        0.47 + 2.25*pow(lam, 0.4)))
    v1 = 200*pow(lam, 0.41)/v1
    x = r/rOpt
    a = 0.72 + 0.44*np.log10(lam)
    b = 1.6*pow(np.e, -0.4*lam)
    c = 2.25*pow(lam, 0.4)
    dv2dx = 2.4034*a*pow(x, 0.22)*pow(x*x + 0.61, 1.43)
    dv2dx = dv2dx - 5.6342*a*pow(x, 2.22)*pow(x*x + 0.61, 0.43)
    dv2dx = dv2dx/pow(x*x + 0.61, 2.86) + 2*b*x*c/pow(x*x + c, 2)
    return dv2dx*v1*v1/2/URC(r*Kpc, lam)/rOpt*1000/Kpc
    #I don't think this derivative method is correctly implemented yet but it
    #isn't used for anything useful in this module so I've neglected it
Q = 1.6
G = 6.67408*10**-11

#This section is used for comparing the results of our model and those of
#the baryonic Tully Fisher relation. mass.txt and vMax3.txt were files I made
#using the testBTFR method of Tester.java
sigmaG = 8000.0
logV = []
logMBar = []
massFile = open("/Users/virinchirallabhandi/Desktop/Data/mass3.txt")
masses = massFile.readlines()
massFile.close()
vMaxFile = open("/Users/virinchirallabhandi/Desktop/Data/vMax3.txt")
vMax = vMaxFile.readlines()
vMaxFile.close()
for i in range(len(masses)):
    logMBar.append(np.log10(float(masses[i])))
    logV.append(np.log10(float(vMax[i])))
fit = np.polyfit(logV, logMBar, 1)
#print(fit)
plt.figure()
plt.plot(logV, logMBar, "bo")
plt.close()

#This section visualises the effect of lambda and v0 in the URC of GalaxyModel4 
x = np.linspace(0.02, 5, 100)
y = []
for lam in x:
    v0 = 0.75*np.e**(-0.4*lam)/(0.47 + 2.25*lam**0.4) + 0.8 + 0.49*np.log10(lam)
    v0 = 200*lam**0.41*np.sqrt(v0)
    y.append(v0)
plt.figure()
plt.plot(x, y)
plt.close()

#This section finds the relation between r25 and r0 (called rd here) for the
#sample of Ianjamasimanana et al. regarding gas velocity dispersions.
rd = [0.9, 1.3, 1.2, 0.7, 2.4, 1.6, 1.3, 2.1, 2.3, 2.9, 2.5, 3.2,
      1.2, 3.2, 2.6, 2.5, 2.5, 4.1]
r25 = [3.8, 2.2, 3.3, 2.9, 12.0, 7.4, 5.9, 7.5, 10.4, 15.2, 10.6, 
       17.3, 5.3, 13.0, 9.4, 9.9, 10.1, 14.3]
linearRegression = np.polyfit(rd, r25, 1)
m = linearRegression[0]
c = linearRegression[1]
#print(linearRegression)
r = 0.0
#The correlation coefficient
rdError = 0.0
r25Error = 0.0
rdMean = 0.0
r25Mean = 0.0
for i in range(len(rd)):
    rdMean = rdMean + rd[i]
    r25Mean = r25Mean + r25[i]
rdMean = rdMean/len(rd)
r25Mean = r25Mean/len(r25)
for i in range(len(rd)):
    rdError = rdError + (rd[i] - rdMean)**2
    r25Error = r25Error + (r25[i] - r25Mean)**2
rdError = np.sqrt(rdError)
r25Error = np.sqrt(r25Error)
for i in range(len(rd)):
    r = r + (rd[i] - rdMean)*(r25[i] - r25Mean)
r = r/rdError/r25Error
#print(rdError)
#print(r)
sigma = 0.0
for i in range(len(rd)):
    sigma = sigma + (r25[i] - m*rd[i] - c)**2
sigma = np.sqrt(sigma/len(rd))
#print(sigma)
fitPoints = []
for value in rd:
    fitPoints.append(value*m + c)
plt.figure()
plt.plot(rd, r25, "bo")
plt.plot(rd, fitPoints, "r-")
plt.xlabel("r0 (kpc)")
plt.ylabel("r25 (kpc)")
plt.title("Relationship between r0 and r25")
plt.close()

#This section creates figure 21 in Studentship report.pdf but plots different
#Q value data in different figures. The least square error for each Q value
#is printed/ gasMassToLight was made using the getGasMasses method of Tester.java.
vMaxObs = []
massToLightObs = []
dataFile = open("/Users/virinchirallabhandi/Desktop/Data/Observed gas data.txt", "r")
data = dataFile.readlines()
dataFile.close()
for i in range(len(data)):
    numbers = data[i].split()
    vMaxObs.append(float(numbers[1]))
    massToLightObs.append(float(numbers[2]) + np.log10(1.3))
dataFile = open("/Users/virinchirallabhandi/Desktop/Data/gasMassToLight.txt", "r")
data = dataFile.readlines()
dataFile.close()
gasMassToLight = []
vMax = []
vMaxLower = []
vMaxHigher = []
gasMassLower = []
gasMassHigher = []
qError = []
numbers = data[0].split()
qValues = -1
Q = float(numbers[2])
for i in range(len(data)):
    numbers = data[i].split()
    if float(numbers[2]) == Q:
        gasMassToLight.append(float(numbers[0]))
        vMax.append(float(numbers[1]))
        if float(numbers[1]) > 2.32:
            vMaxHigher.append(float(numbers[1]))
            gasMassHigher.append(float(numbers[0]))
        else:
            vMaxLower.append(float(numbers[1]))
            gasMassLower.append(float(numbers[0]))
    else:
        qError.append(0)
        qValues = qValues + 1
        lowerFit = np.polyfit(vMaxLower, gasMassLower, 1)
        higherFit = np.polyfit(vMaxHigher, gasMassHigher, 1)
        piecewiseFit = []
        for j in range(len(vMaxObs)):
            if vMaxObs[j] > 2.32:
                qError.append(0)
                del qError[len(qError) - 1]
                fitValue = higherFit[0]*vMaxObs[j] + higherFit[1]
                qError[qValues] = qError[qValues] + (fitValue - massToLightObs[j])**2
                piecewiseFit.append(fitValue)
            else:
                fitValue = lowerFit[0]*vMaxObs[j] + lowerFit[1]
                qError[qValues] = qError[qValues] + (fitValue - massToLightObs[j])**2
                piecewiseFit.append(fitValue)
        plt.figure()
        plt.plot(vMax, gasMassToLight, "bo")
        plt.plot(vMaxObs, massToLightObs, "ko")
        plt.plot(vMaxObs, piecewiseFit, "r-")
        plt.xlabel("log(vMax (km/s))")
        plt.ylabel("log(gas mass/R band luminosity (solar units))")
        plt.legend(["Data from model", "Observed data", "Best fit to model"])
        Q = float(int(Q*10))/10
        plt.title("Gas mass to light as a function of vMax for Q = "+str(Q))
        plt.close()
        for j in range(len(gasMassToLight)):
            del gasMassToLight[len(gasMassToLight) - 1]
            del vMax[len(vMax) - 1]
        for j in range(len(vMaxLower)):
            del vMaxLower[len(vMaxLower) - 1]
            del gasMassLower[len(gasMassLower) - 1]
        for j in range(len(vMaxHigher)):
            del vMaxHigher[len(vMaxHigher) - 1]
            del gasMassHigher[len(gasMassHigher) - 1]
        Q = float(numbers[2])
        gasMassToLight.append(float(numbers[0]))
        vMax.append(float(numbers[1]))
qError.append(0)
qValues = qValues + 1
gasMassToLight.append(float(numbers[0]))
vMax.append(float(numbers[1]))
if float(numbers[1]) > 2.32:
    vMaxHigher.append(float(numbers[1]))
    gasMassHigher.append(float(numbers[0]))
else:
    vMaxLower.append(float(numbers[1]))
    gasMassLower.append(float(numbers[0]))
lowerFit = np.polyfit(vMaxLower, gasMassLower, 1)
higherFit = np.polyfit(vMaxHigher, gasMassHigher, 1)
for j in range(len(vMaxObs)):
    if vMaxObs[j] > 2.32:
        qError[qValues] = qError[qValues] + (higherFit[0]*vMaxObs[j] 
            + higherFit[1] - massToLightObs[j])**2
    else:
        qError[qValues] = qError[qValues] + (lowerFit[0]*vMaxObs[j] 
            + lowerFit[1] - massToLightObs[j])**2
plt.figure()
plt.plot(vMax, gasMassToLight, "bo")
plt.plot(vMaxObs, massToLightObs, "ko")
plt.xlabel("vMax (km/s)")
plt.ylabel("Gas mass (solar masses)")
Q = float(int(Q*10))/10
plt.title("Gas mass as a function of vMax for Q = "+str(Q))
plt.close()
#print(qError)

#This section creates figure 21 in Studentship report.pdf
vMaxObs = []
massToLightObs = []
dataFile = open("/Users/virinchirallabhandi/Desktop/Data/Observed gas data.txt", "r")
data = dataFile.readlines()
dataFile.close()
for i in range(len(data)):
    numbers = data[i].split()
    vMaxObs.append(float(numbers[1]))
    massToLightObs.append(float(numbers[2]) + np.log10(1.3))
dataFile = open("/Users/virinchirallabhandi/Desktop/Data/gasMassToLight.txt", "r")
data = dataFile.readlines()
dataFile.close()
vMax1 = []
vMax2 = []
vMax4 = []
gasML1 = []
gasML2 = []
gasML4 = []
for i in range(len(data)):
    numbers = data[i].split()
    q = float(numbers[2])
    if q > 0.99 and q < 1.01:
        vMax1.append(float(numbers[1]))
        gasML1.append(float(numbers[0]))
    if q > 1.99 and q < 2.01:
        vMax2.append(float(numbers[1]))
        gasML2.append(float(numbers[0]))
    if q > 3.99 and q < 4.01:
        vMax4.append(float(numbers[1]))
        gasML4.append(float(numbers[0]))
plt.figure()
plt.plot(vMaxObs, massToLightObs, "ko")
plt.plot(vMax1, gasML1, "co")
plt.plot(vMax2, gasML2, "go")
plt.plot(vMax4, gasML4, "mo")
plt.xlabel("$V_{max}$ (km/s)")
plt.ylabel("$M_{gas}/L_R$ (solar units)")
plt.title("Gas mass to light ratio as a function of $V_{max}$ and $Q$")
plt.legend(["Observed points", "$Q = 1$ models", "$Q = 2$ models", "$Q = 4$ models"])
#plt.savefig("/Users/virinchirallabhandi/Desktop/Data/Mass to light for varying Q")
plt.close()

#This section fits my own model to the 14 normalised data points from which
#Karukes and Salucci's dwarf disk URC is made.
normalisedRadii = [0.11, 0.22, 0.32, 0.41, 0.52, 0.63, 0.77, 0.91, 1.03, 
    1.18, 1.32, 1.45, 1.65, 1.88]
normalisedVelocity = [0.21, 0.37, 0.49, 0.6, 0.68, 0.78, 0.86, 0.95, 0.99,
    1.05, 1.07, 1.07, 1.12, 1.2]
def exponentialProfile(r, v0, a):
    return v0*(1 - np.e**(-r/a))
def arctanProfile(r, v0, rT):
    return 2*v0/np.pi*np.arctan(r/rT)
expFit, covariance = opt.curve_fit(exponentialProfile, normalisedRadii, normalisedVelocity)
arctanFit, covariance = opt.curve_fit(arctanProfile, normalisedRadii, normalisedVelocity)
expv0Fit = expFit[0]
aFit = expFit[1]
arctanv0Fit = arctanFit[0]
rTFit = arctanFit[1]
#print(arctanv0Fit, rTFit)
expFitCurve = []
arctanFitCurve = []
for i in range(len(normalisedRadii)):
    expFitCurve.append(exponentialProfile(normalisedRadii[i], expv0Fit, aFit))
    arctanFitCurve.append(arctanProfile(normalisedRadii[i], arctanv0Fit, rTFit))
plt.figure()
plt.plot(normalisedRadii, normalisedVelocity, "bo")
plt.plot(normalisedRadii, expFitCurve, "r-")
plt.plot(normalisedRadii, arctanFitCurve, "g-")
plt.close()
#print(expv0Fit, aFit)
#print(arctanv0Fit, rTFit)

#This section plots the vMax - r0 relation for Wong et al (2016)
vMaxValues = []
r0Values = []
for i in range(50):
    vMax = float(i)*5000.0 + 20000.0
    vMaxValues.append(vMax/1000)
    mR = -3.9 - 7.622*np.log10(vMax/1000)
    #The absolute magnitude of the galaxy in the R band
    surfaceBrightness = pow(10, 5.3785 + 1.1757*np.log10(vMax/1000))
    #The effective surface brightness of the disk in Lsun/kpc^2
    r0 = np.sqrt(pow(10, -0.4*(mR - 4.61))/5.647/np.pi/surfaceBrightness)
    #This is in units of Kpc
    r0Values.append(r0)
r0Fit = []
linearFit = np.polyfit(vMaxValues, r0Values, 1)
for i in range(len(vMaxValues)):
    r0Fit.append(vMaxValues[i]*linearFit[0] + linearFit[1])
plt.figure()
plt.plot(vMaxValues, r0Values, "bo")
plt.plot(vMaxValues, r0Fit, "r-")
plt.close()

#This section attempts (unsuccessfully) to reformulate the URC in terms of the arctan
#model of Courteau (1997)
vMaxValues = []
rTValues = []
#for i in range(11):
#    vMax = float(i)*5000.0 + 20000.0
#    vMaxValues.append(vMax/1000)
#    mR = -3.9 - 7.622*np.log10(vMax/1000)
    #The absolute magnitude of the galaxy in the R band
#    surfaceBrightness = pow(10, 5.3785 + 1.1757*np.log10(vMax/1000))
    #The effective surface brightness of the disk in Lsun/kpc^2
#    r0 = np.sqrt(pow(10, -0.4*(mR - 4.61))/5.647/np.pi/surfaceBrightness)
    #This is in units of Kpc
#    print(r0)
#    rTValues.append(rTFit*3.2*r0)
for count in range(37):
    radius = []
    urcCurve = []
    vFinal = float(count)*5000.0 + 75000.0
    rMax = vFinal*pow(10, 9)*365.25*24*3600/2/np.pi
    mB = -10*np.log10(vFinal/1000) + 3;
    #Very rough estimate from expressions in Carrol and Ostlie and assuming
    #vMax = vFinal. These assumptions are only made so that a reasonable
    #initial guess is made for the parameter lambda.
    luminosity = pow(10, -0.4*(mB - 5.31));
    lam = luminosity/6/pow(10, 10)*pow(67.4/50.0, 2)
    #Using H0 = 67 km/s/Mpc, the latest value by the Planck collaboration
    #Iteratively work out what the correct value of lambda is in the URC
    for i in range(2000):
        vFinalPredicted = URC(rMax, lam);
        lam = lam + (vFinal - vFinalPredicted)/2000000
        #v0 in the URC is a strictly increasing function of lambda so
        #this process should work given enough iterations
    for i in range(70):
        radius.append(float(i)*rMax/100)
        urcCurve.append(URC(radius[i], lam)/10000)
        radius[i] = radius[i]/3.086/10**19
        #The units must be adjusted to better suit the capabilities of the
        #fitting algorithm used by scipy
    arctanFit, covariance = opt.curve_fit(arctanProfile, radius, urcCurve)
    vMaxValues.append(arctanFit[0]*10)
    rTValues.append(arctanFit[1])
cubicFit = np.polyfit(vMaxValues, rTValues, 3)
linearFit = np.polyfit(vMaxValues, rTValues, 1)
#print(linearFit)
#print(cubicFit)
cubicCurve = []
linearCurve = []
vMaxRange = []
for i in range(31):
    v = float(i)*10.0 + 20.0
    vMaxRange.append(v)
    cubicCurve.append(cubicFit[0]*v**3 + cubicFit[1]*v**2 + cubicFit[2]*v + cubicFit[3])
    linearCurve.append(linearFit[0]*v + linearFit[1])
#print(vMaxValues)
plt.figure()
plt.plot(vMaxValues, rTValues, "bo")
plt.plot(vMaxRange, cubicCurve, "r-")
plt.plot(vMaxRange, linearCurve, "g-")
plt.xlabel("vMax (km/s)")
plt.ylabel("rT (kpc)")
plt.title("Scatter plot of vMax and rT value for arctan reformulation of URC")
plt.close()

#Again an arctan fit to the URC
radius = []
urcCurve = []
vFinal = 70000.0
rMax = vFinal*pow(10, 9)*365.25*24*3600/2/np.pi;
mB = -10*np.log10(vFinal/1000) + 3;
#Very rough estimate from expressions in Carrol and Ostlie and assuming
#vMax = vFinal. These assumptions are only made so that a reasonable
#initial guess is made for the parameter lambda.
luminosity = pow(10, -0.4*(mB - 5.31));
lam = luminosity/6/pow(10, 10)*pow(67.4/50.0, 2);
#Using H0 = 67 km/s/Mpc, the latest value by the Planck collaboration
#Iteratively work out what the correct value of lambda is in the URC
for i in range(2000):
	vFinalPredicted = URC(rMax, lam);
	lam = lam + (vFinal - vFinalPredicted)/2000000;
	#v0 in the URC is a strictly increasing function of lambda so
	#this process should work given enough iterations
for i in range(100):
    radius.append(float(i)*rMax/100)
    urcCurve.append(URC(radius[i], lam)/10000)
    radius[i] = radius[i]/3.086/10**19
    #The units must be adjusted to better suit the capabilities of the
    #fitting algorithm used by scipy
arctanFit, covariance = opt.curve_fit(arctanProfile, radius, urcCurve)
arctanv0Fit = arctanFit[0]*10
rTFit = arctanFit[1]
#print(arctanv0Fit, rTFit)
arctanFitCurve = []
for i in range(len(radius)):
    arctanFitCurve.append(arctanProfile(radius[i], arctanv0Fit, rTFit))
    urcCurve[i] = urcCurve[i]*10
plt.figure()
plt.plot(radius, urcCurve, "r-")
plt.plot(radius, arctanFitCurve, "b-")
plt.xlabel("Radius (kpc)")
plt.ylabel("Velocity (km/s)")
plt.close()

#Creates figure 19 in studentship.pdf
dwarfURC = []
mR = -3.9 - 7.622*np.log10(vFinal/1000)
#The absolute magnitude of the galaxy in the R band
surfaceBrightness = pow(10, 5.3785 + 1.1757*np.log10(vFinal/1000))
#The effective surface brightness of the disk in Lsun/kpc^2
r0 = np.sqrt(pow(10, -0.4*(mR - 4.61))/5.647/np.pi/surfaceBrightness)
for i in range(len(radius)):
    dwarfURC.append(vFinal/1000*(1 - np.e**(-radius[i]/1.995/r0)))
plt.figure()
plt.plot(radius, urcCurve, "r-")
plt.plot(radius, dwarfURC, "b-")
plt.xlabel("Radius (kpc)")
plt.ylabel("Velocity (km/s)")
plt.title("URC discontinuity")
plt.legend(["Spiral galaxy URC", "Dwarf disk URC"])
plt.close()

#Plots star formation efficiency predicted from my models against vMax
#to compare with the results of Wong et al. (2016)
dataFile = open("/Users/virinchirallabhandi/Desktop/Data/SFE5.txt", "r")
data = dataFile.readlines()
dataFile.close()
logSFE = []
logVMax =[]
for i in range(len(data)):
    numbers = data[i].split()
    logVMax.append(float(numbers[0]))
    logSFE.append(float(numbers[1]))
plt.figure()
plt.plot(logVMax, logSFE, "bo")
plt.xlabel("$\log(V_{max} \, (km/s))$")
plt.ylabel("$\log(SFE (1/yr)})$")
plt.title("Star formation efficieny (HI) as a function of $V_{max}$")
#plt.ylim(top = -8)
#plt.ylim(bottom = -13)
#plt.close()

#A scatter plot looking for any correlations between the parameters of Courteau's
#arctan and multi parameter fit models of rotation curves.
dataFile = open("/Users/virinchirallabhandi/Desktop/Data/Arctan fits.txt", "r")
data = dataFile.readlines()
dataFile.close()
vC = []
rT = []
for i in range(len(data)):
    numbers = data[i].split()
    if np.abs(float(numbers[0])) < 800:
        vC.append(np.abs(float(numbers[0])))
        rT.append(float(numbers[1]))
plt.figure()
plt.plot(vC, rT, "bo")
plt.xlabel("vC (km/s)")
plt.ylabel("rT (kpc)")
plt.title("A scatter plot of vC and rT parameters in Courteau's arctan model")
plt.close()
dataFile = open("/Users/virinchirallabhandi/Desktop/Data/Multi parameter fits.txt", "r")
data = dataFile.readlines()
dataFile.close()
vC = []
rT = []
gamma = []
beta = []
for i in range(len(data)):
    numbers = data[i].split()
    if (np.abs(float(numbers[0])) < 400 and float(numbers[1]) < 50 
        and float(numbers[2]) < 20):
        vC.append(np.abs(float(numbers[0])))
        rT.append(float(numbers[1]))
        gamma.append(float(numbers[2]))
        beta.append(float(numbers[3]))
plt.figure()
ax = plt.axes(projection = "3d")
ax.set_xlabel("gamma")
ax.set_ylabel("rT (kpc)")
ax.set_zlabel("vC (km/s)")
ax.scatter3D(gamma, rT, vC, c = rT)
plt.title("Scatter plot of multi paramter fit values in Courteau's model")
plt.close()