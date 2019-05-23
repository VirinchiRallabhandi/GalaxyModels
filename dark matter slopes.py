# -*- coding: utf-8 -*-
"""
Created on Thu May 23 19:19:19 2019

@author: Virinchi Rallabhandi
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp1
import scipy.interpolate as sp
G = 6.67408*10**-11
gD = 1.2*10**-10
def sech(x):
    return 1.0/np.cosh(x)
def f1(x):
    return 1.0/(1.0 - np.e**(-np.sqrt(x)))
def df1(x):
    return -np.e**(-np.sqrt(x))/2.0/np.sqrt(x)/(1.0 - np.e**(-np.sqrt(x)))**2
def d2f1(x):
    return (np.e**(-x**0.5)*(1.0 - np.e**(-x**0.5))*(
        1.0 + 1.0/x**0.5) + 2.0*np.e**(-2.0*x**0.5))/4.0/x/(1.0 - np.e**(-x**0.5))**3
def f2(x):
    return (1 + np.sqrt(1.0 + 4.0/x))/2.0
def df2(x):
    return -1.0/x/x/np.sqrt(1.0 + 4.0/x)
def d2f2(x):
    return 2.0*(x + 3.0)/x**4/(1.0 + 4.0/x)**1.5
def f3(x):
    return np.tanh(x**(1.0/3.5))**-1.75
def df3(x):
    return -0.5*x**(1.0/3.5 - 1.0)*np.tanh(x**(1.0/3.5))**-2.75*sech(x**(1.0/3.5))**2
def d2f3(x):
    return (-0.5*(1.0/3.5 - 1.0)*x**(1.0/3.5 - 2.0)*np.tanh(x**(1.0/3.5))**-2.75*sech(
        x**(1.0/3.5))**2 + 2.75/7.0*x**(1.0/1.75 - 2.0)*np.tanh(x**(1.0/3.5))**-3.75*sech(
        x**(1.0/3.5))**4 - 1.0/3.5*x**(1.0/1.75 - 2.0)*np.tanh(x**(1.0/3.5))**-1.75*sech(
        x**(1.0/3.5))**2)

#Graphing the exponential disk case with f1
r0Step = 20.0**0.25
mStarStep = 1000.0**0.25
fig, ax = plt.subplots(2, 3, sharex = True, sharey = True, figsize = (11, 14))
fig.subplots_adjust(hspace = 0.2, wspace = 0.8)
for rCount in range(5):
    r0 = 1.0*10**19*r0Step**rCount
    for densityCount in range(5):
        maxDensity = 1.0*10**39*mStarStep**densityCount/r0/r0/2.0/np.pi
        radius = np.zeros(1000)
        logHaloSlope = np.zeros(len(radius))
        for i in range(len(radius)):
            r = float(i + 1)*r0/50
            radius[i] = r
            s = 0.5*r/r0
            i0 = sp1.iv(0, s)
            i1 = sp1.iv(1, s)
            i2 = sp1.iv(2, s)
            i3 = sp1.iv(3, s)
            k0 = sp1.kv(0, s)
            k1 = sp1.kv(1, s)
            k2 = sp1.kv(2, s)
            k3 = sp1.kv(3, s)
            gBar = np.pi*G*maxDensity*r/r0*(i0*k0 - i1*k1)
            x = float(gBar/gD)
            dgBardr = gBar/r + np.pi*G*maxDensity*r/2/r0/r0*(i1*k0 - i0*k1 - 0.5*k1*(
                i0 + i2) + 0.5*i1*(k0 + k2))
            d2gBardr2 = -gBar/r/r + dgBardr/r + np.pi*G*maxDensity/2/r0/r0*(i1*k0 - i0*k1
                - 0.5*k1*(i0 + i2) + 0.5*i1*(k0 + k2)) + np.pi*G*maxDensity*r/4/r0/r0/r0*(
                0.5*(i0 + i2)*k0 - 2.0*i1*k1 + 0.5*i0*(k0 + k2) + 0.5*(k2 + k0)*(i0 
                + i2) - 0.5*k1*(i1 + 0.5*(i1 + i3)) + 0.5*i1*(-k1 -0.5*(k1 + k3)))
            rhoDM = gBar*(f1(x) - 1.0)/2.0/np.pi/G/r + 1.0/4.0/np.pi/G*(
                f1(x) - 1.0 + x*df1(x))*dgBardr
            drhoDMdr = 1.0/2.0/np.pi/G/r/r*(r*(f1(x) - 1.0 + x*df1(x))*dgBardr - gBar*(
                f1(x) - 1.0)) + 1.0/4.0/np.pi/G/gD*(2.0*df1(x) + x*d2f1(x)
                )*dgBardr**2 + 1.0/4.0/np.pi/G*(f1(x) - 1.0 + x*df1(x))*d2gBardr2
            logHaloSlope[i] = r/rhoDM*drhoDMdr
        #print(logHaloSlope[0])
        #print(logHaloSlope[len(logHaloSlope) - 1])
        for i in range(len(radius)):
            radius[i] = radius[i]/r0
        plt.subplot(5, 5, 5*rCount + densityCount + 1)
        plt.plot(radius, logHaloSlope, "c-")
        #plt.title(str(rCount)+"  "+str(densityCount))
fig.add_subplot(111, frameon=False)
plt.grid(False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel("Radius (scale lengths)")
plt.ylabel("logarithmic derivative of dark matter density")
plt.title("Slope of log(dark matter density) as per MDAR")
plt.tight_layout()
#plt.savefig("D:/University/Physics/PHYS3003/exponential1.pdf")

#Graphing the exponential disk case with f2
r0Step = 20.0**0.25
mStarStep = 1000.0**0.25
fig, ax = plt.subplots(2, 3, sharex = True, sharey = True, figsize = (11, 14))
fig.subplots_adjust(hspace = 0.2, wspace = 0.8)
for rCount in range(5):
    r0 = 1.0*10**19*r0Step**rCount
    for densityCount in range(5):
        maxDensity = 1.0*10**39*mStarStep**densityCount/r0/r0/2.0/np.pi
        radius = np.zeros(1000)
        logHaloSlope = np.zeros(len(radius))
        for i in range(len(radius)):
            r = float(i + 1)*r0/50
            radius[i] = r
            s = 0.5*r/r0
            i0 = sp1.iv(0, s)
            i1 = sp1.iv(1, s)
            i2 = sp1.iv(2, s)
            i3 = sp1.iv(3, s)
            k0 = sp1.kv(0, s)
            k1 = sp1.kv(1, s)
            k2 = sp1.kv(2, s)
            k3 = sp1.kv(3, s)
            gBar = np.pi*G*maxDensity*r/r0*(i0*k0 - i1*k1)
            x = float(gBar/gD)
            dgBardr = gBar/r + np.pi*G*maxDensity*r/2/r0/r0*(i1*k0 - i0*k1 - 0.5*k1*(
                i0 + i2) + 0.5*i1*(k0 + k2))
            d2gBardr2 = -gBar/r/r + dgBardr/r + np.pi*G*maxDensity/2/r0/r0*(i1*k0 - i0*k1
                - 0.5*k1*(i0 + i2) + 0.5*i1*(k0 + k2)) + np.pi*G*maxDensity*r/4/r0/r0/r0*(
                0.5*(i0 + i2)*k0 - 2.0*i1*k1 + 0.5*i0*(k0 + k2) + 0.5*(k2 + k0)*(i0 
                + i2) - 0.5*k1*(i1 + 0.5*(i1 + i3)) + 0.5*i1*(-k1 -0.5*(k1 + k3)))
            rhoDM = gBar*(f2(x) - 1.0)/2.0/np.pi/G/r + 1.0/4.0/np.pi/G*(
                f2(x) - 1.0 + x*df2(x))*dgBardr
            drhoDMdr = 1.0/2.0/np.pi/G/r/r*(r*(f2(x) - 1.0 + x*df2(x))*dgBardr - gBar*(
                f2(x) - 1.0)) + 1.0/4.0/np.pi/G/gD*(2.0*df2(x) + x*d2f2(x)
                )*dgBardr**2 + 1.0/4.0/np.pi/G*(f2(x) - 1.0 + x*df2(x))*d2gBardr2
            logHaloSlope[i] = r/rhoDM*drhoDMdr
        #print(logHaloSlope[0])
        #print(logHaloSlope[len(logHaloSlope) - 1])
        for i in range(len(radius)):
            radius[i] = radius[i]/r0
        plt.subplot(5, 5, 5*rCount + densityCount + 1)
        plt.plot(radius, logHaloSlope, "c-")
        #plt.title(str(rCount)+"  "+str(densityCount))
fig.add_subplot(111, frameon=False)
plt.grid(False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel("Radius (scale lengths)")
plt.ylabel("logarithmic derivative of dark matter density")
plt.title("Slope of log(dark matter density) as per MDAR")
plt.tight_layout()
#plt.savefig("D:/University/Physics/PHYS3003/exponential2.pdf")

#Graphing the exponential disk case with f3
r0Step = 20.0**0.25
mStarStep = 1000.0**0.25
fig, ax = plt.subplots(2, 3, sharex = True, sharey = True, figsize = (11, 14))
fig.subplots_adjust(hspace = 0.2, wspace = 0.8)
for rCount in range(5):
    r0 = 1.0*10**19*r0Step**rCount
    for densityCount in range(5):
        maxDensity = 1.0*10**39*mStarStep**densityCount/r0/r0/2.0/np.pi
        radius = np.zeros(1000)
        logHaloSlope = np.zeros(len(radius))
        for i in range(len(radius)):
            r = float(i + 1)*r0/50
            radius[i] = r
            s = 0.5*r/r0
            i0 = sp1.iv(0, s)
            i1 = sp1.iv(1, s)
            i2 = sp1.iv(2, s)
            i3 = sp1.iv(3, s)
            k0 = sp1.kv(0, s)
            k1 = sp1.kv(1, s)
            k2 = sp1.kv(2, s)
            k3 = sp1.kv(3, s)
            gBar = np.pi*G*maxDensity*r/r0*(i0*k0 - i1*k1)
            x = float(gBar/gD)
            dgBardr = gBar/r + np.pi*G*maxDensity*r/2/r0/r0*(i1*k0 - i0*k1 - 0.5*k1*(
                i0 + i2) + 0.5*i1*(k0 + k2))
            d2gBardr2 = -gBar/r/r + dgBardr/r + np.pi*G*maxDensity/2/r0/r0*(i1*k0 - i0*k1
                - 0.5*k1*(i0 + i2) + 0.5*i1*(k0 + k2)) + np.pi*G*maxDensity*r/4/r0/r0/r0*(
                0.5*(i0 + i2)*k0 - 2.0*i1*k1 + 0.5*i0*(k0 + k2) + 0.5*(k2 + k0)*(i0 
                + i2) - 0.5*k1*(i1 + 0.5*(i1 + i3)) + 0.5*i1*(-k1 -0.5*(k1 + k3)))
            rhoDM = gBar*(f3(x) - 1.0)/2.0/np.pi/G/r + 1.0/4.0/np.pi/G*(
                f3(x) - 1.0 + x*df3(x))*dgBardr
            drhoDMdr = 1.0/2.0/np.pi/G/r/r*(r*(f3(x) - 1.0 + x*df3(x))*dgBardr - gBar*(
                f3(x) - 1.0)) + 1.0/4.0/np.pi/G/gD*(2.0*df3(x) + x*d2f3(x)
                )*dgBardr**2 + 1.0/4.0/np.pi/G*(f3(x) - 1.0 + x*df3(x))*d2gBardr2
            logHaloSlope[i] = r/rhoDM*drhoDMdr
        #print(logHaloSlope[0])
        #print(logHaloSlope[len(logHaloSlope) - 1])
        for i in range(len(radius)):
            radius[i] = radius[i]/r0
        plt.subplot(5, 5, 5*rCount + densityCount + 1)
        plt.plot(radius, logHaloSlope, "c-")
        #plt.title(str(rCount)+"  "+str(densityCount))
fig.add_subplot(111, frameon=False)
plt.grid(False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel("Radius (scale lengths)")
plt.ylabel("logarithmic derivative of dark matter density")
plt.title("Slope of log(dark matter density) as per MDAR")
plt.tight_layout()
#plt.savefig("D:/University/Physics/PHYS3003/exponential3.pdf")

#Graphing the data generated by GalaxyModel14.java for the stable disk model
fig, ax = plt.subplots(2, 3, sharex = True, sharey = True, figsize = (11, 14))
fig.subplots_adjust(hspace = 0.2, wspace = 0.8)
for vCount in range(6):
    for fCount in range(3):
        dataFile = open("D:/University/Physics/PHYS3003/Assignment/radius1"+str(
            vCount)+str(fCount)+".txt")
        data = dataFile.readlines()
        dataFile.close()
        radius = []
        for line in data:
            numbers = line.split()
            radius.append(float(numbers[0]))
        dataFile = open("D:/University/Physics/PHYS3003/Assignment/gBar1"+str(
            vCount)+str(fCount)+".txt")
        data = dataFile.readlines()
        dataFile.close()
        gBar = []
        for line in data:
            numbers = line.split()
            gBar.append(float(numbers[0]))
        logHaloSlope = np.zeros(int(float(len(radius))*0.75))
        radiusForFit = []
        gBarForFit = []
        for i in range(len(radius)/35):
            radiusForFit.append(radius[i*35])
            gBarForFit.append(gBar[i*35])
        tck = sp.splrep(radiusForFit, gBarForFit, k=5)
        r0 = (radius[1] - radius[0])*50.0
        for i in range(len(logHaloSlope)):
            if i == 0:
                i = i + 1
            r = radius[i]
            g = gBar[i]
            dg = sp.splev(r, tck, der=1)
            d2g = sp.splev(r, tck, der=2)
            x = float(g/gD)
            rhoDM = 0.0
            drhoDMdr = 0.0
            if fCount == 0:
                rhoDM = g*(f1(x) - 1.0)/2.0/np.pi/G/r + 1.0/4.0/np.pi/G*(
                    f1(x) - 1.0 + x*df1(x))*dg
                drhoDMdr = 1.0/2.0/np.pi/G/r/r*(r*(f1(x) - 1.0 + x*df1(x))*dg - g*(
                    f1(x) - 1.0)) + 1.0/4.0/np.pi/G/gD*(2.0*df1(x) + x*d2f1(x)
                    )*dg**2 + 1.0/4.0/np.pi/G*(f1(x) - 1.0 + x*df1(x))*d2g
            if fCount == 1:
                rhoDM = g*(f2(x) - 1.0)/2.0/np.pi/G/r + 1.0/4.0/np.pi/G*(
                    f2(x) - 1.0 + x*df2(x))*dg
                drhoDMdr = 1.0/2.0/np.pi/G/r/r*(r*(f2(x) - 1.0 + x*df2(x))*dg - g*(
                    f2(x) - 1.0)) + 1.0/4.0/np.pi/G/gD*(2.0*df2(x) + x*d2f2(x)
                    )*dg**2 + 1.0/4.0/np.pi/G*(f2(x) - 1.0 + x*df2(x))*d2g
            if fCount == 2:
                rhoDM = g*(f3(x) - 1.0)/2.0/np.pi/G/r + 1.0/4.0/np.pi/G*(
                    f3(x) - 1.0 + x*df3(x))*dg
                drhoDMdr = 1.0/2.0/np.pi/G/r/r*(r*(f3(x) - 1.0 + x*df3(x))*dg - g*(
                    f3(x) - 1.0)) + 1.0/4.0/np.pi/G/gD*(2.0*df3(x) + x*d2f3(x)
                    )*dg**2 + 1.0/4.0/np.pi/G*(f3(x) - 1.0 + x*df3(x))*d2g
            logHaloSlope[i] = r/rhoDM*drhoDMdr
        logHaloSlope[0] = logHaloSlope[1]
        #print(logHaloSlope[0])
        #print(logHaloSlope[len(logHaloSlope) - 1])
        radiusForPlot = []
        for i in range(len(logHaloSlope)):
            radiusForPlot.append(radius[i]/r0)
        plt.subplot(6, 3, 3*vCount + fCount + 1)
        plt.plot(radiusForPlot, logHaloSlope, "c-")
        #plt.title(str(rCount)+"  "+str(densityCount))
fig.add_subplot(111, frameon=False)
plt.grid(False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel("Radius (scale lengths)")
plt.ylabel("logarithmic derivative of dark matter density")
plt.title("Slope of log(dark matter density) as per MDAR")
plt.tight_layout()
#plt.savefig("D:/University/Physics/PHYS3003/stableDisks2.pdf")