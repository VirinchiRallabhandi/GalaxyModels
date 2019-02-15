# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 20:14:16 2019
A module for plotting the implied dark matter density power law exponent so
that MDAR holds at every point.
@author: Virinchi Rallabhandi
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
G = 6.67408*10**-11
gD = 1.2*10**-10
fig, ax = plt.subplots(2, 3, sharex = True, sharey = True, figsize = (11, 14))
fig.subplots_adjust(hspace = 0.2, wspace = 0.8)
for rCount in range(5):
    r0 = 4.0*10**19 + 3.0*10**19*float(rCount)
    for densityCount in range(5):
        maxDensity = (2.5*10**39 + 1.0*10**41)/r0/r0/2/np.pi
        radius = np.zeros(2500)
        logHaloSlope = np.zeros(len(radius))
        for i in range(len(radius)):
            r = float(i + 1)*r0/50
            radius[i] = r
            x = 0.5*r/r0
            i0 = sp.iv(0, x)
            i1 = sp.iv(1, x)
            i2 = sp.iv(2, x)
            i3 = sp.iv(3, x)
            k0 = sp.kv(0, x)
            k1 = sp.kv(1, x)
            k2 = sp.kv(2, x)
            k3 = sp.kv(3, x)
            gBar = np.pi*G*maxDensity*r/r0*(i0*k0 - i1*k1)
            eg = np.e**(np.sqrt(gBar/gD))
            dgBardr = gBar/r + np.pi*G*maxDensity*r/2/r0/r0*(i1*k0 - i0*k1 - 0.5*k1*(
                i0 + i2) + 0.5*i1*(k0 + k2))
            d2gBardr2 = -gBar/r/r + dgBardr/r + np.pi*G*maxDensity/2/r0/r0*(i1*k0 - i0*k1
                - 0.5*k1*(i0 + i2) + 0.5*i1*(k0 + k2)) + np.pi*G*maxDensity*r/4/r0/r0/r0*(
                0.5*(i0 + i2)*k0 - 2.0*i1*k1 + 0.5*i0*(k0 + k2) + 0.5*(k2 + k0)*(i0 
                + i2) - 0.5*k1*(i1 + 0.5*(i1 + i3)) + 0.5*i1*(-k1 -0.5*(k1 + k3)))
            haloDensity = gBar/2/np.pi/G/r/(eg - 1) + (eg*(1 - 0.5*np.sqrt(gBar/gD
                )) - 1)*0.25/np.pi/G/(eg - 1)**2*dgBardr
            logHaloSlope[i] = r/haloDensity*(-gBar/2/np.pi/G/r/r/(eg - 1) + (eg*(1 - 0.5*np.sqrt(
                gBar/gD)) - 1)/2/np.pi/G/r/(eg - 1)**2*dgBardr + dgBardr**2*((eg - 1)**2*(
                0.5/np.sqrt(gBar*gD)*eg*(1 - 0.5*np.sqrt(gBar/gD)) - 0.25*eg/np.sqrt(
                gBar*gD)) - 2*(eg - 1)*0.5/np.sqrt(gBar*gD)*eg*(eg*(1 - 0.5*np.sqrt(gBar/gD))
                - 1))*0.25/np.pi/G/(eg - 1)**4 + d2gBardr2*(eg*(1 - 0.5*np.sqrt(gBar/gD)) - 1)
                *0.25/np.pi/G/(eg - 1)**2)
        print(logHaloSlope[0])
        print(logHaloSlope[len(logHaloSlope) - 1])
        for i in range(len(radius)):
            radius[i] = radius[i]/r0
        plt.subplot(5, 5, 5*rCount + densityCount + 1)
        plt.plot(radius, logHaloSlope, "c-")
fig.add_subplot(111, frameon=False)
plt.grid(False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel("Radius (scale lengths)")
plt.ylabel("logarithmic derivative of dark matter density")
plt.title("Slope of log(dark matter density) as per MDAR")
plt.tight_layout()
plt.savefig("D:/University/ICRAR Studentship/Dark matter by MDAR compiled.pdf")