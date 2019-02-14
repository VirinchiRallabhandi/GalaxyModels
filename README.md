# GalaxyModels
This repository contains code from my studentship at the International Centre for Radio Astronomy Research (ICRAR) in 2018/2019. My project, supervised by Professor Gerhardt Meurer involved building template galaxies based upon the "constant stability disk" or "constant Q" model. The hope was the constant stability disk model could be the physical basis underlying the mass discrepancy acceleration relation (MDAR). While this proved misguided, my implementation of the constant Q model is still valid and I'm able to present some very well calibrated Java code for building and exploring typical galaxies. I've also included my project report which explains the astrophysics behind these models and presents the results of my studentship.

I have 13 models available overall. All 13 implement the GalaxyModel.java interface I've provided. Whenever there are stars and gas in the galaxy, it is assumed the gas distribution responds to the stellar disk to maintain constant Q. The models' salient features are

GalaxyModel1: 2D pure exponential stellar disk and 2D gas disk, uses Hubble type to calibrate luminosity and the universal rotation curve (URC)

GalaxyModel2: Same as GalaxyModel1 except the stellar disk is now 3D - the volume mass density takes the canonical form of decaying exponentially radially and decaying by sech^2 vertically

GalaxyModel3: Same as GalaxyModel2 but with a better implementation of URC (the "note added in proof" of Persic, Salucci and Stel's 1996 paper) which no longer requires Hubble type

GalaxyModel4: Same as GalaxyModel3 except the gas disk is also now 3D allowing for more realistic forces in the outskirts of galaxies where the gas disk tends to flare. This is by the far the best model out of the 13 and should be prefered for modelling galaxies with vMax > 70 km/s.

GalaxyModel5: This class is the only one using the Romeo and Wiegert expression for Q rather than the Wang and Silk expression. The Romeo and Wiegert expression, being piecewise in definition doesn't turn out well for modelling a gas distribution in these models.

GalaxyModel6: This is the only class where the dark matter distribution is not assumed to adjust itself to reproduce the rotation curve. The dark matter distribution is fit to a cored isothermal sphere which is then iteratively improved with the rotation curve to achieve self consistency. The fitting procedure in this class doesn't yet work and requires much improvement.

GalaxyModel7: Same as GalaxyModel4 except the URC is assumed is very similar to GalaxyModel 1-3 but without relying on Hubble type. This class didn't really resolve any issue with the URC of GalaxyModel1.

GalaxyModel8: This class models a galaxy with a constant Q gas disk alone - no stars. Its purpose is to isolate the influence of the constant stability disk on MDAR. 

GalaxyModel9: Same as GalaxyModel4 but with the Karukes and Salucci URC tailored to dwarf disk galaxies. This model should be used only for galaxies with vMax < 70 km/s, a vMax range not necessarily well modelled by GalaxyModel4. Honestly, I didn't fully understand Karukes and Salucci's paper and the implementation in this class gives eccentric results.

GalaxyModel10: Same as GalaxyModel9, but with a rotation curve model I fit myself to Karukes and Salucci's data. This class is the best of the 13 for modelling galaxies with vMax < 70 km/s.

GalaxyModel11 and GalaxyModel12: Both are unsuccessful attempts to patch the two regimes of the URC (vMax > and < 70 km/s) together using Courteau's tan^-1 model of rotation curves.

GalaxyModel13: Same as GalaxyModel4 except the stellar mass density now follows a Sersic, rather than exponential, profile radially. This model works as well as GalaxyModel4 except the user must choose the Sersic index. 

Tester: This class contains the main method I used for running my models. It is completely customised to me and my computer. If other people want to run my models they should ideally write their own main method. At the very least, all the path names will need to be adjusted to reflect the location of files on the user's directory.
