# GalaxyModels
This repository contains code from my studentship at the International Centre for Radio Astronomy Research (ICRAR) in 2018/2019. My project, supervised by Professor Gerhardt Meurer involved building template galaxies based upon the "constant stability disk" or "constant Q" model. The hope was the constant stability disk model could be the physical basis underlying the mass discrepancy acceleration relation (MDAR). While this proved misguided, my implementation of the constant Q model is still valid and I'm able to present some very well calibrated Java code for building and exploring typical galaxies. I've also included my project report which explains the astrophysics behind these models and presents the results of my studentship. Section 6 of my studentship report has since been superseded by the file titled Dark matter via MDAR.pdf, which I wrote for an assignment for PHYS3003 at UWA. The latter document presents a more holistic and refined approach to that in section 6 of my studentship report.

I pursued a follow-up to the studentship during a research project of an astrophysics unit I studied in the next semester. For that purpose I built GalaxyModel14 which iteratively builds a galaxy within which MDAR holds strictly and the stable disk model also holds.

I have 14 models available overall. All but GalaxyModel14 implement the GalaxyModel.java interface I've provided. Whenever there are stars and gas in the galaxy, it is assumed the gas distribution responds to the stellar disk to maintain constant Q. The models' salient features are

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

GalaxyModel14: Same as GalaxyModel14 but the model's fields/parameters (e.g. scale length, rotation curve etc.) are iteratively improved so that MDAR holds strictly in the galaxy and the gas + stars form a constant Q disk. This model doesn't implement the GalaxyModel interface because its purpose is somewhat different. I should also caution that the iteration process almost certainly doesn't converge for vFinal < 125 km/s. This model is the first where I've update the numerical integration method to use Simpson's rule as opposed to less accurate, but simpler, methods.

Tester: This class contains the main method I used for running my models. It is completely customised to me and my computer. If other people want to run my models they should ideally write their own main method. At the very least, all the path names will need to be adjusted to reflect the location of files on the user's directory.

I've also included 5 python modules I used for visualising the output of methods in Tester.java. These files are currently not very well documented and quite personalised for my computer and the methods in Tester.java and GalaxyModel.java. MDAR plotting.py is used for visualising the output of the calculateAccelerations and MDAR methods of GalaxyModel and Tester respectively. One galaxy plotting.py is designed for recordData and oneGalaxy methods respectively. Dark matter.py was used to create figure 20 in Studentship report.pdf. Miscellaneous plotting.py contains a smorgasbord of different graphs I found worthwhile making (only figures 19 and 21 in the report are made using Miscellaneous plotting.py). The dark matter slopes.py is designed to produce the figures in the strict MDAR paper which also uses GalaxyModel14.java. In general the five python files are very poorly documented and difficult to use as of now.
