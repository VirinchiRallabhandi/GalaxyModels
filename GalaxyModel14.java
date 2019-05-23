import static java.lang.Math.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
/**
 * This class models a galaxy's mass distribution given its rotation velocity at
 * the edge of the disk. It does this by deriving the appropriate URC,
 * a stellar disk which is exponential radially and sech^2 vertically and a 3D gas
 * disk which maintains a constant Toomre stability parameter across the two disks.
 * The gas disk flares towards higher radii as per observations.
 * The dark matter distribution is therefore deduced from the URC and baryonic
 * matter distribution via Poisson's equation.
 * Any time a quantity is documented as surface ..., but refers to a 3D
 * quantity, it means volume density integrated over the z axis.
 * This class differs from GalaxyModel3 in that the gas disk in now 3D.
 * All units are SI unless otherwise stated.
 * @author Virinchi Rallabhandi
 *
 */
public class GalaxyModel14 {
	/**
	 * The constant value of the Toomre stability parameter
	 */
	public static final double Q = 1.6;
	/**
	 * The MOND acceleration constant
	 */
	public static final double gD = 1.2*pow(10, -10);
	public static final double G = 6.67408*pow(10, -11);
	public static final double Kpc = 3.086*pow(10, 19);
	public static final double Msun = 1.989*pow(10, 30);
	/**
	 * The set of radii (in ascending order) at which the evaluations are made.
	 * Each calculation is for a cylindrical rim or spherical shell of that radius.
	 * The values in this list are separated by 0.02 times the original stellar
	 * disk scale length. The 
	 */
	private ArrayList<Double> radius;
	/**
	 * The original maximum radius of the disk set by the requirement that
	 * the orbital period at the maximum radius is one 1 Gyr.
	 */
	private double rMax;
	/**
	 * The scale length of the stellar disk set in the first iteration through
	 * the program
	 */
	private double r0Original;
	/**
	 * The scale length of the stellar disk
	 */
	private double r0;
	/**
	 * The scale height of the stellar disk. By van der Kruit and Freeman 2011's
	 * notation, the vertical mass distribution is proportional to sech(z/z0/2)^2.
	 */
	private double z0;
	/**
	 * The central volume mass density of the stellar disk
	 */
	private double maxDensity;
	/**
	 * The maximum gas velocity dispersion. This is set, along with rg below,
	 * by the requirement that the average gas velocity dispersion is 8 km/s.
	 */
	private double sigma0;
	/**
	 * The scale length of gas velocity dispersion's exponential decay
	 */
	private double rG;
	/**
	 * The value which the stellar velocity dispersions, in any direction,
	 * asymptotically approach
	 */
	private double sigmaSAsymptote;
	/**
	 * rotationCurve.get(i) contains the rotation speed at radius.get(i)
	 */
	private ArrayList<Double> rotationCurve;
	/**
	 * The derivative of the rotation curve with respect to radius
	 */
	private ArrayList<Double> dVdr;
	private boolean function1;
	private boolean function2;
	private boolean function3;
	/**
	 * sechIntegrals is an i by 2 table. sechIntegrals[i][0] contains a value p
	 * while sechIntegrals[i][1] contains an estimate of Integrate[Sech[x]^(2*p),
	 * {x, 0, Infinity}]
	 */
	private double[][] sechIntegrals;
	
	/**
	 * @param vFinal The circular speed at the edge of the galaxy disk (the 
	 * disk will be truncated in a continuously differentiable manner
	 * past the edge)
	 * @param sechFile A full path to a text file which contains a table. The left
	 * column of the table should contain a value p while the right column
	 * should contain an estimate of Integrate[Sech[x]^(2*p), {x, 0, Infinity}] 
	 * @param function1 Whether or not the model should use the first MOND
	 * interpolating function - the one chosen by McGaugh, Lelli and Schombert
	 * (2016) - described in the companion paper
	 * @param function2 Whether or not the model should use the second MOND
	 * interpolating function - 1/2*(1 + sqrt(1 + 4/x)) - described in the
	 * companion paper
	 * @param function3 Whether or not the model should use the third MOND
	 * interpolating function - tanh(x^(1/3.5))^-1.75 - described in the 
	 * companion paper
	 */
	public GalaxyModel14(double vFinal, String sechFile, boolean function1,
		boolean function2, boolean function3) {
		if ((function1 && function2) || (function2 && function3) || (
			function1 && function3) || (!function1 && !function2 && !function3)) {
			throw new IllegalArgumentException("You must select exactly one "
				+ "MOND interpolating function.");
		}
		this.function1 = function1;
		this.function2 = function2;
		this.function3 = function3;
		//Most of the following scalings come from Wong et al. 2016,
		//Meurer et al. 2018 or van der Kruit and Freeman 2011.
		rMax = vFinal*pow(10, 9)*365.25*24*3600/2/PI;
		double mB = -10*log10(vFinal/1000) + 3;
		//Very rough estimate from expressions in Carrol and Ostlie and assuming
		//vMax = vFinal. These assumptions are only made so that a reasonable
		//initial guess is made for the parameter lambda.
		double luminosity = pow(10, -0.4*(mB - 5.31));
		double lambda = luminosity/6/pow(10, 10)*pow(67.4/50.0, 2);
		//Using H0 = 67.4 km/s/Mpc, the latest value by the Planck collaboration
		//Iteratively work out what the correct value of lambda is in the URC
		for (int i = 0; i < 2000; i = i + 1) {
			double vFinalPredicted = URC(rMax, lambda);
			lambda = lambda + (vFinal - vFinalPredicted)/2000000;
			//v0 in the URC is a strictly increasing function of lambda so
			//this process should work given enough iterations.
		}
		double vMax = 0;
		for (double r = 0; r < rMax; r = r + 0.001*rMax) {
			double v = URC(r, lambda);
			if (v > vMax) {
				vMax = v;
			}
		}
		double mR = -3.9 - 7.622*log10(vMax/1000);
		//The absolute magnitude of the galaxy in the R band
		double surfaceBrightness = pow(10, 5.3785 + 1.1757*log10(vMax/1000));
		//The effective surface brightness of the disk in Lsun/kpc^2
		r0 = sqrt(pow(10, -0.4*(mR - 4.61))/5.647/PI/surfaceBrightness)*Kpc;
		r0Original = r0;
		z0 = (0.45*vMax/100000 - 0.14)*Kpc;
		double stellarMass = pow(10, -0.4*(mR - 4.61))*pow(10,
			-1.578 - 0.0856*mR)*Msun;
		maxDensity = stellarMass/8/PI/z0/r0/r0;
		double r25 = (4.42631842*r0/Kpc - 0.55073662)*Kpc;
		//The optical radius for Ianjamasimanana's sample
		if (r25 < 0) {
			r25 = 10*r0;
			//r25 < 0 can only happen if r0 is very small. For dwarf galaxies,
			//the decay of the velocity dispersions is much smaller and
			//so it makes sense to have r25 very high so that rG is very high
			//and hence there is little decay.
		}
		rG = 3.7*r25;
		double sigmaGMean = 8000;
		//The mean gas velocity dispersion
		sigma0 = rMax*sigmaGMean/rG/(1 - pow(E, -rMax/rG));
		sechIntegrals = readSechFile(sechFile);
		sigmaSAsymptote = sigma0*pow(E, -rMax/rG)/sqrt(3);
		radius = new ArrayList<Double>();
		rotationCurve = new ArrayList<Double>();
		dVdr = new ArrayList<Double>();
		for (double r = 0; r < rMax; r = r + r0/50) {
			radius.add(r);
			rotationCurve.add(URC(r, lambda));
			dVdr.add(dURCdr(r, lambda));
		}
		double[] parameters = getTruncations();
		double rFinalS = parameters[4];
		double rFinalG = parameters[5];
		radius = new ArrayList<Double>();
		rotationCurve = new ArrayList<Double>();
		dVdr = new ArrayList<Double>();
		for (double r = 0; r < max(rFinalS, max(rFinalG, rMax)); r = r + r0/50) {
			radius.add(r);
			rotationCurve.add(URC(r, lambda));
			dVdr.add(dURCdr(r, lambda));
		}
	}
	
	/** 
	 * Iteratively calculates the total acceleration due to the baryons and hence
	 * the total observed acceleration as determined by the chosen MOND 
	 * interpolating function. The rotation curve and its derivative are
	 * updated at each iteration.
	 * @return A list of three arrays - the first is the list of radius values,
	 * the second is the rotation speed assumed at each radius and the third
	 * is list of gBar = gradPotentialBaryon at each radius.
	 */
	public double[][] calculateAccelerations() {
		int maxIteration = 4;
		for(int iteration = 0; iteration <= maxIteration;
			iteration = iteration + 1) {
			System.out.println(iteration);
			double[] parameters = getTruncations();
			double aS = parameters[0];
			double bS = parameters[1];
			double aG = parameters[2];
			double bG = parameters[3];
			double[] height = new double[101];
			//The height above the z axis/midplane for the stellar disk
			double[][] stellarDensity = new double[radius.size()][height.length];
			//The stellar disk's volume mass density at each radius and height
			double[][] gasDensity = new double[radius.size()][height.length];
			//The gas disk's volume mass density at each radius and height
			//given constant Q
			double[] gradPotentialBaryon = new double[radius.size()];
			//The acceleration due to the "baryons" alone
			double[] gradPotentialTotal = new double[radius.size()];
			//The acceleration due to all the mass as determined by the URC
			double[] newRotationCurve = processModel(height, stellarDensity,
				gasDensity, gradPotentialBaryon, gradPotentialTotal, aS,
				bS, aG, bG);
			if (iteration == maxIteration) {
				double[] radiusValues = new double[radius.size()];
				double[] rotationSpeeds = new double[rotationCurve.size()];
				for (int i = 0; i < radius.size(); i = i + 1) {
					radiusValues[i] = radius.get(i);
					rotationSpeeds[i] = rotationCurve.get(i);
				}
				return new double[][] {radiusValues, rotationSpeeds,
					gradPotentialBaryon};
			}
			reInitialise(newRotationCurve);
		}
		return null;
		//It is not actually possible to reach this last line. The return statement
		//earlier will always be executed first.
	}
	
	/**
	 * Gives the parameters required to extend the disk beyond rMax in
	 * a continuously differentiable but truncatable way.
	 * @return aS, bS, aG, bG, rFinalS and rFinalG in that order
	 * The stellar disk's surface mass density is extended as aS + bSe^(r/r0)
	 * beyond rMax. The equivalent expression for the gas is aG + bGe^(r/r0)
	 * rFinalS and rFinalG are the radii at which the stellar and gas surface
	 * mass densities drop to zero respectively.
	 */
	private double[] getTruncations() {
		//To truncate the disk in a continuously differentiable manner,
		//I fit a + be^(r/r0) to the function, with b negative, to
		//the calculated values beyond rMax.
		//These densities refer to surface mass densities - i.e. volume
		//density integrated over the z axis.
		int j = 0;
		while (j < radius.size() - 1 && radius.get(j) < rMax) {
			j = j + 1;
		}
		double starDensity2 = 4*z0*maxDensity*pow(E, -1*rMax/r0);
		double starDensity1 = 4*z0*maxDensity*pow(E, -1*(rMax - 0.02*r0)/r0);
		double sigmaG2 = sigma0*pow(E, -1*rMax/rG);
		double sigmaG1 = sigma0*pow(E, -1*(rMax - 0.02*r0)/rG);
		double v2 = rotationCurve.get(j);
		double k2 = v2/rMax*sqrt(2*(1 + rMax/v2*dVdr.get(j)));
		double sigmaSR2 = 1.67*sqrt(2*PI*G*starDensity2*z0) + sigmaSAsymptote;
		double gasDensity2 = sigmaG2*k2/PI/G/Q - sigmaG2/sigmaSR2*starDensity2;
		double v1 = rotationCurve.get(j - 1);
		double k1 = v1/(rMax - 0.02*r0)*sqrt(2*(1 + (rMax - 0.02*r0)/v1*dVdr.get(
			j - 1)));
		double sigmaSR1 = 1.67*sqrt(2*PI*G*starDensity1*z0) + sigmaSAsymptote;
		double gasDensity1 = sigmaG1*k1/PI/G/Q - sigmaG1/sigmaSR1*starDensity1;
		double aS = starDensity2 - (starDensity2 - starDensity1)/0.02;
		double bS = pow(E, -rMax/r0)*(starDensity2 - starDensity1)/0.02;
		double aG = gasDensity2 - (gasDensity2 - gasDensity1)/0.02;
		double bG = pow(E, -rMax/r0)*(gasDensity2 - gasDensity1)/0.02;
		double rFinalS = r0*log(abs(aS/bS));
		double rFinalG = r0*log(abs(aG/bG));
		return new double[] {aS, bS, aG, bG, rFinalS, rFinalG};
	}
	
	/**
	 * Gives a linear extension/tangent to the gas density profiles for radii
	 * less than radius[i]. This method assumes i < radius.length - 1.
	 * @param i the index in the radius array below which the gas density profiles
	 * should be linearly extended.
	 * @param radius The radii in ascending order of all the points used for
	 * the numerical methods of this class
	 * @return a 2 element array. The gas surface density should be 
	 * return[0]*r + return[1] for r <= radius[i]
	 */
	private double[] getInnerDensities(int i) {
		double sigmaG1 = sigma0*pow(E, -1*radius.get(i)/rG);
		double stellarDensity1 = 4*z0*maxDensity*pow(E, -1*radius.get(i)/r0);
		double v1 = rotationCurve.get(i);
		double k1 = v1/radius.get(i)*sqrt(2*(1 + radius.get(i)/v1*dVdr.get(i)));
		double sigmaSR1 = 1.67*sqrt(2*PI*G*stellarDensity1*z0) + sigmaSAsymptote;
		double gasDensity1 = sigmaG1*k1/PI/G/Q - sigmaG1/sigmaSR1*
			stellarDensity1;
		double sigmaG2 = sigma0*pow(E, -1*radius.get(i + 1)/rG);
		double stellarDensity2 = 4*z0*maxDensity*pow(E, -1*radius.get(i + 1)/r0);
		double v2 = rotationCurve.get(i);
		double k2 = v2/radius.get(i + 1)*sqrt(2*(1 + radius.get(i + 1)/v2*dVdr.get(
			i + 1)));
		double sigmaSR2 = 1.67*sqrt(2*PI*G*stellarDensity2*z0) + sigmaSAsymptote;
		double gasDensity2 = sigmaG2*k2/PI/G/Q - sigmaG2/sigmaSR2*
			stellarDensity2;
		double m = (gasDensity2 - gasDensity1)/(radius.get(i + 1) - radius.get(i));
		double c = gasDensity1 - radius.get(i)*m;
		return new double[] {m, c};
	}
	
	/**
	 * Re-calculates the fields based on the new rotation curve calculated
	 * at the end of a given iteration.
	 * @param newRotationCurve The rotation curve that should result from strict
	 * adherence to MOND and the current configuration of stars and gas.
	 */
	private void reInitialise(double[] newRotationCurve) {
		for (int i = 0; i < newRotationCurve.length; i = i + 1) {
			rotationCurve.set(i, 0.5*(rotationCurve.get(i) + newRotationCurve[i]));
			//The rotation curve is adjusted to take the average between its
			//new and old values.
		}
		int j = 0;
		while (j < radius.size() - 1 && radius.get(j) < rMax) {
			j = j + 1;
		}
		double vFinal = rotationCurve.get(j);
		//Most of the following scalings come from Wong et al. 2016,
		//Meurer et al. 2018 or van der Kruit and Freeman 2011.
		rMax = vFinal*pow(10, 9)*365.25*24*3600/2/PI;
		double vMax = 0;
		for (int i = 0; i <= j; i = i + 1) {
			double v = rotationCurve.get(i);
			if (v > vMax) {
				vMax = v;
			}
		}
		double mR = -3.9 - 7.622*log10(vMax/1000);
		//The absolute magnitude of the galaxy in the R band
		double surfaceBrightness = pow(10, 5.3785 + 1.1757*log10(vMax/1000));
		//The effective surface brightness of the disk in Lsun/kpc^2
		r0 = sqrt(pow(10, -0.4*(mR - 4.61))/5.647/PI/surfaceBrightness)*Kpc;
		z0 = (0.45*vMax/100000 - 0.14)*Kpc;
		double stellarMass = pow(10, -0.4*(mR - 4.61))*pow(10,
			-1.578 - 0.0856*mR)*Msun;
		maxDensity = stellarMass/8/PI/z0/r0/r0;
		double r25 = (4.42631842*r0/Kpc - 0.55073662)*Kpc;
		//The optical radius for Ianjamasimanana's sample
		if (r25 < 0) {
			r25 = 10*r0;
			//r25 < 0 can only happen if r0 is very small. For dwarf galaxies,
			//the decay of the velocity dispersions is much smaller and
			//so it makes sense to have r25 very high so that rG is very high
			//and hence there is little decay.
		}
		rG = 3.7*r25;
		double sigmaGMean = 8000;
		//The mean gas velocity dispersion
		sigma0 = rMax*sigmaGMean/rG/(1 - pow(E, -rMax/rG));
		sigmaSAsymptote = sigma0*pow(E, -rMax/rG)/sqrt(3);
		double[] truncations = getTruncations();
		double rFinalS = truncations[4];
		double rFinalG = truncations[5];
		double newRFinal = max(rFinalS, rFinalG);
		while (newRFinal > radius.get(radius.size() - 1)) {
			radius.add(radius.get(radius.size() - 1) + r0Original/50);
			rotationCurve.add(rotationCurve.get(rotationCurve.size() - 1));
			//The actual value of the rotation speed doesn't matter too much this
			//far out because we're in the truncation regime where I'm just
			//adding in an analytic extension to the disk dependent only the
			//rotation curve before truncation. Even if the rotation curve value
			//is required, making the rotation curve flat like this isn't
			//unrealistic; we should have reached the flat part of the rotation
			//curve by now.
		}
		dVdr = new ArrayList<Double>();
		for (int i = 0; i < radius.size() - 1; i = i + 1) {
			dVdr.add((rotationCurve.get(i + 1) - rotationCurve.get(
				i))/(radius.get(i + 1) - radius.get(i)));
		}
		dVdr.add(dVdr.get(dVdr.size() - 1));
		//The last line make radius, rotationCurve and dVdr all the same length.
	}
	
	/**
	 * This method does the majority of the actually important work - namely
	 * calculating the distribution of stars and gas based on the current
	 * values of the fields and the consequent accelerations.
	 * All the arrays will be re-written with their correct values for
	 * this model galaxy so it is irrelevant how one initializes them.
	 * The dimensions of the arrays are assumed to be properly related.
	 * @param stellarDensity The stellar disk's surface mass density
	 * stellarDensity[i][j] = volume mass density at radius[i] and height[j]
	 * @param gasDensity The surface mass density of the gas disk at each
	 * radius given constant Q
	 * @param rotationCurve The URC velocity at each radius
	 * @param gradPotentialGas The acceleration due to the gas alone at each radius
	 * @param gradPotentialStar The acceleration due to the stellar disk alone
	 * @param gradPotentialBaryon The acceleration due to the baryons alone
	 * @param gradPotentialTotal The total acceleration
	 * @param lower gradPotential*[i] should correspond to radius[i + lower]
	 * @return The rotation curve that would result from strict adherence to MOND
	 * and an acceleration due to baryons determined by the current configuration
	 * of stars and gas.
	 */
	private double[] processModel(double[] height,
		double[][] stellarDensity, double[][] gasDensity,
		double[] gradPotentialBaryon, double[] gradPotentialTotal,
		double aS, double bS, double aG, double bG) {
		for (int i = 0; i < height.length; i = i + 1) {
			height[i] = ((double) i)*z0*6/101;
		}
		//First get the parameters of the linear fit to the interior of 
		//the gas density function
		//Assume radius.length >= 52
		//These are all surface mass densities
		double[] linearFit = getInnerDensities(50);
		double m = linearFit[0];
		double c = linearFit[1];
		for (int i = 0; i < radius.size(); i = i + 1) {
			double gasSurfaceDensity = 0;
			if (radius.get(i) <= rMax) {
				double starSurfaceDensity = 4*z0*maxDensity*pow(E,
					-1*radius.get(i)/r0);
				//The surface mass density of the stars at radius.get(i)
				for (int j = 0; j < height.length; j = j + 1) {
					stellarDensity[i][j] = 4*maxDensity*pow(E, -1*radius.get(
						i)/r0)*pow(E, height[j]/z0)/pow(pow(E, height[j]/z0) + 1, 2);
				}
				if (i > 50) {
					double sigmaG = sigma0*pow(E, -1*radius.get(i)/rG);
					double v = rotationCurve.get(i);
					double k = v/radius.get(i)*sqrt(2*(1 + radius.get(
						i)/v*dVdr.get(i)));
					double sigmaSR = 1.67*sqrt(2*PI*G*starSurfaceDensity*z0) 
						+ sigmaSAsymptote;
					//The stellar radial velocity dispersion
					gasSurfaceDensity = sigmaG*k/PI/G/Q - sigmaG/sigmaSR*
						starSurfaceDensity;
				} else {
					gasSurfaceDensity = m*radius.get(i) + c;
				}
			} else {
				for (int j = 0; j < height.length; j = j + 1) {
					stellarDensity[i][j] = (aS + bS*pow(E, radius.get(i)/r0))/z0*pow(
						E, height[j]/z0)/pow(pow(E, height[j]/z0) + 1, 2); 
				}
				gasSurfaceDensity = aG + bG*pow(E, radius.get(i)/r0);
			}
			if (gasSurfaceDensity < 0) {
				gasSurfaceDensity = 0;
		        //Gas is not necessarily found everywhere within a galaxy.
				//If negative gas density is required to maintain constant Q,
				//then it means there is really no gas there at all and the
				//two fluid Q is inappropriate.
			}
			double starSurfaceDensity = 4*z0*maxDensity*pow(E, -1*radius.get(i)/r0);
			double sigmaG = sigma0*pow(E, -1*radius.get(i)/rG);
			double sigmaSZ = sqrt(2*PI*G*starSurfaceDensity*z0) + sigmaSAsymptote;
			//The stellar velocity dispersion in the z direction
			double p = 3*pow(sigmaSZ/sigmaG, 2);
			//The p parameter of van der Kruit (1981)
			double integral = 0;
			for (int j = 0; j < sechIntegrals.length; j = j + 1) {
				//Maybe I'll change this to a binary search later. It seems
				//irrelevant which search procedure is used as the numerical
				//integral later is the one taking the bulk of the time.
				if (j == 0 && p <= sechIntegrals[0][0]) {
					p = sechIntegrals[0][0];
					integral = sechIntegrals[0][1];
				} else {
					if (j == sechIntegrals.length - 1 && p >= sechIntegrals[j][0]) {
						p = sechIntegrals[j][0];
						integral = sechIntegrals[j][1];
					} else {
						if (sechIntegrals[j][0] <= p && p < sechIntegrals[j + 1][0]) {
							p = sechIntegrals[j][0];
							integral = sechIntegrals[j][1];
						}
					}
				}
			}
			double w = 1.7*pow(p, -0.5)*sigmaSZ/pow(2*PI*G*maxDensity*pow(E,
				-1*radius.get(i)/r0), 0.5);
			//The full width at half maximum of the vertical gas distribution
			double z1 = w/2/log(pow(2, 1/2/p) + sqrt(pow(2, 1/p) - 1));
			//The scale height of the gas distribution. While the stellar vertical
			//distribution is sech(z/z0/2)^2, the gas vertical distribution is
			//sech(z/z1)^(2*p)
			double rho0 = gasSurfaceDensity/2/z1/integral;
			//The volume mass density at radius[i] in the galactic midplane
			for (int j = 0; j < height.length; j = j + 1) {
				gasDensity[i][j] = rho0*pow(2*pow(E, height[j]/z1)/(pow(
					E, 2*height[j]/z1) + 1), 2*p);
				if (stellarDensity[i][j] < 0) {
					stellarDensity[i][j] = 0;
					//This can happen if radius[i] is beyond the truncation
					//radius of the stellar disk.
				}
			}
		}
		//The process below is essentially a numerical integration via Simpson's
		//rule. There's three loops as this is a triple integral to find the
		//acceleration due to gas and stars.
		int rMax = radius.size() - 1;
		if (rMax % 2 == 0) {
			rMax = rMax - 1;
		}
		for (int i = 0; i < gradPotentialBaryon.length; i = i + 1) {
			double R = radius.get(i);
			double dr = radius.get(1) - radius.get(0);
			double dtheta = PI/100;
			double dz = height[1] - height[0];
			double totalSum = 0;
			for (int j = 1; j <= rMax; j = j + 1) {
				double r = radius.get(j);
				double angleSum = 0;
				for (int k = 0; k <= 100; k = k + 1) {
					double theta = ((double) k)*100/PI;
					double zSum = 0;
					for (int n = 0; n < height.length; n = n + 1) {
						double z = height[n];
						if (r != R || k != 0 || z != 0) {
							double zAdd = (
								stellarDensity[j][n] + gasDensity[j][n])/pow(
								z*z + r*r + R*R - 2*r*R*cos(theta), 1.5);
							if (n % 2 == 1) {
								zAdd = zAdd*4;
							} else if (n !=0 && n != height.length - 1) {
								zAdd = zAdd*2;
							}
							zSum = zSum + zAdd;
						}
					}
					double thetaAdd = zSum*(R - r*cos(theta));
					if (k % 2 == 1) {
						thetaAdd = thetaAdd*4;
					} else if (k != 0 && k != 100) {
						thetaAdd = thetaAdd*2;
					}
					angleSum = angleSum + thetaAdd;
				}
				double rAdd = 0;
				if (j == 1 || j == rMax) {
					rAdd = angleSum*r;
				} else if (j % 2 == 0) {
					rAdd = 4*angleSum*r;
				} else {
					rAdd = 2*angleSum*r;
				}
				totalSum = totalSum + rAdd;
			}
			gradPotentialBaryon[i] = totalSum*4*G*dz*dtheta*dr/27;
		}
		double[] newRotationCurve = new double[gradPotentialBaryon.length];
		for (int i = 0; i < gradPotentialBaryon.length; i = i + 1) {
			if (radius.get(i) == 0) {
				gradPotentialTotal[i] = 0;
				gradPotentialBaryon[i] = 0;
			} else {
				if (function1) {
					gradPotentialTotal[i] = gradPotentialBaryon[i]*f1(
						gradPotentialBaryon[i]/gD);
				} else if (function2) {
					gradPotentialTotal[i] = gradPotentialBaryon[i]*f2(
						gradPotentialBaryon[i]/gD);
				} else if (function3) {
					gradPotentialTotal[i] = gradPotentialBaryon[i]*f3(
						gradPotentialBaryon[i]/gD);
				}
			}
			newRotationCurve[i] = sqrt(radius.get(i)*gradPotentialTotal[i]);
		}
		return newRotationCurve;
	}
	
	private double f1(double x) {
		return 1/(1 - pow(E, -1*sqrt(x)));
	}
	
	private double f2(double x) {
		return (1 + sqrt(1 + 4/x))/2;
	}
	
	private double f3(double x) {
		return pow(tanh(pow(x, 1/3.5)), -1.75);
	}
	
	/**
	 * @param r The galocentric radius
	 * @param lambda the lambda parameter controlling the URC as formulated in
	 * the "note added in proof" of Persic, Salucci and Stel (1996)
	 * @return The velocity at a radius r for this model galaxy based
	 * on the URC of Persic, Salucci and Stel.
	 */
	private double URC(double r, double lambda) {
		r = r*75/70;
		//A correction for the different Hubble constants used by different papers
		double rOpt = 13*pow(lambda, 0.5);
		//rOpt is now in units of kpc
		r = r/Kpc;
		//r is now also in units of kpc
		double v0 = sqrt(0.8 + 0.49*log10(lambda) + 0.75*pow(E, -0.4*lambda)/(
			0.47 + 2.25*pow(lambda, 0.4)));
		v0 = 200*pow(lambda, 0.41)/v0;
		double v = (0.72 + 0.44*log10(lambda))*1.97*pow(r/rOpt, 1.22)/pow(
			pow(r/rOpt, 2) + 0.61, 1.43);
		v = v + 1.6*pow(E, -0.4*lambda)*pow(r/rOpt, 2)/(pow(r/rOpt, 2)
			+ 2.25*pow(lambda, 0.4));
		return v0*sqrt(v)*1000;
	}
	
	/**
	 * @param r The galocentric radius
	 * @param lambda the lambda parameter controlling the URC as formulated in
	 * the "note added in proof" of Persic, Salucci and Stel (1996)
	 * @return The derivative of the rotation curve with respect to radius
	 */
	private double dURCdr(double r, double lambda) {
		r = r*75/70;
		//A correction for the different Hubble constants used by different papers
		double rOpt = 13*pow(lambda, 0.5);
		//rOpt is now in units of kpc
		r = r/Kpc;
		//r is now also in units of kpc
		double v0 = sqrt(0.8 + 0.49*log10(lambda) + 0.75*pow(E, -0.4*lambda)/(
			0.47 + 2.25*pow(lambda, 0.4)));
		v0 = 200*pow(lambda, 0.41)/v0;
		double x = r/rOpt;
		double a = 0.72 + 0.44*log10(lambda);
		double b = 1.6*pow(E, -0.4*lambda);
		double c = 2.25*pow(lambda, 0.4);
		double dv2dx = 2.4034*a*pow(x, 0.22)*pow(x*x + 0.61, 1.43);
		dv2dx = dv2dx - 5.6342*a*pow(x, 2.22)*pow(x*x + 0.61, 0.43);
		dv2dx = dv2dx/pow(x*x + 0.61, 2.86) + 2*b*x*c/pow(x*x + c, 2);
		return dv2dx*v0*v0/2/URC(r*Kpc, lambda)/rOpt*1000*1000/Kpc;
	}
	
	private double[][] readSechFile(String sechFile) {
		try {
			ArrayList<String> lines = new ArrayList<String>();
			String line = "";
			File f = new File(sechFile);
			FileReader reader = new FileReader(f);
			BufferedReader br = new BufferedReader(reader);
			while (line != null) {
				if (!line.equals("")) {
					lines.add(line);
				}
				line = br.readLine();
			}
			br.close();
			reader.close();
			double[][] values = new double[lines.size()][2];
			for (int i = 0; i < lines.size(); i = i + 1) {
				String[] numbers = lines.get(i).split(" ");
				values[i][0] = Double.parseDouble(numbers[0]);
				values[i][1] = Double.parseDouble(numbers[numbers.length - 1]);
			}
			//Assume the file is sorted by ascending order of p
			//Some sorting may be added later.
			return values;
		} catch (Exception e) {
			System.out.println("A disaster occurred.");
			return null;
		}
	}
}