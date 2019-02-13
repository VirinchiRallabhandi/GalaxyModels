import static java.lang.Math.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
/**
 * This class models a galaxy's mass distribution given its maximum rotation velocity.
 * It does this by deriving a stellar disk which is exponential radially and 
 * sech^2 vertically and a 3D gas disk which maintains a constant Toomre stability
 * parameter across the two disks. The gas disk flares towards higher
 * radii as per observations. The dark matter distribution is therefore deduced
 * from the rotation curve and baryonic matter distribution via Poisson's equation.
 * This class uses an arctan model for the rotation curve. The rotation curve's
 * scale length, rT, is a cubic function of vMax. The cubic arises from trying to
 * fit an arctan model to the spiral galaxy URC. While the fit is quite strong, the
 * extrapolation to vMax lower than the spiral vMax is questionable. 
 * Any time a quantity is documented as surface ..., but refers to a 3D
 * quantity, it means volume density integrated over the z axis.
 * All units are SI unless otherwise stated.
 * @author Virinchi Rallabhandi
 *
 */
public class GalaxyModel12 implements GalaxyModel {
	/**
	 * The constant value of the Toomre stability parameter
	 */
	public static final double Q = 2;
	public static final double G = 6.67408*pow(10, -11);
	public static final double Kpc = 3.086*pow(10, 19);
	public static final double Msun = 1.989*pow(10, 30);
	/**
	 * The maximum radius of the disk set by the requirement that
	 * the orbital period at the maximum radius is one 1 Gyr.
	 */
	private double rMax;
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
	 * The maximum rotational speed of any object in the galaxy
	 */
	private double vMax;
	/**
	 * The scale length of the rotation curve
	 */
	private double rT;
	/**
	 * sechIntegrals is an i by 2 table. sechIntegrals[i][0] contains a value p
	 * while sechIntegrals[i][1] contains an estimate of Integrate[Sech[x]^(2*p),
	 * {x, 0, Infinity}]
	 */
	private double[][] sechIntegrals;
	
	/**
	 * @param vMax The maximum rotational velocity of an object orbiting within
	 * the disk's truncation radius.
	 * @param sechFile A full path to a text file which contains a table. The left
	 * column of the table should contain a value p while the right column
	 * should contain an estimate of Integrate[Sech[x]^(2*p), {x, 0, Infinity}] 
	 */
	public GalaxyModel12(double vMax, String sechFile) {
		this.vMax = vMax;
		double vkm = vMax/1000;
		rT = (2.65704169*pow(10, -7)*pow(vkm, 3) - 1.67503665*pow(10, -4)*pow(
			vkm, 2) + 3.76942491*pow(10, -2)*vkm - 1.05643153)*Kpc;
		//Most of the following scalings come from Wong et al. 2016,
		//Meurer et al. 2018 or van der Kruit and Freeman 2011.
		rMax = vMax*pow(10, 9)*365.25*24*3600/2/PI;
		//As this a dwarf galaxy, the rotation curve will be strictly increasing
		//within observational limits and hence V(rMax) = vMax 
		double mR = -3.9 - 7.622*log10(vMax/1000);
		//The absolute magnitude of the galaxy in the R band
		double surfaceBrightness = pow(10, 5.3785 + 1.1757*log10(vMax/1000));
		//The effective surface brightness of the disk in Lsun/kpc^2
		r0 = sqrt(pow(10, -0.4*(mR - 4.61))/5.647/PI/surfaceBrightness)*Kpc;
		z0 = r0/7.3;
		//Using the "average flattening ratio" as quoted in Zheng et al. 2013.
		//The van der Kruit and Freeman expression which has z0 as a function
		//of vMax only gives positive values for vMax > 31 km/s
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
	}
	
	/** 
	 * Calculates the total accelerations (gradient of the gravitational potential)
	 * due to each galaxy component - gas, stellar disk and dark matter halo -
	 * assuming the URC with parameters a/alpha and beta.
	 * @param start The percentage of the lower radius range whose acceleration
	 * should not be evaluated.
	 * @param end The percentage of the upper radius range whose acceleration should
	 * not be evaluated.
	 * e.g. If lower = 0.25 and upper = 0.75 then only the middle 50% of the
	 * radius range of the galaxy will have its acceleration evaluated.
	 * @return A two row table. return[0] contains log(gBar) at radii defined
	 * by the start and end values while return[1] contains log(gObs) at
	 * radii defined by the start and end values.
	 */
	public double[][] calculateAccelerations(double start, double end) {
		double[] parameters = getTruncations();
		double aS = parameters[0];
		double bS = parameters[1];
		double aG = parameters[2];
		double bG = parameters[3];
		double rFinalS = parameters[4];
		double rFinalG = parameters[5];
		double[] radius = new double[(int) ((max(
			max(rMax, rFinalS), rFinalG))/r0*50)];
		//The radius of each cylindrical rim or spherical shell at which
		//the summations will be performed
		double[] height = new double[100];
		//The height above the z axis/midplane for the stellar disk
		double[][] stellarDensity = new double[radius.length][height.length];
		//The stellar disk's volume mass density at each radius and height
		double[][] gasDensity = new double[radius.length][height.length];
		//The gas disk's volume mass density at each radius and height
		//given constant Q
		double[] rotationCurve = new double[radius.length];
		//The URC velocity at each radius
		int lower = (int) (((double) radius.length)*start);
		int upper = (int) (((double) radius.length)*end - 1);
		double[] gradPotentialGas = new double[upper - lower + 1];
		//The acceleration due to the gas alone at each radius on the midplane
		double[] gradPotentialStar = new double[gradPotentialGas.length];
		//The acceleration due to the stellar disk alone 
		double[] gradPotentialBaryon = new double[gradPotentialGas.length];
		//The acceleration due to the "baryons" alone
		double[] gradPotentialTotal = new double[gradPotentialGas.length];
		//The acceleration due to all the mass as determined by the URC
		processModel(radius, height, stellarDensity, gasDensity, rotationCurve,
			gradPotentialGas, gradPotentialStar, gradPotentialBaryon,
			gradPotentialTotal, lower, aS, bS, aG, bG);
		double[][] results = new double[2][gradPotentialTotal.length];
		for (int i = 0; i < gradPotentialTotal.length; i = i + 1) {
			results[0][i] = log10(gradPotentialBaryon[i]);
			results[1][i] = log10(gradPotentialTotal[i]);
		}
		return results;
	}
	
	/**
	 * Calculates the total accelerations (gradient of the gravitational potential)
	 * due to each galaxy component - gas, stellar disk and dark matter halo -
	 * assuming the URC with parameters a/alpha and beta.
	 * @param pathname The full path and filename to which the this galaxy's
	 * model results are to be logged to.
	 * The results will be in essentially a table format with each row listing
	 * a radius followed by the stellar surface density, gas surface density,
	 * dark matter volume density, gradient of the gravitational potential = 
	 * acceleration for stars, gas, all baryons, dark matter and all matter.
	 */
	public void recordData(String pathname) {
		double[] parameters = getTruncations();
		double aS = parameters[0];
		double bS = parameters[1];
		double aG = parameters[2];
		double bG = parameters[3];
		double rFinalS = parameters[4];
		double rFinalG = parameters[5];
		double[] radius = new double[(int) ((max(
			max(rMax, rFinalS), rFinalG))/r0*50)];
		//The radius of each cylindrical rim or spherical shell at which
		//the summations will be performed
		double[] height = new double[100];
		//The height above the z axis/midplane for the stellar disk
		double[][] stellarDensity = new double[radius.length][height.length];
		//The stellar disk's volume mass density at each radius and height
		double[][] gasDensity = new double[radius.length][height.length];
		//The gas disk's volume mass density at each radius and height
		//given constant Q
		double[] rotationCurve = new double[radius.length];
		//The URC velocity at each radius
		double[] gradPotentialGas = new double[radius.length];
		//The acceleration due to the gas alone at each radius
		double[] gradPotentialStar = new double[radius.length];
		//The acceleration due to the stellar disk alone 
		double[] gradPotentialBaryon = new double[radius.length];
		//The acceleration due to the "baryons" alone
		double[] gradPotentialTotal = new double[radius.length];
		//The acceleration due to all the mass as determined by the URC
		processModel(radius, height, stellarDensity, gasDensity, rotationCurve,
			gradPotentialGas, gradPotentialStar, gradPotentialBaryon,
			gradPotentialTotal, 0, aS, bS, aG, bG);
		double[] gradPotentialHalo = new double[radius.length];
		for (int i = 0; i < gradPotentialHalo.length; i = i + 1) {
			gradPotentialHalo[i] = gradPotentialTotal[i] - gradPotentialBaryon[i];
		}
		double[] haloDensity = new double[radius.length];
		//We calculate the halo density by Poisson's equation
		for (int i = 0; i < haloDensity.length - 1; i = i + 1) {
			if (radius[i] != 0) {
				haloDensity[i] = gradPotentialHalo[i]/2/PI/G/radius[i] + (
					gradPotentialHalo[i + 1] - gradPotentialHalo[i])/(
					radius[i + 1] - radius[i])/4/PI/G;
				if (haloDensity[i] < 0) {
					haloDensity[i] = 0;
				}
			}
		}
		haloDensity[0] = haloDensity[1];
		//The Laplacian in spherical coordinates is indeterminate at r = 0
		haloDensity[haloDensity.length - 1] = haloDensity[haloDensity.length - 2];
		ArrayList<String> data = new ArrayList<String>();
		String space = "   ";
		for (int i = 0; i < radius.length; i = i + 1) {
			double sigmaG = sigma0*pow(E, -1*radius[i]/rG);
			double starSurfaceDensity = 4*z0*stellarDensity[i][0];
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
					-1*radius[i]/r0), 0.5);
			//The full width at half maximum of the vertical gas distribution
			double z1 = w/2/log(pow(2, 1/2/p) + sqrt(pow(2, 1/p) - 1));
			//The scale height of the gas distribution. While the stellar vertical
			//distribution is sech(z/z0/2)^2, the gas vertical distribution is
			//sech(z/z1)^(2*p)
			StringBuilder line = new StringBuilder();
			line.append(radius[i]);
			line.append(space);
			line.append(stellarDensity[i][0]*4*z0);
			line.append(space);
			line.append(gasDensity[i][0]*2*z1*integral);
			line.append(space);
			line.append(haloDensity[i]);
			line.append(space);
			line.append(gradPotentialStar[i]);
			line.append(space);
			line.append(gradPotentialGas[i]);
			line.append(space);
			line.append(gradPotentialBaryon[i]);
			line.append(space);
			line.append(gradPotentialHalo[i]);
			line.append(space);
			line.append(gradPotentialTotal[i]);
			data.add(line.toString());
		}
		printToFile(data, pathname);
	}
	
	/**
	 * @return The total baryonic mass in this model galaxy
	 */
	public double getBaryonicMass() {
		double mass = 8*PI*maxDensity*z0*r0*r0;
		//This is currently just the stellar mass
		double[] parameters = getTruncations();
		double aS = parameters[0];
		double bS = parameters[1];
		double aG = parameters[2];
		double bG = parameters[3];
		double rFinalS = parameters[4];
		double rFinalG = parameters[5];
		double[] radius = new double[(int) ((max(
			max(rMax, rFinalS), rFinalG))/r0*50)];
		for (int i = 0; i < radius.length; i = i + 1) {
			radius[i] = ((double) i)*0.02*r0;
		}
		//First work out the parameters of the linear fit to the
		//interior of the gas density function
		//Assume radius.length >= 52
		//These are all surface mass densities
		double[] linearFit = getInnerDensities(50, radius);
		double m = linearFit[0];
		double c = linearFit[1];
		for (int i = 0; i < radius.length; i = i + 1) {
			double gasDensity = 0;
			double stellarDensity = 0;
			if (radius[i] <= rMax) {
				double starSurfaceDensity = 4*z0*maxDensity*pow(E, -1*radius[i]/r0);
				//The surface mass density of the stars at radius[i]
				if (i > 50) {
					double sigmaG = sigma0*pow(E, -1*radius[i]/rG);
					double v = URC(radius[i]);
					double k = v/radius[i]*sqrt(2*(1 + radius[i]/v*dURCdr(
						radius[i])));
					double sigmaSR = 1.67*sqrt(2*PI*G*starSurfaceDensity*z0) 
						+ sigmaSAsymptote;
					//The stellar radial velocity dispersion
					gasDensity = sigmaG*k/PI/G/Q - sigmaG/sigmaSR*starSurfaceDensity;
				} else {
					gasDensity = m*radius[i] + c;
				}
			} else {
				gasDensity = aG + bG*pow(E, radius[i]/r0);
				stellarDensity = aS + bS*pow(E, radius[i]/r0);
			}
			if (gasDensity > 0) {
				mass = mass + 2*PI*radius[i]*gasDensity*(radius[1] - radius[0]);
			}
			if (stellarDensity > 0) {
				mass = mass + 2*PI*radius[i]*stellarDensity*(radius[1] - radius[0]);
			}
		}
		return mass;
	}
	
	/**
	 * @return The total stellar mass of this galaxy
	 */
	public double getStellarMass() {
		return 8*PI*maxDensity*z0*r0*r0;
	}
	
	/**
	 * @return The total integrated star formation efficiency (HI) of this galaxy,
	 * i.e. the total star formation rate divided by the total HI mass as per the
	 * star formation law of Wong et al. (2016) and the molecular to atomic hydrogen
	 * prescription based on hydrostatic pressure.
	 */
	public double getSFEHI() {
		double kB = 1.38065*pow(10, -23);
		//Botzmann's constant
		double sfr = 0;
		//The net star formation rate of the galaxy
		double hIMass = 0;
		double[] parameters = getTruncations();
		double aS = parameters[0];
		double bS = parameters[1];
		double aG = parameters[2];
		double bG = parameters[3];
		double rFinalS = parameters[4];
		double rFinalG = parameters[5];
		double[] radius = new double[(int) ((max(
			max(rMax, rFinalS), rFinalG))/r0*200)];
		for (int i = 0; i < radius.length; i = i + 1) {
			radius[i] = ((double) i)*r0/200;
		}
		double dr = radius[1] - radius[0];
		double[] linearFit = getInnerDensities(200, radius);
		double m = linearFit[0];
		double c = linearFit[1];
		for (int i = 0; i < radius.length; i = i + 1) {
			double starDensity = 0;
			double gasDensity = 0;
			double r = radius[i];
			double sigmaG = sigma0*pow(E, -1*r/rG);
			if (r > rMax) {
				starDensity = aS + bS*pow(E, r/r0);
				gasDensity = aG + bG*pow(E, r/r0);
			} else {
				starDensity = 4*z0*maxDensity*pow(E, -r/r0);
				if (i < 200) {
					gasDensity = m*r + c;
				} else {
					double v = URC(radius[i]);
					double k = v/radius[i]*sqrt(2*(1 + radius[i]/v*dURCdr(
						radius[i])));
					double sigmaSR = 1.67*sqrt(2*PI*G*starDensity*z0) 
						+ sigmaSAsymptote;
					//The stellar radial velocity dispersion
					gasDensity = sigmaG*k/PI/G/Q - sigmaG/sigmaSR*starDensity;
				}
			}
			double sigmaSZ = sqrt(2*PI*G*starDensity*z0) + sigmaSAsymptote; 
			double pH = PI/2*G*gasDensity*(gasDensity + sigmaG/sigmaSZ*starDensity);
			//The hydrostatic pressure at this radius
			double rMol = pow(pH/1.7/pow(10, 10)/kB, 0.8);
			//The ratio of the H2 and HI surface mass densities at this radius (here
			//H2 can mean anything heavier than HI for some authors).
			double hIDensity = 0;
			double h2Density = 0;
			double sfrDensity = 0;
			if (rMol > 0) {
				hIDensity = gasDensity/(rMol + 1);
				//This is in kg/m^2
				h2Density = gasDensity/(1 + 1/rMol)/Msun*pow(Kpc/1000, 2);
				//This is in Msun/pc^2
				sfrDensity = 8.66*pow(10, -4)*h2Density*pow(10, 0.84*(atan(
					5*(log10(h2Density) - 2.5))/PI + 0.5));
				//This is in Msun/Kpc^2/yr
			}
			sfrDensity = sfrDensity*Msun/Kpc/Kpc/365.25/24/3600;
			//Now the units are kg/m^2/s
			sfr = sfr + 2*PI*sfrDensity*r*dr;
			hIMass = hIMass + 2*PI*hIDensity*r*dr;
		}
		return sfr/hIMass;
	}
	
	/**
	 * @return The galaxy's R band luminosity in solar units
	 */
	public double getRBandLuminosity() {
		double mR = -3.9 - 7.622*log10(vMax/1000);
		return pow(10, -0.4*(mR - 4.61));
	}
	
	/**
	 * @return The maximum rotational speed of any object in the galaxy
	 */
	public double getVMax() {return vMax;}
	
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
		double starDensity2 = 4*z0*maxDensity*pow(E, -1*rMax/r0);
		double starDensity1 = 4*z0*maxDensity*pow(E, -1*(rMax - 0.001*r0)/r0);
		double sigmaG2 = sigma0*pow(E, -1*rMax/rG);
		double sigmaG1 = sigma0*pow(E, -1*(rMax - 0.001*r0)/rG);
		double v2 = URC(rMax);
		double k2 = v2/rMax*sqrt(2*(1 + rMax/v2*dURCdr(rMax)));
		double sigmaSR2 = 1.67*sqrt(2*PI*G*starDensity2*z0) + sigmaSAsymptote;
		double gasDensity2 = sigmaG2*k2/PI/G/Q - sigmaG2/sigmaSR2*starDensity2;
		double v1 = URC(rMax - 0.001*r0);
		double k1 = v1/(rMax - 0.001*r0)*sqrt(2*(1 + (rMax - 0.001*r0)/v1*dURCdr(
			rMax - 0.001*r0)));
		double sigmaSR1 = 1.67*sqrt(2*PI*G*starDensity1*z0) + sigmaSAsymptote;
		double gasDensity1 = sigmaG1*k1/PI/G/Q - sigmaG1/sigmaSR1*starDensity1;
		double aS = starDensity2 - (starDensity2 - starDensity1)/0.001;
		double bS = pow(E, -rMax/r0)*(starDensity2 - starDensity1)/0.001;
		double aG = gasDensity2 - (gasDensity2 - gasDensity1)/0.001;
		double bG = pow(E, -rMax/r0)*(gasDensity2 - gasDensity1)/0.001;
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
	private double[] getInnerDensities(int i, double[] radius) {
		double sigmaG1 = sigma0*pow(E, -1*radius[i]/rG);
		double stellarDensity1 = 4*z0*maxDensity*pow(E, -1*radius[i]/r0);
		double v1 = URC(radius[i]);
		double k1 = v1/radius[i]*sqrt(2*(1 + radius[i]/v1*dURCdr(radius[i])));
		double sigmaSR1 = 1.67*sqrt(2*PI*G*stellarDensity1*z0) + sigmaSAsymptote;
		double gasDensity1 = sigmaG1*k1/PI/G/Q - sigmaG1/sigmaSR1*
			stellarDensity1;
		double sigmaG2 = sigma0*pow(E, -1*radius[i + 1]/rG);
		double stellarDensity2 = 4*z0*maxDensity*pow(E, -1*radius[i + 1]/r0);
		double v2 = URC(radius[i + 1]);
		double k2 = v2/radius[i + 1]*sqrt(2*(1 + radius[i + 1]/v2*dURCdr(
			radius[i + 1])));
		double sigmaSR2 = 1.67*sqrt(2*PI*G*stellarDensity2*z0) + sigmaSAsymptote;
		double gasDensity2 = sigmaG2*k2/PI/G/Q - sigmaG2/sigmaSR2*
			stellarDensity2;
		double m = (gasDensity2 - gasDensity1)/(radius[i + 1] - radius[i]);
		double c = gasDensity1 - radius[i]*m;
		return new double[] {m, c};
	}
	
	/**
	 * All the arrays will be re-written with their correct values for
	 * this model galaxy so it is irrelevant how one initializes them.
	 * The dimensions of the arrays are assumed to be properly related.
	 * @param radius The radius of each circular rim or spherical shell being
	 * used to approximate an exponential stellar disk
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
	 */
	private void processModel(double[] radius, double[] height,
		double[][] stellarDensity, double[][] gasDensity, double[] rotationCurve,
		double[] gradPotentialGas, double[] gradPotentialStar,
		double[] gradPotentialBaryon, double[] gradPotentialTotal, int lower,
		double aS, double bS, double aG, double bG) {
		for (int i = 0; i < radius.length; i = i + 1) {
			radius[i] = ((double) i)*r0/50;
		}
		for (int i = 0; i < height.length; i = i + 1) {
			height[i] = ((double) i)*z0*6/100;
		}
		//First get the parameters of the linear fit to the interior of 
		//the gas density function
		//Assume radius.length >= 52
		//These are all surface mass densities
		double[] linearFit = getInnerDensities(50, radius);
		double m = linearFit[0];
		double c = linearFit[1];
		for (int i = 0; i < radius.length; i = i + 1) {
			rotationCurve[i] = URC(radius[i]);
			double gasSurfaceDensity = 0;
			if (radius[i] <= rMax) {
				double starSurfaceDensity = 4*z0*maxDensity*pow(E, -1*radius[i]/r0);
				//The surface mass density of the stars at radius[i]
				for (int j = 0; j < height.length; j = j + 1) {
					stellarDensity[i][j] = 4*maxDensity*pow(E, -1*radius[i]/r0)*pow(
						E, height[j]/z0)/pow(pow(E, height[j]/z0) + 1, 2);
				}
				if (i > 50) {
					double sigmaG = sigma0*pow(E, -1*radius[i]/rG);
					double v = rotationCurve[i];
					double k = v/radius[i]*sqrt(2*(1 + radius[i]/v*dURCdr(
						radius[i])));
					double sigmaSR = 1.67*sqrt(2*PI*G*starSurfaceDensity*z0) 
						+ sigmaSAsymptote;
					//The stellar radial velocity dispersion
					gasSurfaceDensity = sigmaG*k/PI/G/Q - sigmaG/sigmaSR*
						starSurfaceDensity;
				} else {
					gasSurfaceDensity = m*radius[i] + c;
				}
			} else {
				for (int j = 0; j < height.length; j = j + 1) {
					stellarDensity[i][j] = (aS + bS*pow(E, radius[i]/r0))/z0*pow(
						E, height[j]/z0)/pow(pow(E, height[j]/z0) + 1, 2); 
				}
				gasSurfaceDensity = aG + bG*pow(E, radius[i]/r0);
			}
			if (gasSurfaceDensity < 0) {
				gasSurfaceDensity = 0;
		        //Gas is not necessarily found everywhere within a galaxy.
				//If negative gas density is required to maintain constant Q,
				//then it means there is really no gas there at all and the
				//two fluid Q is inappropriate.
			}
			double starSurfaceDensity = 4*z0*maxDensity*pow(E, -1*radius[i]/r0);
			double sigmaG = sigma0*pow(E, -1*radius[i]/rG);
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
				-1*radius[i]/r0), 0.5);
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
		//The process below is essentially a numerical integration.
		//I'm using a single set of loops to find the accelerations due to both
		//gas and stars.
		for (int i = 0; i < gradPotentialStar.length; i = i + 1) {
			double R = radius[i + lower];
			double dr = radius[1] - radius[0];
			double dtheta = PI/100;
			double dz = height[1] - height[0];
			double totalSumS = 0;
			double totalSumG = 0;
			for (int j = 0; j < radius.length; j = j + 1) {
				double r = radius[j];
				if (r != 0) {
					double angleSumS = 0;
					double angleSumG = 0;
					for (int k = 0; k < 100; k = k + 1) {
						double theta = ((double) k)*100/PI;
						double zSumS = 0;
						double zSumG = 0;
						for (int n = 0; n < height.length; n = n + 1) {
							double z = height[n];
							if (r != R || k != 0 || z != 0) {
								zSumS = zSumS + stellarDensity[j][n]/pow(
									z*z + r*r + R*R - 2*r*R*cos(theta), 1.5);
								zSumG = zSumG + gasDensity[j][n]/pow(
										z*z + r*r + R*R - 2*r*R*cos(theta), 1.5);
							}
						}
						angleSumS = angleSumS + zSumS*(R - r*cos(theta));
						angleSumG = angleSumG + zSumG*(R - r*cos(theta));
					}
					totalSumS = totalSumS + angleSumS*r;
					totalSumG = totalSumG + angleSumG*r;
				}
			}
			gradPotentialStar[i] = totalSumS*4*G*dz*dtheta*dr;
			gradPotentialGas[i] = totalSumG*4*G*dz*dtheta*dr;
		}
		for (int i = 0; i < gradPotentialBaryon.length; i = i + 1) {
			if (radius[i + lower] == 0) {
				gradPotentialTotal[i] = 0;
				gradPotentialBaryon[i] = 0;
			} else {
				gradPotentialTotal[i] = pow(rotationCurve[i + lower],
					2)/radius[i + lower];
				gradPotentialBaryon[i] = gradPotentialStar[i] + gradPotentialGas[i];
			}
		}
	}
	
	/**
	 * This method is called URC to keep consistency with the previous classes,
	 * but it is not an implementation of a URC.
	 * @param r The galocentric radius
	 * @return The velocity at a radius r for this model galaxy based
	 * on data used to create the URC of Karukes and Salucci.
	 */
	private double URC(double r) {return 2*vMax/PI*atan(r/rT);}
	
	/**
	 * This method is called dURCdr to keep consistency with the previous classes,
	 * but it is not an implementation of a URC.
	 * @param r The galactocentric radius
	 * @return The derivative of the rotation curve with respect to radius
	 */
	private double dURCdr(double r) {return 2*vMax/PI/rT/(1 + pow(r/rT, 2));}
	
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
	
	/**
	 * @param data A list of doubles to be appended one after the other
	 * to the file given by pathname
	 * @param pathname The complete path and filename to which the data
	 */
	private <E> void printToFile(ArrayList<E> data, String pathname) {
		try {
			File f = new File(pathname);
			FileWriter fw = new FileWriter(f, true);
			BufferedWriter bw = new BufferedWriter(fw);
			PrintWriter printer = new PrintWriter(bw);
			for (E line : data) {
				printer.println(line);
			}
			printer.close();
		} catch (IOException e) {
			System.out.println("Something terrible happened");
		}
	}
}