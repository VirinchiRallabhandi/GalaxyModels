import static java.lang.Math.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
/**
 * This class models a galaxy's mass distribution given its rotation velocity at
 * the edge of the disk. It does this by deriving the appropriate URC and a 2D gas
 * disk which maintains a constant Toomre stability parameter. There are no
 * stars in this model. This model is mainly for deciding whether a stable disk
 * along can reproduce MDAR. The dark matter distribution is therefore
 * deduced from the URC and baryonic matter distribution via Poisson's equation.
 * All units are SI unless otherwise stated.
 * @author Virinchi Rallabhandi
 *
 */
public class GalaxyModel8 implements GalaxyModel {
	/**
	 * The constant value of the Toomre stability parameter
	 */
	public static final double Q = 1;
	public static final double G = 6.67408*pow(10, -11);
	public static final double Kpc = 3.086*pow(10, 19);
	/**
	 * The maximum radius of the disk set by the requirement that
	 * the orbital period at the maximum radius is one 1 Gyr.
	 */
	private double rMax;
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
	 * The maximum rotational speed of any object in the galaxy
	 */
	private double vMax;
	/**
	 * The lambda parameter of the URC as given in the "note added in proof"
	 */
	private double lambda;
	
	/**
	 * @param vFinal The circular speed at the edge of the galaxy disk (the 
	 * disk will be truncated in a continuously differentiable manner
	 * past the edge)
	 */
	public GalaxyModel8(double vFinal) {
		//Most of the following scalings come from Wong et al. 2016,
		//Meurer et al. 2018 or van der Kruit and Freeman 2011.
		rMax = vFinal*pow(10, 9)*365.25*24*3600/2/PI;
		double mB = -10*log10(vFinal/1000) + 3;
		//Very rough estimate from expressions in Carrol and Ostlie and assuming
		//vMax = vFinal. These assumptions are only made so that a reasonable
		//initial guess is made for the parameter lambda.
		double luminosity = pow(10, -0.4*(mB - 5.31));
		lambda = luminosity/6/pow(10, 10)*pow(67.4/50.0, 2);
		//Using H0 = 67 km/s/Mpc, the latest value by the Planck collaboration
		//Iteratively work out what the correct value of lambda is in the URC
		for (int i = 0; i < 2000; i = i + 1) {
			double vFinalPredicted = URC(rMax);
			lambda = lambda + (vFinal - vFinalPredicted)/2000000;
			//v0 in the URC is a strictly increasing function of lambda so
			//this process should work given enough iterations.
		}
		vMax = 0;
		for (double r = 0; r < rMax; r = r + 0.001*rMax) {
			double v = URC(r);
			if (v > vMax) {
				vMax = v;
			}
		}
		double mR = -3.9 - 7.622*log10(vMax/1000);
		//The absolute magnitude of the galaxy in the R band
		double surfaceBrightness = pow(10, 5.3785 + 1.1757*log10(vMax/1000));
		//The effective surface brightness of the disk in Lsun/kpc^2
		double r0 = sqrt(pow(10, -0.4*(mR - 4.61))/5.647/PI/surfaceBrightness)*Kpc;
		//This would be the scale length of the exponential stellar disk if there
		//actually was one. I still need r0 to get rG, which is required.
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
		double aG = parameters[0];
		double bG = parameters[1];
		double rFinalG = parameters[2];
		double[] radius = new double[1200];;
		//The radius of each cylindrical rim or spherical shell at which
		//the summations will be performed
		double[] gasDensity = new double[radius.length];
		//The gas disk's volume mass density at each radius given constant Q
		double[] rotationCurve = new double[radius.length];
		//The URC velocity at each radius
		int lower = (int) (((double) radius.length)*start);
		int upper = (int) (((double) radius.length)*end - 1);
		double[] gradPotentialGas = new double[upper - lower + 1];
		//The acceleration due to the gas alone at each radius on the midplane
		double[] gradPotentialTotal = new double[gradPotentialGas.length];
		//The acceleration due to all the mass as determined by the URC
		processModel(radius, gasDensity, rotationCurve,
			gradPotentialGas, gradPotentialTotal, lower, aG, bG, rFinalG);
		double[][] results = new double[2][gradPotentialTotal.length];
		for (int i = 0; i < gradPotentialTotal.length; i = i + 1) {
			results[0][i] = log10(gradPotentialGas[i]);
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
		double aG = parameters[0];
		double bG = parameters[1];
		double rFinalG = parameters[2];
		double[] radius = new double[1200];
		//The radius of each cylindrical rim or spherical shell at which
		//the summations will be performed
		double[] gasDensity = new double[radius.length];
		//The gas disk's volume mass density at each radius given constant Q
		double[] rotationCurve = new double[radius.length];
		//The URC velocity at each radius
		double[] gradPotentialGas = new double[radius.length];
		//The acceleration due to the gas alone at each radius
		double[] gradPotentialTotal = new double[radius.length];
		//The acceleration due to all the mass as determined by the URC
		processModel(radius, gasDensity, rotationCurve, gradPotentialGas,
			gradPotentialTotal, 0, aG, bG, rFinalG);
		double[] gradPotentialHalo = new double[radius.length];
		for (int i = 0; i < gradPotentialHalo.length; i = i + 1) {
			gradPotentialHalo[i] = gradPotentialTotal[i] - gradPotentialGas[i];
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
			StringBuilder line = new StringBuilder();
			line.append(radius[i]);
			line.append(space);
			line.append(0);
			line.append(space);
			line.append(gasDensity[i]);
			line.append(space);
			line.append(haloDensity[i]);
			line.append(space);
			line.append(0);
			line.append(space);
			line.append(gradPotentialGas[i]);
			line.append(space);
			line.append(gradPotentialGas[i]);
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
		double mass = 0;
		double[] parameters = getTruncations();
		double aG = parameters[0];
		double bG = parameters[1];
		double rFinalG = parameters[2];
		double[] radius = new double[1200];
		double dr = max(rFinalG, rMax)/radius.length;
		for (int i = 0; i < radius.length; i = i + 1) {
			radius[i] = ((double) i)*dr;
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
			if (radius[i] <= rMax) {
				if (i > 50) {
					double sigmaG = sigma0*pow(E, -1*radius[i]/rG);
					double v = URC(radius[i]);
					double k = v/radius[i]*sqrt(2*(1 + radius[i]/v*dURCdr(
						radius[i])));
					gasDensity = sigmaG*k/PI/G/Q;
				} else {
					gasDensity = m*radius[i] + c;
				}
			} else {
				gasDensity = aG + bG*pow(E, radius[i]/0.1/rMax);
			}
			if (gasDensity > 0) {
				mass = mass + 2*PI*radius[i]*gasDensity*dr;
			}
		}
		return mass;
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
		//I fit a + be^(r/0.1/rMax) to the function, with b negative, to
		//the calculated values beyond rMax.
		double sigmaG2 = sigma0*pow(E, -1*rMax/rG);
		double v2 = URC(rMax);
		double k2 = v2/rMax*sqrt(2*(1 + rMax/v2*dURCdr(rMax)));
		double gasDensity2 = sigmaG2*k2/PI/G/Q;
		double sigmaG1 = sigma0*pow(E, -1*(0.999*rMax)/rG);
		double v1 = URC(0.999*rMax);
		double k1 = v1/(0.999*rMax)*sqrt(2*(1 + (0.999*rMax)/v1*dURCdr(0.999*rMax)));
		double gasDensity1 = sigmaG1*k1/PI/G/Q;
		double aG = gasDensity2 - (gasDensity2 - gasDensity1)/0.01;
		double bG = pow(E, -0.1)*(gasDensity2 - gasDensity1)/0.01;
		double rFinalG = 0.1*rMax*log(abs(aG/bG));
		return new double[] {aG, bG, rFinalG};
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
		double v1 = URC(radius[i]);
		double k1 = v1/radius[i]*sqrt(2*(1 + radius[i]/v1*dURCdr(radius[i])));
		double gasDensity1 = sigmaG1*k1/PI/G/Q;
		double sigmaG2 = sigma0*pow(E, -1*radius[i + 1]/rG);
		double v2 = URC(radius[i + 1]);
		double k2 = v2/radius[i + 1]*sqrt(2*(1 + radius[i + 1]/v2*dURCdr(
			radius[i + 1])));
		double gasDensity2 = sigmaG2*k2/PI/G/Q;
		double m = (gasDensity2 - gasDensity1)/(radius[i + 1] - radius[i]);
		double c = gasDensity1 - radius[i]*m;
		return new double[] {m, c};
	}
	
	/**
	 * All the arrays will be re-written with their correct values for
	 * this model galaxy so it is irrelevant how one initializes them.
	 * The dimensions of the arrays are assumed to be properly related.
	 * @param radius The radius of each circular rim or spherical shell at
	 * which calculations are performed
	 * @param gasDensity The surface mass density of the gas disk at each
	 * radius given constant Q
	 * @param rotationCurve The URC velocity at each radius
	 * @param gradPotentialGas The acceleration due to the gas alone at each radius
	 * @param gradPotentialTotal The total acceleration
	 * @param lower gradPotential*[i] should correspond to radius[i + lower]
	 */
	private void processModel(double[] radius, double[] gasDensity,
		double[] rotationCurve, double[] gradPotentialGas, double[] gradPotentialTotal,
		int lower, double aG, double bG, double rFinalG) {
		double dr = max(rMax, rFinalG)/radius.length;
		for (int i = 0; i < radius.length; i = i + 1) {
			radius[i] = ((double) i)*dr;
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
				if (i > 50) {
					double sigmaG = sigma0*pow(E, -1*radius[i]/rG);
					double v = rotationCurve[i];
					double k = v/radius[i]*sqrt(2*(1 + radius[i]/v*dURCdr(
						radius[i])));
					gasSurfaceDensity = sigmaG*k/PI/G/Q;
				} else {
					gasSurfaceDensity = m*radius[i] + c;
				}
			} else {
				gasSurfaceDensity = aG + bG*pow(E, radius[i]/0.1/rMax);
			}
			if (gasSurfaceDensity < 0) {
				gasSurfaceDensity = 0;
		        //Gas is not necessarily found everywhere within a galaxy.
				//If negative gas density is required to maintain constant Q,
				//then it means there is really no gas there at all and the
				//two fluid Q is inappropriate.
			}
			gasDensity[i] = gasSurfaceDensity;
		}
		//The process below is essentially a numerical integration.
		//I'm using a single set of loops to find the accelerations due to both
		//gas and stars.
		for (int i = 0; i < gradPotentialGas.length; i = i + 1) {
			double R = radius[i + lower];
			dr = radius[1] - radius[0];
			double dtheta = PI/100;
			double totalSum = 0;
			for (int j = 0; j < radius.length; j = j + 1) {
				double r = radius[j];
				if (r != 0) {
					double angleSum = 0;
					for (int k = 0; k < 100; k = k + 1) {
						if (r != R || k != 0) {
							double theta = PI*((double) k)/100;
							angleSum = angleSum + (R - r*cos(theta))/pow(
								r*r + R*R - 2*r*R*cos(theta), 1.5);
						}
					}
					if (j == 1 || j == radius.length - 1) {
						totalSum = totalSum + G*gasDensity[j]*r*angleSum;
					} else {
						totalSum = totalSum + 2*G*gasDensity[j]*r*angleSum;
					}
				}
			}
			gradPotentialGas[i] = totalSum*dr*dtheta;
		}
		for (int i = 0; i < gradPotentialTotal.length; i = i + 1) {
			if (radius[i + lower] > 0) {
				gradPotentialTotal[i] = pow(URC(radius[i + lower]),
					2)/radius[i + lower];
			} else {
				gradPotentialTotal[i] = 0;
			}
		}
	}
	
	/**
	 * @param r The galocentric radius
	 * @return The velocity at a radius r for this model galaxy based
	 * on the URC of Persic, Salucci and Stel.
	 */
	private double URC(double r) {
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
	 * 
	 * @param r The galocentric radius
	 * @return The derivative of the rotation curve with respect to radius
	 */
	private double dURCdr(double r) {
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
		return dv2dx*v0*v0/2/URC(r*Kpc)/rOpt*1000*1000/Kpc;
	}
	
	/**
	 * @param data A list of doubles to be appended one after the other
	 * to the file given by pathname
	 * @param pathname The complete path and filename to which the data
	 */
	private void printToFile(ArrayList<String> data, String pathname) {
		try {
			File f = new File(pathname);
			FileWriter fw = new FileWriter(f, true);
			BufferedWriter bw = new BufferedWriter(fw);
			PrintWriter printer = new PrintWriter(bw);
			for (String line : data) {
				printer.println(line);
			}
			printer.close();
		} catch (IOException e) {
			System.out.println("Something terrible happened");
		}
	}
}