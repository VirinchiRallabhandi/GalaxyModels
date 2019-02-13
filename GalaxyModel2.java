import static java.lang.Math.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
/**
 * This class models a galaxy's mass distribution given its maximum rotation
 * velocity and Hubble type. It does this by deriving the appropriate URC,
 * a stellar disk which is exponential radially and sech^2 vertically,
 * maintain a constant Toomre stability parameter across the two disks. The
 * dark matter distribution is therefore deduced from the URC and baryonic
 * matter distribution via Poisson's equation.
 * Any time a quantity is documented as surface ..., but refers to a 3D
 * quantity, it means volume density integrated over the z axis.
 * All units are SI unless otherwise stated.
 * @author Virinchi Rallabhandi
 *
 */
public class GalaxyModel2 implements GalaxyModel {
	/**
	 * The constant value of the Toomre stability parameter
	 */
	public static final double Q = 1.6;
	public static final double G = 6.67408*pow(10, -11);
	public static final double Kpc = 3.086*pow(10, 19);
	/**
	 * The leading coefficient of the URC rotation curve
	 */
	private double v0;
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
	 * The velocity dispersion of the gaseous disk
	 */
	public static final double sigmaG = 8000;
	/**
	 * The a/alpha parameter of the URC (Persic et al. call it a but I've
	 * re-named it alpha because who calls their parameters a and beta?)
	 */
	private double alpha;
	/**
	 * The beta parameter of the URC
	 */
	private double beta;
	
	/**
	 * @param vMax The maximum rotational speed of any object in the galaxy
	 * @param sA Whether the galaxy is of Sa type
	 * @param sB Whether the galaxy is of Sb type
	 * @param sC Whether the galaxy is of Sc type
	 */
	public GalaxyModel2(double vMax, boolean sA, boolean sB, boolean sC) {
		if ((sA && sB) || (sB && sC) || (sC && sA) || !(sA || sB || sC)) {
			throw new IllegalArgumentException("A galaxy must be of "
				+ "exactly one Hubble type.");
		}
		//Most of the following scalings come from Wong et al. 2016
		//or van der Kruit and Freeman 2011
		rMax = vMax*pow(10, 9)*365.25*24*3600/2/PI;
		double mR = -3.9 - 7.622*log10(vMax/1000);
		//The absolute magnitude of the galaxy in the R band
		double surfaceBrightness = pow(10, 5.3785 + 1.1757*log10(vMax/1000));
		//The effective surface brightness of the disk in Lsun/kpc^2
		r0 = sqrt(pow(10, -0.4*(mR - 4.61))/5.647/PI/surfaceBrightness)*Kpc;
		z0 = (0.45*vMax/100000 - 0.14)*Kpc;
		double stellarMass = pow(10, -0.4*(mR - 4.61))*pow(10,
			-1.578 - 0.0856*mR)*1.989*pow(10, 30);
		maxDensity = stellarMass/8/PI/z0/r0/r0;
		double mB = 0;
		//Absolute magnitude in the blue band by the appropriate
		//Tully-Fisher relation as per Carrol and Ostlie
		if (sA) {
			mB = -9.95*log10(vMax/1000) + 3.15;
		}
		if (sB) {
			mB = -10.2*log10(vMax/1000) + 2.71;
		}
		if (sC) {
			mB = -11.0*log10(vMax/1000) + 3.31;
		}
		double luminosity = pow(10, -0.4*(mB - 5.31));
		alpha = 1.5*pow((luminosity/6/pow(10, 10)*1.4*1.4), 0.2);
		beta = 0.72 + 0.44*log10(luminosity/6/pow(10, 10)*1.4*1.4);
		v0 = vMax;
		//Iteratively work out what the correct value of v0 should be so that
		//URC's predicted rotation curve actually has maximum velocity vMax.
		for (int i = 0; i < 20; i = i + 1) {
			double maxVFound = 0;
			for (double r = 0; r < rMax; r = r + 0.1*r0) {
				//The actual model extends slightly beyond rMax so that the
				//disks are not truncated abruptly, but we don't yet know
				//how far that will be so let r < rMax. The velocity hardly
				//changes that far out anyway.
				double v = URC(r);
				if (v > maxVFound) {
					maxVFound = v;
				}
			}
			v0 = v0 + 0.5*(vMax - maxVFound);
		}
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
		//To truncate the disk in a continuously differentiable manner,
		//I fit a + be^(r/r0) to the function, with b negative, to
		//the calculated values beyond rMax.
		//These densities refer to surface mass densities - i.e. volume
		//density integrated over the z axis.
		double starDensity2 = 4*z0*maxDensity*pow(E, -1*rMax/r0);
		double starDensity1 = 4*z0*maxDensity*pow(E, -1*(rMax - 0.02*r0)/r0);
		double v2 = URC(rMax);
		double k2 = v2/rMax*sqrt(2*(1 + rMax/v2*dURCdr(rMax)));
		double sigmaSR2 = 1.67*sqrt(2*PI*G*starDensity2*z0);
		double gasDensity2 = sigmaG*k2/PI/G/Q - sigmaG/sigmaSR2*starDensity2;
		double v1 = URC(rMax - 0.02*r0);
		double k1 = v1/(rMax - 0.02*r0)*sqrt(2*(1 + (rMax - 0.02*r0)/v1*dURCdr(
			rMax - 0.02*r0)));
		double sigmaSR1 = 1.67*sqrt(2*PI*G*starDensity1*z0);
		double gasDensity1 = sigmaG*k1/PI/G/Q - sigmaG/sigmaSR1*starDensity1;
		double aS = starDensity2 - (starDensity2 - starDensity1)/0.02;
		double bS = pow(E, -rMax/r0)*(starDensity2 - starDensity1)/0.02;
		double aG = gasDensity2 - (gasDensity2 - gasDensity1)/0.02;
		double bG = pow(E, -rMax/r0)*(gasDensity2 - gasDensity1)/0.02;
		double rFinalS = r0*log(abs(aS/bS));
		double rFinalG = r0*log(abs(aG/bG));
		double[] radius = new double[(int) ((max(
			max(rMax, rFinalS), rFinalG))/r0*50)];
		//double[] radius = new double[(int) (rMax/r0*100)];
		//The radius of each circular rim or spherical shell being used to
		//approximate an exponential stellar disk
		double[] height = new double[100];
		//The height above the z axis/midplane
		double[][] stellarDensity = new double[radius.length][height.length];
		//The exponential stellar disk's surface mass density at each radius
		double[] gasDensity = new double[radius.length];
		//The surface mass density of the gas disk at each radius given constant Q
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
		//To truncate the disk in a continuously differentiable manner,
		//I fit a + be^(r/r0) to the function, with b negative, to
		//the calculated values beyond rMax.
		//These densities refer to surface mass densities - i.e. volume
		//density integrated over the z axis.
		double starDensity2 = 4*z0*maxDensity*pow(E, -1*rMax/r0);
		double starDensity1 = 4*z0*maxDensity*pow(E, -1*(rMax - 0.02*r0)/r0);
		double v2 = URC(rMax);
		double k2 = v2/rMax*sqrt(2*(1 + rMax/v2*dURCdr(rMax)));
		double sigmaSR2 = 1.67*sqrt(2*PI*G*starDensity2*z0);
		double gasDensity2 = sigmaG*k2/PI/G/Q - sigmaG/sigmaSR2*starDensity2;
		double v1 = URC(rMax - 0.02*r0);
		double k1 = v1/(rMax - 0.02*r0)*sqrt(2*(1 + (rMax - 0.02*r0)/v1*dURCdr(
			rMax - 0.02*r0)));
		double sigmaSR1 = 1.67*sqrt(2*PI*G*starDensity1*z0);
		double gasDensity1 = sigmaG*k1/PI/G/Q - sigmaG/sigmaSR1*starDensity1;
		double aS = starDensity2 - (starDensity2 - starDensity1)/0.02;
		double bS = pow(E, -rMax/r0)*(starDensity2 - starDensity1)/0.02;
		double aG = gasDensity2 - (gasDensity2 - gasDensity1)/0.02;
		double bG = pow(E, -rMax/r0)*(gasDensity2 - gasDensity1)/0.02;
		double rFinalS = r0*log(abs(aS/bS));
		double rFinalG = r0*log(abs(aG/bG));
		double[] radius = new double[(int) ((max(
			max(rMax, rFinalS), rFinalG))/r0*50)];
		//double[] radius = new double[(int) (rMax/r0*100)];
		//The radius of each circular rim or spherical shell being used to
		//approximate an exponential stellar disk
		double[] height = new double[100];
		//The height above the z axis/midplane
		double[][] stellarDensity = new double[radius.length][height.length];
		//The exponential stellar disk's surface mass density at each radius
		double[] gasDensity = new double[radius.length];
		//The surface mass density of the gas disk at each radius given constant Q
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
			StringBuilder line = new StringBuilder();
			line.append(radius[i]);
			line.append(space);
			line.append(stellarDensity[i][0]*4*z0);
			line.append(space);
			line.append(gasDensity[i]);
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
		double[][] stellarDensity, double[] gasDensity, double[] rotationCurve,
		double[] gradPotentialGas, double[] gradPotentialStar,
		double[] gradPotentialBaryon, double[] gradPotentialTotal, int lower,
		double aS, double bS, double aG, double bG) {
		for (int i = 0; i < radius.length; i = i + 1) {
			radius[i] = ((double) i)*r0/50;
		}
		for (int i = 0; i < height.length; i = i + 1) {
			height[i] = ((double) i)*z0*6/100;
		}
		//First work out the parameters of the linear fit to the
		//interior of the gas density function
		//Assume radius.length >= 27
		//These are all surface mass densities
		double stellarDensity25 = 4*z0*maxDensity*pow(E, -1*radius[50]/r0);
		double v25 = URC(radius[50]);
		double k25 = v25/radius[50]*sqrt(2*(1 + radius[50]/v25*dURCdr(
			radius[50])));
		double sigmaSR25 = 1.67*sqrt(2*PI*G*stellarDensity25*z0);
		double gasDensity25 = sigmaG*k25/PI/G/Q - sigmaG/sigmaSR25*
			stellarDensity25;
		double stellarDensity26 = 4*z0*maxDensity*pow(E, -1*radius[51]/r0);
		double v26 = URC(radius[51]);
		double k26 = v26/radius[51]*sqrt(2*(1 + radius[51]/v26*dURCdr(
			radius[51])));
		double sigmaSR26 = 1.67*sqrt(2*PI*G*stellarDensity26*z0);
		double gasDensity26 = sigmaG*k26/PI/G/Q - sigmaG/sigmaSR26*
			stellarDensity26;
		double m = (gasDensity26 - gasDensity25)/(radius[51] - radius[50]);
		double c = gasDensity25 - radius[50]*m;
		for (int i = 0; i < radius.length; i = i + 1) {
			rotationCurve[i] = URC(radius[i]);
			if (radius[i] <= rMax) {
				double starSurfaceDensity = 4*z0*maxDensity*pow(E, -1*radius[i]/r0);
				//The surface mass density of the stars at radius[i]
				for (int j = 0; j < height.length; j = j + 1) {
					stellarDensity[i][j] = 4*maxDensity*pow(E, -1*radius[i]/r0)*pow(
						E, height[j]/z0)/pow(pow(E, height[j]/z0) + 1, 2);
				}
				if (i > 50) {
					double v = URC(radius[i]);
					double k = v/radius[i]*sqrt(2*(1 + radius[i]/v*dURCdr(
						radius[i])));
					double sigmaSR = 1.67*sqrt(2*PI*G*starSurfaceDensity*z0);
					//The stellar radial velocity dispersion
					gasDensity[i] = sigmaG*k/PI/G/Q - sigmaG/sigmaSR*
						starSurfaceDensity;
				} else {
					gasDensity[i] = m*radius[i] + c;
				}
			} else {
				for (int j = 0; j < height.length; j = j + 1) {
					stellarDensity[i][j] = (aS + bS*pow(E, radius[i]/r0))/z0*pow(
						E, height[j]/z0)/pow(pow(E, height[j]/z0) + 1, 2); 
				}
				gasDensity[i] = aG + bG*pow(E, radius[i]/r0);
			}
			if (gasDensity[i] < 0) {
				gasDensity[i] = 0;
		        //Gas is not necessarily found everywhere within a galaxy.
				//If negative gas density is required to maintain constant Q,
				//then it means there is really no gas there at all and the
				//two fluid Q is inappropriate.
			}
			for (int j = 0; j < height.length; j = j + 1) {
				if (stellarDensity[i][j] < 0) {
					stellarDensity[i][j] = 0;
					//This can happen if radius[i] is beyond the truncation
					//radius of the stellar disk.
				}
			}
		}
		//The process below is essentially a numerical integration.
		//Begin with the gas component.
		for (int i = 0; i < gradPotentialGas.length; i = i + 1) {
			double R = radius[i + lower];
			double dr = radius[1] - radius[0];
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
			//if (gradPotentialGas[i] < 0) {
			//	gradPotentialGas[i] = 0;
			//}
		}
		//Get the acceleration due to the stellar disk alone
		for (int i = 0; i < gradPotentialStar.length; i = i + 1) {
			double R = radius[i + lower];
			double dr = radius[1] - radius[0];
			double dtheta = PI/100;
			double dz = height[1] - height[0];
			double totalSum = 0;
			for (int j = 0; j < radius.length; j = j + 1) {
				double r = radius[j];
				if (r != 0) {
					double angleSum = 0;
					for (int k = 0; k < 100; k = k + 1) {
						double theta = ((double) k)*100/PI;
						double zSum = 0;
						for (int n = 0; n < height.length; n = n + 1) {
							double z = height[n];
							if (r != R || k != 0 || z != 0) {
								zSum = zSum + stellarDensity[j][n]/pow(
									z*z + r*r + R*R - 2*r*R*cos(theta), 1.5);
							}
						}
						angleSum = angleSum + zSum*(R - r*cos(theta));
					}
					totalSum = totalSum + angleSum*r;
				}
			}
			gradPotentialStar[i] = totalSum*4*G*dz*dtheta*dr;
		}
		for (int i = 0; i < gradPotentialBaryon.length; i = i + 1) {
			if (radius[i + lower] == 0) {
				gradPotentialTotal[i] = 0;
				gradPotentialBaryon[i] = 0;
			} else {
				gradPotentialTotal[i] = pow(rotationCurve[i + lower],
					2)/radius[i + lower];
				gradPotentialBaryon[i] = gradPotentialGas[i] + gradPotentialStar[i];
			}
		}
	}
	
	/**
	 * @param r The galocentric radius
	 * @return The velocity at a radius r for this model galaxy based
	 * on the URC of Persic and Salucci.
	 */
	private double URC(double r) {
		double x = r/3.2/r0;
	    return sqrt(v0*v0*beta*1.97*pow(x, 1.22)/pow(x*x + 0.6084, 1.43) + v0*v0*(
	        1 - beta)*(1 + alpha*alpha)*x*x/(x*x + alpha*alpha));
	}
	
	private double dURCdr(double r) {
		double x = r/3.2/r0;
		double dv2 = v0*v0*beta*(2.4*pow(x, 0.22)*pow(x*x + 0.6084,
			1.43) - 5.63*pow(x, 2.22)*pow(x*x + 0.6084, 0.43))/pow(
			x*x + 0.6084, 2.86) + 2*v0*v0*(1 - beta)*(
			1 + alpha*alpha)*alpha*alpha*x/pow(x*x + alpha*alpha, 2);
		return dv2/2/3.2/r0/URC(r);
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
