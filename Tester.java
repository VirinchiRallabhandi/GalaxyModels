import static java.lang.Math.*;
import java.util.ArrayList;
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.IOException;
/**
 * A class for running various simulations with the GalaxyModels and
 * collating/plotting results. The contents of this class are very
 * frequently changing and the user is required to edit the main method
 * to perform whichever task they are interested in.
 * All units are SI unless otherwise stated.
 * @author Virinchi Rallabhandi
 *
 */
public class Tester {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//oneGalaxy(75000, false, false, false, 
		//	"/Users/virinchirallabhandi/Desktop/Data/75w.txt");
		MDAR(0.25, 0.75);
		//testBTFR();
		//getSechIntegrals();
		//getGasMasses();
		//getGasMassToLight(
		//	"/Users/virinchirallabhandi/Desktop/Data/gasMassToLight.txt");
		//getSFE("/Users/virinchirallabhandi/Desktop/Data/SFE6.txt");
	}
	
	private static void MDAR(double start, double end) {
		long begin = System.currentTimeMillis();
		ArrayList<Double> gBar = new ArrayList<Double>();
		ArrayList<Double> gObs = new ArrayList<Double>();
		for (int i = 0; i < 21; i = i + 1) {
			System.out.println(i);
			double vFinal = 75000 + i*10000;
			GalaxyModel galaxy = new GalaxyModel4(vFinal, 
				"/Users/virinchirallabhandi/Desktop/Data/sechIntegrals.txt");
			double[][] accelerations = galaxy.calculateAccelerations(start, end);
			for (int j = 0; j < accelerations[0].length; j = j + 1) {
				gBar.add(accelerations[0][j]);
				gObs.add(accelerations[1][j]);
			}
		}
		printToFile(gBar, "/Users/virinchirallabhandi/Desktop/Data/gBar27.txt");
		printToFile(gObs, "/Users/virinchirallabhandi/Desktop/Data/gObs27.txt");
		long complete = System.currentTimeMillis();
		System.out.println("The program took "+(complete - begin)+" ms.");
	}
	
	/**
	 * @param v0 The v0 parameter of the URC if using GalalxyModel1 or GalaxyModel2
	 * Otherwise v0 can represent the velocity at rMax or simply vMax. This class
	 * adheres to very few OOP principles and hence the best way to find meaning
	 * of parameters is to read this class' code itself.
	 * @param sA Whether the galaxy is of Sa Hubble type
	 * @param sB Whether the galaxy is of Sb Hubble type
	 * @param sC Whether the galaxy is of Sc Hubble type
	 * The Hubble types are only relevant if using GalaxyModel1 or GalaxyModel2
	 * @param pathname The complete path and file to which the data for
	 * this model galaxy is to be stored to
	 */
	private static void oneGalaxy(double v0, boolean sA, boolean sB, boolean sC,
		String pathname) {
		long start = System.currentTimeMillis();
		GalaxyModel galaxy = new GalaxyModel13(0.75, v0,
			"/Users/virinchirallabhandi/Desktop/Data/sechIntegrals.txt");
		galaxy.recordData(pathname);
		long end = System.currentTimeMillis();
		System.out.println("The program took "+(end - start)+" ms.");
	}
	
	/**
	 * Writes the baryonic masses and vMax values to file so that the baryonic
	 * Tully Fisher relation may be tested
	 */
	private static void testBTFR() {
		ArrayList<Double> masses = new ArrayList<Double>();
		ArrayList<Double> speeds = new ArrayList<Double>();
		for (int i = 0; i < 21; i = i + 1) {
			GalaxyModel4 galaxy = new GalaxyModel4(75000 + i*10000,
				"/Users/virinchirallabhandi/Desktop/Data/sechIntegrals.txt");
			speeds.add(galaxy.getVMax());
			masses.add(galaxy.getBaryonicMass());
		}
		printToFile(masses, "/Users/virinchirallabhandi/Desktop/Data/mass3.txt");
		printToFile(speeds, "/Users/virinchirallabhandi/Desktop/Data/vMax3.txt");
	}
	
	/**
	 * Writes the estimated values of Integrate[Sech[x]^(2*p), {x, 0, Infinity}]
	 * to file so that they don't have to be recalculated every time. p is taken
	 * to be in the range 1 to 3000 (inclusive) in steps of 0.1
	 */
	private static void getSechIntegrals() {
		ArrayList<String> values = new ArrayList<String>();
		for (double p = 1; p <= 3000; p = p + 0.1) {
			String line = p+"   ";
			double sum = 0;
			double dx = 0.01;
			//Use Simpson' rule
			for (int i = 0; i < 1001; i = i + 1) {
				double x = ((double) i)*dx;
				if (i == 0 || i == 1000) {
					sum = sum + pow(2*pow(E, x)/(pow(E, 2*x) + 1), 2*p);
				} else {
					if (i % 2 == 1) {
						sum = sum + 4*pow(2*pow(E, x)/(pow(E, 2*x) + 1), 2*p);
					} else {
						sum = sum + 2*pow(2*pow(E, x)/(pow(E, 2*x) + 1), 2*p);
					}
				}
			}
			line = line+(sum*dx/3);
			values.add(line);
		}
		printToFile(values, 
			"/Users/virinchirallabhandi/Desktop/Data/sechIntegrals.txt");
	}
	
	/**
	 * Creates many galaxy models with varying Q values and vMax and gets
	 * their gas mass. To vary the Q values the Q field in the GalaxyModel
	 * class will need to be adjusted to remove the final keyword (this is
	 * very bad object oriented programming).
	 */
	private static void getGasMasses() {
		ArrayList<String> simulations = new ArrayList<String>();
		for (double q = 1; q <= 3.1; q = q + 0.1) {
			for (double vFinal = 75000; vFinal <= 255000; vFinal = vFinal + 5000) {
				//GalaxyModel4.Q = q;
				GalaxyModel4 galaxy = new GalaxyModel4(vFinal,
					"/Users/virinchirallabhandi/Desktop/Data/sechIntegrals.txt");
				double vMax = galaxy.getVMax()/1000;
				double gasMass = (galaxy.getBaryonicMass() 
					- galaxy.getStellarMass())/1.989/pow(10, 30);
				String line = Double.toString(gasMass)+"  "+Double.toString(vMax);
				line = line+"  "+Double.toString(q);
				simulations.add(line);
			}
		}
		printToFile(simulations,
			"/Users/virinchirallabhandi/Desktop/Data/gasMassData.txt");
	}
	
	/**
	 * Prints into pathname a table. The first column is the 
	 * log10(gasMassToLightRatio (solar units)). The second column is log10(vMax).
	 * The third column is the value of Q.
	 * Using this method requires removing the final keyword from the Q field 
	 * in GalaxyModel4 (which is very bad software engineering).
	 * @param pathname The file to which the data is to be output to
	 */
	private static void getGasMassToLight(String pathname) {
		ArrayList<String> simulations = new ArrayList<String>();
		for (double q = 1; q <= 4.1; q = q + 0.1) {
			for (double vFinal = 20000; vFinal <= 255000; vFinal = vFinal + 5000) {
				//GalaxyModel10.Q = q;
				//GalaxyModel4.Q = q;
				double vMax = 0;
				double gasMass = 0;
				double luminosity = 0;
				if (vFinal < 75000) {
					GalaxyModel10 galaxy = new GalaxyModel10(vFinal,
						"/Users/virinchirallabhandi/Desktop/Data/sechIntegrals.txt");
					vMax = galaxy.getVMax()/1000;
					gasMass = (galaxy.getBaryonicMass() 
						- galaxy.getStellarMass())/1.989/pow(10, 30);
					luminosity = galaxy.getRBandLuminosity();
				} else {
					GalaxyModel4 galaxy = new GalaxyModel4(vFinal,
						"/Users/virinchirallabhandi/Desktop/Data/sechIntegrals.txt");
					vMax = galaxy.getVMax()/1000;
					gasMass = (galaxy.getBaryonicMass() 
						- galaxy.getStellarMass())/1.989/pow(10, 30);
					luminosity = galaxy.getRBandLuminosity();
				}
				String line = Double.toString(log10(gasMass/luminosity))+"  ";
				line = line+Double.toString(log10(vMax))+"  "+Double.toString(q);
				simulations.add(line);
			}
		}
		printToFile(simulations, pathname);
	}
	
	/**
	 * Appends data for the variation in star formation efficiency (with respect
	 * to HI) to the file pathname. The first column will have log(vMax (km/s)) and 
	 * the second will have log(SFE (/yr))
	 * @param pathname The file to which the data is to be appended to
	 */
	private static void getSFE(String pathname) {
		ArrayList<String> data = new ArrayList<String>();
		for (double vFinal = 20000; vFinal <= 256000; vFinal = vFinal + 5000) {
			double vMax = 0;
			double sfe = 0;
			if (vFinal < 75000) {
				GalaxyModel11 galaxy = new GalaxyModel11(vFinal,
					"/Users/virinchirallabhandi/Desktop/Data/sechIntegrals.txt");
				vMax = galaxy.getVMax()/1000;
				sfe = galaxy.getSFEHI();
			} else {
				GalaxyModel11 galaxy = new GalaxyModel11(vFinal,
					"/Users/virinchirallabhandi/Desktop/Data/sechIntegrals.txt");
				vMax = galaxy.getVMax()/1000;
				sfe = galaxy.getSFEHI();
			}
			String line = Double.toString(log10(vMax))+"  "+Double.toString(
				log10(sfe*365.25*24*3600));
			data.add(line);
		}
		printToFile(data, pathname);
	}
	
	/**
	 * @param data A list of lines (of a type which can be printed e.g. String,
	 * Double) to be appended one after the other
	 * to the file given by pathname
	 * @param pathname The complete path and filename to which the data
	 */
	private static <E> void printToFile(ArrayList<E> data, String pathname) {
		try {
			File f = new File(pathname);
			FileWriter fw = new FileWriter(f, true);
			BufferedWriter bw = new BufferedWriter(fw);
			PrintWriter printer = new PrintWriter(bw);
			for (E number : data) {
				printer.println(number);
			}
			printer.close();
		} catch (IOException e) {
			System.out.println("Something terrible happened");
		}
	}
}