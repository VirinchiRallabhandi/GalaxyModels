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
		//oneGalaxy(255000, false, false, false, 
		//	"D:/University/ICRAR Studentship/Data/215a.txt");
		//MDAR(0.25, 0.75);
		//testBTFR();
		//getSechIntegrals();
		//getGasMasses();
		//strictMDAR(125000, true, false, false,
		//	"D:/University/ICRAR Studentship/Data/sechIntegrals.txt");
		strictMDARFullSet();
	}
	
	private static void MDAR(double start, double end) {
		long begin = System.currentTimeMillis();
		ArrayList<Double> gBar = new ArrayList<Double>();
		ArrayList<Double> gObs = new ArrayList<Double>();
		for (int i = 0; i < 21; i = i + 1) {
			System.out.println(i);
			double vFinal = 75000 + i*10000;
			GalaxyModel galaxy = new GalaxyModel4(vFinal, 
				"D:/University/ICRAR Studentship/Data/sechIntegrals.txt");
			double[][] accelerations = galaxy.calculateAccelerations(start, end);
			for (int j = 0; j < accelerations[0].length; j = j + 1) {
				gBar.add(accelerations[0][j]);
				gObs.add(accelerations[1][j]);
			}
		}
		printToFile(gBar, "D:/University/ICRAR Studentship/Data/gBar26.txt");
		printToFile(gObs, "D:/University/ICRAR Studentship/Data/gObs26.txt");
		long complete = System.currentTimeMillis();
		System.out.println("The program took "+(complete - begin)+" ms.");
	}
	
	/**
	 * @param v0 The v0 parameter of the URC if using GalalxyModel1 or GalaxyModel2
	 * Otherwise v0 represents the circular speed at rMax - the truncation radius
	 * of the galaxy's disk.
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
		GalaxyModel galaxy = new GalaxyModel4(v0,
			"D:/University/ICRAR Studentship/Data/sechIntegrals.txt");
		galaxy.recordData(pathname);
		long end = System.currentTimeMillis();
		System.out.println("The program took "+(end - start)+" ms.");
	}
	
	/**
	 * Finds the accelerations and rotation curves for a galaxy strictly
	 * adhering to MDAR.
	 * @param vFinal The final velocity used to initialise the galaxy
	 * @param f1 Whether or not to use the interpolating function, f1
	 * @param f2 Whether or not to use the interpolating function, f2
	 * @param f3 Whether or not to use the interpolating function, f3
	 * @param sechFile The location of the sech integrals file
	 */
	private static void strictMDAR(double vFinal, boolean f1, boolean f2,
		boolean f3, String sechFile) {
		long start = System.currentTimeMillis();
		GalaxyModel14 galaxy = new GalaxyModel14(vFinal, sechFile, f1, f2, f3);
		double[][] data = galaxy.calculateAccelerations();
		ArrayList<Double> radius = new ArrayList<Double>();
		ArrayList<Double> rotationCurve = new ArrayList<Double>();
		ArrayList<Double> gBar = new ArrayList<Double>();
		for (int i = 0; i < data[0].length; i = i + 1) {
			radius.add(data[0][i]);
			rotationCurve.add(data[1][i]);
			gBar.add(data[2][i]);
		}
		printToFile(radius,
			"D:/University/Physics/PHYS3003/Assignment/radius1.txt");
		printToFile(rotationCurve,
			"D:/University/Physics/PHYS3003/Assignment/v1.txt");
		printToFile(gBar, "D:/University/Physics/PHYS3003/Assignment/gBar1.txt");
		long end = System.currentTimeMillis();
		System.out.println("The program took "+(end - start)+" ms.");
	}
	
	private static void strictMDARFullSet() {
		String sechFile = "D:/University/ICRAR Studentship/Data/sechIntegrals.txt";
		long start = System.currentTimeMillis();
		for (int k = 0; k < 6; k = k + 1) {
			double vFinal = 125000 + 25000*k;
			GalaxyModel14 galaxy = null;
			for (int j = 0; j < 3; j = j + 1) {
				if (j == 0) {
					galaxy = new GalaxyModel14(vFinal, sechFile,
						true, false, false);
				}
				if (j == 1) {
					galaxy = new GalaxyModel14(vFinal, sechFile,
						false, true, false);
				}
				if (j == 2) {
					galaxy = new GalaxyModel14(vFinal, sechFile,
						false, false, true);
				}
				double[][] data = galaxy.calculateAccelerations();
				ArrayList<Double> radius = new ArrayList<Double>();
				ArrayList<Double> gBar = new ArrayList<Double>();
				for (int i = 0; i < data[0].length; i = i + 1) {
					radius.add(data[0][i]);
					gBar.add(data[2][i]);
				}
				printToFile(radius, "D:/University/Physics/PHYS3003/Assignment/"
					+ "radius1"+k+j+".txt");;
				printToFile(gBar, "D:/University/Physics/PHYS3003/Assignment/"
					+ "gBar1"+k+j+".txt");
			}
		}
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
		printToFile(masses, "D:/University/ICRAR Studentship/Data/mass2.txt");
		printToFile(speeds, "D:/University/ICRAR Studentship/Data/vMax2.txt");
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
			"D:/University/ICRAR Studentship/Data/sechIntegrals.txt");
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
					"D:/University/ICRAR Studentship/Data/sechIntegrals.txt");
				double vMax = galaxy.getVMax()/1000;
				double gasMass = (galaxy.getBaryonicMass() 
					- galaxy.getStellarMass())/1.989/pow(10, 30);
				String line = Double.toString(gasMass)+"  "+Double.toString(vMax);
				line = line+"  "+Double.toString(q);
				simulations.add(line);
			}
		}
		printToFile(simulations,
			"D:/University/ICRAR Studentship/Data/gasMassData.txt");
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