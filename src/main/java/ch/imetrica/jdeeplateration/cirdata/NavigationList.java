package ch.imetrica.jdeeplateration.cirdata;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.math.plot.Plot2DPanel;

import ch.imetrica.jdeeplateration.anchors.Anchors;
import ch.imetrica.jdeeplateration.matrix.Matrix;
import ch.imetrica.jdeeplateration.mstat.GradDescentResult;
import ch.imetrica.jdeeplateration.mstat.Mstat;
import de.micromata.opengis.kml.v_2_2_0.AltitudeMode;
import de.micromata.opengis.kml.v_2_2_0.Boundary;
import de.micromata.opengis.kml.v_2_2_0.Coordinate;
import de.micromata.opengis.kml.v_2_2_0.Document;
import de.micromata.opengis.kml.v_2_2_0.Feature;
import de.micromata.opengis.kml.v_2_2_0.Geometry;
import de.micromata.opengis.kml.v_2_2_0.Kml;
import de.micromata.opengis.kml.v_2_2_0.KmlFactory;
import de.micromata.opengis.kml.v_2_2_0.LineString;
import de.micromata.opengis.kml.v_2_2_0.LinearRing;
import de.micromata.opengis.kml.v_2_2_0.Placemark;
import de.micromata.opengis.kml.v_2_2_0.Point;
import de.micromata.opengis.kml.v_2_2_0.Polygon;
import de.micromata.opengis.kml.v_2_2_0.Style;


public class NavigationList {

	
	public final static double C = 299792458; 
	public final static double drift_bound = 3e-8; 
	double[] Frequencies; 
	double[] localOrigin = null;

	double CenterFreq; 
	double Bandwidth;
	
	long creationTimestamp;
	
	ArrayList<NavigationChannelList> navigationList = null;
	ArrayList<NavigationChannelList> filteredNavigationList = null;
	ArrayList<TimeDiffOnArrival> tdoaFrequency0;
	ArrayList<double[]> simulationCoordinates;
	
	private double[] noiseModel;
	private int maxOutlierThreshold = 55; 
	
	Matrix bounds_in;
	private double sourceLatitude = 0;
	private double sourceLongitude = 0;
	private boolean plotTDOAs = false;
	
	/* Parameters for Stochastic Gradient Solution and Simulations */
	
    private int n_trial = 10000; 
    private double alpha = 0.0001; 
    private double time_threshold = 50000000;
	private double noisePercent = 0.00;
	private double navigationSpeedMetersPerSec = 50.0;
	private int numberEstimatesAvg = 6;
	

	/**
	* Initiates a navigation list which is a dynamic array 
	* of Navigation channels read, parsed, and stored from a 
	* navigation log file
	* @see         NavigationChannelList
	*/
	public NavigationList() {
		navigationList = new ArrayList<NavigationChannelList>();
	}
	
	/**
	* Toggles plotting the filtered TDOA computation results
	* @param  plotTDOA if true will plot TDOAs
	* @see         NavigationChannelList
	*/
	public void setPlotTDOAs(boolean plotTDOA) {
		this.plotTDOAs = plotTDOA; 
	}
	
	public void setKnownSource(double lat, double longi) {
		this.sourceLatitude = lat; 
		this.sourceLongitude = longi; 
	}
	
	/**
	* Sets the number of Trials used in the stochastic gradient solution. This number should 
	* usually be at least 1000 to get reasonable estimates of the source location. 
	* Each trial generates a solution be starting a random
	* coordinate inside the bounds of a region where the source is thought to be located. This 
	* region is a square box about 200 km outside each coordinate in the navigation list. 
	* The default value is 5000.
	* 
	* @param  n_trials Number of trials used in the stochastic gradient solution
	* @throws Exception must be greater than 100
	*/
	
	public void setNumberTrials(int n_trials) throws Exception {
		if(n_trials < 100) { 
			throw new Exception("Number of Trials must be greater than 100");
		}
		this.n_trial = n_trials;
	}
	
	/**
	* Sets the learning rate of the stochastic gradient descent method. This number should 
	* usually be between .001 and .0001, and determines the step size in the gradients. 
	* For very large geographic areas, a larger step size of .001 should most likely be used.
	* Default value is set at .0001
	* 
	* @param  alpha The step-size in the stochastic gradient descent 
	* @throws Exception must be a positive number.
	*/
	public void setAlpha(double alpha) throws Exception {
		if(alpha < 0) { 
			throw new Exception("Alpha must be a positive number");
		}
		this.alpha = alpha;
	}
	
	/**
	* Sets the noise percentage used in computing the noisy TDOAs or FDOAs in the 
	* simulated multilateration environment. The noise is a number between 0 and 1 (inclusive)
	* that adds a random Gaussian number with standard deviation a percentage of
	* the distance to the given source. Thus the further away from the source, the more 
	* noisy the simulated TDOA/FDOA measurement. 
	* 
	* @param  noise Percentage of Gaussian noise relative to the distance to source
	* added to the simulated TDOA/FDOA computation
	* @throws Exception must be between 0 and 1 (inclusive)
	*/
	public void setNoisePercentage(double noise) throws Exception {
		if(noise < 0 || noise > 1) { 
			throw new Exception("Noise percentage must be between 0 and 1");
		}
		this.noisePercent = noise;
	}
	
	/**
	* Sets the constant speed of the simulated navigation for computing the FDOAs in the 
	* simulated multilateration environment. This number should be greater than 0. 
	* 
	* @param  navigationSpeedMetersPerSec Speed in meters per second
	* added to the simulated TDOA/FDOA computation
	* @throws Exception must be greater than 0
	*/
	public void setSimulatedNavigationSpeed(double navigationSpeedMetersPerSec) throws Exception {
		if(navigationSpeedMetersPerSec <= 0) { 
			throw new Exception("Speed must be greater than 0");
		}
		this.navigationSpeedMetersPerSec = navigationSpeedMetersPerSec;
	}
	
	/**
	* Extracts from the filtered navigation list the time-difference on arrival
	* information at each navigation point. Assumes a frameLength of 65536 milliseconds.
	* 
	* @param  freq Filters on this desired frequency channel
	* @throws Exception on filteredNavigationList must not be empty
	* @see    NavigationChannelList
	*/
	public void computeTDOAfromNavigationList00(int freq) throws Exception {
        
	       if(filteredNavigationList == null) {    
	           throw new Exception("Apply filtering to navigation list first");   
	       }
	        
	        double refTime = -1.0;
	        double tdoa = 0.0;
	        double frameLength = .065536;
	        
	        ArrayList<Double> times = new ArrayList<Double>();
	        
	        tdoaFrequency0 = new ArrayList<TimeDiffOnArrival>();
	        
	        System.out.println("tdoa with respect to first point (reference)----------------------------------------");
	        for (NavigationChannelList navList : filteredNavigationList) {

	            double time = navList.getTimeStampAtFreq(freq);

	            if(refTime < 0 && time > 0) { 
	                refTime = time; 
	            
	                tdoa = 0.0;
	                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa, 0);                
	                
	                times.add(time);
	                tdoaFrequency0.add(td);
	            }
	            else if(time > 0) {
	                
	                double timeStamp = (time - refTime) % frameLength;
	                if (timeStamp < frameLength / 2) {
	                    tdoa = timeStamp;
	                }
	                else {
	                    tdoa = - (frameLength - timeStamp);
	                }
	                
	                System.out.println(tdoa);
	                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa, time - refTime);                
	                
	                times.add(time);
	                tdoaFrequency0.add(td);
	            }            
	        }
	        System.out.println("----------------------------------------");
	        System.out.println("----------------------------------------");
	        
	        double[] x = new double[tdoaFrequency0.size()];
	        double[] y = new double[tdoaFrequency0.size()];
	        
	        for(int i = 0; i < x.length; i++) {
	        	y[i] = tdoaFrequency0.get(i).tdoa;
	        	x[i] = times.get(i).doubleValue(); 
	        	x[i] = (double)i; 
	        }
	        
	        Plot2DPanel plot = new Plot2DPanel();
	        		 
	         // add a line plot to the PlotPanel
	        plot.addLinePlot("TDOA plot", x, y);
	        JFrame frame = new JFrame("Plot of TDOAs");
	        frame.setSize(600, 500);
	        frame.setContentPane(plot);
	        frame.setVisible(true);
	          
	}
	
	/**
	* Extracts from the filtered navigation list the time-difference on arrival
	* information at each navigation point at the given the frequency channel index. 
	* Assumes a frameLength of 65536 milliseconds. It also applies difference 
	* filtering where outliers are discarded by taking a difference
	* in time. If the absolute difference is greater than the parameter thresh than the 
	* observation is thrown out.  
	* 
	* @param  freq Filters on this desired frequency channel
	* @param  thresh Given threshhold on which to throw out outliers
	* @throws Exception on filteredNavigationList must not be empty
	* @see    NavigationChannelList
	*/
	
	public void filterTDOAfromNavigationList(int freq, double thresh) throws Exception {
        
	       if(filteredNavigationList == null) {    
	           throw new Exception("Apply filtering to navigation list first");   
	       }
	        
	        double refTime = -1.0;
	        double tdoa = 0.0;
	        double frameLength = .065536;
	        
	        int outlierCount = 0;
	        
	        ArrayList<Double> times = new ArrayList<Double>();
	        
	        tdoaFrequency0 = new ArrayList<TimeDiffOnArrival>();
	        
	        System.out.println("tdoa with respect to first point (reference)----------------------------------------");
	        
	        for (NavigationChannelList navList : filteredNavigationList) {

	            double time = navList.getTimeStampAtFreq(freq);

	            if(refTime < 0 && time > 0) { 
	                refTime = time;  //set reference time for first measurement
	            
	                tdoa = 0.0;
	                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa, 0);                
	                
	                times.add(time);
	                tdoaFrequency0.add(td);
	                
	            }
	            else if(time > 0) {
	                
	                double timeStamp = (time - refTime) % frameLength;
	                if (timeStamp < frameLength / 2) {
	                    tdoa = timeStamp;
	                }
	                else {
	                    tdoa = - (frameLength - timeStamp);
	                }
	                
	                
	                
	                if(Math.abs(tdoa - tdoaFrequency0.get(tdoaFrequency0.size()-1).getTDOA()) < thresh) {
	                	
	                	System.out.println(tdoaFrequency0.size() + " " + tdoa + " and diff " + (tdoa - tdoaFrequency0.get(tdoaFrequency0.size()-1).getTDOA()));
		                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa, time - refTime);                
		                times.add(time);
		                tdoaFrequency0.add(td);
		                outlierCount = 0;
	                }
	                else if(outlierCount > maxOutlierThreshold) {
	                	
	                	System.out.println(tdoaFrequency0.size() + " " + tdoa + " and diff " + (tdoa - tdoaFrequency0.get(tdoaFrequency0.size()-1).getTDOA()));
		                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa, time - refTime);                
		                times.add(time);
		                tdoaFrequency0.add(td);
		                outlierCount = 0;
	                }
	                else {
	                	outlierCount++;
	                }
	            }            
	        }
	        System.out.println("----------------------------------------");
	        System.out.println("----------------------------------------");
	        
	        double[] x = new double[tdoaFrequency0.size()];
	        double[] y = new double[tdoaFrequency0.size()];
	        
	        for(int i = 0; i < x.length; i++) {
	        	y[i] = tdoaFrequency0.get(i).tdoa;
	        	x[i] = times.get(i).doubleValue(); 
	        	x[i] = (double)i; 
	        }
	        
	        Plot2DPanel plot = new Plot2DPanel();
	        plot.addLinePlot("TDOA plot", x, y);
	        JFrame frame = new JFrame("Plot of TDOAs");
	        frame.setSize(900, 700);
	        frame.setContentPane(plot);
	        frame.setVisible(true);
	          
	}
	
	/**
	* Extracts from the filtered navigation list the valid doppler measurements that 
	* satisfy the peak list conditions. It modifies the original peak list by dropping 
	* all peaks that have a power that is a threshhold smaller than the strongest peak. 
	* This rejects correlation peaks artifacts in the time domain that were created
	* by a hard filtering of the transmission channel in the frequency domain. It also only 
	* keeps the channels such that the strongest remaining peak is also the first one arriving.
	* This function is used in filterFDOAfromNavigationList
	* 
	* @param  freq Filters on this desired frequency channel index
	* @param  threshold Given threshhold on which to throw out the peaks 
	*         this much smaller than the strongest
	* @throws Exception on filteredNavigationList must not be empty
	* @see    NavigationChannelList
	* @see    filterFDOAfromNavigationList
	*/
	
	public void validateDopplerMeasurement(int freq, double threshold) throws Exception {
		
		if(filteredNavigationList == null) {    
	           throw new Exception("Apply filtering to navigation list first");   
	    }
	       
	    int i = 0;
	    System.out.println("Navigation readings before fdoa filtering: " + filteredNavigationList.size());	
		
	    while(i < filteredNavigationList.size()) {
	    	
	    	if(!filteredNavigationList.get(i).validateDopplerMeasurementOnChannel(freq, threshold)) {
	    		filteredNavigationList.remove(i);
	    	}
	    	else {
	    		i++;
	    	}	
	    }
	    System.out.println("Navigation readings after fdoa filtering: " + filteredNavigationList.size());
	}
	
	/**
	* Extracts from the filtered navigation list the valid doppler measurements that 
	* satisfy the peak list conditions from validateDopplerMeasurement. It also filters
	* on differences and on observations where the navigation is stationary and yet the fdoa
	* measurment is nonzero. Lastly, all absolute values of fdoa greater that the threshhold parameter
	* are discarded as well.   
	* 
	* @param  freq Filters on this desired frequency channel index
	* @param  threshold Given threshhold used to throw out first-order difference outliers
	* @param  powerThreshold threshhold used in the validateDopplerMeasurement   
	* @throws Exception on filteredNavigationList must not be empty
	* @see    NavigationChannelList
	* @see    validateDopplerMeasurement
	*/
	
	public void filterFDOAfromNavigationList(int freq, double threshold, double powerThreshold) throws Exception {
		
	       if(filteredNavigationList == null) {    
	           throw new Exception("Apply filtering to navigation list first");   
	       }
	       
	       this.validateDopplerMeasurement(freq, powerThreshold);
	       
	              
	       int i = 0;
	       System.out.println("Navigation readings before fdoa filtering: " + filteredNavigationList.size());	       
	       while(i < filteredNavigationList.size()) {
	            
	    	    filteredNavigationList.get(i).setFrequencyError(freq);
	    	    	    	    
	    	    //System.out.println(navList.getFDOA() + " " + navList.isStationary(.5));
	            if(filteredNavigationList.get(i).getFDOA() == 0 && !filteredNavigationList.get(i).isStationary(.5)) {
	            	filteredNavigationList.remove(i);
	            }
	            else if(Math.abs(filteredNavigationList.get(i).getFDOA()) > threshold) {
	            	filteredNavigationList.remove(i);
	            }
	            else if(i > 2 && Math.abs(filteredNavigationList.get(i).getFDOA() - filteredNavigationList.get(i-1).getFDOA()) > 1.0) {
	            	filteredNavigationList.remove(i);
	            }
	            else {i++;}
	       }	
	       System.out.println("Navigation readings after fdoa filtering: " + filteredNavigationList.size());
	       
	       double[] y = new double[filteredNavigationList.size()];
	       double[] x = new double[filteredNavigationList.size()];
	       
	       for (i = 0; i < filteredNavigationList.size(); i++) {
	    	   //System.out.print(filteredNavigationList.get(i).getFDOA()); 
	    	   //filteredNavigationList.get(i).printVelocity();
	    	   y[i] = filteredNavigationList.get(i).getFDOA();
	    	   x[i] = i;
	       }
	       
	        Plot2DPanel plot = new Plot2DPanel();
	        JFrame frame = new JFrame("Plot of FDOAs");
	        plot.addLinePlot("TDOA plot", x, y);	        
	        frame.setSize(900, 700);
	        frame.setContentPane(plot);
	        frame.setVisible(true);
	}
	

	
	/**
	 * Once the navigation log file has been read, this function will 
	 * compute the velocity given the recorded latitude and longitudes 
	 * coordinates along the their local timestamp measurements. 
	 * <p>
	 * The resulting velocities are stored in the navigation list
	 * @see         NavigationChannelList
	 */
	
	public void computeVelocity() {
		
		navigationList.get(0).setVelocityToZero();
		for(int i = 1; i < navigationList.size(); i++) {
		
			NavigationChannelList navList = navigationList.get(i-1);
			long prevTime = navList.getDaySeconds();
			double[] prevlocalECEF = navList.getLocalECEF();
			
			navigationList.get(i).computeVelocityECEF(prevTime, prevlocalECEF);
		}
	}
	
	
	/**
	 * Filters the navigation list on a specific frequency channel 
	 * index at a certain threshhold (in dB).  
	 * <p>
	 * The measurement will be discarded of the final navigation list 
	 * if either the frequency index does not exist or if the RSSI measurement
	 * is less than the given threshhold. 
	 *
	 * 
	 * @param  threshhold  the RSSI threshhold. Values less than this will be
	 * thrown out of the navigation list
	 * @param  freq desired frequency index
	 */
	
	public void filterOnRSSI(double threshhold, int freq) {
		
		int i = 0;
		System.out.println("Navigation readings before RSSI filtering: " + navigationList.size());
		while(i < navigationList.size()) { 
			
			if(navigationList.get(i).getRSSI(freq) < threshhold) {
				navigationList.remove(i);
			}
			else {i++;}
		}
		System.out.println("Navigation readings after RSSI filtering: " + navigationList.size());
	}
	
	
	/**
	 * Filters the out all navigation readings that have duplicate location readings. 
	 * Creates the filteredNavigationList which is used for the rest of the filtering 
	 * procedures.
	 *
	 */
    public void filterUnique() {
		
    	filteredNavigationList = new ArrayList<NavigationChannelList>();
    	
    	double prevLat = 0; double prevLong = 0;
    	
		prevLat = navigationList.get(0).getLatitude(); 
		prevLong = navigationList.get(0).getLongitude();
    	
		filteredNavigationList.add(navigationList.get(0));
		
		for(int i = 1; i < navigationList.size(); i++) {
			
			if(prevLat != navigationList.get(i).getLatitude() || 
					prevLong != navigationList.get(i).getLongitude()) {
				
				filteredNavigationList.add(navigationList.get(i));
				prevLat = navigationList.get(i).getLatitude(); 
				prevLong = navigationList.get(i).getLongitude();
			}
		}		
	}
	
    /**
	 * Filters the out all navigation readings that have duplicate location readings
	 * and on a randomized basis. A filtering procedure used for experiementing with 
	 * navigation location importance. Creates the filteredNavigationList which is used 
	 * for the rest of the filtering procedures.
	 *
	 * @param seed Random number generator seed
	 * @param skip Number of navigation entries to skip
	 */

    public void filterRandomlyUnique(int seed, int skip) {
		
    	Random random = new Random(seed);
    	filteredNavigationList = new ArrayList<NavigationChannelList>();
    	
    	double prevLat = 0; double prevLong = 0;
    	
		prevLat = navigationList.get(0).getLatitude(); 
		prevLong = navigationList.get(0).getLongitude();
    	
		filteredNavigationList.add(navigationList.get(0));
		
		int i = random.nextInt(skip);
		while(i < navigationList.size()) {
			
			if(prevLat != navigationList.get(i).getLatitude() || 
					prevLong != navigationList.get(i).getLongitude()) {
				
				filteredNavigationList.add(navigationList.get(i));
				prevLat = navigationList.get(i).getLatitude(); 
				prevLong = navigationList.get(i).getLongitude();
			}
			i = i + random.nextInt(skip);
		}				
	}    
    
	/**
	 * Once the navigation log file has been read, this function will 
	 * compute the velocity given the recorded latitude and longitudes 
	 * coordinates along the their local timestamp measurements in the 
	 * ECEF coordinate system in meters.
	 * <p>
	 * The resulting velocities are stored in the navigation list 
	 * @see         computeVelocity
	 */
	
    
    public void computeVelocityECEF() {
    	
    	localOrigin = NavigationChannelList.GeodeticToECEF(navigationList.get(0).getLatitude(),
    			navigationList.get(0).getLongitude(), 0);  	
    	for(NavigationChannelList nav : navigationList) {    		
    		nav.GeodeticToLocal(localOrigin);
    	}
    	computeVelocity();   	
    }
    
	/**
	 * Once the filtered navigation list has been computed, this function will compute 
	 * the bounds in which the estimated source might be located. The bounds are found by 
	 * taking the most extreme north/south/east/west coordinates and adding an additional 
	 * 200 kilometers in each direction. These bounds are then used in the stochastic gradient
	 * method to estimate the source for either the TDOA or FDOA multilateration.  
	 * @see estimateSourceSolutionFDOA
	 * @see estimateSourceSolutionTDOA
	 * @throws Exception if matrix dimensions are wrong
	 */
    public void createEstimationBounds() throws Exception {
    	
	   if(filteredNavigationList == null) {    
	           throw new Exception("Apply filtering to navigation list first");   
	   }
	   
	   double minLong = 400;
	   double maxLong = -1.0;
	   double minLat = 100.0;
	   double maxLat = -100.0;
	   
	   
	   localOrigin = NavigationChannelList.GeodeticToECEF(filteredNavigationList.get(0).getLatitude(), 
			   filteredNavigationList.get(0).getLongitude(), 0);
	   
	   for(NavigationChannelList nav : filteredNavigationList) {  
		   
		   double lat = nav.getLatitude();
		   double longitude = nav.getLongitude();
		   
		   if(lat > maxLat) {maxLat = lat;}
		   else if(lat < minLat) {minLat = lat;}
		
		   if(longitude > maxLong) {maxLong = longitude;}
		   else if(longitude < minLong) {minLong = longitude;}		   
	   }
	   
	   double[] boundNW = new double[3];
	   double[] boundNE = new double[3];
	   double[] boundSW = new double[3];
	   double[] boundSE = new double[3];
	   
	   boundNW[0] = maxLat + .02; boundNW[1] = minLong - .02;
	   boundNE[0] = maxLat + .02; boundNE[1] = maxLong + .02;
	   boundSW[0] = minLat - .02; boundSW[1] = minLong - .02;
	   boundSE[0] = minLat - .02; boundSE[1] = maxLong + .02;
	   
       double[] v0 = NavigationChannelList.GeodeticToECEF(boundNW[0], boundNW[1], 0, localOrigin);
       double[] v1 = NavigationChannelList.GeodeticToECEF(boundNE[0], boundNE[1], 0, localOrigin);
       double[] v2 = NavigationChannelList.GeodeticToECEF(boundSW[0], boundSW[1], 0, localOrigin);
       double[] v3 = NavigationChannelList.GeodeticToECEF(boundSE[0], boundSE[1], 0, localOrigin);

       System.out.println(v0[0] + " " + v0[1] + " " + v0[2]);
       System.out.println(v1[0] + " " + v1[1] + " " + v1[2]);
       System.out.println(v2[0] + " " + v2[1] + " " + v2[2]);
       System.out.println(v3[0] + " " + v3[1] + " " + v3[2]);
       
       bounds_in = new Matrix(4,3);
       
       bounds_in.setRow(0, v0);
       bounds_in.setRow(1, v1);
       bounds_in.setRow(2, v2);
       bounds_in.setRow(3, v3); 	   
    	
    }
    
    
	/**
	 * Once the navigation log file has been read, and FDOA filtering has been 
	 * applied using filterFDOAfromNavigationList, this function will estimate the solution using the given FDOAs 
	 * at the given anchors and computed velocities. The source solution is estimated
	 * via means of a stochastic gradient descent approach. All of the given anchors will be
	 * used for finding the gradient descent solution 
	 * The resulting velocities are stored in the navigation list 
	 * @see  filterFDOAfromNavigationList
	 * @throws Exception flteredNavigationList must be nonempty
	 */
      
    public void estimateSourceSolutionFDOA() throws Exception {
    	
    	if(filteredNavigationList.size() == 0) {
    		throw new Exception("Must first apply FDOA filtering to get "
    				+ "filtered Navigation list");
    	}
    	
    	Anchors myAnchors = new Anchors();
	    ArrayList<Double> fdoas = new ArrayList<Double>();
    	
	    for(int i = 0; i < filteredNavigationList.size(); i++) {
	    
	      NavigationChannelList nav = filteredNavigationList.get(i);	
	      
	      double[] s = nav.getLocalECEF();
	      double[] v = nav.getVelocity();
	      
	      myAnchors.setCoordinatesAndVelocity(s, v);
	      fdoas.add(nav.getFDOA());
	    }
	    
	    myAnchors.commitCoordinatesAndVelocity();
	    int num_anchors = myAnchors.getNumberOfAnchors();
	    
        Matrix ranges = new Matrix(num_anchors);        
        for (int i = 0; i < num_anchors; i++) {
        	
            ranges.w[i] = fdoas.get(i).doubleValue();
        }

                    
        double[] source = null;
        if(sourceLatitude != 0.0) {
        	
        	source = NavigationChannelList.GeodeticToECEF(sourceLatitude, sourceLongitude, 0);
    		source[0] -= localOrigin[0];
    		source[1] -= localOrigin[1];
    		source[2] -= localOrigin[2];	
        }
        
        GradDescentResult gdescent_result = GradDescentResult.mlatFdoa(myAnchors.getAnchors(), 
        		myAnchors.getVelocities(), ranges, bounds_in, n_trial, alpha, time_threshold, source);
	    
        System.out.println("\nEstimator");
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
        gdescent_result.estimator.printMatrix();
        
        
        for(int j = 0; j < gdescent_result.estimator_candidate.rows; j++) {
        	gdescent_result.estimator_candidate.transformRowCoord(j,localOrigin); 
        }
    	
        
        ArrayList<double[]> estimates = new ArrayList<double[]>();
        ArrayList<Double> error_est = new ArrayList<Double>();
                     
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
                
        final Kml kml = createSolutionDocument(estimates, error_est);
        kml.marshal(new File("SolutionMarkers.kml"));        
    }
    
    
	/**
	 * Once the navigation log file has been read, and TDOA filtering has been 
	 * applied using filterTDOAfromNavigationList, this function will estimate the 
	 * solution using the given anchors and their TDOA estimates. The source solution is estimated
	 * via means of a stochastic gradient descent approach. All of the given anchors will be
	 * used for finding the gradient descent solution.
	 * @see  filterTDOAfromNavigationList
	 * @throws Exception flteredNavigationList must be nonempty
	 */
      
	
	public void estimateSourceSolutionTDOA() throws Exception {
	
        int num_anchors;
        
        final Kml kml0 = createCIRDocument(filteredNavigationList);
		kml0.marshal(new File("FilteredCIRMarkers.kml"));
        
		Anchors myAnchors = new Anchors();
	    ArrayList<Double> rtdoas = new ArrayList<Double>();
	    
	    localOrigin = NavigationChannelList.GeodeticToECEF(tdoaFrequency0.get(0).getLatitude(), 
	    		tdoaFrequency0.get(0).getLongitude(), 0);
	    
	    
        double[] source = null;
        if(sourceLatitude != 0.0) {
        	
        	source = NavigationChannelList.GeodeticToECEF(sourceLatitude, sourceLongitude, 0);
    		source[0] -= localOrigin[0];
    		source[1] -= localOrigin[1];
    		source[2] -= localOrigin[2];
        }

		for(int i = 0; i < tdoaFrequency0.size(); i++) {
						
			TimeDiffOnArrival tdoa = tdoaFrequency0.get(i);
					
			double[] locs = NavigationChannelList.GeodeticToECEF(tdoa.getLatitude(), tdoa.getLongitude(), 0);
			
			myAnchors.setCoordinates(locs[0] - localOrigin[0], 
					                 locs[1] - localOrigin[1], 
					                 locs[2] - localOrigin[2]);
						
			rtdoas.add(tdoa.getTDOA()*C);
			
		}
		System.out.println("Origin: " + localOrigin[0] + " " + localOrigin[1] + " " + localOrigin[2]);

		myAnchors.commitCoordinates();
		num_anchors = myAnchors.getNumberOfAnchors();	
        Matrix ranges_with_error = new Matrix(num_anchors);
        
        for (int i = 0; i < num_anchors; i++) {
         
            ranges_with_error.w[i] = rtdoas.get(i);
        }

        
        GradDescentResult gdescent_result = GradDescentResult.mlatTdoa(myAnchors, ranges_with_error, bounds_in, 
        		n_trial, alpha, time_threshold, source);

        
        System.out.println("\nEstimator");
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
        gdescent_result.estimator.printMatrix();
        
        gdescent_result.error.printMatrix();
        
        for(int j = 0; j < gdescent_result.estimator_candidate.rows; j++) {
        	gdescent_result.estimator_candidate.transformRowCoord(j,localOrigin); 
        }
        
        
        ArrayList<double[]> estimates = new ArrayList<double[]>();
        ArrayList<Double> error_est = new ArrayList<Double>();
        
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
        
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        
        
        final Kml kml = createSolutionDocument(estimates, error_est);
        kml.marshal(new File("SolutionMarkers.kml"));
        
	}
	
	/**
	 * Once the navigation log file has been read, and TDOA filtering has been 
	 * applied using filterTDOAfromNavigationList, this function will estimate the solution 
	 * using the given anchors and their TDOA estimates in a successive manner. Namely, 
	 * an estimate of the source will be produced for every successive anchor/tdoa added from 
	 * the filtered navigation queue. This an estimate will be produced for each additional 
	 * anchor/tdoa. The solution estimates will be found in a .kml file. 
	 * @see  filterTDOAfromNavigationList
	 * @throws Exception flteredNavigationList must be nonempty
	 */	
	public void estimateAdaptiveSourceTDOA() throws Exception {
		
        int num_anchors;
        
        
        final Kml kml0 = createCIRDocument(filteredNavigationList);
		kml0.marshal(new File("FilteredCIRMarkers.kml"));
        
        
		Anchors myAnchors = new Anchors();
	    ArrayList<Double> rtdoas = new ArrayList<Double>();
	    
	    localOrigin = NavigationChannelList.GeodeticToECEF(tdoaFrequency0.get(0).getLatitude(), 
	    		tdoaFrequency0.get(0).getLongitude(), 0);
	    
        double[] source = null;
        if(sourceLatitude != 0.0) {
        	
        	source = NavigationChannelList.GeodeticToECEF(sourceLatitude, sourceLongitude, 0);
    		source[0] -= localOrigin[0];
    		source[1] -= localOrigin[1];
    		source[2] -= localOrigin[2];
        }
	
		
		for(int i = 0; i < tdoaFrequency0.size(); i++) {
						
			TimeDiffOnArrival tdoa = tdoaFrequency0.get(i);
					
			double[] locs = NavigationChannelList.GeodeticToECEF(tdoa.getLatitude(), tdoa.getLongitude(), 0);
			
			myAnchors.setCoordinates(locs[0] - localOrigin[0], 
					                 locs[1] - localOrigin[1], 
					                 locs[2] - localOrigin[2]);
						
			rtdoas.add(tdoa.getTDOA()*C);
			
		}
		System.out.println("Origin: " + localOrigin[0] + " " + localOrigin[1] + " " + localOrigin[2]);
		
		myAnchors.commitCoordinates();
		num_anchors = myAnchors.getNumberOfAnchors();
					
        Matrix ranges_with_error = new Matrix(num_anchors);
        
        for (int i = 0; i < num_anchors; i++) {
            ranges_with_error.w[i] = rtdoas.get(i);            
        }

       ArrayList<double[]> estimates = new ArrayList<double[]>();
       ArrayList<Double> error_est = new ArrayList<Double>();
       
       
       for(int i = 10; i < num_anchors; i++) {

    	Anchors updateAnchors = myAnchors.subset(i);
    	Matrix updateRanges = ranges_with_error.subset(i);
        
        GradDescentResult gdescent_result = GradDescentResult.mlatTdoa(updateAnchors, updateRanges, bounds_in, 
        		n_trial, alpha, time_threshold, source);

        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
                
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        System.out.print("Error with " + i + " anchor navigation nodes: " + gdescent_result.error.w[0] 
        		+ " " + updateRanges.w[updateRanges.w.length-1] + " ");
        updateAnchors.getRow(updateAnchors.getAnchors().rows - 1).printMatrix(); 
         
       } 
       
       double[] stockArr = new double[error_est.size()];
       double[] x = new double[error_est.size()];
       for(int i = 0; i < error_est.size(); i++) {
    	   stockArr[i] = error_est.get(i).doubleValue();
    	   x[i] = i;
       }
       
       Plot2DPanel plot = new Plot2DPanel();
		 
       // add a line plot to the PlotPanel
       plot.addLinePlot("TDOA plot", x, stockArr);
       JFrame frame = new JFrame("Cost function value");
       frame.setSize(900, 700);
       frame.setContentPane(plot);
       frame.setVisible(true);

        
       final Kml kml = createSolutionDocument(estimates, error_est);
       kml.marshal(new File("SolutionMarkers.kml"));
       
	}
		
	
	
    public void estimateAdaptiveSourceWithDrift() throws Exception {
		
        int num_anchors;
        
        
        final Kml kml0 = createCIRDocument(filteredNavigationList);
		kml0.marshal(new File("FilteredCIRMarkers.kml"));
        
        
		Anchors myAnchors = new Anchors();
	    ArrayList<Double> tdoas = new ArrayList<Double>();
	    ArrayList<Double> rtdoas = new ArrayList<Double>();
	    
	    localOrigin = NavigationChannelList.GeodeticToECEF(tdoaFrequency0.get(0).getLatitude(), 
	    		tdoaFrequency0.get(0).getLongitude(), 0);
	    
		double[] source = NavigationChannelList.GeodeticToECEF(47.184733, 7.707917, 0);
		source[0] -= localOrigin[0];
		source[1] -= localOrigin[1];
		source[2] -= localOrigin[2];
	    

		for(int i = 0; i < tdoaFrequency0.size(); i++) {
						
			TimeDiffOnArrival tdoa = tdoaFrequency0.get(i);
					
			double[] locs = NavigationChannelList.GeodeticToECEF(tdoa.getLatitude(), tdoa.getLongitude(), 0);
			
			myAnchors.setCoordinates(locs[0] - localOrigin[0], 
					                 locs[1] - localOrigin[1], 
					                 locs[2] - localOrigin[2], 
					                 0.0);
						
			rtdoas.add(tdoa.getTDOA()*C);
			tdoas.add(tdoa.timeDiff*C);
			
		}
		System.out.println("Origin: " + localOrigin[0] + " " + localOrigin[1] + " " + localOrigin[2]);
		System.out.println("Source: " + source[0] + " " + source[1] + " " + source[2]);
		
		myAnchors.commitCoordinates();
		num_anchors = myAnchors.getNumberOfAnchors();

			
		
        Matrix ranges_with_error = new Matrix(num_anchors);
        Matrix timeDiffs = new Matrix(num_anchors);
        
        for (int i = 0; i < num_anchors; i++) {
        	
            ranges_with_error.w[i] = rtdoas.get(i).doubleValue();
            timeDiffs.w[i] = tdoas.get(i).doubleValue();
            
        }

     
       ArrayList<double[]> estimates = new ArrayList<double[]>();
       ArrayList<Double> error_est = new ArrayList<Double>();
       
       
       for(int i = 10; i < num_anchors; i++) {

    	Anchors updateAnchors = myAnchors.subset(i);
    	Matrix updateRanges = ranges_with_error.subset(i);
    	Matrix updateTimeDiff = timeDiffs.subset(i);
        
        GradDescentResult gdescent_result = GradDescentResult.mlatTdoaDrift(updateAnchors, updateRanges, updateTimeDiff.w, 
        		bounds_in, n_trial, alpha, time_threshold, source);

        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
                
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        System.out.print("Error with " + i + " anchor navigation nodes: " + gdescent_result.error.w[0] 
        		+ " "); 
        gdescent_result.estimator.printMatrix(); 
       } 
       
       double[] stockArr = new double[error_est.size()];
       double[] x = new double[error_est.size()];
       for(int i = 0; i < error_est.size(); i++) {
    	   stockArr[i] = error_est.get(i).doubleValue();
    	   x[i] = i;
       }
       
       Plot2DPanel plot = new Plot2DPanel();
		 
       // add a line plot to the PlotPanel
       plot.addLinePlot("TDOA plot", x, stockArr);
       JFrame frame = new JFrame("Cost function value");
       frame.setSize(900, 700);
       frame.setContentPane(plot);
       frame.setVisible(true);

        
       final Kml kml = createSolutionDocument(estimates, error_est);
       kml.marshal(new File("SolutionMarkers.kml"));
       
	}
	
	
	
	
	
	/**
	 * Once the navigation log file has been read, this function will estimate the solution 
	 * using the given anchors and their TDOA estimates in a successive manner where the first
	 * reference anchor changes. Namely, an estimate of the source will be produced for every successive anchor/tdoa added from 
	 * the filtered navigation queue. The navigation queue will be split up in at most n_estimates
	 * separate navigation queues, where the first anchor will be the reference anchor for the TDOA
	 * computation. Thus an estimate will be produced for each n_estimate navigation queue. 
	 * This is to study the effects of the reference location with respect to the source and tdoa
	 * computation. The solution estimates will be found in a .kml file called SolutionMarkers_i where
	 * i is the i-th segment of the navigation list.  
	 * @param freq The frequency index to use (0,1,2)  
	 * @param n_estimates The number of sections the navigation path will be split up into. 
	 * The reference point for tdoa computation will be the first in the segmented navigation queue. 
	 * @param thresh The threshhold used for eliminating outliers in first-order differencing 
	 * @see  filterTDOAfromNavigationList
	 * @throws Exception flteredNavigationList must be nonempty
	 */		
	public void estimateDynamicReferenceTDOA(int freq, int n_estimates, double thresh) throws Exception {
		
				
        double time = 0; 
        double refTime = -1.0;
        double tdoa = 0.0;
        double frameLength = .065536;
        double[] x;
        
        int outlierCount = 0;
        localOrigin = NavigationChannelList.GeodeticToECEF(filteredNavigationList.get(0).getLatitude(), 
        		filteredNavigationList.get(0).getLongitude(), 0);
	    
        double[] source = null;
        if(sourceLatitude != 0.0) {
        	
        	source = NavigationChannelList.GeodeticToECEF(sourceLatitude, sourceLongitude, 0);
    		source[0] -= localOrigin[0];
    		source[1] -= localOrigin[1];
    		source[2] -= localOrigin[2];	
        }
        
        int count = 0; 
		int n_samples = (int)Math.floor(filteredNavigationList.size()/n_estimates); 
		tdoaFrequency0 = new ArrayList<TimeDiffOnArrival>();
		x = new double[n_samples];
		
		System.out.println("Number of samples: " + n_samples);
		
        for (NavigationChannelList navList : filteredNavigationList) {
        	
            time = navList.getTimeStampAtFreq(freq);
            if(refTime < 0 && time > 0) { 
            	
                refTime = time;  //set reference time for first measurement            
                tdoa = 0.0;
                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa, 0);                
                tdoaFrequency0.add(td);                
            }
            else if(time > 0) {
                            	
                double timeStamp = (time - refTime) % frameLength;
                if (timeStamp < frameLength / 2) {
                    tdoa = timeStamp;
                }
                else {
                    tdoa = - (frameLength - timeStamp);
                }
                
                if(Math.abs(tdoa - tdoaFrequency0.get(tdoaFrequency0.size()-1).getTDOA()) < thresh) {
                	
       
	                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa, time - refTime);                
	                tdoaFrequency0.add(td);
	                outlierCount = 0;
                }
                else if(outlierCount > maxOutlierThreshold) {
                	
                	
	                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa, time - refTime);                
	                tdoaFrequency0.add(td);
	                outlierCount = 0;
                }
                else {
                	outlierCount++;
                }
                               
            }
            
            
            //If enough samples collected, estimate new source//
            if(tdoaFrequency0.size() == n_samples) {
            	
            	ArrayList<double[]> estimates = new ArrayList<double[]>();
                ArrayList<Double> error_est = new ArrayList<Double>();

                ArrayList<double[]> coordinates = new ArrayList<double[]>();
            	Anchors updateAnchors = new Anchors();
            	Matrix updateRanges = new Matrix(n_samples);

            	for(int i = 0; i < tdoaFrequency0.size(); i++) {
					        					
        			double[] locs = NavigationChannelList.GeodeticToECEF(tdoaFrequency0.get(i).getLatitude(), 
        					tdoaFrequency0.get(i).getLongitude(), 0);
        			
        			updateAnchors.setCoordinates(locs[0] - localOrigin[0], 
        					                 locs[1] - localOrigin[1], 
        					                 locs[2] - localOrigin[2]);
        			
        			double[] coord = {tdoaFrequency0.get(i).getLongitude(), tdoaFrequency0.get(i).getLatitude(), 0}; 
        			coordinates.add(coord);
        			
        			updateRanges.w[i] = tdoaFrequency0.get(i).getTDOA()*C; 
        			x[i] = i;
        		}

        		updateAnchors.commitCoordinates();
        		GradDescentResult gdescent_result = GradDescentResult.mlatTdoa(updateAnchors, updateRanges, bounds_in, 
                		n_trial, alpha, time_threshold, source);

                gdescent_result.estimator.transformRowCoord(0,localOrigin);         
                
    	        estimates.add(gdescent_result.estimator.w);
    	        error_est.add(gdescent_result.error.w[0]);
    	        System.out.print("Error with " + n_samples + " anchor navigation nodes: " + gdescent_result.error.w[0] 
    	        		+ " " + updateRanges.w[updateRanges.w.length-1] + " ");
    	        updateAnchors.getRow(updateAnchors.getAnchors().rows - 1).printMatrix(); 
            	
            	//Output to GoogleEarth the estimate and the navigation// 
    	        final Kml kml = createNavigationAndSolutionDocument(coordinates, gdescent_result.estimator.w, count);
    	        kml.marshal(new File("SolutionMarkers_" + count +".kml"));
            	
    	        
    	        if(plotTDOAs) {
    	        	
    		        Plot2DPanel plot = new Plot2DPanel();
    		        JFrame frame = new JFrame("Plot of TDOAs " + count);
    		        plot.addLinePlot("TDOA plot", x, updateRanges.w);	        
    		        frame.setSize(900, 700);
    		        frame.setContentPane(plot);
    		        frame.setVisible(true);	
    	        }
    	        
    	        tdoaFrequency0.clear();
    	        count++;
    	        refTime = -1.0;
    	        tdoa = 0.0;
            }   
        }
        	
	}
	
	
	/**
	 * Once the navigation log file has been read, and FDOA filtering has been 
	 * applied using filterFDOAfromNavigationList, this function will estimate the 
	 * solution using the given FDOAs at the given anchors and computed velocities. The source solution is estimated
	 * via means of a stochastic gradient descent approach. All of the given anchors will be
	 * used for finding the gradient descent solution 
	 * The resulting velocities are stored in the navigation list 
	 * @see  filterFDOAfromNavigationList
	 * @throws Exception flteredNavigationList must be nonempty
	 */
	
    public void estimateAdaptiveSourceFDOA() throws Exception {
		
    	Anchors myAnchors = new Anchors();
	    ArrayList<Double> fdoas = new ArrayList<Double>();
    	
	    for(int i = 0; i < filteredNavigationList.size(); i++) {
	    
	      NavigationChannelList nav = filteredNavigationList.get(i);	
	      
	      double[] s = nav.getLocalECEF();
	      double[] v = nav.getVelocity();
	      
	      myAnchors.setCoordinatesAndVelocity(s, v);
	      fdoas.add(nav.getFDOA());
	    
	      System.out.println(nav.getFDOA());
	    }
	    
	    myAnchors.commitCoordinatesAndVelocity();
	    int num_anchors = myAnchors.getNumberOfAnchors();
	    
        Matrix ranges = new Matrix(num_anchors);        
        for (int i = 0; i < num_anchors; i++) {
        	
            ranges.w[i] = fdoas.get(i).doubleValue();
        }
        
       
        double[] source = NavigationChannelList.GeodeticToECEF(47.184733, 7.707917, 0);
		source[0] -= localOrigin[0];
		source[1] -= localOrigin[1];
		source[2] -= localOrigin[2];

        
        int n_trial = 3000; 
        double alpha = 0.001; 
        double time_threshold = 50000;
        
     
       ArrayList<double[]> estimates = new ArrayList<double[]>();
       ArrayList<Double> error_est = new ArrayList<Double>();
       
       
       for(int i = 10; i < num_anchors; i=i+20) {

    	Matrix updateAnchors = myAnchors.subsetCoordinates(i);
    	Matrix updateVelocities = myAnchors.subsetVelocity(i);
    	Matrix updateRanges = ranges.subset(i);
        
    	
    	GradDescentResult gdescent_result = GradDescentResult.mlatFdoa(updateAnchors, 
    			updateVelocities, updateRanges, bounds_in, n_trial, alpha, time_threshold, source);
    	
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
                
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        System.out.print("Error with " + i + " anchor navigation nodes: " + gdescent_result.error.w[0] 
        		+ " " + updateRanges.w[updateRanges.w.length-1] + " ");
        updateAnchors.getRow(updateAnchors.rows - 1).printMatrix(); 
         
       } 
       
       double[] stockArr = new double[error_est.size()];
       double[] x = new double[error_est.size()];
       for(int i = 0; i < error_est.size(); i++) {
    	   stockArr[i] = error_est.get(i).doubleValue();
    	   x[i] = i;
       }
       
       Plot2DPanel plot = new Plot2DPanel();       
       JFrame frame = new JFrame("FDOA Cost function value");
       plot.addLinePlot("TDOA plot", x, stockArr);
       frame.setSize(900, 700);
       frame.setContentPane(plot);
       frame.setVisible(true);

        
       final Kml kml = createSolutionDocument(estimates, error_est);
       kml.marshal(new File("SolutionMarkers.kml"));
       
	}
	
    
    
    
    
	
    public void estimateAdaptiveSourceRSSI(double txPower) throws Exception {
		
    	
    	final Kml kml0 = createCIRDocument(filteredNavigationList);
		kml0.marshal(new File("FilteredCIRMarkers.kml"));
    	
    	Anchors myAnchors = new Anchors();
	    ArrayList<Double> distance = new ArrayList<Double>();
    	
	    for(int i = 0; i < filteredNavigationList.size(); i++) {
	    
	      NavigationChannelList nav = filteredNavigationList.get(i);		      
	      double[] s = nav.getLocalECEF();	      
	      myAnchors.setCoordinates(s);

	      distance.add(getDistance(nav.getRSSI(0), txPower));
	    
	      System.out.println(distance.get(distance.size() - 1));
	    }
	    
	    myAnchors.commitCoordinates();
	    int num_anchors = myAnchors.getNumberOfAnchors();
	    
	    double[] xv = new double[num_anchors];
        Matrix ranges = new Matrix(num_anchors);        
        for (int i = 0; i < num_anchors; i++) {        	
            ranges.w[i] = distance.get(i).doubleValue();
            xv[i] = i;
        }
        
        Plot2DPanel plotRssi = new Plot2DPanel();       
        JFrame frameRssi = new JFrame("Distance values");
        plotRssi.addLinePlot(" ", xv, ranges.w);
        frameRssi.setSize(900, 700);
        frameRssi.setContentPane(plotRssi);
        frameRssi.setVisible(true);
        
        
       
        double[] source = NavigationChannelList.GeodeticToECEF(47.184733, 7.707917, 0);
		source[0] -= localOrigin[0];
		source[1] -= localOrigin[1];
		source[2] -= localOrigin[2];

        
        int n_trial = 3000; 
        double alpha = 0.001; 
        double time_threshold = 50000;
        
     
       ArrayList<double[]> estimates = new ArrayList<double[]>();
       ArrayList<Double> error_est = new ArrayList<Double>();
       
       
       for(int i = 10; i < num_anchors; i++) {

    	Anchors updateAnchors = myAnchors.subset(i);
    	Matrix updateRanges = ranges.subset(i);
        
    	GradDescentResult gdescent_result = GradDescentResult.mlat(updateAnchors, updateRanges, bounds_in, n_trial, alpha, time_threshold, source);
    	
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
                
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        System.out.print("Error with " + i + " anchor navigation nodes: " + gdescent_result.error.w[0] 
        		+ " " + updateRanges.w[updateRanges.w.length-1] + " ");
        updateAnchors.getRow(updateAnchors.getAnchors().rows - 1).printMatrix(); 
         
       } 
       
       double[] stockArr = new double[error_est.size()];
       double[] x = new double[error_est.size()];
       for(int i = 0; i < error_est.size(); i++) {
    	   stockArr[i] = error_est.get(i).doubleValue();
    	   x[i] = i;
       }
       
       Plot2DPanel plot = new Plot2DPanel();       
       JFrame frame = new JFrame("Rssi Cost function value");
       plot.addLinePlot("TDOA plot", x, stockArr);
       frame.setSize(900, 700);
       frame.setContentPane(plot);
       frame.setVisible(true);

       
       final Kml kml = createSolutionDocument(estimates, error_est);
       kml.marshal(new File("SolutionMarkers.kml"));
       
	}
    
    double getDistance(double rssi, double txPower) {
        /*
         * RSSI = TxPower - 10 * n * lg(d)
         * n = 2 (in free space)
         * 
         * d = 10 ^ ((TxPower - RSSI) / (10 * n))
         */
     
        return Math.pow(10.0, (txPower - rssi) / (20.0));
    }

    
	/**
	 * Once the navigation log file has been read, and FDOA filtering has been 
	 * applied using filterTDOAfromNavigationList, this function will estimate the solution 
	 * using the given anchors and also compute the exact Doppler shift values to compare with
	 * the ones extracted from the .log file. The source must be known in this case is set using
	 * setKnownSource. The function will plot the difference of exact FDOA and  
	 *  
	 * @see setKnownSource  
	 * @see  filterFDOAfromNavigationList
	 * @throws Exception flteredNavigationList must be nonempty
	 */	
    
	public void testSGDfdoa() throws Exception {
		
		DescriptiveStatistics stats = new DescriptiveStatistics();
		
		Anchors myAnchors = new Anchors();
    	double[] fds = new double[filteredNavigationList.size()];
    	double[] xvs = new double[filteredNavigationList.size()];
    	
    	
	    for(int i = 0; i < filteredNavigationList.size(); i++) {
	    
	      NavigationChannelList nav = filteredNavigationList.get(i);	
	      
	      double[] s = nav.getLocalECEF();
	      double[] v = nav.getVelocity();
	      
	      filteredNavigationList.get(i).setFrequencyError(2);
	      myAnchors.setCoordinatesAndVelocity(s, v);	    
	      System.out.println(nav.getFDOA());
	      fds[i] = nav.getFDOA();
	      xvs[i] = (double)i;
	      stats.addValue(fds[i]);
	    }
	    
	    myAnchors.commitCoordinatesAndVelocity();
	    int num_anchors = myAnchors.getNumberOfAnchors();
	    
	    Matrix anchors_in = myAnchors.getAnchors();
	    Matrix velocities = myAnchors.getVelocities();
	    double[] source = NavigationChannelList.GeodeticToECEF(46.762878, 7.600469, 0);
        
	    //double[] source = NavigationChannelList.GeodeticToECEF(47.184733, 7.707917, 0);
		source[0] -= localOrigin[0];
		source[1] -= localOrigin[1];
		source[2] -= localOrigin[2];
		
		double[] xv = new double[num_anchors];
	    Matrix sol = new Matrix(source, 1);
	    Matrix ranges = new Matrix(num_anchors,1);           
        ranges.w[0] = 0;
        
        // Compute some statistics
        double mean = stats.getMean();
        
        noiseModel = new double[num_anchors];
        for (int j = 0; j < num_anchors-1; j++) {  
        	
        	double fdoa_est = GradDescentResult.fdoaEstimate0(anchors_in.getRow(j+1), 
        			velocities.getRow(j+1), sol);
        	
        	ranges.w[j+1] = fdoa_est/10 + mean;
        	noiseModel[j+1] = Math.abs(ranges.w[j+1] - fds[j+1]);
        	xv[j+1] = j+1;
        }
	    System.out.println(num_anchors + " " + filteredNavigationList.size());
        
        
        Plot2DPanel plotfdoa = new Plot2DPanel();
        JFrame framefdoa = new JFrame("Plot of Simulated FDOAs");
        plotfdoa.addLinePlot("Plot of Simulation FDOAs", xv, ranges.w);	        
        framefdoa.setSize(900, 700);
        framefdoa.setContentPane(plotfdoa);
        framefdoa.setVisible(true);
        
        Plot2DPanel plotfdoa2 = new Plot2DPanel();
        JFrame framefdoa2 = new JFrame("Plot of FDOAs");
        plotfdoa2.addLinePlot("Plot of FDOAs", xvs, fds);	        
        framefdoa2.setSize(900, 700);
        framefdoa2.setContentPane(plotfdoa2);
        framefdoa2.setVisible(true);
        
        Plot2DPanel plotfdoa3 = new Plot2DPanel();
        JFrame framefdoa3 = new JFrame("Plot of FDOAs");
        plotfdoa3.addLinePlot("Plot of Residual error", xvs, noiseModel);	        
        framefdoa3.setSize(900, 700);
        framefdoa3.setContentPane(plotfdoa3);
        framefdoa3.setVisible(true);
        
     
       ArrayList<double[]> estimates = new ArrayList<double[]>();
       ArrayList<Double> error_est = new ArrayList<Double>();
       
       
       for(int i = 10; i < num_anchors; i++) {

    	Matrix updateAnchors = myAnchors.subsetCoordinates(i);
    	Matrix updateVelocities = myAnchors.subsetVelocity(i);
    	Matrix updateRanges = ranges.subset(i);
        
    	
    	GradDescentResult gdescent_result = GradDescentResult.mlatFdoa(updateAnchors, 
    			updateVelocities, updateRanges, bounds_in, n_trial, alpha, time_threshold, source);
    	
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
                
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        System.out.print("Error with " + i + " anchor navigation nodes: " + gdescent_result.error.w[0] 
        		+ " " + updateRanges.w[updateRanges.w.length-1] + " ");
        updateAnchors.getRow(updateAnchors.rows - 1).printMatrix(); 
         
       } 
       
       double[] stockArr = new double[error_est.size()];
       double[] x = new double[error_est.size()];
       for(int i = 0; i < error_est.size(); i++) {
    	   stockArr[i] = error_est.get(i).doubleValue();
    	   x[i] = i;
       }
       
        
       final Kml kml = createSolutionDocument(estimates, error_est);
       kml.marshal(new File("SolutionMarkers.kml"));        		
		
	}
    
    
	
	public void testSGDtdoa(File file) throws Exception {
		
		String strline; 
		String[] tokens; 
		String delims = "\\s+";
		
		ArrayList<double[]> locs = new ArrayList<double[]>();
		FileInputStream fin; DataInputStream din; BufferedReader br;
		fin = new FileInputStream(file);
        din = new DataInputStream(fin);
        br = new BufferedReader(new InputStreamReader(din));
		
        strline = br.readLine();
        tokens = strline.split(delims);
        double[] loc0 = {(new Double(tokens[0])).doubleValue(), 
		        (new Double(tokens[1])).doubleValue(),
		        (new Double(tokens[2])).doubleValue()};
        
        double[] origin = NavigationChannelList.GeodeticToECEF(loc0[0], loc0[1], loc0[2]);        
        double[] orig = {0.0, 0.0, 0.0}; 
        locs.add(orig);
        
        
        while((strline = br.readLine()) != null) { 
        	
        	tokens = strline.split(delims);
        	double[] loc = {(new Double(tokens[0])).doubleValue(), 
        			        (new Double(tokens[1])).doubleValue(),
        			        (new Double(tokens[2])).doubleValue()};
        	
        	//System.out.println(loc[0] + " " + loc[1] + " " + loc[2]);
        	double[] ecef = NavigationChannelList.GeodeticToECEF(loc[0], loc[1], loc[2]);
        	
        	ecef[0] -= origin[0];
        	ecef[1] -= origin[1];
        	ecef[2] -= origin[2];
        	
        	//System.out.println(ecef[0] + " " + ecef[1] + " " + ecef[2]);
        	locs.add(ecef);	
        	
        	
        	ecef[0] += origin[0];
        	ecef[1] += origin[1];
        	ecef[2] += origin[2];
        	
        	double[] back = Coord.ecef_to_geo(ecef);
        	
        	System.out.println(back[0] + " " + back[1] + " " + back[2]);
        }
        br.close(); 	
        
        double[] source = {-2762.18694234, -5068.57649907, 0};
        
        source[0] += origin[0];
        source[1] += origin[1];
        source[2] += origin[2];
        
        double[] backOrig = Coord.ecef_to_geo(source);
        System.out.println("Origin in lat/long");
        System.out.println(backOrig[0] + " " + backOrig[1] + " " + backOrig[2]);
                
	}
	
	
	public void testCaseData(File file) throws Exception {
		

		String strline; 
		String[] tokens; 
		String delims = "\\s+";
		
		FileInputStream fin; DataInputStream din; BufferedReader br;
		fin = new FileInputStream(file);
        din = new DataInputStream(fin);
        br = new BufferedReader(new InputStreamReader(din));
		
        
        int num_anchors;
		Anchors myAnchors = new Anchors();
	    ArrayList<Double> tdoas = new ArrayList<Double>();
        

        while((strline = br.readLine()) != null) {        	
        	
        	tokens = strline.split(delims);
        	myAnchors.setCoordinates((new Double(tokens[1])).doubleValue(), 
        			                 (new Double(tokens[2])).doubleValue(),
        			                 (new Double(tokens[3])).doubleValue());
        	
        	//tdoas.add((new Double(tokens[7])).doubleValue());
        }
        br.close();
        
        myAnchors.commitCoordinates();
		num_anchors = myAnchors.getNumberOfAnchors();
		   
		Matrix node = new Matrix(3);
		double[] n = {-2762.18694234, -5068.57649907, 2522.74598415};
		node.setPointer(n);
		
		double[] ctdoas = new double[num_anchors];
		
		tdoas.add(0.0);
		for(int i = 0; i < num_anchors-1; i++) {
			
			double ctdoa = Mstat.distance(node, myAnchors.getRow(i+1)) - 
		                   Mstat.distance(node, myAnchors.getRow(i));
			
			ctdoa = ctdoa/C; 	
			ctdoas[i+1] = ctdoa; 
			tdoas.add(ctdoa);		
		}
		
				
        Matrix ranges = new Matrix(num_anchors);
        Matrix ranges_with_error = new Matrix(num_anchors);
        
        for (int i = 0; i < num_anchors; i++) {
        	
            ranges.w[i] = tdoas.get(i).doubleValue()*C;
            ranges_with_error.w[i] = ranges.w[i];
        }
        ranges_with_error.w[0] = 0; 
        
        int n_trial = 5000; 
        double alpha = 0.0001; 
        double time_threshold = 500000;
        
        Matrix bounds_in = new Matrix(2,3);
        bounds_in.set(0, 0, myAnchors.getColumnMin(0) - 100.0);
        bounds_in.set(0, 1, myAnchors.getColumnMin(1) - 100.0);
        bounds_in.set(0, 2, myAnchors.getColumnMin(2));
        
        bounds_in.set(1, 0, myAnchors.getColumnMax(0) + 100.0);
        bounds_in.set(1, 1, myAnchors.getColumnMax(1) + 100.0);
        bounds_in.set(1, 2, myAnchors.getColumnMax(2));
        
        GradDescentResult gdescent_result = GradDescentResult.mlatTdoa(myAnchors, ranges_with_error, bounds_in, 
        		n_trial, alpha, time_threshold, node.w);

        System.out.println("Anchors");
        myAnchors.printMatrix();
        
        System.out.println("Ranges");
        ranges.printMatrix();
        
        
        System.out.println("\nEstimator");
        gdescent_result.estimator.printMatrix();
        
        System.out.println("\nTrue location: " + "-2762.18694234 -5068.57649907  2522.74598415");
        
	}
	
	
	/**
	 * Function to print out the coordinates of the navigation file, if non empty
	 *
	 * @throws Exception navigation list must be non-empty
	 */
	
	public void printNavigationHistory() throws Exception {
		
		if(navigationList == null) {
			throw new Exception("NavigationList must be nonempty");
		}
		
		System.out.println("My size: " + navigationList.size());
		
		for(NavigationChannelList nav : navigationList) {
			nav.printCoordinates();
		}
	}
	
	
	/**
	 * Main function that will parse a navigation log file (*.cir or *.log)
	 * into a navigation list that holds all recordings of longitude and latitude
	 * in the form of a NavigationChannelList. This Channel list for each coordinate
	 * will then give all the frequency index information and store the peak lists for each 
	 * frequency.   
	 * <p>
	 * The output of the function will be a list of NavigationChannelList which can then be
	 * filtered on using the different filtering functions
	 *
	 * @param file  A navigation file 
	 * @throws Exception if the file is not in the correct format
	 * @see filterUnique NavigationChannelList
	 */
	
	
	public void createNavigationLog(File file) throws Exception {
		
		String strline; 
		String[] tokens; 		
		String delims = "\\s+";
		boolean channelValue = false;
        ArrayList<String> navList = new ArrayList<String>();
		
		FileInputStream fin; DataInputStream din; BufferedReader br;
		fin = new FileInputStream(file);
        din = new DataInputStream(fin);
        br = new BufferedReader(new InputStreamReader(din));
		
        //----------------- get header ---
        strline = br.readLine(); 
        System.out.println(strline);
        
        //----------------- get creation time stamp
        strline = br.readLine(); 
        
        tokens = strline.split("\\s+");
        System.out.println(strline + " " + tokens.length);
        
        if(!tokens[3].equals("CREATION")) {
        	br.close();
        	throw new Exception("File not in correct format " + tokens[3]);
        }
        this.creationTimestamp = (new Long(tokens[0])).longValue();
        		
        //----------------- get frequencies
        strline = br.readLine(); 
        System.out.println(strline);
        tokens = strline.split("\\s+");
        if(!tokens[3].equals("FREQUENCY")) {
        	br.close();
        	throw new Exception("File not in correct format " + tokens[3]);
        }
        
        
        Frequencies = new double[3];
        Frequencies[0] = (new Double(tokens[4])).doubleValue();
        Frequencies[1] = (new Double(tokens[5])).doubleValue();
        Frequencies[2] = (new Double(tokens[6])).doubleValue();
        
        //----------------- start measurement 
        strline = br.readLine(); 
        System.out.println(strline);
        //----------------- start measurement 
        strline = br.readLine(); 
        tokens = strline.split("\\s+");
        System.out.println(strline + " " + tokens.length);
        
        if(!tokens[3].equals("SETUP")) {
        	br.close();
        	throw new Exception("File not in correct format " + tokens[3]);
        }
        String[] centerToks = tokens[4].split("[=]+");
        CenterFreq = (new Double(centerToks[1])).doubleValue();
        
        centerToks = tokens[5].split("[=]+");
        Bandwidth = (new Double(centerToks[1])).doubleValue();
        
        //--- now read rest of file to extract navigation and channel logs 
        while((strline = br.readLine()) != null) {        	
        	
        	tokens = strline.split("\\s+");
        	if(!tokens[3].equals("STOP")) {
        		navList.add(strline);        	
        	}
        }
        br.close();         
        
        //--- now read rest of file to extract navigation and channel logs 
        int count = 0;
        while(count < navList.size()) {
        	
          tokens = navList.get(count).split(delims);
          
          if(tokens[3].equals("NAVIGATION")) {
          
        	  long time = (new Long(tokens[0])).longValue();
        	  
        	  String[] toks = tokens[4].split("[=]+");
        	  int seconds = (new Integer(toks[1])).intValue();
             
        	  
              toks = tokens[5].split("[=]+");
              double longitude = (new Double(toks[1])).doubleValue();
 
              toks = tokens[6].split("[=]+");
              double lat = (new Double(toks[1])).doubleValue();              

        	  NavigationChannelList channelList = 
        			  new NavigationChannelList(time, seconds, longitude, lat);
        	  
        	  count++;
        	  channelValue = true;
        	  
        	  while(channelValue) {
        		  
            	  tokens = navList.get(count).split(delims);
            	  if(tokens[3].equals("CHANNEL")) {            		  
            		  
            		  channelList.addChannel(navList.get(count));  
            		  count++;
            	  }
            	  else {
            		  
            		  if(channelList.ChannelList.size() > 0) {
            		   navigationList.add(channelList);
            		  }
            		  channelValue = false;
            	  }
        	  }    
           }
           else {
        	  count++;
           }
        }
	 }

	
     public void createNavigationLog_Interpolation(File file) throws Exception {
		
		String strline; 
		String[] tokens; 		
		String delims = "\\s+";
		int channel_count = 0;
		
		boolean channelValue = false;
        ArrayList<String> navList = new ArrayList<String>();
		
		FileInputStream fin; DataInputStream din; BufferedReader br;
		fin = new FileInputStream(file);
        din = new DataInputStream(fin);
        br = new BufferedReader(new InputStreamReader(din));
		
        //----------------- get header ---
        strline = br.readLine(); 
        System.out.println(strline);
        
        //----------------- get creation time stamp
        strline = br.readLine(); 
        
        tokens = strline.split("\\s+");
        System.out.println(strline + " " + tokens.length);
        
        if(!tokens[3].equals("CREATION")) {
        	br.close();
        	throw new Exception("File not in correct format " + tokens[3]);
        }
        this.creationTimestamp = (new Long(tokens[0])).longValue();
        		
        //----------------- get start measurement
        strline = br.readLine(); 
        System.out.println(strline);

      
        
        strline = br.readLine(); 
        tokens = strline.split("\\s+");
        System.out.println(strline + " " + tokens.length);
        
        if(!tokens[3].equals("SETUP")) {
        	br.close();
        	throw new Exception("File not in correct format " + tokens[3]);
        }
        String[] centerToks = tokens[4].split("[=]+");
        CenterFreq = (new Double(centerToks[1])).doubleValue();
        
        centerToks = tokens[5].split("[=]+");
        Bandwidth = (new Double(centerToks[1])).doubleValue();
        
        //--- now read rest of file to extract navigation and channel logs 
        while((strline = br.readLine()) != null) {        	
        	
        	tokens = strline.split("\\s+");
        	if(!tokens[3].equals("STOP")) {
        		navList.add(strline);        	
        	}
        }
        br.close();         
        
        //--- now read rest of file to extract navigation and channel logs 
        int count = 0;
        while(count < navList.size()) {
        	
          tokens = navList.get(count).split(delims);
          
          if(tokens[3].equals("NAVIGATION")) {
          
        	  long time = (new Long(tokens[0])).longValue();
        	  
        	  String[] toks = tokens[4].split("[=]+");
        	  int seconds = (new Integer(toks[1])).intValue();
             
        	  
              toks = tokens[5].split("[=]+");
              double longitude = (new Double(toks[1])).doubleValue();
 
              toks = tokens[6].split("[=]+");
              double lat = (new Double(toks[1])).doubleValue();              

        	  NavigationChannelList channelList = 
        			  new NavigationChannelList(time, seconds, longitude, lat);
        	  
        	  count++;
        	  channelValue = true;
        	  channel_count = 0;
        	  while(channelValue && count < navList.size()) {
        		  
            	  tokens = navList.get(count).split(delims);
            	  if(tokens[3].equals("CHANNEL")) {            		  
            		  
            		  channelList.addChannel_interpolation(navList.get(count), channel_count);  
            		  count++;
            		  channel_count++;
            	  }
            	  else {
            		  
            		  if(channelList.ChannelList.size() > 0) {
            		   navigationList.add(channelList);
            		  }
            		  channelValue = false;
            	  }
        	  }    
           }
           else {
        	  count++;
           }
        }
	 }
	
	
	

	
	
	
	
	
	

	
	
	public static void main(String[] args) throws Exception {
		
		NavigationList navigation = new NavigationList();
		String sourceFile = "data/mySource.kml";
		
		String[] anchorFiles = new String[15];
		anchorFiles[0] = "zigzagOne.kml";
		navigation.testHybridPath(anchorFiles[0], sourceFile);
		
//		anchorFiles[0] = "straightOne.kml";
//		anchorFiles[1] = "straightTwo.kml";
//		anchorFiles[2] = "straightThree.kml";
//		anchorFiles[3] = "zigzagOne.kml";
//		anchorFiles[4] = "zigzagTwo.kml";
//		anchorFiles[5] = "zigzagThree.kml";		
//		anchorFiles[6] = "incirOne.kml";
//		anchorFiles[7] = "incirTwo.kml";
//		anchorFiles[8] = "incirThree.kml";
//		anchorFiles[9] = "outcirOne.kml";
//		anchorFiles[10] = "outcirTwo.kml";
//		anchorFiles[11] = "outcirThree.kml";
//		anchorFiles[12] = "linesOne.kml";
//		anchorFiles[13] = "linesTwo.kml";
//		anchorFiles[14] = "linesThree.kml";		
//		
//		for(int i = 0; i < 15; i++) {
//			navigation.testFDOAPath(anchorFiles[i], sourceFile);
//		}

	}
	
	
 	/**
 	 * Creates a .kml file to be used in Google Earth that will plot with placemarks all
 	 * the navigation points in the given NavigationChannelList.
 	 *  
 	 * @param navList  The NavigationChannelList to plot
 	 * @return Placemark handle on the placemark 
 	 */	 
     
	public static Placemark createCIRPlaceMark(NavigationChannelList navList) {
		
		Placemark placemark = KmlFactory.createPlacemark();
		placemark.setName("CIR measure at " + navList.getTimeStamp());
		placemark.setVisibility(true);
		placemark.setOpen(false);
		placemark.setDescription(navList.getDescription());
		placemark.setStyleUrl("styles.kml#exampleBalloonStyle");
		
		Point point = KmlFactory.createPoint();
		point.setExtrude(false);
		point.setAltitudeMode(AltitudeMode.CLAMP_TO_GROUND);
		point.getCoordinates().add(new Coordinate(navList.getPointCoordinates()));
				
		placemark.setGeometry(point);

		return placemark; 
	}

 	/**
 	 * Creates a .kml file to be used in Google Earth that will plot with placemarks all
 	 * the navigation points in the given ArrayList of NavigationChannelList.
 	 *  
 	 * @param navigationList The NavigationChannelList to plot
 	 * @return Kml the handle on the Kml file 
 	 */	
	
	public static Kml createCIRDocument(ArrayList<NavigationChannelList> navigationList) {
		
		Kml kml = KmlFactory.createKml();
		Document document = kml.createAndSetDocument().withName("CIRMarkers.kml");
		
		for (NavigationChannelList navList : navigationList){
            document.createAndAddPlacemark()
            //.withName("CIR measure at " + navList.getTimeStamp())
            .withDescription(navList.getDescription())
            .withVisibility(true)
            .createAndSetPoint().addToCoordinates(navList.getLongitude(), navList.getLatitude());
        }
        kml.setFeature(document);
		
		return kml; 
	}
	
	public static Kml createCIRDocumentArrayList(ArrayList<double[]> navigationList) {
		
		Kml kml = KmlFactory.createKml();
		Document document = kml.createAndSetDocument().withName("CIRMarkers.kml");
		
		for (double[] navList : navigationList){
            document.createAndAddPlacemark()
            .withVisibility(true)
            .createAndSetPoint().addToCoordinates(navList[1], navList[0]);
        }
        kml.setFeature(document);		
		return kml; 
	}
	
 	/**
 	 * Creates a .kml file to be used in Google Earth that will plot with placemarks all
 	 * the source estimates and label these source estimates with their respective errors 
 	 * in the given ArrayList. The placemarks will be colored according to their proximity 
 	 * to the source given by the error terms.
 	 *  
 	 * @param estimates A collection of source estimates
 	 * @param errors Their respective errors
 	 * @return Kml the handle on the Kml file 
 	 */	
	
    public static Kml createSolutionDocument(ArrayList<double[]> estimates, ArrayList<Double> errors) {
		
		Kml kml = KmlFactory.createKml();
		Document document = kml.createAndSetDocument().withName("SolutionMarkers.kml");

		for(int j = 0; j < estimates.size() - 1; j++) {
		  
		  double[] myEst = estimates.get(j);	
			
		  document.createAndAddPlacemark()
          .withDescription("Estimate No" + j + ": " + errors.get(j))
          .withVisibility(true)
          .createAndSetPoint().addToCoordinates(myEst[1], myEst[0]);
		  		  
		}
		
		double[] myEst = estimates.get(estimates.size() - 1);	
		document.createAndAddPlacemark()
        .withDescription("Final averaged estimate: " + errors.get(estimates.size() - 1))
        .withVisibility(true)
        .createAndSetPoint().addToCoordinates(myEst[1], myEst[0]);
		
		
		System.out.println("Number of features: " + document.getFeature().size());
		
		List<Color> colors = pick(document.getFeature().size()+2);
		
		for (int i = 0; i < document.getFeature().size()-1; i++) {
			
			Color color = colors.get(i);
			Style style = document.getFeature().get(i).createAndAddStyle();			
			String hex = String.format("%02x%02x%02x", color.getRed(), color.getGreen(), color.getBlue());

			style.createAndSetIconStyle().withColor("e5"+hex).withScale(1.0);
			
			//style.setIconStyle((new IconStyle()).withColor(hex).withHeading(1.0).withScale(1.0));
			document.getFeature().get(i).getStyleSelector().add(style);
		}
		
		Style style = document.getFeature().get(document.getFeature().size() - 1).createAndAddStyle();			
		String hex = String.format("%02x%02x%02x", 0, 0, 0);
		style.createAndSetIconStyle().withColor("e5"+hex).withScale(1.0);
		document.getFeature().get(document.getFeature().size() - 1).getStyleSelector().add(style);
		
		
        kml.setFeature(document);
		
		return kml; 		
	}
	
 	/**
 	 * Creates a .kml file to be used in Google Earth that will plot with placemarks all
 	 * the source estimates and label these source estimates with their respective errors 
 	 * in the given ArrayList. The placemarks will be colored according to their proximity 
 	 * to the source given by the error terms. These estimates are derived from a segmentation 
 	 * solution, so the segments will be colored as well.
 	 *  
 	 * @param coordinates A collection of anchors node coordinates
 	 * @param estimate Their respective estimates
 	 * @param n Number [0, number of total segments-1] of the segment being computed
 	 * @return Kml the handle on the Kml file 
 	 */	
    public static Kml createNavigationAndSolutionDocument(ArrayList<double[]> coordinates, double[] estimate, int n) {
		
		Kml kml = KmlFactory.createKml();
		Document document = kml.createAndSetDocument().withName("SolutionMarkers_" + n + ".kml");
		
		for(int j = 0; j < coordinates.size(); j++) {
		  
		  double[] myEst = coordinates.get(j);	
			
		  document.createAndAddPlacemark()
          //.withName("Error " + errors.get(j))
          .withDescription("Navigation Loc" + j)
          .withVisibility(true)
          .createAndSetPoint().addToCoordinates(myEst[0], myEst[1]);		  		  
		}
		document.createAndAddPlacemark()
        .withName("Estimate " + n)
        .withDescription("Estimate Loc" + n)
        .withVisibility(true)
        .createAndSetPoint().addToCoordinates(estimate[1], estimate[0]);	
				
		System.out.println("Number of features: " + document.getFeature().size());
		
		List<Color> colors = pick(document.getFeature().size()+2);
		Color color = colors.get((int)(colors.size()/(n+1)) - 1);
		
		for (int i = 0; i < document.getFeature().size(); i++) {
						
			Style style = document.getFeature().get(i).createAndAddStyle();			
			String hex = String.format("%02x%02x%02x", color.getRed(), color.getGreen(), color.getBlue());

			style.createAndSetIconStyle().withColor("e5"+hex).withScale(1.0);
			document.getFeature().get(i).getStyleSelector().add(style);			
		}
		
        kml.setFeature(document);
		
		return kml; 		
	}
    
    
    public void parseKml(String src) {

        Kml kml = Kml.unmarshal(new File(src));
        Feature feature = kml.getFeature();
        parseFeature(feature);
    }

    public ArrayList<double[]> parseKmlAnchors(String src) {
    	
    	final Kml kml = Kml.unmarshal(new File(src));
    	final Document document = (Document) kml.getFeature();
    	List<Feature> featureList = document.getFeature();   	
    	
    	ArrayList<double[]> myListAnchors = new ArrayList<double[]>();
    
    	if(featureList.get(0) instanceof Placemark) {

    		final Placemark placemark = (Placemark) featureList.get(0);

    		LineString linestring = (LineString) placemark.getGeometry();
    		List<Coordinate> coordinates = linestring.getCoordinates();
    		
    		for (Coordinate coordinate : coordinates) {
        	    
                double[] locs = {coordinate.getLatitude(), coordinate.getLongitude()};
                myListAnchors.add(locs);
    		}	
    	}   
    	
    	return myListAnchors; 
    }
    
    public double[] parseKmlSource(String src) {
    	
    	double[] mySource = new double[2];
    	final Kml kml = Kml.unmarshal(new File(src));
    	final Document document = (Document) kml.getFeature();
    	List<Feature> featureList = document.getFeature();
    	
    	if(featureList.get(0) instanceof Placemark) {
    		
    		final Placemark placemark = (Placemark) featureList.get(0);
    		Point point = (Point) placemark.getGeometry();
    		List<Coordinate> coordinates = point.getCoordinates();
    		mySource[0] = coordinates.get(0).getLatitude();
    		mySource[1] = coordinates.get(0).getLongitude();
    	}
    	
    	return mySource;
    }
    

    public void testTDOAPath(String anchors, String sourceFile) throws Exception {
    	
    	Random random = new Random(); 
        int num_anchors;
 	    
        double minLong = 400;
 	    double maxLong = -400.0;
 	    double minLat = 400.0;
 	    double maxLat = -400.0;
        
        System.out.println("\nFile: " + anchors);
    	ArrayList<double[]> myListAnchors = parseKmlAnchors("data/" + anchors); 
    	double[] mySource = parseKmlSource(sourceFile);
    	
    	
        final Kml kml0 = createCIRDocumentArrayList(myListAnchors);
		kml0.marshal(new File("FilteredCIRMarkers_" + anchors));
        
		Anchors myAnchors = new Anchors();	    
	    double[] origin = myListAnchors.get(0);	    
	    localOrigin = NavigationChannelList.GeodeticToECEF(origin[0], origin[1], 0);

		double[] source = NavigationChannelList.GeodeticToECEF(mySource[0], mySource[1], 0);
		source[0] -= localOrigin[0];
		source[1] -= localOrigin[1];
		source[2] -= localOrigin[2];
	   		
		
		Matrix sourceMat = new Matrix(source,1);
		
		for(int i = 0; i < myListAnchors.size(); i++) {
											
			double latitude = myListAnchors.get(i)[0];
			double longitude = myListAnchors.get(i)[1];
			
			double[] locs = NavigationChannelList.GeodeticToECEF(latitude, longitude, 0);
			
			myAnchors.setCoordinates(locs[0] - localOrigin[0], 
					                 locs[1] - localOrigin[1], 
					                 locs[2] - localOrigin[2]);
			
			if(latitude > maxLat) {maxLat = latitude;}
			else if(latitude < minLat) {minLat = latitude;}
			
			if(longitude > maxLong) {maxLong = longitude;}
			else if(longitude < minLong) {minLong = longitude;}	
			
		}
	
		double[] boundNW = new double[3];
		double[] boundNE = new double[3];
		double[] boundSW = new double[3];
		double[] boundSE = new double[3];
		   
		boundNW[0] = maxLat + .02; boundNW[1] = minLong - .02;
		boundNE[0] = maxLat + .02; boundNE[1] = maxLong + .02;
		boundSW[0] = minLat - .02; boundSW[1] = minLong - .02;
		boundSE[0] = minLat - .02; boundSE[1] = maxLong + .02;
		   
	    double[] v0 = NavigationChannelList.GeodeticToECEF(boundNW[0], boundNW[1], 0, localOrigin);
	    double[] v1 = NavigationChannelList.GeodeticToECEF(boundNE[0], boundNE[1], 0, localOrigin);
	    double[] v2 = NavigationChannelList.GeodeticToECEF(boundSW[0], boundSW[1], 0, localOrigin);
	    double[] v3 = NavigationChannelList.GeodeticToECEF(boundSE[0], boundSE[1], 0, localOrigin);


	    bounds_in = new Matrix(4,3);
	       
	    bounds_in.setRow(0, v0);
	    bounds_in.setRow(1, v1);
	    bounds_in.setRow(2, v2);
	    bounds_in.setRow(3, v3); 		
		
		
		myAnchors.commitCoordinates();
		num_anchors = myAnchors.getNumberOfAnchors();
		Matrix anchors_in = myAnchors.getAnchors(); 
		Matrix ranges_with_error = new Matrix(num_anchors);
	    
		ranges_with_error.w[0] = 0;
		for(int j = 0; j < num_anchors-1; j++) {
		
			double tdoa_est = GradDescentResult.tdoaEstimate(anchors_in.getRow(0), 
					anchors_in.getRow(j+1), sourceMat);
			
			ranges_with_error.w[j+1] = tdoa_est + noisePercent*tdoa_est*random.nextGaussian(); 
		}
		
	
       ArrayList<double[]> estimates = new ArrayList<double[]>();
       ArrayList<Double> error_est = new ArrayList<Double>();
       
       
       for(int i = 10; i < num_anchors; i=i+20) {

    	Anchors updateAnchors = myAnchors.subset(i);
    	Matrix updateRanges = ranges_with_error.subset(i);
        
        GradDescentResult gdescent_result = GradDescentResult.mlatTdoa(updateAnchors, updateRanges, bounds_in, 
        		n_trial, alpha, time_threshold, source);

        double distance = Mstat.distance(gdescent_result.estimator, sourceMat);
        
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
                
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        
        System.out.print("Error with " + i + " anchor navigation nodes: " + distance + ", estimate: ");
        gdescent_result.estimator.printMatrix();
        
       } 
       
       
       double avgLat = 0; 
       double avgLong = 0; 
       for(int i = 0; i < numberEstimatesAvg ; i++) {
    	   avgLat  += estimates.get(estimates.size() - 1 - i)[0];
    	   avgLong += estimates.get(estimates.size() - 1 - i)[1];
       }
       avgLat = avgLat/numberEstimatesAvg; 
       avgLong = avgLong/numberEstimatesAvg; 
       
       System.out.println(avgLat + ", " + avgLong + ", source = " + mySource[0] + ", " + mySource[1]);
       double finalError = Math.abs(avgLat - mySource[0])*Math.abs(avgLat - mySource[0]) 
    		   + Math.abs(avgLong - mySource[1])*Math.abs(avgLong - mySource[1]);
       
       double[] finalEstimate = new double[3];
       finalEstimate[0] = avgLat; 
       finalEstimate[1] = avgLong; 
       
       finalError = Math.sqrt(finalError); 
       
       System.out.println("Final average error: " + finalError); 
       estimates.add(finalEstimate);
       error_est.add(finalError);
       
       double[] stockArr = new double[ranges_with_error.w.length];
       double[] x = new double[ranges_with_error.w.length];
       for(int i = 0; i < ranges_with_error.w.length; i++) {
    	   stockArr[i] = ranges_with_error.w[i];
    	   x[i] = i;
       }
       
       Plot2DPanel plot = new Plot2DPanel();
		 
       // add a line plot to the PlotPanel
       plot.addLinePlot("TDOA plot", x, stockArr);
       JFrame frame = new JFrame("TDOA simulation");
       frame.setSize(900, 700);
       frame.setContentPane(plot);
       frame.setVisible(true);

       final Kml kml = createSolutionDocument(estimates, error_est);
       kml.marshal(new File("SolutionMarkersTDOA_" + anchors));
    }
    
    
    
	/**
	 * Tests an FDOA path drawn in Google Earth.  
	 * <p>
	 * The measurement will be discarded of the final navigation list 
	 * if either the frequency index does not exist or if the RSSI measurement
	 * is less than the given threshhold. 
	 *
	 * @param  anchors  A string that gives the file name of the kml file drawn in 
	 * google earth that points to the simulated navigation anchors
	 * @param  sourceFile A string that gives the file name of the kml file pointing 
	 * to the source
	 * @throws Exception if dimensions to Matrix are wrong
	 */
    
    
    public void testFDOAPath(String anchors, String sourceFile) throws Exception {
    	
    	Random random = new Random(); 
        int num_anchors;
 	    
        double minLong = 400;
 	    double maxLong = -400.0;
 	    double minLat = 400.0;
 	    double maxLat = -400.0;
        double timeDiff = 1.0;
        double dThreshold = 3.0;
        
        System.out.println("\nFile: " + anchors);
    	ArrayList<double[]> myListAnchors = parseKmlAnchors("data/" + anchors);
    	double[] mySource = parseKmlSource(sourceFile);
    	
    	
        final Kml kml0 = createCIRDocumentArrayList(myListAnchors);
        kml0.marshal(new File("FilteredCIRMarkers_" + anchors));
        
		Anchors myAnchors = new Anchors();	    
	    double[] origin = myListAnchors.get(0);	    
	    localOrigin = NavigationChannelList.GeodeticToECEF(origin[0], origin[1], 0);

		double[] source = NavigationChannelList.GeodeticToECEF(mySource[0], mySource[1], 0);
		source[0] -= localOrigin[0];
		source[1] -= localOrigin[1];
		source[2] -= localOrigin[2];
	   				
		Matrix sourceMat = new Matrix(source,1);
		double[] prevLocation = new double[3]; 
		double[] velocity = new double[3];
		
		
		myAnchors.setCoordinatesAndVelocity(prevLocation, velocity); 
		
		for(int i = 1; i < myListAnchors.size(); i++) {
											
			double latitude = myListAnchors.get(i)[0];
			double longitude = myListAnchors.get(i)[1];
			
			double[] locs = NavigationChannelList.GeodeticToECEF(latitude, longitude, 0);
			double[] anchor = {locs[0] - localOrigin[0], 
	                           locs[1] - localOrigin[1], 
	                           locs[2] - localOrigin[2]};
			
			velocity = new double[3];
			
			timeDiff = Mstat.distance(anchor, prevLocation)/navigationSpeedMetersPerSec;
		
			
			for(int k = 0; k < 3; k++) {
				velocity[k] = (anchor[k] - prevLocation[k])/timeDiff;			
			}
			
			myAnchors.setCoordinatesAndVelocity(anchor, velocity); 
			   
			if(latitude > maxLat) {maxLat = latitude;}
			else if(latitude < minLat) {minLat = latitude;}
			
			if(longitude > maxLong) {maxLong = longitude;}
			else if(longitude < minLong) {minLong = longitude;}	
			
			prevLocation = anchor; 
			
		}
		System.out.println("Origin: " + localOrigin[0] + " " + localOrigin[1] + " " + localOrigin[2]);
		System.out.println("Source: " + source[0] + " " + source[1] + " " + source[2]);
		
		double[] boundNW = new double[3];
		double[] boundNE = new double[3];
		double[] boundSW = new double[3];
		double[] boundSE = new double[3];
		   
		boundNW[0] = maxLat + .02; boundNW[1] = minLong - .02;
		boundNE[0] = maxLat + .02; boundNE[1] = maxLong + .02;
		boundSW[0] = minLat - .02; boundSW[1] = minLong - .02;
		boundSE[0] = minLat - .02; boundSE[1] = maxLong + .02;
		   
	    double[] v0 = NavigationChannelList.GeodeticToECEF(boundNW[0], boundNW[1], 0, localOrigin);
	    double[] v1 = NavigationChannelList.GeodeticToECEF(boundNE[0], boundNE[1], 0, localOrigin);
	    double[] v2 = NavigationChannelList.GeodeticToECEF(boundSW[0], boundSW[1], 0, localOrigin);
	    double[] v3 = NavigationChannelList.GeodeticToECEF(boundSE[0], boundSE[1], 0, localOrigin);

	    System.out.println(v0[0] + " " + v0[1] + " " + v0[2]);
	    System.out.println(v1[0] + " " + v1[1] + " " + v1[2]);
	    System.out.println(v2[0] + " " + v2[1] + " " + v2[2]);
	    System.out.println(v3[0] + " " + v3[1] + " " + v3[2]);
	       
	    bounds_in = new Matrix(4,3);
	       
	    bounds_in.setRow(0, v0);
	    bounds_in.setRow(1, v1);
	    bounds_in.setRow(2, v2);
	    bounds_in.setRow(3, v3); 		
		
		
		myAnchors.commitCoordinatesAndVelocity();
		num_anchors = myAnchors.getNumberOfAnchors();
		Matrix anchors_in = myAnchors.getAnchors(); 
		Matrix velocities = myAnchors.getVelocities(); 
		Matrix ranges_with_error = new Matrix(num_anchors);
	    
		ranges_with_error.w[0] = 0;
		for(int j = 1; j < num_anchors; j++) {
		
			double fdoa_est = GradDescentResult.fdoaEstimate0(anchors_in.getRow(j), 
        			velocities.getRow(j), sourceMat);
			
			double d1 = fdoa_est - ranges_with_error.w[j-1];
			ranges_with_error.w[j] = fdoa_est; 
			
			if(Math.abs(d1) < dThreshold) {
				ranges_with_error.w[j] += fdoa_est*noisePercent*random.nextGaussian(); 
			}
		}
		
	
       ArrayList<double[]> estimates = new ArrayList<double[]>();
       ArrayList<Double> error_est = new ArrayList<Double>();
       
       
       for(int i = 10; i < num_anchors; i=i+20) {

    	Matrix updateAnchors = myAnchors.subsetCoordinates(i);
       	Matrix updateVelocities = myAnchors.subsetVelocity(i);
       	Matrix updateRanges = ranges_with_error.subset(i);
           
       	
       	GradDescentResult gdescent_result = GradDescentResult.mlatFdoa(updateAnchors, 
       			updateVelocities, updateRanges, bounds_in, n_trial, alpha, time_threshold, source);
       	

        double distance = Mstat.distance(gdescent_result.estimator, sourceMat);
        
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
                
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        
        System.out.print("Error with " + i + " anchor navigation nodes: " + distance + ", estimate: ");
        gdescent_result.estimator.printMatrix();
        
       } 
       
       int numberEstimatesAvg = 5;
       double avgLat = 0; 
       double avgLong = 0; 
       for(int i = 0; i < numberEstimatesAvg; i++) {
    	   avgLat  += estimates.get(estimates.size() - 1 - i)[0];
    	   avgLong += estimates.get(estimates.size() - 1 - i)[1];
       }
       avgLat = avgLat/numberEstimatesAvg; 
       avgLong = avgLong/numberEstimatesAvg; 
       
       System.out.println(avgLat + ", " + avgLong + ", source = " + mySource[0] + ", " + mySource[1]);
       double finalError = Math.abs(avgLat - mySource[0])*Math.abs(avgLat - mySource[0]) 
    		   + Math.abs(avgLong - mySource[1])*Math.abs(avgLong - mySource[1]);
       
       double[] finalEstimate = new double[3];
       finalEstimate[0] = avgLat; 
       finalEstimate[1] = avgLong; 
       
       finalError = Math.sqrt(finalError); 
       
       System.out.println("Final average error: " + finalError); 
       estimates.add(finalEstimate);
       error_est.add(finalError);
       
       double[] stockArr = new double[ranges_with_error.w.length];
       double[] x = new double[ranges_with_error.w.length];
       for(int i = 0; i < ranges_with_error.w.length; i++) {
    	   stockArr[i] = ranges_with_error.w[i];
    	   x[i] = i;
       }
       
       Plot2DPanel plot = new Plot2DPanel();
		 
       // add a line plot to the PlotPanel
       plot.addLinePlot("FDOA plot", x, stockArr);
       JFrame frame = new JFrame("FDOA simulation");
       frame.setSize(900, 700);
       frame.setContentPane(plot);
       frame.setVisible(true);

        
       final Kml kml = createSolutionDocument(estimates, error_est);
       kml.marshal(new File("SolutionMarkersFDOA_" + anchors));
		
    	
    }
    
    
    
    /**
	 * Tests an FDOA+TDOA path drawn in Google Earth.  
	 * <p>
	 *
	 * @param  anchors  A string that gives the file name of the kml file drawn in 
	 * google earth that points to the simulated navigation anchors
	 * @param  sourceFile A string that gives the file name of the kml file pointing 
	 * to the source
	 * @throws Exception if dimensions to Matrix are wrong
	 */
    
    
    public void testHybridPath(String anchors, String sourceFile) throws Exception {
    	
    	Random random = new Random(); 
        int num_anchors;
 	    
        double minLong = 400;
 	    double maxLong = -400.0;
 	    double minLat = 400.0;
 	    double maxLat = -400.0;
        double timeDiff = 1.0;
        double dThreshold = 3.0;
        
        System.out.println("\nFile: " + anchors);
    	ArrayList<double[]> myListAnchors = parseKmlAnchors("data/" + anchors);
    	double[] mySource = parseKmlSource(sourceFile);
    	
    	
        final Kml kml0 = createCIRDocumentArrayList(myListAnchors);
        kml0.marshal(new File("FilteredCIRMarkers_" + anchors));
        
		Anchors myAnchors = new Anchors();	    
	    double[] origin = myListAnchors.get(0);	    
	    localOrigin = NavigationChannelList.GeodeticToECEF(origin[0], origin[1], 0);

		double[] source = NavigationChannelList.GeodeticToECEF(mySource[0], mySource[1], 0);
		source[0] -= localOrigin[0];
		source[1] -= localOrigin[1];
		source[2] -= localOrigin[2];
	   				
		Matrix sourceMat = new Matrix(source,1);
		double[] prevLocation = new double[3]; 
		double[] velocity = new double[3];
		
		
		myAnchors.setCoordinatesAndVelocity(prevLocation, velocity); 
		
		for(int i = 1; i < myListAnchors.size(); i++) {
											
			double latitude = myListAnchors.get(i)[0];
			double longitude = myListAnchors.get(i)[1];
			
			double[] locs = NavigationChannelList.GeodeticToECEF(latitude, longitude, 0);
			double[] anchor = {locs[0] - localOrigin[0], 
	                           locs[1] - localOrigin[1], 
	                           locs[2] - localOrigin[2]};
			
			velocity = new double[3];
			
			timeDiff = Mstat.distance(anchor, prevLocation)/navigationSpeedMetersPerSec;
		
			
			for(int k = 0; k < 3; k++) {
				velocity[k] = (anchor[k] - prevLocation[k])/timeDiff;			
			}
			
			myAnchors.setCoordinatesAndVelocity(anchor, velocity); 
			   
			if(latitude > maxLat) {maxLat = latitude;}
			else if(latitude < minLat) {minLat = latitude;}
			
			if(longitude > maxLong) {maxLong = longitude;}
			else if(longitude < minLong) {minLong = longitude;}	
			
			prevLocation = anchor; 
			
		}
		System.out.println("Origin: " + localOrigin[0] + " " + localOrigin[1] + " " + localOrigin[2]);
		System.out.println("Source: " + source[0] + " " + source[1] + " " + source[2]);
		
		double[] boundNW = new double[3];
		double[] boundNE = new double[3];
		double[] boundSW = new double[3];
		double[] boundSE = new double[3];
		   
		boundNW[0] = maxLat + .02; boundNW[1] = minLong - .02;
		boundNE[0] = maxLat + .02; boundNE[1] = maxLong + .02;
		boundSW[0] = minLat - .02; boundSW[1] = minLong - .02;
		boundSE[0] = minLat - .02; boundSE[1] = maxLong + .02;
		   
	    double[] v0 = NavigationChannelList.GeodeticToECEF(boundNW[0], boundNW[1], 0, localOrigin);
	    double[] v1 = NavigationChannelList.GeodeticToECEF(boundNE[0], boundNE[1], 0, localOrigin);
	    double[] v2 = NavigationChannelList.GeodeticToECEF(boundSW[0], boundSW[1], 0, localOrigin);
	    double[] v3 = NavigationChannelList.GeodeticToECEF(boundSE[0], boundSE[1], 0, localOrigin);

	    System.out.println(v0[0] + " " + v0[1] + " " + v0[2]);
	    System.out.println(v1[0] + " " + v1[1] + " " + v1[2]);
	    System.out.println(v2[0] + " " + v2[1] + " " + v2[2]);
	    System.out.println(v3[0] + " " + v3[1] + " " + v3[2]);
	       
	    bounds_in = new Matrix(4,3);
	       
	    bounds_in.setRow(0, v0);
	    bounds_in.setRow(1, v1);
	    bounds_in.setRow(2, v2);
	    bounds_in.setRow(3, v3); 		
		
		
		myAnchors.commitCoordinatesAndVelocity();
		num_anchors = myAnchors.getNumberOfAnchors();
		Matrix anchors_in = myAnchors.getAnchors(); 
		Matrix velocities = myAnchors.getVelocities(); 
		Matrix ranges_with_error = new Matrix(num_anchors);
	    
		ranges_with_error.w[0] = 0;
		for(int j = 1; j < num_anchors; j++) {
		
			double fdoa_est = GradDescentResult.fdoaEstimate0(anchors_in.getRow(j), 
        			velocities.getRow(j), sourceMat);
			
			double d1 = fdoa_est - ranges_with_error.w[j-1];
			ranges_with_error.w[j] = fdoa_est; 
			
			if(Math.abs(d1) < dThreshold) {
				ranges_with_error.w[j] += fdoa_est*noisePercent*random.nextGaussian(); 
			}
		}
		
		Matrix tdoas_with_error = new Matrix(num_anchors);
		ranges_with_error.w[0] = 0;
		for(int j = 0; j < num_anchors-1; j++) {
		
			double tdoa_est = GradDescentResult.tdoaEstimate(anchors_in.getRow(0), 
					anchors_in.getRow(j+1), sourceMat);
			
			tdoas_with_error.w[j+1] = tdoa_est + noisePercent*tdoa_est*random.nextGaussian(); 
		}
		
	
       ArrayList<double[]> estimates = new ArrayList<double[]>();
       ArrayList<Double> error_est = new ArrayList<Double>();
       
       
       for(int i = 10; i < num_anchors; i=i+20) {

    	Matrix updateAnchors = myAnchors.subsetCoordinates(i);
       	Matrix updateVelocities = myAnchors.subsetVelocity(i);
       	Matrix updateRanges = ranges_with_error.subset(i);
       	Matrix updateTdoa = tdoas_with_error.subset(i);
           
       	
       	GradDescentResult gdescent_result = GradDescentResult.mlatHybrid(updateAnchors, 
       			updateVelocities, updateRanges, updateTdoa, bounds_in, n_trial, alpha, time_threshold, source);
       	

        double distance = Mstat.distance(gdescent_result.estimator, sourceMat);
        
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
                
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        
        System.out.print("Error with " + i + " anchor navigation nodes: " + distance + ", estimate: ");
        gdescent_result.estimator.printMatrix();
        
       } 
       
       int numberEstimatesAvg = 5;
       double avgLat = 0; 
       double avgLong = 0; 
       for(int i = 0; i < numberEstimatesAvg; i++) {
    	   avgLat  += estimates.get(estimates.size() - 1 - i)[0];
    	   avgLong += estimates.get(estimates.size() - 1 - i)[1];
       }
       avgLat = avgLat/numberEstimatesAvg; 
       avgLong = avgLong/numberEstimatesAvg; 
       
       System.out.println(avgLat + ", " + avgLong + ", source = " + mySource[0] + ", " + mySource[1]);
       double finalError = Math.abs(avgLat - mySource[0])*Math.abs(avgLat - mySource[0]) 
    		   + Math.abs(avgLong - mySource[1])*Math.abs(avgLong - mySource[1]);
       
       double[] finalEstimate = new double[3];
       finalEstimate[0] = avgLat; 
       finalEstimate[1] = avgLong; 
       
       finalError = Math.sqrt(finalError); 
       
       System.out.println("Final average error: " + finalError); 
       estimates.add(finalEstimate);
       error_est.add(finalError);
       
       double[] stockArr = new double[ranges_with_error.w.length];
       double[] x = new double[ranges_with_error.w.length];
       for(int i = 0; i < ranges_with_error.w.length; i++) {
    	   stockArr[i] = ranges_with_error.w[i];
    	   x[i] = i;
       }
       
       Plot2DPanel plot = new Plot2DPanel();
		 
       // add a line plot to the PlotPanel
       plot.addLinePlot("FDOA plot", x, stockArr);
       JFrame frame = new JFrame("FDOA simulation");
       frame.setSize(900, 700);
       frame.setContentPane(plot);
       frame.setVisible(true);

        
       final Kml kml = createSolutionDocument(estimates, error_est);
       kml.marshal(new File("SolutionMarkersFDOA_" + anchors));
		
    	
    }
    
    
    
    
   
    
    
    private void parseFeature(Feature feature) {
        
    	
    	if(feature != null) {
            if(feature instanceof Document) {
                Document document = (Document) feature;
                List<Feature> featureList = document.getFeature();
                for(Feature documentFeature : featureList) {
                    if(documentFeature instanceof Placemark) {
                        Placemark placemark = (Placemark) documentFeature;
                        System.out.println(placemark.getDescription());
                        Geometry geometry = placemark.getGeometry();
                        parseGeometry(geometry);
                    }
                }
            }
        }
    }

    private void parseGeometry(Geometry geometry) {
        if(geometry != null) {
            if(geometry instanceof Polygon) {
                Polygon polygon = (Polygon) geometry;
                Boundary outerBoundaryIs = polygon.getOuterBoundaryIs();
                if(outerBoundaryIs != null) {
                    LinearRing linearRing = outerBoundaryIs.getLinearRing();
                    if(linearRing != null) {
                        List<Coordinate> coordinates = linearRing.getCoordinates();
                        if(coordinates != null) {
                            for(Coordinate coordinate : coordinates) {
                                parseCoordinate(coordinate);
                            }
                        }
                    }
                }
            }
        }
    }

    private void parseCoordinate(Coordinate coordinate) {
        if(coordinate != null) {
            System.out.println("Longitude: " +  coordinate.getLongitude());
            System.out.println("Latitude : " +  coordinate.getLatitude());
            System.out.println("Altitude : " +  coordinate.getAltitude());
            System.out.println("");
        }
    }
    
    
    
    
    
    
    public static List<Color> pick(int num) {
        List<Color> colors = new ArrayList<Color>();
        if (num < 2)
            return colors;
        float dx = 1.0f / (float) (num - 1);
        for (int i = 0; i < num; i++) {
            colors.add(get(i * dx));
        }
        return colors;
    }

    public static Color get(float x) {
        float r = 0.0f;
        float g = 0.0f;
        float b = 1.0f;
        if (x >= 0.0f && x < 0.2f) {
            x = x / 0.2f;
            r = 0.0f;
            g = x;
            b = 1.0f;
        } else if (x >= 0.2f && x < 0.4f) {
            x = (x - 0.2f) / 0.2f;
            r = 0.0f;
            g = 1.0f;
            b = 1.0f - x;
        } else if (x >= 0.4f && x < 0.6f) {
            x = (x - 0.4f) / 0.2f;
            r = x;
            g = 1.0f;
            b = 0.0f;
        } else if (x >= 0.6f && x < 0.8f) {
            x = (x - 0.6f) / 0.2f;
            r = 1.0f;
            g = 1.0f - x;
            b = 0.0f;
        } else if (x >= 0.8f && x <= 1.0f) {
            x = (x - 0.8f) / 0.2f;
            r = 1.0f;
            g = 0.0f;
            b = x;
        }
        return new Color(r, g, b);
    }
    
    
    
	
//	public static void main(String[] args) throws Exception {
//		
//		NavigationList navigation = new NavigationList();
//		
//		//------- Choose log file --------------------
//		String file = "data/ChannelLog_ETZIKEN.log"; //"data/ChannelLog.log" "data/ChannelLog_ZWEI.log"
//		navigation.createNavigationLog(new File(file));
//
//		final Kml kml = createCIRDocument(navigation.navigationList);
//		kml.marshal(new File("CIRMarkers.kml"));
//		
//		
//		int n_estimates = 8;
//		int freqIndex = 2; 
//		double threshold_tdoaDiff = .000001;
//		double threshold_dynamicRange = 20.0;
//		double threshold_FDOAerror = 3.0;
//		
//		navigation.computeVelocityECEF();
//		navigation.filterOnRSSI(-90.0, freqIndex);
//		navigation.filterUnique();
//		navigation.createEstimationBounds();		
//		navigation.setPlotTDOAs(true);
//		navigation.filterFDOAfromNavigationList(freqIndex, threshold_FDOAerror, threshold_dynamicRange);
//		navigation.estimateAdaptiveSourceFDOA();
//
//		//navigation.testSGDfdoa();
//        //navigation.filterTDOAfromNavigationList(freqIndex, threshold_tdoaDiff);
//		//navigation.estimateDynamicReferenceTDOA(freqIndex, n_estimates, threshold_tdoaDiff);
//		//navigation.estimateAdaptiveSourceTDOA();
//		//navigation.estimateAdaptiveSourceWithDrift();
//	}
    
    
}
