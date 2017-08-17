package ch.imetrica.jdeeplateration.cirdata;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Random;

import ch.imetrica.jdeeplateration.anchors.Anchors;
import ch.imetrica.jdeeplateration.matrix.Matrix;
import ch.imetrica.jdeeplateration.mstat.GradDescentResult;
import ch.imetrica.jdeeplateration.mstat.Mstat;
import de.micromata.opengis.kml.v_2_2_0.AltitudeMode;
import de.micromata.opengis.kml.v_2_2_0.Coordinate;
import de.micromata.opengis.kml.v_2_2_0.Document;
import de.micromata.opengis.kml.v_2_2_0.Kml;
import de.micromata.opengis.kml.v_2_2_0.KmlFactory;
import de.micromata.opengis.kml.v_2_2_0.Placemark;
import de.micromata.opengis.kml.v_2_2_0.Point;


public class NavigationList {

	double[] Frequencies; 

	double CenterFreq; 
	double Bandwidth;
	
	long creationTimestamp;
	
	ArrayList<NavigationChannelList> navigationList;
	ArrayList<TimeDiffOnArrival> tdoaFrequency0;
	
	
	public NavigationList() {
		navigationList = new ArrayList<NavigationChannelList>();
	}
	
	
	public void computeTDOAfromNavigationList(int freq) {
		
		double refTime = -1.0;
		double residual = 0;
		
		tdoaFrequency0 = new ArrayList<TimeDiffOnArrival>();
		
		for (NavigationChannelList navList : navigationList) {
			
			double time = navList.getTimeStampAtFreq(freq);
			
			if(refTime < 0 && time > 0) { 
				refTime = time;  //set reference time for first measurement
			
				residual = refTime*100.0 - Math.floor(refTime*100.0);
				residual = residual/100.0;
			}
			else if(time > 0) {
				
				double currentRefTime = Math.floor(time*100.0)/100.0 + residual;				
				double tdoa = time - currentRefTime;
				
				TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa);
				td.printTDOA();
				
				tdoaFrequency0.add(td);
			}
		}	
	}
	
	
	public void filter() {
		//Filters and prunes entire raw navigation list
		filterOnRSSI(-70, 0);
	}
	
	
	public void computeVelocity() {
		
		for(int i = 1; i < navigationList.size(); i++) {
		
			NavigationChannelList navList = navigationList.get(i-1);
			long prevTime = navList.getTimeStamp();
			double prevLongitude = navList.getLongitude();
			double prevLatitude = navList.getLatitude();
			
			navigationList.get(i).computeVelocity(prevTime, prevLongitude, prevLatitude);
			navigationList.get(i).printVelocity();
		}
		
	}
	
	public void filterOnRSSI(double threshhold, int freq) {
		
		for(int i = 0; i < navigationList.size(); i++) {
			
			if(navigationList.get(i).getRSSI(freq) < threshhold) {
				navigationList.remove(i);
			}	
		}		
	}
	
	
	
	/*
	 * 
	 * Use stochastic gradient method to find solution given 
	 * set of long/lat and tdoa 
	 * 
	 */
	
	public void estimateSourceSolution() throws Exception {
		
		
        
        int num_anchors;
        
		Anchors myAnchors = new Anchors();
	    ArrayList<Double> tdoas = new ArrayList<Double>();
	    
		
		for(int i = 53; i < tdoaFrequency0.size() - 20; i++) {
						
			TimeDiffOnArrival tdoa = tdoaFrequency0.get(i);
					
			myAnchors.setCoordinates(tdoa.getLongitude(), tdoa.getLatitude(), 0);
			tdoas.add(tdoa.getTDOA());
		}
		
		myAnchors.commitCoordinates();
		num_anchors = myAnchors.getNumberOfAnchors();
		   
        Matrix ranges = new Matrix(num_anchors);
        Matrix ranges_with_error = new Matrix(num_anchors);
        
        for (int i = 0; i < num_anchors; i++) {
        	
            ranges.w[i] = tdoas.get(i).doubleValue();
            ranges_with_error.w[i] = ranges.w[i];
        }

        
        int n_trial = 5000; 
        double alpha = 0.0001; 
        double time_threshold = 50000;
        
        // Set the bounds in which source will be located  7.362370 46.45439
        Matrix bounds_in = new Matrix(2,3);
        bounds_in.set(0, 0, myAnchors.getColumnMin(0) - 1.0);
        bounds_in.set(0, 1, myAnchors.getColumnMin(1) - 1.0);
        
        bounds_in.set(1, 0, myAnchors.getColumnMax(0) + 1.0);
        bounds_in.set(1, 1, myAnchors.getColumnMax(1) + 1.0);
        
        GradDescentResult gdescent_result = GradDescentResult.mlat(myAnchors, ranges_with_error, bounds_in, n_trial, alpha, time_threshold);

        System.out.println("Anchors");
        myAnchors.printMatrix();
        
        System.out.println("Ranges");
        ranges.printMatrix();
        
//        System.out.println("Ranges with error");
//        ranges_with_error.printMatrix();
        
        System.out.println("\nEstimator");
        gdescent_result.estimator.printMatrix();
        
//        System.out.println("\nFull result");
//        gdescent_result.estimator_candidate.printMatrix();
//        gdescent_result.error.printMatrix();
		
        final Kml kml = createSolutionDocument(gdescent_result.estimator_candidate, 
        		                               gdescent_result.error, 
        		                               gdescent_result.estimator);
        
		kml.marshal(new File("SolutionMarkers.kml"));

   
	}
	
	
	public void testCaseData(File file) throws Exception {
		
		//x y z vx vyvz toa(time of arrival) foa(frequency of arrival) absolute time(not really used)
		
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
        	
        	tdoas.add((new Double(tokens[7])).doubleValue());
        
        }
        br.close();
        
        myAnchors.commitCoordinates();
		num_anchors = myAnchors.getNumberOfAnchors();
		   
        Matrix ranges = new Matrix(num_anchors);
        Matrix ranges_with_error = new Matrix(num_anchors);
        
        for (int i = 0; i < num_anchors; i++) {
        	
            ranges.w[i] = tdoas.get(i).doubleValue();
            ranges_with_error.w[i] = ranges.w[i];
        }

        
        int n_trial = 5000; 
        double alpha = 0.00001; 
        double time_threshold = 500000;
        
        // Set the bounds in which source will be located  7.362370 46.45439
        Matrix bounds_in = new Matrix(2,3);
        bounds_in.set(0, 0, myAnchors.getColumnMin(0) - 100.0);
        bounds_in.set(0, 1, myAnchors.getColumnMin(1) - 100.0);
        bounds_in.set(0, 2, myAnchors.getColumnMin(2));
        
        bounds_in.set(1, 0, myAnchors.getColumnMax(0) + 100.0);
        bounds_in.set(1, 1, myAnchors.getColumnMax(1) + 100.0);
        bounds_in.set(1, 2, myAnchors.getColumnMax(2));
        
        GradDescentResult gdescent_result = GradDescentResult.mlatTdoa(myAnchors, ranges_with_error, bounds_in, n_trial, alpha, time_threshold);

        System.out.println("Anchors");
        myAnchors.printMatrix();
        
        System.out.println("Ranges");
        ranges.printMatrix();
        
//        System.out.println("Ranges with error");
//        ranges_with_error.printMatrix();
        
        System.out.println("\nEstimator");
        gdescent_result.estimator.printMatrix();
        
        System.out.println("\nTrue location: " + "-2762.18694234 -5068.57649907  2522.74598415");
        
        
        
	}
	
	
	
	
	public void printNavigationHistory() {
		
		System.out.println("My size: " + navigationList.size());
		
		for(NavigationChannelList nav : navigationList) {
			nav.printCoordinates();
		}
	}
	
	
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
            		  
            		  navigationList.add(channelList);
            		  channelValue = false;
            	  }
        	  }    
           }
           else {
        	  count++;
           }
        }
	 }

	
	public void createLocalCoordinateSystem() {
		
		double longitude = navigationList.get(0).getLongitude();
		double latitude = navigationList.get(0).getLatitude();
		double altitude = 0;
		
		double[] localOrigin = NavigationChannelList.GeodeticToECEF(longitude, latitude, altitude);
		
		for (NavigationChannelList navList : navigationList) {			
			navList.GeodeticToLocal(localOrigin);			
		}
	}
	
	public static void main(String[] args) throws Exception {
		
		NavigationList navigation = new NavigationList();
		navigation.createNavigationLog(new File("data/ChannelLog.log"));
				
		final Kml kml = createCIRDocument(navigation.navigationList);
		kml.marshal(new File("CIRMarkers.kml"));
		
		navigation.computeTDOAfromNavigationList(0);
		navigation.computeVelocity();
		
		//navigation.estimateSourceSolution();
		
		navigation.testCaseData(new File("data/test_data.txt"));
		
	}
	
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
	
	public static Kml createCIRDocument(ArrayList<NavigationChannelList> navigationList) {
		
		Kml kml = KmlFactory.createKml();
		Document document = kml.createAndSetDocument().withName("CIRMarkers.kml");
		
		for (NavigationChannelList navList : navigationList){
            document.createAndAddPlacemark()
            .withName("CIR measure at " + navList.getTimeStamp())
            .withDescription(navList.getDescription())
            .withVisibility(true)
            .createAndSetPoint().addToCoordinates(navList.getLongitude(), navList.getLatitude());
        }
        kml.setFeature(document);
		
		return kml; 
		
	}
	
    public static Kml createSolutionDocument(Matrix locs, Matrix error, Matrix solution) {
		
		Kml kml = KmlFactory.createKml();
		Document document = kml.createAndSetDocument().withName("SolutionMarkers.kml");
				
		document.createAndAddPlacemark()
        .withName("Error " + error.w[error.rows-1])
        .withDescription("Solution")
        .withVisibility(true)
        .createAndSetPoint().addToCoordinates(solution.getW(0, 0), solution.getW(0, 1));
		
        kml.setFeature(document);
		
		return kml; 
		
	}
	
	
}
