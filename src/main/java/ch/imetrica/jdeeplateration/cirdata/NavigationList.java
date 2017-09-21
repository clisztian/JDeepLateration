package ch.imetrica.jdeeplateration.cirdata;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;

import org.math.plot.Plot2DPanel;

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
import de.micromata.opengis.kml.v_2_2_0.Style;


public class NavigationList {

	
	public final static double C = 299792458;
	double[] Frequencies; 
	double[] localOrigin = null;

	double CenterFreq; 
	double Bandwidth;
	
	long creationTimestamp;
	
	ArrayList<NavigationChannelList> navigationList = null;
	ArrayList<NavigationChannelList> filteredNavigationList = null;
	ArrayList<TimeDiffOnArrival> tdoaFrequency0;
	
	Matrix bounds_in;
	private double sourceLatitude = 0;
	private double sourceLongitude = 0;
	private boolean plotTDOAs = false;
	
	public NavigationList() {
		navigationList = new ArrayList<NavigationChannelList>();
	}
	
	public void setPlotTDOAs(boolean f) {
		this.plotTDOAs = f; 
	}
	
	public void setKnownSource(double lat, double longi) {
		this.sourceLatitude = lat; 
		this.sourceLongitude = longi; 
	}
	
	public void computeTDOAfromNavigationList(int freq) throws Exception {
		
	   if(filteredNavigationList == null) {	
		   throw new Exception("Apply filtering to navigation list first");   
	   }
		
		double refTime = -1.0;
		double residual = 0;
		
		tdoaFrequency0 = new ArrayList<TimeDiffOnArrival>();
		
		for (NavigationChannelList navList : filteredNavigationList) {

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
	
	
	
	
	public void computeTDOAfromNavigationList0(int freq) throws Exception {
		
	   if(filteredNavigationList == null) {	
		   throw new Exception("Apply filtering to navigation list first");   
	   }
		
		double pulse = .065536;
		
		tdoaFrequency0 = new ArrayList<TimeDiffOnArrival>();
		
		//for (NavigationChannelList navList : filteredNavigationList) {
      
		double time = filteredNavigationList.get(0).getTimeStampAtFreq(freq);
				
		for(int i = 1; i < filteredNavigationList.size(); i++) {
			
			double currentRefTime = filteredNavigationList.get(i).getTimeStampAtFreq(freq);
						
			if(currentRefTime > 0) {
							
				double tdoa = (currentRefTime - time)%pulse;
				
				if(tdoa > .06) {
					tdoa -= pulse;
				}
				
				TimeDiffOnArrival td = new TimeDiffOnArrival(filteredNavigationList.get(i).getLongitude(),
						filteredNavigationList.get(i).getLatitude(), tdoa);
				td.printTDOA();
				
				tdoaFrequency0.add(td);
			}
		}
	}
	
	
	public void filterPeakList() {
		
		
		
	}
	
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
	                refTime = time;  //set reference time for first measurement
	            
	                tdoa = 0.0;
	                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa);                
	                
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
	                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa);                
	                
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
	
	
	public void filterTDOAfromNavigationList(int freq, double thresh) throws Exception {
        
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
	                refTime = time;  //set reference time for first measurement
	            
	                tdoa = 0.0;
	                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa);                
	                
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
		                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa);                
		                times.add(time);
		                tdoaFrequency0.add(td);   
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
	
	
	public void filterFDOAfromNavigationList(int freq, double threshold) throws Exception {
		
	       if(filteredNavigationList == null) {    
	           throw new Exception("Apply filtering to navigation list first");   
	       }
	       
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
	            	System.out.println("Eliminating " + i);
	            	filteredNavigationList.remove(i);
	            }
	            else {i++;}
	       }	
	       System.out.println("Navigation readings after fdoa filtering: " + filteredNavigationList.size());
	       
	       double[] y = new double[filteredNavigationList.size()];
	       double[] x = new double[filteredNavigationList.size()];
	       
	       for (i = 0; i < filteredNavigationList.size(); i++) {
	    	   System.out.print(filteredNavigationList.get(i).getFDOA()); 
	    	   filteredNavigationList.get(i).printVelocity();
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
	
	
	public void filter() {
		//Filters and prunes entire raw navigation list
		filterOnRSSI(-70, 0);
	}
	
	
	public void computeVelocity() {
		
		navigationList.get(0).setVelocityToZero();
		for(int i = 1; i < navigationList.size(); i++) {
		
			NavigationChannelList navList = navigationList.get(i-1);
			long prevTime = navList.getDaySeconds();
			double[] prevlocalECEF = navList.getLocalECEF();
			
			navigationList.get(i).computeVelocityECEF(prevTime, prevlocalECEF);
			//navigationList.get(i).printVelocity();
		}
	}
	
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
    
    
    public void computeVelocityECEF() {
    	
    	localOrigin = NavigationChannelList.GeodeticToECEF(navigationList.get(0).getLatitude(),
    			navigationList.get(0).getLongitude(), 0);  	
    	for(NavigationChannelList nav : navigationList) {    		
    		nav.GeodeticToLocal(localOrigin);
    	}
    	computeVelocity();   	
    }
    
    
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
    
    
    public void estimateSourceSolutionFDOA() throws Exception {
    	
    	
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

        
        int n_trial = 500; 
        double alpha = 0.0001; 
        double time_threshold = 50000;
        
        
       
        double[] source = NavigationChannelList.GeodeticToECEF(46.762606, 7.600533, 0);
		source[0] -= localOrigin[0];
		source[1] -= localOrigin[1];
		source[2] -= localOrigin[2];
        
        GradDescentResult gdescent_result = GradDescentResult.mlatFdoa(myAnchors.getAnchors(), 
        		myAnchors.getVelocities(), ranges, bounds_in, n_trial, alpha, time_threshold, source);
	    
        System.out.println("\nEstimator");
        gdescent_result.estimator.transformRowCoord(0,localOrigin);         
        gdescent_result.estimator.printMatrix();
        
        //gdescent_result.error.printMatrix();
        
        for(int j = 0; j < gdescent_result.estimator_candidate.rows; j++) {
        	gdescent_result.estimator_candidate.transformRowCoord(j,localOrigin); 
        }
    	
        
        ArrayList<double[]> estimates = new ArrayList<double[]>();
        ArrayList<Double> error_est = new ArrayList<Double>();
        
        //gdescent_result.estimator.transformRowCoord(0,localOrigin);         
        
        estimates.add(gdescent_result.estimator.w);
        error_est.add(gdescent_result.error.w[0]);
        
        
        final Kml kml = createSolutionDocument(estimates, error_est);
        kml.marshal(new File("SolutionMarkers.kml"));
        
        
    }
    
    


	
	public void estimateSourceSolution() throws Exception {
		
		
        Random random = new Random(); 
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
	    

		Matrix sourceMat = new Matrix(source,1);
		
		for(int i = 0; i < tdoaFrequency0.size(); i++) {
						
			TimeDiffOnArrival tdoa = tdoaFrequency0.get(i);
					
			double[] locs = NavigationChannelList.GeodeticToECEF(tdoa.getLatitude(), tdoa.getLongitude(), 0);
			
			myAnchors.setCoordinates(locs[0] - localOrigin[0], 
					                 locs[1] - localOrigin[1], 
					                 locs[2] - localOrigin[2]);
						
			rtdoas.add(tdoa.getTDOA()*C);
			
		}
		System.out.println("Origin: " + localOrigin[0] + " " + localOrigin[1] + " " + localOrigin[2]);
		System.out.println("Source: " + source[0] + " " + source[1] + " " + source[2]);
		
		myAnchors.commitCoordinates();
		num_anchors = myAnchors.getNumberOfAnchors();
		Matrix anchors_in = myAnchors.getAnchors(); 
		
		
		tdoas.add(0.0);
		for(int j = 0; j < num_anchors-1; j++) {
		
			double tdoa_est = GradDescentResult.tdoaEstimate(anchors_in.getRow(j), 
					anchors_in.getRow(j+1), sourceMat);
			
			tdoas.add(tdoa_est); 
		}
				 
		for(int j = 0; j < num_anchors; j++) {
			
			System.out.println(rtdoas.get(j) + " " + tdoas.get(j));
			
		}
		
        Matrix ranges = new Matrix(num_anchors);
        Matrix ranges_with_error = new Matrix(num_anchors);
        
        for (int i = 0; i < num_anchors; i++) {
        	
            ranges.w[i] = tdoas.get(i).doubleValue();
            ranges_with_error.w[i] = ranges.w[i] + 4.0*random.nextGaussian();
        }

        
        int n_trial = 500; 
        double alpha = 0.0001; 
        double time_threshold = 50000;
        
        
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
	
	
	
	public void estimateAdaptiveSource() throws Exception {
		
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
	    

		Matrix sourceMat = new Matrix(source,1);
		
		for(int i = 0; i < tdoaFrequency0.size(); i++) {
						
			TimeDiffOnArrival tdoa = tdoaFrequency0.get(i);
					
			double[] locs = NavigationChannelList.GeodeticToECEF(tdoa.getLatitude(), tdoa.getLongitude(), 0);
			
			myAnchors.setCoordinates(locs[0] - localOrigin[0], 
					                 locs[1] - localOrigin[1], 
					                 locs[2] - localOrigin[2]);
						
			rtdoas.add(tdoa.getTDOA()*C);
			
		}
		System.out.println("Origin: " + localOrigin[0] + " " + localOrigin[1] + " " + localOrigin[2]);
		System.out.println("Source: " + source[0] + " " + source[1] + " " + source[2]);
		
		myAnchors.commitCoordinates();
		num_anchors = myAnchors.getNumberOfAnchors();
		Matrix anchors_in = myAnchors.getAnchors(); 
		
		
		tdoas.add(0.0);
		for(int j = 0; j < num_anchors-1; j++) {
		
			double tdoa_est = GradDescentResult.tdoaEstimate(anchors_in.getRow(0), 
					anchors_in.getRow(j+1), sourceMat);
			
			tdoas.add(tdoa_est); 
		}
			
		
        Matrix ranges = new Matrix(num_anchors);
        Matrix ranges_with_error = new Matrix(num_anchors);
        
        for (int i = 0; i < num_anchors; i++) {
        	
            ranges.w[i] = rtdoas.get(i).doubleValue();
            ranges_with_error.w[i] = ranges.w[i];
            
        }

        
        int n_trial = 3000; 
        double alpha = 0.001; 
        double time_threshold = 50000;
        
     
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
		
	
	
	
	
	
	
	public void estimateDynamicSourceTDOA(int freq, int n_estimates, double thresh) throws Exception {
		
	    //Define default SGD parameters//
		int n_trial = 1000; 
        double alpha = 0.0001; 
        double time_threshold = 50000;
		
        double time = 0; 
        double refTime = -1.0;
        double tdoa = 0.0;
        double frameLength = .065536;
        double[] x;
        
        //Compute original source location//
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
		
        for (NavigationChannelList navList : filteredNavigationList) {
        	
            time = navList.getTimeStampAtFreq(freq);
            if(refTime < 0 && time > 0) { 
            	
                refTime = time;  //set reference time for first measurement            
                tdoa = 0.0;
                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa);                
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
                
                System.out.println(time + " " + refTime + " " + frameLength + " " + tdoa);
                
                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa);                
                tdoaFrequency0.add(td);
                
//                System.out.println(tdoaFrequency0.size() + " " + tdoa + " and diff " + (tdoa - tdoaFrequency0.get(tdoaFrequency0.size()-1).getTDOA()));
//                
//                if(Math.abs(tdoa - tdoaFrequency0.get(tdoaFrequency0.size()-1).getTDOA()) < thresh) {
//                	
//	                TimeDiffOnArrival td = new TimeDiffOnArrival(navList.getLongitude(),navList.getLatitude(), tdoa);                
//	                tdoaFrequency0.add(td);   
//                }	                
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
    	        
    	        //Now clear the list
    	        tdoaFrequency0.clear();
    	        count++;
    	        refTime = -1.0;
    	        tdoa = 0.0;
            }   
        }
        	
	}
	
	
	
	
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

    
	public void testSGDfdoa() throws Exception {
		
		Anchors myAnchors = new Anchors();
    	
	    for(int i = 0; i < filteredNavigationList.size(); i++) {
	    
	      NavigationChannelList nav = filteredNavigationList.get(i);	
	      
	      double[] s = nav.getLocalECEF();
	      double[] v = nav.getVelocity();
	      
	      myAnchors.setCoordinatesAndVelocity(s, v);	    
	      System.out.println(nav.getFDOA());
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
        for (int j = 0; j < num_anchors-1; j++) {                    	
        	ranges.w[j+1] = GradDescentResult.fdoaEstimate(anchors_in.getRow(j), anchors_in.getRow(j+1), 
        			velocities.getRow(j), velocities.getRow(j+1), sol);
            
        	xv[j+1] = j+1;
        }
	    
        
        Plot2DPanel plotfdoa = new Plot2DPanel();
        JFrame framefdoa = new JFrame("Plot of FDOAs");
        plotfdoa.addLinePlot("TDOA plot", xv, ranges.w);	        
        framefdoa.setSize(900, 700);
        framefdoa.setContentPane(plotfdoa);
        framefdoa.setVisible(true);
        

        
        int n_trial = 500; 
        double alpha = 0.0001; 
        double time_threshold = 50000;
        
     
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
       
       Plot2DPanel plot = new Plot2DPanel();       
       JFrame frame = new JFrame("FDOA Cost function value");
       plot.addLinePlot("TDOA plot", x, stockArr);
       frame.setSize(900, 700);
       frame.setContentPane(plot);
       frame.setVisible(true);

        
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
        
        //System.out.println(loc0[0] + " " + loc0[1] + " " + loc0[2]);
        
        double[] origin = NavigationChannelList.GeodeticToECEF(loc0[0], loc0[1], loc0[2]);        
        //System.out.println(origin[0] + " " + origin[1] + " " + origin[2] + "\n");
        
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
	
	
	
	
	
	public void createLocalCoordinateSystem() {
		
		double longitude = navigationList.get(0).getLongitude();
		double latitude = navigationList.get(0).getLatitude();
		double altitude = 0;
		
		localOrigin = NavigationChannelList.GeodeticToECEF(latitude, longitude, altitude);
		
		for (NavigationChannelList navList : navigationList) {			
			navList.GeodeticToLocal(localOrigin);			
		}
	}
	
	
	
	
	
	
	
	public static void main(String[] args) throws Exception {
		
		NavigationList navigation = new NavigationList();
		//navigation.createNavigationLog(new File("data/ChannelLog_ZWEI.log"));
		//navigation.createNavigationLog(new File("data/ChannelLog.log"));
		//navigation.createNavigationLog(new File("data/ChannelLog_SIVIRIEZ.log"));
		//navigation.createNavigationLog(new File("data/ChannelLog_ETZIKEN.log"));
		
		navigation.createNavigationLog_Interpolation(new File("data/ChannelLog_Lausanne.log"));
		
		final Kml kml = createCIRDocument(navigation.navigationList);
		kml.marshal(new File("CIRMarkers.kml"));
		
		
		int n_estimates = 5;
		int freqIndex = 0; 
		double threshold = .000001;
		
		navigation.computeVelocityECEF();
		navigation.filterOnRSSI(-150.0, freqIndex);
		navigation.filterUnique();
		navigation.createEstimationBounds();		
		navigation.setPlotTDOAs(true);
		navigation.estimateDynamicSourceTDOA(freqIndex, n_estimates, threshold);
		
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
            //.withName("CIR measure at " + navList.getTimeStamp())
            .withDescription(navList.getDescription())
            .withVisibility(true)
            .createAndSetPoint().addToCoordinates(navList.getLongitude(), navList.getLatitude());
        }
        kml.setFeature(document);
		
		return kml; 
	}
	
    public static Kml createSolutionDocument(ArrayList<double[]> estimates, ArrayList<Double> errors) {
		
		Kml kml = KmlFactory.createKml();
		Document document = kml.createAndSetDocument().withName("SolutionMarkers.kml");
		
		
		
		for(int j = 0; j < estimates.size(); j++) {
		  
		  double[] myEst = estimates.get(j);	
			
		  document.createAndAddPlacemark()
          //.withName("Error " + errors.get(j))
          .withDescription("Estimate No" + j + ": " + errors.get(j))
          .withVisibility(true)
          .createAndSetPoint().addToCoordinates(myEst[1], myEst[0]);
		  		  
		}
		
		System.out.println("Number of features: " + document.getFeature().size());
		
		List<Color> colors = pick(document.getFeature().size()+2);
		
		for (int i = 0; i < document.getFeature().size(); i++) {
			
			Color color = colors.get(i);
			Style style = document.getFeature().get(i).createAndAddStyle();			
			String hex = String.format("%02x%02x%02x", color.getRed(), color.getGreen(), color.getBlue());

			style.createAndSetIconStyle().withColor("e5"+hex).withScale(1.0);
			
			//style.setIconStyle((new IconStyle()).withColor(hex).withHeading(1.0).withScale(1.0));
			document.getFeature().get(i).getStyleSelector().add(style);
			
		}
		
        kml.setFeature(document);
		
		return kml; 		
	}
	
    
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
    
    
	//navigation.estimateAdaptiveSourceRSSI(txPower);
	//navigation.estimateSourceSolutionFDOA();
	//navigation.estimateAdaptiveSourceFDOA();
	//navigation.estimateSourceSolution();
	//navigation.estimateAdaptiveSource();
	//navigation.testSGDfdoa();
	//navigation.filterRandomlyUnique(3, 4);
	
	//navigation.filterFDOAfromNavigationList(0, 3.0);
	//navigation.filterTDOAfromNavigationList(0,.00001);
}
