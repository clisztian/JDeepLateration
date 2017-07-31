package ch.imetrica.jdeeplateration.cirdata;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;


public class NavigationList {

	double[] Frequencies; 
	double CenterFreq; 
	double Bandwidth;
	
	long creationTimestamp;
	
	ArrayList<NavigationChannelList> navigationList;
	
	public NavigationList() {
		navigationList = new ArrayList<NavigationChannelList>();
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
        	throw new Exception("File not in correct format " + tokens[3]);
        }
        this.creationTimestamp = (new Long(tokens[0])).longValue();
        		
        //----------------- get frequencies
        strline = br.readLine(); 
        System.out.println(strline);
        tokens = strline.split("\\s+");
        if(!tokens[3].equals("FREQUENCY")) {
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
            		  navigationList.add(channelList);
            		  count++;
            	  }
            	  else {           		  
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
		navigation.createNavigationLog(new File("data/ChannelLog.log"));
		System.out.println("Printing history");
		navigation.printNavigationHistory();
		
		
		
		static Layer addPoint(double latitude, double longitude) {
		    SimpleFeatureTypeBuilder b = new SimpleFeatureTypeBuilder();

		    b.setName("MyFeatureType");
		    b.setCRS(DefaultGeographicCRS.WGS84);
		    b.add("location", Point.class);
		    // building the type
		    final SimpleFeatureType TYPE = b.buildFeatureType();

		    SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(TYPE);
		    GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory();
		    Point point = geometryFactory.createPoint(new Coordinate( latitude, longitude));
		    featureBuilder.add(point);
		    SimpleFeature feature = featureBuilder.buildFeature(null);
		    DefaultFeatureCollection featureCollection = new DefaultFeatureCollection("internal", TYPE);
		    featureCollection.add(feature);
		    Style style = SLD.createSimpleStyle(TYPE,Color.red);

		    Layer layer = new FeatureLayer(featureCollection, style);
		    return layer;
		  }
		
		
		
		
		
		
		
		
		
	}
	
}
