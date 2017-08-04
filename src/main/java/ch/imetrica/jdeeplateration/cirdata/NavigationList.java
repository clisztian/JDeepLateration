package ch.imetrica.jdeeplateration.cirdata;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

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
				
		final Kml kml = createCIRDocument(navigation.navigationList);
		kml.marshal(new File("CIRMarkers.kml"));
		
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
	
}
