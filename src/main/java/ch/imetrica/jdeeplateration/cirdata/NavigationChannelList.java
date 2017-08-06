package ch.imetrica.jdeeplateration.cirdata;

import java.util.ArrayList;

public class NavigationChannelList {

	
	ArrayList<Channel> ChannelList; 
	
	long TimeStamp;
	int DaySeconds; 
	
	double longitude;
	double latitude;
	
	
	public NavigationChannelList(long time, int seconds, double longitude, double lat) {
		
		this.TimeStamp = time; 
		this.DaySeconds = seconds; 
		this.latitude = lat; 
		this.longitude = longitude;
		
		ChannelList = new ArrayList<Channel>();		
	}
	
	public long getTimeStamp() {
		return TimeStamp;
	}
	
	public double[] getCoordinates() {
		double[] coords = new double[2]; 
		
		coords[0] = longitude; 
		coords[1] = latitude;
		
		return coords;   
	}
	
	public double getTimeStampAtFreq(int i) {
	   
		double time = -1.0;
		for(Channel chan : ChannelList) {
			
			if(chan.getFreqIndex() == i) {
				time = chan.ABStime;
				break;
			}			
		}
		return time;		
	}
	
	
	public String getPointCoordinates() {
		
		return new String("" + longitude + "," + latitude + ",0");
	}
	
	public double getLongitude() {
		return longitude; 
	}
	
	public double getLatitude() {
		return latitude;
	}
	
	public void printCoordinates() {
		
		System.out.println("Longitude: " + longitude + ", latitude: " + latitude);
	}
	
	public void addChannel(String channel) {
		
		String[] tokens = channel.split("\\s+");
		
		long time = (new Long(tokens[0])).longValue();
		
		String[] paramtoks = tokens[4].split("[=]+");
		int freq = (new Integer(paramtoks[1])).intValue();
		
		paramtoks = tokens[5].split("[=]+");
		double Abs = (new Double(paramtoks[1])).doubleValue();
		
		paramtoks = tokens[6].split("[=]+");
		double Rssi = (new Double(paramtoks[1])).doubleValue();
		
		paramtoks = tokens[7].split("[=]+");
		double snr = (new Double(paramtoks[1])).doubleValue();
		
		paramtoks = tokens[8].split("[=]+");
		int nbpeaks = (new Integer(paramtoks[1])).intValue();
		
		paramtoks = tokens[9].split("[=]+");
		String peaklist = paramtoks[1];

		Channel myChannel = new Channel(time, freq, Abs, Rssi, snr, nbpeaks, peaklist);
		
		ChannelList.add(myChannel);
	}

	public String getDescription() {
		
		 String description = "";
		 for(Channel ch : ChannelList) {
			 
			 description += ch.toString();
		 }
		 return description;
	}
	
	
	
}
