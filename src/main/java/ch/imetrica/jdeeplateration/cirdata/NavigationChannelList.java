package ch.imetrica.jdeeplateration.cirdata;

import java.util.ArrayList;

import org.apache.commons.math3.stat.descriptive.moment.VectorialCovariance;

public class NavigationChannelList {

	
	final static double semiMajorAxis = 6378137.0; 
	final static double firstEccentricitySquared=.00669437999014;
	
	ArrayList<Channel> ChannelList; 
	
	long TimeStamp;
	int DaySeconds; 
	
	double fdoa; 
	double longitude;
	double latitude;	
	double altitude = 0;
	
	double[] velocity = null;
	double[] localECEF = null;
	double[] localOrigin = null;
	
	public NavigationChannelList(long time, int seconds, double longitude, double lat) {
		
		this.TimeStamp = time; 
		this.DaySeconds = seconds; 
		this.latitude = lat; 
		this.longitude = longitude;		
		
		ChannelList = new ArrayList<Channel>();		
	}
	
	public void setFrequencyError(int freq) throws Exception {

		for(Channel ch : ChannelList) {
		
			if(ch.getFreqIndex() == freq) {
				fdoa =  ch.getFirstPeakFrequencyError();
				break;
			}
		}
	}
	
	public NavigationChannelList(long time, int seconds, double longitude, double lat, double alt) {
		
		this.TimeStamp = time; 
		this.DaySeconds = seconds; 
		this.latitude = lat; 
		this.longitude = longitude;		
		this.altitude = alt; 
		
		ChannelList = new ArrayList<Channel>();		
	}
	
	public void GeodeticToLocal(double[] localOrigin) {
		
		if(localECEF == null) {
		   	
			localECEF = new double[3];
			
			double f = 1.0/298.257224;
			double f1 = (1.0 - f)*(1.0 - f);
			
			double loc_longitude = (Math.PI/180.0)*longitude;
		    double loc_latitude =  (Math.PI/180.0)*latitude;
		           
		    double Cn = 1.0/Math.sqrt(Math.cos(loc_latitude)*Math.cos(loc_latitude) +
		    		          f1*Math.sin(loc_latitude)*Math.sin(loc_latitude));      
		    		
		    double S = Cn*f1;
		    		    
		    localECEF[0] = (semiMajorAxis*Cn+altitude)*Math.cos(loc_latitude)*Math.cos(loc_longitude);
		    localECEF[1] = (semiMajorAxis*Cn+altitude)*Math.cos(loc_latitude)*Math.sin(loc_longitude);
		    localECEF[2] = (semiMajorAxis*S+altitude)*Math.sin(loc_latitude);
		    
		    localECEF[0] = localECEF[0] - localOrigin[0];
		    localECEF[1] = localECEF[1] - localOrigin[1];
		    localECEF[2] = localECEF[2] - localOrigin[2];
		}
	}
	

	
	public void computeVelocity(long prevTime, double prevLongitude, double prevLatitude) {
		
		double meterPerDegreeLat = 111166.01;
		double meterPerDegreeLong = 76404.09;
		velocity = new double[3];

		double timeDiff = (double)(TimeStamp - prevTime)/1000.0;
		double longDiffmeters = (longitude - prevLongitude)*meterPerDegreeLong;
		double latDiffmeters = (latitude - prevLatitude)*meterPerDegreeLat;
		
		velocity[0] = longDiffmeters/timeDiff;
		velocity[1] = latDiffmeters/timeDiff;		
	}
	
	public void computeVelocityECEF(long prevTime, double[] prevPos) {
		
		velocity = new double[3];

		double timeDiff = (double)(DaySeconds - prevTime);
		
		for(int i = 0; i < 3; i++) {
			velocity[i] = (localECEF[i] - prevPos[i])/timeDiff;			
		}	
	}
	
	
	public long getTimeStamp() {
		return TimeStamp;
	}
	
	public long getDaySeconds() {
		return DaySeconds;
	}
	
	public double[] getCoordinates() {
		double[] coords = new double[2]; 
		
		coords[0] = latitude; 
		coords[1] = longitude;
		
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
	
	
	public void addChannel_interpolation(String channel, int f) {
		
		String[] tokens = channel.split("\\s+");
		
		long time = (new Long(tokens[0])).longValue();
		int freq = f;
		
		String[] paramtoks = tokens[4].split("[=]+");
		double Abs = (new Double(paramtoks[1])).doubleValue();
		
		paramtoks = tokens[5].split("[=]+");
		double Rssi = (new Double(paramtoks[1])).doubleValue();
		
		paramtoks = tokens[6].split("[=]+");
		double snr = (new Double(paramtoks[1])).doubleValue();
		
		paramtoks = tokens[7].split("[=]+");
		int nbpeaks = (new Integer(paramtoks[1])).intValue();
		
		paramtoks = tokens[8].split("[=]+");
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

	public void printVelocity() {
		if(velocity != null) {
			System.out.println("<" + velocity[0] + ", " + velocity[1] + ", 0>");
		}		
	}
	
	static private double cbrt(double x) {
		
        if (x >= 0)
            return Math.pow(x, 1.0/3.0);
        else
            return -Math.pow(Math.abs(x), 1.0/3.0);		
	}

	
	static public double[] GeodeticToECEF(double latitude, double longitude, double altitude) {
		 	
			double[] local = new double[3];
	
			double f = 1.0/298.257224;
			double f1 = (1.0 - f)*(1.0 - f);
			
			double loc_longitude = (Math.PI/180.0)*longitude;
		    double loc_latitude =  (Math.PI/180.0)*latitude;
		           
		    //double N = semiMajorAxis/Math.sqrt(1.0 - firstEccentricitySquared*Math.pow(Math.sin(loc_latitude),2));
		    
		    double Cn = 1.0/Math.sqrt(Math.cos(loc_latitude)*Math.cos(loc_latitude) +
		    		          f1*Math.sin(loc_latitude)*Math.sin(loc_latitude));      
		    		
		    double S = Cn*f1;
		    		    
		    local[0] = (semiMajorAxis*Cn+altitude)*Math.cos(loc_latitude)*Math.cos(loc_longitude);
		    local[1] = (semiMajorAxis*Cn+altitude)*Math.cos(loc_latitude)*Math.sin(loc_longitude);
		    local[2] = (semiMajorAxis*S+altitude)*Math.sin(loc_latitude);
		    
		    return local; 
	}
	
	
	static public double[] ECEFToGeodetic(double x, double y, double z) {
		
	
		 double[] coords = new double[2];
	     double a = 6378137.0;
	     double b = 6356752.3142;
	     double esq = 6.69437999014 * 0.001;
	     double e1sq = 6.73949674228 * 0.001;
	     double f = 1.0 / 298.257223563;   

	     double r = Math.sqrt(x * x + y * y);
	     double Esq = a * a - b * b;
	     double F = 54.0 * b * b * z * z;
	     double G = r * r + (1.0 - esq) * z * z - esq * Esq;
	     double C = (esq * esq * F * r * r) / (Math.pow(G, 3.0));
	     double S = cbrt(1.0 + C + Math.sqrt(C * C + 2 * C));
	     double P = F / (3 * Math.pow((S + 1 / S + 1), 2) * G * G);
	     double Q = Math.sqrt(1 + 2 * esq * esq * P);
	     double r_0 =  -(P * esq * r) / (1 + Q) + Math.sqrt(0.5 * a * a*(1 + 1.0 / Q) -
	               P * (1 - esq) * z * z / (Q * (1 + Q)) - 0.5 * P * r * r);
	     double U = Math.sqrt(Math.pow((r - esq * r_0), 2) + z * z);
	     double V = Math.sqrt(Math.pow((r - esq * r_0), 2) + (1 - esq) * z * z);
	     double Z_0 = b * b * z / (a * V);
	     double h = U * (1.0 - b * b / (a * V));
	     double lat = Math.atan((z + e1sq * Z_0) / r);
	     double lon = Math.atan2(y, x);
	
	     coords[0] = 180.0*lat/Math.PI;
	     coords[1] = 180.0*lon/Math.PI;
	
	     return coords; 
    }

	public double getRSSI(int freq) {
		
		 double rssi = -Double.MAX_VALUE;
         for(Channel chan : ChannelList) {
			
			if(chan.getFreqIndex() == freq) {
			   rssi = chan.getRSSI();
			   break;
			}
         }	
		
		return rssi;
	}

	public static double[] GeodeticToECEF(double latitude, double longitude, double altitude, double[] localOrig) {
		
		double[] local = new double[3];
		
		double f = 1.0/298.257224;
		double f1 = (1.0 - f)*(1.0 - f);
		
		double loc_longitude = (Math.PI/180.0)*longitude;
	    double loc_latitude =  (Math.PI/180.0)*latitude;
	           
	    double Cn = 1.0/Math.sqrt(Math.cos(loc_latitude)*Math.cos(loc_latitude) +
	    		          f1*Math.sin(loc_latitude)*Math.sin(loc_latitude));      
	    		
	    double S = Cn*f1;
	    		    
	    local[0] = (semiMajorAxis*Cn+altitude)*Math.cos(loc_latitude)*Math.cos(loc_longitude);
	    local[1] = (semiMajorAxis*Cn+altitude)*Math.cos(loc_latitude)*Math.sin(loc_longitude);
	    local[2] = (semiMajorAxis*S+altitude)*Math.sin(loc_latitude);
	    
	    local[0] -= localOrig[0];
	    local[1] -= localOrig[1];
	    local[2] -= localOrig[2];
	    
	    return local; 
	}

	public double[] getLocalECEF() {
		return localECEF;
	}

	public double[] getVelocity() throws Exception {
		
		if(velocity == null) {
			throw new Exception("velocity not yet defined");
		}
		return velocity;
	}
	
	public double getFDOA() {
		return fdoa;
	}
	
	public void setVelocityToZero() {
		velocity = new double[3];
	}

	public boolean isStationary(double thresh) {
		
		boolean it_is = true;		
		it_is = (closeToZero(velocity[0], thresh) && closeToZero(velocity[1], thresh));		
		return it_is;
	}
	
	
	public boolean closeToZero(double val, double thresh) {
		
		if(Math.abs(val) < thresh) {return true;}
		else {return false;}
	}
	
	
}
