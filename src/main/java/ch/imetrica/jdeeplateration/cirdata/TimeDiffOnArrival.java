package ch.imetrica.jdeeplateration.cirdata;

public class TimeDiffOnArrival {

	
	double longitude;
	double latitude;
	
	double[] sphericalCoords;
	
	double tdoa;
	double fdoa;
	
	public TimeDiffOnArrival(double longitude, double latitude, double tdoa) {
		
		this.longitude = longitude;
		this.latitude = latitude;
		this.tdoa = tdoa;
	}
	
	
	public void printTDOA() {
		System.out.println("Longitude: " + longitude + ", latitude: " + latitude + ", TDOA:" + tdoa);
	}


	public double getLongitude() {	
		return this.longitude;
	}


	public double getLatitude() {
		return this.latitude;
	}


	public Double getTDOA() {
		return this.tdoa;
	}
}
