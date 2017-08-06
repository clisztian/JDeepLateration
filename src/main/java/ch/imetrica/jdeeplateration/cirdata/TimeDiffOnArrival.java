package ch.imetrica.jdeeplateration.cirdata;

public class TimeDiffOnArrival {

	
	double longitude;
	double latitude;
	
	double[] sphericalCoords;
	
	double tdoa;
	
	public TimeDiffOnArrival(double longitude, double latitude, double tdoa) {
		
		this.longitude = longitude;
		this.latitude = latitude;
		this.tdoa = tdoa;
	}
	
	
	public void printTDOA() {
		System.out.println("Longitude: " + longitude + ", latitude: " + latitude + ", TDOA:" + tdoa);
	}
}
