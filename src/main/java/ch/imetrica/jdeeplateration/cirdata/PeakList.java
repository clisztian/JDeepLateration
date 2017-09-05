package ch.imetrica.jdeeplateration.cirdata;

public class PeakList {

	
	int frequencyNumber;  //0,1,2
	double frequency;
	
	double power;
	double timeDelay; //micro from Abstime, 
	double frequency_error;
	
	
	PeakList(int freqNum, double power, double timeDelay, double freq_error) {
		
		this.frequencyNumber = freqNum;
		this.power = power;
		this.timeDelay = timeDelay; 
		this.frequency_error = freq_error;
		
	}
	
	public double getFrequencyError() {
		return frequency_error;
	}
	
}
