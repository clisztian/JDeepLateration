package ch.imetrica.jdeeplateration.cirdata;

public class Peak {

	
	int frequencyNumber;  //0,1,2
	double frequency;
	
	double power;
	double timeDelay; //micro from Abstime, 
	double frequency_error;
	
	
	
	Peak(int freqNum, double freq, double power, double timeDelay, double freq_error) {
		
		this.frequencyNumber = freqNum;
		this.frequency = freq;
		this.power = power;
		this.timeDelay = timeDelay; 
		this.frequency_error = freq_error;
		
	}
}
