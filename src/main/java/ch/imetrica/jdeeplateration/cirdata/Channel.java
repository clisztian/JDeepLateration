package ch.imetrica.jdeeplateration.cirdata;

import java.util.ArrayList;

public class Channel {

	
	   int FreqIndex;
	   double RSSI; 
	   double SNR;
	   double ABStime; 
	   
	   int NbPeaks;
	   long TimeStamp;
	 
	   ArrayList<PeakList> peaks;
	   
	   
	   double[] power;
	   double[] timeDelay;
	   double[] freq_error;
	   double maxPower = -Double.MAX_VALUE;
 
	   
	   final static String delims = "|";
	   String[] tokens; 
	   
	   public Channel(long time, int freq, double Abs, double Rssi, double snr, 
			   int nbpeaks, String peaklist) {
		   
		   this.TimeStamp = time; 
		   this.FreqIndex = freq;
		   this.ABStime = Abs; 
		   this.SNR = snr; 
		   this.RSSI = Rssi; 
		   this.NbPeaks = nbpeaks; 
		   peaks = new ArrayList<PeakList>();	   
		  
		   power = new double[nbpeaks];
		   timeDelay = new double[nbpeaks];
		   freq_error = new double[nbpeaks];
		   
		   tokens = peaklist.split("[|]+");
		   
		   for(int i = 0; i < nbpeaks; i++) {	   
			   power[i] = (new Double(tokens[3*i + 0])).doubleValue();
			   timeDelay[i] = (new Double(tokens[3*i + 1])).doubleValue();
			   freq_error[i] = (new Double(tokens[3*i + 2])).doubleValue();
		   }

		   sortPeaks();		   
		   for(int i = 0; i < nbpeaks; i++) {
			   
			   PeakList peak = new PeakList(freq, power[i], timeDelay[i], freq_error[i]);		   
			   peaks.add(peak);
			   
			   if(power[i] > maxPower) {
				   maxPower = power[i];
			   }
		   }	   
	   }

	   public void printPeaks(String message) {
		   for(int i = 0; i < this.NbPeaks; i++) {	   
			   System.out.println("\n" + message + ": " + power[i] + ", " + timeDelay[i] + ", " + freq_error[i]);
		   }
	   }
	   

	   public void sortPeaks() {
		   //bubble sort: orders the peak list with respect to time; 
		   //first arriving peak on the top;
		      boolean swapped = true;
		      int j = 0;
		      double tmp_power;
		      double tmp_timeDelay;
		      double tmp_freq_error;
		      while (swapped) {
		            swapped = false;
		            j++;
		            for (int i = 0; i < NbPeaks - j; i++) {                                       
		                  if (timeDelay[i] > timeDelay[i + 1]) {         
		                        tmp_timeDelay = timeDelay[i];
		                        tmp_power = power[i];
		                        tmp_freq_error = freq_error[i];
		                        
		                        timeDelay[i] = timeDelay[i + 1];
		                        power[i] = power[i+1];
		                        freq_error[i] = freq_error[i+1];
		                        
		                        timeDelay[i + 1] = tmp_timeDelay;
		                        power[i+1] = tmp_power;
		                        freq_error[i+1] = tmp_freq_error;

		                        swapped = true;
		                  }
		            }                
		      }
	   }
	   
	   public double getFirstPeakFrequencyError() throws Exception {
		   if(peaks.size() == 0) {
			   throw new Exception("No peaks recorded");
		   }
		   return peaks.get(0).getFrequencyError();
	   }
	   
	   public double getFirstPeakTime() {
		   return peaks.get(0).getTimeDelay() * 0.000001; //offset time in .log-file is in Microseconds
	   }
	   
	   public int getFreqIndex() {
		   return FreqIndex; 
	   }
	   
	   public String toString() {
		   	  
		   String output = "FreqIndex=" + FreqIndex + ", ABStime=" + ABStime + ", NbPeaks=" + NbPeaks + "\n";	   
		   return output; 
	   }


	   public double getRSSI() {
		  return RSSI;	
	   }
	   
	   
	   public boolean isValidforFDOA() {
		   
		   int i = 1;  
		   double topPower = peaks.get(0).getPower();
		   if(peaks.size() > 1) {
			   
			 while(i < peaks.size()) {  
			  
			  if(peaks.get(i).getPower() > topPower) {
				  System.out.println("Power at " + i + " " + peaks.get(i).getPower() + ", topPower: " +  topPower);
				  return false;
			  }		  
			  i++;
			 }
		   }
		   return true;
	   }
	   
	   
	   public void dropFalseCorrelationPeaks(double threshold) {
		   
		   int i = 0;
		   if(peaks.size() > 1) {			  
			   while(i < peaks.size()) { 			
			    if(maxPower - peaks.get(i).getPower() > threshold) {
			    	peaks.remove(i);
			    }
			    else {i++;}
		       }		   
	       }
	   }
	   
	   
	   //ABStime - when peak was recieved
	   //RSSI (recieved signal strength indicator) 
	   //SignalNoise Ratio 
	   //PeakList (power of peak, time (micro) from Abstime, frequency error in Hz used 
	   //to doppler shift (towards the target if positive)
	   //34*10^6 - .4
	   
	   //Dopp shift Delta (f) = Delta (v)/c *fc (fc = carrier frequency = 2.GHz 
	   
	   
	   //Compute Delta(t)
	   //At frequency f0, compute ABStime_t/.01 and 
}
