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
	  
	   tokens = peaklist.split("[|]+");
	   
	   for(int i = 0; i < nbpeaks; i++) {
		   
		   double power = (new Double(tokens[3*i + 0])).doubleValue();
		   double timeDelay = (new Double(tokens[3*i + 1])).doubleValue();
		   double freq_error = (new Double(tokens[3*i + 2])).doubleValue();
		   
		   PeakList peak = new PeakList(freq, power, timeDelay, freq_error);		   
		   peaks.add(peak);
	   }
	   
   }
   
   public double getFirstPeakFrequencyError() throws Exception {
	   if(peaks.size() == 0) {
		   throw new Exception("No peaks recorded");
	   }
	   return peaks.get(0).getFrequencyError();
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
