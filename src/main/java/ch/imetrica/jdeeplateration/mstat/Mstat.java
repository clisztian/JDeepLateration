package ch.imetrica.jdeeplateration.mstat;

import java.io.Serializable;

import ch.imetrica.jdeeplateration.matrix.Matrix;

public class Mstat implements Serializable {

	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public static double distance(Matrix m1, Matrix m2) throws Exception
    {
		if (m1.size != m2.size) {
			throw new Exception("Expected same size matrices");
		}		
		double distance = 0;
	
		for(int i = 0; i < 3; i++) {
			distance += Math.abs(m1.w[i] - m2.w[i])*Math.abs(m1.w[i] - m2.w[i]); 
		}
		distance = Math.sqrt(distance);		
        return distance;
    }
	
	public static double distance(double[] m1, double[] m2) throws Exception
    {
		if (m1.length != m2.length) {
			throw new Exception("Expected same size matrices");
		}		
		double distance = 0;
	
		for(int i = 0; i < 3; i++) {
			distance += Math.abs(m1[i] - m2[i])*Math.abs(m1[i] - m2[i]); 
		}
		distance = Math.sqrt(distance);		
        return distance;
    }
	
	public static double norm(Matrix m1, Matrix m2) throws Exception
    {
		if (m1.size != m2.size) {
			throw new Exception("Expected same size matrices");
		}		
		double distance = 0;
	
		for(int i = 0; i < m1.size; i++) {
			distance += Math.abs(m1.w[i] - m2.w[i])*Math.abs(m1.w[i] - m2.w[i]); 
		}
		distance = Math.sqrt(distance)/m1.size;		
        return distance;
    }
	
	public static double fdoanorm(Matrix m1, Matrix m2) throws Exception
    {
		if (m1.size != m2.size) {
			throw new Exception("Expected same size matrices");
		}		
		double distance = 0;
	
		for(int i = 1; i < m1.size; i++) {
			distance += Math.abs((m1.w[i] - m1.w[i-1]) - m2.w[i])*Math.abs((m1.w[i] - m1.w[i-1]) - m2.w[i]); 
		}
		distance = Math.sqrt(distance)/m1.size;		
        return distance;
    }
	
	
	
}
