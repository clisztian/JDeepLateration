package ch.imetrica.jdeeplateration.anchors;

import java.util.ArrayList;

import ch.imetrica.jdeeplateration.matrix.Matrix;


public class Anchors {

	/**
	 * 
	 */
	//-- Locations of anchors 
	ArrayList<double[]> coordinates = null;
    Matrix anchors = null;
	
	int number_anchors = 0;
	
	private static final long serialVersionUID = 1L;


	public Anchors(int n) {
		this.number_anchors = n;		
		coordinates = new ArrayList<double[]>();
	}
	
	public Anchors() {
		number_anchors = 0;
		coordinates = new ArrayList<double[]>();
	}
	
    public Matrix getAnchors() throws Exception {
    	if(anchors == null) {
    		throw new Exception("Must have 3 coordinates");
    	}
    	return anchors; 
    }
    
    public int getNumberOfAnchors() {
    	return number_anchors;
    }
	
    public Matrix getRow(int i) throws Exception {
    	return anchors.getRow(i);
    }
    
	public void setCoordinates(double[] a) throws Exception {
		if(a.length != 3 ) {
			throw new Exception("Must have 3 coordinates");
		}
		if(coordinates == null) {
			coordinates = new ArrayList<double[]>();
		}
		coordinates.add(a);
		number_anchors++;
	}
	
	
	public void setCoordinates(double a, double b, double c) {

		if(coordinates == null) {
			coordinates = new ArrayList<double[]>();
			number_anchors = 0;
		}
		
		double[] loc = {a, b, c};
		coordinates.add(loc);
		number_anchors++;
	}
	
	public void commitCoordinates() throws Exception {
		if(coordinates.size() < 4) {
			throw new Exception("Must have at least 4 locations");
		}
		
		double[] vec_coordinates = new double[coordinates.size()*3];
		for(int i = 0; i < coordinates.size(); i++) {
			double[] temp = coordinates.get(i);
			for(int j = 0; j < 3; j++) {
				vec_coordinates[i*3 + j] = temp[j];
			}	
		}
		this.number_anchors = coordinates.size();
		anchors = new Matrix(this.number_anchors, 3);
		anchors.setPointer(vec_coordinates);
	}

	public void printMatrix() {
		anchors.printMatrix();		
	}
	
	public double getColumnMin(int i) throws Exception {
		return anchors.ColumnMin(i);
	}
	
	public double getColumnMax(int i) throws Exception {
		return anchors.ColumnMax(i);
	}
	
	
}
