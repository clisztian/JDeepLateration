package ch.imetrica.jdeeplateration.anchors;

import java.util.ArrayList;

import ch.imetrica.jdeeplateration.matrix.Matrix;


public class Anchors {

	/**
	 * 
	 */
	//-- Locations of anchors 
	ArrayList<double[]> coordinates = null;
	ArrayList<double[]> velocity = null;
	
    Matrix anchors = null;
	Matrix anchor_velocity = null;
    
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
    		throw new Exception("Must define velocities");
    	}
    	return anchors; 
    }
    
    public Matrix getVelocities() throws Exception {
    	
    	if(anchor_velocity == null) {
    		throw new Exception("Must define velocities");
    	}
    	return anchor_velocity; 
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
	
	public void setCoordinatesAndVelocity(double[] a, double[] v) throws Exception {
		if(a.length != 3 ) {
			throw new Exception("Must have 3 coordinates");
		}
		if(coordinates == null) {
			coordinates = new ArrayList<double[]>();
		}
		if(velocity == null) {
			velocity = new ArrayList<double[]>();
		}
		
		coordinates.add(a);
		velocity.add(v);
		
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


	public void commitCoordinatesAndVelocity() throws Exception {
		
		if(coordinates.size() < 4) {
			throw new Exception("Must have at least 4 locations");
		}
		
		double[] vec_coordinates = new double[coordinates.size()*3];
		double[] vec_velocity = new double[velocity.size()*3];
		
		for(int i = 0; i < coordinates.size(); i++) {
		
			double[] temp = coordinates.get(i);
			for(int j = 0; j < 3; j++) {
				vec_coordinates[i*3 + j] = temp[j];
			}	
			
			double[] temp_v = velocity.get(i);
			for(int j = 0; j < 3; j++) {
				vec_velocity[i*3 + j] = temp_v[j];
			}						
		}
		this.number_anchors = coordinates.size();
		anchors = new Matrix(this.number_anchors, 3);
		anchors.setPointer(vec_coordinates);
		
		anchor_velocity = new Matrix(this.number_anchors, 3);
		anchor_velocity.setPointer(vec_velocity);
	
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

	public Anchors subset(int nrows) throws Exception {
		
        Anchors subset = new Anchors();
        
        for(int j = 0; j < nrows; j++) {
        	
        	Matrix row = this.getRow(j);
        	subset.setCoordinates(row.w);        	
        
        }
        subset.commitCoordinates();        
		return subset; 
	}
	
	
}
