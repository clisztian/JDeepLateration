package ch.imetrica.jdeeplateration.matrix;
import java.io.Serializable;
import java.util.Random;


public class Matrix implements Serializable {
	
	private static final long serialVersionUID = 1L;
	
	public int size;
	public int rows;
	public int cols;
	
	public double[] w;
	public double[] dw;

	@Override
	public String toString() {
		String result = "";
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				result += String.format("%.4f",getW(r, c)) + "\t";
			}
			result += "\n";
		}
		return result;
	}
	
	public Matrix clone() {
		Matrix result = new Matrix(rows, cols);
		for (int i = 0; i < w.length; i++) {
			result.w[i] = w[i];
			result.dw[i] = dw[i];
		}
		return result;
	}

	public Matrix concat(double[][] data) throws Exception {
		
		int rows = data.length; 
		int cols = data[0].length;
		Matrix out = new Matrix(rows, cols);
						
	    for(int i = 0; i < rows; i++) {
	        for(int j = 0; j < cols; j++){
	            out.w[i*cols + j] = data[i][j];
	        }
	    }
	     	
		return out; 
	}
	
	public Matrix stack(Matrix data) throws Exception {
		
		if (this.cols != data.cols) {
			throw new Exception("Expected same number of columns");
		}
		
		int rows = data.rows; 
		int cols = data.cols;
		
		Matrix out = new Matrix(this.rows + rows, this.cols);
			
		for(int i = 0; i < this.rows; i++) {
	        for(int j = 0; j < this.cols; j++){
	            out.w[i*this.cols + j] = this.w[i*this.cols + j];
	        }
	    }
		
	    for(int i = 0; i < rows; i++) {
	        for(int j = 0; j < cols; j++){
	            out.w[this.size + i*cols + j] = data.w[i*cols + j];
	        }
	    }
	     	
		return out; 
	}
	
	//--- returns the minimum value in the ith column
	public double ColumnMin(int i) throws Exception {
		
		if (this.cols <= i) {
			throw new Exception("Exceeds number of columns");
		}
		
		double min = Double.MAX_VALUE;
		
		for(int j = 0; j < this.rows; j++) {
			if(this.w[j*this.cols + i] < min) {
				min = this.w[j*this.cols + i];
			}
		}
		
		return min;		
	}
	
	//--- returns the minimum value in the ith column
	public double ColumnMax(int i) throws Exception {
		
		if (this.cols <= i) {
			throw new Exception("Exceeds number of columns");
		}
		
		double max = -Double.MAX_VALUE;		
		for(int j = 0; j < this.rows; j++) {
			if(this.w[j*this.cols + i] > max) {
				max = this.w[j*this.cols + i];
			}
		}		
		return max;		
	}
	
	
	public void normalizeBatch() {
		
		for(int i = 0; i < this.rows; i++) {
			
			double sumcol = 0; 
			double var = 0;
			
			for(int j = 0; j < this.cols; j++) {
				sumcol += this.getW(i, j)/this.cols;
			}
			
			for(int j = 0; j < this.cols; j++) {
				var += (this.getW(i, j) - sumcol)*(this.getW(i, j) - sumcol);
			}
			
			for(int j = 0; j < this.cols; j++) {
				this.w[i*this.cols + j] = (this.getW(i, j) - sumcol)/Math.sqrt(var);
			}			
		}		
	}
	
	public void resetDw() {
		for (int i = 0; i < dw.length; i++) {
			dw[i] = 0;
		}
	}

	
	public static Matrix transpose(Matrix m) {
		Matrix result = new Matrix(m.cols, m.rows);
		for (int r = 0; r < m.rows; r++) {
			for (int c = 0; c < m.cols; c++) {
				result.setW(c, r, m.getW(r, c));
			}
		}
		return result;
	}
	
	public static Matrix rand(int rows, int cols, double initParamsStdDev, Random rng) {
		Matrix result = new Matrix(rows, cols);
		for (int i = 0; i < result.w.length; i++) {
			result.w[i] = rng.nextDouble() * 2*initParamsStdDev - initParamsStdDev;
			//result.w[i] = .2 * initParamsStdDev;
		}
		return result;
	}
	
	public static Matrix ident(int dim) {
		Matrix result = new Matrix(dim, dim);
		for (int i = 0; i < dim; i++) {
			result.setW(i, i, 1.0);
		}
		return result;
	}
	
	public static Matrix uniform(int rows, int cols, double s) {
		Matrix result = new Matrix(rows, cols);
		for (int i = 0; i < result.w.length; i++) {
			result.w[i] = s;
		}
		return result;
	}
	
	public void printMatrix()
	{
		for(int i = 0; i < this.rows; i++) {
			for(int j = 0; j < this.cols; j++) {
				System.out.print(w[this.cols*i + j] + ", ");
			}
			System.out.println("");
		}		
	}
	
	public void printMatrixDW()
	{
		for(int i = 0; i < this.w.length; i++)
		{System.out.print(dw[i] + ", ");}
		System.out.println("");
	}
	
	
	public static Matrix ones(int rows, int cols) {
		return uniform(rows, cols, 1.0);
	}
	
	public static Matrix negones(int rows, int cols) {
		return uniform(rows, cols, -1.0);
	}
	
	public Matrix(int dim) {
		this.rows = dim;
		this.cols = 1;
		this.size = rows*cols;
		this.w = new double[rows * cols];
		this.dw = new double[rows * cols];
	}
	
	public Matrix(int rows, int cols) {
		this.rows = rows;
		this.cols = cols;
		this.size = rows*cols;
		this.w = new double[rows * cols];
		this.dw = new double[rows * cols];
	}
	
	public Matrix(double[] vector) {
		this.rows = vector.length;
		this.cols = 1;
		this.size = rows*cols;
		this.w = vector;
		this.dw = new double[vector.length];
	}
	
	public Matrix(double[] vector, int n, int batchsize) throws Exception 
	{		
		if (n*batchsize != vector.length) {
			throw new Exception("matrix dimension mismatch");
		}
		
		this.rows = n;
		this.cols = batchsize;
		this.size = rows*cols;
		this.w = vector;
		this.dw = new double[vector.length];
	}
	
	private int index(int row, int col) {
		int ix = cols * row + col;
		return ix;
	}
	
	public void set(int row, int col, double val) {
		w[index(row, col)] = val;
	}
	
	public double getW(int row, int col) {
		return w[index(row, col)];
	}
	
	public double[] returnColumnW(int column) throws Exception {
		
		if (column >= this.cols) {
			throw new Exception("matrix dimension mismatch");
		}
		
		double[] myCol = new double[rows];
		for(int i = 0; i < rows; i++) myCol[i] = this.getW(i, column);		
		return myCol;
		
	}
	
	public void add(Matrix rhs) throws Exception {
		
		if (this.size != rhs.size) {
			throw new Exception("matrix dimension mismatch");
		}
		
		for(int i = 0; i < this.size; i++) {
			this.w[i] += rhs.w[i];
		}
		
	}
	
	public Matrix getRow(int row) throws Exception {
		
		if (row >= this.rows) {
			throw new Exception("matrix dimension mismatch");
		}
		
		double[] myRow = new double[this.cols];
		for(int i = 0; i < this.cols; i++) {
			  myRow[i] = this.getW(row, i);		
		}
		return new Matrix(myRow);	
	}
	
	
	private void setW(int row, int col, double val) {
		w[index(row, col)] = val;
	}

	public void scale(double d) {
		for(int i = 0; i < w.length; i++) w[i] *= d;		
	}

	public void setRow(int i, Matrix estimator) throws Exception {
		
		if (estimator.size != this.cols) {
			throw new Exception("matrix dimension mismatch");
		}
		
		for(int j = 0; j < this.cols; j++) {
			setW(i, j, estimator.w[j]);
		}
	}
	
	public void setPointer(double[] wi) throws Exception {
		
		if(wi.length != size) {
			throw new Exception("matrix dimension mismatch");
		}
		
		w = wi;
	}
	
	
	
}
