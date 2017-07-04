package ch.imetrica.jdeeplateration.mstat;

import java.util.Random;

import ch.imetrica.jdeeplateration.matrix.Matrix;

public class GradDescentResult {

	
	public Matrix estimator;
	public Matrix estimator_candidate;
	public Matrix error; 
	
	
	
	public void setEstimator(Matrix estimator) {
		this.estimator = estimator;
	}
	
	public void setEstimatorCandidate(Matrix estimator) {
		this.estimator = estimator;
	}
	

	public GradDescentResult(int n, int dim)
    {
        estimator_candidate = new Matrix(n, dim);
        error = new Matrix(n, 1);
    }
	
	
	public static GradDescentResult GradDescent(Matrix anchors_in, Matrix ranges_in,
            Matrix bounds_in, int n_trial, double alpha, double time_threshold) throws Exception {
		
	
            Random random = new Random();

            int n = anchors_in.rows;
            int dim = anchors_in.cols;
            
            GradDescentResult gdescent_result = new GradDescentResult(n_trial, dim);

            if (bounds_in == null) {
                bounds_in = new Matrix(1, dim);
            }
            
            Matrix bounds_temp = anchors_in.stack(bounds_in);
            Matrix bounds = new Matrix(2, dim);
            
            for(int i = 0; i < dim; i++) {
                bounds.set(0, i, bounds_temp.ColumnMin(i));
                bounds.set(1, i, bounds_temp.ColumnMax(i));
            }
            
            
            Matrix ranges = new Matrix(n,1);
            for (int i = 0; i < n_trial; i++)
            {
            	Matrix estimator0 = new Matrix(dim);
            	
                for (int j = 0; j < dim; j++) {
                    estimator0.w[j] = random.nextDouble() * (bounds.getW(1, j) - bounds.getW(0, j)) + bounds.getW(0, j);
                }
                Matrix estimator = new Matrix(estimator0.w);

            
                long startTime = System.nanoTime();
                while (true)
                {
                    for (int j = 0; j < n; j++) {
                        ranges.w[j] = Mstat.distance(anchors_in.getRow(j), estimator);
                    }
                    
                    double error = Mstat.distance(ranges_in, ranges);
                    Matrix delta = new Matrix(dim);
                    
                    for (int j = 0; j < n; j++) {
                   	
                        delta.add(sub(estimator, anchors_in.getRow(j)));
                        delta.scale((ranges_in.w[j] - ranges.w[j]) / ranges.w[j]);
                    }
                    
                    delta.scale(2*alpha);
                    Matrix estimator_next = sub(estimator,delta);
                    
                    for (int j = 0; j < n; j++) {
                        ranges.w[j] = Mstat.distance(anchors_in.getRow(j), estimator_next);
                    }
                    
                    double error_next = Mstat.distance(ranges_in, ranges);
                    if (error_next < error) {
                        estimator = estimator_next;
                    }
                    else
                    {
                        gdescent_result.estimator_candidate.setRow(i, estimator);
                        gdescent_result.error.w[i] = error;
                        break;
                    }
                    
                   
                    if (System.nanoTime() - startTime > time_threshold)
                    {
                        gdescent_result.error.w[i] = Double.MAX_VALUE;
                        break;
                    }
                }
            }
            
            return gdescent_result;
            
	}

	
	
	public static GradDescentResult mlat(Matrix anchors_in, Matrix ranges_in, Matrix bounds_in, int n_trial, double alpha, double time_threshold) throws Exception {
		
            GradDescentResult gdescent_result = GradDescent(anchors_in, ranges_in, bounds_in, n_trial, alpha, time_threshold);

            int idx = -1;
            double error = Double.MAX_VALUE;
            
            for (int i = 0; i < gdescent_result.error.size; i++)
            {
                if (gdescent_result.error.w[i] < error)
                {
                    idx = i;
                    error = gdescent_result.error.w[i];
                }
            }
            gdescent_result.estimator = gdescent_result.estimator_candidate.getRow(idx);
            return gdescent_result;
    }	
	
	
	
	private static Matrix sub(Matrix estimator2, Matrix delta) throws Exception {
		
		if (estimator2.size != delta.size) {
			throw new Exception("matrix dimension mismatch");
		}
		
		Matrix diff = new Matrix(delta.size);
		for(int i = 0; i < delta.size; i++) {
			diff.w[i] = estimator2.w[i] - delta.w[i];
		}
		
		return diff;
	}

	
	
	
	public static void main(String[] args) throws Exception {
		
		
		Random random = new Random();
		
		double error = 0.5;
        double W = 9, L = 9, H = 3;
        
        Matrix anchors = new Matrix(4,3);
        double[] anch = {0,0,H,W,0,H,W,L,H,0,L,H};
        anchors.setPointer(anch);
        
        Matrix node = new Matrix(3);
        //double[] n = {W*random.nextDouble(), L*random.nextDouble(), H*random.nextDouble()};
        double[] n = {7.46429794,  1.26322155,  2.95852001};
        node.setPointer(n);       
 
        Matrix ranges = new Matrix(anchors.rows);
           
        Matrix ranges_with_error = new Matrix(anchors.rows);
        for (int i = 0; i < anchors.rows; i++)
        {
            ranges.w[i] = Mstat.distance(anchors.getRow(i), node);
            ranges_with_error.w[i] = ranges.w[i];  //+ 2 * error * (random.nextDouble() - 0.5);
        }

        
        int n_trial = 100; 
        double alpha = 0.001; 
        double time_threshold = 10000;
        Matrix bounds_in = null;
        
        GradDescentResult gdescent_result = GradDescentResult.mlat(anchors, ranges_with_error, bounds_in, n_trial, alpha, time_threshold);

        System.out.println("Anchors");
        anchors.printMatrix();
        
        System.out.println("Node");
        node.printMatrix();
        
        System.out.println("Ranges");
        ranges.printMatrix();
        
        System.out.println("Ranges with error");
        ranges_with_error.printMatrix();
        
        System.out.println("Estimator");
        gdescent_result.estimator.printMatrix();
        
        System.out.println("Full result");
        gdescent_result.estimator_candidate.printMatrix();
        gdescent_result.error.printMatrix();


    }
		
			
	
	
	
	
	
}
