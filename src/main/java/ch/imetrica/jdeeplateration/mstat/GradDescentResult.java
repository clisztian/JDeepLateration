package ch.imetrica.jdeeplateration.mstat;

import java.util.Random;

import ch.imetrica.jdeeplateration.anchors.Anchors;
import ch.imetrica.jdeeplateration.cirdata.NavigationChannelList;
import ch.imetrica.jdeeplateration.matrix.Matrix;

public class GradDescentResult {

	
	public final static double C = 299792.458;
	public final static double f0 = 2147000000; //in Hz
	public final static double f_0 = 1;
	
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
	
	
	public static GradDescentResult GradDescent(Anchors anchors, Matrix ranges_in,
            Matrix bounds_in, int n_trial, double alpha, double time_threshold) throws Exception {
		
	
            Random random = new Random();
            Matrix anchors_in = anchors.getAnchors(); 
            
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
                    
                    //System.out.println(error_next + " " + error);
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

	
	
	public static Matrix triDistanceGradient0(Matrix s0, Matrix s1, Matrix estimator, double tdoa_i) throws Exception {
		 
		Matrix num0 = sub(estimator, s0);
		Matrix num1 = sub(estimator, s1);
		
		double diff1 = Mstat.distance(estimator, s1);
		double diff0 = Mstat.distance(estimator, s0);
		double D_i = diff1 - diff0;
		
		double tdoaDiff = D_i - tdoa_i; 
		
		num0.scale(tdoaDiff/(diff0));
		num1.scale(tdoaDiff/(diff1));
				
	    Matrix grad = sub(num1, num0);
	    
		return grad; 		
	}
	
	  
	public static double tdoaEstimate(Matrix s_i, Matrix s_i1, Matrix estimator) throws Exception {
		
		double D_i = Mstat.distance(estimator, s_i1) - Mstat.distance(estimator, s_i);
		return D_i; 
	}
	
	public static double fdoaEstimate(Matrix s_i, Matrix s_i1, Matrix v_i, Matrix v_i1, Matrix estimator) throws Exception {
		
		double diff = 0;
		
		Matrix num0 = sub(estimator, s_i);
		Matrix num1 = sub(estimator, s_i1);
		
		double diffi = Mstat.distance(estimator, s_i);
		double diffi1 = Mstat.distance(estimator, s_i1);
		
		num0.scale(1.0/diffi);
		num1.scale(1.0/diffi1);
		
		double angle_i = num0.dot(v_i);
		double angle_i1 = num1.dot(v_i1);
		
		diff = (angle_i1 - angle_i); 

		return diff; 
	}
	
	
	public static Matrix fdoaGradient0(Matrix s_i, Matrix s_i1, Matrix v_i, Matrix v_i1, 
			Matrix estimator, double fdoa_i) throws Exception {
			
		double D_i = fdoaEstimate(s_i, s_i1, v_i, v_i1, estimator);
		double fdoaDiff = D_i - fdoa_i; 
		
		double diffi = Mstat.distance(estimator, s_i);
		double diffi1 = Mstat.distance(estimator, s_i1);
		
		Matrix numi = sub(estimator, s_i);
		Matrix numi1 = sub(estimator, s_i1);
			
		double angle_i = numi.dot(v_i);
		double angle_i1 = numi1.dot(v_i1);
		
		Matrix gradient = new Matrix(1,3);
		
		for(int i = 0; i < estimator.cols; i++) {
					
			double grad1 = v_i1.w[i]/diffi1 - angle_i1*numi1.w[i]/Math.pow(diffi1, 3);			
			double grad0 = v_i.w[i]/diffi - angle_i*numi.w[i]/Math.pow(diffi, 3);					
			gradient.w[i] = grad1 - grad0;
		}
		
		gradient.scale(fdoaDiff);
		return gradient; 			
	}	
	
	
	
	public static GradDescentResult GradDescentTDOA(Anchors anchors, Matrix tdoas_in,
            Matrix bounds_in, int n_trial, double alpha, double time_threshold, double[] source) throws Exception {
		
	
            Random random = new Random();
            Matrix anchors_in = anchors.getAnchors(); 
            
            int n = anchors_in.rows;
            int dim = anchors_in.cols;
            
            
            GradDescentResult gdescent_result = new GradDescentResult(n_trial, dim);
            
            if (bounds_in == null) {
                bounds_in = new Matrix(1, dim);
            }
            
            Matrix bounds = new Matrix(2, dim);
            
            for(int i = 0; i < dim; i++) {
                bounds.set(0, i, bounds_in.ColumnMin(i));
                bounds.set(1, i, bounds_in.ColumnMax(i));
            }

    
            Matrix sol = new Matrix(source, 1);            
            Matrix ranges = new Matrix(n,1);           
            ranges.w[0] = 0;
            for (int j = 0; j < n-1; j++) {                    	
            	ranges.w[j+1] = tdoaEstimate(anchors_in.getRow(j), anchors_in.getRow(j+1), sol);
            }
            
            
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
                	
                	ranges.w[0] = 0;
                    for (int j = 0; j < n-1; j++) {                    	
                    	ranges.w[j+1] = tdoaEstimate(anchors_in.getRow(j), anchors_in.getRow(j+1), estimator);
                    }
                    
                    double error = Mstat.norm(tdoas_in, ranges);
                    Matrix delta = new Matrix(dim);
                    
                    for (int j = 0; j < n-1; j++) {
                   	
                    	Matrix grad = 
                    			triDistanceGradient0(anchors_in.getRow(j+1), anchors_in.getRow(j), estimator, tdoas_in.w[j+1]);
                    	
                    	delta.add(grad);	
                    }
                    
                    delta.scale(2.0*alpha);
                    Matrix estimator_next = sub(estimator,delta);
                    
                    for (int j = 0; j < n-1; j++) {
                        ranges.w[j+1] = tdoaEstimate(anchors_in.getRow(j), anchors_in.getRow(j+1), estimator);
                    }
                    
                    double error_next = Mstat.norm(tdoas_in, ranges);

                    if (error_next < error) {
                        estimator = estimator_next;
                    }
                    else {
                        gdescent_result.estimator_candidate.setRow(i, estimator);
                        gdescent_result.error.w[i] = error;
                        break;
                    }
                   
                    if (System.nanoTime() - startTime > time_threshold) {
                        gdescent_result.error.w[i] = Double.MAX_VALUE;
                        break;
                    }
                }
            }
            return gdescent_result;
            
	}
	

	public static GradDescentResult GradDescentTDOA0(Anchors anchors, Matrix tdoas_in,
            Matrix bounds_in, int n_trial, double alpha, double time_threshold, double[] source) throws Exception {
		
	
            Random random = new Random();
            Matrix anchors_in = anchors.getAnchors(); 
            
            int n = anchors_in.rows;
            int dim = anchors_in.cols;
            
            
            GradDescentResult gdescent_result = new GradDescentResult(n_trial, dim);
            
            if (bounds_in == null) {
                bounds_in = new Matrix(1, dim);
            }
            
            Matrix bounds = new Matrix(2, dim);
            
            for(int i = 0; i < dim; i++) {
                bounds.set(0, i, bounds_in.ColumnMin(i));
                bounds.set(1, i, bounds_in.ColumnMax(i));
            }

            Matrix ranges = new Matrix(n,1); 
            
//            if(source != null) {
//            	
//            	Matrix sol = new Matrix(source, 1);            
//                
//                ranges.w[0] = 0;
//                for (int j = 0; j < n-1; j++) {                    	
//                	ranges.w[j+1] = tdoaEstimate(anchors_in.getRow(0), anchors_in.getRow(j+1), sol);
//                }
//                
//               // System.out.println("Error at source " + Mstat.norm(tdoas_in, ranges));
//            }
            
            
            
            
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
                	
                	ranges.w[0] = 0;
                    for (int j = 0; j < n-1; j++) {                    	
                    	ranges.w[j+1] = tdoaEstimate(anchors_in.getRow(0), anchors_in.getRow(j+1), estimator);
                    }
                    
                    double error = Mstat.norm(tdoas_in, ranges);
                    Matrix delta = new Matrix(dim);
                    
                    for (int j = 0; j < n-1; j++) {
                   	
                    	Matrix grad = 
                    			triDistanceGradient0(anchors_in.getRow(0), anchors_in.getRow(j+1), estimator, tdoas_in.w[j+1]);
                    	
                    	delta.add(grad);	
                    }
                    
                    delta.scale(2.0*alpha);
                    Matrix estimator_next = sub(estimator,delta);
                    
                    for (int j = 0; j < n-1; j++) {
                        ranges.w[j+1] = tdoaEstimate(anchors_in.getRow(0), anchors_in.getRow(j+1), estimator);
                    }
                    
                    double error_next = Mstat.norm(tdoas_in, ranges);

                    if (error_next < error) {
                        estimator = estimator_next;
                    }
                    else {
                        gdescent_result.estimator_candidate.setRow(i, estimator);
                        gdescent_result.error.w[i] = error;
                        break;
                    }
                   
                    if (System.nanoTime() - startTime > time_threshold) {
                        gdescent_result.error.w[i] = Double.MAX_VALUE;
                        break;
                    }
                }
            }
            return gdescent_result;
            
	}	
	
	
	public static GradDescentResult GradDescentFDOA(Matrix anchors_in, Matrix v, Matrix fdoas_in,
            Matrix bounds_in, int n_trial, double alpha, double time_threshold, double[] source) throws Exception {
		
	
            Random random = new Random();
            
            int n = anchors_in.rows;
            int dim = anchors_in.cols;
            
            
            GradDescentResult gdescent_result = new GradDescentResult(n_trial, dim);
            
            if (bounds_in == null) {
                bounds_in = new Matrix(1, dim);
            }
            
            Matrix bounds = new Matrix(2, dim);
            
            for(int i = 0; i < dim; i++) {
                bounds.set(0, i, bounds_in.ColumnMin(i));
                bounds.set(1, i, bounds_in.ColumnMax(i));
            }

    
            Matrix sol = new Matrix(source, 1);            
            Matrix ranges = new Matrix(n,1);           
            ranges.w[0] = 0;
            for (int j = 0; j < n-1; j++) {                    	
            	ranges.w[j+1] = fdoaEstimate(anchors_in.getRow(j), anchors_in.getRow(j+1), 
            			                     v.getRow(j), v.getRow(j+1), sol);
            }
            
            System.out.println("Error at source: " + Mstat.norm(fdoas_in, ranges));
            
            
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
                	
                	ranges.w[0] = 0;
                    for (int j = 0; j < n-1; j++) {                    	
                    	ranges.w[j+1] = fdoaEstimate(anchors_in.getRow(j), anchors_in.getRow(j+1), 
			                                         v.getRow(j), v.getRow(j+1), estimator);
                    }
                    
                    double error = Mstat.fdoanorm(fdoas_in, ranges);
                    Matrix delta = new Matrix(dim);
                    
                    for (int j = 0; j < n-1; j++) {
                   	
                    	Matrix grad = 
                    			fdoaGradient0(anchors_in.getRow(j), anchors_in.getRow(j+1), 
       			                     v.getRow(j), v.getRow(j+1), estimator, 
       			                     fdoas_in.w[j+1] - fdoas_in.w[j]);
                    	
                    	delta.add(grad);	
                    }
                    
                    delta.scale(2.0*alpha);
                    Matrix estimator_next = sub(estimator,delta);
                    
                    for (int j = 0; j < n-1; j++) {
                        ranges.w[j+1] = fdoaEstimate(anchors_in.getRow(j), anchors_in.getRow(j+1), 
			                                         v.getRow(j), v.getRow(j+1), estimator);
                    }
                    
                    double error_next = Mstat.fdoanorm(fdoas_in, ranges);

                    if (error_next < error) {
                        estimator = estimator_next;
                    }
                    else {
                        gdescent_result.estimator_candidate.setRow(i, estimator);
                        gdescent_result.error.w[i] = error;
                        break;
                    }
                   
                    if (System.nanoTime() - startTime > time_threshold) {
                        gdescent_result.error.w[i] = Double.MAX_VALUE;
                        break;
                    }
                }
            }
            return gdescent_result;
            
	}
	
	
	
	
	
	public static GradDescentResult mlat(Anchors anchors_in, Matrix ranges_in, Matrix bounds_in, int n_trial, double alpha, double time_threshold) throws Exception {
		
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
	
	
	public static GradDescentResult mlatTdoa(Anchors anchors_in, Matrix ranges_in, Matrix bounds_in, 
			int n_trial, double alpha, double time_threshold, double[] source) throws Exception {
		
        GradDescentResult gdescent_result = GradDescentTDOA0(anchors_in, ranges_in, bounds_in, 
        		n_trial, alpha, time_threshold, source);

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
        gdescent_result.error.w[0] = error; 
        gdescent_result.estimator = gdescent_result.estimator_candidate.getRow(idx);
        return gdescent_result;
    }	
	
	
	public static GradDescentResult mlatFdoa(Matrix anchors_in, Matrix v, Matrix fdoas_in, Matrix bounds_in, 
			int n_trial, double alpha, double time_threshold, double[] source) throws Exception {
		
        GradDescentResult gdescent_result = GradDescentFDOA(anchors_in, v, fdoas_in, bounds_in, 
        		n_trial, alpha, time_threshold, source);

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
        gdescent_result.error.w[0] = error; 
        gdescent_result.estimator = gdescent_result.estimator_candidate.getRow(idx);
        return gdescent_result;
    }
	
	
	private static Matrix sub(Matrix estimator2, Matrix delta) throws Exception {
		
		if (estimator2.size != delta.size) {
			throw new Exception("matrix dimension mismatch");
		}
		
		Matrix diff = new Matrix(delta.rows, delta.cols);
		for(int i = 0; i < delta.size; i++) {
			diff.w[i] = estimator2.w[i] - delta.w[i];
		}
		
		return diff;
	}

	public static void testExample() throws Exception {
		
        Random random = new Random();
		
		double error = 0.5;
        double W = 100, L = 100, H = 1;
		
        double W_min = -100;
        double L_min = -100;
        double H_min = 0; 
        
        int n_anchors = 50; 
        Anchors myAnchors = new Anchors(n_anchors);
        for(int i = 0; i < n_anchors; i++) {
        	
        	double cur_W = random.nextDouble() * (W - W_min) + W_min;
        	double cur_L = random.nextDouble() * (L - L_min) + L_min;
        	double cur_H = random.nextDouble() * (H - H_min) + H_min;
        	myAnchors.setCoordinates(cur_W, cur_L, cur_H);     	
        }
        myAnchors.commitCoordinates();
        
        Matrix node = new Matrix(3);
        double[] n = {W*random.nextDouble(), L*random.nextDouble(), H*random.nextDouble()};      
        node.setPointer(n); 
        
        int num_anchors = myAnchors.getNumberOfAnchors();
        Matrix ranges = new Matrix(num_anchors);
           
        
        Matrix ranges_with_error = new Matrix(num_anchors);
        for (int i = 0; i < num_anchors; i++) {
            ranges.w[i] = Mstat.distance(myAnchors.getRow(i), node);
            ranges_with_error.w[i] = ranges.w[i];// + 2 * error * (random.nextDouble() - 0.5);
        }
        
        int n_trial = 10000; 
        double alpha = 0.0001; 
        double time_threshold = 10000;
        Matrix bounds_in = null;
        
        GradDescentResult gdescent_result = GradDescentResult.mlat(myAnchors, ranges_with_error, bounds_in, n_trial, alpha, time_threshold);

        System.out.println("Anchors");
        myAnchors.printMatrix();
        
        System.out.println("Node");
        node.printMatrix();
        
//        System.out.println("Ranges");
//        ranges.printMatrix();
        
        System.out.println("Ranges with error");
        ranges_with_error.printMatrix();
        
        System.out.println("Estimator");
        gdescent_result.estimator.printMatrix();
        
//        System.out.println("Full result");
//        gdescent_result.estimator_candidate.printMatrix();
//        gdescent_result.error.printMatrix();
        
	}
	
	
	public static void main(String[] args) throws Exception {
				
		testExample();
	
    }
		
		
	/*  Algorithm:
	 * 
	 * 
	 *  1) Create NavigationList from log file
	 *  
	 *  2) Compute local coordinate system on navigation list
	 *  
	 *  3) Apply filtering rules to anchors to update NavigationList 
	 *  
	 *  4) Compute tdoas from filtered NavigationList
	 *  
	 *  5) Compute boundaries of local coordinate system where estimate exists
	 *  
	 *  6) Compute estimated solution for sequential additions to navigationList using
	 *     stochastic gradient tdoa
	 *     
	 *  7) Plot trajectory of estimated solutions every sequential addition
	 *  
	 * 
	 * 
	 * 
	 */
	
	
	
	
	
	
	
}
