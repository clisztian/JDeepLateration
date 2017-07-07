package ch.imetrica.jdeeplateration.fixedpoint.util;

import ch.imetrica.jdeeplateration.fixedpoint.FPMath;

public class Array {
	private int[] arr;
	private int length;
	
	public Array(int length){
		arr = new int[length];
		this.length = length;
	}
	
	public Array(int[] a){
		arr = a;
		length = a.length;
	}
	
	public Array(double[] a){
		length = a.length;
		arr = new int[length];
		for(int i = 0 ; i < length; i++){
			arr[i] = FPMath.FloatToFP(a[i]);
		}
	}
		
	public Array copy(){
		return new Array(arr.clone());
	}
	
	public int[] getArray(){
		return arr.clone();
	}
	
	public int getLength(){
		return length;
	}
	
	public int getElement(int m){
		return arr[m];
	}
	
	public Array matrixMultiplyFromLeft(Matrix A){
		int[][] a = A.getArray();
		int[] b = new int[a.length];
		for(int i = 0 ; i < a.length; i++){
			b[i] = 0;
			for(int j = 0; j < a[0].length; j++){
				b[i] += FPMath.FPMul(a[i][j] , arr[j]);
			}
		}
		return new Array(b);
	}
	
	public int norm(int type) {
		int rtn = 0;
		if (type == 1){
			for(int i = 0 ; i < length; i++){
				rtn += Math.abs(arr[i]);
			}
		}
		else if (type == 2){
			for(int i = 0 ; i < length; i++){
				rtn += FPMath.FPMul(arr[i],arr[i]);
			}
		}
		else if (type == 0){
			for(int i = 0 ; i < length; i++){
				if (arr[i] != 0)
					rtn += FPMath.IntToFP(1);
			}
		}
		return rtn;
	}
	
	public int[] minArray(){
		int[] rtn  = new int[2];
		rtn[0] = arr[0];
		rtn[1] = 0;
		for(int i = 0 ; i < length; i++){
			if (arr[i] < rtn[0]){
				rtn[0] = arr[i];
				rtn[1] = i;
			}
		}
		return rtn;
	}
	
	public int[] maxArray(){
		int[] rtn  = new int[2];
		rtn[0] = arr[0];
		rtn[1] = 0;
		for(int i = 0 ; i < length; i++){
			if (arr[i] > rtn[0]){
				rtn[0] = arr[i];
				rtn[1] = i;
			}
		}
		return rtn;
	}
	
	public Array absArray(){
		int[] a = new int[length];
		for(int i = 0 ; i < length; i++){
			a[i] = Math.abs(arr[i]);
		}
		return new Array(a);
	}
	
	public Array addConstant(int s){
		int[] b = new int[length];
		for(int i = 0 ; i < length; i++){
			b[i] = s + arr[i];
		}
		return new Array(b);
	}
	
	public Array addArray(Array a){
		int[] b = new int[length];
		for(int i = 0 ; i < length; i++){
			b[i] = arr[i] + a.getElement(i);
		}
		return new Array(b);
	}
	
	public Array subtractConstantFromLeft(int s){
		int[] b = new int[length];
		for(int i = 0 ; i < length; i++){
			b[i] = s - arr[i];
		}
		return new Array(b);
	}
	
	public Array subtractConstantFromRight(int s){
		int[] b = new int[length];
		for(int i = 0 ; i < length; i++){
			b[i] = arr[i] - s;
		}
		return new Array(b);
	}
	
	public Array subtractArray(Array a){
		int[] b = new int[length];
		for(int i = 0 ; i < length; i++){
			b[i] = arr[i] - a.getElement(i);
		}
		return new Array(b);
	}
	
	public Array scalarMultArray(int a){
		int[] b = new int[length];
		for(int i = 0 ; i < length; i++){
			b[i] = FPMath.FPMul(arr[i],a);
		}
		return new Array(b);
	}
	
	public Array getElementsFromIndices(int[] indices){
		if (indices == null)
			return null;
		int[] rtn = new int[indices.length];
		for(int i = 0 ; i < indices.length; i++){
			rtn[i] = arr[indices[i]];
		}
		return new Array(rtn);
	}
	
	public Array setIndices(int[] indices, int B){
		if (indices  == null)
			return this;
		int[] rtn = arr.clone();
		for(int i = 0 ; i < indices.length; i++){
			rtn[indices[i]] = B;
		}
		return new Array(rtn);
	}
	
	public Array setIndices(int[] indices, int[] B){
		if (indices  == null)
			return this;
		int[] rtn = arr.clone();
		for(int i = 0 ; i < indices.length; i++){
			rtn[indices[i]] = B[i];
		}
		return new Array(rtn);
	}
	
	public void setIndex(int index, int element){
		arr[index] = element;
	}
	
	public int[] find(int a, int type){
		SupportSet rtn = new SupportSet(length);
		if (type == 1){
			for(int i = 0 ; i < length; i++){
				if (arr[i] > a){
					rtn.addElement(i);
				}
			}
		} else if (type == 2){
			for(int i = 0 ; i < length; i++){
				if (arr[i] == a){
					rtn.addElement(i);
				}
			}
		} else if (type == 3){
			for(int i = 0 ; i < length; i++){
				if (arr[i] < a){
					rtn.addElement(i);
				}
			}
		}
		return rtn.getSet();
	}
	
	public Array uminus(){
		int[] rtn = new int[length];
		for(int i = 0 ; i < length; i++){
			rtn[i] = -arr[i];
		}
		return new Array(rtn);
	}
	
	public Array sign(){
		int[] rtn = new int[length];
		for(int i = 0 ; i < length; i++){
			rtn[i] = Helper.sign(arr[i]);
		}
		return new Array(rtn);
	}
	
	public void printArray(){
		for(int i = 0 ; i < length; i++){
			System.out.print(FPMath.FPToFloat(arr[i]) + " ");
		}
		System.out.println();
	}
	
	// specific to line 142 in SolveHomotopy
	public void controlMachinePrecision(int epsilon, int eps){
		int min = Math.min(epsilon, eps);
		for(int i = 0 ; i < arr.length; i++){
			if (Math.abs(Math.abs(arr[i]) - epsilon) < min){
				arr[i] = FPMath.FPMul(Helper.sign(arr[i]) , epsilon);
			}
		}
	}

	// both vectors have the same length
	public static Array dotMult(Array A, Array B){
		int[] rtn = new int[A.length];
		for(int i = 0; i < A.length; i++){
			rtn[i] = FPMath.FPMul(A.getElement(i) , B.getElement(i));
		}
		return new Array(rtn);
	}
	
	// both vectors have the same length
	public static Array dotDivide(Array A, Array B){
		int[] rtn = new int[A.length];
		for(int i = 0; i < A.length; i++){
			rtn[i] = FPMath.FPDiv(A.getElement(i) , B.getElement(i));
		}
		return new Array(rtn);
	}
	
	// normalizing the vector
	public void normalize(){
		int total = 0;
		for(int i = 0 ; i < length; i++){
			total += FPMath.FPMul(arr[i],arr[i]);
		}
		//System.out.println("total: " + total);
		total = FPMath.FPSqrt(total);
		for(int i = 0 ; i < length; i++){
			arr[i] = FPMath.FPDiv(arr[i], total);
		}
		//System.out.println("total: " + total);
	}
	
	/*
	public static int[] findGreatestPositiveMinimum(Array A){
		int[] rtn = new int[2];
		int min = 10000000;
		int index = -1;
		for(int i = 0 ; i < A.length; i++){
			if (A.getElement(i) > 0){
				if (A.getElement(i) < min){
					min = A.getElement(i);
					index = i;
				}
			}
		}
		if (min == 10000000)
			return null;
		rtn[0] = min;
		rtn[1] = (int)index;
		return rtn;
	}
	*/
}

