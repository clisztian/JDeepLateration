package ch.imetrica.jdeeplateration.fixedpoint.util;




public class Helper {

	public static int sign(int a){
		if (a > 0)
			return 65536;
		else if (a < 0)
			return -65536;
		return 0;
	}

	public static int[] getIndices(int[] A, int[] indices){
		int[] rtn = new int[indices.length];
		for(int i = 0 ; i < indices.length; i++){
			rtn[i] = A[indices[i]];
		}
		return rtn;
	}

	public static int[] setIndices(int[] A, int[] indices, int[] B){
		int[] rtn = A.clone();
		for(int i = 0 ; i < indices.length; i++){
			rtn[indices[i]] = B[i];
		}
		return rtn;
	}

	public static int[] unionIntArray(int[] A, int[] B, int maxInt){
		boolean[] temp = new boolean[maxInt];
		int repeats = 0;

		if (A == null)
			A = new int[0];
		if (B == null)
			B = new int[0];

		for(int i = 0 ; i < A.length; i++){
			temp[A[i]] = true;
		}
		for(int i = 0 ; i < B.length; i++){
			if (B[i] >= 0 && B[i] < maxInt){
				if (temp[B[i]] == true){
					repeats += 1;
					B[i] = -1;
				}
				if (B[i] != -1)
					temp[B[i]] = true;
			}
		}
		int[] rtn = new int[A.length + B.length - repeats];
		for(int i = 0 ; i < A.length; i++){
			rtn[i] = A[i];
		}
		for(int i = 0 ; i < B.length; i++){
			if(B[i] != -1)
				rtn[i] = B[i];
		}
		return rtn;
	}

	public static int[] find(int a, int[] arr, int type){
		int length = arr.length;
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


	public static int[] line381(Array A1){
		int[] a1 = A1.getArray(); 
		Helper.Quicksort(a1, 0, a1.length); // a1 is sorted
		int[] a2 = A1.getArray(); // a2 is unsorted
		int[] a3 = new int[a2.length];
		
		for(int i = 0 ; i < a2.length; i++){
			for(int j = 0; j < a1.length; j++){
				if (a2[i] == a1[j]){
					a3[i] = j;
					continue;
				}
			}
		}
		return a3;
	}
	
	public static void swap (int A[], int x, int y)
	{
		int temp = A[x];
		A[x] = A[y];
		A[y] = temp;
	}

	// Reorganizes the given list so all elements less than the first are 
	// before it and all greater elements are after it.                   
	public static int partition(int A[], int f, int l)
	{
		double pivot = A[f];
		while (f < l)
		{
			if (A[f] == pivot || A[l] == pivot) 
			{
				//System.out.print("Only distinct doubleegers allowed - C321");
				//System.out.print("students should ignore this if statement");
				break;
			}
			while (A[f] < pivot) f++;
			while (A[l] > pivot) l--;
			swap (A, f, l);
		}
		return f;
	}

	public static void Quicksort(int A[], int f, int l)
	{
		if (f >= l) return;
		int pivot_index = Helper.partition(A, f, l);
		Quicksort(A, f, pivot_index);
		Quicksort(A, pivot_index+1, l);
	}
}
