package ch.imetrica.jdeeplateration.fixedpoint.util;



public class SupportSet {

	private int tail;
	private int length;
	private int[] set;

	private boolean changed = true;
	private int[] rtn;

	public SupportSet(int length){
		this.set = new int[length];
		this.length = length;
		this.tail = 0;
	}

	public SupportSet(int[] set, int tail){
		this.set    = set;
		this.length = set.length;
		this.tail   = tail;
	}

	public SupportSet copy(){
		return new SupportSet(this.set.clone(), this.tail);
	}

	public int getLength(){
		return tail;
	}

	public int[] getSet(){
		if (changed == false)
			return rtn;
		else{
			rtn = new int[tail];
			for(int i = 0; i < tail; i++){
				rtn[i] = set[i];
			}
			changed = false;
			return rtn;
		}
	}

	public void addElement(int x){
		set[tail] = x;
		tail = tail + 1;
		changed = true;
	}

	public void addElements(int[] x){
		for(int i = 0 ; i < x.length; i++){
			set[tail + i] = x[i];
		}
		tail = tail + x.length;
		changed = true;
	}

	public int getElement(int index){
		if (index < tail){
			return set[index];
		}
		return -1;
	}

	public int[] getElementsFromIndex(int[] index){
		int[] rtn = new int[index.length];
		for(int i = 0; i < index.length; i++){
			rtn[i] = set[index[i]];
		}
		return rtn;
	}

	public void setElement(int index, int value){
		set[index] = value;
		changed = true;
	}

	public void deleteIndex(int index){
		int t       = set[index];
		set[index]  = set[tail-1];
		set[tail-1] = t;
		tail        = tail - 1;
		changed = true;
	}

	public int[] find(int a, int type){
		SupportSet rtn = new SupportSet(length);
		if (type == 1){
			for(int i = 0 ; i < length; i++){
				if (set[i] > a){
					rtn.addElement(i);
				}
			}
		} else if (type == 2){
			for(int i = 0 ; i < length; i++){
				if (set[i] == a){
					rtn.addElement(i);
				}
			}
		} else if (type == 3){
			for(int i = 0 ; i < length; i++){
				if (set[i] < a){
					rtn.addElement(i);
				}
			}
		}
		return rtn.getSet();
	}


	public static SupportSet setDiffAndUnion(SupportSet X, int[] B, int N){
		boolean[] temp = new boolean[N];
		for(int i = 0 ; i < X.getLength(); i++){
			temp[X.getElement(i)] = true;
		}
		if(B != null){
			for(int i = 0 ; i < B.length; i++){
				if (B[i] >= 0 && B[i] < N){
					if (temp[B[i]] == true) {
					}
					temp[B[i]] = true;
				}
			}
		}
		SupportSet rtn;
		rtn = new SupportSet(N - X.getLength());
		for(int i = 0 ; i < N ; i++){
			if (temp[i] == false)
				rtn.addElement(i);
		}
		return rtn;
	}

	public void printSet(){
		for(int i = 0 ; i < tail; i++)
			System.out.print(set[i] + " ");
		System.out.println();
	}
}
