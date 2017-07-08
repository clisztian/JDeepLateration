package ch.imetrica.jdeeplateration.fixedpoint;

public class Primal_ReturnFP {

	private int index;
	private int delta;
	private int[] out_x;
	public Primal_ReturnFP(int index, int delta, int[] out_x){
		this.index = index;
		this.delta = delta;
		this.out_x = out_x;
	}

	public int get_index(){
		return index;
	}

	public int get_delta(){
		return delta;
	}

	public int[] get_out_x(){
		return out_x;
	}
}
