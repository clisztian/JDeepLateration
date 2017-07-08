package ch.imetrica.jdeeplateration.fixedpoint;

import ch.imetrica.jdeeplateration.fixedpoint.util.Array;
import ch.imetrica.jdeeplateration.fixedpoint.util.Helper;
import ch.imetrica.jdeeplateration.fixedpoint.util.Matrice;
import ch.imetrica.jdeeplateration.fixedpoint.util.SupportSet;

public class HomotopyFP {
	static int N;
	static int K;
	static SupportSet gamma_x;
	static int[] z_x;
	static Array xk_temp;
	static Array del_x_vec;
	static Array pk_temp;
	static Array dk;
	static int epsilon;
	static boolean isNonnegative = true;
	static final int eps = 1;
	static int tolerance = FPMath.FloatToFP(0.001);

	public static Array SolveHomotopy(Matrice A, Array y){
		int lambda = FPMath.FloatToFP(1e-3);
		int maxiter = 100;

		Array xk_1= null;

		K = A.getRowDimension();
		N = A.getColumnDimension();

		z_x = new int[N];
		gamma_x = new SupportSet(N);
		SupportSet gamma_xk = new SupportSet(N);
		Array Primal_constrk = A.tranposeAndNegate().times(y);

		int[] temp;
		int c;
		int i;
		if (isNonnegative){
			temp = Primal_constrk.minArray();
			c = Math.max(-temp[0], 0);
			i = (int)temp[1];
		}
		else{
			temp = Primal_constrk.absArray().maxArray();
			c = temp[0];
			i = (int)temp[1];
		}

		epsilon = c;
		xk_1 = new Array(N);
		gamma_xk.addElement(i);

		int f = FPMath.FPMul(epsilon , xk_1.norm(1)) + FPMath.FPMul( FPMath.FloatToFP(0.5) , y.subtractArray(A.times(xk_1)).norm(2));
		for(int t = 0 ; t < gamma_xk.getLength(); t++){
			z_x[gamma_xk.getElement(t)] = - Helper.sign(Primal_constrk.getElement(gamma_xk.getElement(t)));
		}

		int[] z_xk = z_x;

		int iter = 0;
		int i_delta;
		int delta; 
		int[] out_x = new int[1];
		out_x[0] = -1;
		int old_delta = 0;
		int count_delta_stop = 0;


		Matrice m1       = A.getMatrice(0, A.getRowDimension()-1, gamma_xk.getSet());
		Matrice AtgxAgx  = m1.transpose().times(m1);
		Matrice iAtgxAgx = AtgxAgx.inverse(); 

		Array del_x = null;
		Array x_k = null;
		Matrice Asupported = null;
		Array Agdelx = null;
		int epsilon_old;
		int prev_f;

		Primal_ReturnFP temp1;

		while (iter < maxiter){
			iter += 1;
			gamma_x = gamma_xk;
			z_x = z_xk;
			x_k = xk_1;

			// update on x

			// update direction
			del_x = iAtgxAgx.times(Helper.getIndices(z_x, gamma_x.getSet()));
			del_x_vec = new Array(N);
			del_x_vec = del_x_vec.setIndices(gamma_x.getSet(), del_x.getArray());

			// %dk = A'*(A*del_x_vec);
			Asupported =  A.getMatrice(0, A.getRowDimension()-1, gamma_x.getSet());
			Agdelx = Asupported.times(del_x);
			dk = (A.transpose()).times(Agdelx);

			// CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW.
			pk_temp = Primal_constrk.copy();
			pk_temp.controlMachinePrecision(epsilon, eps*2);

			xk_temp = x_k.copy();
			xk_temp.setIndices(x_k.absArray().find(2*eps, 3), 0);

			temp1 = update_primal(out_x);
			i_delta = temp1.get_index();
			delta   = temp1.get_delta();
			out_x   = temp1.get_out_x();

			if (old_delta < 4 * eps && delta < 4 * eps){
				count_delta_stop += 1;
				if (count_delta_stop >= 500){
					System.out.println("Stuck in a corner");
					break;
				}
			}
			else{
				count_delta_stop = 0;
			}

			old_delta = delta;

			xk_1 = x_k.addArray(del_x_vec.scalarMultArray(delta));
			Primal_constrk = Primal_constrk.addArray(dk.scalarMultArray(delta));
			epsilon_old = epsilon;
			epsilon = epsilon - delta;

			if (epsilon <= lambda){
				xk_1 = x_k.addArray(del_x_vec.scalarMultArray(epsilon_old - lambda));
				break;
			}

			boolean keep_going = true;
			if (delta != 0){
				prev_f = f;
				f = FPMath.FPMul(lambda , xk_1.norm(1)) + FPMath.FPMul(FPMath.FloatToFP(0.5) , y.subtractArray(Asupported.times(xk_1.getElementsFromIndices(gamma_x.getSet()))).norm(2));
				keep_going = Math.abs(FPMath.FPDiv((prev_f-f),prev_f)) > tolerance;
			}
			if (!keep_going)
				break;

			int rowi, colj, n;
			Matrice AtgxAgx_ij, temp_row, temp_col, iAtgxAgx_ij, AtgxAnx;
			Matrice Q11, Q12, Q21, Q12Q21_Q22, iA11, iA11A12, A21iA11, Q11_right;
			int Q22, S;
			int new_x;

			if (out_x != null && out_x[0] != -1){
				// if an element is removed from gamma_x
				int len_gamma = gamma_x.getLength();
				int[] outx_index = gamma_x.find(out_x[0], 2);
				gamma_x.setElement(outx_index[0], gamma_x.getElement(len_gamma-1));
				gamma_x.setElement(len_gamma-1, out_x[0]);
				gamma_x.deleteIndex(len_gamma-1);
				gamma_xk = gamma_x;

				rowi = outx_index[0];
				colj = outx_index[0];

				AtgxAgx_ij =  AtgxAgx.copy();
				temp_row =  AtgxAgx_ij.getMatrice(rowi, rowi, 0, AtgxAgx_ij.getColumnDimension()-1);
				AtgxAgx_ij.setMatrice(rowi, rowi, 0, AtgxAgx_ij.getColumnDimension()-1, AtgxAgx_ij.getMatrice(len_gamma-1, len_gamma-1, 0, AtgxAgx_ij.getColumnDimension()-1));
				AtgxAgx_ij.setMatrice(len_gamma-1, len_gamma-1, 0, AtgxAgx_ij.getColumnDimension()-1, temp_row);

				temp_col =  AtgxAgx_ij.getMatrice(0, AtgxAgx_ij.getRowDimension()-1, colj, colj);
				AtgxAgx_ij.setMatrice(0, AtgxAgx_ij.getRowDimension()-1, colj, colj, AtgxAgx_ij.getMatrice(0, AtgxAgx_ij.getRowDimension()-1, len_gamma-1, len_gamma-1));
				AtgxAgx_ij.setMatrice(0, AtgxAgx_ij.getRowDimension()-1, len_gamma-1, len_gamma-1, temp_col);

				iAtgxAgx_ij =  iAtgxAgx.copy();
				temp_row =  iAtgxAgx_ij.getMatrice(colj, colj, 0, iAtgxAgx_ij.getColumnDimension()-1);
				iAtgxAgx_ij.setMatrice(colj, colj, 0, iAtgxAgx_ij.getColumnDimension()-1, iAtgxAgx_ij.getMatrice(len_gamma-1, len_gamma-1, 0, iAtgxAgx_ij.getColumnDimension()-1));
				iAtgxAgx_ij.setMatrice(len_gamma-1, len_gamma-1, 0, iAtgxAgx_ij.getColumnDimension()-1, temp_row);

				temp_col =  iAtgxAgx_ij.getMatrice(0, iAtgxAgx_ij.getRowDimension()-1, rowi, rowi);
				iAtgxAgx_ij.setMatrice(0, iAtgxAgx_ij.getRowDimension()-1, rowi, rowi, iAtgxAgx_ij.getMatrice(0, iAtgxAgx_ij.getRowDimension()-1, len_gamma-1, len_gamma-1));
				iAtgxAgx_ij.setMatrice(0, iAtgxAgx_ij.getRowDimension()-1, len_gamma-1, len_gamma-1, temp_col);

				AtgxAgx =  AtgxAgx_ij.getMatrice(0, len_gamma-2, 0, len_gamma-2);

				n = AtgxAgx_ij.getRowDimension();
				Q11 =  iAtgxAgx_ij.getMatrice(0, n-2, 0, n-2);
				Q12 =  iAtgxAgx_ij.getMatrice(0, n-2, n-1, n-1);
				Q21 =  iAtgxAgx_ij.getMatrice(n-1, n-1, 0, n-2);
				Q22 = iAtgxAgx_ij.get(n-1, n-1);
				Q12Q21_Q22 =  Q12.times(Q21.divide(Q22));
				iAtgxAgx = Q11.subtract(Q12Q21_Q22);

				xk_1.setIndex(out_x[0], 0);
			}else{
				gamma_xk = gamma_x.copy();
				gamma_xk.addElement(i_delta);
				new_x = i_delta;

				AtgxAnx = A.getMatrice(0, A.getRowDimension()-1, gamma_x.getSet()).transpose().times(A.getMatrice(0, A.getRowDimension()-1, new_x, new_x));
				AtgxAgx = Matrice.append(AtgxAgx, AtgxAnx, AtgxAnx.transpose(), A.getMatrice(0, A.getRowDimension()-1, new_x, new_x).transpose().times(A.getMatrice(0, A.getRowDimension()-1, i_delta, i_delta)));

				n = AtgxAgx.getRowDimension();
				iA11 = iAtgxAgx;
				iA11A12 = iA11.times(AtgxAgx.getMatrice(0, n-2, n-1, n-1));
				A21iA11 = AtgxAgx.getMatrice(n-1, n-1, 0, n-2).times(iA11);
				//int temp2 = AtgxAgx.get(n-1, n-1);
				//int temp3 = AtgxAgx.getMatrice(n-1, n-1, 0, n-2).times(iA11A12).get(0, 0);
				S = AtgxAgx.get(n-1, n-1) - AtgxAgx.getMatrice(n-1, n-1, 0, n-2).times(iA11A12).get(0, 0);
				Q11_right = iA11A12.times(A21iA11.divide(S));

				int[][] tt = new int[1][1];
				tt[0][0] = FPMath.FPDiv(FPMath.IntToFP(1),S);
				iAtgxAgx = Matrice.append(iA11.plus(Q11_right), iA11A12.divide(-1*S), A21iA11.divide(-1*S), new Matrice(tt));

				xk_1.setIndex(i_delta, 0);
			}

			z_xk = new int[N];
			for(int i1 = 0 ; i1 < gamma_xk.getLength(); i1++){
				z_xk[gamma_xk.getElement(i1)] = -1 * Helper.sign(Primal_constrk.getElement(gamma_xk.getElement(i1)));
			}
			for(int i1 = 0; i1 < gamma_x.getLength(); i1++){
				Primal_constrk.setIndex(gamma_x.getElement(i1), FPMath.FPMul(Helper.sign(Primal_constrk.getElement(gamma_x.getElement(i1))),epsilon));
			}
		}
		Array x_out = xk_1;
		//System.out.println("\t Total iterations: " + total_iter);
		return x_out;
	}

	// out_x may need to be an int[] array
	public static Primal_ReturnFP update_primal(int[] out_x){
		SupportSet gamma_lc = SupportSet.setDiffAndUnion(gamma_x, out_x, HomotopyFP.N);

		Array temp = null;

		int delta = 0;
		int idelta = 0;
		int delta1[] = null;
		int[] delta1_pos_ind = null;
		int[] delta2_pos_ind = null;
		int[] delta3_pos_ind = null;
		int delta2[] = null;
		int delta3[] = null;
		int out_x_index;
		if (!isNonnegative){
			//    delta1_constr = (epsilon-pk_temp(gamma_lc))./(1+dk(gamma_lc));
			temp = Array.dotDivide(pk_temp.getElementsFromIndices(gamma_lc.getSet()).subtractConstantFromLeft(epsilon), dk.getElementsFromIndices(gamma_lc.getSet()).addConstant(1));
			//     delta1_pos_ind = find(delta1_constr>0);
			//     delta1_pos = delta1_constr(delta1_pos_ind);
			//     [delta1 i_delta1] = min(delta1_pos);
			//     if isempty(delta1)
			//	       delta1 = inf;
			//     end
			delta1_pos_ind = temp.find(0, 1);
			delta1 = temp.getElementsFromIndices(delta1_pos_ind).minArray();
		}
		/*  delta2_constr = (epsilon+pk_temp(gamma_lc))./(1-dk(gamma_lc));
				delta2_pos_ind = find(delta2_constr>0);
				delta2_pos = delta2_constr(delta2_pos_ind);
				[delta2 i_delta2] = min(delta2_pos);
				if isempty(delta2)
				    delta2 = inf;
				end
		 */
		temp = Array.dotDivide(pk_temp.getElementsFromIndices(gamma_lc.getSet()).addConstant(epsilon), dk.getElementsFromIndices(gamma_lc.getSet()).subtractConstantFromLeft(FPMath.IntToFP(1)));
		delta2_pos_ind = temp.find(0, 1);
		delta2 = temp.getElementsFromIndices(delta2_pos_ind).minArray();

		if (delta1 == null){
			delta = delta2[0];
			idelta = gamma_lc.getElement(delta2_pos_ind[(int)delta2[1]]);
		}
		else if (delta2 == null){
			delta = delta1[0];
			idelta = gamma_lc.getElement(delta1_pos_ind[(int)delta1[1]]);
		}
		else if(delta1[0] > delta2[0]){
			delta = delta2[0];
			idelta = gamma_lc.getElement(delta2_pos_ind[(int)delta2[1]]);
		}
		else{
			delta = delta1[0];
			idelta = gamma_lc.getElement(delta1_pos_ind[(int)delta1[1]]);
		}

		temp = Array.dotDivide(xk_temp.getElementsFromIndices(gamma_x.getSet()).uminus(), del_x_vec.getElementsFromIndices(gamma_x.getSet()));
		delta3_pos_ind = temp.find(0, 1);
		
		Array t_delta = temp.getElementsFromIndices(delta3_pos_ind);
		out_x = null;
		
		if (t_delta != null && t_delta.getLength() != 0){
			delta3 = t_delta.minArray();
			out_x_index = gamma_x.getElement(delta3_pos_ind[(int)delta3[1]]);

			if (delta3 != null && delta3[0] > 0 && delta3[0] <= delta){
				out_x = new int[1];
				delta = delta3[0];
				out_x[0] = out_x_index;
			}
		}
		/*
		xk_1 = xk_temp.addArray(del_x_vec.scalarMultArray(delta));
		xk_1 = xk_1.setIndices(out_x, 0);

		// fix this
		int[] wrong_sign = Array.dotMult(xk_1.getElementsFromIndices(gamma_x.getSet()).sign(), new Array(Helper.getIndices(z_x, gamma_x.getSet()))).find(-1, 2);
		if (isNonnegative){
			wrong_sign = Helper.unionIntArray(wrong_sign, xk_1.getElementsFromIndices(gamma_x.getSet()).find(0, 3), N);
		}
		if (wrong_sign.length > 0){
			delta = 0;
			int[] ind_wrong_x = Helper.line381(del_x_vec.getElementsFromIndices(gamma_x.getElementsFromIndex(wrong_sign)).absArray());
			out_x = gamma_x.getElementsFromIndex(Helper.getIndices(wrong_sign, ind_wrong_x));
		}

		Array temp1 = pk_temp.getElementsFromIndices(gamma_lc.getSet()).addArray(dk.getElementsFromIndices(gamma_lc.getSet()).scalarMultArray(delta)).absArray().subtractConstantFromRight(epsilon - delta);
		int[] i_delta_temp = gamma_lc.getElementsFromIndex(pk_temp.getElementsFromIndices(gamma_lc.getSet()).addArray(dk.getElementsFromIndices(gamma_lc.getSet()).scalarMultArray(delta)).absArray().subtractConstantFromRight(epsilon - delta).find(10*eps, 1));
		if (i_delta_temp.length != 0){
			if (i_delta_temp.length >= 1 && Helper.find(idelta, i_delta_temp, 2).length == 0){
				int i_temp = (int) Array.dotDivide(pk_temp.getElementsFromIndices(i_delta_temp).uminus(), dk.getElementsFromIndices(i_delta_temp)).maxArray()[1];
				idelta = i_delta_temp[i_temp];
				delta = 0;
				out_x = new int[1];
				out_x[0] = -1;
			}
		}
		*/
		return new Primal_ReturnFP(idelta, delta, out_x);
		//return new Homotopy.Primal_Return(idelta, delta, out_x);
	}

}
