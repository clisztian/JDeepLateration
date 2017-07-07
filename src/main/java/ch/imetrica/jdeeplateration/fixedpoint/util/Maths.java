package ch.imetrica.jdeeplateration.fixedpoint.util;

import ch.imetrica.jdeeplateration.fixedpoint.FPMath;

public class Maths {


		public static int hypot(int a, int b) {
			int r;
			if (Math.abs(a) > Math.abs(b)){
				r = FPMath.FPDiv(b, a);
				r = FPMath.FPMul(Math.abs(a), FPMath.FPSqrt( FPMath.IntToFP(1) + FPMath.FPMul(r, r)));
			}
			else if (b != 0){
				r = FPMath.FPDiv(a, b);
				r = FPMath.FPMul(Math.abs(b), FPMath.FPSqrt( FPMath.IntToFP(1) + FPMath.FPMul(r, r)));
			}
			else{
				r = FPMath.IntToFP(0);
			}
			return r;
		}


		/** sqrt(a^2 + b^2) without under/overflow. **/



		public static double hypot(double a, double b) {
			double r;
			if (Math.abs(a) > Math.abs(b)) {
				r = b/a;
				r = Math.abs(a)*Math.sqrt(1+r*r);
			} else if (b != 0) {
				r = a/b;
				r = Math.abs(b)*Math.sqrt(1+r*r);
			} else {
				r = 0.0;
			}
			return r;
		}
		
		public static void main(String[] args){
			int a = FPMath.FloatToFP(6.5);
			int b = FPMath.FloatToFP(8.5);
			double c = FPMath.FPToFloat(Maths.hypot(a, b));
			System.out.println(c);
		}
	
	
}
