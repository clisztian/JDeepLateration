package ch.imetrica.jdeeplateration.fixedpoint;



/**
 *
 * 16:16 fixed point math routines, for IAppli/CLDC platform.
 * A fixed point number is a 32 bit int containing 16 bits of integer and 16 bits of fraction.
 *<p>
 * (C) 2001 Beartronics 
 * Author: Henry Minsky (hqm@alum.mit.edu)
 *<p>
 * Licensed under terms "Artistic License"<br>
 * <a href="http://www.opensource.org/licenses/artistic-license.html">http://www.opensource.org/licenses/artistic-license.html</a><br>
 *
 *<p>
 * Numerical algorithms based on
 * http://www.cs.clemson.edu/html_docs/SUNWspro/common-tools/numerical_comp_guide/ncg_examples.doc.html
 * <p>
 * Trig routines based on numerical algorithms described in 
 * http://www.magic-software.com/MgcNumerics.html
 *
 * http://www.dattalo.com/technical/theory/logs.html
 */

public class FPMath {

	/** Convert a 16:16 fixed-point to an int
	 */
	public static int FPToInt (int x) {
		return x>>16;
	}

	/** Convert an int to a 16:16 fixed-point 
	 */

	public static int IntToFP (int x) {
		return x<<16;
	}

	/** Multiply two fixed-point numbers */
	public static int FPMul (int x, int y) {
		long z = (long) x * (long) y;
		return ((int) (z >> 16));
	}

	/** Divides two fixed-point numbers */
	public static int FPDiv (int x, int y) {
		if (y == 0)
			y = 1;
		long z = (((long) x) << 32);
		return (int) ((z / y) >> 16);
	}


	/** Compute square-root of a 16:16 fixed point number */
	public static int FPSqrt (int n) {
		int s = (n + 65536) >> 1;
		for (int i = 0; i < 8; i++) {
			//converge six times
			s = (s + FPDiv(n, s)) >> 1;
		}
		return s;
	}

	/** Round to nearest fixed poitn integer */
	public static int FPRound (int n) {
		if (n > 0) {
			if ((n & 0x8000) != 0) {
				return (((n+0x10000)>>16)<<16);
			} else {
				return (((n)>>16)<<16);
			}
		} else {
			int k;
			n = -n;
			if ((n & 0x8000) != 0) {
				k = (((n+0x10000)>>16)<<16);
			} else {
				k = (((n)>>16)<<16);
			}
			return -k;
		}
	}

	public static final int PI = 205887;
	public static final int E = FloatToFP(2.718281828459, "E");
	public static final int HALF = FloatToFP(0.5, "HALF");
	public static final int HALF_PI = PI/2;
	public static final int QUARTER_PI = PI/4;

	/**
	 * For the inverse tangent calls, all approximations are valid for |t| <= 1.
	 * To compute ATAN(t) for t > 1, use ATAN(t) = PI/2 - ATAN(1/t).  For t < -1,
	 * use ATAN(t) = -PI/2 - ATAN(1/t).
	 */

	static final int SK1 = 498;
	static final int SK2 = 10882;

	/** Computes SIN(f), f is a fixed point number in radians.
	 * 0 <= f <= PI/2
	 */
	public static int FPSin (int f) {
		int sqr = FPMul(f,f);
		int result = SK1;
		result = FPMul(result, sqr);
		result -= SK2;
		result = FPMul(result, sqr);
		result += (1<<16);
		result = FPMul(result, f);
		return result;
	}

	static final int CK1 = 2428;
	static final int CK2 = 32551;

	/** Computes COS(f), f is a fixed point number in radians.
	 * 0 <= f <= PI/2
	 */
	public static int FPCos (int f) {
		int sqr = FPMul(f,f);
		int result = CK1;
		result = FPMul(result, sqr);
		result -= SK2;
		result = FPMul(result, sqr);
		result += (1<<16);
		return result;
	}


	/** Computes Tan(f), f is a fixed point number in radians.
	 * 0 <= f <= PI/4
	 */

	static final int TK1 = 13323;
	static final int TK2 = 20810;

	public static int FPTan (int f) {
		int sqr = FPMul(f,f);
		int result = TK1;
		result = FPMul(result, sqr);
		result += TK2;
		result = FPMul(result, sqr);
		result += (1<<16);
		result = FPMul(result, f);
		return result;
	}


	/** Computes ArcTan(f), f is a fixed point number in radians.
	 * |f| <= 1
	 */

	static final int ATK1 = 18350;

	public static int FPArcTan (int f) {
		int sqr = FPMul(f,f);
		return FPDiv(f, (FPMul(ATK1, sqr) + (1<<16)));
	}

	static final int AS1 = -1228;
	static final int AS2 = 4866;
	static final int AS3 = 13901;
	static final int AS4 = 102939;

	/** Compute ArcSin(f), 0 <= f <= 1
	 */

	public static int ArcSin (int f) {
		int fRoot = FPSqrt((1<<16)-f);
		int result = AS1;
		result = FPMul(result, f);
		result += AS2;
		result = FPMul(result, f);
		result -= AS3;
		result = FPMul(result, f);
		result += AS4;
		result = HALF_PI - (FPMul(fRoot,result));
		return result;
	}


	/** Compute ArcCos(f), 0 <= f <= 1
	 */

	public static int ArcCos (int f) {
		int fRoot = FPSqrt((1<<16)-f);
		int result = AS1;
		result = FPMul(result, f);
		result += AS2;
		result = FPMul(result, f);
		result -= AS3;
		result = FPMul(result, f);
		result += AS4;
		result = FPMul(fRoot,result);
		return result;
	}


	/** Exponential

	e^x = 1 + x/1! + x^2/2! + x^3/3!

	 */

	static int fpfact[] = { 1<<16,
		1<<16,
		2<<16,
		6<<16,
		24<<16,
		120<<16,
		720<<16,
		5040<<16,
		40320<<16
	};

	public static int FPExp (int x) {
		int result = 1<<16;
		int x2 = FPMul(x,x);
		int x3 = FPMul(x2,x);
		int x4 = FPMul(x2,x2);
		int x5 = FPMul(x4,x);
		int x6 = FPMul(x4,x2);
		int x7 = FPMul(x6,x);
		int x8 = FPMul(x4,x4);
		return result + x 
		+ FPDiv(x2,fpfact[2]) 
		+ FPDiv(x3,fpfact[3]) 
		+ FPDiv(x4,fpfact[4])
		+ FPDiv(x5,fpfact[5]) 
		+ FPDiv(x6,fpfact[6]) 
		+ FPDiv(x7,fpfact[7])
		+ FPDiv(x8,fpfact[8]);
	}


	/** Logarithms: 
	 * 
	 * (2) Knuth, Donald E., "The Art of Computer Programming Vol 1",
	 * Addison-Wesley Publishing Company, ISBN 0-201-03822-6 ( this
	 * comes from Knuth (2), section 1.2.3, exercise 25).
	 *
	 * http://www.dattalo.com/technical/theory/logs.html
	 *
           (loop for k from 0 to 31 do
           (setq z (log (+ 1 (expt 2.0 (- (+ k 1)))))) 
           (insert (format "FloatToFP(%f,\"%f\"),\n"  z z)))
	 */

	/*
	static int log2arr[] = {
		FloatToFP(0.405465,"0.405465"),
		FloatToFP(0.223144,"0.223144"),
		FloatToFP(0.117783,"0.117783"),
		FloatToFP(0.060625,"0.060625"),
		FloatToFP(0.030772,"0.030772"),
		FloatToFP(0.015504,"0.015504"),
		FloatToFP(0.007782,"0.007782"),
		FloatToFP(0.003899,"0.003899"),
		FloatToFP(0.001951,"0.001951"),
		FloatToFP(0.000976,"0.000976"),
		FloatToFP(0.000488,"0.000488"),
		FloatToFP(0.000244,"0.000244"),
		FloatToFP(0.000122,"0.000122"),
		FloatToFP(0.000061,"0.000061"),
		FloatToFP(0.000031,"0.000031"),
		FloatToFP(0.000015,"0.000015"),
		FloatToFP(0.000008,"0.000008"),
		FloatToFP(0.000004,"0.000004"),
	};
*/
	/*
      Binary Logarithm:

      case is very similar to the previous one. The only difference is
      in how the input is factored. Like before we are given:

      Input: 16 bit unsigned integer x; 0 < x < 65536

      (or 8 bit unsigned integer...)

      Output: g, lg(x) the logarithm of x with respect to base 2.

      3(b).i) Create a table of logarithms of the following constants:
      log2arr[i] = lg(1 + 2^(-(i+1))) 

      i = 0..M, M == desired size of the table.

      The first few values of the array are

      lg(3/2), lg(5/4), lg(9/8), lg(17/16),...

      Recall that in the previous case the factors were

      lg(2/1), lg(4/3), lg(8/7), lg(16/15),...

      Again, if you wish to compute logarithms to a different base,
      then substitute the lg() function with the appropriate based
      logarithm function.

      3(b).ii)Scale y to a value between 1 and 2.
      This is identical to 3(a).ii.

      3(b).iii) Changing Perspective
      Again, this is identical to 3(a).iii.

      3(b).iv) Factor y.

      This is very similar to step (iv) above. However, we now have
      different factors. Using the same example, x = 1.9, we can find
      the factors for this case.

      a) 1.9 > 1.5,
      x = x/1.5 ==> 1.266666
      b) 1.26 >  1.25
      x ==> 1.0133333
      c) 1.0133 < 1.125 so don't divide
      d) 1.0133 < 1.0625  "     "
      e) etc.


      So, x ~= 1.5 * 1.25 * etc.


      Like the previouse case, these factors are not perfect. Also,
      they're somewhat redundant in the sense that
      1.5*1.25*1.125*1.0625*... spans a range that is larger than 2 (
      ~2.38423 for i<=22). So unlike the previous factoring method,
      this one will not have repeated factors. Here's some psuedo
      code:

      for(i=1,d=0.5; i<M; i++, d/=2)
      if( x > 1+d)
      {
      x /= (1+d);
      g += log2arr[i-1];   // log2arr[i-1] = log2(1+d);
      }


      Here, d takes on the values of 0.5, 0.25, 0.125, ... , 2^(-i). Then
      1+d is the trial factor at each step. If x is greater than this trial
      factor, then we divide the trial factor out and add to g (ultimately
      the logarithm of x) the partial logarithm of the factor.

	 */
/*
	public static int FPLn (int x) {
		int g = 0;
		int d = HALF;
		for (int i = 1; i < 16; i++) {
			if (x > ((1<<16) + d)) {
				x = FPDiv(x, ( (1<<16) + d));
				g += log2arr[i-1];   // log2arr[i-1] = log2(1+d);
			}
			d >>= 1;
		}
		return g;
	}

*/
	// The x,y point where two lines intersect
	static int xint;
	static int yint;

	/**
	 * Does line segment A intersection line segment B?
	 *
	 * Assumes 16 bit fixed point numbers with 16 bits of fraction.
	 *
	 * For debugging, side effect xint, yint, the intersection point.
	 *
	 * <pre>
	 * Algorithm 
	 * 
	 * As an example of algorithm development, consider the intersection of
	 * two line segments.  Given line segment A goes from point XA1 and YA1
	 * to point XA2 and YA2 and given line segment B goes from point XB1 and
	 * YB1 to point XB2 and YB2.  Find whether there is zero, one, or an
	 * infinite number of points of intersection (the line segments overlap)
	 * and the values of the points of intersection.  Assume all numbers are
	 * double.
	 * 
	 * For case 1 where line segment A is not vertical, line segment B is not
	 * vertical, and line segment A is not parallel to line segment B, the
	 * equations for line segment A and B are:
	 * 
	 * 
	 * XMA = (YA2-YA1)/(XA2-XA1) = slope of line segment A
	 * XBA = YA1 - XA1*XMA = Y-intercept for line segment A
	 * YA = XMA*XA + XBA
	 * 
	 * XMB = (YB2-YB1)/(XB2-XB1) = slope of line segment B
	 * XBB = YB1 - XB1*XMB = Y-intercept for line segment B
	 * YB = XMB*XB + XBB
	 * 
	 * At the intersection of line segment A and B, XA=XB=XINT and YA=YB=YINT.
	 * YINT = XMA*XINT + XBA
	 * YINT = XMB*XINT + XBB
	 * XMA*XINT + XBA = XMB*XINT + XBB
	 * XMA*XINT - XMB*XINT = XBB - XBA
	 * XINT*(XMA-XMB) = XBB - XBA
	 * XINT = (XBB-XBA)/(XMA-XMB)
	 * YINT = XMA*XINT + XBA
	 * There is one point of intersection.
	 * 
	 * For case 2 where line segment A is vertical (XA1 is close to XA2) and
	 * line segment B is not vertical, the equations for line segment A and B
	 * are:
	 * 
	 * XA = 0.5*(XA1+XA2)
	 * 
	 * XMB = (YB2-YB1)/(XB2-XB1) = slope of line segment B
	 * XBB = YB1 - XB1*XMB = Y-intercept for line segment B
	 * YB = XMB*XB + XBB
	 * 
	 * At the intersection of line segment A and B, XA=XB=XINT and YA=YB=YINT.
	 * XINT = XA
	 * YINT = XMB*XINT + XBB
	 * There is one point of intersection.
	 * 
	 * For case 3 where line segment A is not vertical and line segment B is
	 * vertical (XB1 is close to XB2), the equations for line segment A and B
	 * are:
	 * 
	 * XMA = (YA2-YA1)/(XA2-XA1) = slope of line segment A
	 * XBA = YA1 - XA1*XMA = Y-intercept for line segment A
	 * YA = XMA*XA + XBA
	 * 
	 * XB= 0.5*(XB1+XB2)
	 * 
	 * At the intersection of line segment A and B, XA=XB=XINT and YA=YB=YINT.
	 * XINT = XB
	 * YINT = XMA*XINT + XBA
	 * There is one point of intersection.
	 * 
	 * For case 4 where line segment A is vertical (XA1 is close to XA2) and
	 * line segment B is vertical (XB1 is close to XB2), the distance between
	 * the parallel line segments is:
	 * 
	 * DIST = ABS ( 0.5*(XA1+XA2) - 0.5*(XB1+XB2) )
	 * 
	 * If DIST is close to zero, then:
	 * 
	 * XINT1 = 0.5*(0.5*(XA1+XA2)+0.5*(XB1+XB2))
	 * YINT1 = MAX(MIN(YA1,YA2),MIN(YB1,YB2))
	 * XINT2 = XINT1
	 * YINT2 = MIN(MAX(YA1,YA2),MAX(YB1,YB2))
	 * There are two points of intersection.
	 * 
	 * For case 5 where line segment A is not vertical, line segment B is not
	 * vertical, and line segment A is parallel to line segment B (XMA is
	 * close to XMB), the equations for line segment A and B are:
	 * 
	 * XMA = (YA2-YA1)/(XA2-XA1) = slope of line segment A
	 * XBA = YA1 - XA1*XMA = Y-intercept for line segment A
	 * YA = XMA*XA + XBA
	 * 
	 * XMB = (YB2-YB1)/(XB2-XB1) = slope of line segment B
	 * XBB = YB1 - XB1*XMB = Y-intercept for line segment B
	 * YB = XMB*XB + XBB
	 * 
	 * The distance between the parallel line segments is:
	 * 
	 * DIST = ABS(XBA-XBB)*COS(ATAN(0.5*(XMA+XMB)))
	 * 
	 * If DIST is close to zero, then:
	 * 
	 * XINT1 = MAX(MIN(XA1,XA2),MIN(XB1,XB2))
	 * YINT1 = MAX(MIN(YA1,YA2),MIN(YB1,YB2))
	 * XINT2 = MIN(MAX(XA1,XA2),MAX(XB1,XB2))
	 * YINT2 = MIN(MAX(YA1,YA2),MAX(YB1,YB2))
	 * There are two points of intersection.
	 * 
	 * After the point or points of intersection are calculated, each
	 * solution must be checked to ensure that the point of intersection lies
	 * on line segment A and B by checking if XINT >= MIN(XA1,XA2) and XINT
	 * <= MAX(XA1,XA2) and YINT >= MIN(YA1,YA2) and YINT <= MAX(YA1,YA2) and
	 * checking if XINT >= MIN(XB1,XB2) and XINT <= MAX(XB1,XB2) and YINT >=
	 * MIN(XB1,XB2) and YINT <= MAX(YB1,YB2).
	 * 
	 * Note that case 2, 3, 4, and 5 are all special instances of case 1
	 * where a division by zero would have caused the creation of an infinite
	 * number and thus a program error.
	 * 
	 * </pre>
	 */
	public static boolean intersects (int ax0, int ay0, int ax1, int ay1,
			int bx0, int by0, int bx1, int by1) {

		ax0 <<= 16;
		ay0 <<= 16;
		ax1 <<= 16;
		ay1 <<= 16;

		bx0 <<= 16;
		by0 <<= 16;
		bx1 <<= 16;
		by1 <<= 16;

		int adx = (ax1 - ax0);
		int ady = (ay1 - ay0);
		int bdx = (bx1 - bx0);
		int bdy = (by1 - by0);

		int xma;
		int xba;

		int xmb;
		int xbb;	
		int TWO = (2 << 16);

		if ((adx == 0) && (bdx == 0)) { // both vertical lines
			int dist = Math.abs(FPDiv((ax0+ax1)-(bx0+bx1), TWO));
			return (dist == 0);
		} else if (adx == 0) { // A  vertical
			int xa = FPDiv((ax0 + ax1), TWO);
			xmb = FPDiv(bdy,bdx);           // slope segment B
			xbb = by0 - FPMul(bx0, xmb); // y intercept of segment B
			xint = xa;
			yint = (FPMul(xmb,xint)) + xbb;
		} else if ( bdx == 0) { // B vertical
			int xb = FPDiv((bx0+bx1), TWO);
			xma = FPDiv(ady,adx);           // slope segment A
			xba = ay0 - (FPMul(ax0,xma)); // y intercept of segment A
			xint = xb;
			yint = (FPMul(xma,xint)) + xba;
		} else {
			xma = FPDiv(ady,adx);           // slope segment A
			xba = ay0 - (FPMul(ax0, xma)); // y intercept of segment A

			xmb = FPDiv(bdy,bdx);           // slope segment B
			xbb = by0 - (FPMul(bx0,xmb)); // y intercept of segment B

			// parallel lines? 
			if (xma == xmb) {
				// Need trig functions
				int dist = Math.abs(FPMul((xba-xbb),
						(FPCos(FPArcTan(FPDiv((xma+xmb), TWO))))));
				if (dist < (1<<16) ) {
					return true;
				} else {
					return false;
				}
			} else {
				// Calculate points of intersection
				// At the intersection of line segment A and B, XA=XB=XINT and YA=YB=YINT.
				if ((xma-xmb) == 0) {
					return false;
				}
				xint = FPDiv((xbb-xba),(xma-xmb));
				yint = (FPMul(xma,xint)) + xba;
			}
		}

		// After the point or points of intersection are calculated, each
		// solution must be checked to ensure that the point of intersection lies
		// on line segment A and B.

		int minxa = Math.min(ax0, ax1);
		int maxxa = Math.max(ax0, ax1);

		int minya = Math.min(ay0, ay1);
		int maxya = Math.max(ay0, ay1);

		int minxb = Math.min(bx0, bx1);
		int maxxb = Math.max(bx0, bx1);

		int minyb = Math.min(by0, by1);
		int maxyb = Math.max(by0, by1);

		return ((xint >= minxa) && (xint <= maxxa) && (yint >= minya) && (yint <= maxya) 
				&& 
				(xint >= minxb) && (xint <= maxxb) && (yint >= minyb) && (yint <= maxyb));
	}


	////////////////////////////////////////////////////////////////
	static double pi = 3.141592;

	/*
	public static void main (String args[]) {
		debug("0 * 1 = " + FPToInt(FPMul(IntToFP(0), IntToFP(1))));
		debug("10 * 1 = " + FPToInt(FPMul(IntToFP(10), IntToFP(1))));
		debug("100 * 100 = " + FPToInt(FPMul(IntToFP(100), IntToFP(100))));
		debug("-2 * 3 = " + FPToInt(FPMul(IntToFP(-2), IntToFP(3))));
		debug("-6 * -7 = " + FPToInt(FPMul(IntToFP(-6), IntToFP(-7))));
		debug("-6 * 7 = " + FPToInt(FPMul(IntToFP(-6), IntToFP(7))));
		debug("259 / 22 = " + FPToInt(FPMul(IntToFP(259), IntToFP(22))));
		debug("259 / -22 = " + FPToInt(FPMul(IntToFP(259), IntToFP(-22))));
		debug("sqrt(2537) = " + FPToFloat(FPSqrt(IntToFP(2537))));
		debug("round(20.7) = " + FPToFloat(FPRound(FloatToFP(20.7))));
		debug("round(200.5) = " + FPToFloat(FPRound(FloatToFP(200.5))));
		debug("round(827.3) = " + FPToFloat(FPRound(FloatToFP(827.3))));
		debug("round(-827.3) = " + FPToFloat(FPRound(FloatToFP(-827.3))));
		debug("round(-827.7) = " + FPToFloat(FPRound(FloatToFP(-827.7))));
		debug("sin(30) = " + FPToFloat(FPSin(FloatToFP(30*(pi/180)))));
		debug("sin(60) = " + FPToFloat(FPSin(FloatToFP(60*(pi/180)))));
		debug("sin(45) = " + FPToFloat(FPSin(FloatToFP(45*(pi/180)))));
		debug("exp(1) = " + FPToFloat(FPExp(FloatToFP(1))));
		debug("exp(.5) = " + FPToFloat(FPExp(FloatToFP(0.5))));
		debug("exp(2) = " + FPToFloat(FPExp(FloatToFP(2))));
		debug("exp(5) = " + FPToFloat(FPExp(FloatToFP(5))));
		debug("exp(7.27) = " + FPToFloat(FPExp(FloatToFP(7.27))));
		debug("exp(-3) = " + FPToFloat(FPExp(FloatToFP(-3))));
		debug("ln(e) = " + FPToFloat(FPLn(FloatToFP(E))));
		debug("ln(e^2) = " + FPToFloat(FPLn(FloatToFP(2.718281828459 * 2.718281828459))));
		debug("ln(10) = " + FPToFloat(FPLn(FloatToFP(10))));
		debug("ln(1.4) = " + FPToFloat(FPLn(FloatToFP(1.4))));
		debug("ln(1.9) = " + FPToFloat(FPLn(FloatToFP(1.9))));
		debug("ln(1.2319) = " + FPToFloat(FPLn(FloatToFP(1.2319))));
		debug("ln(2.1) = " + FPToFloat(FPLn(FloatToFP(2.1))));
		debug("ln(6.1) = " + FPToFloat(FPLn(FloatToFP(6.1))));
		debug("ln(180000) = " + FPToFloat(FPLn(FloatToFP(180000))));

	}*/
	////////////////////////////////////////////////////////////////


	/** Convert a fixed point number to a float */
	public static double FPToFloat (int x) {
		return (x >> 16) + ( ((double) (x & 0xFFFF)) / 65536.0);
	}

	public static int FloatToFP (double z, String desc) {
		int n = FloatToFP(z);
		debug(desc+" = "+n);
		return n;
	}

	public static int FloatToFP (double z) {
		int n = ((((int) Math.floor (z))) << 16) + 
		((int)( ((z - Math.floor (z))) * 65536));
		return n;
	}


	static void debug (Object s) {
		System.out.println(s);
	}

	static String binary (int n) {
		StringBuffer buf = new StringBuffer();
		for (int i = 31; i >= 0; i--) {
			buf.append(""+ ((n>>i)&1));
			if (i%8 == 0) {
				buf.append(" ");
			}
		}
		return buf.toString();
	}


}
