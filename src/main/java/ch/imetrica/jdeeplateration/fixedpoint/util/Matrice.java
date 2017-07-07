package ch.imetrica.jdeeplateration.fixedpoint.util;


import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import ch.imetrica.jdeeplateration.fixedpoint.FPMath;

import java.text.FieldPosition;
import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.StreamTokenizer;



public class Matrice implements Cloneable, java.io.Serializable {

	/* ------------------------
   Class variables
	 * ------------------------ */

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/** Array for internal storage of elements.
   @serial internal array storage.
	 */
	private int[][] A;

	/** Row and column dimensions.
   @serial row dimension.
   @serial column dimension.
	 */
	private int m, n;

	/* ------------------------
   Constructors
	 * ------------------------ */

	/** Construct an m-by-n Matrice of zeros. 
   @param m    Number of rows.
   @param n    Number of colums.
	 */

	public Matrice (int m, int n) {
		this.m = m;
		this.n = n;
		A = new int[m][n];
	}

	/** Construct an m-by-n constant Matrice.
   @param m    Number of rows.
   @param n    Number of colums.
   @param s    Fill the Matrice with this scalar value.
	 */

	public Matrice (int m, int n, int s) {
		this.m = m;
		this.n = n;
		A = new int[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = s;
			}
		}
	}
	
	public Matrice (double[][] A) {
		m = A.length;
		n = A[0].length;
		this.A = new int[m][n];
		for (int i = 0; i < m; i++) {
			if (A[i].length != n) {
				throw new IllegalArgumentException("All rows must have the same length.");
			}
			for(int j = 0 ; j < n; j++){
				this.A[i][j] = FPMath.FloatToFP(A[i][j]);
			}
		}
	}
	
	/** Construct a Matrice from a 2-D array.
   @param A    Two-dimensional array of doubles.
   @exception  IllegalArgumentException All rows must have the same length
   @see        #constructWithCopy
	 */

	public Matrice (int[][] A) {
		m = A.length;
		n = A[0].length;
		for (int i = 0; i < m; i++) {
			if (A[i].length != n) {
				throw new IllegalArgumentException("All rows must have the same length.");
			}
		}
		this.A = A;
	}

	/** Construct a Matrice quickly without checking arguments.
   @param A    Two-dimensional array of doubles.
   @param m    Number of rows.
   @param n    Number of colums.
	 */

	public Matrice (int[][] A, int m, int n) {
		this.A = A;
		this.m = m;
		this.n = n;
	}

	/** Construct a Matrice from a one-dimensional packed array
   @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
   @param m    Number of rows.
   @exception  IllegalArgumentException Array length must be a multiple of m.
	 */

	public Matrice (int vals[], int m) {
		this.m = m;
		n = (m != 0 ? vals.length/m : 0);
		if (m*n != vals.length) {
			throw new IllegalArgumentException("Array length must be a multiple of m.");
		}
		A = new int[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = vals[i+j*m];
			}
		}
	}

	/* ------------------------
   Public Methods
	 * ------------------------ */

	/** Construct a Matrice from a copy of a 2-D array.
   @param A    Two-dimensional array of doubles.
   @exception  IllegalArgumentException All rows must have the same length
	 */

	public static Matrice constructWithCopy(int[][] A) {
		int m = A.length;
		int n = A[0].length;
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			if (A[i].length != n) {
				throw new IllegalArgumentException
				("All rows must have the same length.");
			}
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j];
			}
		}
		return X;
	}

	/** Make a deep copy of a Matrice
	 */

	public Matrice copy () {
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j];
			}
		}
		return X;
	}

	/** Clone the Matrice object.
	 */

	public Object clone () {
		return this.copy();
	}

	/** Access the internal two-dimensional array.
   @return     Pointer to the two-dimensional array of Matrice elements.
	 */

	public int[][] getArray () {
		return A;
	}

	/** Copy the internal two-dimensional array.
   @return     Two-dimensional array copy of Matrice elements.
	 */

	public int[][] getArrayCopy () {
		int[][] C = new int[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j];
			}
		}
		return C;
	}

	/** Make a one-dimensional column packed copy of the internal array.
   @return     Matrice elements packed in a one-dimensional array by columns.
	 */

	public int[] getColumnPackedCopy () {
		int[] vals = new int[m*n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				vals[i+j*m] = A[i][j];
			}
		}
		return vals;
	}

	/** Make a one-dimensional row packed copy of the internal array.
   @return     Matrice elements packed in a one-dimensional array by rows.
	 */

	public int[] getRowPackedCopy () {
		int[] vals = new int[m*n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				vals[i*n+j] = A[i][j];
			}
		}
		return vals;
	}

	/** Get row dimension.
   @return     m, the number of rows.
	 */

	public int getRowDimension () {
		return m;
	}

	/** Get column dimension.
   @return     n, the number of columns.
	 */

	public int getColumnDimension () {
		return n;
	}

	/** Get a single element.
   @param i    Row index.
   @param j    Column index.
   @return     A(i,j)
   @exception  ArrayIndexOutOfBoundsException
	 */

	public int get (int i, int j) {
		return A[i][j];
	}

	/** Get a subMatrice.
   @param i0   Initial row index
   @param i1   Final row index
   @param j0   Initial column index
   @param j1   Final column index
   @return     A(i0:i1,j0:j1)
   @exception  ArrayIndexOutOfBoundsException SubMatrice indices
	 */

	public Matrice getMatrice (int i0, int i1, int j0, int j1) {
		Matrice X = new Matrice(i1-i0+1,j1-j0+1);
		int[][] B = X.getArray();
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = j0; j <= j1; j++) {
					B[i-i0][j-j0] = A[i][j];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("SubMatrice indices");
		}
		return X;
	}

	/** Get a subMatrice.
   @param r    Array of row indices.
   @param c    Array of column indices.
   @return     A(r(:),c(:))
   @exception  ArrayIndexOutOfBoundsException SubMatrice indices
	 */

	public Matrice getMatrice (int[] r, int[] c) {
		Matrice X = new Matrice(r.length,c.length);
		int[][] B = X.getArray();
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < c.length; j++) {
					B[i][j] = A[r[i]][c[j]];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("SubMatrice indices");
		}
		return X;
	}

	/** Get a subMatrice.
   @param i0   Initial row index
   @param i1   Final row index
   @param c    Array of column indices.
   @return     A(i0:i1,c(:))
   @exception  ArrayIndexOutOfBoundsException SubMatrice indices
	 */

	public Matrice getMatrice (int i0, int i1, int[] c) {
		Matrice X = new Matrice(i1-i0+1,c.length);
		int[][] B = X.getArray();
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = 0; j < c.length; j++) {
					B[i-i0][j] = A[i][c[j]];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("SubMatrice indices");
		}
		return X;
	}

	/** Get a subMatrice.
   @param r    Array of row indices.
   @param i0   Initial column index
   @param i1   Final column index
   @return     A(r(:),j0:j1)
   @exception  ArrayIndexOutOfBoundsException SubMatrice indices
	 */

	public Matrice getMatrice (int[] r, int j0, int j1) {
		Matrice X = new Matrice(r.length,j1-j0+1);
		int[][] B = X.getArray();
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = j0; j <= j1; j++) {
					B[i][j-j0] = A[r[i]][j];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("SubMatrice indices");
		}
		return X;
	}

	/** Set a single element.
   @param i    Row index.
   @param j    Column index.
   @param s    A(i,j).
   @exception  ArrayIndexOutOfBoundsException
	 */

	public void set (int i, int j, int s) {
		A[i][j] = s;
	}

	/** Set a subMatrice.
   @param i0   Initial row index
   @param i1   Final row index
   @param j0   Initial column index
   @param j1   Final column index
   @param X    A(i0:i1,j0:j1)
   @exception  ArrayIndexOutOfBoundsException SubMatrice indices
	 */

	public void setMatrice (int i0, int i1, int j0, int j1, Matrice X) {
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = j0; j <= j1; j++) {
					A[i][j] = X.get(i-i0,j-j0);
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("SubMatrice indices");
		}
	}

	/** Set a subMatrice.
   @param r    Array of row indices.
   @param c    Array of column indices.
   @param X    A(r(:),c(:))
   @exception  ArrayIndexOutOfBoundsException SubMatrice indices
	 */

	public void setMatrice (int[] r, int[] c, Matrice X) {
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < c.length; j++) {
					A[r[i]][c[j]] = X.get(i,j);
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("SubMatrice indices");
		}
	}

	/** Set a subMatrice.
   @param r    Array of row indices.
   @param j0   Initial column index
   @param j1   Final column index
   @param X    A(r(:),j0:j1)
   @exception  ArrayIndexOutOfBoundsException SubMatrice indices
	 */

	public void setMatrice (int[] r, int j0, int j1, Matrice X) {
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = j0; j <= j1; j++) {
					A[r[i]][j] = X.get(i,j-j0);
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("SubMatrice indices");
		}
	}

	/** Set a subMatrice.
   @param i0   Initial row index
   @param i1   Final row index
   @param c    Array of column indices.
   @param X    A(i0:i1,c(:))
   @exception  ArrayIndexOutOfBoundsException SubMatrice indices
	 */

	public void setMatrice (int i0, int i1, int[] c, Matrice X) {
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = 0; j < c.length; j++) {
					A[i][c[j]] = X.get(i-i0,j);
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("SubMatrice indices");
		}
	}

	/** Matrice transpose.
   @return    A'
	 */

	public Matrice transpose () {
		Matrice X = new Matrice(n,m);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[j][i] = A[i][j];
			}
		}
		return X;
	}

	/** One norm
   @return    maximum column sum.
	 */

	public int norm1 () {
		int f = 0;
		for (int j = 0; j < n; j++) {
			int s = 0;
			for (int i = 0; i < m; i++) {
				s += Math.abs(A[i][j]);
			}
			f = Math.max(f,s);
		}
		return f;
	}

	/** Two norm
   @return    maximum singular value.
	 */

	public int norm2 () {
		return (new SingularValueDecomposition(this).norm2());
	}

	/** Infinity norm
   @return    maximum row sum.
	 */

	public int normInf () {
		int f = 0;
		for (int i = 0; i < m; i++) {
			int s = 0;
			for (int j = 0; j < n; j++) {
				s += Math.abs(A[i][j]);
			}
			f = Math.max(f,s);
		}
		return f;
	}

	/** Frobenius norm
   @return    sqrt of sum of squares of all elements.
	 */

	public int normF () {
		int f = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				f = Maths.hypot(f,A[i][j]);
			}
		}
		return f;
	}

	/**  Unary minus
   @return    -A
	 */

	public Matrice uminus () {
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = -A[i][j];
			}
		}
		return X;
	}

	/** C = A + B
   @param B    another Matrice
   @return     A + B
	 */

	public Matrice plus (Matrice B) {
		checkMatriceDimensions(B);
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] + B.A[i][j];
			}
		}
		return X;
	}

	/** A = A + B
   @param B    another Matrice
   @return     A + B
	 */

	public Matrice plusEquals (Matrice B) {
		checkMatriceDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] + B.A[i][j];
			}
		}
		return this;
	}

	/** C = A - B
   @param B    another Matrice
   @return     A - B
	 */

	public Matrice minus (Matrice B) {
		checkMatriceDimensions(B);
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] - B.A[i][j];
			}
		}
		return X;
	}

	/** A = A - B
   @param B    another Matrice
   @return     A - B
	 */

	public Matrice minusEquals (Matrice B) {
		checkMatriceDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] - B.A[i][j];
			}
		}
		return this;
	}

	/** Element-by-element multiplication, C = A.*B
   @param B    another Matrice
   @return     A.*B
	 */

	public Matrice arrayTimes (Matrice B) {
		checkMatriceDimensions(B);
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = FPMath.FPMul(A[i][j] , B.A[i][j]);
			}
		}
		return X;
	}

	/** Element-by-element multiplication in place, A = A.*B
   @param B    another Matrice
   @return     A.*B
	 */

	public Matrice arrayTimesEquals (Matrice B) {
		checkMatriceDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = FPMath.FPMul(A[i][j] , B.A[i][j]);
			}
		}
		return this;
	}

	/** Element-by-element right division, C = A./B
   @param B    another Matrice
   @return     A./B
	 */

	public Matrice arrayRightDivide (Matrice B) {
		checkMatriceDimensions(B);
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = FPMath.FPDiv(A[i][j] , B.A[i][j]);
			}
		}
		return X;
	}

	/** Element-by-element right division in place, A = A./B
   @param B    another Matrice
   @return     A./B
	 */

	public Matrice arrayRightDivideEquals (Matrice B) {
		checkMatriceDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = FPMath.FPDiv(A[i][j] , B.A[i][j]);
			}
		}
		return this;
	}

	/** Element-by-element left division, C = A.\B
   @param B    another Matrice
   @return     A.\B
	 */

	public Matrice arrayLeftDivide (Matrice B) {
		checkMatriceDimensions(B);
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = FPMath.FPDiv(B.A[i][j] , A[i][j]);
			}
		}
		return X;
	}

	/** Element-by-element left division in place, A = A.\B
   @param B    another Matrice
   @return     A.\B
	 */

	public Matrice arrayLeftDivideEquals (Matrice B) {
		checkMatriceDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = FPMath.FPDiv(B.A[i][j] , A[i][j]);
			}
		}
		return this;
	}

	/** Multiply a Matrice by a scalar, C = s*A
   @param s    scalar
   @return     s*A
	 */

	public Matrice times (int s) {
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = FPMath.FPMul(s,A[i][j]);
			}
		}
		return X;
	}


	/** Multiply a Matrice by a scalar in place, A = s*A
   @param s    scalar
   @return     replace A by s*A
	 */

	public Matrice timesEquals (int s) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = FPMath.FPMul(s,A[i][j]);
			}
		}
		return this;
	}

	/** Linear algebraic Matrice multiplication, A * B
   @param B    another Matrice
   @return     Matrice product, A * B
   @exception  IllegalArgumentException Matrice inner dimensions must agree.
	 */

	public Matrice times (Matrice B) {
		if (B.m != n) {
			throw new IllegalArgumentException("Matrice inner dimensions must agree.");
		}
		Matrice X = new Matrice(m,B.n);
		int[][] C = X.getArray();
		int[] Bcolj = new int[n];
		for (int j = 0; j < B.n; j++) {
			for (int k = 0; k < n; k++) {
				Bcolj[k] = B.A[k][j];
			}
			for (int i = 0; i < m; i++) {
				int[] Arowi = A[i];
				int s = 0;
				for (int k = 0; k < n; k++) {
					s += FPMath.FPMul(Arowi[k],Bcolj[k]);
				}
				C[i][j] = s;
			}
		}
		return X;
	}

	/** LU Decomposition
   @return     LUDecomposition
   @see LUDecomposition
	 */

	public LUDecomposition lu () {
		return new LUDecomposition(this);
	}

	/** QR Decomposition
   @return     QRDecomposition
   @see QRDecomposition
	 */

	public QRDecomposition qr () {
		return new QRDecomposition(this);
	}

	/** Cholesky Decomposition
   @return     CholeskyDecomposition
   @see CholeskyDecomposition
	 */

	public CholeskyDecomposition chol () {
		return new CholeskyDecomposition(this);
	}

	/** Singular Value Decomposition
   @return     SingularValueDecomposition
   @see SingularValueDecomposition
	 */

	public SingularValueDecomposition svd () {
		return new SingularValueDecomposition(this);
	}

	/** Eigenvalue Decomposition
   @return     EigenvalueDecomposition
   @see EigenvalueDecomposition
	 */
/*
	public EigenvalueDecomposition eig () {
		return new EigenvalueDecomposition(this);
	}
*/
	/** Solve A*X = B
   @param B    right hand side
   @return     solution if A is square, least squares solution otherwise
	 */

	public Matrice solve (Matrice B) {
		return (m == n ? (new LUDecomposition(this)).solve(B) :
			(new QRDecomposition(this)).solve(B));
	}

	/** Solve X*A = B, which is also A'*X' = B'
   @param B    right hand side
   @return     solution if A is square, least squares solution otherwise.
	 */

	public Matrice solveTranspose (Matrice B) {
		return transpose().solve(B.transpose());
	}

	/** Matrice inverse or pseudoinverse
   @return     inverse(A) if A is square, pseudoinverse otherwise.
	 */

	public Matrice inverse () {
		return solve(identity(m,m));
	}

	/** Matrice determinant
   @return     determinant
	 */

	public int det () {
		return new LUDecomposition(this).det();
	}

	/** Matrice rank
   @return     effective numerical rank, obtained from SVD.
	 */

	public int rank () {
		return new SingularValueDecomposition(this).rank();
	}

	/** Matrice condition (2 norm)
   @return     ratio of largest to smallest singular value.
	 */

	public int cond () {
		return new SingularValueDecomposition(this).cond();
	}

	/** Matrice trace.
   @return     sum of the diagonal elements.
	 */

	public int trace () {
		int t = 0;
		for (int i = 0; i < Math.min(m,n); i++) {
			t += A[i][i];
		}
		return t;
	}

	/** Generate Matrice with random elements
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n Matrice with uniformly distributed random elements.
	 */

	public static Matrice random (int m, int n) {
		Matrice A = new Matrice(m,n);
		int[][] X = A.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				X[i][j] = (int)(Math.random()*100);
			}
		}
		return A;
	}

	/** Generate identity Matrice
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n Matrice with ones on the diagonal and zeros elsewhere.
	 */

	public static Matrice identity (int m, int n) {
		Matrice A = new Matrice(m,n);
		int[][] X = A.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				X[i][j] = (i == j ? FPMath.IntToFP(1) : 0);
			}
		}
		return A;
	}


	/** Print the Matrice to stdout.   Line the elements up in columns
	 * with a Fortran-like 'Fw.d' style format.
   @param w    Column width.
   @param d    Number of digits after the decimal.
	 */

	public void print (int w, int d) {
		print(new PrintWriter(System.out,true),w,d); }

	/** Print the Matrice to the output stream.   Line the elements up in
	 * columns with a Fortran-like 'Fw.d' style format.
   @param output Output stream.
   @param w      Column width.
   @param d      Number of digits after the decimal.
	 */

	public void print (PrintWriter output, int w, int d) {
		DecimalFormat format = new DecimalFormat();
		format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
		format.setMinimumIntegerDigits(1);
		format.setMaximumFractionDigits(d);
		format.setMinimumFractionDigits(d);
		format.setGroupingUsed(false);
		print(output,format,w+2);
	}

	/** Print the Matrice to stdout.  Line the elements up in columns.
	 * Use the format object, and right justify within columns of width
	 * characters.
	 * Note that is the Matrice is to be read back in, you probably will want
	 * to use a NumberFormat that is set to US Locale.
   @param format A  Formatting object for individual elements.
   @param width     Field width for each column.
   @see java.text.DecimalFormat#setDecimalFormatSymbols
	 */

	public void print (NumberFormat format, int width) {
		print(new PrintWriter(System.out,true),format,width); }

	// DecimalFormat is a little disappointing coming from Fortran or C's printf.
	// Since it doesn't pad on the left, the elements will come out different
	// widths.  Consequently, we'll pass the desired column width in as an
	// argument and do the extra padding ourselves.

	/** Print the Matrice to the output stream.  Line the elements up in columns.
	 * Use the format object, and right justify within columns of width
	 * characters.
	 * Note that is the Matrice is to be read back in, you probably will want
	 * to use a NumberFormat that is set to US Locale.
   @param output the output stream.
   @param format A formatting object to format the Matrice elements 
   @param width  Column width.
   @see java.text.DecimalFormat#setDecimalFormatSymbols
	 */

	public void print (PrintWriter output, NumberFormat format, int width) {
		output.println();  // start on new line.
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				String s = format.format(A[i][j]); // format the number
				int padding = Math.max(1,width-s.length()); // At _least_ 1 space
				for (int k = 0; k < padding; k++)
					output.print(' ');
				output.print(s);
			}
			output.println();
		}
		output.println();   // end with blank line.
	}

	/** Read a Matrice from a stream.  The format is the same the print method,
	 * so printed matrices can be read back in (provided they were printed using
	 * US Locale).  Elements are separated by
	 * whitespace, all the elements for each row appear on a single line,
	 * the last row is followed by a blank line.
   @param input the input stream.
	 */

	public static Matrice read (BufferedReader input) throws java.io.IOException {
		StreamTokenizer tokenizer= new StreamTokenizer(input);

		// Although StreamTokenizer will parse numbers, it doesn't recognize
		// scientific notation (E or D); however, Double.valueOf does.
		// The strategy here is to disable StreamTokenizer's number parsing.
		// We'll only get whitespace delimited words, EOL's and EOF's.
		// These words should all be numbers, for Double.valueOf to parse.

		tokenizer.resetSyntax();
		tokenizer.wordChars(0,255);
		tokenizer.whitespaceChars(0, ' ');
		tokenizer.eolIsSignificant(true);
		java.util.Vector v = new java.util.Vector();

		// Ignore initial empty lines
		while (tokenizer.nextToken() == StreamTokenizer.TT_EOL);
		if (tokenizer.ttype == StreamTokenizer.TT_EOF)
			throw new java.io.IOException("Unexpected EOF on Matrice read.");
		do {
			v.addElement(Integer.valueOf(tokenizer.sval)); // Read & store 1st row.
		} while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);

		int n = v.size();  // Now we've got the number of columns!
		int row[] = new int[n];
		for (int j=0; j<n; j++)  // extract the elements of the 1st row.
			row[j]=((Integer)v.elementAt(j)).intValue();
		v.removeAllElements();
		v.addElement(row);  // Start storing rows instead of columns.
		while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
			// While non-empty lines
			v.addElement(row = new int[n]);
			int j = 0;
			do {
				if (j >= n) throw new java.io.IOException
				("Row " + v.size() + " is too long.");
				row[j++] = Integer.valueOf(tokenizer.sval).intValue();
			} while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
			if (j < n) throw new java.io.IOException
			("Row " + v.size() + " is too short.");
		}
		int m = v.size();  // Now we've got the number of rows.
		int[][] A = new int[m][];
		v.copyInto(A);  // copy the rows out of the vector
		return new Matrice(A);
	}


	/* ------------------------
   Private Methods
	 * ------------------------ */

	/** Check if size(A) == size(B) **/

	private void checkMatriceDimensions (Matrice B) {
		if (B.m != m || B.n != n) {
			throw new IllegalArgumentException("Matrice dimensions must agree.");
		}
	}

	/*
	 * My methods
	 */
	public Matrice tranposeAndNegate(){
		Matrice X = new Matrice(n,m);
		int[][] A = getArray();
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[j][i] = -1 * A[i][j];
			}
		}
		return X;
	}

	public Array times (int[] a) {
		return this.times(new Array(a));
	}

	public Array times (Array A) {
		int[][] a = this.getArray();
		int[] b = new int[a.length];
		for(int i = 0 ; i < a.length; i++){
			b[i] = 0;
			for(int j = 0; j < a[0].length; j++){
				b[i] += FPMath.FPMul(a[i][j] , A.getElement(j));
			}
		}
		return new Array(b);
	}

	public Matrice subtract (Matrice B) {
		int[][] A = this.getArray();
		int[][] D = B.getArray();
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] - D[i][j];
			}
		}
		return X;
	}

	public static Matrice append(Matrice topleft, Matrice topright, Matrice botleft, Matrice botright){
		int m = topleft.getRowDimension() + botright.getRowDimension();
		int n = topleft.getColumnDimension() + botright.getColumnDimension();
		Matrice rtn = new Matrice(m,n);
		rtn.setMatrice(0, topleft.getRowDimension()-1  , 0, topleft.getColumnDimension()-1, topleft);
		rtn.setMatrice(topleft.getRowDimension(), m-1  , 0, botleft.getColumnDimension()-1, botleft);
		rtn.setMatrice(0, topright.getRowDimension()-1 , topleft.getColumnDimension(), n-1, topright);
		rtn.setMatrice(topleft.getRowDimension(), m-1  , topleft.getColumnDimension(), n-1, botright);
		return rtn;
	}


	public Matrice divide (int s) {
		Matrice X = new Matrice(m,n);
		int[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = FPMath.FPDiv(A[i][j],s);
			}
		}
		return X;
	}
}

