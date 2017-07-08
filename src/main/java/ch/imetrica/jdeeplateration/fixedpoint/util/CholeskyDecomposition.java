package ch.imetrica.jdeeplateration.fixedpoint.util;

import ch.imetrica.jdeeplateration.fixedpoint.FPMath;

/** Cholesky Decomposition.
   <P>
   For a symmetric, positive definite Matrice A, the Cholesky decomposition
   is an lower triangular Matrice L so that A = L*L'.
   <P>
   If the Matrice is not symmetric or positive definite, the constructor
   returns a partial decomposition and sets an internal flag that may
   be queried by the isSPD() method.
 */

public class CholeskyDecomposition implements java.io.Serializable {

	/* ------------------------
   Class variables
	 * ------------------------ */

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/** Array for internal storage of decomposition.
   @serial internal array storage.
	 */
	private int[][] L;

	/** Row and column dimension (square Matrice).
   @serial Matrice dimension.
	 */
	private int n;

	/** Symmetric and positive definite flag.
   @serial is symmetric and positive definite flag.
	 */
	private boolean isspd;

	/* ------------------------
   Constructor
	 * ------------------------ */

	/** Cholesky algorithm for symmetric and positive definite Matrice.
   @param  A   Square, symmetric Matrice.
   @return     Structure to access L and isspd flag.
	 */

	public CholeskyDecomposition (Matrice Arg) {


		// Initialize.
		int[][] A = Arg.getArray();
		n = Arg.getRowDimension();
		L = new int[n][n];
		isspd = (Arg.getColumnDimension() == n);
		// Main loop.
		for (int j = 0; j < n; j++) {
			int[] Lrowj = L[j];
			int d = 0;
			for (int k = 0; k < j; k++) {
				int[] Lrowk = L[k];
				int s = 0;
				for (int i = 0; i < k; i++) {
					s += FPMath.FPMul(Lrowk[i],Lrowj[i]);
				}
				Lrowj[k] = s = FPMath.FPDiv((A[j][k] - s),L[k][k]);
				d = d + FPMath.FPMul(s,s);
				isspd = isspd & (A[k][j] == A[j][k]); 
			}
			d = A[j][j] - d;
			isspd = isspd & (d > 0);  
			L[j][j] = FPMath.FPSqrt(Math.max(d,0));
			for (int k = j+1; k < n; k++) {
				L[j][k] = 0;
			}
		}
	}

	/* ------------------------
   Temporary, experimental code.
	 * ------------------------ *\

   \** Right Triangular Cholesky Decomposition.
   <P>
   For a symmetric, positive definite Matrice A, the Right Cholesky
   decomposition is an upper triangular Matrice R so that A = R'*R.
   This constructor computes R with the Fortran inspired column oriented
   algorithm used in LINPACK and MATLAB.  In Java, we suspect a row oriented,
   lower triangular decomposition is faster.  We have temporarily included
   this constructor here until timing experiments confirm this suspicion.
	 *\

   \** Array for internal storage of right triangular decomposition. **\
   private transient double[][] R;

   \** Cholesky algorithm for symmetric and positive definite Matrice.
   @param  A           Square, symmetric Matrice.
   @param  rightflag   Actual value ignored.
   @return             Structure to access R and isspd flag.
	 *\

   public CholeskyDecomposition (Matrice Arg, int rightflag) {
      // Initialize.
      double[][] A = Arg.getArray();
      n = Arg.getColumnDimension();
      R = new double[n][n];
      isspd = (Arg.getColumnDimension() == n);
      // Main loop.
      for (int j = 0; j < n; j++) {
         double d = 0.0;
         for (int k = 0; k < j; k++) {
            double s = A[k][j];
            for (int i = 0; i < k; i++) {
               s = s - R[i][k]*R[i][j];
            }
            R[k][j] = s = s/R[k][k];
            d = d + s*s;
            isspd = isspd & (A[k][j] == A[j][k]); 
         }
         d = A[j][j] - d;
         isspd = isspd & (d > 0.0);
         R[j][j] = Math.sqrt(Math.max(d,0.0));
         for (int k = j+1; k < n; k++) {
            R[k][j] = 0.0;
         }
      }
   }

   \** Return upper triangular factor.
   @return     R
	 *\

   public Matrice getR () {
      return new Matrice(R,n,n);
   }

\* ------------------------
   End of temporary code.
	 * ------------------------ */

	/* ------------------------
   Public Methods
	 * ------------------------ */

	/** Is the Matrice symmetric and positive definite?
   @return     true if A is symmetric and positive definite.
	 */

	public boolean isSPD () {
		return isspd;
	}

	/** Return triangular factor.
   @return     L
	 */

	public Matrice getL () {
		return new Matrice(L,n,n);
	}

	/** Solve A*X = B
   @param  B   A Matrice with as many rows as A and any number of columns.
   @return     X so that L*L'*X = B
   @exception  IllegalArgumentException  Matrice row dimensions must agree.
   @exception  RuntimeException  Matrice is not symmetric positive definite.
	 */

	public Matrice solve (Matrice B) {
		if (B.getRowDimension() != n) {
			throw new IllegalArgumentException("Matrice row dimensions must agree.");
		}
		if (!isspd) {
			throw new RuntimeException("Matrice is not symmetric positive definite.");
		}

		// Copy right hand side.
		int[][] X = B.getArrayCopy();
		int nx = B.getColumnDimension();

		// Solve L*Y = B;
		for (int k = 0; k < n; k++) {
			for (int j = 0; j < nx; j++) {
				for (int i = 0; i < k ; i++) {
					X[k][j] -= FPMath.FPMul(X[i][j],L[k][i]);
				}
				X[k][j] = FPMath.FPDiv(X[k][j], L[k][k]);
			}
		}

		// Solve L'*X = Y;
		for (int k = n-1; k >= 0; k--) {
			for (int j = 0; j < nx; j++) {
				for (int i = k+1; i < n ; i++) {
					X[k][j] -= FPMath.FPMul(X[i][j],L[i][k]);
				}
				X[k][j] = FPMath.FPDiv(X[k][j], L[k][k]);
			}
		}


		return new Matrice(X,n,nx);
	}
}


