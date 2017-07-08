package ch.imetrica.jdeeplateration.fixedpoint.util;

import ch.imetrica.jdeeplateration.fixedpoint.FPMath;

/** QR Decomposition.
<P>
   For an m-by-n Matrice A with m >= n, the QR decomposition is an m-by-n
   orthogonal Matrice Q and an n-by-n upper triangular Matrice R so that
   A = Q*R.
<P>
   The QR decompostion always exists, even if the Matrice does not have
   full rank, so the constructor will never fail.  The primary use of the
   QR decomposition is in the least squares solution of nonsquare systems
   of simultaneous linear equations.  This will fail if isFullRank()
   returns false.
 */

public class QRDecomposition implements java.io.Serializable {

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
	private int[][] QR;

	/** Row and column dimensions.
   @serial column dimension.
   @serial row dimension.
	 */
	private int m, n;

	/** Array for internal storage of diagonal of R.
   @serial diagonal of R.
	 */
	private int[] Rdiag;

	/* ------------------------
   Constructor
	 * ------------------------ */

	/** QR Decomposition, computed by Householder reflections.
   @param A    Rectangular Matrice
   @return     Structure to access R and the Householder vectors and compute Q.
	 */

	public QRDecomposition (Matrice A) {
		// Initialize.
		QR = A.getArrayCopy();
		m = A.getRowDimension();
		n = A.getColumnDimension();
		Rdiag = new int[n];

		// Main loop.
		for (int k = 0; k < n; k++) {
			// Compute 2-norm of k-th column without under/overflow.
			int nrm = 0;
			for (int i = k; i < m; i++) {
				nrm = Maths.hypot(nrm,QR[i][k]);
			}

			if (nrm != 0) {
				// Form k-th Householder vector.
				if (QR[k][k] < 0) {
					nrm = -nrm;
				}
				for (int i = k; i < m; i++) {
					QR[i][k] = FPMath.FPDiv(QR[i][k],nrm);
				}
				QR[k][k] += FPMath.IntToFP(1);

				// Apply transformation to remaining columns.
				for (int j = k+1; j < n; j++) {
					int s = 0; 
					for (int i = k; i < m; i++) {
						s += FPMath.FPMul(QR[i][k],QR[i][j]);
					}
					s = FPMath.FPDiv(-s,QR[k][k]);
					for (int i = k; i < m; i++) {
						QR[i][j] += FPMath.FPMul(s,QR[i][k]);
					}
				}
			}
			Rdiag[k] = -nrm;
		}
	}

	/* ------------------------
   Public Methods
	 * ------------------------ */

	/** Is the Matrice full rank?
   @return     true if R, and hence A, has full rank.
	 */

	public boolean isFullRank () {
		for (int j = 0; j < n; j++) {
			if (Rdiag[j] == 0)
				return false;
		}
		return true;
	}

	/** Return the Householder vectors
   @return     Lower trapezoidal Matrice whose columns define the reflections
	 */

	public Matrice getH () {
		Matrice X = new Matrice(m,n);
		int[][] H = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (i >= j) {
					H[i][j] = QR[i][j];
				} else {
					H[i][j] = 0;
				}
			}
		}
		return X;
	}

	/** Return the upper triangular factor
   @return     R
	 */

	public Matrice getR () {
		Matrice X = new Matrice(n,n);
		int[][] R = X.getArray();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i < j) {
					R[i][j] = QR[i][j];
				} else if (i == j) {
					R[i][j] = Rdiag[i];
				} else {
					R[i][j] = 0;
				}
			}
		}
		return X;
	}

	/** Generate and return the (economy-sized) orthogonal factor
   @return     Q
	 */

	public Matrice getQ () {
		Matrice X = new Matrice(m,n);
		int[][] Q = X.getArray();
		for (int k = n-1; k >= 0; k--) {
			for (int i = 0; i < m; i++) {
				Q[i][k] = 0;
			}
			Q[k][k] = FPMath.IntToFP(1);
			for (int j = k; j < n; j++) {
				if (QR[k][k] != 0) {
					int s = 0;
					for (int i = k; i < m; i++) {
						s += FPMath.FPMul(QR[i][k],Q[i][j]);
					}
					s = FPMath.FPDiv(-s,QR[k][k]);
					for (int i = k; i < m; i++) {
						Q[i][j] += FPMath.FPMul(s,QR[i][k]);
					}
				}
			}
		}
		return X;
	}

	/** Least squares solution of A*X = B
   @param B    A Matrice with as many rows as A and any number of columns.
   @return     X that minimizes the two norm of Q*R*X-B.
   @exception  IllegalArgumentException  Matrice row dimensions must agree.
   @exception  RuntimeException  Matrice is rank deficient.
	 */

	public Matrice solve (Matrice B) {
		if (B.getRowDimension() != m) {
			throw new IllegalArgumentException("Matrice row dimensions must agree.");
		}
		if (!this.isFullRank()) {
			throw new RuntimeException("Matrice is rank deficient.");
		}

		// Copy right hand side
		int nx = B.getColumnDimension();
		int[][] X = B.getArrayCopy();

		// Compute Y = transpose(Q)*B
		for (int k = 0; k < n; k++) {
			for (int j = 0; j < nx; j++) {
				int s = 0; 
				for (int i = k; i < m; i++) {
					s += FPMath.FPMul(QR[i][k],X[i][j]);
				}
				s = FPMath.FPDiv(-s,QR[k][k]);
				for (int i = k; i < m; i++) {
					X[i][j] += FPMath.FPMul(s,QR[i][k]);
				}
			}
		}
		// Solve R*X = Y;
		for (int k = n-1; k >= 0; k--) {
			for (int j = 0; j < nx; j++) {
				X[k][j] = FPMath.FPMul(X[k][j],Rdiag[k]);
			}
			for (int i = 0; i < k; i++) {
				for (int j = 0; j < nx; j++) {
					X[i][j] -= FPMath.FPMul(X[k][j],QR[i][k]);
				}
			}
		}
		return (new Matrice(X,n,nx).getMatrice(0,n-1,0,nx-1));
	}
}

