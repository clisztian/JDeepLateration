package ch.imetrica.jdeeplateration.tdoafoa;

import java.util.ArrayList;

import ch.imetrica.jdeeplateration.matrix.Matrix;



public class NavigationAlgorithms {

	
	int m_nbDataSets;
	boolean m_localOriginSet = false;
	
	Matrix Hn; 
	Matrix W = null;
	Matrix rho;
	
	ArrayList<Double> m_tdoa;
	ArrayList<Double> m_foa;
	ArrayList<Double> m_t;
	
	ArrayList<double[]> m_v;
	ArrayList<double[]> m_x;
	
	double[] m_uHat;
	double[] m_localOrigin;
	private int m_startPoint; 
	
	final static double m_cspeed = 300000000.0;
	final static double m_fc = 2147000000.0;
	
    public NavigationAlgorithms() {
    	
    	m_localOriginSet = false;
    	m_localOrigin = new double[3];   	
    }
	

	public int AddPoint(double[] x, double[] v, double tdoa, double foa, double time)
	{
		m_nbDataSets++;
		m_tdoa.add(tdoa);
		m_foa.add(foa);
		m_t.add(time);
		m_v.add(v);

		if (m_localOriginSet)
			CoordinateTranslation(x);
		else
			return -1;

		m_x.add(x);
		return 0;
	}

	public int CoordinateTranslation(double[] point)
	{ 
		for (int i = 0; i < 3; i++) {
			point[i] = point[i] - m_localOrigin[i];
		}

		return 0;
	}

	public int SetTOAWeights()
	{
		for (int i = 0; i < m_nbDataSets - 1; i++) {
			W.set(i, i, 1);
		}
		return 0;
	}

	public int SetLocalOrigin(double[] localOrigin)
	{	
		
		for (int i = 0; i < 3; i++) {
			m_localOrigin[i] = localOrigin[i];
		}
		
		m_localOriginSet = true;
		return 0;
	}

	public Matrix ComputeSphericalInterpolation()
	{
		int i;
		double tmp;
		double r = 0;
		
		double den;
		Matrix x_est_1, x_est_2;
		Matrix solution;
		
		Matrix R = new Matrix(m_nbDataSets);
		Matrix d = new Matrix(m_nbDataSets);
		Matrix Delta = new Matrix(m_nbDataSets);
		Matrix V = Matrix.ident(m_nbDataSets-1);
		Matrix W = Matrix.ident(m_nbDataSets-1);
		Matrix S = new Matrix(m_nbDataSets-1,3);



			//for each navigation location other than the first point (origin)
		for (i = 0; i < m_nbDataSets - 1; i++)
		{
			//translate each vector to local origin and compute the length
			r = 0;
			double[] x = m_x.get(i + 1);
			for (int j = 0; j < 3; j++) {
								
				tmp = x[j]; 
				S.set(i, j,tmp);
				r += tmp*tmp;
			}
			R.set(i, 0, Math.sqrt(r));
				
			d.set(i, 0, (m_tdoa.get(i + 1) - m_tdoa.get(0))*m_cspeed); 
			Delta.set(i, 0, R.getW(i, 0)*R.getW(i, 0) - d.getW(i, 0)*d.getW(i, 0));
		}

		Matrix Strans = S.trans();
		Matrix StransW = Strans.mul(W);
		
		
		
		H = inv(S.t()*W*S)*S.t()*W;
		
		
		Ps.set_size(m_nbDataSets - 1, m_nbDataSets - 1);
		Ps = S*H;
		Ps_v = eye(m_nbDataSets - 1, m_nbDataSets - 1) - Ps;
			

			mat num(m_nbDataSets - 1, m_nbDataSets - 1);


			num = (d*(d.t()*(Ps_v.t()*(V*Ps_v))));
			den = as_scalar(d.t()*Ps_v.t()*V*Ps_v*d);
			

			I.eye(m_nbDataSets - 1, m_nbDataSets - 1);

			solution = 0.5*(H*(I - num / den)*DELTA);

			num = d.t()*Ps_v.t()*V*Ps_v;
			den = 2 * as_scalar(d.t()*Ps_v.t()*V*Ps_v*d);
			mat Rs = 1 / den*num*DELTA;
			x_est_2 = 0.5*H*(DELTA - 2 * as_scalar(Rs)*d);
		
	}
	
	
	
	

	public void SetLinearisationPoint(double[] u)
	{
		m_uHat = new double[3];
		for (int i = 0; i < 3; i++) {
			m_uHat[i] = u[i];
		}
	}

	double dfdx(double[] xj, double[] vj)
	{
		double den, num;
		den = Math.sqrt(Math.pow((xj[0] - m_uHat[0]), 2) + Math.pow((xj[1] - m_uHat[1]), 2) + Math.pow((xj[2] - m_uHat[2]), 2));
		num = (xj[0] - m_uHat[0])*(vj[0] * (xj[0] - m_uHat[0]) + vj[1] * (xj[1] - m_uHat[1]) + vj[2] * (xj[2] - m_uHat[2]));
		return -m_fc / m_cspeed*(num / Math.pow(den, 1.5) - vj[0] / Math.pow(den, 0.5));
	}

	double dfdy(double[] xj, double[] vj)
	{
		double den, num;
		den = Math.sqrt(Math.pow((xj[0] - m_uHat[0]), 2) + Math.pow((xj[1] - m_uHat[1]), 2) + Math.pow((xj[2] - m_uHat[2]), 2));
		num = (xj[1] - m_uHat[1])*(vj[0] * (xj[0] - m_uHat[0]) + vj[1] * (xj[1] - m_uHat[1]) + vj[2] * (xj[2] - m_uHat[2]));
		return -m_fc / m_cspeed*(num / Math.pow(den, 1.5) - vj[1] / Math.pow(den, 0.5));
	}

	double dfdz(double[] xj, double[] vj)
	{
		double den, num;
		den = Math.sqrt(Math.pow((xj[0] - m_uHat[0]), 2) + Math.pow((xj[1] - m_uHat[1]), 2) + Math.pow((xj[2] - m_uHat[2]), 2));
		num = (xj[2] - m_uHat[2])*(vj[0] * (xj[0] - m_uHat[0]) + vj[1] * (xj[1] - m_uHat[1]) + vj[2] * (xj[2] - m_uHat[2]));
		return -m_fc / m_cspeed*(num / Math.pow(den, 1.5) - vj[2] / Math.pow(den, 0.5));
	}

	double dpdx(double[] xj)
	{
		double den, num;
		num = -(xj[0] - m_uHat[0]);
		den = Math.sqrt(Math.pow((xj[0] - m_uHat[0]), 2) + Math.pow((xj[1] - m_uHat[1]), 2) + Math.pow((xj[2] - m_uHat[2]), 2));
		return num / den;
	}

	double dpdy(double[] xj)
	{
		double den, num;
		num = -(xj[1] - m_uHat[1]);
		den = Math.sqrt(Math.pow((xj[0] - m_uHat[0]), 2) + Math.pow((xj[1] - m_uHat[1]), 2) + Math.pow((xj[2] - m_uHat[2]), 2));
		return num / den;
	}

	double dpdz(double[] xj)
	{
		double den, num;
		num = -(xj[2] - m_uHat[2]);
		den = Math.sqrt(Math.pow((xj[0] - m_uHat[0]), 2) + Math.pow((xj[1] - m_uHat[1]), 2) + Math.pow((xj[2] - m_uHat[2]), 2));
		return num / den;
	}

	double p(double[] xj)
	{
		return Math.sqrt(Math.pow((xj[0] - m_uHat[0]), 2) + Math.pow((xj[1] - m_uHat[1]), 2) + Math.pow((xj[2] - m_uHat[2]), 2));
	}
	
	double f(double[] xj, double[] vj)
	{
		double den, num;
		den = Math.sqrt(Math.pow((xj[0] - m_uHat[0]), 2) + Math.pow((xj[1] - m_uHat[1]), 2) + Math.pow((xj[2] - m_uHat[2]), 2))*m_cspeed;
		num = (vj[0] * (xj[0] - m_uHat[0]) + vj[1] * (xj[1] - m_uHat[1]) + vj[2] * (xj[2] - m_uHat[2]))*m_fc;
		return m_fc - num / den;
	}

	public void CreateHn_TDOAFDOA()
	{
		double dp1dx = dpdx(m_x.get(m_startPoint));
		double dp1dy = dpdy(m_x.get(m_startPoint));
		double dp1dz = dpdz(m_x.get(m_startPoint));

		
		int j = 0;
		Hn = new Matrix(2*(m_nbDataSets - m_startPoint - 1), 4);
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++) {
			
			double[] x = m_x.get(i + 1);
			double[] v = m_v.get(i + 1);
			
			Hn.set(2*j, 0, dp1dx - dpdx(x));
			Hn.set(2*j, 1, dp1dy - dpdy(x));
			Hn.set(2*j, 2, dp1dz - dpdz(x));
			Hn.set(2*j, 3, m_t.get(m_startPoint) - m_t.get(i + 1)*m_cspeed);
			Hn.set(2 * j + 1, 0, dfdx(x, v));
			Hn.set(2 * j + 1, 1, dfdy(x, v));
			Hn.set(2 * j + 1, 2 , dfdz(x, v));
			Hn.set(2 * j + 1, 3, 0);
			j++;
		}
	}

	public void CreateHn_TDOA()
	{
		
		double dp1dx = dpdx(m_x.get(m_startPoint));
		double dp1dy = dpdy(m_x.get(m_startPoint));
		double dp1dz = dpdz(m_x.get(m_startPoint));
		
		int j = 0;
		Hn = new Matrix(m_nbDataSets - m_startPoint - 1, 4);
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{
			double[] x = m_x.get(i + 1);
			Hn.set(j, 0, dp1dx - dpdx(x));
			Hn.set(j, 1, dp1dy - dpdy(x));
			Hn.set(j, 2, dp1dz - dpdz(x));
			Hn.set(j, 3, (m_t.get(m_startPoint) - m_t.get(i + 1))*m_cspeed);
			j++;
		}
	}

	public void CreateHn_FDOA()
	{	
		int j = 0;
		double[] x;
		double[] v;
		double a, b, c;
		
		Hn = new Matrix(m_nbDataSets - m_startPoint - 1, 3);
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{
			x = m_x.get(i + 1);
			v = m_v.get(i + 1);
			
			a = dfdx(x, v);
			b = dfdy(x, v);
			c = dfdz(x, v);
			Hn.set(j, 0, a);
			Hn.set(j, 1, b);
			Hn.set(j, 2, c);
			j++;
		}
	}

	public void CreateRho_TDOAFDOA()
	{
		rho = new Matrix(2*(m_nbDataSets - m_startPoint - 1), 1);
		double p1 = p(m_x.get(m_startPoint));
		
		int j = 0;
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{

			rho.set(2*j, 0, (m_tdoa.get(m_startPoint) - m_tdoa.get(i + 1))*m_cspeed - (p1 - p(m_x.get(i + 1))));
			rho.set(2 * j + 1, 0, m_foa.get(i + 1) - f(m_x.get(i + 1), m_v.get(i + 1)));
			j++;
		}
	}

	public void CreateRho_TDOA()
	{
		int j = 0;
		rho = new Matrix(m_nbDataSets - m_startPoint- 1, 1);
		double p1 = p(m_x.get(m_startPoint));
		
		
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++) {

			rho.set(j, 0, (m_tdoa.get(m_startPoint) - m_tdoa.get(i + 1))*m_cspeed - (p1 - p(m_x.get(i + 1))));
			j++;
		}
	}

	public void CreateRho_FDOA()
	{
		int j = 0;
		rho = new Matrix(m_nbDataSets - m_startPoint - 1, 1);
		
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{
			rho.set(j, 0, m_foa.get(i + 1) - f(m_x.get(i + 1), m_v.get(i + 1)));
			j++;
		}
	}

	public void CreateHnMinus()
	{
		mat A(Hn.t()*Hn);
		mat B(inv(A));
		HnMinus = B*Hn.t();
		return 0;
	}

	int CNavigationAlgorithms::ComputeTDOATaylorApprox(int lastNbPoints)
	{
		vector<double> solution(4,0);
		m_navSolution.set_size(4,1);
		m_startPoint=0;

		if (m_nbDataSets < 3)
			return -1;

		if (m_nbDataSets>lastNbPoints)
			m_startPoint = m_nbDataSets - lastNbPoints-1;
			
		CreateHn_TDOA();
		CreateRho_TDOA();
		CreateHnMinus();

		m_navSolution = HnMinus*rho;

		//add new solution to list of all computed solutions
		for (int i = 0; i < 4; i++)
		{
			solution[i]= m_navSolution(i, 0);
		}

		m_allNavSolutions.push_back(solution); 
		//************************************

		return 0;
	}

	int CNavigationAlgorithms::ComputeFDOATaylorApprox(int lastNbPoints)
	{
		vector<double> solution(3, 0);
		m_navSolution.set_size(3, 1);
		m_startPoint = 0;

		if (m_nbDataSets < 3)
			return -1;

		if (m_nbDataSets>lastNbPoints)
			m_startPoint = m_nbDataSets - lastNbPoints - 1;

		CreateHn_FDOA();
		CreateRho_FDOA();
		CreateHnMinus();

		m_navSolution = HnMinus*rho;

		//add new solution to list of all computed solutions
		for (int i = 0; i < 3; i++)
		{
			solution[i] = m_navSolution(i, 0);
		}

		m_allNavSolutions.push_back(solution);
		//************************************
		return 0;
	}

	int CNavigationAlgorithms::ComputeTDOAFDOATaylorApprox(int lastNbPoints)
	{
		vector<double> solution(4, 0);
		m_navSolution.set_size(4, 1);
		m_startPoint = 0;

		if (m_nbDataSets < 3)
			return -1;

		if (m_nbDataSets>lastNbPoints)
			m_startPoint = m_nbDataSets - lastNbPoints - 1;

		CreateHn_TDOAFDOA();
		CreateRho_TDOAFDOA();
		CreateHnMinus();

		m_navSolution = HnMinus*rho;

		cout << m_navSolution << "\n";
		//add new solution to list of all computed solutions
		for (int i = 0; i < 4; i++)
		{
			solution[i] = m_navSolution(i, 0);
		}

		m_allNavSolutions.push_back(solution);
		//************************************
		return 0;
	}

	int CNavigationAlgorithms::ReadDataFromFile(string filename, int& nbPoints, CNavigationDataRecord& d)
	{
		ifstream myFile(filename, ifstream::in);
		string line;

		nbPoints = 0;

		////stackoverflow.com/questions/1894886/parsing-a-comma-delimited-stdstring
		//myFile.open(filename);//read from a text file with I and Q in two columns, seperated by a white space
		int i = 0;
		double value;
		stringstream ss;
		vector<double> tmp;

		std::cout << "\n" << "Start reading file" << '\n';
		if (myFile.is_open())
		{
			while (getline(myFile, line))
			{
				ss << line;

				nbPoints++;

				tmp.clear();
				//parse file

				//coordinates x,y,z
				for (int i = 0; i < 3; i++)
				{
					ss >> value;
					tmp.push_back(value);
				}
				d.x.push_back(tmp);
				tmp.clear();
				//speeds vx,vy,vz
				for (int i = 0; i < 3; i++)
				{
					ss >> value;
					tmp.push_back(value);
				}
				d.v.push_back(tmp);
				tmp.clear();

				//delay to set the channel simulator
				ss >> value;
				d.tdoa.push_back(value);
				//doppler frequency to set the channel simulator
				ss.clear();
				ss >> value;
				d.foa.push_back(value);
				//
				ss.clear();
				ss >> value;
				d.time.push_back(value);

				//clear sstring
				ss.clear();
				ss.str("");
			}
		}
		else
			return -1;

		myFile.close();

		std::cout << "Stop reading file" << '\n';
		return 0;
	}
	
	
	
	
}
