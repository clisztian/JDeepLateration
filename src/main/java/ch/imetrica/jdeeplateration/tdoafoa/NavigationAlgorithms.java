package ch.imetrica.jdeeplateration.tdoafoa;

public class NavigationAlgorithms {

	
	int m_nbDataSets;
	int m_localOriginSet;
	
	ArrayList<Double> m_tdoa;
	ArrayList<Double> m_foa;
	
	double[] m_localOrigin; 
	
	final static double m_cspeed = 300000000.0;
	final static double m_fc = 2147000000.0;
	
    public NavigationAlgorithms() {
    	
    	m_localOriginSet = 0;
    	m_localOrigin = {0, 0, 0};
    	
    	
    }
	


	int CNavigationAlgorithms::AddPoint(vector<double>& x, vector<double>& v, double tdoa, double foa, double time)
	{
		m_nbDataSets++;
		m_tdoa.push_back(tdoa);
		m_foa.push_back(foa);
		m_t.push_back(time);
		m_v.push_back(v);

		if (m_localOriginSet)
			CoordinateTranslation(x);
		else
			return -1;

		m_x.push_back(x);
		return 0;
	}

	int CNavigationAlgorithms::CoordinateTranslation(vector<double>& point)
	{ 
		for (int i = 0; i < 3; i++)
			point[i] = point[i] - m_localOrigin[i];

		return 0;
	}

	int CNavigationAlgorithms::SetTOAWeigths(void)
	{
		for (int i = 0; i < m_nbDataSets - 1; i++)
		{
			W(i, i) = 1;// m_weights[i];
		}
		return 0;
	}

	int CNavigationAlgorithms::SetLocalOrigin(vector<double>& localOrigin)
	{	
		//we must copy x[0] into localOrigin, if we only assign addresses, content will be changed as soon as we act on x[0] in "AddPoint"-Method
		for (int i = 0; i < 3; i++)
			m_localOrigin[i] = localOrigin[i];

		m_localOriginSet = 1;
		return 0;
	}

	int CNavigationAlgorithms::ComputeSphericalInterpolation(mat& solution)
	{
		double tmp;
		double r = 0;
		int i;
		double den;
		mat x_est_1, x_est_2;

		//set the first acquired point as the local origin
		if (m_localOriginSet)
		{

			//create vectors R, d and identity matrices M, V and W; set TOA weights to diagonal elements
			R.set_size(m_nbDataSets - 1, 1);
			d.set_size(m_nbDataSets - 1, 1);
			DELTA.set_size(m_nbDataSets - 1, 1);
			V.eye(m_nbDataSets - 1, m_nbDataSets - 1);//generate identity matrix
			W.eye(m_nbDataSets - 1, m_nbDataSets - 1);
			S.set_size(m_nbDataSets - 1, 3);
			SetTOAWeigths();

			//for each navigation location other than the first point (origin)
			for (i = 0; i < m_nbDataSets - 1; i++)
			{
				//translate each vector to local origin and compute the length
				r = 0;
				for (int j = 0; j < 3; j++)
				{
					tmp = m_x[i + 1][j];// -m_localOrigin[j];
					m_x[i + 1][j] = tmp;
					S(i, j) = tmp;
					r += pow(tmp, 2);
				}
				R(i, 0) = sqrt(r);
				//TDOA between current navigation point and origin (first point); tdoa[0] is inherently 0
				d(i, 0) = (m_tdoa[i + 1] - m_tdoa[0])*m_cspeed; //and that is why the equation here is the same as m_toa[i+1]-m_toa[0]
				DELTA(i, 0) = R(i, 0)*R(i, 0) - d(i, 0)*d(i, 0);
			}

			H = inv(S.t()*W*S)*S.t()*W;
			Ps.set_size(m_nbDataSets - 1, m_nbDataSets - 1);
			Ps = S*H;
			Ps_v = eye(m_nbDataSets - 1, m_nbDataSets - 1) - Ps;
			

			mat num(m_nbDataSets - 1, m_nbDataSets - 1);
			//d.print("d:");
			//Ps_v.print("Ps_v:");

			//mat a; mat b; mat c; mat e; mat f;
			//a = Ps_v;
			//a.print("a:");
			//b = V*a;
			//b.print("b:");
			//c = Ps_v.t()*b;
			//c.print("c:");
			//e = d.t()*c;
			//e.print("e:");
			//f = d*e;
			//f.print("f: ");
			//V.print("V: ");

			num = (d*(d.t()*(Ps_v.t()*(V*Ps_v))));
			den = as_scalar(d.t()*Ps_v.t()*V*Ps_v*d);
			/*num.print("num :");*/

			I.eye(m_nbDataSets - 1, m_nbDataSets - 1);

			solution = 0.5*(H*(I - num / den)*DELTA);//x_est_1 = 0.5*(H*(I - num / den)*DELTA);

			num = d.t()*Ps_v.t()*V*Ps_v;
			den = 2 * as_scalar(d.t()*Ps_v.t()*V*Ps_v*d);
			mat Rs = 1 / den*num*DELTA;
			x_est_2 = 0.5*H*(DELTA - 2 * as_scalar(Rs)*d);
			return 0;
		}
		else
			return -1;
	}

	int CNavigationAlgorithms::SetLinearisationPoint(vector<double> u)
	{
		m_uHat.set_size(3, 1);
		for (int i = 0; i < 3; i++)
			m_uHat(i) = u[i];
		return 0;
	}

	double CNavigationAlgorithms::dfdx(vector<double> xj, vector<double> vj)
	{
		double den, num;
		den = pow((xj[0] - m_uHat(0)), 2) + pow((xj[1] - m_uHat(1)), 2) + pow((xj[2] - m_uHat(2)), 2);
		num = (xj[0] - m_uHat(0))*(vj[0] * (xj[0] - m_uHat(0)) + vj[1] * (xj[1] - m_uHat(1)) + vj[2] * (xj[2] - m_uHat(2)));
		return -m_fc / m_cspeed*(num / pow(den, 1.5) - vj[0] / pow(den, 0.5));
	}

	double CNavigationAlgorithms::dfdy(vector<double> xj, vector<double> vj)
	{
		double den, num;
		den = pow((xj[0] - m_uHat(0)), 2) + pow((xj[1] - m_uHat(1)), 2) + pow((xj[2] - m_uHat(2)), 2);
		num = (xj[1] - m_uHat(1))*(vj[0] * (xj[0] - m_uHat(0)) + vj[1] * (xj[1] - m_uHat(1)) + vj[2] * (xj[2] - m_uHat(2)));
		return -m_fc / m_cspeed*(num / pow(den, 1.5) - vj[1] / pow(den, 0.5));
	}

	double CNavigationAlgorithms::dfdz(vector<double> xj, vector<double> vj)
	{
		double den, num;
		den = pow((xj[0] - m_uHat(0)), 2) + pow((xj[1] - m_uHat(1)), 2) + pow((xj[2] - m_uHat(2)), 2);
		num = (xj[2] - m_uHat(2))*(vj[0] * (xj[0] - m_uHat(0)) + vj[1] * (xj[1] - m_uHat(1)) + vj[2] * (xj[2] - m_uHat(2)));
		return -m_fc / m_cspeed*(num / pow(den, 1.5) - vj[2] / pow(den, 0.5));
	}

	double CNavigationAlgorithms::dpdx(vector<double> xj)
	{
		double den, num;
		num = -(xj[0] - m_uHat(0));
		den = sqrt(pow((xj[0] - m_uHat(0)), 2) + pow((xj[1] - m_uHat(1)), 2) + pow((xj[2] - m_uHat(2)), 2));
		return num / den;
	}

	double CNavigationAlgorithms::dpdy(vector<double> xj)
	{
		double den, num;
		num = -(xj[1] - m_uHat(1));
		den = sqrt(pow((xj[0] - m_uHat(0)), 2) + pow((xj[1] - m_uHat(1)), 2) + pow((xj[2] - m_uHat(2)), 2));
		return num / den;
	}

	double CNavigationAlgorithms::dpdz(vector<double> xj)
	{
		double den, num;
		num = -(xj[2] - m_uHat(2));
		den = sqrt(pow((xj[0] - m_uHat(0)), 2) + pow((xj[1] - m_uHat(1)), 2) + pow((xj[2] - m_uHat(2)), 2));
		return num / den;
	}

	double CNavigationAlgorithms::p(vector<double> xj)
	{
		return sqrt(pow((xj[0] - m_uHat(0)), 2) + pow((xj[1] - m_uHat(1)), 2) + pow((xj[2] - m_uHat(2)), 2));
	}

	double CNavigationAlgorithms::f(vector<double> xj, vector<double> vj)
	{
		double den, num;
		den = sqrt(pow((xj[0] - m_uHat(0)), 2) + pow((xj[1] - m_uHat(1)), 2) + pow((xj[2] - m_uHat(2)), 2))*m_cspeed;
		num = (vj[0] * (xj[0] - m_uHat(0)) + vj[1] * (xj[1] - m_uHat(1)) + vj[2] * (xj[2] - m_uHat(2)))*m_fc;
		return m_fc - num / den;
	}

	int CNavigationAlgorithms::CreateHn_TDOAFDOA(void)
	{
		double dp1dx = dpdx(m_x[m_startPoint]);
		double dp1dy = dpdy(m_x[m_startPoint]);
		double dp1dz = dpdz(m_x[m_startPoint]);
		vector<double> x;
		vector<double> v;
		int j = 0;
		Hn.set_size(2*(m_nbDataSets - m_startPoint - 1), 4);
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{
			x = m_x[i + 1];
			v = m_v[i + 1];
			Hn(2*j, 0) = dp1dx - dpdx(x);
			Hn(2*j, 1) = dp1dy - dpdy(x);
			Hn(2*j, 2) = dp1dz - dpdz(x);
			Hn(2*j, 3) = (m_t[m_startPoint] - m_t[i + 1])*m_cspeed;
			Hn(2 * j + 1, 0) = dfdx(x, v);
			Hn(2 * j + 1, 1) = dfdy(x, v);
			Hn(2 * j + 1, 2) = dfdz(x, v);
			Hn(2 * j + 1, 3) = 0;
			j++;
		}

		return 0;
	}

	int CNavigationAlgorithms::CreateHn_TDOA(void)
	{
		double dp1dx = dpdx(m_x[m_startPoint]);
		double dp1dy = dpdy(m_x[m_startPoint]);
		double dp1dz = dpdz(m_x[m_startPoint]);
		vector<double> x;
		int j = 0;
		Hn.set_size(m_nbDataSets - m_startPoint - 1, 4);
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{
			x = m_x[i + 1];
			Hn(j, 0) = dp1dx - dpdx(x);
			Hn(j, 1) = dp1dy - dpdy(x);
			Hn(j, 2) = dp1dz - dpdz(x);
			Hn(j, 3) = (m_t[m_startPoint] - m_t[i + 1])*m_cspeed;
			j++;
		}

		return 0;
	}

	int CNavigationAlgorithms::CreateHn_FDOA(void)
	{	
		int j = 0;
		vector<double> x;
		vector<double> v;
		double a, b, c;
		Hn.set_size(m_nbDataSets - m_startPoint - 1, 3);
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{
			x = m_x[i + 1];
			v = m_v[i + 1];
			a= dfdx(x, v);
			b = dfdy(x, v);
			c= dfdz(x, v);
			Hn(j, 0) = a;
			Hn(j, 1) = b;
			Hn(j, 2) = c;
			j++;
		}
		return 0;
	}

	int CNavigationAlgorithms::CreateRho_TDOAFDOA(void)
	{
		rho.set_size(2*(m_nbDataSets - m_startPoint - 1), 1);
		double p1 = p(m_x[m_startPoint]);
		int j = 0;
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{
			//m_tdoa[i] refers to m_toa[0], which is 0 (see 'InitialFix'),
			//the equation below is the same as m_toa[0]-m_toa[i+1]
			rho(2*j, 0) = (m_tdoa[m_startPoint] - m_tdoa[i + 1])*m_cspeed - (p1 - p(m_x[i + 1]));
			rho(2 * j + 1, 0) = m_foa[i + 1] - f(m_x[i + 1], m_v[i + 1]);
			j++;
		}

		return 0;
	}

	int CNavigationAlgorithms::CreateRho_TDOA(void)
	{
		rho.set_size(m_nbDataSets - m_startPoint- 1, 1);
		double p1 = p(m_x[m_startPoint]);
		int j = 0;
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{
			//m_tdoa[i] refers to m_toa[0], which is 0 (see 'InitialFix'),
			//the equation below is the same as m_toa[0]-m_toa[i+1]
			rho(j, 0) = (m_tdoa[m_startPoint] - m_tdoa[i + 1])*m_cspeed - (p1 - p(m_x[i + 1]));
			j++;
		}

		return 0;
	}

	int CNavigationAlgorithms::CreateRho_FDOA(void)
	{
		rho.set_size(m_nbDataSets - m_startPoint - 1, 1);
		int j = 0;
		for (int i = m_startPoint; i < m_nbDataSets - 1; i++)
		{
			rho(j, 0) = m_foa[i + 1] - f(m_x[i + 1], m_v[i + 1]);
			j++;
		}
		return 0;
	}

	int CNavigationAlgorithms::CreateHnMinus(void)
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
