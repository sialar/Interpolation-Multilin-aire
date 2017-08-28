#include "../include/RealDataFunctions.hpp"

RealDataFunctions::RealDataFunctions(int d, int n, string filePath) : Functions(d,n)
{
	m_n = n;
	m_d = d;
	m_filePath = filePath;
	readDataFromFile();
	m_tuckerApprox = make_shared<TuckerApproximation>(m_coreName,m_csName);
}

vector<double> RealDataFunctions::evaluate(MultiVariatePoint<double> x)
{
		vector<double> res;
		for (int i=0; i<m_d; i++)
				x(i) = Utils::convertToFunctionDomain(m_parametersDomain[i][0], m_parametersDomain[i][1], x(i));
		res.push_back(m_tuckerApprox->evaluate(x,m_csName));
    return res;
}

void RealDataFunctions::readDataFromFile()
{
    ifstream file(m_filePath, ios::in);
    if(file)
  	{
	    	string line;
	      vector<double> data;
	      vector<double> max(m_d,-numeric_limits<double>::max());
	      vector<double> min(m_d,numeric_limits<double>::max());
	      MultiVariatePoint<double> p(m_d,0,0);
	      while (getline(file, line))
	      {
	          data = Utils::str2vector(line);
						if (int(data.size()) != m_d + m_n + 1)
						{
								cerr << "Choix des dimensions (n et d) et contenu du fichier incompatibles" << endl;
								cout << data.size() << " " << m_d + m_n + 1 << endl;
								exit(1);
						}

        		for (int i=0; i<m_d; i++)
	              p(i) = data[i+1];
        		m_refPoints.push_back(p);
			      for (int i=0; i<m_d; i++)
			      {
			      		if (p(i) > max[i]) max[i] = p(i);
			          if (p(i) < min[i]) min[i] = p(i);
			      }
						for (int i=0; i<m_d+1; i++)
								data.erase(data.begin());
						m_exactValues.push_back(data);
	      }
	      for (int i=0; i<m_d; i++)
	      {
						m_parametersDomain[i][0] = min[i];
						m_parametersDomain[i][1] = max[i];
	      }
	      for (int i=0; i<int(m_refPoints.size()); i++)
	          for (int j=0; j<m_d; j++)
	              m_refPoints[i](j) = Utils::adaptCoordsToFunctionDomain(min[j], max[j], m_refPoints[i](j));
	      file.close();
	  }
	  else cerr << "Error while opening file!" << endl;
}
