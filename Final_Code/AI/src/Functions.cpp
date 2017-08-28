#include "../include/Functions.hpp"

Functions::Functions(int d, int n)
{
	m_d = d;
	m_n = n;
	vector<double> defaultDomain(2,1);
	defaultDomain[0] = -1;
	for (int i=0; i<m_d; i++)
		m_parametersDomain.push_back(defaultDomain);
}

void Functions::displayParametersDomain()
{
    cout << " - Parameters domain of function : ";
    for (int i=0; i<m_d-1; i++)
      cout << "[" << m_parametersDomain[i][0] << "," << m_parametersDomain[i][1] << "] x ";
    cout << "[" << m_parametersDomain[m_d-1][0] << "," << m_parametersDomain[m_d-1][1] << "]" << endl;
}
