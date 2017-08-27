#include "../include/AnalyticalFunctions.hpp"

AnalyticalFunctions::AnalyticalFunctions(int d, int n)
{
	m_n = n;
	m_d = d;
}

vector<double> AnalyticalFunctions::evaluate(MultiVariatePoint<double> x)
{
	double res = 0;
	for (int i=0; i<m_d; i++)
		res += 	pow(x(i), 2);
	vector<double> vecres(m_n,sqrt(res));
	return vecres;
}
