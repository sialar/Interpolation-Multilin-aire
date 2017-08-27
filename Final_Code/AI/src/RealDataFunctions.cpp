#include "../include/RealDataFunctions.hpp"

RealDataFunctions::RealDataFunctions(int d, int n, string filePath)
{
	m_n = n;
	m_d = d;
	m_filePath = filePath;
}

vector<double> RealDataFunctions::evaluate(MultiVariatePoint<double> x)
{
	double res = 0;
	for (int i=0; i<m_d; i++)
		res += 	pow(x(i), i);
	vector<double> vecres(m_n,res);
	return vecres;
}
