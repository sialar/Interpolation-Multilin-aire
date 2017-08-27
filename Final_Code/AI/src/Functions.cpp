#include "../include/Functions.hpp"

// Exemple de fonctions analytiques à intepoler
// Fonction polynomiale
vector<double> p(MultiVariatePoint<double> x, int d, int n)
{
	double res = 0;
	for (int i=0; i<d; i++)
		res += 	pow(x(i), i);
	vector<double> vecres(n,res);
	return vecres;
}

// sin( norm_2(x) )
vector<double> f(MultiVariatePoint<double> x, int d, int n)
{
	double res = 0;
	for (int i=0; i<d; i++)
		res += 	pow(x(i), 2);
	vector<double> vecres(n,sqrt(res));
	return vecres;
}

// sqrt(1-x**2) * sin( sum(xi)_{0<j<d} )
vector<double> g(MultiVariatePoint<double> x, int d, int n)
{
	double res = 0;
	for (int i=0; i<d; i++)
		res += 	pow(x(i), 2);
	vector<double> vecres(n,sqrt(res));
	return vecres;
}


Functions::Functions(int d, int n)
{
	m_d = d;
	m_n = n;
}
