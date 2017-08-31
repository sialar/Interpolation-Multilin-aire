#include "../include/AnalyticalFunctions.hpp"

/******************************************************************************/
/****************  Exemple de fonctions analytiques à intepoler ***************/

// Polynome P(x) = (p(x), ..., p(x)) : n fois
// p(x) = x(0) + x(1)^2 + x(2)^3 + ... + x(d-1)^d
vector<double> p(MultiVariatePoint<double> x, int d, int n)
{
	double res = 0;
	for (int i=0; i<d; i++)
		res += 	pow(x(i), i+1);
	vector<double> vecres(n,res);
	return vecres;
}

// (sin( norm_2(x) ), .., sin( norm_2(x) )) : n fois
vector<double> f(MultiVariatePoint<double> x, int d, int n)
{
	double res = 0;
	for (int i=0; i<d; i++)
		res += 	pow(x(i), 2);
	vector<double> vecres(n,sqrt(res));
	return vecres;
}

// (sqrt(1-x**2) * sin( sum(xi)_{0<j<d} ), .., sqrt(1-x**2) * sin( sum(xi)_{0<j<d} )) : n fois
vector<double> g(MultiVariatePoint<double> x, int d, int n)
{
	double res = 1;
	for (int i=1; i<d; i++)
		res += pow(x(i),i);
	vector<double> vecres(n,sqrt(1-x(0)*x(0)) * res);
	return vecres;
}

/******************************************************************************/
/******************************************************************************/

AnalyticalFunctions::AnalyticalFunctions(int d, int n) : Functions(d,n)
{
	m_n = n;
	m_d = d;
}

vector<double> AnalyticalFunctions::evaluate(MultiVariatePoint<double> x)
{
		// Ici on choisi la fonction analytique qu'on veut approcher
		return g(x,m_d,m_n);
}
