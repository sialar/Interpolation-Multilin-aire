#include "../include/Functions.hpp"

vector<vector<vector<double>>> Functions::m_coefs;
int Functions::m_polynomialDegree;

double Functions::changeFunctionDomain(double a, double b, double x)
{
    return (2*x)/(b-a) + 1 - (2*b)/(b-a);
}

vector<double> Functions::f(MultiVariatePoint<double> x, int n)
{
    // write here the interpolated function
    return autoPolynomialFunction(x,n);
}

vector<double> Functions::toAlternatingVector(double x, int n)
{
    vector<double> res(n,x);
    for (int i=0; i<n; i++)
        if (i%2) res[i] *= -1;
    return res;
}

vector<double> Functions::functionToPlot(MultiVariatePoint<double> x, int n)
{
    double temp = 1;
    for (int i=1; i<x.getD(); i++) temp *= exp(-x(i));
    return toAlternatingVector(sqrt(1-x(0)*x(0)) * temp, n);
}

vector<double> Functions::h(MultiVariatePoint<double> x, int n)
{
    double temp = 1, f = 2;
    for (int i=1; i<x.getD(); i++) temp *= exp(-x(i));
    return toAlternatingVector(cos(2*M_PI*f*x(0)) * temp,n);
}

vector<double> Functions::function1D(MultiVariatePoint<double> x, int n)
{
    if (x.getD()==1)
    {
        double f = 10;
        if (x(0)<0) return toAlternatingVector(cos(2*M_PI*f*x(0)),n);
        else return toAlternatingVector(cos(0.1*M_PI*f*x(0)),n);
    }
    else return toAlternatingVector(0,n);
}

vector<double> Functions::sinOfNorm2(MultiVariatePoint<double> x, int n)
{
    double temp = 0;
    for (int i=0; i<x.getD(); i++)
        temp += pow(x(i),2);
    return toAlternatingVector(sin(sqrt(temp)),n);
}

vector<double> Functions::autoPolynomialFunction(MultiVariatePoint<double> x, int n)
{
    vector<double> res;
    res.resize(n,0);
    for (int i=0; i<n; i++)
        for (int j=0; j<m_polynomialDegree; j++)
            for (int k=0; k<x.getD(); k++)
                res[i] += m_coefs[i][j][k] * pow(x(k),j+1);
    return res;
}

void Functions::setCoefs(int degree, int d, int n)
{
    m_polynomialDegree = degree;
    m_coefs.resize(n);
    for (int i=0; i<n; i++)
        m_coefs[i].resize(degree);
    for (int i=0; i<n; i++)
        for (int j=0; j<degree; j++)
            m_coefs[i][j].resize(d);

    for (int i=0; i<n; i++)
        for (int j=0; j<degree; j++)
            for (int k=0; k<d; k++)
                m_coefs[i][j][k] = Utils::randomValue(-1,1);

    saveCoefsInFile(d,n);
}

void Functions::saveCoefsInFile(int d, int n)
{
    ofstream file(Utils::projectPath + "data/polynomial_coefs.txt", ios::out | ios::trunc);
    if(file)
    {
        file << n << " " << m_polynomialDegree << " " << d << endl;
        for (int i=0; i<n; ++i)
        {
            for (int j=0; j<m_polynomialDegree; ++j)
            {
                for (int k=0; k<d; k++)
                {
                    file << m_coefs[i][j][k] << " ";
                }
                file << endl;
            }
            file << endl;
        }
        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}
