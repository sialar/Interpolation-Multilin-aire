#include "../include/Functions.hpp"

vector<vector<vector<double>>> Functions::m_coefs;
int Functions::m_polynomialDegree;
int Functions::m_nbPerturbations;
vector<double> Functions::m_perturbations;
vector<double> Functions::m_perturbationsWidth;

double Functions::getPointInPerturbationNeighborhood()
{
    if (!m_perturbations.size()) return 0;
    return m_perturbations[0] + m_perturbationsWidth[0]/10;
}

double Functions::hat(double a, double b, double fa, double fb, double t)
{
    double c = (a+b)/2;
    if (t <= c && t >= a) return ((t-a)/(c-a)+fa);
    else if (t >= c && t <= b) return ((t-b)/(c-b)+fb);
    else return 0;
}

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
    double temp = 0;
    for (int i=0; i<x.getD(); i++) temp += x(i);
    if (x.getD()>1) return toAlternatingVector(sin(temp),n);
    else return toAlternatingVector(2*exp(-x(0)*x(0) + sin(2*M_PI*x(0))),n);
}

int Functions::inOneHat(double x)
{
    double xc, xe;
    for (int i=0; i<m_nbPerturbations; i++)
    {
        xc = m_perturbations[i];
        xe = m_perturbationsWidth[i];
        if (x>=xc-xe/2 && x<=xc+xe/2)
            return i;
    }
    return -1;
}
vector<double> Functions::phi(MultiVariatePoint<double> x, int n)
{
    double temp = 1, f = 1;
    for (int i=1; i<x.getD(); i++) temp *= exp(-x(i));
    int hatPos = inOneHat(x(0));
    double xc, xe, a, b;
    if (hatPos>-1)
    {
        xc = m_perturbations[hatPos];
        xe = m_perturbationsWidth[hatPos];
        a = xc - xe/2;
        b = xc + xe/2;
        return toAlternatingVector(hat(a,b,sin(2*M_PI*f*a),sin(2*M_PI*f*b),x(0)) * temp,n);
    }
    else
        return toAlternatingVector(sin(2*M_PI*f*x(0)) * temp,n);
}

vector<double> Functions::cosinus(MultiVariatePoint<double> x, int n)
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
    m_nbPerturbations = Utils::randomValue(1,3);
    for (int i=0; i<m_nbPerturbations; i++)
    {
        m_perturbations.push_back(Utils::randomValue(-0.9,0.9));
        m_perturbationsWidth.push_back(Utils::randomValue(0.01,0.05));
    }
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

void Functions::savePerturbationsInFile()
{
    ofstream file(Utils::projectPath + "data/perturbations.txt", ios::out | ios::trunc);
    if(file)
    {
        file << m_nbPerturbations << endl;

        for (int i=0; i<m_nbPerturbations; ++i)
            file << m_perturbations[i] << " ";
        file << endl;

        for (int i=0; i<m_nbPerturbations; ++i)
            file << m_perturbationsWidth[i] << " ";
        file << endl;

        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}
