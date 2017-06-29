#include "../../include/Tucker/LagrangePolynomial.hpp"


LagrangePolynomial::LagrangePolynomial(vector<double> x, vector<double> fx)
{
    if (x.size() != fx.size())
    {
        cout << x.size() << " " << fx.size() << endl;
        cout << "Error : x and fx must be the same size" << endl;
        exit(0);
    }
    m_degree = x.size();
    m_axisValues = x;
    m_polynomialValues = fx;
}

double LagrangePolynomial::getInterpolation(double x)
{
    double f = 0.0;
    int m = 0;
    while (m < int(m_axisValues.size()))
    {
        double l = 1;
        int n = 0;
        while (n < int(m_axisValues.size()))
        {
            if (n != m)
                l *= (x - m_axisValues[n])/(m_axisValues[m] - m_axisValues[n]);
            n++;
        }
        f += m_polynomialValues[m]*l;
        m ++;
    }
    return f;
}

double LagrangePolynomial::getInterpolation_bis(double x)
{
    int n = m_axisValues.size();
    double valueOfPolynome = 0.0, li;
    for (int i=0; i<n; i++)
    {
        li = 1.0;
        for (int j=0; j<n; j++)
            if (i!=j)
                li = li*(x - m_axisValues[j])/(m_axisValues[i]-m_axisValues[j]);
        valueOfPolynome += m_polynomialValues[i]*li;
    }
    return valueOfPolynome;
}

vector<double> LagrangePolynomial::getInterpolationArr_bis(vector<double> x)
{
    vector<double> farr;
    for (double xi : x)
        farr.push_back(getInterpolation_bis(xi));
    return farr;
}

vector<double> LagrangePolynomial::getInterpolationArr(vector<double> x)
{
    vector<double> farr;
    for (double xi : x)
        farr.push_back(getInterpolation(xi));
    return farr;
}

double LagrangePolynomial::getDerivative(double x)
{
    int n = m_axisValues.size();
    double polynome_derivative = 0.0;
    for (int i=0; i<n; i++)
    {
        double denominator_i = 1.0;
        for (int j=0; j<n; j++)
            if (i!=j)
                denominator_i *=  m_axisValues[i]-m_axisValues[j];

        double numerator_i = 0.0;
        for (int j=0; j<n; j++)
        {
            if (i!=j)
            {
                double numerator_i_derivative_at_xj = 1.0;
                for (int k=0; k<n; k++)
                    if (k!=i && k!=j)
                        numerator_i_derivative_at_xj *= x - m_axisValues[k];
                numerator_i += numerator_i_derivative_at_xj;
            }
        }
        polynome_derivative += m_polynomialValues[i]*numerator_i/denominator_i;
    }
    return  polynome_derivative;
}
