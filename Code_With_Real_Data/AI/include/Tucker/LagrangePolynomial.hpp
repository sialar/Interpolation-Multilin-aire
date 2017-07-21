#ifndef LAGRANGEPOLYNOMIAL
#define LAGRANGEPOLYNOMIAL

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <memory>
#include <limits>

#include "../Utils.hpp"
using namespace std;

class LagrangePolynomial
{
    private:
        int m_degree;
        vector<double> m_axisValues;
        vector<double> m_polynomialValues;

    public:
        ~LagrangePolynomial() {};
        LagrangePolynomial() {};
        LagrangePolynomial(vector<double> x, vector<double> fx);

        const int degree() { return m_degree; };
        const vector<double> axisValues() { return m_axisValues; };
        const vector<double> polynomialValues() { return m_polynomialValues; };

        double getInterpolation(double x);
        double getDerivative(double x);
        double getInterpolation_bis(double x);
        vector<double> getInterpolationArr_bis(vector<double> x);
        vector<double> getInterpolationArr(vector<double> x);
};

#endif
