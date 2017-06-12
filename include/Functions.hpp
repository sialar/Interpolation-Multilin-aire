#ifndef FUNCTIONS
#define FUNCTIONS

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>

#include "Utils.hpp"
#include "MultiVariatePoint.hpp"
using namespace std;

typedef vector<double>(*Function)(MultiVariatePoint<double>, int n);

class Functions
{
    private:
      static vector<vector<vector<double>>> m_coefs;
      static int m_polynomialDegree;

    public:
        Functions()  {};
        ~Functions() {};

        static void setCoefs(int degree, int d, int n);
        static vector<double> toAlternatingVector(double x, int n);
        static vector<double> functionToPlot(MultiVariatePoint<double> x, int n);
        static vector<double> function1D(MultiVariatePoint<double> x, int n);
        static vector<double> sinOfNorm2(MultiVariatePoint<double> x, int n);
        static vector<double> autoPolynomialFunction(MultiVariatePoint<double> x, int n);
        static void saveCoefsInFile(int d, int n);

        // Function to approximate
        static vector<double> f(MultiVariatePoint<double> x, int n);
};

#endif
