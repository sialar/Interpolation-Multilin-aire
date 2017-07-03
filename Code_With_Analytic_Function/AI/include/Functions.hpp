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
      static int m_nbPerturbations;
      static vector<double> m_perturbations;
      static vector<double> m_perturbationsWidth;

    public:
        Functions() {};
        ~Functions() {};

        static void savePerturbationsInFile();
        static void saveCoefsInFile(int d, int n);
        static void setCoefs(int degree, int d, int n);
        static double changeFunctionDomain(double a, double b, double x);
        static double hat(double a, double b, double fa, double fb, double t);
        static double getPointInPerturbationNeighborhood();
        static int inOneHat(double x);

        static vector<double> toAlternatingVector(double x, int n);
        static vector<double> functionToPlot(MultiVariatePoint<double> x, int n);
        static vector<double> function1D(MultiVariatePoint<double> x, int n);
        static vector<double> sinOfNorm2(MultiVariatePoint<double> x, int n);
        static vector<double> autoPolynomialFunction(MultiVariatePoint<double> x, int n);
        static vector<double> cosinus(MultiVariatePoint<double> x, int n);
        static vector<double> phi(MultiVariatePoint<double> x, int n);
        static vector<double> h(MultiVariatePoint<double> x, int n);

        // Function to approximate
        static vector<double> f(MultiVariatePoint<double> x, int n);
};

#endif