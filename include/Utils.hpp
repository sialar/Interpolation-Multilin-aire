#ifndef UTILS
#define UTILS

#include <iostream>
#include <vector>
#include <array>
#include <cmath>

using namespace std;

typedef array<double,2> point2D;
class Utils
{
    public:
        Utils()  {};
        ~Utils() {};

        static vector<double> createPointsForSimpleTest(int nbPoints);
        static vector<double> createDataPoints(int nbPoints);
        static void binaryDecomposition(int number, vector<double>& binary_decomp);
        static double randomValue(double a, double b);
        static double squareError(vector<double> realValue, vector<double> estimate);
};

#endif
