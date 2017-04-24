#ifndef UTILS
#define UTILS

#include <fstream>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <limits>
using namespace std;

class Utils
{
    private:
        static vector<double> m_1dGrid;

    public:
        Utils()  {};
        ~Utils() {};

        static void separateur();
        static void displayPoints(vector<vector<double>> v, int d);
        static vector<double> displayGRealValues(vector<double> vx, vector<double> vy, int d, bool debug);
        static void displayApproximation(vector<double> approx, int nx, int ny, int d, bool debug);

        static double g1d(double x);
        static double g2d(double x, double y);
        static double gNd(vector<double> x);
        static double randomValue(double a, double b);
        static double squareError(vector<double> realValue, vector<double> estimate);
        static void binaryDecomposition(int number, vector<double>& binary_decomp);
        static bool isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold);
        static double computeNewLejaPointFromSequence(vector<double> seq);
        static void storeResult1D(vector<double> x, vector<double> y, vector<double> real_y);
        static void storeResult2D(vector<double> x, vector<double> y, vector<double> z, vector<double> real_z);
        static void store2DLejaSequenceInFile(vector<double> x, vector<double> y);
        static void store3DLejaSequenceInFile(vector<double> x, vector<double> y, vector<double> z);

        // Uniform sequence of data points
        static vector<double> createUniformSequence(int nbPoints);
        // Sequence of Leja points
        static vector<double> createLejaSequence(int nbPoints);
        // Sequence of Tchebychev zeros
        static vector<double> createChebychevSequence(int nbPoints);

};

#endif
