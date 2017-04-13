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
        static double randomValue(double a, double b);
        static double squareError(vector<double> realValue, vector<double> estimate);
        static void binaryDecomposition(int number, vector<double>& binary_decomp);
        static bool isTooCloseToOneLejiPoint(double y, vector<double> seq, double threshold);
        static double computeNewLejaPointFromSequence(vector<double> seq);
        static void storeLejaSequenceInFile(vector<double> seq);

        // For uniform sequence of data points
        static vector<double> createUniformSequence(int nbPoints);
        // For Leja sequence (L)
        static vector<double> createLejaSequence(int nbPoints);
        // For Tchebychev zeros
        static vector<double> createChebychevSequence(int nbPoints);

};

#endif
