#ifndef UTILS
#define UTILS

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>

#include "MultiVariatePoint.hpp"
#include "BinaryTree.hpp"

using namespace std;

class Utils
{
    private:
        static vector<double> m_1dGrid;

    public:
        Utils()  {};
        ~Utils() {};

        // Display data and results
        static void separateur();
        static void displayPoints(vector<double> points);
        static void displayPoints(vector<MultiVariatePoint<double>> v);

        // Useful functions
        static double randomValue(double a, double b);
        static MultiVariatePoint<double> createRandomMultiVariatePoint(int d);
        static double interpolationError(vector<double> realValue, vector<double> estimate);

        // Intermediate function for Leja points computation
        static bool isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold);
        static double computeNewLejaPointFromSequence(vector<double> seq);

        // Write data in external file
        static void storeResult(vector<MultiVariatePoint<double>> x, vector<double> approx, vector<double> realValue);
        static void storeFunction(vector<double> x, vector<double> y, vector<double> z);
        static void storeLejaSequenceInFile(vector<vector<double>> x);
        static void storeLejaSequenceInFile(int length);
        static vector<double> loadLejaSequenceFromFile(int length);
        static vector<double> createSequenceByDichotomy(int length);

        // Function to approximate
        static double gNd(MultiVariatePoint<double> x);

        // Useful for uniform points and Leja sequence creation
        static void binaryDecomposition(int number, vector<double>& binary_decomp);

        // Create uniform sequence of data points
        static vector<double> createUniformSequence(int nbPoints);
        // Create sequence of Tchebychev zeros to construct Leja sequence
        static vector<double> createChebychevSequence(int nbPoints);
        // Create sequence of Leja points
        static vector<double> createLejaSequence(int nbPoints);
        // Create sequence of middle points
        static vector<double> createSequenceOfMiddles(int nbPoints);

};

#endif