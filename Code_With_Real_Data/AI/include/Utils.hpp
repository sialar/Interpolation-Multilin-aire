#ifndef UTILS
#define UTILS

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>

#include "MultiVariatePoint.hpp"

using namespace std;

class Utils
{
    private:
        static vector<double> m_1dGrid;

    public:
        static string projectPath;

        Utils()  {};
        ~Utils() {};

        // Display data and results
        static void separateur();
        static bool displayResults();
        static void displayValues(vector<double> v);
        static void displayPoints(vector<vector<double>> v);
        static void displayPoints(vector<MultiVariatePoint<double>> v);

        // Useful functions
        static double relativeInterpolationError(vector<vector<double>> realValue, vector<vector<double>> estimate);
        static double mseInterpolationError(vector<vector<double>> realValue, vector<vector<double>> estimate);
        static MultiVariatePoint<double> createRandomMultiVariatePoint(int d);
        static vector<double> diff(vector<double> x,vector<double> y);
        static double randomValue(double a, double b);
        static string vector2str(vector<double> x);
        static double norm(vector<double> x, int p);
        static double maxAbsValue(vector<double> vec);
        static double adaptCoordsToFunctionDomain(double a, double b, double x);
        static double convertToFunctionDomain(double a, double b, double x);
        static MultiVariatePoint<double> getCoordsFromString(string s);
        static string replace(string strs, string str_old, string str_new);
        static string eraseExtraSpaces(string strs);

        // Intermediate function for Leja points computation
        static bool isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold);
        static double computeNewLejaPointFromSequence(vector<double> seq);

        // Write data in external file
        static vector<double> loadLejaSequenceFromFile(int length);
        static vector<double> createSequenceByDichotomy(int length);
        static void storeDichotomySequenceInFile(int length);
        static void storeLejaSequenceInFile(int length);

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

        // Comparaison function
        static bool equals(MultiVariatePoint<string> nu1, MultiVariatePoint<string> nu2);
};

#endif
