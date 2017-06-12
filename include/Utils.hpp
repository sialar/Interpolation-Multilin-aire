#ifndef UTILS
#define UTILS

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>

#include "MultiVariatePoint.hpp"

using namespace std;

class Utils
{
    private:
        static vector<double> m_1dGrid;

    public:
      static vector<vector<vector<double>>> m_coefs;
        Utils()  {};
        ~Utils() {};

        // Display data and results
        static void separateur();
        static void displayPoints(vector<vector<double>> v);
        static void displayPoints(vector<MultiVariatePoint<double>> v);

        // Useful functions
        static double interpolationError(vector<vector<double>> realValue, vector<vector<double>> estimate);
        static MultiVariatePoint<double> createRandomMultiVariatePoint(int d);
        static vector<double> toVector(double x, int n);
        static vector<double> diff(vector<double> x,vector<double> y);
        static double randomValue(double a, double b);
        static string vector2str(vector<double> x);
        static double norm(vector<double> x, int p);

        // Intermediate function for Leja points computation
        static bool isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold);
        static double computeNewLejaPointFromSequence(vector<double> seq);

        // Write data in external file
        static vector<double> loadLejaSequenceFromFile(int length);
        static vector<double> createSequenceByDichotomy(int length);
        static void storeDichotomySequenceInFile(int length);
        static void storeLejaSequenceInFile(int length);

        // Function to approximate
        static void setCoefs(int degree, int d, int n);
        static vector<double> g(MultiVariatePoint<double> x, int n);
        static vector<double> f(MultiVariatePoint<double> x, int n);
        static vector<double> polynomial(MultiVariatePoint<double> x, int degree, int n);

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

        // Read arguments
        static int chooseDimensionD(int argc, char* argv[], int argNum);
        static int chooseDimensionN(int argc, char* argv[], int argNum);
        static int chooseNbTestPoints(int argc, char* argv[], int argNum);
        static int chooseMaxIteration(int argc, char* argv[], int argNum);
        static int chooseMethod(int argc, char* argv[], int argNum);
        static bool saveResults(int argc, char* argv[], int argNum);
        static bool saveError(int argc, char* argv[], int argNum);
        static bool plotPath(int argc, char* argv[], int argNum);
        static bool displayResults();

};

#endif
