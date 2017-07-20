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

enum method { Apollo, Cocagne, Tucker, Tucker_bis, AI };

class Utils
{
    private:
        static vector<double> m_1dGrid;

    public:
        static string projectPath;
        static void separateur();
        static void displayValues(vector<double> v);
        static vector<double> diff(vector<double> x,vector<double> y);
        static double randomValue(double a, double b);
        static string vector2str(vector<double> x);
        static vector<double> str2vector(string line);
        static double norm(vector<double> x, int p);
        static double maxAbsValue(vector<double> vec);
        static double computeMseError(vector<double> f,vector<double> f_tilde);
        static bool strInVector(string required, vector<string> vec);
        static MultiVariatePoint<double> createRandomMultiVariatePoint(int d);
        static double adaptCoordsToFunctionDomain(double a, double b, double x);
        static double convertToFunctionDomain(double a, double b, double x);
        static string replace(string strs, string str_old, string str_new);
        static string eraseExtraSpaces(string strs);
        static vector<double> createUniformSequence(int nbPoints);
        static vector<double> createChebychevSequence(int nbPoints);
        static vector<double> createLejaSequence(int nbPoints);
        static vector<double> createSequenceOfMiddles(int nbPoints);
        static vector<double> loadLejaSequenceFromFile(int length);
        static vector<double> createSequenceByDichotomy(int length);
        static void storeDichotomySequenceInFile(int length);
        static void storeLejaSequenceInFile(int length);
        static void binaryDecomposition(int number, vector<double>& binary_decomp);
        static bool isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold);
        static bool equals(MultiVariatePoint<string> nu1, MultiVariatePoint<string> nu2);
        static double computeNewLejaPointFromSequence(vector<double> seq);
};

#endif
