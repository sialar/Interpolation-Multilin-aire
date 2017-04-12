#ifndef UTILS
#define UTILS

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <cmath>

using namespace std;

class Function
{
    private:
        string expression;
        int m_d, m_n;
    public:
        Function(int d, int n) : m_d(d), m_n(n) {};
        ~Function() {};

        const string& getExpression() { return expression; };
        void setExpression(string expr) { expression = expr; };
};

class Utils
{
    public:
        Utils()  {};
        ~Utils() {};

        static double randomValue(double a, double b);
        static double squareError(vector<double> realValue, vector<double> estimate);

        // For uniform sequence of data points
        static vector<double> createUniformSequence(int nbPoints);
        static void binaryDecomposition(int number, vector<double>& binary_decomp);

        // For Leja sequence (L)
        static vector<double> createLejaSequence(int nbPoints);

        // For Tchebychev zeros
        static vector<double> createChebychevSequence(int nbPoints);



};

#endif
