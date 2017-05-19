#ifndef DICHOTOMY
#define DICHOTOMY

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>

using namespace std;

class Dichotomy
{

    public:
        Dichotomy()  {};
        ~Dichotomy() {};

        static int getIndice(double l);
        static double getValue(int i);
        static int getParentIndice(int i);
        static double getParentValue(int i);
        static vector<double> computeChildrenValue(int i);
        static bool isAncestorOf(int i, int j);
};

#endif
