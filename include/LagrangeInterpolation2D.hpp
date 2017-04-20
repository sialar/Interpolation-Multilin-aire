#ifndef LAGRANGEINTERPOLATION2D
#define LAGRANGEINTERPOLATION2D

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <time.h>

#include "Utils.hpp"

using namespace std;

typedef array<int,2> indice2D;

class LagrangeInterpolation2D
{
    private:
        double m_alphaInitVal;
        vector<double> m_pointsX;
        vector<double> m_pointsY;
        vector<indice2D> m_path;
        vector<indice2D> m_curSetInAIAlgo;
        vector<vector<double>> m_alphaTab;

    public:

        LagrangeInterpolation2D(int sizeX, int sizeY,int path);
        ~LagrangeInterpolation2D();

        void clear();
        static void displayVector(vector<indice2D> vect);

        const vector<double>& pointsX() { return m_pointsX; };
        void setPointsX(vector<double> points) { m_pointsX = points; };
        const vector<double>& pointsY() { return m_pointsY; };
        void setPointsY(vector<double> points) { m_pointsY = points; };
        void setAlphaInitVal(double val) { m_alphaInitVal = val; };

        const vector<indice2D>& path() { return m_path; };
        void choosePath(int n, int m, int v /* 0, 1 ou 2 */);
        int getIndiceInPath(int maxI, int maxJ);
        void initAlphaTab(double initVal);
        void showAlphaTab();
        void savePathInFile();
        void showPath();

        void buildPathWithAIAlgo(int k);
        double testPathBuilt(int nbIteration);
        vector<indice2D> getCurentNeighbours();
        bool indexInPath(indice2D index);

        double lagrangeBasisFunction_1D(int j, int k, double y, int axis);
        double lagrangeInterpolation_2D_simple(double x, double y);

        void computeAllAlphaNu();
        double computeOneAlphaNu(indice2D nu);
        double computeLastAlphaNu(indice2D nu);
        double lagrangeInterpolation_2D_iterative(double x, double y);
        double computeExecTimeOfOneApprox();

};

#endif
