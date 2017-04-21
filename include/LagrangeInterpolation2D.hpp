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
#include "IndiceND.hpp"

using namespace std;


class LagrangeInterpolation2D
{
    private:
        double m_alphaInitVal;
        vector<double> m_pointsX;
        vector<double> m_pointsY;
        vector<IndiceND> m_path;
        vector<IndiceND> m_curSetInAIAlgo;
        vector<vector<double>> m_alphaTab;

    public:

        LagrangeInterpolation2D(int sizeX, int sizeY,int path);
        ~LagrangeInterpolation2D();

        void clear();

        const vector<double>& pointsX() { return m_pointsX; };
        void setPointsX(vector<double> points) { m_pointsX = points; };
        const vector<double>& pointsY() { return m_pointsY; };
        void setPointsY(vector<double> points) { m_pointsY = points; };
        void setAlphaInitVal(double val) { m_alphaInitVal = val; };

        const vector<IndiceND>& path() { return m_path; };
        void choosePath(int n, int m, int v /* 0, 1 ou 2 */);
        int getIndiceInPath(int maxI, int maxJ);
        void initAlphaTab(double initVal);
        void showAlphaTab();
        void savePathInFile();
        void showPath();

        void buildPathWithAIAlgo(int k);
        double testPathBuilt(int nbIteration);
        vector<IndiceND> getCurentNeighbours();
        bool indexInPath(IndiceND& index);

        double lagrangeBasisFunction_1D(int j, int k, double y, int axis);
        double lagrangeInterpolation_2D_simple(double x, double y);

        void computeAllAlphaNu();
        double computeOneAlphaNu(IndiceND& nu);
        double computeLastAlphaNu(IndiceND& nu);
        double lagrangeInterpolation_2D_iterative(double x, double y);
        double computeExecTimeOfOneApprox();

};

#endif
