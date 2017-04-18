#ifndef LAGRANGEINTERPOLATION2D
#define LAGRANGEINTERPOLATION2D

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>

using namespace std;

typedef array<int,2> indice2D;

class LagrangeInterpolation2D
{
    private:
        vector<double> m_pointsX;
        vector<double> m_pointsY;
        vector<indice2D> m_path;
        vector<vector<double>> m_alphaTab;

        vector<indice2D> m_curSetInAIAlgo;

    public:

        LagrangeInterpolation2D(int sizeX, int sizeY,int path);
        ~LagrangeInterpolation2D();

        const vector<double>& pointsX() { return m_pointsX; };
        void setPointsX(vector<double> points) { m_pointsX = points; };
        const vector<double>& pointsY() { return m_pointsY; };
        void setPointsY(vector<double> points) { m_pointsY = points; };

        const vector<indice2D>& path() { return m_path; };
        void choosePath(int n, int m, int v /* 0, 1 ou 2 */);
        int getCurentIndice(int maxI, int maxJ);
        void showPath();


        void buildPathWithAIAlgo(int n, int m, bool lejaSeq);
        vector<indice2D> getCurentNeighbours();
        bool indexInPath(indice2D index);
        void showAlphaTab();

        double g(double x, double y);
        double lagrangeBasisFunction_1D(int j, int k, double y, int axis);
        double lagrangeInterpolation_2D_simple(double x, double y);

        void computeAllAlphaNu(int k1, int k2);
        double lagrangeInterpolation_2D_iterative(double x, double y, int k1, int k2);

};

#endif
