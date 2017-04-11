#ifndef LAGRANGEINTERPOLATION2D
#define LAGRANGEINTERPOLATION2D

#include <iostream>
#include <vector>
#include <array>
#include <cmath>

using namespace std;

typedef array<int,2> indice2D;

class LagrangeInterpolation2D
{
    private:
        vector<double> m_pointsX;
        vector<double> m_pointsY;

    public:
        vector<indice2D> m_indices;
        vector<vector<double>> allLiMoins1Fi;

        LagrangeInterpolation2D(int sizeX, int sizeY);
        ~LagrangeInterpolation2D();

        const vector<double>& pointsX() { return m_pointsX; };
        void setPointsX(vector<double> points) { m_pointsX = points; };

        const vector<double>& pointsY() { return m_pointsY; };
        void setPointsY(vector<double> points) { m_pointsY = points; };

        void fillIndicesWithoutMax(int maxI, int maxJ);

        double g(double x, double y);
        double lagrangeBasisFunction_1D(int j, int k, double y, int axis);

        double lagrangeInterpolation_2D_simple(double x, double y);
        double lagrangeInterpolation_2D_iterative(double x, double y, int k1, int k2);
        void computeLiMinus1Fi(int k1, int k2);

};

#endif
