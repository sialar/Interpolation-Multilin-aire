#ifndef LAGRANGEINTERPOLATION2D
#define LAGRANGEINTERPOLATION2D

#include <iostream>
#include <vector>
#include <array>
#include <cmath>

using namespace std;

class LagrangeInterpolation2D
{
    private:
        vector<double> m_pointsX;
        vector<double> m_pointsY;
        vector<double> allLiMoins1Fi;

    public:

        LagrangeInterpolation2D(int sizeX, int sizeY);
        ~LagrangeInterpolation2D();

        const vector<double>& pointsX() { return m_pointsX; };
        void setPointsX(vector<double> points) { m_pointsX = points; };

        const vector<double>& pointsY() { return m_pointsY; };
        void setPointsY(vector<double> points) { m_pointsY = points; };

        double g(double x, double y);
        double lagrangeBasisFunction_1D(int j, int k, double y, int axis);

        double lagrangeInterpolation_2D_simple(double x, double y);
};

#endif
