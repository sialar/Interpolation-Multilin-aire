#ifndef LAGRANGEINTERPOLATION1D
#define LAGRANGEINTERPOLATION1D

#include <iostream>
#include <vector>
#include <array>
#include <cmath>

using namespace std;

class LagrangeInterpolation1D
{
    private:
        vector<double> m_points;
        vector<double> m_alphaTab;

    public:

        LagrangeInterpolation1D(int size);
        ~LagrangeInterpolation1D();

        const vector<double>& points() { return m_points; };
        void setPoints(vector<double> points) { m_points = points; };

        double g(double y);
        double lagrangeBasisFunction_1D(int j, int k, double y);
        void computeAllAlphaI(int k);

        double lagrangeInterpolation_1D_simple(double y);
        double lagrangeInterpolation_1D_recursive(double y, int k);
        double lagrangeInterpolation_1D_iterative(double y, int k);
};
#endif
