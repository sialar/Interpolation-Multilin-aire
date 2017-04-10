#ifndef LAGRANGEINTERPOLATIONND
#define LAGRANGEINTERPOLATIONND

#include <iostream>
#include <vector>
#include <array>
#include <cmath>

using namespace std;

class LagrangeInterpolationND
{
    private:
        vector<double> m_points;

    public:
        LagrangeInterpolationND(int size);
        ~LagrangeInterpolationND();

        const vector<double>& points() { return m_points; };
        void setPoints(vector<double> points) { m_points = points; };
};

#endif
