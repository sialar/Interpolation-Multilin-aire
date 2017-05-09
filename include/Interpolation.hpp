#ifndef INTERPOLATIONND
#define INTERPOLATIONND

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <sys/time.h>
#include <ctime>
#include <chrono>

#include "MultiVariatePoint.hpp"
#include "Utils.hpp"

using namespace std;

class Interpolation
{
    private:

        int m_d;
        map<MultiVariatePoint<int>, double> m_alphaMap;
        vector<vector<double>> m_points;
        vector<MultiVariatePoint<double>> m_testPoints;
        vector<MultiVariatePoint<int>> m_path;
        list<MultiVariatePoint<int>> m_curentNeighbours;

    public:

        Interpolation(vector<int> sizes);
        ~Interpolation();

        /************************* Data points ********************************/
        const vector<vector<double>>& points() { return m_points; };
        MultiVariatePoint<double> getPoint(MultiVariatePoint<int> nu);
        void setDirPoints(int i, vector<double> pointsI) { m_points[i] = pointsI; };
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };

        /************************* AI algo ************************************/
        const vector<MultiVariatePoint<int>>& path() { return m_path; };
        int buildPathWithAIAlgo(int maxIteration, auto start_time, double threshold);
        bool indiceInPath(MultiVariatePoint<int>& index);
        double computeLastAlphaNu(MultiVariatePoint<int>& nu);
        void updateCurentNeighbours(MultiVariatePoint<int>& nu);
        bool isCorrectNeighbourToCurentPath(MultiVariatePoint<int>& nu);

        /*********************** Interpolation ********************************/
        double piecewiseFunction(int k, double t, int axis);
        double lagrangeBasisFunction_1D(int j, int k, double t, int axis);
        double lagrangeInterpolation_ND(MultiVariatePoint<double> x);

        /********************** Test functions ********************************/
        void testPathBuilt(int maxIteration, double threshold);
        double testAlphaNuComputation(MultiVariatePoint<int>& nu);
        double tryWithCurentPath();

        /********************** Display functions *****************************/
        void displayPath();
        void displayAlphaTab();
        void displayCurentNeighbours();
        void displayInterpolationPoints();
        void savePathInFile();
};

#endif
