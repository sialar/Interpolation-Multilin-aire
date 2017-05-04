#ifndef INTERPOLATIONND
#define INTERPOLATIONND

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <sys/time.h>
#include <ctime>
#include <omp.h>

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

        void clear();

        /************************* Data points ********************************/
        const vector<vector<double>>& points() { return m_points; };
        MultiVariatePoint<double> getPoint(MultiVariatePoint<int> nu);
        void setDirPoints(int i, vector<double> pointsI) { m_points[i] = pointsI; };
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };
        void smartDiscretization();
        void displayPoints();

        /************************* Path ***************************************/
        const vector<MultiVariatePoint<int>>& path() { return m_path; };
        int getLastIndiceInPath(MultiVariatePoint<int> max);
        int buildPathWithAIAlgo(int maxIteration, double start_time, double threshold);
        double testPathBuilt(int maxIteration, double threshold);
        bool indiceInPath(MultiVariatePoint<int>& index);
        double tryWithCurentPath();
        void savePathInFile();
        void displayPath();

        /************************* Alpha **************************************/
        double computeLastAlphaNu(MultiVariatePoint<int>& nu);
        double testAlphaNuComputation(MultiVariatePoint<int>& nu);
        void displayAlphaTab();

        /************************ Neighbours **********************************/
        void updateCurentNeighbours(MultiVariatePoint<int>& nu);
        bool isCorrectNeighbourToCurentPath(MultiVariatePoint<int>& nu);
        void displayCurentNeighbours();

        /*********************** Interpolation ********************************/
        double lagrangeBasisFunction_1D(int j, int k, double t, int axis);
        double lagrangeInterpolation_ND_iterative(MultiVariatePoint<double> x);
        double computeExecTimeOfOneApprox();
        void displayInterpolationPoints();

};

#endif
