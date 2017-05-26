#ifndef LAGRANGEINTERPOLATION
#define LAGRANGEINTERPOLATION

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <memory>
#include <cmath>
#include <chrono>

#include "MultiVariatePoint.hpp"
#include "BinaryTree.hpp"
#include "Utils.hpp"

using namespace std;

class LagrangeInterpolation
{
    private:
        int m_method;
        int m_maxIteration;
        int m_d;

        vector<double> m_lejaSequence;
        vector<vector<double>> m_interpolationPoints;
        vector<MultiVariatePoint<double>> m_interpolationNodes;
        vector<MultiVariatePoint<double>> m_testPoints;

        vector<MultiVariatePointPtr<int>> m_path;
        map<MultiVariatePointPtr<int>, double> m_alphaMap;
        list<MultiVariatePointPtr<int>> m_curentNeighbours;

    public:
        LagrangeInterpolation(int d, int nIter);
        ~LagrangeInterpolation() {};

        /************************* Data points ********************************/
        const vector<vector<double>>& points() { return m_interpolationPoints; };
        const vector<MultiVariatePoint<double>>& interpolationPoints() { return m_interpolationNodes; };
        MultiVariatePoint<double> getPoint(MultiVariatePointPtr<int> nu);
        void addInterpolationPoint(MultiVariatePoint<double> p);
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };

        /************************* AI algo ************************************/
        const vector<MultiVariatePointPtr<int>>& path() { return m_path; };
        int buildPathWithAIAlgo(auto start_time, double threshold, bool debug);
        bool indiceInPath(MultiVariatePoint<int> index);
        bool indiceInNeighborhood(MultiVariatePoint<int> index);
        double computeLastAlphaNu(MultiVariatePointPtr<int> nu);
        void updateCurentNeighbours(MultiVariatePointPtr<int> nu);
        bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<int> nu);

        /*********************** Interpolation ********************************/
        double lagrangeBasisFunction_1D(int k, double t, int axis);
        double interpolation_ND(MultiVariatePoint<double>& x, int end);

        /********************** Test functions ********************************/
        void testPathBuilt(double threshold, bool debug);
        double tryWithCurentPath();

        /********************** Display functions *****************************/
        void displayPath();
        void displayAlphaTab();
        void displayCurentNeighbours();
        void displayInterpolationPointsInEachDirection();
        void displayInterpolationMultiVariatePoints();
        void storeInterpolationBasisFunctions();
        void storeInterpolationProgression();
        void savePathInFile();

};

typedef std::unique_ptr<LagrangeInterpolation> LagrangeInterpolationPtr;

#endif
