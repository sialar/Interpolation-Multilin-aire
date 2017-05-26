#ifndef PIECEWISEINTERPOLATION
#define PIECEWISEINTERPOLATION

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

class PiecewiseInterpolation
{
    private:
        int m_method; // 1 ou 2
        int m_maxIteration;
        int m_d;

        vector<double> m_lejaSequence;
        vector<vector<double>> m_interpolationPoints;
        vector<MultiVariatePoint<double>> m_interpolationNodes;
        vector<MultiVariatePoint<double>> m_testPoints;

        vector<MultiVariatePointPtr<string>> m_path;
        map<MultiVariatePointPtr<string>, double> m_alphaMap;
        list<MultiVariatePointPtr<string>> m_curentNeighbours;

        vector<BinaryTreePtr> m_trees;

    public:
        PiecewiseInterpolation(int d, int nIter, int method);
        ~PiecewiseInterpolation() {};

        /************************* Data points ********************************/
        const vector<vector<double>>& points() { return m_interpolationPoints; };
        const vector<MultiVariatePoint<double>>& interpolationPoints() { return m_interpolationNodes; };
        MultiVariatePoint<double> getPoint(MultiVariatePointPtr<string> nu);
        void addInterpolationPoint(MultiVariatePoint<double> p);
        void setTestPoints(vector<MultiVariatePoint<double>> points) { m_testPoints = points; };
        void computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis);

        /************************* AI algo ************************************/
        const vector<MultiVariatePointPtr<string>>& path() { return m_path; };
        int buildPathWithAIAlgo(auto start_time, double threshold, bool debug);
        bool indiceInPath(MultiVariatePoint<string> index);
        bool indiceInNeighborhood(MultiVariatePoint<string> index);
        double computeLastAlphaNu(MultiVariatePointPtr<string> nu);
        void updateCurentNeighbours(MultiVariatePointPtr<string> nu);
        bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<string> nu);

        /*********************** Interpolation ********************************/
        double piecewiseFunction_1D(string code, double t, int axis);
        double quadraticFunction_1D(string code, double t, int axis);
        double interpolation_ND(MultiVariatePoint<double>& x, int end);
        void setMethod(int method) { m_method = method;};

        /********************** Test functions ********************************/
        void testPathBuilt(double threshold, bool debug);
        double tryWithCurentPath();

        /********************** Display functions *****************************/
        void displayPath();
        void displayAlphaTab();
        void displayCurentNeighbours();
        void displayInterpolationPointsInEachDirection();
        void displayInterpolationMultiVariatePoints();
        void displayTrees();
        void storeInterpolationBasisFunctions();
        void storeInterpolationProgression();
        void savePathInFile();
};

typedef std::unique_ptr<PiecewiseInterpolation> PiecewiseInterpolationPtr;

#endif
