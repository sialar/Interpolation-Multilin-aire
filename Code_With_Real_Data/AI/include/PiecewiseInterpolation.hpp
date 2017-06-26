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
#include "Interpolation.hpp"
#include "BinaryTree.hpp"
#include "Utils.hpp"

using namespace std;

class PiecewiseInterpolation : public Interpolation<string>
{
    private:
        int m_method; // 1 ou 2
        vector<BinaryTreePtr> m_trees;

    public:
        PiecewiseInterpolation(int d, int n, int nIter, int method);
        ~PiecewiseInterpolation() {};
        void clearAllTrees();

        /************************* Data points ********************************/
        MultiVariatePoint<double> getPoint(MultiVariatePointPtr<string> nu);
        void addInterpolationPoint(MultiVariatePoint<double> p);
        void computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis);

        /************************* AI algo ************************************/
        MultiVariatePointPtr<string> getFirstMultivariatePoint();
        MultiVariatePointPtr<string> maxElement(int iteration);
        bool indiceInPath(MultiVariatePoint<string> index);
        void updateCurentNeighbours(MultiVariatePointPtr<string> nu);
        bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<string> nu);

        /*********************** Interpolation ********************************/
        double basisFunction_1D(string code, double t, int axis);

        /********************** Display functions *****************************/
        void displayTrees();
        void saveInterpolationBasisFunctions();
};

typedef std::unique_ptr<PiecewiseInterpolation> PiecewiseInterpolationPtr;

#endif
