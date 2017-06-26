#ifndef MIXEDINTERPOLATION
#define MIXEDINTERPOLATION

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

class MixedInterpolation : public Interpolation<string>
{
    private:
      vector<double> m_lejaSequence;
      MultiVariatePoint<int> m_methods; // 0 ou 1 ou 2 sur chaque variable
      map<MultiVariatePoint<int>,vector<double>> m_methods_errors;
      vector<BinaryTreePtr> m_trees;

    public:
      MixedInterpolation(int d, int n, int nIter, MultiVariatePoint<int> methods);
      ~MixedInterpolation() {};
      void clearAllTrees();

      /************************* Data points ********************************/
      MultiVariatePoint<double> getPoint(MultiVariatePointPtr<string> nu);
      void addInterpolationPoint(MultiVariatePoint<double> p);
      void computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis);
      void setMethods(MultiVariatePoint<int> methods) { m_methods = methods; }

      /************************* AI algo ************************************/
      MultiVariatePointPtr<string> getFirstMultivariatePoint();
      MultiVariatePointPtr<string> maxElement(int iteration);
      bool indiceInPath(MultiVariatePoint<string> index);
      void updateCurentNeighbours(MultiVariatePointPtr<string> nu);
      bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<string> nu);

      /*********************** Interpolation ********************************/
      double basisFunction_1D(string code, double t, int axis);
      vector<double> tryWithDifferentMethods(MultiVariatePoint<int> methods, double threshold);
      MultiVariatePoint<int> tryAllCases(double threshold);

      /********************** Display functions *****************************/
      void savePathInFile(string fileName);
      void saveInterpolationBasisFunctions();
};

typedef std::unique_ptr<MixedInterpolation> MixedInterpolationPtr;

#endif
