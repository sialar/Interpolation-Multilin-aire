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
      map<MultiVariatePoint<int>,double> m_methods_errors;
      vector<BinaryTreePtr> m_trees;

    public:
      MixedInterpolation(int d, int nIter, MultiVariatePoint<int> methods);
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
      bool indiceInNeighborhood(MultiVariatePoint<string> index);
      void updateCurentNeighbours(MultiVariatePointPtr<string> nu);
      bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<string> nu);

      /*********************** Interpolation ********************************/
      double basisFunction_1D(string code, double t, int axis);
      double tryWithDifferentMethods(MultiVariatePoint<int> methods, double threshold, int function);
      MultiVariatePoint<int> tryAllCases(double threshold, int function);

      /********************** Display functions *****************************/
      void savePathInFile();
      void storeInterpolationBasisFunctions();
};

typedef std::unique_ptr<MixedInterpolation> MixedInterpolationPtr;

#endif
