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
#include "Interpolation.hpp"
#include "BinaryTree.hpp"
#include "Utils.hpp"

using namespace std;

class LagrangeInterpolation : public Interpolation<int>
{
    private:
      vector<double> m_lejaSequence;

    public:
      LagrangeInterpolation(int d, int nIter);
      ~LagrangeInterpolation() {};


      /************************* Data points ********************************/
      MultiVariatePoint<double> getPoint(MultiVariatePointPtr<int> nu);
      void addInterpolationPoint(MultiVariatePoint<double> p);

      /************************* AI algo ************************************/
      MultiVariatePointPtr<int> getFirstMultivariatePoint();
      MultiVariatePointPtr<int> maxElement(int iteration);
      bool indiceInPath(MultiVariatePoint<int> index);
      bool indiceInNeighborhood(MultiVariatePoint<int> index);
      void updateCurentNeighbours(MultiVariatePointPtr<int> nu);
      bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<int> nu);

      /*********************** Interpolation ********************************/
      double basisFunction_1D(int code, double t, int axis);

      /********************** Display functions *****************************/
      void storeInterpolationBasisFunctions();
};

typedef std::unique_ptr<LagrangeInterpolation> LagrangeInterpolationPtr;

#endif
