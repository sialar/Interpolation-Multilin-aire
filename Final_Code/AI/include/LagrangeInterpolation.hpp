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

// Classe qui hérite de Interpolation
// Interpolation utilisant des polynomes de Lagrange définis globalement et les points de Leja
// Le type d'ordre est un entier (multindice en grande dimension) correspond à l'indice du point dans la séquence de Leja
class LagrangeInterpolation : public Interpolation<int>
{
    private:
      // Vecteur Contenant des points de Leja 
      vector<double> m_lejaSequence;

    public:
      LagrangeInterpolation(int d, int n, int nIter);
      ~LagrangeInterpolation() {};


      /************************* Data points ********************************/
      MultiVariatePoint<double> getPoint(MultiVariatePointPtr<int> nu);
      void addInterpolationPoint(MultiVariatePoint<double> p);

      /************************* AI algo ************************************/
      MultiVariatePointPtr<int> getFirstMultivariatePoint();
      MultiVariatePointPtr<int> maxElement(int iteration);
      bool indiceInPath(MultiVariatePoint<int> index);
      void updateCurentNeighbours(MultiVariatePointPtr<int> nu);
      bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<int> nu);

      /*********************** Interpolation ********************************/
      double basisFunction_1D(int code, double t, int axis);
};

typedef std::unique_ptr<LagrangeInterpolation> LagrangeInterpolationPtr;

#endif
