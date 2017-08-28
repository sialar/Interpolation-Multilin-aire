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

// Classe qui hérite de Interpolation
// Interpolation utilisant des fonctions définies par morceaux et une construction des points d'interpolation par dichotomie.
// Le type d'ordre est une chaine de caractère (code de Huffman) correspond (en 1d) au chemin du nœud dans l'arbre binaire.

class PiecewiseInterpolation : public Interpolation<string>
{
    private:
	// m_method = 1 ou 2 (version)
	// 1 : Fonctions affines par morceaux + construction des points d'interpolation par dichotomie
	// 2 : Fonctions quadratiques par morceaux + construction des points d'interpolation par dichotomie
        int m_method;
	// Vecteur d'arbre: 1 arbre par direction
        vector<BinaryTreePtr> m_trees;

    public:
        PiecewiseInterpolation(FunctionsPtr f, int nIter, int method);
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
};

typedef std::unique_ptr<PiecewiseInterpolation> PiecewiseInterpolationPtr;

#endif
