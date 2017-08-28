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

// Classe qui hérite de Interpolation
// Interpolation utilisant une combinaison des deux versions (LagrangeInterpolation, PiecewiseInterpolation). Chaque direction correspond à une des 3 version (voir ligne 27).
// Le type d'ordre est une chaine de caractère (code de Huffman) correspond (en 1d) au chemin du nœud dans l'arbre binaire (si version 1 ou 2)
// ou a la conversion de l'indice qui correspond à l'ordre du point de Leja (si version 0)
class MixedInterpolation : public Interpolation<string>
{
    private:
      	vector<double> m_lejaSequence;
 	// 0 ou 1 ou 2 sur chaque variable
	// 1 : Polynômes de Lgrange définis globalement + points de Leja
	// 1 : Fonctions affines par morceaux + construction des points d'interpolation par dichotomie
	// 2 : Fonctions quadratiques par morceaux + construction des points d'interpolation par dichotomie
      	MultiVariatePoint<int> m_methods;
	// Vecteur d'arbre: 1 arbre par direction (si version 1 ou 2)
      	vector<BinaryTreePtr> m_trees;

    public:
      	MixedInterpolation(FunctionsPtr f, int nIter, MultiVariatePoint<int> methods);
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
};

typedef std::unique_ptr<MixedInterpolation> MixedInterpolationPtr;

#endif
