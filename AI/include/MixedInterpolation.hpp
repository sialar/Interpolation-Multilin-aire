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

/**
 *  \file MixedInterpolation.hpp
 *  \brief Classe dérivée de Interpolation<string>: interpolation utilisant une combinaison des deux versions (LagrangeInterpolation, PiecewiseInterpolation). Chaque direction correspond à une des 3 version (voir README ligne 4).
 *  Le type d'ordre est une chaine de caractère (code de Huffman) qui correspond (en 1d), soit au chemin du nœud dans l'arbre binaire (si version 1 ou 2),
 *  soit à l'indice du point de Leja (si version 0)
 *  \author SIALA Rafik
 *  \date 08/16
*/
class MixedInterpolation : public Interpolation<string>
{
    private:
        // Séquence de points de Leja lu a partir du fichier leja_sequence.dat
      	vector<double> m_lejaSequence;
       	// 0 ou 1 ou 2 sur chaque variable (README ligne 4)
      	// 1 : Polynômes de Lgrange définis globalement + points de Leja
      	// 1 : Fonctions affines par morceaux + construction des points d'interpolation par dichotomie
      	// 2 : Fonctions quadratiques par morceaux + construction des points d'interpolation par dichotomie
      	MultiVariatePoint<int> m_methods;
	      // Liste d'arbres: 1 arbre par direction (si version 1 ou 2)
      	vector<BinaryTreePtr> m_trees;

    public:
        /**
          * Constructeur
          * \param f : fonction à interpoler
          * \param nIter : nombre d'itérations dans l'algorithme AI
          * \param methods : indique la version utilisée par l'algorithme AI dans chaque direction
        */
      	MixedInterpolation(FunctionsPtr f, int nIter, MultiVariatePoint<int> methods);
      	~MixedInterpolation() {};
        /**
          * Mutateur de l'attribut m_methods
          * param methods : liste des versions utilisées sur chaque directions
        */
        void setMethods(MultiVariatePoint<int> methods) { m_methods = methods; }

/******************************************************************************/
/************************ Points d'interpolation ******************************/
        /**
          * Supprimer tous les arbres de m_trees
        */
        void clearAllTrees();
        /**
          * Retourne le point multivarié correspondant au vecteur d'ordres nu
          * \param nu : vecteur contenant des codes de Huffman (chemin du point 1d dans l'arbre binaire) ou des indices (indice du point 1d de Leja) sur chaque direction
          * \return point multivarié correspondant à nu
        */
      	MultiVariatePoint<double> getPoint(MultiVariatePointPtr<string> nu);
        /**
          * Ajoute le nouveau point d'interpolation p dans m_interpolationPoints, m_interpolationNodes et m_trees
          * \param p : nouveau point d'interpolation
        */
      	void addInterpolationPoint(MultiVariatePoint<double> p);

/******************************************************************************/
/*************************** Algorithme AI ************************************/
        /**
          * Retourner le premier vecteur d'ordres dans l'algo AI
          * \return le premier vecteur d'ordres (codes de Huffman ou indice de point de Leja) dans l'algo AI
        */
        MultiVariatePointPtr<string> getFirstMultivariatePoint();
        /**
          * Retourner le vecteur d'ordres du prochain point d'interpolation sélectionné par l'algo AI
          * 3 fois sur 4: on choisi celui ayant la plus grande erreur (algo AI)
          * 1 fois sur 4: on choisi celui qui a attendu le plus longtemps pour éviter de bloquer une direction dans le cas où l'erreur vaut 0
          * \param iteration : numéro de l'itération courante
          * \return teration : vecteur d'ordres du prochain point d'interpolation
        */
      	MultiVariatePointPtr<string> maxElement(int iteration);
        /**
          * Vérifier si le vecteur d'ordre (codes de Huffman, ou indices) index existe dans m_path (correspond à un point d'interpolation)
          * \param index : vecteur d'ordres
          * \return true si et seulement si index est dans m_path
        */
      	bool indiceInPath(MultiVariatePoint<string> index);
        /**
          * Mettre à jour, en fonction du nouveau point d'interpolation correspondant au vecteur d'ordre nu, la liste des candidats courant
          * \param nu : nouveau vecteur d'ordre choisi par l'algo AI
        */
      	void updateCurentNeighbours(MultiVariatePointPtr<string> nu);
        /**
          * Vérifier si nu est bien un voisin de la séquence de points d'interpolation courante (vérifier la monotonie de l'ensemble de points choisis)
          * \param nu : vecteur d'ordre
          * \return true si et seulement si nu est voisin à l'ensemble de points de m_path
        */
      	bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<string> nu);

/******************************************************************************/
/************************** Fonctions de base *********************************/
        /**
          * Implémente la fonction de base (polynome de Lagrange global, fonction affine (ou quadratique) par morceaux (dépend de la version de l'algo AI)
          * \param code : ordre (indice ou code de Huffman) correspondant à une coordonnée d'un point d'interpolation
          * \param t : point d'évaluation (1d)
          * \param axis : direction concidéré (0 <= axis <= d-1)
          * \return la valeur au point t de la fonction de base correspondant au point 1d d'ordre code sur la direction axis
        */
        double basisFunction_1D(string code, double t, int axis);
        /**
          * Stocker les valeurs des fonctions de bases dans basis_functions.dat
        */
        void saveInterpolationBasisFunctions();
        /**
          * Rechercher les valeurs des points d'interpolation 1d (inf et sup) les plus proches de t (de part et d'autre de t) sur la direction axis pour construire la fonction de base
          * \param t : point d'interpolation 1d
          * \param inf : plus grand point d'interpolation 1d inférieur à t
          * \param sup : plus petit point d'interpolation 1d supérieur à t
          * \param axis : direction concidérée
        */
        void computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis);
};

typedef std::unique_ptr<MixedInterpolation> MixedInterpolationPtr;

#endif
