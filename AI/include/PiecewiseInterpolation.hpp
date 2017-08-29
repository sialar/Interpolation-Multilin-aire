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

/**
 *  \file PiecewiseInterpolation.hpp
 *  \brief Classe dérivée de Interpolation<string>: interpolation utilisant des fonctions définies par morceaux et une construction des points d'interpolation par dichotomie.
 *  Le type d'ordre est une chaine de caractère (code de Huffman) qui correspond (en 1d) au chemin du nœud dans l'arbre binaire.
 *  \author SIALA Rafik
 *  \date 08/16
*/
class PiecewiseInterpolation : public Interpolation<string>
{
    private:
      	// m_method = 1 ou 2 (version)
      	// 1 : Fonctions affines par morceaux + construction des points d'interpolation par dichotomie
      	// 2 : Fonctions quadratiques par morceaux + construction des points d'interpolation par dichotomie
        int m_method;
	      // Liste d'arbres: 1 arbre par direction
        vector<BinaryTreePtr> m_trees;

    public:
        /**
          * Constructeur
          * \param f : fonction à interpoler
          * \param nIter : nombre d'itérations dans l'algorithme AI
          * \param method : indique la version utilisée par l'algorithme AI
        */
        PiecewiseInterpolation(FunctionsPtr f, int nIter, int method);
        ~PiecewiseInterpolation() {};

/******************************************************************************/
/************************ Points d'interpolation ******************************/
        /**
          * Supprimer tous les arbres de m_trees
        */
        void clearAllTrees();
        /**
          * Retourne le point multivarié correspondant au vecteur d'ordres nu
          * \param nu : vecteur de codes de Huffman (chemin du point 1d dans l'arbre binaire, sur chaque direction)
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
        * \return le premier vecteur d'ordres (codes de Huffman) dans l'algo AI = ("",..,"") correspont au point multivarié (0.0, .., 0.0) (car la racine de l'arbre 0.0)
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
          * Vérifier si le vecteur d'ordre (codes de Huffman) index existe dans m_path (correspond à un point d'interpolation)
          * \param index : vecteur de codes de Huffman
          * \return true si et seulement si index est dans m_path
        */
        bool indiceInPath(MultiVariatePoint<string> index);
        /**
          * Mettre à jour, en fonction du nouveau point d'interpolation correspondant au vecteur de code de Huffman nu, la liste des candidats courant
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
          * Implémente la fonction de base (fonction affine (ou quadratique) par morceaux)
          * \param code : code de Huffman correspondant à une coordonnée d'un point d'interpolation
          * \param t : point d'évaluation (1d)
          * \param axis : direction concidéré (0 <= axis <= d-1)
          * \return la valeur au point t de la fonction de base correspondant au point 1d d'ordre code sur la direction axis
        */
        double basisFunction_1D(string code, double t, int axis);
        /**
          * Rechercher les valeurs des points d'interpolation 1d (inf et sup) les plus proches de t (de part et d'autre de t) sur la direction axis pour construire la fonction de base
          * \param t : point d'interpolation 1d
          * \param inf : plus grand point d'interpolation 1d inférieur à t
          * \param sup : plus petit point d'interpolation 1d supérieur à t
          * \param axis : direction concidérée
        */
        void computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis);
};

typedef std::unique_ptr<PiecewiseInterpolation> PiecewiseInterpolationPtr;

#endif
