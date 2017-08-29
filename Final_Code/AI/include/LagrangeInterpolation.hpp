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

/**
 *  \file LagrangeInterpolation.hpp
 *  \brief Classe dérivée de Interpolation<int>: interpolation utilisant des polynomes de Lagrange définis globalement et les points de Leja.
 *  Le type d'ordre est un entier (multindice en grande dimension) correspond à l'indice du point dans la séquence de Leja.
 *  \author SIALA Rafik
 *  \date 08/16
*/
class LagrangeInterpolation : public Interpolation<int>
{
    private:
      // Séquence de points de Leja lu a partir du fichier leja_sequence.dat
      vector<double> m_lejaSequence;

    public:
      /**
        * Constructeur
        * \param f : fonction à interpoler
        * \param nIter : nombre d'itérations dans l'algorithme AI
      */
      LagrangeInterpolation(FunctionsPtr f, int nIter);
      ~LagrangeInterpolation() {};

/******************************************************************************/
/************************ Points d'interpolation ******************************/
      /**
        * Retourne le point multivarié correspondant au multi-indice nu
        * \param nu : mult-indice (indice du point 1d de Leja sur chaque direction)
        * \return point multivarié correspondant à nu
      */
      MultiVariatePoint<double> getPoint(MultiVariatePointPtr<int> nu);
      /**
        * Ajoute le nouveau point d'interpolation p dans m_interpolationPoints et m_interpolationNodes
        * \param p : nouveau point d'interpolation
      */
      void addInterpolationPoint(MultiVariatePoint<double> p);

/******************************************************************************/
/*************************** Algorithme AI ************************************/
      /**
        * Retourner le premier multi-indice dans l'algo AI
        * \return le premier multi-indice dans l'algo AI = (0,..,0) correspont au point multivarié (1.0, .., 1.0)
      */
      MultiVariatePointPtr<int> getFirstMultivariatePoint();
      /**
        * Retourner le vecteur d'ordres du prochain point d'interpolation sélectionné par l'algo AI
        * 3 fois sur 4: on choisi celui ayant la plus grande erreur (algo AI)
        * 1 fois sur 4: on choisi celui qui a attendu le plus longtemps pour éviter de bloquer une direction dans le cas où l'erreur vaut 0
        * \param iteration : numéro de l'itération courante
        * \return teration : vecteur d'ordres du prochain point d'interpolation
      */
      MultiVariatePointPtr<int> maxElement(int iteration);
      /**
        * Vérifier si le multi-indice index existe dans m_path (correspond à un point d'interpolation)
        * \param index multi-indice
        * \return true si et seulement si index est dans m_path
      */
      bool indiceInPath(MultiVariatePoint<int> index);
      /**
        * Mettre à jour, en fonction du nouveau point d'interpolation correspondant au multi-indice nu, la liste des candidats courant
        * \param nu nouveau multi-indices choisi par l'algo AI
      */
      void updateCurentNeighbours(MultiVariatePointPtr<int> nu);
      /**
        * Vérifier si nu est bien un voisin de la séquence de points d'interpolation courante (vérifier la monotonie de l'ensemble de points choisis)
        * \param nu vecteur d'ordre
        * \return true si et seulement si nu est voisin à l'ensemble de points de m_path
      */
      bool isCorrectNeighbourToCurentPath(MultiVariatePointPtr<int> nu);

/******************************************************************************/
/************************** Fonctions de base *********************************/
      /**
        * Implémente la fonction de base (polynome de Lagrange global)
        * \param code : indice correspondant à une coordonnée d'un point d'interpolation
        * \param t : point d'évaluation (1d)
        * \param axis : direction concidéré (0 <= axis <= d-1)
        * \return la valeur au point t de la fonction de base correspondant au point 1d d'indice code sur la direction axis
      */
      double basisFunction_1D(int code, double t, int axis);
};

typedef std::unique_ptr<LagrangeInterpolation> LagrangeInterpolationPtr;

#endif
