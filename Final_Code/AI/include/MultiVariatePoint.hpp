#ifndef MULTIVARIATEPOINT
#define MULTIVARIATEPOINT

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <memory>
#include <limits>
#include <sstream>

#include "BinaryTree.hpp"

using namespace std;

/**
 *  \file MultiVariatePoint.hpp
 *  \brief Classe générique qui implémente un point multivarié
 *  \author SIALA Rafik
 *  \date 08/16
*/

// La classe est générique car le contenu du point multivarié peut être :
//    - des double si on modélise des points cartésiens
//    - des indices si on modélise des multi-indices (qui correspondent aux ordre des points cartésien dans
//    la version 0 (voir README, ligne 4) de l'algorithme AI)
//    - des chaine de caractères si on modélise des code de Huffman (qui correspondent aux ordre des points cartésien dans
//    la version 1 ou 2 (voir README, ligne 4) de l'algorithme AI)
template <typename T>
class MultiVariatePoint
{
    private:
      	// m_d : taille du point multivarié = dimension de l'éspace de départ de f
      	// m_n : taille de alpha (erreur d'interpolation en un point multivarié dans l'algorithme AI) = dimension de l'éspace d'arrivé de f
      	// f: (x_1, .., x_{m_d}) --> (y_1, .., y_{m_n}) est l'interpolé
        int m_d, m_n;
	      // Les coordonées du point multivarié
        T* m_nu = NULL;
        // erreur d'interpolation en un point multivarié dans l'algorithme AI
        // il s'agit d'un vecteur car l'erreur est la difference entre f et f_tilde (approximation de f) qui sont à valeurs vectorielles
        vector<double> m_alpha;
      	// Valeur initiale arbitraire: permet de vérifier si l'erreur a déja été calculée au point en question (pour éviter les recalculs)
      	double m_initialAlpha;
      	// Nombre d'itération au cours desquelles ce point est un candidat pour devenir un nouveau point d'interplation
        int m_nbWaitingIter;

   public:

        ~MultiVariatePoint();
        MultiVariatePoint() : m_d(0), m_n(0) {};

        /**
          *  Constructeur
          *  \param d : dimension de l'éspace de départ de f, (taille du point multivarié)
          *  \param n : dimension de l'éspace d'arrivé de f, (taille de l'erreur au point multivarié)
          *  \param T : valeur par défaut des coordonnées du point multivariés
        */
        MultiVariatePoint(int d, int n, T val);
        /**
          *  Constructeur par copie
          *  \param nu : point multivarié
        */
        MultiVariatePoint(const MultiVariatePoint<T>& nu);

        /**
          *  Accesseur à l'attribut m_d
          *  \return valeur de m_d
        */
        int getD() const { return m_d; };
        /**
          *  Accesseur à l'attribut m_n
          *  \return valeur de m_n
        */
        int getN() const { return m_n; };

        /**
          *  Accesseur à l'attribut m_nbWaitingIter
          *  \return valeur de m_nbWaitingIter
        */
        int getNbWaitingIter() { return m_nbWaitingIter; };
        /**
          *  Incrément la valeur de m_nbWaitingIter
        */
        void incrNbWaitingIter() { m_nbWaitingIter++; };

        /**
          *  Accesseur à l'attribut m_alpha
          *  \return valeur de m_alpha
        */
        vector<double> getAlpha() const { return m_alpha; };
        /**
          *  Mutateur de l'attribut m_alpha
          *  \param alpha : nouvelle valeur de m_alpha
        */
        void setAlpha(vector<double> alpha);

        /**
          *  Vérifier si m_alpha a déjà été calculé (en le comparant à m_initialAlpha)
          *  \return True si et seulement m_alpha != m_initialAlpha
        */
        bool alphaAlreadyComputed();

        /**
          * Vérifier si m_alpha est nulle
          *  \return True si et seulement si m_alpha vaut (0,..,0)
        */
        bool alphaIsNull();
        void reinit();


/******************************************************************************/
/************************ Surcharge d'opérateurs ******************************/
        // Surcharge de l'opérateur d'acces
        T &operator()(int d) const { return m_nu[d]; };
        // Surcharge de l'opérateur +=
        MultiVariatePoint& operator+=(const MultiVariatePoint & v);
        // Surcharge de l'opérateur -=
        MultiVariatePoint& operator-=(const MultiVariatePoint & v);
        // Surcharge de l'opérateur *=
        MultiVariatePoint& operator*=(const double & r);
        // Surcharge de l'opérateur =
        MultiVariatePoint& operator=(MultiVariatePoint const &nu);
};

template <typename T>
using MultiVariatePointPtr = shared_ptr<MultiVariatePoint<T>>;

template <typename T>
MultiVariatePoint<T>::~MultiVariatePoint()
{
    if (m_nu != NULL)
        delete[] m_nu;
}

template <typename T>
MultiVariatePoint<T>::MultiVariatePoint(int d, int n, T val) : m_d(d), m_n(n)
{
    // Initialisation des coordonnées du point multivarié
    if (m_d != 0)
    {
        m_nu = new T[m_d];
        for (int i=0; i<m_d; i++)
            m_nu[i] = val;
    }
    // Initialisation la valeur par défaut de m_alpha
    m_initialAlpha = numeric_limits<int>::max();
    // Initialisation de m_alpha (= (m_initialAlpha, .., m_initialAlpha))
    m_alpha.resize(m_n,m_initialAlpha);
    // Initialisation de m_nbWaitingIter
    m_nbWaitingIter = 0;
}

template <typename T>
MultiVariatePoint<T>::MultiVariatePoint(const MultiVariatePoint<T>& nu)
{
    // Construction par copie des informations de nu
    m_d = nu.m_d;
    m_n = nu.m_n;
    if (m_d != 0)
    {
        m_nu = new T[m_d];
        for (int i=0; i<m_d; i++)
            m_nu[i] = nu.m_nu[i];
    }
    m_initialAlpha = numeric_limits<int>::max();
    m_alpha.resize(m_n,m_initialAlpha);
    m_nbWaitingIter = 0;
};

template <typename T>
void MultiVariatePoint<T>::setAlpha(vector<double> alpha)
{
    for (int i=0; i<m_n; i++)
        m_alpha[i] = alpha[i];
}

template <typename T>
bool MultiVariatePoint<T>::alphaAlreadyComputed()
{
    // On compare les elements de m_alpha à m_initialAlpha
    // Si une coordonnée à changé alors m_alpha a déjà été calculé
    for (int i=0; i<m_n; i++)
        if (m_alpha[i]!=m_initialAlpha)
            return true;
    return false;
}

template <typename T>
bool MultiVariatePoint<T>::alphaIsNull()
{
    // On compare les m_alpha au point multivarié (0, .., 0)
    for (int i=0; i<m_n; i++)
        if (m_alpha[i])
            return false;
    return true;
}

template <typename T>
void MultiVariatePoint<T>::reinit()
{
    for (int i=0; i<m_n; i++)
        m_alpha[i] = m_initialAlpha;
}


template <typename T>
MultiVariatePoint<T>& MultiVariatePoint<T>::operator=(const MultiVariatePoint<T> & nu)
{
    // Construction d'un point multivarié par copie
    if (this != &nu)
    {
        m_d = nu.getD();
        m_n = nu.getN();
        if (m_d != 0)
        {
            delete[] m_nu;
            m_nu = new T[m_d];
            m_nbWaitingIter = 0;
            m_alpha.resize(m_n,0);
            memcpy(m_nu,nu.m_nu,sizeof(T)*nu.getD());
        }
    }
    return *this;
}

template <typename T>
MultiVariatePoint<T>& MultiVariatePoint<T>::operator+=(const MultiVariatePoint<T> & nu)
{
    // Ajout de nu à un point multivarié
    for (int i=0; i<nu.getD(); i++)
        (*this)(i) += nu(i);
    return *this;
}

template <typename T>
MultiVariatePoint<T>& MultiVariatePoint<T>::operator-=(const MultiVariatePoint<T> & nu)
{
  // Soustraction de nu à un point multivarié
    for (int i=0; i<nu.getD(); i++)
        (*this)(i) -= nu(i);
    return *this;
}

template <typename T>
MultiVariatePoint<T>& MultiVariatePoint<T>::operator*=(const double & r)
{
    // Multiplication du scalaire r par un point multivarié
    for (int i=0; i<getD(); i++)
        (*this)(i) *= r;
    return *this;
}


template <typename T>
MultiVariatePoint<T> operator+(const MultiVariatePoint<T> & p1 ,const MultiVariatePoint<T> & p2)
{
    // Retourne la somme de deux points multivariés p1 et p2
    // Résultat = p1 + p2
    MultiVariatePoint<T> p(p1);
    for (int i=0;i<p.getD();i++)
        p(i) += p2(i);
    return p;
}

template <typename T>
MultiVariatePoint<T> operator-(const MultiVariatePoint<T> & p1 ,const MultiVariatePoint<T> & p2)
{
    // Retourne la difference entre deux points multivariés p1 et p2
    // Résultat = p1 - p2
    MultiVariatePoint<T> p(p1);
    for (int i=0;i<p.getD();i++)
        p(i) -= p2(i);
    return p;
}

template <typename T>
MultiVariatePoint<T> operator*(MultiVariatePoint<T> p1 , const double & r)
{
    // Retourne le produit d'un point multivarié p1 par un scalaire r
    // Résultat = r * p1 = (r*p1(0), .., r*p1(d-1))
    MultiVariatePoint<T> p(p1);
    for (int i=0;i<p.getD();i++)
        p(i) *= r;
    return p;
}

template <typename T>
std::ostream & operator<<(std::ostream &out, const MultiVariatePoint<T> &nu)
{
    // Surcharge de l'opérateur de redirection de flux
    // Utile pour l'affichage des points multivariés
    out << "(";
    for (int k=0; k<int(nu.getD()-1); k++)
        cout << nu(k) << ",";
    out << nu(nu.getD()-1) << ")";
    return out;
}

template <typename T>
bool operator==(MultiVariatePoint<T> const &nu1 , MultiVariatePoint<T> const &nu2)
{
    // Comparaison de 2 points multivariés nu1 et nu2
    // nu1 != nu2 s'ils n'ont pas la même taille
    if (nu1.getD() != nu2.getD()) return false;
    int i = 0;
    // On parcours les coordonnées de nu1 et nu2
    while (i < nu1.getD())
    {
        // Si on trouve une différence alors les points ne sont pas égaux
        if (nu1(i) != nu2(i)) return false ;
        i++ ;
    }
    // Sinon ils sont égaux
    return true ;
}

#endif
